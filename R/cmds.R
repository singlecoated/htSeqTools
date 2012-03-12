setMethod("cmdsFit", signature=c(d='matrix'), function(d,k=2,type='classic',add=FALSE,cor.method='pearson') {
  if (type=='classic') {
    ans <- cmdscale(as.dist(d), k=k, add=add)
    if (add) ans <- ans$points
  } else if (type=='isoMDS') {
    ans <- isoMDS(d, k=k, maxit=100, trace=FALSE)$points
  }
  dapprox <- dist(ans, method='euclidean')
  dapprox <- as.matrix(dapprox)
  R.square <- ifelse(nrow(d)==2, 1, cor(d[upper.tri(d)], dapprox[upper.tri(dapprox)])^2)
  new("cmdsFit", points=ans, d=d, dapprox=dapprox, R.square=R.square)
}
)

setMethod("cmds", signature(x='RangedDataList'), function(x, k=2, logscale=TRUE, mc.cores=1, cor.method='pearson') {
  cat('Computing coverage...\n')
  if (mc.cores>1) {
    if ('multicore' %in% loadedNamespaces()) {
      cover <- multicore::mclapply(as.list(x), coverage, mc.cores=mc.cores, mc.preschedule=FALSE)
      if (logscale) cover <- multicore::mclapply(cover, function(z) log(z+1), mc.cores=mc.cores, mc.preschedule=FALSE)
    } else stop('multicore library has not been loaded!')
  } else {
    cover <- lapply(x, coverage)
    if (logscale) cover <- lapply(cover, function(z) log(z+1))
  }
  # Remove chrs without length
  cover <- lapply(cover, function(x) x[unlist(lapply(x,length))!=0])
  #
  cat('Computing correlations...\n')
  index <- expand.grid(1:length(cover),1:length(cover))
  index <- index[index[,1]<index[,2],]
  index <- as.list(data.frame(t(index)))
  if (mc.cores>1) {
    d <- multicore::mclapply(index, function(z) corRleList(cover[[z[1]]], cover[[z[2]]], cor.method=cor.method), mc.cores=mc.cores, mc.preschedule=FALSE)
  } else {
    d <- lapply(index, function(z) corRleList(cover[[z[1]]], cover[[z[2]]], cor.method=cor.method))
  }
  #
  r <- diag(length(cover))
  for (i in 1:length(index)) {
    r[index[[i]][1],index[[i]][2]] <- r[index[[i]][2],index[[i]][1]] <- d[[i]]
  }
  rownames(r) <- colnames(r) <- names(cover)
  d <- as.dist((1 - r)/2)
  if (k>=nrow(r)) k <- nrow(r)-1
  ans <- cmdscale(d, k=k)
  dapprox <- dist(ans, method='euclidean')
  #
  d <- as.matrix(d); dapprox <- as.matrix(dapprox)
  R.square <- ifelse(nrow(d)==2, 1, cor(d[upper.tri(d)], dapprox[upper.tri(dapprox)])^2)
  new("cmdsFit", points=ans, d=d, dapprox=dapprox, R.square=R.square)
}
)

corRleList <- function(z1, z2, cor.method='pearson') {
  n <- names(z1)[names(z1) %in% names(z2)]
  ans <- mapply(function(x,y) corRle(x,y,cor.method=cor.method), z1[n], z2[n])
  #Correlation is 0 for elements present only in z1 or only in z2
  notinz2 <- names(z1)[!names(z1) %in% names(z2)]
  if (length(notinz2)>0) {
    ans <- c(ans, rep(0,length(notinz2)))
    names(ans)[(length(ans)-length(notinz2)+1):length(ans)] <- notinz2
  }
  notinz1 <- names(z2)[!names(z2) %in% names(z1)]
  if (length(notinz1)>0) {
    ans <- c(ans, rep(0,length(notinz1)))
    names(ans)[(length(ans)-length(notinz1)+1):length(ans)] <- notinz1
  }
  #Weight correlations based on length
  l <- c(sapply(z1,length),sapply(z2,length))
  l <- tapply(l, FUN=max, INDEX=names(l))/max(l)  #divide by max length to avoid numerical overflow
  l <- l[names(ans)]
  sum(ans*l/sum(l))
}

corRle <- function(z1, z2, cor.method='pearson') {
  l1 <- length(z1); l2 <- length(z2)
  if (l1>l2) {
    z2 <- c(z2, Rle(0,l1-l2))
  } else {
    z1 <- c(z1, Rle(0,l2-l1))
  }
  if (cor.method=='spearman') {
    z1@values <- rank(z1@values)
    z2@values <- rank(z2@values)    
  }
  cor(z1,z2)
}

setMethod("cmds", signature(x='GRangesList'),
  function(x, k=2, logscale=TRUE, mc.cores=1, cor.method='pearson') {
    x <- RangedDataList(lapply(x,function(y) as(y,'RangedData')))
    ans <- cmds(x,k=k,logscale=logscale,mc.cores=mc.cores,cor.method=cor.method)
    #ans <- as(ans,'GRangesList')
    return(ans)
  }
)
