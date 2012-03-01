setMethod("plotMeanCoverage",signature(cover='RleList', x='RangedData'), function(cover, x, upstreambp=1000, downstreambp=5000, startpos='start_position', endpos='end_position', normalize=FALSE, main='', xlab='(bp)', ylab='Average coverage', ...) {
  sp <- x[['space']]; st <- x[[startpos]]; en <- x[[endpos]]
  #Remove missings
  sel <- !is.na(st) & !is.na(en)
  sp <- sp[sel]; st <- st[sel]; en <- en[sel]
  #Add upstream bases
  sel <- st<=en
  st[sel] <- st[sel] - upstreambp; st[!sel] <- st[!sel] + upstreambp
  l <- upstreambp+downstreambp
  fpos <- function(z) { if (length(z)<l) { ans <- c(as.vector(z),rep(0,l-length(z))) } else { ans <- as.vector(z)[1:l] } }
  fneg <- function(z) { z <- rev(z); fpos(z) }
  if (any(sel)) {
    rcovpos <- regionsCoverage(chr=sp[sel],start=st[sel],end=en[sel],cover=cover)
    matpos <- viewApply(rcovpos$views, FUN=fpos, simplify=TRUE)
    matpos <- matpos[sapply(matpos,length)>0]
    matpos <- do.call(cbind,as.list(matpos))
  } else {
    matpos <- matrix(ncol=0,nrow=upstreamdb+downstreamdb)
  }
  if (any(!sel)) {
    rcovneg <- regionsCoverage(chr=sp[!sel],start=st[!sel],end=en[!sel],cover=cover)
    matneg <- viewApply(rcovneg$views, FUN=fneg, simplify=TRUE)
    matneg <- matneg[sapply(matneg,length)>0]
    matneg <- do.call(cbind,as.list(matneg))
  } else {
    matneg <- matrix(ncol=1,nrow=upstreamdb+downstreamdb)
  }
  #Format coverage as matrix
  x2plot <- rowMeans(cbind(matpos,matneg))
  if (normalize) x2plot <- x2plot/mean(x2plot)
  plot(-upstreambp:(downstreambp-1),x2plot,type='l',main=main,xlab=xlab,ylab=ylab,...)
}
)
