setGeneric("findPeakHeight", function(regions, sample1, sample2, hmin=5, hmax=200, myfdr=0.01, gridSize=25, space, mc.cores=1) standardGeneric("findPeakHeight"))

setMethod("findPeakHeight", signature(regions='RangedData', sample1='RangedData', sample2='RangedData'),
function(regions, sample1, sample2, hmin=5, hmax=200, myfdr=0.01, gridSize=25, space, mc.cores=1) {
  findPeakHeight(regions,ranges(sample1),ranges(sample2),hmin,hmax,myfdr,gridSize,space,mc.cores)
})

setMethod("findPeakHeight", signature(regions='RangedData', sample1='IRangesList', sample2='IRangesList'),
function(regions, sample1, sample2, hmin=5, hmax=200, myfdr=0.01, gridSize=25, space, mc.cores=1) {
  # Function for Peak Calling FDR
  n <- names(regions)
  sample1 <- sample1[n]; sample2 <- sample2[n]
  # Computing global coverage and height
  cat('\nComputing coverage...')
  cov1 <- coverage(sample1)
  cov2 <- coverage(sample2)
  options(warn=-1)
  hmin <- max(hmin,mean(mean(cov1-cov2)))
  hmax <- min(hmax,mean(max(cov1-cov2)))
  options(warn=1)
  # Call peakCallFDR function to estimate best minHeight value
  f <- function(idx) {
    ans <- vector("list",length(idx))
    #for (j in 1:length(idx)) ans[[j]] <- findPeakHeight(regions=regions,sample1=sample1[[idx[j]]],sample2=sample2[[idx[j]]],hmin=hmin,hmax=hmax,myfdr=myfdr,space=n[idx[j]]) # With rangedDatas
    for (j in 1:length(idx)) ans[[j]] <- findPeakHeight(regions=regions,sample1=cov1[[idx[j]]],sample2=cov2[[idx[j]]],hmin=hmin,hmax=hmax,myfdr=myfdr,space=n[idx[j]]) # With coverages directly
    return(ans)
  }
  cat('\nEstimating peak calling threshold for FDR',myfdr,'...\n')
  if (mc.cores>1) {
    if ('parallel' %in% loadedNamespaces()) {
      fdr <- pvec(1:length(n), f, mc.cores=mc.cores)
    } else stop('parallel library has not been loaded!')
  } else {
    fdr <- f(1:length(n))
  }
  npeaks <- colSums(do.call(rbind,lapply(fdr,function(x) x$sel)),na.rm=TRUE)
  fdr.ans <- fdr <- colSums(do.call(rbind,lapply(fdr,function(x) x$sel.fdr)),na.rm=TRUE) / colSums(do.call(rbind,lapply(fdr,function(x) x$sel)),na.rm=TRUE)
  # No FDR can be higher than 1
  fdr[fdr>1 | is.na(fdr)] <- 1
  fdr.ans[fdr.ans>1 | is.na(fdr.ans)] <- 1
  # Transforming this fdr into isotonic to ensure decreasing
  isofdr <- isoreg(fdr)
  # plot(isofdr$yf,type='l')
  # Now use approxfun
  fdrfun <- tryCatch(expr=(approxfun(x=rev(as.numeric(names(fdr))),y=rev(isofdr$yf))),error=function(x) stop('\nError in approxfun (no peaks found?), consider using a lower hmin value\n'))
  fdrfun2 <- tryCatch(expr=(approxfun(y=rev(as.numeric(names(fdr))),x=rev(isofdr$yf))),error=function(x) stop('\nError in approxfun (no peaks found?), consider using a lower hmin value\n'))
  fdr <- fdrfun2(myfdr)
  return(list(fdr=fdr.ans,opt=fdr,cut=myfdr,npeaks=npeaks))
}
)

plotminHeight <- function(x,...)
  {
    fdr <- x$fdr
    isofdr <- isoreg(fdr)
    # Now use approxfun
    fdrfun <- approxfun(x=rev(as.numeric(names(fdr))),y=rev(isofdr$yf))
    fdrfun2 <- approxfun(y=rev(as.numeric(names(fdr))),x=rev(isofdr$yf))
    plot(fdrfun,xlim=c(0,max(as.numeric(names(fdr)))),xlab='minHeight',ylab='FDR / # of peaks',ylim=c(0,max(x$fdr)),...)
    abline(v=fdrfun2(x$cut),lty=2,col='red',...)
    par(new=TRUE)
    plot(rev(x$npeaks),type='l',xaxt='n',yaxt='n',lty=3,xlab='',ylab='',ylim=c(0,max(x$npeaks)))
    points(rev(x$npeaks),ylim=c(0,max(x$npeaks)))
    axis(4)
    legend('topright',legend=c(sprintf('FDR: %.3f',x$cut),sprintf('minHeight: %.2f',x$opt)))
    legend('bottomleft',legend=c('FDR','# of peaks'),lty=c(1,3))
  }

setMethod("findPeakHeight", signature(regions='RangedData', sample1='Rle', sample2='Rle'),
function(regions, sample1, sample2, hmin=5, hmax=200, myfdr=0.01, gridSize=25, space) {
  if (missing(space)) stop("Argument space must be specified for signature (regions='RangedData', sample1='Rle', sample2='Rle')")
  # Estimate optimum minHeight based on peak calling FDR
  peakCall <- function(cov1,cov2,minHeight)
    {
      m1 <- length(cov1); m2 <- length(cov2)
      if (m1<m2) {
        cov1 <- c(cov1-window(cov2,start=1,end=m1),window(cov2,start=m1+1,end=m2))
      } else {
        cov1 <- c(window(cov1,start=1,end=m2)-cov2,window(cov1,start=m2+1,end=m1))
      }
      islands1 <- slice(cov1, lower=minHeight)
      m <- viewMaxs(islands1)
      islands1 <- IRanges(start=start(islands1), end=end(islands1))
      counts1 <- findOverlaps(islands1, ranges(regions)[[space]], select='first')
      sel <- !is.na(counts1)
      sel
    } 
  # Straight peak calling, IP vs Input
  h <- seq(hmax,hmin,length.out=gridSize); names(h) <- as.character(round(h,2))
  sel <- unlist(lapply(h,function(x) sum(peakCall(sample1,sample2,minHeight=x))))
  # Reversed peak calling, Input vs IP
  sel.fdr <- unlist(lapply(h,function(x) sum(peakCall(sample2,sample1,minHeight=x))))
  # Compute FDR
  sel <- unlist(lapply(sel,function(x) sum(unlist(x))))
  sel.fdr <- unlist(lapply(sel.fdr,function(x) sum(unlist(x))))
  return(list(sel=sel,sel.fdr=sel.fdr))
})
