setMethod("enrichedPeaks", signature(regions='RangedData', sample1='RangedData', sample2='missing'),
function(regions, sample1, sample2, minHeight=100, space, mc.cores=1) {
  enrichedPeaks(regions=regions, sample1=ranges(sample1), minHeight=minHeight, space=space, mc.cores=mc.cores)
}
)

setMethod("enrichedPeaks", signature(regions='RangedData', sample1='IRangesList', sample2='missing'),
function(regions, sample1, sample2, minHeight=100, space, mc.cores=1) {
  f <- function(y,z) { enrichedPeaks(regions=regions, sample1=y, minHeight=minHeight, space=z) }
  n <- names(regions)
  sample1 <- sample1[n]
  if (mc.cores>1) {
    if ('multicore' %in% loadedNamespaces()) {
      for (i in 1:length(n)) { multicore::parallel(f(sample1[[i]],n[i])) }
      ans <- multicore::collect()
    } else stop('multicore library has not been loaded!')
  } else {
    ans <- vector("list",length(n))
    for (i in 1:length(n)) { ans[[i]] <- f(sample1[[i]],n[i]) }
  }
  names(ans) <- rep('',length(ans))
  ans <- do.call(c,ans)
  return(ans)
}
)


setMethod("enrichedPeaks", signature(regions='RangedData', sample1='IRanges', sample2='missing'),
function(regions, sample1, sample2, minHeight=100, space, mc.cores=1) {
  if (missing(space)) stop("Argument space must be specified for signature (regions='RangedData', sample1='IRanges', sample2='missing')")
  if (space %in% names(regions)) {
    cov1 <- coverage(sample1)
    islands1 <- slice(cov1, lower=minHeight)
    m <- viewMaxs(islands1)
    islands1 <- IRanges(start=start(islands1), end=end(islands1))
    counts1 <- findOverlaps(islands1, ranges(regions)[[space]], select='first')
    sel <- !is.na(counts1)
    ans <- RangedData(islands1[sel], height=m[sel], region.pvalue=values(regions)[[space]][counts1[sel],'pvalue'], space=space)
  } else {
    ans <- RangedData(IRanges(start=integer(0),end=integer(0)), space=space)
  }
  return(ans)
}
)


setMethod("enrichedPeaks", signature(regions='RangedData', sample1='RangedData', sample2='RangedData'),
function(regions, sample1, sample2, minHeight=100, space, mc.cores=1) {
  enrichedPeaks(regions=regions, sample1=ranges(sample1), sample2=ranges(sample2), minHeight=minHeight, space=space, mc.cores=mc.cores)
}
)

setMethod("enrichedPeaks", signature(regions='RangedData', sample1='IRangesList', sample2='IRangesList'),
function(regions, sample1, sample2, minHeight=100, space, mc.cores=1) {
  f <- function(y1,y2,z) { enrichedPeaks(regions=regions, sample1=y1, sample2=y2, minHeight=minHeight, space=z) }
  n <- names(regions)
  sample1 <- sample1[n]; sample2 <- sample2[n]
  if (mc.cores>1) {
    if ('multicore' %in% loadedNamespaces()) {
      for (i in 1:length(n)) { multicore::parallel(f(sample1[[i]],sample2[[i]],n[i])) }
      ans <- multicore::collect()
    } else stop('multicore library has not been loaded!')
  } else {
    ans <- vector("list",length(n))
    for (i in 1:length(n)) { ans[[i]] <- f(sample1[[i]],sample2[[i]],n[i]) }
  }
  names(ans) <- rep('',length(ans))
  ans <- do.call(c,ans)
  return(ans)
}
)


setMethod("enrichedPeaks", signature(regions='RangedData', sample1='IRanges', sample2='IRanges'),
function(regions, sample1, sample2, minHeight=5, space, mc.cores=1) {
  if (missing(space)) stop("Argument space must be specified for signature (regions='RangedData', sample1='IRanges', sample2='missing')")
  if (space %in% names(regions)) {
    cov1 <- coverage(sample1)
    cov2 <- coverage(sample2)
    m1 <- length(cov1); m2 <- length(cov2)
    if (m1<m2) {
      cov1 <- c(cov1-seqselect(cov2,start=1,end=m1),seqselect(cov2,start=m1+1,end=m2))
    } else {
      cov1 <- c(seqselect(cov1,start=1,end=m2)-cov2,seqselect(cov1,start=m2+1,end=m1))
    }
    islands1 <- slice(cov1, lower=minHeight)
    m <- viewMaxs(islands1)
    islands1 <- IRanges(start=start(islands1), end=end(islands1))
    counts1 <- findOverlaps(islands1, ranges(regions)[[space]], select='first')
    sel <- !is.na(counts1)
    ans <- RangedData(islands1[sel], height=m[sel], region.pvalue=values(regions)[[space]][counts1[sel],'pvalue'], space=space)
  } else {
    ans <- RangedData(IRanges(start=integer(0),end=integer(0)), space=space)
  }
  return(ans)
}
)
