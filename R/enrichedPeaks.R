setMethod("enrichedPeaks", signature(regions='RangedData', sample1='RangedData', sample2='missing'),
function(regions, sample1, sample2, minHeight=100, space, mc.cores=1) {
  enrichedPeaks(regions=regions, sample1=ranges(sample1), minHeight=minHeight, space=space, mc.cores=mc.cores)
}
)

setMethod("enrichedPeaks", signature(regions='RangedData', sample1='IRangesList', sample2='missing'),
function(regions, sample1, sample2, minHeight=100, space, mc.cores=1) {
  n <- names(regions)
  sample1 <- sample1[n]
  f <- function(idx) {
    ans <- vector("list",length(idx))
    for (j in 1:length(idx)) ans[[j]] <- enrichedPeaks(regions=regions,sample1=sample1[[idx[j]]],minHeight=minHeight,space=n[idx[j]])
    return(ans)
  }
  if (mc.cores>1) {
    if ('parallel' %in% loadedNamespaces()) {
      ans <- pvec(1:length(n), f, mc.cores=mc.cores)
    } else stop('parallel library has not been loaded!')
  } else {
    ans <- f(1:length(n))
  }
  ans <- do.call(rbind,ans)
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
  n <- names(regions)
  sample1 <- sample1[n]; sample2 <- sample2[n]
  f <- function(idx) {
    ans <- vector("list",length(idx))
    for (j in 1:length(idx)) ans[[j]] <- enrichedPeaks(regions=regions,sample1=sample1[[idx[j]]],sample2=sample2[[idx[j]]],minHeight=minHeight,space=n[idx[j]])
    return(ans)
  }
  if (mc.cores>1) {
    if ('parallel' %in% loadedNamespaces()) {
       ans <- pvec(1:length(n), f, mc.cores=mc.cores)
    } else stop('parallel library has not been loaded!')
  } else {
    ans <- f(1:length(n))
  }
  ans <- do.call(rbind,ans)
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
      cov1 <- c(cov1-window(cov2,start=1,end=m1),window(cov2,start=m1+1,end=m2))
    } else {
      cov1 <- c(window(cov1,start=1,end=m2)-cov2,window(cov1,start=m2+1,end=m1))
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

setMethod("enrichedPeaks", signature(regions='GRanges', sample1='missing', sample2='missing'),
  function(regions, sample1, sample2, minHeight=100, space, mc.cores=1) {
    regions <- as(regions,'RangedData')
    ans <- enrichedPeaks(regions=regions,minHeight=minHeight,space=space,mc.cores=mc.cores)
    ans <- as(ans,'GRanges')
    return(ans)
  }
)

setMethod("enrichedPeaks", signature(regions='GRanges', sample1='GRanges', sample2='missing'),
  function(regions, sample1, sample2, minHeight=100, space, mc.cores=1) {
    regions <- as(regions,'RangedData')
    sample1 <- as(sample1,'RangedData')
    ans <- enrichedPeaks(regions=regions,sample1=sample1,minHeight=minHeight,space=space,mc.cores=mc.cores)
    ans <- as(ans,'GRanges')
    return(ans)
  }
)

setMethod("enrichedPeaks", signature(regions='GRanges', sample1='GRanges', sample2='GRanges'),
  function(regions, sample1, sample2, minHeight=100, space, mc.cores=1) {
    regions <- as(regions,'RangedData')
    sample1 <- as(sample1,'RangedData')
    sample2 <- as(sample2,'RangedData')
    ans <- enrichedPeaks(regions,sample1=sample1,sample2=sample2,minHeight=minHeight,space=space,mc.cores=mc.cores)
    ans <- as(ans,'GRanges')
    return(ans)
  }
)
