setMethod("ssdCoverage", signature(x="IRangesList"), function(x, mc.cores=1) {
  covx <- coverage(x)
  #cvs6 <- sd(covx)/mean(covx)
  cvs6 <- weighted.mean(sd(covx),w=sapply(covx,length))
  cvs6 <- 1000*cvs6/sqrt(sum(sapply(x,length)))
  return(cvs6)
}
)

setMethod("ssdCoverage", signature(x='RangedData'), function(x, mc.cores=1) { ssdCoverage(ranges(x)) } )

setMethod("ssdCoverage", signature(x='RangedDataList'), function(x, mc.cores=1) {
  if (mc.cores>1) {
    if ('parallel' %in% loadedNamespaces()) {
      ans <- parallel::mclapply(as.list(x), ssdCoverage, mc.cores=mc.cores, mc.preschedule=FALSE)
    } else stop('parallel library has not been loaded!')
  } else {
    ans <- lapply(x, ssdCoverage)
  }
  return(unlist(ans))
}
)

setMethod("ssdCoverage", signature(x='GRanges'),
  function(x, mc.cores=1) {
    x <- as(x,'RangedData')
    ans <- ssdCoverage(x,mc.cores=mc.cores)
    return(ans)
  }
)

setMethod("ssdCoverage", signature(x='GRangesList'),
  function(x, mc.cores=1) {
    x <- RangedDataList(lapply(x,function(y) as(y,'RangedData')))
    ans <- ssdCoverage(x,mc.cores=mc.cores)
    return(ans)
  }
)
