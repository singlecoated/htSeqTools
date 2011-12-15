setMethod("filterDuplReads",signature(x='RangedData'),
  function(x, maxRepeats, fdrOverAmp=.01, negBinomUse=.999, components=0, mc.cores=1) {
    nrepeats <- countRepeats(x,mc.cores=mc.cores)
    seqid <- unlist(lapply(nrepeats,function(x) x[['readReps']]))
    nrepeats <- unlist(lapply(nrepeats,function(x) x[['reps']]))
    if (missing(maxRepeats)) {
      counts <- table(nrepeats)
      use <- 1:sum(cumsum(counts)/sum(counts)<negBinomUse)
      if (length(use)<=6) {
        components <- 1
      } else if (length(use)<=8 & components>2) {
        components <- 2
      } else if (length(use)<=12 & components==4) {
        components <- 3
      }
      fdr <- fdrEnrichedCounts(counts,use=use,components=components,mc.cores=mc.cores)$fdrEnriched
      maxRepeats <- max(5,match(TRUE,fdr<fdrOverAmp)-1)
      if (is.na(maxRepeats)) maxRepeats <- max(as.numeric(counts))+1
    }
    x <- x[seqid<= maxRepeats,]
    return(x)
  }
)

setMethod("filterDuplReads",signature(x='RangedDataList'),
  function(x, maxRepeats, fdrOverAmp=.01, negBinomUse=.999, components=0, mc.cores=1) {
    x <- as.list(x)
    mc.cores <- ifelse(missing(mc.cores),1,mc.cores)
    if (missing(maxRepeats)) {
      ans <- lapply(x,function(x) filterDuplReads(x,fdrOverAmp=fdrOverAmp,negBinomUse=negBinomUse,components=components,mc.cores=mc.cores))
    } else {
      ans <- lapply(x,function(x) filterDuplReads(x,maxRepeats=maxRepeats,fdrOverAmp=fdrOverAmp,negBinomUse=negBinomUse,components=components,mc.cores=mc.cores))
    }
    ans <- RangedDataList(ans)
    return(ans)
  }
)

setMethod("filterDuplReads",signature(x='GRanges'),
  function(x, maxRepeats, fdrOverAmp=.01, negBinomUse=.999, components=0, mc.cores=1) {
    x <- as(x,'RangedData')
    ans <- filterDuplReads(x,maxRepeats=maxRepeats,fdrOverAmp=fdrOverAmp,negBinomUse=negBinomUse,components=components,mc.cores=mc.cores)
    ans <- as(ans,'GRanges')
    return(ans)
  }
)

setMethod("filterDuplReads",signature(x='GRangesList'),
  function(x, maxRepeats, fdrOverAmp=.01, negBinomUse=.999, components=0, mc.cores=1) {
    x <- RangedDataList(lapply(x,function(y) as(y,'RangedData')))
    ans <- filterDuplReads(x,maxRepeats=maxRepeats,fdrOverAmp=fdrOverAmp,negBinomUse=negBinomUse,components=components,mc.cores=mc.cores)
    ans <- as(ans,'GRangesList')
    return(ans)
  }
)
