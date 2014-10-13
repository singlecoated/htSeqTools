setMethod("tabDuplReads", signature(x='RangedData'),
  function(x, minRepeats=1, mc.cores=1) {
    ans <- countRepeats(x, mc.cores=mc.cores)
    ans <- do.call(c,lapply(ans,'[[','reps'))
    table(ans)
  }
)

setMethod("tabDuplReads", signature(x='RangedDataList'),
  function(x, minRepeats=1, mc.cores=1) {
    x <- as.list(x)
    ans <- lapply(x, function(z) tabDuplReads(z, minRepeats=minRepeats, mc.cores=mc.cores))
    names(ans) <- names(x)
    return(ans)
  }
)

setMethod("tabDuplReads",signature(x='GRanges'),
  function(x, minRepeats=1, mc.cores=1) {
    x <- as(x,'RangedData')
    ans <- tabDuplReads(x,minRepeats=minRepeats,mc.cores=mc.cores)
    return(ans)
  }
)

setMethod("tabDuplReads",signature(x='GRangesList'),
  function(x, minRepeats=1, mc.cores=1) {
    x <- RangedDataList(lapply(x,function(y) as(y,'RangedData')))
    ans <- tabDuplReads(x,,minRepeats=minRepeats,mc.cores=mc.cores)
    return(ans)
  }
)
