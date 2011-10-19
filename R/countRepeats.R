setMethod(countRepeats, signature(reads='RangedData'), function(reads,mc.cores=1) {
  if (mc.cores>1) {
    if ('multicore' %in% loadedNamespaces()) {
      myfun <- function(idx) {
        ans <- vector('list',length(idx))
        for (i in 1:length(idx)) ans[[i]] <- countRepeats(ranges(reads)[[idx[i]]])
        return(ans)
      }
      ans <- multicore::pvec(1:length(reads),myfun,mc.cores=ifelse(length(reads)<=mc.cores,round(length(reads)/2),mc.cores))
    } else stop('multicore library has not been loaded!')
  } else {
    ans <- lapply(as.list(ranges(reads)),function(x) countRepeats(x))
  }
  ans
}
)

setMethod(countRepeats, signature(reads='IRangesList'), function(reads,mc.cores=1) {
  if (mc.cores>1) {
    myfun <- function(idx) {
      ans <- vector('list',length(idx))
      for (i in 1:length(idx)) ans[[i]] <- countRepeats(as.list(reads)[[idx[i]]])
      return(ans)
    }
    ans <- multicore::pvec(1:length(reads),myfun,mc.cores=ifelse(length(reads)<=mc.cores,length(reads)/2,mc.cores))
  } else {
    ans <- lapply(as.list(reads),countRepeats)
  }
  ans
}
)

setMethod(countRepeats, signature(reads='IRanges'), function(reads) {
  tdf = data.frame(oldOrder=1:length(reads),pos=start(reads), width=width(reads), strand="+")
  tdf = tdf[order(tdf$pos, tdf$width),]
  if(length(unique(tdf$width)) == 1)
  {
    reps = Rle(tdf$pos)@lengths
  } else {
    reps = Rle(paste(tdf$pos, tdf$width, sep="."))@lengths
  }
  cc = as.list(reps)
  cc[reps>1] = lapply(reps[reps>1], function(x) as.numeric(rep(x,x)))  
  readReps = unlist(cc)[order(tdf$oldOrder)]
  ans <- list(reps=reps,readReps=readReps)
  ans
}
)
