mergeRegionsCore <- function(start, end, chromosome, score, annot, aggregateFUN='median', maxDist=300) {
if (missing(start)) stop('Argument start must be specified')
if (missing(end)) stop('Argument end must be specified')
if (any(is.na(start))) stop('start contains NA values')
if (any(is.na(end))) stop('end contains NA values')
if (!missing(chromosome)) { 
  if (length(chromosome)!=length(start)) stop('chromosome and start must have the same length')
  if (length(chromosome)!=length(end)) stop('chromosome and end must have the same length')
}
if (!missing(annot)) { if (length(start)!=length(annot)) stop('start and annot must have the same length') }
if (!missing(score)) { if (length(start)!=length(score)) stop('start and score must have the same length') }
if (length(start)>0) {
  #sort data
  if (!missing(chromosome)) {
    o <- order(chromosome,start)
    chromosome <- chromosome[o]
  } else {
    o <- order(start)
  }
  start <- start[o]; end <- end[o]
  if (!missing(annot)) annot <- annot[o]
  if (!missing(score)) score <- score[o]
  #group neighbouring regions
  if (!missing(chromosome)) {
    myCh <- c(TRUE,(start[-1]-end[-length(start)])>maxDist | chromosome[-1]!=chromosome[-length(start)])
  } else {
    myCh <- c(TRUE,(start[-1]-end[-length(start)])>maxDist)
  }
  myGr <- cumsum(abs(myCh))
  s <- aggregate(x = start, by = list(myGr), FUN = "min")[,2]
  e <- aggregate(x = end, by = list(myGr), FUN = "max")[,2]
  xout <- data.frame(start=s,end=e)
  if (!missing(chromosome)) xout$chromosome <- as.character(chromosome[myCh])
  if (!missing(annot)) xout$annot <- as.character(annot[myCh])
  if (!missing(score)) xout$score <- aggregate(x = score, by = list(myGr), FUN = aggregateFUN)[,2]
  return(xout)
} else {
  return(NULL)
}
}


setMethod("mergeRegions",signature(intervals='IRanges'),
function(intervals, chromosome, score, annot, aggregateFUN='median', maxDist=300) {
  mergeRegionsCore(start=start(intervals), end=end(intervals), chromosome=chromosome, score=score, annot=annot, aggregateFUN=aggregateFUN, maxDist=maxDist)
}
)

setMethod("mergeRegions",signature(intervals='RleViews'),
function(intervals, chromosome, score, annot, aggregateFUN='median', maxDist=300) {
  mergeRegionsCore(start=start(intervals), end=end(intervals), chromosome=chromosome, score=score, annot=annot, aggregateFUN=aggregateFUN, maxDist=maxDist)
}
)

setMethod("mergeRegions",signature(intervals='data.frame'),
function(intervals, chromosome, score, annot, aggregateFUN='median', maxDist=300) {
  mergeRegionsCore(start=intervals$start, end=intervals$end, chromosome=chromosome, score=score, annot=annot, aggregateFUN=aggregateFUN, maxDist=maxDist)
}
)

setMethod("mergeRegions",signature(intervals='matrix'),
function(intervals, chromosome, score, annot, aggregateFUN='median', maxDist=300) {
  mergeRegionsCore(start=intervals[,'start'], end=intervals[,'end'], chromosome=chromosome, score=score, annot=annot, aggregateFUN=aggregateFUN, maxDist=maxDist)
}
)

setMethod("mergeRegions",signature(intervals='RangedData'),
function(intervals, chromosome, score, annot, aggregateFUN='median', maxDist=300) {
  if (!missing(score)) {
    if (!is.character(score)) stop('Argument score must be a character vector of length 1 indicating a variable name in values(intervals)')
    if (length(score)>1) stop("Argument score must be of length 1 for signature(intervals='RangedData')")
    myscore <- intervals[[score]]
    ans <- mergeRegionsCore(start=start(intervals), end=end(intervals), chromosome=space(intervals), score=myscore, aggregateFUN=aggregateFUN, maxDist=maxDist)
    eval(parse(text=paste("ans <- RangedData(IRanges(start=ans$start,end=ans$end), space=ans$chromosome,",score,"=ans$score)",sep='')))
  } else {
    ans <- mergeRegionsCore(start=start(intervals), end=end(intervals), chromosome=space(intervals), aggregateFUN=aggregateFUN, maxDist=maxDist)
    eval(parse(text=paste("ans <- RangedData(IRanges(start=ans$start,end=ans$end), space=ans$chromosome)",sep='')))
  }
  return(ans)
}
)
