setMethod("stdPeakLocation", signature(x='data.frame'),
function(x, startpos='start_position', endpos='end_position', strand='strand', distance, main='', xlab='Distance relative to feature length', xaxt='n', xlim=c(-.25,1.5), densityType='kernel', nbreaks=10, ...) {
  if (!missing(distance)) distance <- x[,distance]
  stdPeakLocationBase(st= x$start, end=x$end, startpos= x[,startpos], endpos=x[,endpos], strand=x[,strand], distance=distance, main=main, xlab=xlab, xaxt=xaxt, xlim=xlim, densityType=densityType, nbreaks=nbreaks, ...)  
}
)

setMethod("stdPeakLocation", signature(x='RangedData'),
function(x, startpos='start_position', endpos='end_position', strand='strand', distance, main='', xlab='Distance relative to feature length', xaxt='n', xlim=c(-.25,1.5), densityType='kernel', nbreaks=10, ...) {
  if (!missing(distance)) distance <- x[[distance]]
  stdPeakLocationBase(st= start(x), end=end(x), startpos=x[[startpos]], endpos=x[[endpos]], strand=x[[strand]], distance=distance, main=main, xlab=xlab, xaxt=xaxt, xlim=xlim, densityType=densityType, nbreaks=nbreaks, ...)
}
)

dist2midpnt <- function(st,end,startpos,endpos,strand) {
  peakloc <- .5*(st+end)
  if ((is.character(strand)) | (is.factor(strand))) {
    distance <- ifelse(strand=='+', peakloc - startpos, endpos - peakloc)
  } else if (is.numeric(strand)) {
    distance <- ifelse(strand>=0, peakloc - startpos, endpos - peakloc)
  } else {
    stop('Invalid strand information')
  }
  return(distance)
}


stdPeakLocationBase <- function(st, end, startpos, endpos, strand, distance, main='', xlab='Distance relative to feature length', xaxt='n', xlim=c(-.25,1.5), densityType='kernel', nbreaks=10, ...) {
  if (missing(distance)) {
    distance <- dist2midpnt(st=st,end=end,startpos=startpos,endpos=endpos,strand=strand)
  }
  distance[(distance< -1*peakDistance) | (distance>3*peakDistance)] <- NA
  geneLength <- abs(endpos-startpos)
  distance <- distance/geneLength
  distance <- distance[!is.na(distance)]
  distance <- distance[(distance> 0) & (distance<1)]
  if (densityType=='kernel') {
    plot(density(distance),xlim=xlim,main=main,xlab=xlab,xaxt=xaxt,...)
  } else if (densityType=='hist') {
    hist(distance,xlim=xlim,main=main,xlab=xlab,xaxt=xaxt,breaks=nbreaks,...)
  } else {
    stop("densityType has to be 'kernel' or 'hist'")
  }
  y0 <- rep(par('usr')[3],4); y1 <- rep(par('usr')[3]+.025*diff(par('usr')[3:4]),4)
  segments(x0=-1:2,x1=-1:2,y0=y0,y1=y1)
  text(0:1,y1,c('Start','End'),pos=3)
}


setMethod("PeakLocation", signature(x='data.frame'),
function(x, peakDistance=10^4, startpos='start_position', endpos='end_position', strand='strand', distance, main='', xlab='Distance (bp)', densityType='kernel',breaks,...) {
  if (!missing(distance)) distance <- x[[distance]]
  PeakLocationBase(st= x$start, end=x$end, startpos= x[,startpos], endpos=x[,endpos], strand=x[,strand], distance=distance, peakDistance=peakDistance, main=main, xlab=xlab, densityType=densityType, breaks=breaks, ...)
}
)

setMethod("PeakLocation", signature(x='RangedData'),
function(x, peakDistance=10^4, startpos='start_position', endpos='end_position', strand='strand', distance, main='', xlab='Distance (bp)', densityType='kernel', breaks, ...) {
  if (!missing(distance)) distance <- x[[distance]]
  PeakLocationBase(st= start(x), end=end(x), startpos=x[[startpos]], endpos=x[[endpos]], strand=x[[strand]], distance=distance, peakDistance=peakDistance, main=main, xlab=xlab, densityType=densityType, breaks=breaks, ...)
}
)


PeakLocationBase <- function(st, end, startpos, endpos, strand, distance, peakDistance=10^4, main='', xlab='Distance (bp)',densityType='Kernel',breaks,...) {
  if (missing(distance)) {
    distance <- dist2midpnt(st=st,end=end,startpos=startpos,endpos=endpos,strand=strand)
  }
  distance <- distance[!is.na(distance)]
  distance <- distance[(distance> -1*peakDistance) & (distance<3*peakDistance)]
  if (densityType=='kernel') {
    plot(density(distance),main=main,xlab=xlab,...)
  } else if (densityType=='hist') {
    if (missing(breaks)) breaks <- seq(-peakDistance,max(distance),length=nclass.Sturges(distance))
    hist(distance,breaks=breaks,main=main,xlab=xlab,...)
  } else {
    stop("densityType has to be 'kernel' or 'hist'")
  }
  y0 <- rep(par('usr')[3],4); y1 <- rep(par('usr')[3]+.025*diff(par('usr')[3:4]),4)
  segments(x0=-1:2,x1=-1:2,y0=y0,y1=y1)
  text(0,y1,'Start')
}

setMethod("stdPeakLocation", signature(x='GRanges'),
  function(x, startpos='start_position', endpos='end_position', strand='strand', distance, main='', xlab='Distance relative to feature length', xaxt='n', xlim=c(-.25,1.5), densityType='kernel', nbreaks=10, ...) {
    x <- as(x,'RangedData')
    ans <- stdPeakLocation(x,startpos=startpos,endpos=endpos,strand=strand,distance=distance,main=main,xlab=xlab,xaxt=xaxt,xlim=xlim,densityType=densityType,nbreaks=nbreaks,...)
    return(ans)
  }
)

setMethod("PeakLocation", signature(x='GRanges'),
  function(x, peakDistance=10^4, startpos='start_position', endpos='end_position', strand='strand', distance, main='', xlab='Distance (bp)', densityType='kernel', breaks, ...) {
    x <- as(x,'RangedData')
    ans <- PeakLocation(x,peakDistance=peakDistance,startpos=startpos,endpos=endpos,strand=strand,distance=distance,main=main,xlab=xlab,densityType=densityType,breaks=breaks,...)
    return(ans)
  }
)
