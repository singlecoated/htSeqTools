##
## This file contains code for enrichedChrRegions and countHitsWindow functions
##

## enrichedChrRegions

setMethod("enrichedChrRegions", signature(hits1='RangedData', hits2='missing'), function(hits1, hits2, chrLength, windowSize=10^4-1, fdr=0.05, nSims=10, mc.cores=1) {
  smoothCount <- countHitsWindow(hits1, chrLength=chrLength, windowSize=windowSize)  #smooth window
  cutseq <- 1:max(max(smoothCount))
  obshits <- hitsPerCut(smoothCount,cutseq=cutseq)
  rhits <- meanRandomHits(nSims=nSims, cutseq=cutseq, nhits=nrow(hits1), chrLength=chrLength, windowSize=windowSize,mc.cores=mc.cores)
  pobs <- rev(cumsum(rev(obshits)))/sum(obshits)
  pran <- rev(cumsum(rev(rhits)))/sum(rhits)
  fdrest <- pran/pobs
  cutoff <- cutseq[which(fdrest<fdr)[1]]
  if (!is.na(cutoff)) {
    regions <- slice(smoothCount,lower=cutoff,includeLower=TRUE,rangesOnly=TRUE)
    RangedData(regions)
  } else {
    RangedData(start=integer(0),end=integer(0))
  }
}
)


setMethod("enrichedChrRegions", signature(hits1='RangedData', hits2='RangedData'), function(hits1, hits2, chrLength, windowSize=10^4-1, fdr=0.05, nSims=10, mc.cores=1) {
  #Book keeping: ensure same chromosomes are in hits1, hits2
  hits1 <- RangedData(ranges(hits1)); hits2 <- RangedData(ranges(hits2))
  n <- names(hits1)[!names(hits1) %in% names(hits2)]
  if (length(n)>0) {
    temp <- IRangesList(lapply(as.list(n), function(z) IRanges(integer(0),integer(0))))
    names(temp) <- n
    hits2 <- RangedData(ranges=c(ranges(hits2),temp))
  }
  n <- names(hits2)[!names(hits2) %in% names(hits1)]
  if (length(n)>0) {
    temp <- IRangesList(lapply(as.list(n), function(z) IRanges(integer(0),integer(0))))
    names(temp) <- n
    hits1 <- RangedData(ranges=c(ranges(hits1),temp))
  }
  hits2 <- hits2[names(hits1)]
  # Smoothed number of hits
  smoothCount1 <- countHitsWindow(hits1, chrLength=chrLength, windowSize=windowSize)/nrow(hits1)  #smooth window
  smoothCount2 <- countHitsWindow(hits2, chrLength=chrLength, windowSize=windowSize)/nrow(hits2)
  countDif <- smoothCount1-smoothCount2
  cutseq <- seq(0,max(max(abs(countDif))),length=50)
  #obshits <- hitsPerCut2dir(countDif,cutseq=cutseq)
  obshits <- basesPerCut2dir(countDif,cutseq=cutseq)
  # FDR control
  rhits <- randomHits2groups(hits1=hits1, hits2=hits2, cutseq=cutseq, windowSize=windowSize, chrLength=chrLength)
  fdrest <- rhits/obshits; fdrest <- ifelse(fdrest>1,1,fdrest)
  cutoff <- cutseq[which(fdrest<fdr)[1]]
  if (!is.na(cutoff)) {
    regionsPos <- RangedData(slice(countDif,lower=cutoff,includeLower=TRUE),rangesOnly=TRUE)
    countDifNeg <- -1*countDif; names(countDifNeg) <- names(countDif)
    regionsNeg <- RangedData(slice(countDifNeg,lower=cutoff,includeLower=TRUE,rangesOnly=TRUE))
    regionsPos[['direction']] <- rep(1,nrow(regionsPos))
    regionsNeg[['direction']] <- rep(-1,nrow(regionsNeg))
    ans <- rbind(regionsPos,regionsNeg)
  } else {
    ans <- RangedData(IRanges(integer(0),integer(0)))
  }
  return(ans)
}
)

## countHitsWindow
#Count nb of hits (elements in x) in a window of specified size. Return RleList object with one elem for each chromosome
setMethod("countHitsWindow", signature(x="RangedData"), function(x, chrLength, windowSize=10^4-1) {
  if ((windowSize %% 2)==0) windowSize <- windowSize-1
  peaksSinglebase <- .5*(start(x) + end(x))
  peaksSinglebase <- RangedData(IRanges(peaksSinglebase,peaksSinglebase),space=space(x)) 
  peaksRle <- coverage(peaksSinglebase, width=as.list(chrLength[names(peaksSinglebase)]))
  peaksSmooth <- runsum(peaksRle, k=windowSize, endrule='constant')
  return(peaksSmooth)
}
)


#### AUXILIARY INTERNAL FUNCTIONS

countDifHitsWindow <- function(hits1, hits2, cutseq, windowSize, chrLength) {
  smoothCount1 <- countHitsWindow(hits1, chrLength=chrLength, windowSize=windowSize)/nrow(hits1)
  smoothCount2 <- countHitsWindow(hits2, chrLength=chrLength, windowSize=windowSize)/nrow(hits2)
  countDif <- smoothCount1-smoothCount2
  #obshits <- hitsPerCut2dir(countDif,cutseq=cutseq)
  obshits <- basesPerCut2dir(countDif,cutseq=cutseq)
  obshits
}

#Counts nb of regions with smoothed hit count (stored in hits) above threshold cutseq
hitsPerCut <- function(hits, cutseq) sapply(as.list(cutseq), function(z) sum(sapply(slice(hits, lower=z),length)))
#Count nb regions with absolute score (stored in hits) above threshold cutseq
hitsPerCut2dir <- function(hits, cutseq) sapply(as.list(cutseq), function(z) sum(sapply(slice(-1*hits,lower=z),length)) + sum(sapply(slice(hits,lower=z),length)))
#Count nb bases with absolute score (stored in hits) above threshold cutseq
basesPerCut2dir <- function(hits, cutseq) {
  sapply(as.list(cutseq), function(z) sum(sum(width(slice(abs(hits),lower=z)))))
  #sapply(as.list(cutseq), function(z) sum(sum(width(slice(-1*hits,lower=z))) + sum(width(slice(hits,lower=z)))))
}

randomHits2groups <- function(hits1, hits2, cutseq, windowSize, chrLength) {
  #Randomly scramble hits in hits1 and hits2 and compute number of regions with score above cutseq threshold
  hitid1 <- 1:nrow(hits1); hitid2 <- (nrow(hits1)+1):(nrow(hits1)+nrow(hits2))
  allhits <- c(hitid1,hitid2)
  sel1 <- allhits %in% sample(allhits,size=nrow(hits1),replace=FALSE) #randomize groups
  rhits1 <- rbind(hits1[sel1[hitid1],], hits2[sel1[hitid2],])
  rhits2 <- rbind(hits1[!sel1[hitid1],], hits2[!sel1[hitid2],])
  countDifHitsWindow(hits1=rhits1, hits2=rhits2, cutseq=cutseq, windowSize=windowSize, chrLength=chrLength) #count
}


meanRandomHits <- function(nSims=10, cutseq, nhits, chrLength, windowSize=10^4-1, mc.cores=1) {
  #Expected number of regions with smoothed hit count above threshold cutseq
  # - nSims: number of simulations to estimate the expectation
  # - cutseq: vector with threshold sequence (e.g. for cutseq=1:5 the expected nb of regions with more than 1,2,3,4 and 5 random hits in a window of size windowSize is computed)
  # - nhits: number of random hits to be simulated (usually set equal to number of hits in observed data)
  # - chrLength: 
  # - windowSize: window size
  # - mc.cores: number of cores to be used for parallel computation. Passed on to mclapply
  # Output: named vector with expected number of regions witih smoothed hit count above specified thresholds
  if (mc.cores>1) {
    if ('multicore' %in% loadedNamespaces()) {
      rhits <- multicore::mclapply(as.list(1:nSims), function(z) randomHitsWindow(nhits=nhits, chrLength=chrLength, windowSize=windowSize), mc.cores=mc.cores, mc.preschedule=FALSE)
      rhits <- matrix(unlist(multicore::mclapply(rhits, hitsPerCut, cutseq=cutseq, mc.cores=mc.cores, mc.preschedule=FALSE)),nrow=length(cutseq))
    } else stop('multicore library has not been loaded!')
  } else {
    rhits <- lapply(as.list(1:nSims), function(z) randomHitsWindow(nhits=nhits, chrLength=chrLength, windowSize=windowSize))
    rhits <- matrix(unlist(lapply(rhits, hitsPerCut, cutseq=cutseq)),nrow=length(cutseq))
  }
  ans <- rowMeans(rhits)
  names(ans) <- cutseq
  return(ans)
}

randomHitsWindow <- function(nhits, chrLength, windowSize=10^4-1) {
#Could nb of hits in a window (as in countHitsWindow) when hits are generated uniformly distributed
  probs <- chrLength/(length(chrLength)*mean(chrLength))
  chrRand <- as.vector(rmultinom(n=1, size=nhits, prob=probs))
  chrRand <- rep(names(chrLength),chrRand)
  posRand <- floor(runif(nhits,min=1,max=chrLength[chrRand]+1))
  xRand <- RangedData(IRanges(posRand,posRand),space=chrRand)
  covRand <- coverage(xRand, width=as.list(chrLength[names(xRand)]))
  covSmooth <- runsum(covRand, k=windowSize, endrule='constant')
  return(covSmooth)
}
