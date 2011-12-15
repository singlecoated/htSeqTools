############################################################################
## AUXILIARY ROUTINES
############################################################################

fillRleList <- function(z, l) {
  #Append 0's at the end of each elem in RleList z so that the resulting lengths are equal to l
  #Note: if l has more elements than z, new elements are added to z (all filled with 0's)
  ans <- vector("list",length(l))
  names(ans) <- names(l)
  for (i in 1:length(ans)) {
    n <- names(ans)[i]
    if (n %in% names(z)) {
      ans[[n]] <- c(z[[n]],Rle(0,l[n]-length(z[[n]])))
    } else {
      ans[[n]] <- Rle(0,l[n]-length(z[[n]]))
    }
  }
  return(RleList(ans))
}


#tabIslandCounts <- function(z, nislands) {
#  #Tabulate overlap counts as returned by findOverlaps
#  #Input: -z: overlaps (RleList); -nislands: named vector indicating the number of islands in each chromosome
#  #Note: - all chromosomes in nislands are included in the output (NULL elements are added to z if necessary)
#  #      - tables count frequencies for all islands, i.e. output length is equal to that indicated in nislands
#  n <- names(nislands)[!(names(nislands) %in% names(z))]
#  if (length(n)>0) {
#    znew <- vector("list",length(n))
#    names(znew) <- n
#    z <- c(as.list(z),znew)
#  } else {
#    z <- as.list(z)
#  }
#  mapply(function(z1,z2) if(z2>0) {table(factor(z1,levels=1:z2))}else{NULL}, z1=z[names(nislands)], z2=nislands)
#}


############################################################################
## MAIN ROUTINES
############################################################################

setMethod(islandCounts, signature=(x='RangedDataList'), function(x, minReads=10, mc.cores=1) {
  #Add chromosomes missing in some samples
  n <- lapply(x,names)
  alln <- unique(unlist(n))
  missx <- lapply(n, function(z) { miss <- alln[!(alln %in% z)]; RangedData(IRanges(integer(0),integer(0)),space=miss) })
  sel <- sapply(missx,length)>0
  if (any(sel)) x[sel] <- RangedDataList(mapply(function(z,zmiss) { c(z,zmiss) }, z=x[sel], zmiss=missx[sel]))
  x <- lapply(x,function(z) RangedData(ranges(z)[alln]))  #sort by chromosome name           
  #Overall coverage
  if (mc.cores>1) {
    if ('multicore' %in% loadedNamespaces()) {
      cov1 <- multicore::mclapply(as.list(x), coverage, mc.cores=mc.cores, mc.preschedule=FALSE)
    } else stop('multicore library has not been loaded!')
  } else {
    cov1 <- lapply(as.list(x), coverage)
  }
  l <- lapply(cov1, function(z) sapply(z, length))
  l <- tapply(unlist(l), INDEX=unlist(sapply(l,names)), FUN=max) #max chr lengths
  l <- l[alln]
  if (mc.cores>1) {
    if ('multicore' %in% loadedNamespaces()) {
      cov1 <- multicore::mclapply(cov1, fillRleList, l=l, mc.cores=mc.cores, mc.preschedule=FALSE)
    } else stop('multicore library has not been loaded!')
  } else cov1 <- lapply(cov1, fillRleList, l=l)
  txt <- paste(paste('cov1[[',1:length(cov1),']]',sep=''),collapse='+')
  cov1 <- eval(parse(text=txt))
  #Islands
  islands1 <- slice(cov1, lower=minReads, rangesOnly=TRUE)
  islands1 <- RangedData(islands1)
  #Count nb reads in each island
  if (nrow(islands1)>0) {
    if (mc.cores>1) {
      if ('multicore' %in% loadedNamespaces()) {
        counts <- multicore::mclapply(as.list(x),function(z) countOverlaps(query=islands1,subject=z,minoverlap=1,type='any'), mc.cores=mc.cores, mc.preschedule=FALSE)
      } else stop('multicore library has not been loaded!')
    } else {
      counts <- lapply(x,function(z) countOverlaps(islands1,z,minoverlap=1,type='any'))
    }
    counts <- do.call(cbind,lapply(counts, unlist))
  } else {
    counts <- matrix(nrow=0,ncol=length(x))
  }
  ans <- RangedData(ranges=ranges(islands1), values=counts)
  if (!is.null(names(x))) colnames(values(ans)) <- names(x) else colnames(values(ans)) <- paste('counts',1:length(x),sep='')
  return(ans)
}
)

setMethod(islandCounts, signature=(x='RangedData'), function(x, minReads=10, mc.cores=1) {
  islandCounts(RangedDataList(x), minReads=minReads, mc.cores=mc.cores)
}
)

setMethod(islandCounts, signature=(x='GRanges'),
  function(x, minReads=10, mc.cores=1) {
    x <- as(x,'RangedData')
    ans <- islandCounts(x,minReads=minReads,mc.cores=mc.cores)
    ans <- as(ans,'GRanges')
    return(ans)
  }
)

setMethod(islandCounts, signature=(x='GRangesList'),
  function(x, minReads=10, mc.cores=1) {
    x <- RangedDataList(lapply(x,function(y) as(y,'RangedData')))
    ans <- islandCounts(x,minReads=minReads,mc.cores=mc.cores)
    ans <- as(ans,'GRanges')
    return(ans)
  }
)
