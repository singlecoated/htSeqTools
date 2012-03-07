setMethod("alignPeaks",signature(x='RangedDataList', strand='character'),
  function(x, strand, npeaks=1000, bandwidth=150, mc.cores=1) {
    x <- as.list(x)
    ans <- multicore::mclapply(x, FUN=function(z) alignPeaks(z, strand=strand, npeaks=npeaks, bandwidth=bandwidth), mc.cores=mc.cores, mc.preschedule=FALSE)
    ans <- RangedDataList(ans)
    names(ans) <- names(x)
    return(ans)
  }
)

setMethod("alignPeaks",signature(x='RangedData',strand='character'),
  function(x, strand, npeaks=1000, bandwidth=150, mc.cores=1) {
    strandnew <- lapply(values(x),function(z) as.factor(z[,strand]))
    ans <- alignPeaks(x=ranges(x), strand=strandnew, npeaks=npeaks, bandwidth=bandwidth)
    ans <- RangedData(ans, values=values(x))
    colnames(ans) <- sub('values.','',colnames(ans))
    return(ans)
  }
)

setMethod("alignPeaks",signature(x='IRangesList',strand='list'),
  function(x, strand, npeaks=1000, bandwidth=150, mc.cores=1) {
    if (length(x) != length(strand)) stop('x and strand must have the same length')
    strand <- strand[names(x)]
    #Find highest peaks
    xcov <- vector("list",length(x)); names(xcov) <- names(x)
    xcov <- coverage(x)
    #for (i in 1:length(xcov)) xcov[[i]] <- coverage(x[[i]])
    xislands <- lapply(xcov,function(x) slice(x,lower=1))
    xmax <- lapply(xislands,function(x) viewMaxs(x))
    xwhichmax <- lapply(xislands,function(x) viewWhichMaxs(x))
    chr <- rep(names(xmax),lapply(xmax,function(x) length(x)))
    xmax <- unlist(xmax)
    xwhichmax <- unlist(xwhichmax)
    chr <- chr[order(xmax,decreasing=TRUE)][1:npeaks]
    xwhichmax <- xwhichmax[order(xmax,decreasing=TRUE)][1:npeaks]
    xmax <- xmax[order(xmax,decreasing=TRUE)][1:npeaks]
     
    #Format peaks +/- bandwidth as an IRangesList
    peakranges <- vector("list",length(x))
    names(peakranges) <- names(x)
    for (i in 1:length(peakranges)) {
      sel <- chr %in% names(peakranges)[i]
      if (any(sel)) peakranges[[i]] <- IRanges(start=xwhichmax[sel]-bandwidth,end=xwhichmax[sel]+bandwidth)
    }
    class(peakranges) <- 'IRangesList'
     
    #Find center of reads overlapping with peaks (+/- bandwidth)
    f <- function(x,strand,y) {
      if (!is.null(y)) {
        # o <- findOverlaps(y,query=x,multiple=TRUE)
        o <- findOverlaps(y,query=x, select='all') # Removed by Oscar on may 19. Check complained with unused argument(s) (multiple = TRUE)
        midpoint <- start(x)[queryHits(o)] - .5*(start(y)[subjectHits(o)]+end(y)[subjectHits(o)])
        # midpoint <- start(x)[as.matrix(o)[,'query']] - .5*(start(y)[as.matrix(o)[,'subject']]+end(y)[as.matrix(o)[,'subject']])
        strandSel <- strand[queryHits(o)]
        #strandSel <- strand[as.matrix(o)[,'query']]
      } else {
        midpoint <- strandSel <- NULL
      }
      return(list(midpoint=midpoint,strand=strandSel))
    }
    readct <- vector("list",length(x))
    for (i in 1:length(x)) readct[[i]] <- f(x=x[[i]],strand=strand[[i]],y=peakranges[[i]])
    #readct <- mapply(f,x,strand,peakranges)
    readstrand <- unlist(sapply(readct,function(x) as.character(x$strand)))
    readct <- unlist(sapply(readct,'[[','midpoint'))

    #Distance between peaks
    d <- mean(readct[readstrand=='-']) - mean(readct[readstrand=='+'])
    if (d<0) {
      d <- 0
      warning('The estimated shift size was below zero. Set to zero instead.')
    } else if (d>300) {
      d <- 300
      warning('The estimated shift size was > 300. Set to 300 instead.')
    }
    #Adjust reads
    if (d!=0) {
      adj <- ifelse(unlist(strand)=='+',d,-d)
      s <- unlist(start(x)) + adj
      e <- unlist(end(x)) + adj
      negs <- s<0
      e[negs] <- e[negs] - adj[negs]
      s[negs] <- s[negs] - adj[negs]
      space <- rep(names(x),sapply(x,length))
      x <- ranges(RangedData(IRanges(start=s,end=e),space=space))
      #
      #for (i in 1:length(x)) {
      #  sel <- ifelse(strand[[i]]=='+',1,-1)
      #  adj <- sel*d
      #  e <- end(x[[i]]) + adj
      #  s <- start(x[[i]]) + adj
      #  negs <- s<0
      #  e[negs] <- e[negs] - adj[negs]
      #  s[negs] <- s[negs] - adj[negs]
      #  x[[i]] <- IRanges(start=s,end=e)
      #}
    }
    cat('Estimated shift size is',d,'\n')
    return(x)
  }
)

setMethod("alignPeaks",signature(x='GRanges',strand='character'),
  function(x, strand, npeaks=1000, bandwidth=150, mc.cores=1) {
    x <- as(x,'RangedData')
    ans <- alignPeaks(x,strand=strand,npeaks=npeaks,bandwidth=bandwidth,mc.cores=mc.cores)
    ans <- as(ans,'GRanges')
    return(ans)
  }
)

setMethod("alignPeaks",signature(x='GRangesList'),
  function(x, strand, npeaks=1000, bandwidth=150, mc.cores=1) {
    x <- RangedDataList(lapply(x,function(y) as(y,'RangedData')))
    ans <- alignPeaks(x,strand=strand,npeaks=npeaks,bandwidth=bandwidth,mc.cores=mc.cores)
    ans <- as(ans,'GRangesList')
    return(ans)
  }
)

