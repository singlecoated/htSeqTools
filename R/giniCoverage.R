setMethod("giniCoverage", signature(sample='RangedDataList', species='missing', chrLengths='missing'),
function(sample, mc.cores=1, mk.plot=FALSE, seqName="", species="", chrLengths, numSim=1) {
sample<-lapply(sample, function(x) x[unique(as.character(space(x)))])
  chrNames<-unique(unlist(lapply(sample, names)))

  chrLengths<-lapply(sample, function(x){
		if (mc.cores>1) require(parallel)
	   	if(mc.cores>1) {
		        if ('parallel' %in% loadedNamespaces()) {
				len<-unlist(parallel::mclapply(as.list(x), function(y) max(end(ranges(y))), mc.cores=mc.cores, mc.preschedule=FALSE))
				names(len)<-names(x)
				len
			} else stop('parallel library has not been loaded!')
		} else {
				len<-unlist(lapply(x, function(y) max(end(ranges(y)))))
				names(len)<-names(x)
				len
		}
	}
 )
  
  chrLengths<-lapply(chrNames, function(x) unname(unlist(lapply(chrLengths, function(y) ifelse(x %in% names(y), y[[x]], 0)))))
  chrLengths<-do.call('cbind', chrLengths)
	chrLengths<-apply(chrLengths, 2, max)
  names(chrLengths)<-chrNames
  mode(chrLengths)<-"integer"
  if (is.null(names(sample))) names(sample) <- paste('sample',1:length(sample),sep='')
 	giniList<-lapply(names(sample), function(x) giniCoverage(sample[[x]], seqName=x, mc.cores=mc.cores, mk.plot=mk.plot, chrLengths=chrLengths, numSim=numSim))
  	names(giniList) <- names(sample)
	giniList<-do.call('cbind', giniList)
  return(t(giniList))
}
)

setMethod("giniCoverage", signature(sample='RangedDataList', species='character', chrLengths='integer'),
function(sample, mc.cores=1, mk.plot=FALSE, seqName="", species="", chrLengths=1, numSim=1) {
	warning("Both \'species\' and \'chrLengths\' provided, chrLengths will be used")
        if (is.null(names(sample))) names(sample) <- paste('sample',1:length(sample),sep='')
	giniList<-lapply(names(sample), function(x) giniCoverage(sample[[x]], seqName=x, mc.cores=mc.cores, mk.plot=mk.plot, chrLengths=chrLengths, numSim=numSim))
	names(giniList)<-names(sample)
	giniList<-do.call('cbind', giniList)
	return(t(giniList))
}
)

setMethod("giniCoverage", signature(sample='RangedDataList', species='missing', chrLengths='integer'),
function(sample, mc.cores=1, mk.plot=FALSE, seqName="", species, chrLengths=1, numSim=1) {
        if (is.null(names(sample))) names(sample) <- paste('sample',1:length(sample),sep='')
	giniList<-lapply(names(sample), function(x) giniCoverage(sample[[x]], seqName=x, mc.cores=mc.cores, mk.plot=mk.plot, chrLengths=chrLengths, numSim=numSim))
	names(giniList)<-names(sample)
	giniList<-do.call('cbind', giniList)
	return(t(giniList))
}
)

setMethod("giniCoverage", signature(sample='RangedDataList', species='character', chrLengths='missing'),
function(sample, mc.cores=1, mk.plot=FALSE, seqName="", species="", chrLengths, numSim=1) {
	if(!(species %in% available.genomes())){ stop('\'species\' must be a valid BSgenome character string')}
	else library(species, character.only=TRUE, logical.return=TRUE)
	simpleSp<-sub("BSgenome.", "", species)
	simpleSp<-sub('.UCSC.+', "", simpleSp)
	chrLengths<-seqlengths(get(simpleSp))
        if (is.null(names(sample))) names(sample) <- paste('sample',1:length(sample),sep='')
	giniList<-lapply(names(sample), function(x) giniCoverage(sample[[x]], seqName=x, mc.cores=mc.cores, mk.plot=mk.plot, chrLengths=chrLengths, numSim=numSim))
  	names(giniList) <- names(sample)
	giniList<-do.call('cbind', giniList)
  	return(t(giniList))
}
)

setMethod("giniCoverage", signature(sample='RangedData', species='missing', chrLengths='missing'),
function(sample, mc.cores=1, mk.plot=FALSE, seqName="", species, chrLengths, numSim=1) {
 sample<-sample[unique(as.character(space(sample)))]
    chrLengths<-unlist(lapply(sample, function(y) max(end(ranges(y)))))
    names(chrLengths)<-names(sample)								           
    gini<- giniCoverage(sample, seqName=seqName, mc.cores=mc.cores, mk.plot=mk.plot, chrLengths=chrLengths, numSim=numSim)
    return(gini)
}
)

setMethod("giniCoverage", signature(sample='RangedData', species='character', chrLengths='integer'),
function(sample, mc.cores=1, mk.plot=FALSE, seqName="", species="", chrLengths=1, numSim=1) {
	warning("Both species and chrLengths provided, chrLengths will be used")
	gini<-giniCoverage(sample, mc.cores=mc.cores, mk.plot=mk.plot, seqName=seqName, chrLengths=chrLengths, numSim=numSim)
    return(gini)
}
)

setMethod("giniCoverage", signature(sample='RangedData', species='character', chrLengths='missing'),
function(sample, mc.cores=1, mk.plot=FALSE, seqName="", species="", chrLengths, numSim=1) {
    if(!(species %in% available.genomes())){ stop('\'species\' must be a valid BSgenome character string')}
    else library(species, character.only=TRUE, logical.return=TRUE)
    simpleSp<-sub("BSgenome.", "", species)
    simpleSp<-sub('.UCSC.+', "", simpleSp)
    chrLengths<-seqlengths(get(simpleSp))
	giniIdx<-giniCoverage(sample, seqName=seqName, mc.cores=mc.cores, mk.plot=mk.plot, chrLengths=chrLengths, numSim=numSim)
	return(giniIdx)
}
)

setMethod("giniCoverage", signature(sample='RangedData', species='missing', chrLengths='integer'),
function(sample, mc.cores=1, mk.plot=FALSE, seqName="", species, chrLengths=1, numSim=1) {

  lorenzC <- function (x, n = rep(1, length(x)), plot = FALSE) 
	{
	    k <- length(x)
	    o <- order(x)
	    x <- x[o]
	    n <- n[o]
	    x <- n * x
	    p <- cumsum(n)/sum(n)
	    L <- cumsum(x)/sum(x)
	    p <- c(0, p)
	    L <- c(0, L)
	    L2 <- L * mean(x)
	    Lc <- list(p, L, L2)
	    names(Lc) <- c("p", "L", "L.general")
	    class(Lc) <- "Lc"
	    if (plot)  {
            plot.Lc <- function (x, general = FALSE, lwd = 2, xlab = "p", ylab = "L(p)", 
                                        main = "Lorenz curve", las = 1, ...) 
        	{
            	if (!general) 
            	L <- x$L
            	else L <- x$L.general
            	plot(x$p, L, type = "l", main = main, lwd = lwd, xlab = xlab, 
                	 ylab = ylab, xaxs = "i", yaxs = "i", las = las, ...)
            	abline(0, max(L))
        	}
	    	plot.Lc(Lc)
        }
	    Lc
	}

  ginifun<-function(x, n) {
	lorenz<-lorenzC(x,n)
	gini<-1-sum(diff(lorenz$p) * (lorenz$L[-1] + lorenz$L[-length(lorenz$L)]))
	return(gini)
   }

  sortFun <- function(x, decreasing = FALSE, na.last = NA, ...)
          {
              ord <- order(runValue(x))	
              Rle(values = runValue(x)[ord], lengths = runLength(x)[ord])
          }


  tableFun<-function(x){
	x<-sortFun(x)
	structure(array(runLength(x), dim=nrun(x), dimnames=structure(list(as.character(runValue(x))), names="")), class="table")
      }


  generateCounts<-function(sample, mc.cores){
   if (mc.cores>1) require(parallel)
   if(mc.cores>1) {
     if ('parallel' %in% loadedNamespaces()) {
       cove<-parallel::mclapply(as.list(ranges(sample)), coverage, mc.cores=mc.cores, mc.preschedule=FALSE)
     } else stop('parallel library has not been loaded!')
     #if (!mk.plot) cat("Done with coverage. ")
   } else {
     cove<-coverage(ranges(sample))
     #if (!mk.plot) cat("Done with coverage. ")
   }
     counts<-lapply(as.list(cove), tableFun)
   counts
 }

  merge.chromo.table<-function(x){
     xnames<-unique(unlist(lapply(x, names)))
     xtable<-vector(length=length(xnames), mode="numeric")
     names(xtable)<-xnames
                #cat(paste("Done sampling ", length(reads), " reads from chromosome ", names(x)[i], " of length ", chrLen[i], "\n", sep=""))
     for(i in xnames){
       for(j in names(x)){
         if(i %in% names(x[[j]])){
           xtable[i]<-xtable[i]+x[[j]][i]
         }
       }
     }
     return(xtable)
   }


 #Simulate data with same number of reads and genome length
   sampleRange<-function(x, chrLen, mc.cores){
	#chrReads<-unlist(lapply(x, nrow))
	#chrReads[chrReads==0]<-1
	chrLen<-chrLen[names(x)]
	#totReads<-sum(unlist(chrReads))
	totReads<-nrow(x)
	chrReads<-totReads*(chrLen/sum(as.numeric(chrLen)))
        if (mc.cores>1) require(parallel)
        if(mc.cores>1) {
          if ('parallel' %in% loadedNamespaces()) {
            rangesl<-parallel::mclapply(1:length(chrReads), function(i) {
                reads<-sample.int(chrLen[i], as.integer(chrReads[i]), replace=T)
		len<-floor(mean(width(x[1:min(nrow(x), 10000, chrLen[i]),])))
                ranges<-IRanges(start=reads, width=rep(len, length(reads)))
                ranges
        }, mc.cores=mc.cores
	)
          }
          } else{
             rangesl<-lapply(1:length(chrReads), function(i) {
                reads<-sample.int(chrLen[i], as.integer(chrReads[i]), replace=T)
                len<-floor(mean(width(x[1:min(nrow(x), 10000, chrLen[i]),])))
                ranges<-IRanges(start=reads, width=rep(len, length(reads)))
                ranges
              })
           }
   RD<-RangedData(RangesList(rangesl))
   RD   
      }

  if(numSim>1) {
    cat(paste("Simulating uniformily distributed data ", numSim, " times\n", sep=""))
    simRange<-lapply(1:numSim, function(x) sampleRange(sample, chrLengths, mc.cores))
    simCounts<-lapply(simRange, generateCounts, mc.cores)
    simCounts<-lapply(simCounts, merge.chromo.table)
    simGini<-lapply(simCounts, function(x) ginifun(as.numeric(names(x)), x))
    cat(paste("Average simulated gini: ", mean(unlist(simGini)), "\nStandard deviation: ", sd(unlist(simGini)), "\n\n", sep="")) 
    simGini<-mean(unlist(simGini))
    } else {
      cat("Simulating uniformily distributed data\n")
      simRange<-sampleRange(sample, chrLengths, mc.cores=mc.cores)
      simCounts<-generateCounts(simRange, mc.cores)
      simCounts<-merge.chromo.table(simCounts)
      simGini<-ginifun(as.numeric(names(simCounts)), simCounts)
  }
    cat("Calculating gini index of original data\n")
 
  counts<-generateCounts(sample, mc.cores)
  counts<-merge.chromo.table(counts)
  gini=ginifun(as.numeric(names(counts)), counts)
  giniIdx<-list(gini=gini, gini.adjust=gini-simGini, counts=counts)


    if(mk.plot){
       
     plotRes<-function(x, seqName, lcurve){
       if(lcurve) par(mfrow=c(2,1))
       if (seqName!="") {
         main=paste("sample: ", seqName, "\nGini index: ",round(x[['gini.adjust']],4), sep="")
       } else {
         main=paste("Gini index: ",round(x[['gini.adjust']],4), sep="") 
       }
         plot(log(x[['counts']]/sum(x[['counts']])), type="h", main=main, ylab='Proportion of bases (log)')
       if(lcurve) {

           lorenzC(as.numeric(names(x[['counts']])), x[['counts']], plot=TRUE)
           par(mfrow=c(1,1))
        } 
     }
     plotRes(giniIdx, seqName, lcurve=T)
   } else {
       return(unlist(giniIdx[c('gini', 'gini.adjust')]))
   }
 }
)

setMethod("giniCoverage", signature(sample='GRanges', species='character', chrLengths='integer'),
  function(sample, mc.cores=1, mk.plot=FALSE, seqName, species, chrLengths, numSim=1) {
    sample <- as(sample,'RangedData')
    ans <- giniCoverage(sample=sample,species=species,chrLengths=chrLengths,mc.cores=mc.cores,mk.plot=mk.plot,seqName=seqName,numSim=numSim)
    return(ans)
  }
)

setMethod("giniCoverage", signature(sample='GRanges', species='character', chrLengths='missing'),
  function(sample, mc.cores=1, mk.plot=FALSE, seqName, species, chrLengths, numSim=1) {
    sample <- as(sample,'RangedData')
    ans <- giniCoverage(sample=sample,species=species,mc.cores=mc.cores,mk.plot=mk.plot,seqName=seqName,numSim=numSim)
    return(ans)
  }
)

setMethod("giniCoverage", signature(sample='GRanges', species='missing', chrLengths='integer'),
  function(sample, mc.cores=1, mk.plot=FALSE, seqName, species, chrLengths, numSim=1) {
    sample <- as(sample,'RangedData')
    ans <- giniCoverage(sample=sample,chrLengths=chrLengths,mc.cores=mc.cores,mk.plot=mk.plot,seqName=seqName,numSim=numSim)
    return(ans)
  }
)

setMethod("giniCoverage", signature(sample='GRanges', species='missing', chrLengths='missing'),
  function(sample, mc.cores=1, mk.plot=FALSE, seqName, species, chrLengths, numSim=1) {
    sample <- as(sample,'RangedData')
    ans <- giniCoverage(sample=sample,mc.cores=mc.cores,mk.plot=mk.plot,seqName=seqName,numSim=numSim)
    return(ans)
  }
)

setMethod("giniCoverage", signature(sample='GRangesList', species='character', chrLengths='integer'),
  function(sample, mc.cores=1, mk.plot=FALSE, seqName, species, chrLengths, numSim=1) {
    sample <- RangedDataList(lapply(sample,function(y) as(y,'RangedData')))
    ans <- giniCoverage(sample=sample,species=species,chrLengths=chrLengths,mc.cores=mc.cores,mk.plot=mk.plot,seqName=seqName,numSim=numSim)
    return(ans)
  }
)

setMethod("giniCoverage", signature(sample='GRangesList', species='character', chrLengths='missing'),
  function(sample, mc.cores=1, mk.plot=FALSE, seqName, species, chrLengths, numSim=1) {
    sample <- RangedDataList(lapply(sample,function(y) as(y,'RangedData')))
    ans <- giniCoverage(sample=sample,species=species,mc.cores=mc.cores,mk.plot=mk.plot,seqName=seqName,numSim=numSim)
    return(ans)
  }
)

setMethod("giniCoverage", signature(sample='GRangesList', species='missing', chrLengths='integer'),
  function(sample, mc.cores=1, mk.plot=FALSE, seqName, species, chrLengths, numSim=1) {
    sample <- RangedDataList(lapply(sample,function(y) as(y,'RangedData')))
    ans <- giniCoverage(sample=sample,chrLengths=chrLengths,mc.cores=mc.cores,mk.plot=mk.plot,seqName=seqName,numSim=numSim)
    return(ans)
  }
)

setMethod("giniCoverage", signature(sample='GRangesList', species='missing', chrLengths='missing'),
  function(sample, mc.cores=1, mk.plot=FALSE, seqName, species, chrLengths, numSim=1) {
    sample <- RangedDataList(lapply(sample,function(y) as(y,'RangedData')))
    ans <- giniCoverage(sample=sample,mc.cores=mc.cores,mk.plot=mk.plot,seqName=seqName,numSim=numSim)
    return(ans)
  }
)
