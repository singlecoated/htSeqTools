#Obtain coverage in specified regions. Returns as RleViewsList and viewsInfo (SplitDataFrameList) indicating chromosome and strand 
setMethod("regionsCoverage", signature(cover='RleList'), function(chr, start, end, cover) {
  #Select unique, non-missing regions
  sel <- !is.na(start) & !is.na(end)
  pos <- unique(data.frame(chr=chr[sel],start=start[sel],end=end[sel]))
  plusStrand <- pos$start<=pos$end
  #Determine start & end position
  st <- ifelse(plusStrand, pos$start, pos$end)
  en <- ifelse(plusStrand, pos$end, pos$start)
  st[st<1] <- 1
  maxlen <- sapply(cover,length)[as.character(pos$chr)]
  en[en>maxlen] <- maxlen[en>maxlen]
  #Obtain views
  myranges <- RangedData(ranges=IRanges(start=st, end=en), space=pos$chr)
  n <- names(myranges)[names(myranges) %in% names(cover)]
  myranges <- myranges[n]; cover <- cover[n]
  views <- Views(cover, ranges(myranges))
  selViews <- sapply(views,length)>0
  views <- views[selViews]
  n <- n[selViews]
  #Obtain viewsInfo
  viewsInfo <- data.frame(chr=rep(n,sapply(views,length)), strand=ifelse(plusStrand,'+','-'), meanCov=unlist(viewMeans(views)), maxCov=unlist(viewMaxs(views)))
  viewsInfo <- by(viewsInfo[,-1],viewsInfo$chr,DataFrame, simplify=FALSE)
  viewsInfo <- SplitDataFrameList(viewsInfo)
  #viewsInfo <- by(viewsInfo[,-1],viewsInfo$chr,FUN= function(z) z, simplify=FALSE)
  #viewsInfo <- SplitDataFrameList(do.call(list,viewsInfo))
  #colnames(viewsInfo) <- c('strand','meanCov','maxCov')
  ans <- list(views=views,viewsInfo=viewsInfo)
  return(ans)
}
)


gridCoverage <- function(cover) {
  strand <- lapply(cover$viewsInfo,function(z) z$strand)
  ans <- mapply(getGridList, as.list(cover$views), strand, SIMPLIFY=FALSE)
  ans <- do.call(rbind,ans)
  ans <- new("gridCover", cover=ans, viewsInfo=unlist(cover$viewsInfo))
  return(ans)
}

getGridList <- function(views, strand) {
  ans <- mapply(getGrid, as.list(views), as.list(strand), SIMPLIFY=FALSE)
  ans <- do.call(rbind,ans)
}

getGrid <- function(views, strand) {
  #Evaluate views on a grid of 500 values (100 values for the promoter region, 400 for the gene region)
  ans <- as.integer(if (strand=='+') views else rev(views))
  geneLength <- length(ans)
  if (geneLength>500) {
    geneGrid <- ans[round(seq(1,geneLength,length=500))]
  } else if (geneLength>100 & geneLength<=500) {
    f <- approxfun(x=seq(0,1,length=geneLength), y= ans)
    geneGrid <- f(seq(0,1,length=500))
  } else {
    geneGrid <- rep(NA,length=500)
  }
  return(geneGrid)
}
