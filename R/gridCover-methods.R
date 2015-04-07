########## SHOW METHOD ########

setMethod("show", signature(object="gridCover"), function(object) {
  cat("Object of class gridCover\n")
  cat("  - Slot cover: coverage for",nrow(object@cover),"regions\n")
  cat("  - Slot viewsInfo: SplitDataFrameList with information about each region\n")
}
)


########## PLOT ##########

setMethod("plot", signature(x="gridCover"), function(x, y, ...) {
  coverage <- colMeans(x@cover,na.rm=TRUE)
  plot(coverage, type='l', xaxt='n',xlab='',...)
  text(1,par('usr')[3],'Start',pos=3)
  text(length(coverage),par('usr')[3],'End', pos=3)
}
)

setMethod("lines", signature(x="gridCover"), function(x, ...) {
  coverage <- colMeans(x@cover)
  lines(coverage,...)
}
)


########## SUBSETTING ##########

setMethod("[",signature(x="gridCover"), function(x, i, j, ..., drop=FALSE) {
  x@cover <- x@cover[i,,drop=FALSE]
  x@viewsInfo <- x@viewsInfo[i,,drop=FALSE]
  return(x)
}
)



########## TRANSFORMATION #########

setGeneric("stdGrid", function(cover, colname='maxCov') standardGeneric("stdGrid"))
setMethod("stdGrid", signature(cover="gridCover"), function(cover, colname='maxCov') {
  m <- unlist(cover@viewsInfo[,colname])
  cover@cover <- cover@cover/m
  return(cover)
}
)



########## ACCESSOR for viewsInfo slot #########

setGeneric("getViewsInfo", function(x="gridCover") standardGeneric("getViewsInfo"))
setMethod("getViewsInfo", signature(x="gridCover"), function(x) {
  ans <- x@viewsInfo
  return(ans)
}
)



########## ACCESSOR for cover slot #########

setGeneric("getCover", function(x="gridCover") standardGeneric("getCover"))
setMethod("getCover", signature(x="gridCover"), function(x) {
  ans <- x@cover
  return(ans)
}
)
