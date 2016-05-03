setClass("cmdsFit", representation(points="matrix", d="matrix", dapprox="matrix", R.square="numeric"))

## SHOW
setMethod("show","cmdsFit",function(object) {
cat("Object of class cmdsFit approximating distances between", nrow(object@d), "objects \n")
cat("R-squared=",round(object@R.square,4), "\n")
}
)


## PLOT
setMethod("plot", signature(x="cmdsFit"), function(x, y, xlim, labels, cex.text=1, ...) {
  if (missing(xlim)) xlim <- 1.2*range(x@points)
  if (missing(labels)) labels <- rownames(x@points)
  if (ncol(x@points)==1) {
    plot(x@points, rep(0,2), xlim=xlim, ylim=xlim, xlab='', ylab='', ...)
    text(x@points, rep(0,2), labels, pos=3)
  } else {
    plot(x@points, xlim=xlim,ylim=xlim,xlab='',ylab='',...)
    text(x@points[,1],x@points[,2],labels,pos=3,cex=cex.text)
  }
}
)


