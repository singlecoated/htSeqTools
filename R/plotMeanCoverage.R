# introduced na.rm=TRUE in rowMeans on y2plot
setMethod("plotMeanCoverage",signature(cover='RleList', x='RangedData'), function(cover, x, upstreambp=1000, downstreambp=5000, startpos='start_position', endpos='end_position', normalize=FALSE, smooth=FALSE, span=0.05, main='', xlab='(bp)', ylab='Average coverage', ...) {
  sp <- x[['space']]; st <- x[[startpos]]; en <- x[[endpos]]
  #Remove missings
  sel <- !is.na(st) & !is.na(en)
  sp <- sp[sel]; st <- st[sel]; en <- en[sel]
  #Add upstream bases
  sel <- st<=en
  st[sel] <- st[sel] - upstreambp; st[!sel] <- st[!sel] + upstreambp
  l <- upstreambp+downstreambp
  fpos <- function(z) { if (length(z)<l) { ans <- c(as.vector(z),rep(0,l-length(z))) } else { ans <- as.vector(z)[1:l] } }
  fneg <- function(z) { z <- rev(z); fpos(z) }
  if (any(sel)) {
    rcovpos <- regionsCoverage(chr=sp[sel],start=st[sel],end=en[sel],cover=cover)
    matpos <- viewApply(rcovpos$views, FUN=fpos, simplify=TRUE)
    matpos <- as.list(matpos)
    matpos <- matpos[sapply(matpos,ncol)>0]
    matpos <- do.call(cbind,as.list(matpos))
  } else {
    matpos <- matrix(ncol=0,nrow=upstreambp+downstreambp)
  }
  if (any(!sel)) {
    rcovneg <- regionsCoverage(chr=sp[!sel],start=st[!sel],end=en[!sel],cover=cover)
    matneg <- viewApply(rcovneg$views, FUN=fneg, simplify=TRUE)
    matneg <- as.list(matneg)
    matneg <- matneg[sapply(matneg,length)>0]
    matneg <- do.call(cbind,as.list(matneg))
  } else {
    matneg <- matrix(ncol=1,nrow=upstreambp+downstreambp)
  }
  #Format coverage as matrix
  x2plot <- -upstreambp:(downstreambp-1)
  y2plot <- rowMeans(cbind(matpos,matneg),na.rm=TRUE)
  if (normalize) y2plot <- y2plot/mean(y2plot)
  if (smooth) {
    fit <- loess(y2plot~x2plot, span=span, degree=2)
    x2plot <- fit$x
    y2plot <- fit$fitted
  }                      
  plot(x2plot,y2plot,type='l',main=main,xlab=xlab,ylab=ylab,...)
}
)
