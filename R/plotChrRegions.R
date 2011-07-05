plotChrRegions <- function(regions, chrLength, markColor='red', ...) {
  if (!is(regions,'RangedData')) stop("Argument regions must be of class 'RangedData'")
  nchrom <- length(chrLength)
  ypos <- nchrom:1; names(ypos) <- names(chrLength)
  plot(NA,NA,xaxt='n',yaxt='n',xlim=c(-.15,1),ylim=c(0,nchrom+1),xlab='',ylab='',...)
  segments(x0=0,x1=chrLength/max(chrLength),y0=ypos,y1=ypos)
  text(0, ypos, names(ypos),pos=2)
  #
  chr <- as.character(space(regions))
  col <- markColor
  rect(xleft=start(regions)/max(chrLength),xright=end(regions)/max(chrLength),ybottom=ypos[chr]-.25,ytop=ypos[chr]+.25,col=col,border=col)
}
