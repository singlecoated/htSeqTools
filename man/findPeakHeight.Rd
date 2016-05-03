\name{findPeakHeight}
\alias{findPeakHeight}
\alias{findPeakHeight-methods}
\alias{findPeakHeight,RangedData,RangedData,RangedData-method}
\alias{findPeakHeight,RangedData,IRangesList,IRangesList-method}
\alias{findPeakHeight,RangedData,Rle,Rle-method}
\alias{findPeakHeight,GRanges,GRanges,GRanges-method}
\alias{plotminHeight}
\title{FDR and optimal minHeight value estimation for ChIP-Seq peak
  calling with enrichedPeaks.}
\description{
Uses False Discovery Rate to estimate the optimal value of the minHeight
parameter in the call to enrichedPeaks. False Discovery Rate is computed
by swapping IP and Input samples and calculating the ratio of 'false' peaks
identified for a given set of minHeight values.
}
\section{Methods}{
\describe{

\item{\code{signature(regions = "RangedData", sample1 = "RangedData",
    sample2 = "RangedData")}}{ \code{sample1} indicates the start/end of
    reads in sample 1 (IP Sample), and similarly for \code{sample2} (Control sample). Only the subset of
    \code{regions} indicated by the argument \code{space} will be used. }
  
\item{\code{signature(regions = "RangedData", sample1 = "IRangesList",
    sample2 = "IRangesList")}}{ \code{regions} contains the regions of
    interest, \code{sample1} and \code{sample2} the reads in sample 1
    and sample 2, respectively. \code{names(sample1)} and
    \code{names(sample2)} must correspond to the space names used in \code{regions}.}
  
  \item{\code{signature(regions = "RangedData", sample1 = "Rle",
    sample2 = "Rle")}}{ \code{regions} contains the regions of
    interest, \code{sample1} and \code{sample2} the coverage for the reads in sample 1
    and sample 2, respectively. \code{names(sample1)} and
    \code{names(sample2)} must correspond to the space names used in \code{regions}.}

}}
\usage{
findPeakHeight(regions, sample1, sample2, hmin=5, hmax=200, myfdr=0.01, gridSize=25, space, mc.cores=1)
plotminHeight(x,...)
}
\arguments{
  \item{regions}{\code{RangedData} indicating
    the regions in which we wish to find peaks.}
  \item{sample1}{\code{IRangesList} or \code{Rle} object containing
    start and end of sequences in sample 1 (IP sample) or their coverage
    respectively.}
  \item{sample2}{Same for sample 2 (Control sample).}
  \item{hmin}{Minimum minHeight value to be considered for FDR
    estimation. Defaults to 5.}
  \item{hmax}{Maximum minHeight value to be considered for FDR
    estimation. Max coverage difference between sample 1 and sample 2 will
    be also calculated, and the minimum of the two will be used as hmax. This
    is done to avoid a skewed distribution of minHeight values to test
    for FDR estimation. }
  \item{myfdr}{Desired FDR cut-off.}
  \item{gridSize}{Number of intermediate steps of minHeight threshold to
  consider between hmin and hmax. Since FDR and optimal minHeight
  estimation is done by actually performing peak calls, selecting a high
value for gridSize can come at a big computational cost. Default value
of 25 or close is recommended. }
  \item{space}{Character text giving the name of the space for the
    \code{RangedData} object. Only used if \code{sample1} and
    \code{sample2} are of class \code{Rle}, for
    \code{RangedData} and {IRangesList} this is set up automatically.}
  \item{mc.cores}{If \code{mc.cores}>1 computations for each element in the
    \code{IRangesList} objects are performed in parallel (using
  the \code{parallel} function from package
  \code{parallel}). Notice: this option launches as many parallel
  processes as there are elements in \code{x}, which can place strong
  demands on the processor and memory.}
\item{x}{An object of class \code{list} as returned by the call to \code{findPeakHeight}.}
\item{...}{Other graphical arguments passed to function \code{plot}.}
}
\value{
  Object of class \code{list} with slots \code{fdr, npeaks} indicating peak calling FDR and number of peaks identified in the IP
  sample for each considered minHeight value and \code{cut, opt} with
  the desired FDR cut-off and the correponding minHeight value to be
  used in the call to \code{enrichedPeaks}.
}
\seealso{
  \code{enrichedPeaks}
}
\examples{
set.seed(1)
st <- round(rnorm(1000,500,100))
strand <- rep(c('+','-'),each=500)
space <- rep('chr1',length(st))
sample1 <- RangedData(IRanges(st,st+38),strand=strand,space=space)
st <- round(runif(1000,1,1000))
sample2 <- RangedData(IRanges(st,st+38),strand=strand,space=space)

#Find enriched regions and call peaks
mappedreads <- c(sample1=nrow(sample1),sample2=nrow(sample2))
regions <- enrichedRegions(sample1,sample2,mappedreads=mappedreads,minReads=50)
minHeight <- findPeakHeight(regions,sample1=sample1,sample2=sample2)
plotminHeight(minHeight)
peaks <- enrichedPeaks(regions,sample1=sample1,sample2=sample2,minHeight=minHeight$opt)
peaks <- peaks[width(peaks)>10,]
peaks

#Compute coverage in peaks
cover <- coverage(sample1)
coverinpeaks <- regionsCoverage(chr=peaks$space,start=start(peaks),end=end(peaks),cover=cover)

#Evaluate coverage in regular grid and plot
#Can be helpful fo clustering of peak profiles
coveringrid <- gridCoverage(coverinpeaks)
coveringrid
plot(coveringrid)

#Standardize peak profiles dividing by max coverage
stdcoveringrid <- stdGrid(coveringrid, colname='maxCov')
stdcoveringrid
plot(stdcoveringrid)
}
\keyword{datasets}
