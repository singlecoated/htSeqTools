\name{enrichedPeaks}
\alias{enrichedPeaks}
\alias{enrichedPeaks-methods}
\alias{enrichedPeaks,RangedData,IRanges,IRanges-method}
\alias{enrichedPeaks,RangedData,IRanges,missing-method}
\alias{enrichedPeaks,RangedData,IRangesList,IRangesList-method}
\alias{enrichedPeaks,RangedData,IRangesList,missing-method}
\alias{enrichedPeaks,RangedData,RangedData,missing-method}
\alias{enrichedPeaks,RangedData,RangedData,RangedData-method}
\title{Find peaks in sequencing experiments.}
\description{
Find peaks in 
significantly enriched regions found via \code{enrichedRegions}.
}
\section{Methods}{
\describe{

\item{\code{signature(regions = "RangedData", sample1 = "IRanges",
    sample2 = "IRanges")}}{ \code{sample1} indicates the start/end of
    reads in sample 1, and similarly for \code{sample2}. Only the subset of
    \code{regions} indicated by the argument \code{space} will be used. }

\item{\code{signature(regions = "RangedData", sample1 = "IRanges", sample2 = "missing")}}{\code{sample1} indicates the start/end of
    reads in sample 1, and similarly for \code{sample2}. Only the subset of
    \code{regions} indicated by the argument \code{space} will be used. }

\item{\code{signature(regions = "RangedData", sample1 = "IRangesList",
    sample2 = "IRangesList")}}{ \code{regions} contains the regions of
    interest, \code{sample1} and \code{sample2} the reads in sample 1
    and sample 2, respectively. \code{names(sample1)} and
    \code{names(sample2)} must correspond to the space names used in \code{regions}.}

\item{\code{signature(regions = "RangedData", sample1 = "IRangesList", sample2 = "missing")}}{ \code{regions} contains the regions of
    interest, \code{sample1} the reads in sample 1.
    \code{names(sample1)} must correspond to the space names used in \code{regions}. }

\item{\code{signature(regions = "RangedData", sample1 = "RangedData",
    sample2 = "missing")}}{ \code{space(sample1)} indicates the
    chromosome, and \code{start(sample1)} and
    \code{end(sample1)} indicate the start/end of the reads in sample 1.  }

\item{\code{signature(regions = "RangedData", sample1 = "RangedData", sample2 = "RangedData")}}{ \code{space(sample1)} indicates the
    chromosome, and \code{start(sample1)} and
    \code{end(sample1)} indicate the start/end of the reads in sample 1.
    Similarly for \code{sample2}.}
}}
\usage{
enrichedPeaks(regions, sample1, sample2, minHeight=100, space, mc.cores=1)
}
\arguments{
  \item{regions}{\code{RangedDataList} or \code{RangedData} indicating
    the regions in which we wish to find peaks.}
  \item{sample1}{\code{IRangesList} or \code{IRanges} object containing
    start and end of sequences in sample 1.}
  \item{sample2}{Same for sample 2. May be left missing, in which case
    only \code{sample1} is used to find peaks.}
  \item{minHeight}{If \code{sample2} is missing, peaks are defined as regions where the coverage in
    \code{sample1} is greater or equal than \code{minHeight}. If
    \code{sample2} is specified, the difference of coverage in sample1
    minus sample2 must be greater or equal than \code{minHeight}.}
  \item{space}{Character text giving the name of the space for the
    \code{RangedData} object. Only used if \code{sample1} and
    \code{sample2} are of class \code{RangedData}, for
    \code{RangedDataList} this is set up automatically.}
  \item{mc.cores}{If \code{mc.cores}>1 computations for each element in the
    \code{IRangesList} objects are performed in parallel (using
  the \code{parallel} function from package
  \code{multicore}). Notice: this option launches as many parallel
  processes as there are elements in \code{x}, which can place strong
  demands on the processor and memory.}
}
\value{
  Object of class \code{RangedData}
  indicating peaks higher than \code{minHeight}.
  Only peaks overlapping with \code{regions} are reported.
  The maximum of the coverage in each selected peak is reported in the
    column height (coverage in sample1 - sample2 when sample2 is
    specified).
  The column \code{region.pvalue} returns the p-value associated to the
    region that the peak belongs to (i.e. it is inherited from
    \code{regions}).
    Therefore, notice that all peaks corresponding to a single
    region will present the same \code{region.pvalue}.
}
\seealso{
  \code{enrichedRegions}
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
peaks <- enrichedPeaks(regions,sample1=sample1,sample2=sample2,minHeight=50)
peaks <- peaks[width(peaks)>10,]
peaks

#Compute coverage in peaks
cover <- coverage(sample1)
coverinpeaks <- regionsCoverage(chr=space(peaks),start=start(peaks),end=end(peaks),cover=cover)

#Evaluate coverage in regular grid and plot
#Can be helpful fo clustering of peak profiles
coveringrid <- gridCoverage(coverinpeaks)
coveringrid
plot(coveringrid)

#Standardize peak profiles dividing by max coverage
stdcoveringrid <- stdGrid(coveringrid, colname='maxCov')
stdcoveringrid
}
\keyword{datasets}
