\name{stdPeakLocation}
\alias{stdPeakLocation}
\alias{PeakLocation}
\alias{stdPeakLocation-methods}
\alias{stdPeakLocation,data.frame-method}
\alias{stdPeakLocation,RangedData-method}
\alias{stdPeakLocation,GRanges-method}
\alias{PeakLocation-methods}
\alias{PeakLocation,data.frame-method}
\alias{PeakLocation,RangedData-method}
\alias{PeakLocation,GRanges-method}
\alias{plotMeanCoverage}
\alias{plotMeanCoverage-methods}
\alias{plotMeanCoverage,RleList,RangedData-method}
\title{Peak density with respect to closest gene.}
\description{
\code{stdPeakLocation} plots the density of peaks with respect to the
genomic feature (e.g. gene) in
standardized gene coordinates so that genes with different lengths are
comparable.

\code{PeakLocation} produces the same plot in non-standardized coordinates (i.e. distances
are measured in base pairs).

\code{plotMeanCoverage} plots the mean coverage in a series of selected
genomic regions.
}
\section{Methods}{
Methods for \code{stdPeakLocation}, \code{PeakLocation}

\describe{
\item{\code{signature(x = "data.frame")}}{ The data frame should contain
  columns named \code{start} and \code{end} indicating the peak location,
  \code{txStart}, \code{txEnd} indicating transcription start/end of the
  closest gene and \code{strand} indicating the strand. }

\item{\code{signature(x = "RangedData")}}{ \code{start(x)} and
  \code{end(x)} indicate the peak location. \code{x} should contain
  variables \code{x[['txStart']]}, \code{x[['txEnd']]} indicating the
  transcription start/end of the closest gene and
  \code{x[['strand']]} indicating the strand.}
}

Methods for \code{plotMeanCoverage}

\describe{
\item{\code{signature(cover="RleList", x="RangedData")}}{ \code{cover}
  contains the coverage and \code{x} the genomic regions of interest. }
}

}
\usage{
stdPeakLocation(x, peakDistance=1000, startpos='start_position', endpos='end_position',
strand='strand', distance, main='', xlab='Distance relative to feature length', xaxt='n',
  xlim=c(-1,2), densityType="kernel", nbreaks=10, ...)

PeakLocation(x, peakDistance=1000, startpos='start_position', endpos='end_position',
strand='strand', distance, main='', xlab='Distance (bp)',
  densityType="kernel", breaks, ...)

plotMeanCoverage(cover, x, upstreambp=1000, downstreambp=5000,
startpos='start_position', endpos='end_position', normalize=FALSE,
smooth=FALSE, span=0.05, main='', xlab='(bp)', ylab='Average coverage', ...)
}
\arguments{
\item{x}{A \code{RangedData} or \code{data.frame} indicating peak start and end in
  \code{start} and \code{end}, and start and end of
  the closest genomic feature (e.g. gene) in \code{startpos} and
  \code{endpos}.}
\item{cover}{An \code{RleList} object containing the coverage, as
  returned by the function \code{coverage}.}
\item{peakDistance}{Peaks more than \code{peakDistance} bases upstream
  or more than \code{3*peakDistance} downstream of the closest feature are discarded.}
\item{startpos}{Name of the variable storing the start position of the closest genomic feature.}
\item{endpos}{Name of the variable storing the end position of the closest genomic feature.}
\item{strand}{Name of the variable storing the strand for the closest genomic feature.}
\item{distance}{Name of the variable indicating the distance between the
  peak and the closest genomic feature. If left missing the distance
  between the feature start and the mid-point of the peak is computed.}
\item{main}{Graphical parameter passed on to \code{plot}.}
\item{xlab}{Graphical parameter passed on to \code{plot}.}
\item{ylab}{Graphical parameter passed on to \code{plot}.}
\item{xaxt}{Graphical parameter passed on to \code{plot}.}
\item{xlim}{In \code{stdPeakLocation} the x-axis limit is set to
  \code{xlim*peakDistance}. }
\item{densityType}{If we eant a density plot or a histogram. Has to be
  one of "kernel" (for the density plot) or "hist" for the histogram.}
\item{nbreaks}{Number of breaks to be used. It will not be used if
  \code{densityType} is different from "hist".}
\item{breaks}{This parameter will be passed to the \code{hist} plotting
  function. It will not be used if \code{densityType} is different from
  "hist".}
\item{upstreambp}{Number of bp upstream of the TSS where the coverage
  should be computed}
\item{downstreambp}{Number of bp downstream of the TSS where the
  coverage should be computed}
\item{normalize}{When set to \code{TRUE} the average coverage in each
  position is divided by the average across all positions. This is
  useful when trying to super-impose data from several experiments that
  had different read coverage.}
\item{smooth}{If set to \code{TRUE}, the average coverage is smooth by
  calling \code{loess}.}
\item{span}{Parameter controlling smoothing, passed on to
  \code{loess}. Larger values indicate more smoothing.}
\item{...}{Further parameters passed on to \code{plot}.}
}
\value{
  This function produces a density plot.
}
\examples{
#Generate synthetic peaks
set.seed(1)
st <- runif(100,1,1000)
en <- st+runif(length(st),25,100)
peaks <- RangedData(IRanges(st,en),space='chr1')

#Assign distance to closest gene
#(typically one would call annotatePeakInBatch
#from package ChIPpeakAnno to do this)
peaks[['start_position']] <- start(peaks) + runif(nrow(peaks),-500,1000)
peaks[['end_position']] <- peaks[['start_position']] + 500
peaks[['distance']] <- peaks[['start_position']] - start(peaks)
peaks[['strand']] <- sample(c('+','-'),nrow(peaks),replace=TRUE)
PeakLocation(peaks,peakDistance=1000)
}
\keyword{ graphs }
