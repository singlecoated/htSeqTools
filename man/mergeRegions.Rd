\name{mergeRegions}
\Rdversion{1.1}
\alias{mergeRegions}
\alias{mergeRegions-methods}
\alias{mergeRegions,data.frame-method}
\alias{mergeRegions,IRanges-method}
\alias{mergeRegions,matrix-method}
\alias{mergeRegions,RangedData-method}
\alias{mergeRegions,RleViews-method}
\alias{mergeRegions,GRanges-method}
\title{
Merge nearby chromosomal regions.
}
\description{
Merges regions that are less than \code{maxDist} bases apart.
}
\section{Methods}{
\describe{

\item{\code{signature(intervals = "data.frame")}}{ \code{intervals$start}
  and \code{intervals$end} give the interval start/end positions. }

\item{\code{signature(intervals = "IRanges")}}{ \code{start(intervals)}
  and \code{end(intervals)} give the interval start/end positions. }

\item{\code{signature(intervals = "matrix")}}{ The columns \code{start}
  and \code{end} in \code{intervals} give the interval start/end positions }

\item{\code{signature(intervals = "RangedData")}}{ \code{start(intervals)}
  and \code{end(intervals)} give the interval start/end positions. }

\item{\code{signature(intervals = "RleViews")}}{ \code{start(intervals)}
  and \code{end(intervals)} give the interval start/end positions. }
}}
\usage{
mergeRegions(intervals, chromosome, score, annot, aggregateFUN='median', maxDist=300)
}
\arguments{
  \item{intervals}{ Object indicating start and end of each region. It
    can either be a \code{matrix}, \code{data.frame}, \code{IRanges}, \code{RangedData} or
    an \code{RleViews} object. If a \code{matrix} or \code{data.frame},
    it must have columns named \code{start} and \code{end}.}
  \item{chromosome}{ Chromosome that the region belongs to (optional). If supplied, must be of the same length as \code{start} and \code{end}. }
  \item{score}{ Numerical score for each interval. Scores in merged
    intervals are aggregated using function \code{aggregateFUN}. If \code{intervals} is of class \code{RangedData}, this should be a character vector of length 1 indicating the name of the variable in \code{values(x)} containing the score.}
  \item{annot}{ Character indicating annotation information for each
    interval. Annotations in merged intervals are pasted in a single
    string (annotations appearing in more than one interval are only
    reported once in the merged interval).}
  \item{aggregateFUN}{Function to aggregate \code{score}.}
  \item{maxDist}{Regions less than \code{maxDist} apart are merged into a single region}
}
\value{
  The result is returned in a \code{data.frame} indicating the start and
  end of each merged interval. If the arguments were provided, the information in
  \code{chromosome}, \code{score} and \code{annot} is provided in
  additional columns.
  If the input argument \code{intervals} was of class \code{RangedData}, the results are returned in a \code{RangedData} object. 
}
\author{ David Rossell }
\examples{
st <- c(10,20,1000)
intervals <- RangedData(IRanges(st,st+10),space='chr1')

intervals
mergeRegions(intervals,maxDist=300)
}
\keyword{ manip }
