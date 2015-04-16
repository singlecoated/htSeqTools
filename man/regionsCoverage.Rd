\name{regionsCoverage}
\alias{regionsCoverage}
\alias{gridCoverage}
\alias{stdGrid}
\alias{regionsCoverage-methods}
\alias{regionsCoverage,ANY,ANY,ANY,RleList-method}
\alias{stdGrid-methods}
\alias{stdGrid,gridCover-method}
\title{
Compute coverage on user specified genomic regions.
}
\description{
  \code{regionsCoverage} computes coverage for user specified genomic
  regions.
  
  \code{gridCoverage} evaluates the coverage on a regular grid with the
  same number of points for each region (facilitating further plotting,
  clustering etc).

  \code{stdGrid} standardized the coverage by diviging by the average or
  maximum coverage at each region.
}
\section{Methods}{
Methods for \code{regionsCoverage}:

\describe{
\item{\code{signature(chr = "ANY", start = "ANY", end = "ANY", cover = "RleList")}}{
Evaluates the coverage \code{cover} at the genomic positions specified by \code{chr},
\code{start}, \code{end}.
}}

Methods for \code{stdGrid}:

\describe{
\item{\code{signature(cover = "gridCover")}}{
Standardizes the coverage evaluated on a grid (typically, as returned by
\code{gridCoverage}) by dividing by the mean or maximum coverage.
}}

}
\usage{
regionsCoverage(chr, start, end, cover)

gridCoverage(cover)

stdGrid(cover, colname="maxCov")
}
\arguments{
  \item{chr}{Vector with chromosome names.}
  \item{start}{Vector with start position. \code{start>end} indicates that
    region is on the negative strand.}
  \item{end}{Vector with end position. \code{start>end} indicates that
    region is on the negative strand.}
  \item{cover}{For \code{regionsCoverage}, \code{cover} is an object of
    class \code{RleList} with the genome-wide coverage 
    (typically obtained by a previous call to \code{coverage}).
    For \code{gridCoverage} this is the coverage evaluated at user-specified regions, as returned
    by \code{regionsCoverage}. For \code{stdGrid} this is the coverage
    evaluated on a grid, as returned by \code{gridCoverage}.}
  \item{colname}{Name of the column in \code{cover@viewsInfo} to be used
  for the standardizing. Currently only \code{"meanCov"} and
  \code{"maxCov"} are implemented.}
}
\value{
  \code{regionsCoverage} returns a list with two components
  \item{views}{\code{RleViewsList} with coverage evaluated at specified
    regions. Orientation is always so that start<end, i.e. For most
    practical purposes, regions on
    the reverse strand will need to be inverted.}
  \item{viewsInfo}{\code{SplitDataFrameList} containing information
    about each peak (chromosome, strand, mean and maximum coverage).}

  \code{gridCoverage} and \code{stdGrid} return an object of class
  \code{gridCover}.
  The slot \code{cover} is a matrix with the coverage evaluated on a
  grid of 500 equi-spaced points, whereas the slot \code{viewsInfo} is the same as that returned
  by \code{regionsCoverage} (see above).
  For regions between 100bp and 500bp long, a linear interpolation is
  used to evaluate the coverage on the 500 points grid. For regions less
  than 100bp long, NAs are returned.
}
\seealso{
  \code{\link{gridCover-class}}
}
\examples{
#See help(enrichedPeaks)
}
\keyword{manip}
