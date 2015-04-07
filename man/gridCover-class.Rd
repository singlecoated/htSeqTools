\name{gridCover-class}
\Rdversion{1.1}
\docType{class}
\alias{gridCover-class}
\alias{[,gridCover-method}
\alias{[,gridCover-method,ANY-method}
\alias{[,gridCover-method,ANY-method,ANY-method}
\alias{lines,gridCover-method}
\alias{plot,gridCover,ANY-method}
\alias{show,gridCover-method}
\alias{getViewsInfo,gridCover-method}
\alias{getCover,gridCover-method}
\alias{getViewsInfo}
\alias{getCover}

\title{Class "gridCover"}
\description{
Objects of class \code{gridCover} store coverage information evaluated
on a grid on pre-specified genomic regions.
}
\section{Objects from the Class}{
Objects of this class are returned by call to the function \code{gridCoverage}.
}
\section{Slots}{
  \describe{
    \item{\code{cover}:}{Object of class \code{"matrix"} with one row
      for each genomic region of interest, and 500 columns. Columns
      1-100 contain the coverage in the promoter region (as specified in
    argument \code{promoterDistance} to \code{gridCoverage}). Columns
    101-500 contain the coverage between \code{start} and \code{end} as
    indicated to \code{promoterDistance}.}
    \item{\code{viewsInfo}:}{Object of class \code{"DataFrame"} with
      information relative to each region (strand, mean and maximum coverage). }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "gridCover", i = "ANY", j = "ANY")}:
      Select a subset of peaks. }
    \item{plot}{\code{signature(x = "gridCover", y = "ANY")}: Plot the coverage. }
    \item{lines}{\code{signature(x = "gridCover")}: Add lines to an
      existing plot. }
    \item{show}{\code{signature(object = "gridCover")}: Show method. }
    \item{stdGrid}{\code{signature(x = "gridCover")}: Standardize the
      coverage by dividing by either the mean or the maximum coverage
      in each region. }
    \item{getViewsinfo}{\code{signature(x = "gridCover")}: Accessor for
      the \code{viewsInfo} slot.}
    \item{getCover}{\code{signature(x = "gridCover")}: Accessor for
      the \code{cover} slot.}
  }
}
\author{David Rossell}
\seealso{
\code{\link{regionsCoverage}} to compute coverage on pre-specified regions,
\code{\link{gridCoverage}} to compute coverage on a grid.
}
\examples{
##See help(gridCoverage)
}
\keyword{classes}
