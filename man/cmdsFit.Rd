\name{cmdsFit}
\alias{cmdsFit}
\alias{cmdsFit-methods}
\alias{cmdsFit,matrix-method}
\title{Classical Multi-Dimensional Scaling for a distance matrix}
\description{
\code{cmdsFit} obtains coordinates in a \code{k} dimensional space
which best approximate the given distances between objects.
}
\section{Methods}{
\describe{

\item{\code{signature(d = "matrix")}}{ Use Classical
  Multi-Dimensional Scaling to represent points 
  in a k-dimensional space.}

}}
\usage{
cmdsFit(d, k=2, type='classic')
}
\arguments{
\item{d}{ Distances between objects}
\item{k}{ Dimensionality of the reconstructed space, typically set to 2
  or 3.}
\item{type}{ Set to \code{"classic"} to perform classical MDS (uses
  function \code{cmdscale} from package \code{stats}). Set to
  \code{"isoMDS"} to use Kruskal's non-metric MDS (uses function
  \code{isoMDS} from package \code{MASS}). }
}
\value{
  The function returns a \code{cmdsFit} object.
  See help("cmdsFit-class") for details.
}
\examples{
d <- matrix(c(0,5,10,5,0,15,10,15,0),byrow=TRUE,ncol=3)
cmdsFit(d)
}
\keyword{ graphs }
