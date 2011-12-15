\name{cmds}
\alias{cmds}
\alias{cmds-methods}
\alias{cmds,RangedDataList-method}
\alias{cmds,GRangesList-method}
\title{Classical Multi-Dimensional Scaling}
\description{
\code{cmds} obtain the coordinates of the elements in \code{x} in a
\code{k} dimensional space
which best approximate the distances between objects.
For high-throughput sequencing data we define the distance between two
samples as 1 - correlation between their respective coverages.
This provides PCA analog for sequencing data.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "RangedDataList")}}{ Use Classical
  Multi-Dimensional Scaling to plot each element of the
  \code{RangedDataList} object in a k-dimensional space. The coverage is
  computed for each element in \code{x}, and the pairwise correlations
  between elements is used to define distances. }
}}
\usage{
cmds(x, k=2, logscale=TRUE, mc.cores=1, cor.method='pearson')
}
\arguments{
\item{x}{ A \code{RangedDataList} object, e.g. each element containing the
  output of a sequencing run.}
\item{k}{ Dimensionality of the reconstructed space, typically set to 2
  or 3.}
\item{logscale}{ If set to \code{TRUE} correlations are computed for \code{log(x+1)}.}
\item{mc.cores}{Number of cores. Setting \code{mc.cores>1} allows
  running computations in parallel. Setting \code{mc.cores} to too large
  a value may require a lot of memory.}
\item{cor.method}{ A character string indicating which correlation coefficient
          (or covariance) is to be computed.  One of "pearson"
          (default), "kendall", or "spearman", can be abbreviated.}
}
\value{
  The function returns a \code{mdsFit} object, with slots
  \code{points} containing the coordinates, \code{d} with the distances
  between elements, \code{dapprox} with the distances between objects in
  the approximated space, and \code{R.square} indicating the percentage
  of variability in \code{d} accounted for by \code{dapprox}.

  Since the coverage distribution is typically highly asymetric, setting
  \code{logscale=TRUE} reduces the influence of the highest coverage
  regions in the distance computation, as this is based on the Pearson
  correlation coefficient.
}
\examples{
data(htSample)
cmds1 <- cmds(htSample)

cmds1
plot(cmds1)
}
\keyword{ graphs }
