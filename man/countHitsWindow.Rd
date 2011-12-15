\name{countHitsWindow}
\alias{countHitsWindow}
\alias{countHitsWindow-methods}
\alias{countHitsWindow,RangedData-method}
\alias{countHitsWindow,GRanges-method}
\title{
Compute number of hits in a moving window along the chromosome.
}
\description{
Computes a smoothed number of hits along the chromosome by using moving
windows of user specified size.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "RangedData")}}{
\code{x} contains chromosome, start and end positions for each hit.
}
}}
\usage{
countHitsWindow(x, chrLength, windowSize = 10^4 - 1)
}
\arguments{
  \item{x}{Object containing hits (start, end and chromosome). Currently only \code{RangedData}
  objects are accepted.}
  \item{chrLength}{Named vector indicating the length of each chromosome
    in base pairs.}
  \item{windowSize}{Size of the window used to smooth the hit count.}
}
\examples{
set.seed(1)
st <- round(rnorm(1000,500,100))
st[st>=10000] <- 10000
strand <- rep(c('+','-'),each=500)
space <- rep('chr1',length(st))
x <- RangedData(IRanges(st,st+38),strand=strand,space=space)
countHitsWindow(x, chrLength=c(chr1=10000), windowSize=99)
}
\keyword{stats}
