\name{extendRanges}
\Rdversion{1.1}
\alias{extendRanges}
\alias{extendRanges-methods}
\alias{extendRanges,RangedData-method}
\alias{extendRanges,RangedDataList-method}
\title{
Extend reads or sequences by a user-specified number of bases.
}
\description{
  This function allows to extend ranges up to a user-specified length,
  which can be helpful in ChIP-seq analysis.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "RangedData")}}{ \code{space(x)} indicates the
chromosome, \code{start(x)} and \code{end(x)} the start/end positions of
each read. }

\item{\code{signature(x = "RangedDataList")}}{Each element in \code{x}
  is assumed to correspond to a different sample. }
}}
\usage{
extendRanges(x, seqLen = 200, chrlength, mc.cores=1)
}
\arguments{
  \item{x}{ Object containing reads. }
  \item{seqLen}{ Desired sequence length after extension. }
  \item{chrlength}{ Integer vector indicating the length of each
    chromosome. \code{names(chrlength)} must match those in
    \code{x}. This argument is used to ensure that no reads are extended
  beyond the maximum chromosome length.}
  \item{mc.cores}{ Number of cores to use in parallel computations
    (passed on to \code{mclapply}).}
}
\value{
A list of \code{IRanges} objects with extended sequence length.
}
\author{ David Rossell }
\examples{
set.seed(1)
st <- round(rnorm(1000,500,100))
st[st>2000] <- 2000
strand <- rep(c('+','-'),each=500)
space <- rep('chr1',length(st))
sample1 <- RangedData(IRanges(st,st+38),strand=strand,space=space)
extendRanges(sample1, seqLen=200, chrlength=c(chr1=2000))
}
\keyword{ manip }
