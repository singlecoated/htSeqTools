\name{coverageDiff}
\alias{coverageDiff}
\title{Compute the difference in coverage between two objects}
\description{
Computes coverage of sample1 minus coverage of sample2,
taking into account that the chromosomes in sample1 and sample2 are
not necessarily the same.
}
\usage{
coverageDiff(sample1, sample2, chrLength)
}
\arguments{
\item{sample1}{Object with reads from sample 1. Typically, a \code{RangedData} object.}
\item{sample2}{Object with reads from sample 2. Typically, a \code{RangedData} object.}
\item{chrLength}{Named vector with chromosome lengths. This can be obtained from the Bioconductor annotation packages,
e.g. \code{BSgenome.Dmelanogaster.UCSC.dm3} for drosophila melanogaster, etc.}
}
\details{
Computation is restricted to chromosomes in \code{names(chrLength)}.
}
\value{
\code{SimpleRleList} with differences in coverage.
}
\examples{
sample1 <- RangedData(IRanges(1:10,11:20),space='chr1')
sample2 <- RangedData(IRanges(1:10,11:20),space=rep(c('chr1','chr2'),each=5))
chrLength <- c(50,25); names(chrLength) <- c('chr1','chr2')
coverageDiff(sample1,sample2,chrLength)
}
\keyword{stats}
