\name{islandCounts}
\alias{islandCounts}
\alias{islandCounts-methods}
\alias{islandCounts,RangedData-method}
\alias{islandCounts,RangedDataList-method}
\title{Find genomic regions with high coverage and count number of reads
overlapping each region in each sample}
\description{
  Finds genomics regions where the coverage is above a user-specified
  threshold
  and counts the number of ranges in each sample overlapping each region.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "RangedData")}}{
  \code{x} is assumed to contain the reads from a single sample.
  Genomic regions with high coverage will be detected and the number of
reads overlapping these regions will be computed.
}

\item{\code{signature(x = "RangedDataList")}}{
\code{x} is assumed to contain the reads for several samples, one sample
in each element of the list.
The overall coverage across all samples is computed by adding the
coverage in the individual samples, and the regions with overall
coverage above the user-specified threshold are selected.
Then the number of reads overlapping each region is computed.
}
}}
\usage{
islandCounts(x, minReads=10, mc.cores=1)
}
\arguments{
  \item{x}{\code{RangedData} or \code{RangedDataList} containing the
    reads. If a \code{RangedDataList} is provided, the overall coverage
    across all its elements is used to find the regions of interest, but
  individual counts are computed for each element in the list.}
  \item{minReads}{Only regions with coverage above \code{minReads} are
    considered.}
  \item{mc.cores}{If \code{mc.cores>1} computations are performed in
    parallel, using function \code{mclapply} from package \code{multicore}.}
}
\value{
  Object of class \code{RangedData}
  indicating the regions with coverage above \code{minReads} and the number of reads overlapping
  each sample for those regions.
}
\details{
The output of \code{islandCounts} can be the input data for a number of downstream analysis methods.
Although for a simple-minded analysis one could use \code{enrichedRegions},
one will usually want to use more sofisticated analyses (e.g. from packages \code{DEseq}, \code{BayesPeak}, \code{limma} etc.)
}
\examples{
set.seed(1)
st <- round(rnorm(1000,500,100))
strand <- rep(c('+','-'),each=500)
space <- rep('chr1',length(st))
sample1 <- RangedData(IRanges(st,st+38),strand=strand,space=space)
st <- round(rnorm(1000,1000,100))
sample2 <- RangedData(IRanges(st,st+38),strand=strand,space=space)

regions <- islandCounts(RangedDataList(sample1,sample2),minReads=50)
regions

#Plot coverage
plot(coverage(sample1)[[1]],type='l',xlim=c(0,2000))
lines(coverage(sample2)[[1]],col=2)
}
\keyword{manip}
