\name{plotChrRegions}
\alias{plotChrRegions}
\title{Plot chromosomal regions of interest}
\description{
Produces a plot with all chromosomes for a given organism,
marking regions of interest in a user-specified color.
}
\usage{
plotChrRegions(regions, chrLength, markColor='red', ...)
}
\arguments{
\item{regions}{\code{RangedData} object with chromosome, start and end positions (chromosome must be stored in \code{space(regions)}.}
\item{chrLength}{Named integer vector with chromosome lengths in base pairs.}
\item{markColor}{Color to be used to mark the regions in the chromosome.}
\item{...}{Further parameters passed on to \code{plot}.}
}
\value{
  This function produces a plot.
}
\examples{
set.seed(1)
chr <- rep(c('chr1','chr2'),each=10)
chrLength <- c(chr1=10000,chr2=5000)
st <- c(runif(10,1,10000),runif(10,1,5000))
regions <- RangedData(IRanges(st,st+50),space=chr)

plotChrRegions(regions,chrLength=chrLength)
}
\keyword{ graphs }
