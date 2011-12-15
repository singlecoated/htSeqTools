\name{ssdCoverage}
\alias{ssdCoverage}
\alias{ssdCoverage-methods}
\alias{ssdCoverage,IRangesList-method}
\alias{ssdCoverage,RangedData-method}
\alias{ssdCoverage,RangedDataList-method}
\alias{ssdCoverage,GRanges-method}
\alias{ssdCoverage,GRangesList-method}
\title{Standardized SD of the genomic coverage}
\description{
  Compute variability of the genomic coverage,
  measured as standardized SD per thousand sequences (see details).
  For instance, this can measure how pronounced are the peaks
  in a ChIP-Seq experiments,
  which can serve as a quality control to detect inefficient
  immuno-precipitation.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "IRangesList")}}{
A single coefficient of variation is returned,
as a weighted average of the coefficients of variation for each chromosome
(weighted according to the chromosome length).
}

\item{\code{signature(x = "RangedData")}}{
The method for \code{IRangesList} is used on \code{ranges(x)}.
}

\item{\code{signature(x = "RangedDataList")}}{
A vector with coefficients of variation for each element in \code{x} are
returned, by repeatedly calling the method for \code{RangedData} objects.
Use \code{mc.cores} to speed up computations with \code{mclapply},
but be careful as this requires more memory.
}
}}
\usage{
ssdCoverage(x, mc.cores=1)
}
\arguments{
\item{x}{Object with ranges indicating the start and end of each read. Currently, \code{x} can be of class
  \code{RangedDataList}, \code{RangedData} and \code{IRangesList}.}
\item{mc.cores}{Set \code{mc.cores} to a value greater than 1 to perform
computations in parallel, using package \code{mclapply}.}
}
\value{
  Numeric vector with coefficients of variation.
}
\details{
  ssdCoverage first computes the coverage for each sample and computes
  the standard deviation (SD) of the coverage.
  However, SD is not an appropriate measure of coverage unevenness, as
  its expected value is proportional to sqrt(n), where n is the number of reads
  (this can be seen with simple algebra).
 
  ssdCoverage therefore reports 1000*SD/sqrt(n), which can be interpreted
  as the standardized SD per thousand sequences.
  
}
\examples{
set.seed(1)
#Simulate IP data
peak1 <- round(rnorm(500,100,10))
peak1 <- RangedData(IRanges(peak1,peak1+38),space='chr1')
peak2 <- round(rnorm(500,200,10))
peak2 <- RangedData(IRanges(peak2,peak2+38),space='chr1')
ip <- rbind(peak1,peak2)

#Generate uniform background
bg <- runif(1000,1,300)
bg <- RangedData(IRanges(bg,bg+38),space='chr1')

rdl <- RangedDataList(ip,bg)
ssdCoverage(rdl)
giniCoverage(rdl)
}
\keyword{ univar }
