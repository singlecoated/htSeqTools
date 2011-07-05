\name{giniCoverage}
\alias{giniCoverage}
\alias{giniCoverage-methods}
\alias{giniCoverage,RangedData,ANY,ANY,ANY,character,integer-method}
\alias{giniCoverage,RangedData,ANY,ANY,ANY,character,missing-method}
\alias{giniCoverage,RangedData,ANY,ANY,ANY,missing,integer-method}
\alias{giniCoverage,RangedData,ANY,ANY,ANY,missing,missing-method}
\alias{giniCoverage,RangedDataList,ANY,ANY,ANY,character,integer-method}
\alias{giniCoverage,RangedDataList,ANY,ANY,ANY,character,missing-method}
\alias{giniCoverage,RangedDataList,ANY,ANY,ANY,missing,integer-method}
\alias{giniCoverage,RangedDataList,ANY,ANY,ANY,missing,missing-method}
\title{
Compute Gini coefficient.
}
\description{
  Calculate Gini coefficient of High-throughput Sequencing aligned
  reads. The index provides a measure of "inequality" in 
  read coverage which can be used for quality control purposes (see details).
}
\section{Methods}{
\describe{

\item{\code{signature(sample = "RangedData", mc.cores = "ANY", mk.plot = "ANY", seqName = "ANY", species = "character", chrLengths = "integer")}}{
%%  ~~describe this method here~~
Analize a single RangeData object with 'chrLengths' used for simulations ('Species' is ignored).

}

\item{\code{signature(sample = "RangedData", mc.cores = "ANY", mk.plot = "ANY", seqName = "ANY", species = "character", chrLengths = "missing")}}{
%%  ~~describe this method here~~
Analize a single RangeData object with chromosome lengths for simulations taken from BSgenome 'species' (package must be installed).
}

\item{\code{signature(sample = "RangedData", mc.cores = "ANY", mk.plot = "ANY", seqName = "ANY", species = "missing", chrLengths = "integer")}}{
%%  ~~describe this method here~~
Analize a single RangeData object with 'chrLengths' used as chromosome lengths in simulations.

}
\item{\code{signature(sample = "RangedData", mc.cores = "ANY", mk.plot = "ANY", seqName = "ANY", species = "missing", chrLengths = "missing")}}{
%%  ~~describe this method here~~
Analize all RangeData objects from sample (RangedDataList) with hromosome lengths for simulations taken as the largest end position of reads in each chromosome of all samples.

}

\item{\code{signature(sample = "RangedDataList", mc.cores = "ANY", mk.plot = "ANY", seqName = "ANY", species = "character", chrLengths = "integer")}}{
%%  ~~describe this method here~~
 Analize all RangeData objects from sample (RangedDataList) with 'chrLengths' used as chromosome lengths in simulations ('Species' is ignored).
}

\item{\code{signature(sample = "RangedDataList", mc.cores = "ANY", mk.plot = "ANY", seqName = "ANY", species = "character", chrLengths = "missing")}}{
%%  ~~describe this method here~~
 Analize all RangeData objects from sample (RangedDataList) with chromosome lengths for simulations taken from BSgenome 'species' (package must be installed).
}

\item{\code{signature(sample = "RangedDataList", mc.cores = "ANY", mk.plot = "ANY", seqName = "ANY", species = "missing", chrLengths = "integer")}}{
%%  ~~describe this method here~~
 Analize all RangeData objects from sample (RangedDataList) with 'chrLengths' used as chromosome lengths in simulations.
}

\item{\code{signature(sample = "RangedDataList", mc.cores = "ANY", mk.plot = "ANY", seqName = "ANY", species = "missing", chrLengths = "missing")}}{
%%  ~~describe this method here~~
 Analize all RangeData objects from sample (RangedDataList) with chromosome lengths for simulations taken as the largest end position of reads in each chromosome of sample.
}
}}
\usage{
giniCoverage(sample, mc.cores = 1, mk.plot = FALSE, seqName = "missing", species="missing", chrLengths="missing")
}
\arguments{
\item{sample}{A RangedData or RangedDataList object}
\item{seqName}{ If sample is a RangedData, name of sequence to use in
  plots }
\item{mk.plot}{ Logical. If TRUE, logarithm of coverage values' histogram and Lorenz
  Curve plot are plotted. }
\item{mc.cores}{If \code{mc.cores} is greater than 1, computations are
  performed in parallel for each element in the \code{IRangesList} object.}
\item{chrLengths}{An integer array with lengths of chromosomes in \code{sample} for simluations of uniformily distributed reads.}
\item{species}{A \code{BSgenome} species to obtain chromosome lengths for simluations of uniformily distributed reads.}
}
\details{
  The Gini coefficient provides a measure of "inequality" in 
  read coverage.
  This can be used in any sequencing experiment where the goal is to
  find peaks, i.e. unusual accumulation of reads in some genomic
  regions. For instance, Chip-Seq etc.
  Typically these experiments will consist of samples of interest
  (e.g. immuno-precipitated) and controls. The samples of interest
  should exhibit higher peaks, whereas reads in the controls should show a more
  uniform distribution.
  Since the Gini coefficient can be seen as a measure of departure from
  uniformity,
  the coefficient should present smaller values in the control samples.
  Since the Gini coefficient depends on the number of reads per sample, a correction is performed by substracting the Gini index from a sample with uniformily distributed reads. 
}
\value{
  If \code{mk.plot==FALSE}, the Gini index and adjusted Gini index for each element in the \code{RangedDataList} or \code{RangedData}
  object.
  
  If \code{mk.plot==TRUE}, a plot is produced showing the
  logarithm of coverage values' histogram and Lorenz Curve plot.
}
\references{
  See the definition of the Gini coefficient and Lorenz curve at
  http://en.wikipedia.org/wiki/Gini_coefficient
}
\author{
Camille Stephan-Otto
}
\seealso{
\code{\link{ssdCoverage}} for another measure of inequality in coverage.
}
\examples{
set.seed(1)
peak1 <- round(rnorm(500,100,10))
peak1 <- RangedData(IRanges(peak1,peak1+38),space='chr1')
peak2 <- round(rnorm(500,200,10))
peak2 <- RangedData(IRanges(peak2,peak2+38),space='chr1')
ip <- rbind(peak1,peak2)
bg <- runif(1000,1,300)
bg <- RangedData(IRanges(bg,bg+38),space='chr1')

rdl <- RangedDataList(ip,bg)
ssdCoverage(rdl)
giniCoverage(rdl)
}
\keyword{ univar }
