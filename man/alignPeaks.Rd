\name{alignPeaks}
\alias{alignPeaks}
\alias{alignPeaks-methods}
\alias{alignPeaks,IRangesList,list-method}
\alias{alignPeaks,RangedData,character-method}
\alias{alignPeaks,RangedDataList,character-method}
\alias{alignPeaks,GRanges,character-method}
\alias{alignPeaks,GRangesList-method}
\alias{alignPeaks,GRangesList,ANY-method}
\title{Align peaks in a ChIP-Seq experiment by removing the strand
  specific bias.}
\description{
Align peaks in a ChIP-Seq experiment by removing the shift between reads
aligned
to the plus and the minus strands.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "IRangesList", strand = "list")}}{ Each
  element in \code{x} corresponds to a chromosome, and each range gives
  the start/end of a sequence. \code{strand} indicates the strand for
  the ranges in \code{x}.}

\item{\code{signature(x = "RangedData", strand = "character")}}{
  \code{x} gives read start and end positions, and \code{strand} gives
  the name of the variable in \code{values(x)} containing the strand information. }

\item{\code{signature(x = "RangedDataList", strand = "character")}}{
  The method for \code{RangedData} is applied to each element in \code{x} separately, as each element may
  have a different strand-specific bias. }
}}
\usage{
alignPeaks(x, strand, npeaks = 1000, bandwidth = 150, mc.cores=1)
}
\arguments{
\item{x}{A \code{RangedDataList}, \code{RangedData} or an \code{IRangesList} object
  containing the aligned reads in each chromosome.}
\item{strand}{Strand that each read was aligned to.
  If \code{x} is of class \code{RangedDataList}, \code{strand}
  can be a character vector of length 1 indicating the name
  of the field in \code{x} indicating the strand,
  i.e. \code{x[[1]][[strand]]} contains the strand information.}
\item{npeaks}{Number of peaks to be used to estimate the shift size.}
\item{bandwidth}{Only reads with distance less than \code{bandwidth}
  between them and their closest gene are used to estimate the shift
  size.}
\item{mc.cores}{Number of cores to be used for parallel computing
  (passed on to \code{mclapply}). Only used if \code{x} is of class \code{RangedDataList}.}
}
\value{
  A \code{CompressedIRangesList} object with all reads shifted so that
  the strand specific bias is no longer present.
}
\details{
  The procedure detects the \code{npeaks} highest peaks (using reads
  from both strands simultaneously).
  Then it selects reads which are less than \code{bandwidth} base pairs away from
  any of the peaks.
  Then it computes (a) the average distance between reads on the plus
  strand and the closest peak, (b) the same distance for reads on the
  minus strand. The mean difference between (a) and (b) is the estimated
  shift size.
  Reads on the plus strand are shifted to the right, whereas reads on
  the minus strands are shifted to the left.
}
\examples{
#Generate 1000 reads containing strand-specific bias
st <- runif(1000,1,250)
strand <- rep(c('+','-'),each=500)
st[strand=='-'] <- st[strand=='-'] + runif(500,50,100)
x <- RangedData(IRanges(st,st+38),strand=strand)
#Estimate and remove the bias
xalign <- alignPeaks(x, strand='strand', npeaks=1)
}
\keyword{manip}
