\name{filterDuplReads}
\Rdversion{1.1}
\alias{filterDuplReads}
\alias{tabDuplReads}
\alias{filterDuplReads-methods}
\alias{filterDuplReads,RangedData-method}
\alias{filterDuplReads,RangedDataList-method}
\alias{tabDuplReads-methods}
\alias{tabDuplReads,RangedData-method}
\alias{tabDuplReads,RangedDataList-method}
\title{ Detect and filter duplicated reads/sequences. }
\description{
\code{filterDuplReads} filters highly repeated sequences, i.e. with the same chromosome, start and
end positions.
As many such sequences are likely due to over-amplification artifacts, this
can be a useful pre-processing step for ultra high-throughput sequencing
data.
A false discovery rate is computed for each number of repeats being
unusually high. The reads with a higher false discovery rate will be
removed. For more information on the false discovery rate calculation
please read the \code{fdrEnrichment} manual.

\code{tabDuplReads} counts the number reads with no duplications, duplicated once, twice etc.
}
\section{Methods}{
Methods for \code{filterDuplReads} and \code{tabDuplReads}

\describe{
\item{\code{signature(x = "RangedData")}}{ Two reads are duplicated if
  they have the same space, start and end position. }

\item{\code{signature(x = "RangedDataList")}}{ The method is applied
  separately to each \code{RangedData} element in the list. }
}}
\usage{
filterDuplReads(x, maxRepeats, fdrOverAmp=0.01, negBinomUse=.999,components=0, mc.cores=1)

tabDuplReads(x,  minRepeats=1, mc.cores=1)
}
\arguments{
  \item{x}{ Object containing read locations.
    Currently methods for \code{RangedData} and \code{RangedDataList}.
    Duplication is assessed based only on the space,
    start, end and \code{x[['strand']]}, i.e. even if they are different
    based on other variables stored in \code{values(x)}, the reads are
    considered duplicated and only the first appearance is returned.}
  \item{maxRepeats}{Reads appearing \code{maxRepeats} or more times will
    be excluded. If not specified, this is setup automatically based on
    \code{fdrOverAmp}.} 
  \item{fdrOverAmp}{Reads with false discovery rate of being
    over-amplified greater than \code{fdrOverAmp} are excluded.}
  \item{negBinomUse}{Number of counts that will be used to compute the
    null distribution. Using 1 - 1/1000 would mean that 99.9\% of the
    reads will be used. The ones with higher number of repetitions are
    the excluded ones.}
  \item{components}{number of negative binomials that will be used to fit
    null distribution. The default value is 1. This value hase to be
    between 0 and 4. If 0 is given the optimal number of negative biomials
    is choosen using the Bayesian information criterion (BIC)}
  \item{mc.cores}{Number of cores to be used in parallel computing
    (passed on to \code{mclapply}).}
  \item{minRepeats}{The table is only produced for reads with at least
    \code{minRepeats} repeats.}
}
\value{ \code{filterDuplReads} returns \code{x} without highly
  repetitive sequencesas, determined by
  \code{maxRepeats} or \code{ppOverAmp}.
  
\code{tabDuplReads} returns a table counting the number of sequences
repeating 1 times, 2 times, 3 times etc.  
}
\seealso{ \code{fdrEnrichedCounts} to compute the posterior probability
  that a certain number of repeats is due to over-amplification. }
\author{ Evarist Planet, David Rossell, Oscar Flores }
\examples{
set.seed(1)
st <- round(rnorm(1000,500,100))
strand <- rep(c('+','-'),each=500)
space <- sample(c('chr1','chr2'),size=length(st),replace=TRUE)
sample1 <- RangedData(IRanges(st,st+38),strand=strand,space=space)

#Add artificial repeats
st <- rep(400,20)
repeats <- RangedData(IRanges(st,st+38),strand='+',space='chr1')
sample1 <- rbind(sample1,repeats)

filterDuplReads(sample1)
}
\keyword{ manip }
