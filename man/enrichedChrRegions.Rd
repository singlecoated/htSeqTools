\name{enrichedChrRegions}
\alias{enrichedChrRegions}
\alias{enrichedChrRegions-methods}
\alias{enrichedChrRegions,RangedData,missing-method}
\alias{enrichedChrRegions,RangedData,RangedData-method}
\title{
Find chromosomal regions with a high concentration of hits.
}
\description{
This function looks for chromosomal regions where there is a large
accumulation of hits, e.g. significant peaks in a chip-seq experiment or
differentially expressed genes in an rna-seq or microarray experiment.
Regions are found by computing number of hits in a moving window
and selecting regions based on a FDR cutoff.
}
\section{Methods}{
\describe{

\item{\code{signature(hits1 = "RangedData", hits2 = "missing")}}{
Look for chromosome zones with a large number of hits reported in \code{hits1}.
}

\item{\code{signature(hits1 = "RangedData", hits2 = "RangedData")}}{
Look for chromosomal zones with a different density of hits in
\code{hits1} vs \code{hits2}.
}

}}
\usage{
enrichedChrRegions(hits1, hits2, chrLength, windowSize=10^4-1, fdr=0.05, nSims=10, mc.cores=1)
}
\arguments{
  \item{hits1}{Object containing hits (start, end and chromosome). Currently only \code{RangedData}
    objects are accepted.}
  \item{hits2}{Optionally, another object containing hits. If specified,
  regions will be defined by comparing hits1 vs hits2.}
  \item{chrLength}{Named vector indicating the length of each chromosome
    in base pairs}
  \item{windowSize}{Size of the window used to smooth the hit count (see
    details)}
  \item{fdr}{Desired FDR level (see details)}
  \item{nSims}{Number of simulations to be used to estimate the FDR}
  \item{mc.cores}{Number of processors to be used in parallel
    computations (passed on to mclapply)}
}
\value{
  Object of class \code{RangedData} containing regions with smoothed hit
  count above the specified FDR level.
}
\details{
A smoothed number of hits is computed by counting the number of
hits in a moving window of size \code{windowSize}.
Notice that only the mid-point of each hit in \code{hits1} (and
\code{hits2} if specified) is used. That is,
hits are not treated as intervals but as being located at a single base
pair.

If \code{hits2} is missing, regions with large smoothed number of hits
are selected.
To assess statistical significance, we generate hits (also 1
base pair long) randomly distributed along the genome and
compute the smoothed number of hits.
The number of simulated hits
is set equal to \code{nrow(hits1)}.
The process is repeated
\code{nSims} times, resulting in several independent simulations.
To estimate the FDR, several thresholds to define enriched chromosomal
regions are considered. For each threshold, we count the 
number of regions above the threshold in the observed data and in the
simulations. For each threshold t, the FDR is estimated as
the average number of regions with score >=t in the simulations
over the number of regions with score >=t in the observed data.

If \code{hits2} is not missing, the difference in smoothed proportion of
hits (i.e. the number of hits in the window divided by the overall number of hits)
between the two groups is used as a test statistic.
To assess statistical significance, we generate randomly scramble hits
between sample 1 and sample 2 (maintaining the original number of hits
in each sample), and we re-compute the test statistic.
The FDR for a given threshold t is estimated as the
number of bases in the simulated data with test statistic>t divided by
number of bases in observed data with test statistic>t.

The lowest t with estimated FDR below \code{fdr} is used to define
enriched chromosomal regions.
}
\examples{
set.seed(1)
st <- round(rnorm(100,500,100))
st[st>10000] <- 10000
strand <- rep(c('+','-'),each=50)
space <- rep('chr1',length(st))
hits1 <- RangedData(IRanges(st,st+38),strand=strand,space=space)
chrLength <- c(chr1=10000)
enrichedChrRegions(hits1,chrLength=chrLength, windowSize=99, nSims=1)
}
\keyword{stats}
