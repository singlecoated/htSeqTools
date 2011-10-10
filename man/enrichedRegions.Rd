\name{enrichedRegions}
\alias{enrichedRegions}
\alias{enrichedRegions-methods}
\alias{enrichedRegions,RangedData,RangedData,missing,ANY,ANY-method}
\alias{enrichedRegions,RangedData,missing,missing,ANY,ANY-method}
\alias{enrichedRegions,RangedDataList,missing,missing,ANY,ANY-method}
\alias{enrichedRegions,missing,missing,RangedData,ANY,missing-method}
\alias{enrichedRegions,missing,missing,RangedData,ANY,numeric-method}
\title{Find significantly enriched regions in sequencing experiments.}
\description{
  Find regions with a significant accumulation of reads in a sequencing experiment.
}
\section{Methods}{
  \describe{
    
    \item{\code{signature(sample1 = "missing", sample2 =
	"missing", regions = "RangedData")}}{\code{ranges(regions)}
      indicates the chromosome, start and end of genomic regions, while \code{values{regions}} should
      indicate the observed number of reads for each group in each
      region. \code{enrichedRegions} tests the null hypothesis that the
      proportion of reads in the region is equal across all groups via a
      likelihood-ratio test (or permutation-based chi-square for regions
      where the expected counts are below 5 for some group). }

    \item{\code{signature(sample1 = "RangedDataList", sample2 =
	"missing", regions = "missing")}}{Each element in \code{sample1} contains the read
      start/end of an individual sample. \code{enrichedRegions} identifies
      regions with high concentration of reads (across all samples) and
      then compares the counts across groups using a likelihood-ratio test
      (or permutation-based chi-square for regions
      where the expected counts are below 5 for some group).}
    
    \item{\code{signature(sample1 = "RangedData", sample2 = "RangedData",
	regions = "missing")}}{ \code{space(sample1)}
      indicates the chromosome, \code{start(sample1)} and
      \code{end(sample1)} the start/end position of the reads. Similarly for
      \code{sample2}. \code{enrichedRegions} identifies regions with high
      concentration of reads (across all samples) and then compares the
      counts across groups using a likelihood-ratio test (or
      permutation-based chi-square for regions where the expected counts are
      below 5 for some group).}

    \item{\code{signature(sample1 = "RangedData", sample2 = "missing",
	regions = "missing")}}{\code{space(sample1)} indicates the chromosome, \code{start(sample1)}
      and \code{end(sample1)} the start/end position of the reads.
      \code{enrichedRegions} tests the null hypothesis that an unusually high proportion of reads has been
      observed in the region using an exact binomial test.}

}}
\usage{
enrichedRegions(sample1, sample2, regions, minReads=10, mappedreads,
pvalFilter=0.05, exact=FALSE, p.adjust.method='none', twoTailed=FALSE,
mc.cores=1) 
}
\arguments{
  \item{sample1}{Either start and end of sequences in sample 1 (\code{IRangesList}, \code{RangedData} or \code{IRanges} object),
    of \code{RangedDataList} with sequences for all samples
    (\code{sample2} must be left missing in this case) .}
  \item{sample2}{Same for sample 2. Can be left missing.}
  \item{regions}{If specified, the analysis is restricted to the regions
    indicated in \code{regions}. If not specified, the regions are
    automatically defined using the argument \code{minReads}.}
  \item{minReads}{This argument is only used when \code{regions} is not
    specified. The regions to be tested for enrichment are those with coverage
    greater or equal than \code{minReads}.
    If \code{sample1} is a \code{RangedDataList}, the overall coverage
    adding all samples is used.
    Otherwise, if \code{twoTailed} is FALSE, only the reads in
    sample 1 are counted. If \code{twoTailed} is TRUE, the sum of reads
    in samples 1 and 2 are counted.}
  \item{mappedreads}{Number of mapped reads for the sample. Has to be of
    class integer. Will be used to compute RPKM.}
  \item{pvalFilter}{Only regions with P-value below \code{pvalFilter}
    are reported as being enriched.}
  \item{exact}{ If set to TRUE, an exact test is used whenever
    some expected cell counts are 5 or less
    (chi-square test based on permutations if \code{sample1} is a
    \code{RangedDataList} object, Fisher's exact test otherwise),
    i.e. when the asymptotic chi-square/likelihood-ratio test calculations break down. Ignored
    if \code{sample2} is missing, as in this case calculations are
    always exact.}
  \item{p.adjust.method}{P-value adjustment method, passed on to
    \code{p.adjust}.}
  \item{twoTailed}{If set to FALSE, only regions with a higher
    concentration of reads in sample 1 than in sample 2 are reported. If
    set to TRUE, regions with higher concentration of sample 2 reads are
    also reported. Ignored if \code{sample2} is missing.}
  \item{mc.cores}{If \code{mc.cores} is greater than 1, computations are
    performed in parallel for each element in the
    \code{IRangesList} objects. Whenever possible the \code{mclapply}
    function is used, therefore exactly \code{mc.cores} are used. For
    some signatures \code{mclapply} cannot be used, in which case the \code{parallel} function from package
    \code{multicore} is used. Note: the latter option launches as many parallel
    processes as there are elements in \code{x}, which can place strong
    demands on the processor and memory.}
}
\value{
  Object of class \code{RangedData}
  indicating the significantly enriched regions, the number of reads in
  each sample for those regions, the fold changes (adjusted considering
  the overall number of sequences in each sample) and the chi-square
  test P-values.
}
\details{
  The calculations depend on whether \code{sample2} is missing or not.
  Non-missing \code{sample2} case.
  First, regions with coverage above \code{minReads} are
  selected. Second, the number of reads falling in the selected regions
  are computed for sample 1 and sample 2.
  Third, the counts are compared via a chi-square test (with Yates
  continuity correction), which takes into
  account the total number of sequences in each sample.
  Finally, statistically significant regions are selected and returned
  in \code{RangedData} or \code{RangedDataList} objects.

  Missing \code{sample2}. First, regions with coverage above \code{minReads} are selected.
  Second, the number of reads in sample 1 falling in the selected regions is computed.
  Third, the proportion of reads in each region is tested for enrichment via a one-tailed Binomial exact test.
}
\examples{
set.seed(1)
st <- round(rnorm(1000,500,100))
strand <- rep(c('+','-'),each=500)
space <- rep('chr1',length(st))
sample1 <- RangedData(IRanges(st,st+38),strand=strand,space=space)
st <- round(rnorm(1000,1000,100))
sample2 <- RangedData(IRanges(st,st+38),strand=strand,space=space)
enrichedRegions(sample1,sample2,twoTailed=TRUE)
}
\keyword{datasets}
