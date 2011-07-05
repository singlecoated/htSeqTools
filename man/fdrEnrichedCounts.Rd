\name{fdrEnrichedCounts}
\alias{fdrEnrichedCounts}
\title{Posterior probability that a certain number of repeats are higher than expected by chance.}
\description{
Given a vector of number of repeats
(e.g. there are 100 sequences appearing once, 50 sequences appearing twice etc.)
the function computes the false discovery rate that each number of repeats is unusually high.
}
\usage{
fdrEnrichedCounts(counts,use=1:10,components=0,mc.cores=1)
}
\arguments{
\item{counts}{vector with observed frequencies. The vector must have
  names. \code{tabDuplReads} function can be used for this purpose.}
\item{use}{number of repeats to be used when estimating the null
  distribution. The number of repeats expected if no unusually high 
  repeats are present. The first 10 are used by default.}
\item{components}{number of negative binomials that will be used to fit
  the null distribution. The default value is 1. This value has to be
  between 0 and 4. If 0 is given the optimal number of negative biomials
  is chosen using the Bayesian information criterion (BIC)}
\item{mc.cores}{number of cores to be used to compute calculations. This
  parameter will be passed bt to \code{mclappply}}
}
\value{
  \code{data.frame} with the following columns:
  \item{pdfH0}{vector with pdf under the null hypothesis of no enrichment}
  \item{pdfOverall}{vector with pdf for mixture distribution}
  \item{fdrEnriched}{vector with false discovery rate that each count is significantly enriched}
}
\details{
The null distribution is a combination of n negative binomials where. n
is assigned  through the  \code{components} parameter.
If \code{components} is equal to 0 the optimal number of negative
binomials is choosen using the Bayesian information criterion (BIC).
The parameters of the null distribution are estimated from the number of
observations with as many repeats as told in the \code{use} parameter.
If use is 1:10 the null distribution will be estimated using repeats that
appear 1 time, 2 times, ... or 10 times.

False discovery rate for usually high number of repeats is done
following an empirical Bayes scheme similar to that in Efron et al. 
Let f0(x) be the null distribution, f(x) be the overall distribution and (1-pi0)
the proportion of unusually high repeats.
We assume the two component mixture f(x)= pi0 f0(x) + (1-pi0)f1(x).
Essentially, f(x) is estimated from the data
(imposing that f(x) must be monotone decreasing after its mode using \code{isoreg}
from packabe \code{base}, 
to improve the estimate in the tails).
pi0 is set to 1, i.e. its maximum possible value.
The false discovery rate for enrichment is computed as \code{1-pi0*(1-cumsum(f0(x)))/(1-cumsum(f(x)))}.
}
\references{
Ji et al. An integrated software system for analyzing ChIP-chip and ChIP-seq data. Nature Biotechnology, 2008, 26, 1293-1300.

Efron et al. Empirical Bayes Analysis of a Microarray Experiment, Journal of the American Statistical Association, 2001, 96, 1151-1160.
}
\examples{
#Generate 1000 sequences repeated once, on the average
nrepeats <- c(rpois(10^4,1),rpois(10,10))
nrepeats <- nrepeats[nrepeats>0]
counts <- table(nrepeats)
barplot(counts) -> xaxis #observe bimodality around 10
fdrest <- fdrEnrichedCounts(counts,use=1:5,components=1)
cutoff <- xaxis[which(fdrest$fdrEnriched<0.95)[1]]
abline(v=cutoff,col=2)
text(cutoff,counts[1]/2,'cut-off',col=2)
head(fdrest)
}
\keyword{ univar }
