\name{htSample}
\alias{htSample}
\docType{data}
\title{
Example ChIP-sequencing data with 2 replicates per group obtained in two
different dates.
}
\description{
  This \code{RangedDataList} contains a subset of drosophila
  melanogaster ChIP-sequencing data
  obtained with the Illumina sequencer.
  An immuno-precipitated and a control input sample were obtained at two
  experimental dates (details not provided as this is still unpublished data).
  In order to save space and let the examples run quicker,
  only reads mapping to the first 500kb of chr2L are included.
}
\usage{data(htSample)}
\format{
  \code{RangedDataList} where each element contains reads from a
  different sample.
  \code{names(htSample)} indicate the group and batch (experimental
  date) that each sample corresponds to.
}
\details{
Data was pre-processed using the Illumina pipeline and mapped to the
drosophila melanogaster dm3 genome using Bowtie.
Only uniquely mapping sequences with at most 2 mismatches in the first
28 bases were kept.
See the package vignette for some more details on this dataset.
}
\examples{
data(htSample)
htSample
}
\keyword{datasets}
