\name{listOverlap}
\alias{listOverlap}
\alias{listOverlap-methods}
\alias{listOverlap,character,character,character,character-method}
\alias{listOverlap,character,character,missing,character-method}
\alias{listOverlap,factor,character,missing,missing-method}
\title{Assess the overlap between two or three lists.}
\description{
Assess the overlap between two or three lists, e.g. ChIP-Seq peaks
vs. genes selected from a microarray,
or peaks obtained in different experiments.
}
\section{Methods}{
\describe{

\item{\code{signature(list1 = "character", list2 = "character", list3 = "character", univ = "character")}}{ Studies 3-way associations.  }

\item{\code{signature(list1 = "character", list2 = "character", list3 = "missing", univ = "character")}}{ Studies bivariate associations. }

\item{\code{signature(list1 = "factor", list2 = "character", list3 = "missing", univ = "missing")}}{ Studies bivariate associations. }
}}
\usage{
listOverlap(list1, list2, list3, univ, ...)
}
\arguments{
\item{list1}{Vector with elements in the first list. This can either be
  a character vector indicating the element names, or a named factor
  vector indicating some classification for the elements in the first list.}
\item{list2}{Vector with elements in the second list. This should be a
  character vector indicating the element names.}
\item{list3}{Vector with elements in the third list. This should be a
  character vector indicating the element names. The overlap assesment method used
depends on whether this argument is specified or not. See details.}
\item{univ}{character vector indicating the universe of all elements
  from which \code{list1} and \code{list2} were obtained. The overlap
  assessment depends on whether this argument is specified or not. See details.}
\item{...}{Further arguments to be passed on to \code{chisq.test} in 2
  list overlapping.}
}
\value{
  For comparison of 2 lists, an \code{htest} object from a chi-square test that evaluates if the
  two lists are statistically independent from each other.
  This is a named list: the observed overlap is stored in
  \code{observed} and the P-value in \code{p.value}.
  
  For 3 list comparison, a \code{list} object containing the occurrence
  and frequency tables (\code{xtab}, \code{ftable}), the fitted linear models
  (\code{glm1}, \code{glm2}), and the anova P-value (\code{pvalue}).
}
\details{
  For signature(list1='character', list2='character', list3='missing', univ='character')
  the overlap is assessed with respect to the universe of all possible elements \code{univ}.
  That is, we count the number of elements that are common to
  list1 and list2, those appearing only in either list1 or list2, and
  those not appearing in either (but appearing in \code{univ}).
  A typical example: \code{list1} contains names of genes with a peak in ChIP-Seq
  experiment 1, \code{list2} names of genes with a peak in ChIP-Seq
  experiment 2, and \code{univ} the names of all genes in the organism.
  
  For signature(list1='character', list2='character', list='character', univ='character')
  the overlap is assessed by fitting and anova comparison
  of linear models. This is done to test whether 3-way overlap is significant with
  respect to the universe of all possible elements \code{univ} when
  compared to a model considering just the combination of 2-way overlapping.
  A typical example: \code{list1}, \code{list2} and \code{list3} contain
  names of genes with peaks in three different ChIP-Seq experiments, and \code{univ} the names of all genes in the organism.

  For signature(list1='factor', list2='character', univ='missing') the
  distribution of \code{list1} is compared between elements appearing
  and not appearing in \code{list2}.
  A typical example: \code{list1} indicates the differential expression
  status for a number of genes, and \code{list2} contains the names of
  the genes which had a peak in a ChIP-Seq experiment.
}
\examples{
#Overlap between diff expression and chip-seq peaks
deStatus <- factor(c(0,0,0,0,1,1,1))
names(deStatus) <- paste('Gene',1:7)
peaks <- c('Gene 6','Gene 7')
ans <- listOverlap(list1=deStatus,list2=peaks)
ans$observed
ans$p.value

#Overlap between peaks obtained from two different experiments
peaks2 <- c('Gene 1','Gene 2','Gene 7')
univ <- paste('Gene',1:7)
ans <- listOverlap(list1=peaks,list2=peaks2,univ=univ)
ans$observed
ans$p.value
}
\keyword{datasets}
