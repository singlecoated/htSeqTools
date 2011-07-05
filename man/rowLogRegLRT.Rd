\name{rowLogRegLRT}
\alias{rowLogRegLRT}
\title{
Row-wise logistic regression
}
\description{
  Row-wise logistic regressions are applied to a matrix with counts.
  For each row, an overall test comparing the column counts across
  columns is performed. Optionally, chi-square permutation tests are
  used when the expected counts are below 5 for some column.
}
\usage{ rowLogRegLRT(counts, exact = TRUE, p.adjust.method = "none") }
\arguments{
  \item{counts}{Matrix with counts}
  \item{exact}{ If set to TRUE, an exact test is used whenever
    some expected cell counts are 5 or less}
  \item{p.adjust.method}{p-value adjustment method, passed on to \code{p.adjust}}
}
\details{
For each column, the proportion of counts in each row (with respect to
the overall counts in that column) is computed. Then a statistical
comparison of these proportions across groups is performed via a
likelihood-ratio test (if \code{exact==TRUE} a permutation based
chi-square test is used whenever the expected counts in some column is
below 5).

Notice that data from column \code{j} can be viewed as a multinomial
distribution with probabilities pj, where pj is a vector of length
\code{nrow(x)}.
\code{rowLogRegLRT} tests the null hypothesis p1[i]=...pc[i] for
i=1...\code{nrow(x)},
where c is \code{ncol(x)}.
This actually ignores the multinomial sampling model and focuses on its
binomial margins, which is a reasonable approximation when the number
\code{nrow(x)} is large and substantially improves computation speed.
}
\examples{
#The first two rows present different counts across columns
#The last two columns do not
x <- matrix(c(70,10,10,10,35,35,10,10),ncol=2)
x
rowLogRegLRT(x)
}
\keyword{htest}
