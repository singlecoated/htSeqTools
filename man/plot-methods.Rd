\name{plot-methods}
\docType{methods}
\alias{plot,cmdsFit,ANY-method}
\alias{plot,gridCover,ANY-method}
\title{ Methods for Function plot in Package `htSeq' }
\description{
 Methods for function \code{plot} in Package `htSeq'
}
\section{Methods}{
\describe{

\item{\code{signature(x = "cmdsFit", y = "ANY")}}{
Produces a Multi-Dimensional scaling plot. See \code{cmds} for details.
}

\item{\code{signature(x = "gridCover", y = "ANY")}}{
Plots the average coverage for each point in the grid. See
\code{gridCover} for details.
}
}}
\examples{
### Not run
#d <- matrix(c(0,5,10,5,0,15,10,15,0),byrow=TRUE,ncol=3)
#rownames(d) <- colnames(d) <- letters[1:3]
#fit1 <- cmdsFit(d,add=TRUE)
#plot(fit1)
}
\keyword{methods}

