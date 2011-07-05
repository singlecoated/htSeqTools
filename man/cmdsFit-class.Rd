\name{cmdsFit-class}
\Rdversion{1.1}
\docType{class}
\alias{cmdsFit-class}

\title{Class "cmdsFit" }
\description{ Classical Multi-Dimensional Scaling Fit. Function
  \code{cmds} creates object of this class.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("cmdsFit", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{points}:}{Object of class \code{"matrix"} with (x,y)
      coordinates in the approximated space. }
    \item{\code{d}:}{Object of class \code{"matrix"} with original
      distances between individuals. }
    \item{\code{dapprox}:}{Object of class \code{"matrix"} with
      distances between individuals in the approximated space. }
    \item{\code{R.square}:}{Percentage of variability in \code{d}
      explained by \code{dapprox} (object of class \code{"numeric"}) }
  }
}
\section{Methods}{
There are \code{show} and \code{plot} methods defined for this class.
}
\author{ David Rossell }
\seealso{
	\code{cmdscale} from package base.
}
\examples{
showClass("cmdsFit")
}
\keyword{classes}
