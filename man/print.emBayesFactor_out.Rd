% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emBayesFactor_out.R
\name{print.emBayesFactor_out}
\alias{print.emBayesFactor_out}
\title{Print Methods for Effect-Size & Moment Bayes Factors}
\usage{
\method{print}{emBayesFactor_out}(x, ...)
}
\arguments{
\item{x}{object of class \code{"emBayesFactor_out"}.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
the argument \code{x}, invisibly,
as for all \code{\link{print}} methods.
}
\description{
Printing objects of class \code{"emBayesFactor_out"} by simple print method.
}
\details{
An \code{emBayesFactor_out} object is just a named list of numbers and character strings,
supplemented with \code{method} and \code{note} elements.
The \code{method} is displayed as a title, the \code{note} as a footnote.
}
\examples{
## Generate object 'out' of class 'emBayesFactor_out'
out <- ettest_tStat(
  tStat = 2.03,
  N1 = 80,
  Cohen_d = 0.2
)

## call method 'print.emBayesFactor_out' on 'out'
print(out)
}
\seealso{
\link[stats]{print.power.htest}
}
\author{
Constantin G. Meyer-Grant
(\email{constantin.meyer-grant@psychologie.uni-freiburg.de})
}
\keyword{internal}
