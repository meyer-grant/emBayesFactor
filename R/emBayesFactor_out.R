#' @title Print Methods for Effect-Size & Moment Bayes Factors
#'
#' @description
#' Printing objects of class \code{"emBayesFactor_out"} by simple print method.
#'
#'
#' @param x object of class \code{"emBayesFactor_out"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' An \code{emBayesFactor_out} object is just a named list of numbers and character strings,
#' supplemented with \code{method} and \code{note} elements.
#' The \code{method} is displayed as a title, the \code{note} as a footnote.
#'
#'
#' @return the argument \code{x}, invisibly,
#' as for all \code{\link{print}} methods.
#'
#' @author Constantin G. Meyer-Grant
#' (\email{constantin.meyer-grant@psychologie.uni-freiburg.de})
#'
#' @seealso \link[stats]{print.power.htest}
#'
#' @examples
#' ## Generate object 'out' of class 'emBayesFactor_out'
#' out <- ettest_tStat(
#'   tStat = 2.03,
#'   N1 = 80,
#'   Cohen_d = 0.2
#' )
#'
#' ## call method 'print.emBayesFactor_out' on 'out'
#' print(out)



#' @keywords internal
#' @method print emBayesFactor_out
#' @export
print.emBayesFactor_out <- function(x, ...) {
  cat("\n")
  cat(strwrap(x$method, prefix = "\t"), sep = "\n")
  cat("\n")
  cat(paste("BayesFactor:  ", format(x$BayesFactor,format="e"), sep=""),
      "with absolute error <", formatC(x$accuracy,format="e",digits=4), "\n",
      sep=" ")
  cat("  [",x$note,"]",sep="")
  cat("\n--------------\n")
  cat("prior specifications: \n")
  if(is.null(x$rscale)){
    cat(paste(" ",x$Cohen,"=",format(x$Cohen_size,format="e",digits=4),
              sep=" "),
        paste("nu =", format(x$nu,format="e",digits=4), sep=" "), sep=", ")
  }else{
    cat(paste(" ",x$Cohen,"=",format(x$Cohen_size,format="e",digits=4),
              sep=" "),
        paste("rscale =", format(x$rscale,format="e",digits=4), sep=" "),
        paste("nu =", format(x$nu,format="e",digits=4), sep=" "), sep=", ")
  }
  cat("\n--------------\n")
  cat("observed ", x$Stat, ":  ", format(x$Stat_size,format="e"), sep="")
  cat("\n")
  if(!is.null(x$alternative)){
    alt_char <-
      switch(x$alternative,
             two.sided = "not equal to",
             less = "less than",
             greater = "greater than")
    cat("alternative hypothesis:  \n  \"true effect", alt_char, "0\"")
  }
  cat("\n")
}
