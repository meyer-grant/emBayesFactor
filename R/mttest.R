#' @name mttest
#' @rdname mttest
#'
#' @title Moment Bayes Factors for \eqn{t} Tests
#'
#' @description
#' Computes Bayes factors (alternative over null hypothesis)
#' with moment priors for one- and two-sample \eqn{t} tests
#' based on the observed test statistic \eqn{t}.
#'
#'
#' @param tStat observed value of \eqn{t} statistic.
#' @param N1 sample size of the first group
#' (or only group for one-sample tests).
#' @param N2 sample size of the second group
#' (only specified for two-sample tests).
#' @param Cohen_d focused-on effect size (Cohen's \eqn{d}).
#' Must be either numeric, or alternatively set to one of the options
#' \code{"medium" = .50} (default), \code{"small" = .20}, \code{"large" = .80}.
#' @param nu prior degrees of freedom.
#' If not specified, it is set internally to the recommended value (\code{nu = 5}).
#' @param alternative a character string specifying the alternative hypothesis.
#' Must be one of \code{"two.sided"} (default),
#' \code{"greater"}, or \code{"less"}.
#'
#' @details
#' The Bayes factor provided by \code{mttest}
#' compares the alternative hypothesis with the null hypothesis
#' that there is no true effect.
#'
#' Note that \code{alternative = "greater"} denotes the alternative that
#' the true effect is positive, whereas \code{alternative = "less"}
#' denotes the alternative that the true effect is negative.
#'
#' \code{Cohen_d} should only be negative if \code{alternative = "less"} and must never be 0.
#'
#' \code{nu} must be at least 3.
#'
#' @return An object of class \code{mBf} containing the following components:
#'
#' \item{method}{the focussed-on analysis.}
#' \item{BayesFactor}{value of the Bayes factor.}
#' \item{accuracy}{estimate of the modulus of the absolute error.}
#' \item{Cohen}{the focussed-on effect size.}
#' \item{Cohen_size}{value of the focussed-on effect size.}
#' \item{nu}{prior degrees of freedom.}
#' \item{Stat}{provided statistic.}
#' \item{Stat_size}{value of the provided statistic.}
#' \item{alternative}{the specified alternative hypothesis}
#'
#' @references
#' Klauer, K. C., Meyer-Grant, C. G., & Kellen, D (2023).
#' \emph{On Bayes Factors for Hypotheses Tests}.
#' Manuscript submitted for publication.
#'
#' Pramanik, S., & Johnson, V. E. (2022).
#' Efficient alternatives for bayesian hypothesis tests in psychology.
#' Advance online publication. \emph{Psychological Methods}.
#'
#' @author Constantin G. Meyer-Grant
#' (\email{constantin.meyer-grant@psychologie.uni-freiburg.de})
#'
#' @seealso [integrate], [t.test], \link[BayesFactor]{ttestBF}
#'
#' @export
#'
#' @examples
#' ## Example from Klauer et al. (2023)
#'
#' ## Compute moment Bayes factor based on 'tStat'
#' mttest_tStat(
#'   tStat = 2.03,
#'   N1 = 80,
#'   Cohen_d = 0.2
#' )


#' @importFrom stats integrate dt pt
mttest_tStat <- function(tStat, N1, N2 = NULL,
                         Cohen_d = c("medium", "small", "large"),
                         nu = 5,
                         alternative = c("two.sided", "less", "greater")
) {

  alternative <- match.arg(alternative)

  if (!is.numeric(tStat) || is.na(tStat)) {
    stop("'tStat' must be numeric")
  }

  if (length(Cohen_d)>1 && !missing(Cohen_d)){
    warning("'Cohen_d' only takes one value, the first one was used")
    Cohen_d <- Cohen_d[1]
  }
  if(missing(Cohen_d) || is.null(Cohen_d)){
    d <- 0.5
  } else {
    if (is.numeric(Cohen_d) && !is.na(Cohen_d)){
      d <- Cohen_d
    } else if (is.character(Cohen_d) && !is.na(Cohen_d)) {
      d <- switch(Cohen_d,
                medium=0.5,
                small=0.2,
                large=0.8,
                stop(paste("'Cohen_d' must be either numeric",
                           "or alternatively set to one of the options",
                           "\"medium\", \"small\", or \"large\""))
      )
    } else {
      stop(paste("'Cohen_d' must be either numeric",
                 "or alternatively set to one of the options",
                 "\"medium\", \"small\", or \"large\""))
    }
  }

  if (d==0) {
    stop("'Cohen_d' must not be 0")
  }

  if (!is.numeric(nu) || is.na(nu)) {
    stop("'nu' must be numeric")
  } else if (nu < 3) {
    stop("'nu' must at least be 3")
  }

  is.whole <- function(x,
                       tol = .Machine$double.eps^0.5) {
    if(is.numeric(x)){
      return(abs(x - round(x)) < tol)
    }else{
      return(FALSE)
    }
  }

  if (!is.numeric(N1) || is.na(N1)) {
    stop("'N1' must be an integer")
  } else {
    if (!is.whole(N1)) {
      suppressWarnings(N1 <- as.integer(N1))
      if(is.na(N1)){stop("'N1' must be an integer")}
      warning(paste("provided value for 'N1' not an integer, 'N1' =",
                    as.character(N1),
                    "was used instead"))
    }
  }
  if (!is.numeric(N2) || is.na(N2)) {
    if(is.null(N2)) {
      # one-sample t test
      method <- paste("Moment Bayes Factor:",
                      "t Test (One-Sample)",
                      sep=" ")
      N2 <- 0
      nu_t <- N1 - 1
      N <- N1
    } else {
      stop("'N2' must be an integer")
    }
  } else{
    if (N2 == 0) {
      # one-sample t test
      method <- paste("Moment Bayes Factor:",
                      "t Test (One-Sample)",
                      sep=" ")
      nu_t <- N1 - 1
      N <- N1
    } else if (!is.whole(N2)){
      # two-sample t test
      method <- paste("Moment Bayes Factor:",
                      "t Test (Two-Sample)",
                      sep=" ")
      suppressWarnings(N2 <- as.integer(N2))
      if(is.na(N2)){stop("'N2' must be an integer")}
      nu_t <- N1 + N2 - 2
      N = N1 * N2 / (N1 + N2)
      warning(paste("provided value for 'N2' not an integer, 'N2' =",
                    as.character(N2),
                    "was used instead"))
    } else {
      # two-sample t test
      method <- paste("Moment Bayes Factor:",
                      "t Test (Two-Sample)",
                      sep=" ")
      nu_t <- N1 + N2 - 2
      N <- N1 * N2 / (N1 + N2)
    }
  }

  if (alternative == "two.sided") {

    if(d < 0){
      warning(paste("if 'alternative' is \"two.sided\",",
                    "the sign of 'Cohen_d' is ignored"))
      d = abs(d)
    }

    if(!is.infinite(nu)){
      moment_dt <- Vectorize(function(d_, d, nu) {

        tausq <- d^2 * (nu-1) / 2 / nu

        out <- d_^2 / tausq^1.5 * dt(d_ / sqrt(tausq), nu) * (nu-2) / nu

        return(out)

      }, "d_")
    } else {
      moment_dt <- Vectorize(function(d_, d, nu) {

        tausq <- d^2 / 2

        out <- d_^2 / tausq^1.5 * dt(d_ / sqrt(tausq), nu)

        return(out)

      }, "d_")
    }

    pseudoBayest <- Vectorize(function(d_, N, tStat, nu_t, d, nu) {

      ncp <- d_ * sqrt(N)

      out <- (dt(tStat, nu_t, ncp) / dt(tStat, nu_t, 0)) *
        moment_dt(d_, d, nu)

      #out <- ifelse(is.infinite(out) | is.na(out) | is.nan(out), 0, out)

      return(out)

    }, "d_")

    suppressWarnings({
      z <- integrate(f=pseudoBayest,
                     N=N, tStat=tStat, nu_t=nu_t, d=d, nu=nu,
                     lower=-Inf, upper=Inf)
    })

  } else if (alternative == "less") {

    if(d > 0){
      warning("if 'alternative' is \"less\", 'Cohen_d' should not be positive")
    }

    if(!is.infinite(nu)){
      moment_dt <- Vectorize(function(d_, d, nu) {

        tausq <- d^2 * (nu-1) / 2 / nu

        out <- d_^2 / tausq^1.5 * dt(d_ / sqrt(tausq), nu) * (nu-2) / nu

        return(out)

      }, "d_")
    } else {
      moment_dt <- Vectorize(function(d_, d, nu) {

        tausq <- d^2 / 2

        out <- d_^2 / tausq^1.5 * dt(d_ / sqrt(tausq), nu)

        return(out)

      }, "d_")
    }

    pseudoBayest <- Vectorize(function(d_, N, tStat, nu_t, d, nu) {

      ncp <- d_ * sqrt(N)

      out <- (dt(tStat, nu_t, ncp) / dt(tStat, nu_t, 0)) *
        2 * moment_dt(d_, d, nu)

      #out <- ifelse(is.infinite(out) | is.na(out) | is.nan(out), 0, out)

      return(out)

    }, "d_")

    suppressWarnings({
      z <- integrate(f=pseudoBayest,
                     N=N, tStat=tStat, nu_t=nu_t, d=d, nu=nu,
                     lower=-Inf, upper=0)
    })

  } else if (alternative == "greater") {

    if(d < 0){
      warning("if 'alternative' is \"greater\", 'Cohen_d' should not be negative")
    }

    if(!is.infinite(nu)){
      moment_dt <- Vectorize(function(d_, d, nu) {

        tausq <- d^2 * (nu-1) / 2 / nu

        out <- d_^2 / tausq^1.5 * dt(d_ / sqrt(tausq), nu) * (nu-2) / nu

        return(out)

      }, "d_")
    } else {
      moment_dt <- Vectorize(function(d_, d, nu) {

        tausq <- d^2 / 2

        out <- d_^2 / tausq^1.5 * dt(d_ / sqrt(tausq), nu)

        return(out)

      }, "d_")
    }

    pseudoBayest <- Vectorize(function(d_, N, tStat, nu_t, d, nu) {

      ncp <- d_ * sqrt(N)

      out <- (dt(tStat, nu_t, ncp) / dt(tStat, nu_t, 0)) *
        2 * moment_dt(d_, d, nu)

      #out <- ifelse(is.infinite(out) | is.na(out) | is.nan(out), 0, out)

      return(out)

    }, "d_")

    suppressWarnings({
      z <- integrate(f=pseudoBayest,
                     N=N, tStat=tStat, nu_t=nu_t, d=d, nu=nu,
                     lower=0, upper=Inf)
    })

  } else {
    stop("'arg' should be one of \"two.sided\", \"less\", \"greater\"")
  }

  emBayesFactor_out <- list(method = method,
                  BayesFactor = z$value, accuracy = z$abs.error,
                  Cohen = "Cohen's d", Cohen_size = d,
                  rscale = NULL, nu = nu,
                  Stat = "t statistic", Stat_size = tStat,
                  alternative = alternative,
                  note = "alternative over null")

  class(emBayesFactor_out) <- "emBayesFactor_out"

  return(emBayesFactor_out)

}
