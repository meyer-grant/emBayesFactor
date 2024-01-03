#' @name mreg
#' @title Moment Bayes Factors for (Multiple) Linear Regression
#'
#' @description
#' Computes Bayes factors (alternative over null hypothesis)
#' with moment priors for (multiple) linear regression analyses
#' based on observed partial coefficient of determination (\eqn{R^2_p}),
#' the observed test statistic \eqn{F}, or the observed test statistic \eqn{t}
#' for testing a single partial regression coefficient.
#'
#'
#'
#' @param R2p observed partial coefficient of determination \eqn{R^2_p}.
#' @param FStat observed value of \eqn{F} statistic.
#' @param tStat observed value of \eqn{t} statistic.
#'
#' @param N number of observations.
#' @param k number of parameters of the full model.
#' @param p number of parameters of the reduced model.
#' @param Cohen_f2 focused-on effect size (Cohen's \eqn{f^2}).
#' Must be either numeric, or alternatively set to one of the options
#' \code{"medium" = .15} (default), \code{"small" = .02}, \code{"large" = .35}.
#' @param nu prior degrees of freedom.
#' If not specified, it is set internally to the recommended value (\code{nu = 5}).
#' @param alternative a character string specifying the alternative hypothesis.
#' Must be one of \code{"two.sided"} (default),
#' \code{"greater"}, or \code{"less"}.
#'
#' @details
#' The Bayes factor provided by \code{mreg} compares
#' the alternative hypothesis with the null hypothesis
#' that there is no true effect.
#' For \code{mreg_tStat}, \code{alternative = "greater"} denotes
#' the alternative that the true effect is positive,
#' whereas \code{alternative = "less"} denotes the alternative
#' that the true effect is negative.
#'
#' \code{Cohen_f2} must be positive, but if \code{alternative = "less"}
#' is specified for \code{mreg_tStat},
#' the presumed relationship between predictor and criterion is negative.
#'
#' \code{nu} must be at least 3.
#'
#' @return An object of class \code{emBayesFactor_out} containing the following components:
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
#' Klauer, K. C., Meyer-Grant, C. G., & Kellen, D (2024).
#' \emph{On Bayes Factors for Hypotheses Tests}.
#' Manuscript submitted for publication.
#'
#' @author Constantin G. Meyer-Grant
#' (\email{constantin.meyer-grant@psychologie.uni-freiburg.de})
#'
#' @seealso [integrate], [lm], \link[BayesFactor]{regressionBF}
#'
#' @examples
#' ## Example from Klauer et al. (2024)
#'
#' ## Compute partial R^2 from the two models' R^2
#' R2p <- (0.7109 - 0.5670) / (1.0 - 0.5670)
#'
#' ## Compute moment Bayes factor based on 'R2p'
#' mreg_R2p(
#'   R2p = R2p,
#'   N = 175,
#'   p = 4,
#'   k = 5,
#'   Cohen_f2 = "large"
#' )


#' @importFrom stats df integrate dt pt
#' @rdname mreg
#' @export
mreg_R2p <- function(
    R2p, N, k, p,
    Cohen_f2 = c("medium", "small", "large"),
    nu = 5 + (k - p)) {

  if (!is.numeric(R2p) || is.na(R2p)) {
    stop("'R2p' must be numeric")
  } else if (R2p < 0 || R2p > 1) {
    stop("'R2p' must be greater than zero and smaller than one")
  }

  if (length(Cohen_f2)>1 && !missing(Cohen_f2)){
    warning("'Cohen_f2' only takes one value, the first one was used")
    Cohen_f2 <- Cohen_f2[1]
  }
  if(missing(Cohen_f2) || is.null(Cohen_f2)){
    f <- sqrt(0.15)
  } else {
    if (is.numeric(Cohen_f2) && !is.na(Cohen_f2)){
      f2 <- as.numeric(Cohen_f2)
      if (f2 <= 0) {
        warning(paste("'Cohen_f2' must be positive,",
                      "'Cohen_f2' = \"medium\" used instead"))
        f <- sqrt(0.15)
      } else {
        f <- sqrt(f2)
      }
    } else if (is.character(Cohen_f2) && !is.na(Cohen_f2)) {
      f <- switch(Cohen_f2,
                medium=sqrt(0.15),
                small=sqrt(0.02),
                large=sqrt(0.35),
                stop(paste("'Cohen_f2' must be either numeric",
                           "or alternatively set to one of the options",
                           "\"medium\", \"small\", or \"large\""))
      )
    } else {
      stop(paste("'Cohen_f2' must be either numeric",
                 "or alternatively set to one of the options",
                 "\"medium\", \"small\", or \"large\""))
    }
  }

  is.whole <- function(x,
                       tol = .Machine$double.eps^0.5) {
    if(is.numeric(x)){
      return(abs(x - round(x)) < tol)
    }else{
      return(FALSE)
    }
  }

  if (!is.numeric(k) || is.na(k)) {
    stop("'k' must be an integer")
  } else {
    if (!is.whole(k)) {
      suppressWarnings(k <- as.integer(k))
      if(is.na(k)){stop("'k' must be an integer")}
      warning(paste("provided value for 'k' not an integer, 'k' =",
                    as.character(k),
                    "was used instead"))
    }
  }
  if (!is.numeric(p) || is.na(p)) {
    stop("'p' must be an integer")
  } else {
    if (!is.whole(p)) {
      suppressWarnings(p <- as.integer(p))
      if(is.na(p)){stop("'p' must be an integer")}
      warning(paste("provided value for 'p' not an integer, 'p' =",
                    as.character(p),
                    "was used instead"))
      if (p >= k) {stop("'p' must be smaller than 'k'")}
    } else {
      if (p >= k) {stop("'p' must be smaller than 'k'")}
    }
  }

  q <- k - p

  if (!is.numeric(N) || is.na(N)) {
    stop("'N' must be an integer")
  } else {
    if (!is.whole(N)) {
      suppressWarnings(N <- as.integer(N))
      if(is.na(N)){stop("'N' must be an integer")}
      warning(paste("provided value for 'N' not an integer, 'N' =",
                    as.character(N),
                    "was used instead"))
      if (k >= N) {stop("'N' must be larger than 'k'")}
    } else {
      if (k >= N) {stop("'N' must be larger than 'k'")}
    }
  }

  frsq <- (N-p-q) / q * R2p / (1-R2p)

  if (!is.numeric(nu) || is.na(nu)) {
    stop("'nu' must be numeric")
  } else if (nu < 3) {
    stop("'nu' must at least be 3")
  }

  if(!is.infinite(nu)){
    hyper <- function(fsq,
                      q,
                      nu,
                      f) {

      temp <- f^2 * (nu + q - 2)/2

      gamma((q + nu) / 2) / gamma(nu / 2) / gamma(q / 2) *
        2 * (nu - 2) / q / (nu-2 + q) / f^2 *
        fsq^(q/2) * temp^(nu/2) * (temp + fsq)^(-(nu+q)/2)

    }
  } else {
    hyper <- function(fsq,
                      q,
                      nu,
                      f) {

      (2 * exp(-(fsq / f^2)) * (f^2)^(-1 - q / 2) * fsq^(q / 2)) /
        (q * gamma(q / 2))

    }
  }

  pseudoBayes <- Vectorize(function(fsq, N, p, q, frsq, nu, f_) {

    f <- f_

    ncp <- (N-p) * fsq

    out <- (df(frsq, q, N-p-q, ncp) / df(frsq, q, N-p-q, 0)) *
      hyper(fsq, q=q, nu=nu, f=f)

    #out <- ifelse(is.infinite(out) | is.na(out) | is.nan(out), 0, out)

    return(out)

  }, "fsq")

  suppressWarnings({
    z <- integrate(f=pseudoBayes,
                   N=N, p=p, q=q, frsq=frsq, nu=nu, f_=f,
                   lower=0, upper=Inf)
  })

  emBayesFactor_out <- list(method = paste("Moment Bayes Factor:",
                                  "Multiple Linear Regression",
                                  sep=" "),
                   BayesFactor = z$value, accuracy = z$abs.error,
                   Cohen = "Cohen's f^2", Cohen_size = f^2,
                   rscale = NULL, nu = nu,
                   Stat = "partial R^2", Stat_size = R2p,
                   alternative = NULL,
                   note = "alternative over null")

  class(emBayesFactor_out) <- "emBayesFactor_out"

  return(emBayesFactor_out)

}

#' @rdname mreg
#' @export
mreg_FStat <- function(
    FStat, N, k, p,
    Cohen_f2 = c("medium", "small", "large"),
    nu = 5 + (k - p)) {

  if (!is.numeric(FStat) || is.na(FStat)) {
    stop("'FStat' must be numeric")
  } else if (FStat < 0) {
    stop("'FStat' must not be negative")
  }

  if (missing(Cohen_f2)) {
    Cohen_f2 = NULL
  }

  is.whole <- function(x,
                       tol = .Machine$double.eps^0.5) {
    if(is.numeric(x)){
      return(abs(x - round(x)) < tol)
    }else{
      return(FALSE)
    }
  }

  if (!is.numeric(k) || is.na(k)) {
    stop("'k' must be an integer")
  } else {
    if (!is.whole(k)) {
      suppressWarnings(k <- as.integer(k))
      if(is.na(k)){stop("'k' must be an integer")}
      warning(paste("provided value for 'k' not an integer, 'k' =",
                    as.character(k),
                    "was used instead"))
    }
  }
  if (!is.numeric(p) || is.na(p)) {
    stop("'p' must be an integer")
  } else {
    if (!is.whole(p)) {
      suppressWarnings(p <- as.integer(p))
      if(is.na(p)){stop("'p' must be an integer")}
      warning(paste("provided value for 'p' not an integer, 'p' =",
                    as.character(p),
                    "was used instead"))
      if (p >= k) {stop("'p' must be smaller than 'k'")}
    } else {
      if (p >= k) {stop("'p' must be smaller than 'k'")}
    }
  }

  if (!is.numeric(N) || is.na(N)) {
    stop("'N' must be an integer")
  } else {
    if (!is.whole(N)) {
      suppressWarnings(N <- as.integer(N))
      if(is.na(N)){stop("'N' must be an integer")}
      warning(paste("provided value for 'N' not an integer, 'N' =",
                    as.character(N),
                    "was used instead"))
      if (k >= N) {stop("'N' must be larger than 'k'")}
    } else {
      if (k >= N) {stop("'N' must be larger than 'k'")}
    }
  }

  df_numer <- k - p
  df_denom <- N - k
  R2p <- FStat * df_numer / (FStat * df_numer + df_denom)

  z <- mreg_R2p(R2p=R2p, N=N, k=k, p=p, Cohen_f2=Cohen_f2,
                nu=nu)

  emBayesFactor_out <- list(method = paste("Moment Bayes Factor:",
                                  "F Test (Multiple Linear Regression)",
                                  sep=" "),
                   BayesFactor = z$BayesFactor, accuracy = z$accuracy,
                   Cohen = "Cohen's f^2", Cohen_size = z$Cohen_size,
                   rscale = NULL, nu = z$nu,
                   Stat = "F statistic", Stat_size = FStat,
                   alternative = NULL,
                   note = "alternative over null")

  class(emBayesFactor_out) <- "emBayesFactor_out"

  return(emBayesFactor_out)

}

#' @rdname mreg
#' @export
mreg_tStat <- function(tStat, N, k,
                       Cohen_f2 = c("medium", "small", "large"),
                       nu = 5,
                       alternative = c("two.sided", "less", "greater")) {

  alternative <- match.arg(alternative)

  if (!is.numeric(tStat) || is.na(tStat)) {
    stop("'tStat' must be numeric")
  }

  if (length(Cohen_f2)>1 && !missing(Cohen_f2)){
    warning("'Cohen_f2' only takes one value, the first one was used")
    Cohen_f2 <- Cohen_f2[1]
  }
  if(missing(Cohen_f2) || is.null(Cohen_f2)){
    f <- sqrt(0.15)
  } else {
    if (is.numeric(Cohen_f2) && !is.na(Cohen_f2)){
      f2 <- as.numeric(Cohen_f2)
      if (f2 <= 0) {
        warning(paste("'Cohen_f2' must be positive,",
                      "'Cohen_f2' = \"medium\" used instead"))
        f <- sqrt(0.15)
      } else {
        f <- sqrt(f2)
      }
    } else if (is.character(Cohen_f2) && !is.na(Cohen_f2)) {
      f <- switch(Cohen_f2,
                  medium=sqrt(0.15),
                  small=sqrt(0.02),
                  large=sqrt(0.35),
                  stop(paste("'Cohen_f2' must be either numeric",
                             "or alternatively set to one of the options",
                             "\"medium\", \"small\", or \"large\""))
      )
    } else {
      stop(paste("'Cohen_f2' must be either numeric",
                 "or alternatively set to one of the options",
                 "\"medium\", \"small\", or \"large\""))
    }
  }

  if (alternative == "less") {
    Cohen.d <- -f
  } else {
    Cohen.d <- f
  }

  is.whole <- function(x,
                       tol = .Machine$double.eps^0.5) {
    if(is.numeric(x)){
      return(abs(x - round(x)) < tol)
    }else{
      return(FALSE)
    }
  }

  if (!is.numeric(k) || is.na(k)) {
    stop("'k' must be an integer")
  } else {
    if (!is.whole(k)) {
      suppressWarnings(k <- as.integer(k))
      if(is.na(k)){stop("'k' must be an integer")}
      warning(paste("provided value for 'k' not an integer, 'k' =",
                    as.character(k),
                    "was used instead"))
    }
  }

  if (!is.numeric(N) || is.na(N)) {
    stop("'N' must be an integer")
  } else {
    if (!is.whole(N)) {
      suppressWarnings(N <- as.integer(N))
      if(is.na(N)){stop("'N' must be an integer")}
      warning(paste("provided value for 'N' not an integer, 'N' =",
                    as.character(N),
                    "was used instead"))
      if (k >= N) {stop("'N' must be larger than 'k'")}
    } else {
      if (k >= N) {stop("'N' must be larger than 'k'")}
    }
  }

  N1 <- N - k + 1
  N2 <- 0
  nu_t <- N1 - 1
  N <- N1

  d <- Cohen.d

  if (!is.numeric(nu) || is.na(nu)) {
    stop("'nu' must be numeric")
  } else if (nu < 3) {
    stop("'nu' must be at least 3")
  }

  if (alternative == "two.sided") {

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

  emBayesFactor_out <- list(method = paste("Moment Bayes Factor:",
                                  "t Test (Single Regression Coefficient)",
                                  sep=" "),
                   BayesFactor = z$value, accuracy = z$abs.error,
                   Cohen = "Cohen's f^2", Cohen_size = d^2,
                   rscale = NULL, nu = nu,
                   Stat = "t statistic", Stat_size = tStat,
                   alternative = alternative,
                   note = "alternative over null")

  class(emBayesFactor_out) <- "emBayesFactor_out"

  return(emBayesFactor_out)

}

#' @rdname mreg
#' @usage NULL
#' @order 0
mreg_dummyFunctionName <-
  function(R2p, FStat, tStat, N, k, p,
           Cohen_f2 = c("medium", "small", "large"),
           nu = 5,
           alternative = c("two.sided", "less", "greater")) { NULL }
