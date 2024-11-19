#' @name eANOVA
#' @title Effect-Size Bayes Factors for Analysis of Variance (ANOVA)
#' @description
#' Computes Bayes factors (alternative over null hypothesis)
#' with effect-size priors for analysis of variance (ANOVA)
#' based on observed partial \eqn{{\eta}^2} (\eqn{{\eta}^2_p})
#' or the observed test statistic \eqn{F}.
#'
#'
#' @param eta2p observed partial coefficient of determination \eqn{{\eta}^2_p}.
#' @param FStat observed value of \eqn{F} statistic.
#' @param df_numer numerator degrees of freedom.
#' @param df_denom denominator degrees of freedom.
#' @param Cohen_f focused-on effect size (Cohen's \eqn{f}).
#' Must be either numeric, or alternatively set to one of the options
#' \code{"medium" = .25} (default), \code{"small" = .10}, \code{"large" = .40}.
#' @param rscale prior scale.
#' If not specified, it is set internally to the recommended value.
#' @param nu prior degrees of freedom.
#' If not specified, it is set internally to the recommended value (\code{nu = 3}).
#'
#'
#' @details
#' The Bayes factor provided by \code{eANOVA} compares
#' the alternative hypothesis with the null hypothesis
#' that there is no true effect.
#'
#' \code{Cohen_f} must not be negative.
#'
#' If \code{rscale} is not specified, \code{nu} must be greater than 2.
#'
#' If \code{Cohen_f = 0}, \code{rscale} must be explicitly specified.
#'
#'
#' @return An object of class \code{emBayesFactor_out} containing the following components:
#'
#' \item{method}{the focussed-on analysis.}
#' \item{BayesFactor}{value of the Bayes factor.}
#' \item{accuracy}{estimate of the modulus of the absolute error.}
#' \item{Cohen}{the focussed-on effect size.}
#' \item{Cohen_size}{value of the focussed-on effect size.}
#' \item{rscale}{prior scale.}
#' \item{nu}{prior degrees of freedom.}
#' \item{Stat}{provided statistic.}
#' \item{Stat_size}{value of the provided statistic.}
#' \item{alternative}{the specified alternative hypothesis}
#'
#'
#' @references
#' Klauer, K. C., Meyer-Grant, C. G., & Kellen, D (in press).
#' On Bayes Factors for Hypotheses Tests
#' \emph{Psychonomic Bulletin & Review}.
#' \[PsyArXiv preprint: \href{https://doi.org/10.31234/osf.io/ykp29}{https://doi.org/10.31234/osf.io/ykp29}\]
#'
#'
#' @author Constantin G. Meyer-Grant
#' (\email{constantin.meyer-grant@psychologie.uni-freiburg.de})
#'
#'
#' @seealso [integrate], [aov], \link[BayesFactor]{anovaBF}
#'
#'
#' @examples
#' ## Example from Klauer et al. (in press)
#'
#' ## Compute effect-size Bayes factor based on 'Fstat'
#' eANOVA_FStat(
#'   FStat = 16.73,
#'   df_numer = 1,
#'   df_denom = 36,
#'   Cohen_f = "small"
#' )


#' @rdname eANOVA
#' @export
eANOVA_eta2p <- function(eta2p, df_numer, df_denom,
                         Cohen_f = c("medium", "small", "large"),
                         rscale = NULL, nu = 3) {

  if (!is.numeric(eta2p) || is.na(eta2p)) {
    stop("'R2p' must be numeric")
  } else if (eta2p < 0 || eta2p > 1) {
    stop("'R2p' must be greater than zero and smaller than one")
  }

  is.whole <- function(x,
                       tol = .Machine$double.eps^0.5) {
    if(is.numeric(x)){
      return(abs(x - round(x)) < tol)
    }else{
      return(FALSE)
    }
  }

  if (!is.numeric(df_numer) || is.na(df_numer)) {
    stop("'df_numer' must be an integer")
  } else {
    if (!is.whole(df_numer)) {
      suppressWarnings(df_numer <- as.integer(df_numer))
      if(is.na(df_numer)){stop("'df_numer' must be an integer")}
      warning(paste("provided value for 'df_numer' not an integer,",
                    "'df_numer' =",
                    as.character(df_numer),
                    "was used instead"))
    }
  }
  if (!is.numeric(df_denom) || is.na(df_denom)) {
    stop("'df_denom' must be an integer")
  } else {
    if (!is.whole(df_denom)) {
      suppressWarnings(df_denom <- as.integer(df_denom))
      if(is.na(df_denom)){stop("'df_denom' must be an integer")}
      warning(paste("provided value for 'df_denom' not an integer,",
                    "'df_denom' =",
                    as.character(df_denom),
                    "was used instead"))
    }
  }

  FStat <- df_denom * eta2p / (df_numer - df_numer * eta2p)

  z <- eANOVA_FStat(FStat=FStat, df_numer=df_numer, df_denom=df_denom,
                    Cohen_f=Cohen_f,
                    rscale=rscale, nu=nu)

  emBayesFactor_out <- list(method = paste("Effect-Size Bayes Factor:",
                                 "ANOVA",
                                 sep=" "),
                  BayesFactor = z$BayesFactor, accuracy = z$accuracy,
                  Cohen = "Cohen's f", Cohen_size = z$Cohen_size,
                  rscale = z$rscale, nu = z$nu,
                  Stat = "partial eta^2", Stat_size = eta2p,
                  alternative = NULL,
                  note = "alternative over null")

  class(emBayesFactor_out) <- "emBayesFactor_out"

  return(emBayesFactor_out)

}

#' @rdname eANOVA
#' @export
eANOVA_FStat <- function(FStat, df_numer, df_denom,
                         Cohen_f = c("medium", "small", "large"),
                         rscale = NULL, nu = 3
) {

  if (!is.numeric(FStat) || is.na(FStat)) {
    stop("'FStat' must be numeric")
  } else if (FStat < 0) {
    stop("'FStat' must not be negative")
  }

  if (length(Cohen_f)>1 && !missing(Cohen_f)){
    warning("'Cohen_f' only takes one value, the first one was used")
    Cohen_f <- Cohen_f[1]
  }
  if(missing(Cohen_f) || is.null(Cohen_f)){
    f <- 0.25
  } else {
    if (is.numeric(Cohen_f) && !is.na(Cohen_f)){
      f <- as.numeric(Cohen_f)
      if (f < 0) {
        warning(paste("'Cohen_f' must not be negative,",
                      "'Cohen_f' = \"medium\" used instead"))
        f <- 0.25
      }
      if (is.null(rscale) && f == 0) {
        stop("if 'Cohen_f' is set to zero, 'rscale' must be set by hand")
      }
    } else if (is.character(Cohen_f) && !is.na(Cohen_f)) {
      f<-switch(Cohen_f,
                medium=0.25,
                small=0.10,
                large=0.40,
                stop(paste("'Cohen_f' must be either numeric",
                           "or alternatively set to one of the options",
                           "\"medium\", \"small\", or \"large\""))
      )
    } else {
      stop(paste("'Cohen_f' must be either numeric",
                 "or alternatively set to one of the options",
                 "\"medium\", \"small\", or \"large\""))
    }
  }

  if (!is.numeric(nu) || is.na(nu)) {
    stop("'nu' must be numeric")
  } else if (nu <= 0) {
    stop("'nu' must be greater than zero")
  }

  is.whole <- function(x,
                       tol = .Machine$double.eps^0.5) {
    if(is.numeric(x)){
      return(abs(x - round(x)) < tol)
    }else{
      return(FALSE)
    }
  }

  if (!is.numeric(df_numer) || is.na(df_numer)) {
    stop("'df_numer' must be an integer")
  } else {
    if (!is.whole(df_numer)) {
      suppressWarnings(df_numer <- as.integer(df_numer))
      if(is.na(df_numer)){stop("'df_numer' must be an integer")}
      warning(paste("provided value for 'df_numer' not an integer,",
                    "'df_numer' =",
                    as.character(df_numer),
                    "was used instead"))
    }
  }
  if (!is.numeric(df_denom) || is.na(df_denom)) {
    stop("'df_denom' must be an integer")
  } else {
    if (!is.whole(df_denom)) {
      suppressWarnings(df_denom <- as.integer(df_denom))
      if(is.na(df_denom)){stop("'df_denom' must be an integer")}
      warning(paste("provided value for 'df_denom' not an integer,",
                    "'df_denom' =",
                    as.character(df_denom),
                    "was used instead"))
    }
  }

  q <- df_numer
  N_minus_p <- df_numer + df_denom

  if (is.null(rscale)) {
    if (nu <= 2) {
      stop("'nu' must be greater than two if 'rscale' is not specified")
    }
    if (is.infinite(nu)) {
      rscale <- f / sqrt(q)
    }else{
      rscale <- sqrt(nu-2) / sqrt(nu * q) * f
    }
  } else if (!is.numeric(rscale) || is.na(rscale)) {
    stop("'rscale' must be numeric")
  } else if (rscale <= 0) {
    stop("'rscale' must be positive")
  }

  R2p <- FStat / ((N_minus_p - q) / q + FStat)

  frsq <- (N_minus_p - q) / q * R2p / (1 - R2p)

  if(!is.infinite(nu)){
    hyper <- function(fsq,
                      q,
                      nu,
                      rscale,
                      f) {

      gamma((q + nu) / 2) / gamma(nu / 2) /gamma(q / 2) *
        (nu * rscale^2)^(nu / 2) * fsq^(q / 2 - 1) *
        (nu * rscale^2 + f^2 + fsq)^(-nu / 2 - q / 2) *
        hypergeo::genhypergeo(c((nu + q) / 4, (2 + nu + q) / 4),
                              q / 2, 4 * f^2 * fsq / (nu * rscale^2 + f^2 + fsq)^2)

    }
  }else{
    hyper <- function(fsq,
                      q,
                      nu,
                      rscale,
                      f) {

      (2 * rscale^2)^(-(q / 2)) *
        exp(-((fsq + f^2)/(2 * rscale^2)))* fsq^(-1 + q / 2) *
        hypergeo::genhypergeo(U=NULL, L=q / 2, z=(fsq * f^2)/(4 *rscale^4)) /
        gamma(q / 2)

    }
  }

  pseudoBayes <- Vectorize(function(fsq, N_minus_p, q, frsq, nu, rscale, f_) {

    f <- f_

    ncp <- N_minus_p * fsq

    out <- df(frsq, q, N_minus_p-q, ncp) / df(frsq, q, N_minus_p-q, 0) *
      hyper(fsq, q=q, nu=nu, rscale=rscale, f=f)

    #out <- ifelse(is.infinite(out) | is.na(out) | is.nan(out), 0, out)

    return(out)

  }, "fsq")

  suppressWarnings({
    z <- integrate(f=pseudoBayes,
                   N_minus_p=N_minus_p, q=q, frsq=frsq, nu=nu, rscale=rscale, f_=f,
                   lower=0, upper=Inf)
  })

  emBayesFactor_out <- list(method = paste("Effect-Size Bayes Factor:",
                                 "F Test (ANOVA)",
                                 sep=" "),
                  BayesFactor = z$value, accuracy = z$abs.error,
                  Cohen = "Cohen's f", Cohen_size = f,
                  rscale = rscale, nu = nu,
                  Stat = "F statistic", Stat_size = FStat,
                  alternative = NULL,
                  note = "alternative over null")

  class(emBayesFactor_out) <- "emBayesFactor_out"

  return(emBayesFactor_out)

}

#' @rdname eANOVA
#' @usage NULL
#' @order 0
eANOVA_dummyFunctionName <-
  function(eta2p, FStat, df_numer, df_denom,
           Cohen_f = c("medium", "small", "large"),
           rscale = NULL, nu = 3) { NULL }
