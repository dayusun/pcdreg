#' @keywords internal
#' @aliases panelrate-package
#'
#' @description
#' Panel count data arise when a recurrent event process is observed only at
#' intermittent examination times, so that the number of events between
#' consecutive examinations is known but the event times themselves are not.
#' `panelrate()` fits the semiparametric proportional rate model
#' \deqn{E[dN(t) \mid X(t)] = \exp(\beta' X(t)) \, d\Lambda(t)}
#' to such data, allowing the covariates \eqn{X(t)} to vary over time.
#'
#' Estimation treats the baseline cumulative rate \eqn{\Lambda} nonparametrically
#' and maximises the likelihood that a nonhomogeneous Poisson process would
#' imply, using an EM algorithm that augments the observed counts with latent
#' per-examination-time Poisson counts. The Poisson assumption is a working
#' device only: the estimator stays consistent and asymptotically normal when it
#' fails, and the default covariance estimator is robust to that failure.
#'
#' @references
#' Sun, D., Guo, Y., Li, Y., Tu, W. and Sun, J. (2024).
#' A robust approach for regression analysis of panel count data.
#' *Bernoulli* **30**(4), 3251--3275. \doi{10.3150/23-BEJ1713}
#'
#' @useDynLib panelrate, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' The S3 generics below are imported so that the methods can be registered
#' against them when the package is loaded without stats attached.
#'
#' @importFrom stats coef confint logLik nobs predict vcov
#' @importFrom stats delete.response model.frame model.matrix model.response
#' @importFrom stats na.pass pnorm printCoefmat qnorm terms
"_PACKAGE"
