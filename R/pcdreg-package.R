#' @keywords internal
#' @aliases pcdreg-package
#'
#' @description
#' Panel count data arise when a recurrent event process is observed only at
#' intermittent examination times, so that the number of events between
#' consecutive examinations is known but the event times themselves are not.
#' [pcdreg()] fits either of the two standard semiparametric models to such
#' data, allowing the covariates \eqn{X(t)} to vary over time:
#' \deqn{\textrm{rate:} \quad E[dN(t) \mid X(t)] = \exp(\beta' X(t)) \, d\Lambda(t)}
#' \deqn{\textrm{mean:} \quad E[N(t) \mid X(t)] = \mu(t) \exp(\beta' X(t))}
#'
#' The rate model treats the baseline cumulative rate \eqn{\Lambda}
#' nonparametrically and maximises the likelihood that a nonhomogeneous Poisson
#' process would imply, using an EM algorithm that augments the observed counts
#' with latent per-examination-time Poisson counts. The Poisson assumption is a
#' working device only: the estimator stays consistent and asymptotically normal
#' when it fails, and the default covariance estimator is robust to that
#' failure. The means model is fitted by an estimating equation and is provided
#' as the comparator.
#'
#' @references
#' Sun, D., Guo, Y., Li, Y., Tu, W. and Sun, J. (2024).
#' A robust approach for regression analysis of panel count data.
#' *Bernoulli* **30**(4), 3251--3275. \doi{10.3150/23-BEJ1713}
#'
#' @useDynLib pcdreg, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom ggplot2 autoplot
#' @importFrom stats coef confint logLik nobs predict vcov
#' @importFrom stats delete.response model.frame model.matrix model.response
#' @importFrom stats na.pass pnorm printCoefmat qnorm terms
"_PACKAGE"

# The first stats import above brings in the S3 generics themselves, so that the
# methods can be registered against them when the package is loaded without
# stats attached.  Without it the package fails to load with only its stated
# dependencies, which no session with stats already attached will reveal.

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

# Column names used inside ggplot2::aes(), which R CMD check cannot see are
# columns rather than stray globals.
utils::globalVariables(c("tstart", "tstop", "y", "bin", "time", "value"))
