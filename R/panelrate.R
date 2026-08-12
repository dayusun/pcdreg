#' Tuning parameters for the EM algorithm
#'
#' @param maxit Maximum number of EM iterations.
#' @param reltol Convergence tolerance. Iteration stops once the largest
#'   relative change in the coefficients and in the baseline jump sizes both
#'   fall below this value.
#' @param profile_maxit,profile_reltol Corresponding limits for the inner runs
#'   that maximise over the baseline with the coefficients held fixed, used only
#'   by the profile likelihood covariance. These can be looser than `reltol`
#'   without loss: the profile log likelihood is stationary in the baseline at
#'   its maximum, so an error of size \eqn{\epsilon} in the inner solution moves
#'   it by only \eqn{O(\epsilon^2)}.
#' @param profile_h Step size for the numerical derivative of the profile log
#'   likelihood. The default, `NULL`, uses \eqn{5 n^{-1/2}} as recommended by
#'   Zeng et al. (2016) and used in the paper.
#' @param accelerate Whether to extrapolate the EM iterations by the SQUAREM
#'   scheme of Varadhan and Roland (2008). This reaches the same fixed point in
#'   far fewer passes; set it to `FALSE` to run plain EM.
#'
#' @return A list of control values.
#'
#' @references
#' Varadhan, R. and Roland, C. (2008). Simple and globally convergent methods
#' for accelerating the convergence of any EM algorithm.
#' *Scandinavian Journal of Statistics* **35**, 335--353.
#'
#' @examples
#' panelrate_control(reltol = 1e-9)
#' @export
panelrate_control <- function(maxit = 2000L, reltol = 1e-7,
                              profile_maxit = 20000L, profile_reltol = 1e-6,
                              profile_h = NULL, accelerate = TRUE) {
  stopifnot(
    length(maxit) == 1L, maxit >= 1,
    length(reltol) == 1L, reltol > 0,
    length(profile_maxit) == 1L, profile_maxit >= 1,
    length(profile_reltol) == 1L, profile_reltol > 0,
    is.null(profile_h) || (length(profile_h) == 1L && profile_h > 0),
    length(accelerate) == 1L, is.logical(accelerate), !is.na(accelerate)
  )
  list(maxit = as.integer(maxit), reltol = as.numeric(reltol),
       profile_maxit = as.integer(profile_maxit),
       profile_reltol = as.numeric(profile_reltol),
       profile_h = profile_h, accelerate = accelerate)
}

#' Fit the proportional rate model to panel count data
#'
#' Fits
#' \deqn{E[dN(t) \mid X(t)] = \exp(\beta' X(t)) \, d\Lambda(t)}
#' by nonparametric maximum likelihood, treating \eqn{\Lambda} as a step
#' function with a jump at each observed examination time. Covariates may vary
#' over time.
#'
#' @param formula A formula whose left hand side is a [PanelCount()] object,
#'   for example `PanelCount(id, tstart, tstop, count) ~ x1 + x2`. No intercept
#'   is fitted: it is absorbed into the baseline \eqn{\Lambda}.
#' @param data A data frame containing the variables in `formula`.
#' @param subset,na.action Handled as in [stats::lm()]. Note that missing
#'   covariate values delete rows, which can break the contiguity of a subject's
#'   intervals; resolve them before fitting.
#' @param control A list from [panelrate_control()].
#' @param profile Whether to also compute the profile likelihood covariance.
#'   This is provided for comparison with the literature and costs one extra
#'   baseline-only fit per coefficient, so it is off by default.
#' @param init Optional starting values for the coefficients.
#'
#' @details
#' The likelihood maximised is the one a nonhomogeneous Poisson process would
#' give. That assumption is a working device: it makes the EM algorithm
#' available, but the estimator remains consistent and asymptotically normal
#' when it fails, and the default covariance estimator remains valid.
#'
#' Three covariance estimators are available through [vcov.pcdfit()].
#' `"robust"` is the sandwich \eqn{\Omega^{-1} S \Omega^{-1} / n}, which does not
#' rely on the Poisson assumption and is the default. `"information"` is
#' \eqn{S^{-1} / n}, valid only under that assumption but a fast substitute for
#' the profile likelihood when it holds. `"profile"` is the numerical profile
#' likelihood estimator of Murphy and van der Vaart (2000), computed only when
#' `profile = TRUE`. When the counts are overdispersed relative to Poisson, the
#' latter two understate the standard errors, sometimes severely.
#'
#' The profile estimator is included for comparison with the literature rather
#' than for routine use. It rests on a numerical derivative taken with a step of
#' order \eqn{n^{-1/2}}, which is coarse at small sample sizes, and it costs one
#' baseline-only fit per coefficient. `"information"` estimates the same matrix
#' far more cheaply and more stably.
#'
#' @return An object of class `"panelrate"`, inheriting from `"pcdfit"`, with
#'   components including `coefficients`, `baseline` (a data frame of jump sizes
#'   and the cumulative baseline rate), `Omega` and `S` (the two matrices of
#'   Section 4 of the paper), `vcov` (a list of the available covariance
#'   estimates), `loglik`, and convergence information.
#'
#' @references
#' Sun, D., Guo, Y., Li, Y., Tu, W. and Sun, J. (2024).
#' A robust approach for regression analysis of panel count data.
#' *Bernoulli* **30**(4), 3251--3275. \doi{10.3150/23-BEJ1713}
#'
#' Murphy, S. A. and van der Vaart, A. W. (2000). On profile likelihood.
#' *Journal of the American Statistical Association* **95**, 449--465.
#'
#' @examples
#' set.seed(1)
#' d <- r_panel_count(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
#' fit <- panelrate(PanelCount(id, tstart, tstop, count) ~ x1 + x2, data = d)
#' fit
#' summary(fit)
#'
#' @seealso [PanelCount()], [panelmean()], [vcov.pcdfit()], [r_panel_count()]
#' @export
panelrate <- function(formula, data, subset, na.action,
                      control = panelrate_control(), profile = FALSE,
                      init = NULL) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  md <- pcd_model_data(mf)
  X <- md$X
  p <- ncol(X)
  init <- pcd_init(init, p)

  d <- prepare_panel(md$y, X)
  lambda0 <- rep.int(1 / max(d$K, 1L), d$K)

  fit <- em_fit_cpp(d$X, d$subj, d$grid, d$panel, d$dN, d$panelsubj,
                    d$n, d$K, init, lambda0,
                    control$maxit, control$reltol, control$accelerate)
  if (!fit$converged) {
    warning("The EM algorithm did not converge in ", control$maxit,
            " iterations (relative change ", format(fit$criterion, digits = 3),
            "). Consider raising `maxit` in `panelrate_control()`.",
            call. = FALSE)
  }

  beta <- drop(fit$beta)
  names(beta) <- colnames(X)
  lambda <- drop(fit$lambda)

  cov_parts <- covariance_cpp(d$X, d$subj, d$grid, d$panel, d$dN, d$panelsubj,
                              d$n, d$K, beta, lambda)
  Omega <- cov_parts$Omega
  S <- cov_parts$S
  dimnames(Omega) <- dimnames(S) <- list(names(beta), names(beta))

  # With no covariates there is nothing to make inference about.
  vcovs <- if (p > 0) {
    list(robust = robust_vcov(Omega, S, d$n),
         information = information_vcov(S, d$n))
  } else {
    list(robust = Omega, information = S)
  }
  if (isTRUE(profile) && p > 0) {
    h <- control$profile_h
    if (is.null(h)) h <- 5 / sqrt(d$n)
    prof <- profile_gradient_cpp(d$X, d$subj, d$grid, d$panel, d$dN,
                                 d$panelsubj, d$n, d$K, beta, lambda, h,
                                 control$profile_maxit, control$profile_reltol)
    if (any(prof$converged == 0L)) {
      warning("Some baseline-only fits behind the profile likelihood ",
              "covariance did not converge; treat `type = \"profile\"` with ",
              "caution.", call. = FALSE)
    }
    vcovs$profile <- profile_vcov(prof$gradient)
  }
  for (nm in names(vcovs)) {
    if (!is.null(vcovs[[nm]])) dimnames(vcovs[[nm]]) <- dimnames(S)
  }

  structure(
    list(
      coefficients = beta,
      vcov = vcovs,
      Omega = Omega,
      S = S,
      scores = cov_parts$scores,
      baseline = data.frame(time = d$times, jump = lambda,
                            cumrate = cumsum(lambda)),
      loglik = fit$loglik,
      n = d$n,
      nexam = d$nexam,
      nevent = sum(d$dN),
      ngrid = d$K,
      iterations = fit$iterations,
      passes = fit$passes,
      converged = fit$converged,
      criterion = fit$criterion,
      control = control,
      call = cl,
      terms = md$terms,
      xlevels = md$xlevels,
      contrasts = md$contrasts,
      na.action = md$na.action
    ),
    class = c("panelrate", "pcdfit")
  )
}

# Var(betahat) = Omega^{-1} S Omega^{-1} / n.  Omega is negative definite, and
# the two sign changes cancel, so the positive definite -Omega is inverted.
robust_vcov <- function(Omega, S, n) {
  inv <- if (all(is.finite(Omega))) {
    tryCatch(solve(-Omega), error = function(e) NULL)
  }
  if (is.null(inv) || !all(is.finite(inv))) {
    warning("Omega is singular or not finite, so the robust covariance is ",
            "unavailable.", call. = FALSE)
    return(NULL)
  }
  symmetrise(inv %*% S %*% inv / n)
}

# Var(betahat) = S^{-1} / n, valid under the Poisson assumption.
information_vcov <- function(S, n) {
  inv <- if (all(is.finite(S))) tryCatch(solve(S), error = function(e) NULL)
  if (is.null(inv) || !all(is.finite(inv))) {
    warning("S is singular or not finite, so the information covariance is ",
            "unavailable.", call. = FALSE)
    return(NULL)
  }
  symmetrise(inv / n)
}

# Var(betahat) = (n Vhat)^{-1} with Vhat the empirical covariance of the profile
# score, so the n factors cancel against crossprod of the per-subject gradients.
profile_vcov <- function(gradient) {
  inv <- tryCatch(solve(crossprod(gradient)), error = function(e) NULL)
  if (is.null(inv)) {
    warning("The profile likelihood covariance is singular and unavailable.",
            call. = FALSE)
    return(NULL)
  }
  symmetrise(inv)
}

symmetrise <- function(m) (m + t(m)) / 2
