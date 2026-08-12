# Pull the response and the model matrix out of an evaluated model frame.
# Neither model has an intercept, because a constant in X cannot be told apart
# from a rescaling of the baseline function.
pcd_model_data <- function(mf) {
  mt <- attr(mf, "terms")
  y <- stats::model.response(mf)
  if (!is.pcd(y)) {
    stop("The left hand side of `formula` must be a `pcd()` object.",
         call. = FALSE)
  }
  X <- stats::model.matrix(mt, mf)
  # Keep the contrasts before subsetting: dropping a column drops the attribute
  # with it, and predict() needs the same coding the fit used.
  contrasts <- attr(X, "contrasts")
  intercept <- match("(Intercept)", colnames(X), 0L)
  if (intercept > 0L) X <- X[, -intercept, drop = FALSE]
  list(y = y, X = X, terms = mt, contrasts = contrasts,
       xlevels = stats::.getXlevels(mt, mf), na.action = attr(mf, "na.action"))
}

pcd_init <- function(init, p) {
  if (is.null(init)) return(numeric(p))
  if (length(init) != p) {
    stop("`init` has length ", length(init), " but the model has ", p,
         " coefficient(s).", call. = FALSE)
  }
  as.numeric(init)
}

#' Tuning parameters for the fitting algorithms
#'
#' @param maxit Maximum number of iterations. For the rate model these are EM
#'   iterations; for the means model, Newton iterations of the estimating
#'   equation, which converges in a handful.
#' @param reltol Convergence tolerance on the largest relative change in the
#'   coefficients, and for the rate model in the baseline jump sizes too.
#' @param accelerate Rate model only: whether to extrapolate the EM iterations
#'   by the SQUAREM scheme of Varadhan and Roland (2008). This reaches the same
#'   fixed point in far fewer passes; set it to `FALSE` to run plain EM.
#' @param profile_maxit,profile_reltol Limits for the inner runs that maximise
#'   over the baseline with the coefficients held fixed, used only by the
#'   profile likelihood covariance. These can be looser than `reltol` without
#'   loss: the profile log likelihood is stationary in the baseline at its
#'   maximum, so an error of size \eqn{\epsilon} in the inner solution moves it
#'   by only \eqn{O(\epsilon^2)}.
#' @param profile_h Step size for the numerical derivative of the profile log
#'   likelihood. The default, `NULL`, uses \eqn{5 n^{-1/2}} as recommended by
#'   Zeng et al. (2016) and used in the paper.
#'
#' @return A list of control values.
#'
#' @references
#' Varadhan, R. and Roland, C. (2008). Simple and globally convergent methods
#' for accelerating the convergence of any EM algorithm.
#' *Scandinavian Journal of Statistics* **35**, 335--353.
#'
#' @examples
#' pcdreg_control(reltol = 1e-9)
#' @export
pcdreg_control <- function(maxit = 2000L, reltol = 1e-7, accelerate = TRUE,
                           profile_maxit = 20000L, profile_reltol = 1e-6,
                           profile_h = NULL) {
  stopifnot(
    length(maxit) == 1L, maxit >= 1,
    length(reltol) == 1L, reltol > 0,
    length(accelerate) == 1L, is.logical(accelerate), !is.na(accelerate),
    length(profile_maxit) == 1L, profile_maxit >= 1,
    length(profile_reltol) == 1L, profile_reltol > 0,
    is.null(profile_h) || (length(profile_h) == 1L && profile_h > 0)
  )
  list(maxit = as.integer(maxit), reltol = as.numeric(reltol),
       accelerate = accelerate,
       profile_maxit = as.integer(profile_maxit),
       profile_reltol = as.numeric(profile_reltol),
       profile_h = profile_h)
}

#' Semiparametric regression for panel count data
#'
#' Fits either of the two standard semiparametric models to panel count data,
#' that is, to recurrent events observed only at intermittent examination times.
#' Covariates may vary over time.
#'
#' @param formula A formula whose left hand side is a [pcd()] object, for
#'   example `pcd(id, tstart, tstop, count) ~ x1 + x2`.
#' @param data A data frame containing the variables in `formula`.
#' @param model Which model to fit. `"rate"` is the proportional rate model and
#'   the default; `"mean"` is the proportional means model.
#' @param subset,na.action Handled as in [stats::lm()]. Note that missing
#'   covariate values delete rows, which can break the contiguity of a subject's
#'   intervals; resolve them before fitting.
#' @param control A list from [pcdreg_control()].
#' @param profile Rate model only: whether to also compute the profile
#'   likelihood covariance. It is provided for comparison with the literature
#'   and costs one extra baseline-only fit per coefficient, so it is off by
#'   default.
#' @param init Optional starting values for the coefficients.
#'
#' @section The two models:
#' \deqn{\textrm{rate:} \quad E[dN(t) \mid X(t)] = \exp(\beta' X(t)) \, d\Lambda(t)}
#' \deqn{\textrm{mean:} \quad E[N(t) \mid X(t)] = \mu(t) \exp(\beta' X(t))}
#'
#' The **rate model** is fitted by nonparametric maximum likelihood, treating
#' \eqn{\Lambda} as a step function with a jump at each observed examination
#' time. The likelihood maximised is the one a nonhomogeneous Poisson process
#' would give, but that is a working device rather than an assumption about the
#' data: the estimator stays consistent and asymptotically normal when it fails,
#' and the default covariance estimator stays valid.
#'
#' The **means model** is fitted by the estimating equation of Hu, Sun and Wei
#' (2003). Only the examination times enter, so no covariate values between them
#' are needed and it is much cheaper. There is no likelihood behind it, so it
#' reports no log likelihood, `profile = TRUE` is unavailable, and the only
#' covariance is the sandwich.
#'
#' The two are not reparametrisations of each other when covariates vary over
#' time, so their coefficients answer different questions rather than estimating
#' a common quantity: \eqn{\beta} acts on the instantaneous rate in one and on
#' the cumulative mean in the other. A further practical difference is that
#' nothing constrains the fitted \eqn{\hat\mu} to increase, and with fluctuating
#' covariates it generally does not, which is the drawback of the means model
#' that motivates the rate model.
#'
#' @section Specifying the model:
#' The left hand side is always a [pcd()] object, which carries the subject, the
#' interval and the event count. The right hand side is an ordinary model
#' formula, so interactions, transformations, factors and `.` all behave as
#' usual:
#'
#' ```
#' pcd(id, tstart, tstop, count) ~ x1 + x2
#' pcd(id, tstart, tstop, count) ~ x1 * x2 + log(x3)
#' pcd(id, time, count)          ~ .
#' ```
#'
#' The dot expands to every column that is not part of the response, so the
#' identifier, times and counts are excluded automatically.
#'
#' No intercept is fitted. A constant column cannot be distinguished from a
#' rescaling of the baseline, so it is dropped and the baseline is reported
#' separately by [baseline()]. Adding `- 1` changes nothing. One consequence is
#' that shifting a covariate by a constant \eqn{c} leaves \eqn{\beta} alone and
#' rescales the baseline by \eqn{e^{-c\beta}}: centring covariates moves the
#' baseline, not the coefficients.
#'
#' `vignette("data-preparation")` covers the data layouts in full, including the
#' conversion from cumulative counts, which is the mistake most worth avoiding.
#'
#' @section Covariance estimation:
#' Available through [vcov.pcdfit()]. `"robust"` is the sandwich
#' \eqn{\Omega^{-1} S \Omega^{-1} / n}, which does not rely on the Poisson
#' assumption, is the default, and is the only option for the means model.
#'
#' For the rate model two more are available. `"information"` is
#' \eqn{S^{-1}/n}, valid only under the Poisson assumption but a fast substitute
#' for the profile likelihood when it holds. `"profile"` is the numerical
#' profile likelihood estimator of Murphy and van der Vaart (2000), computed
#' only when `profile = TRUE`; it is there for comparison with the literature
#' rather than for routine use, since it rests on a numerical derivative with a
#' step of order \eqn{n^{-1/2}} that is coarse at small sample sizes.
#'
#' When the counts are overdispersed relative to Poisson, which is the usual
#' situation, the latter two understate the standard errors, sometimes severely.
#'
#' @return An object of class `"pcdfit"`, with components including
#'   `coefficients`, `model` (which model was fitted), `baseline` (a data frame
#'   of the estimated baseline function), `Omega` and `S` (the two matrices
#'   behind the sandwich), `vcov` (a list of the available covariance
#'   estimates), and convergence information. Rate model fits also carry
#'   `loglik`.
#'
#' @references
#' Sun, D., Guo, Y., Li, Y., Tu, W. and Sun, J. (2024).
#' A robust approach for regression analysis of panel count data.
#' *Bernoulli* **30**(4), 3251--3275. \doi{10.3150/23-BEJ1713}
#'
#' Hu, X. J., Sun, J. and Wei, L. J. (2003). Regression parameter estimation
#' from panel counts. *Scandinavian Journal of Statistics* **30**(1), 25--43.
#' \doi{10.1111/1467-9469.00316}
#'
#' Murphy, S. A. and van der Vaart, A. W. (2000). On profile likelihood.
#' *Journal of the American Statistical Association* **95**, 449--465.
#'
#' @examples
#' set.seed(1)
#' d <- r_panel_count(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
#'
#' fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
#' summary(fit)
#'
#' # The same data under the means model. The coefficients are not comparable
#' # term by term, because the two models act on different quantities.
#' mfit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d,
#'                model = "mean")
#' cbind(rate = coef(fit), mean = coef(mfit))
#'
#' @seealso [pcd()], [vcov.pcdfit()], [baseline()], [r_panel_count()]
#' @export
pcdreg <- function(formula, data, model = c("rate", "mean"), subset, na.action,
                   control = pcdreg_control(), profile = FALSE, init = NULL) {
  model <- match.arg(model)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  md <- pcd_model_data(mf)
  p <- ncol(md$X)
  init <- pcd_init(init, p)

  if (model == "mean" && isTRUE(profile)) {
    stop("The means model is fitted by an estimating equation rather than a ",
         "likelihood, so `profile = TRUE` is only available for ",
         "`model = \"rate\"`.", call. = FALSE)
  }

  # Only the rate model needs the covariate trajectory expanded onto the pooled
  # grid, and that expansion is the expensive part of preparing the data.
  d <- prepare_panel(md$y, md$X, expand = (model == "rate"))
  part <- if (model == "rate") fit_rate(d, init, control, profile)
          else fit_mean(d, init, control)

  beta <- drop(part$beta)
  names(beta) <- colnames(md$X)
  Omega <- part$Omega
  S <- part$S
  dimnames(Omega) <- dimnames(S) <- list(names(beta), names(beta))

  # With no covariates there is nothing to make inference about.
  vcovs <- if (p > 0) {
    c(list(robust = robust_vcov(Omega, S, d$n)), part$extra_vcov)
  } else {
    list(robust = Omega)
  }
  for (nm in names(vcovs)) {
    if (!is.null(vcovs[[nm]])) dimnames(vcovs[[nm]]) <- dimnames(S)
  }

  structure(
    c(list(
      coefficients = beta,
      model = model,
      vcov = vcovs,
      Omega = Omega,
      S = S,
      scores = part$scores,
      baseline = part$baseline,
      n = d$n,
      nexam = d$nexam,
      nevent = sum(d$dN),
      ngrid = d$K,
      control = control,
      call = cl,
      terms = md$terms,
      xlevels = md$xlevels,
      contrasts = md$contrasts,
      na.action = md$na.action
    ), part$info),
    class = "pcdfit"
  )
}

# The rate model: EM, then the two covariance matrices of Section 4, and
# optionally the profile likelihood covariance.
fit_rate <- function(d, init, control, profile) {
  lambda0 <- rep.int(1 / max(d$K, 1L), d$K)
  fit <- em_fit_cpp(d$X, d$subj, d$grid, d$panel, d$dN, d$panelsubj,
                    d$n, d$K, init, lambda0,
                    control$maxit, control$reltol, control$accelerate)
  if (!fit$converged) {
    warning("The EM algorithm did not converge in ", control$maxit,
            " iterations (relative change ", format(fit$criterion, digits = 3),
            "). Consider raising `maxit` in `pcdreg_control()`.",
            call. = FALSE)
  }

  beta <- drop(fit$beta)
  lambda <- drop(fit$lambda)
  cov_parts <- covariance_cpp(d$X, d$subj, d$grid, d$panel, d$dN, d$panelsubj,
                              d$n, d$K, beta, lambda)

  extra <- list(information = if (length(beta))
    information_vcov(cov_parts$S, d$n))
  if (isTRUE(profile) && length(beta)) {
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
    extra$profile <- profile_vcov(prof$gradient)
  }

  list(beta = beta, Omega = cov_parts$Omega, S = cov_parts$S,
       scores = cov_parts$scores, extra_vcov = extra,
       baseline = data.frame(time = d$times, jump = lambda,
                             cumrate = cumsum(lambda)),
       info = list(loglik = fit$loglik, iterations = fit$iterations,
                   passes = fit$passes, converged = fit$converged,
                   criterion = fit$criterion))
}

# The means model: Newton on the estimating equation, then the sandwich.
fit_mean <- function(d, init, control) {
  fit <- mean_fit_cpp(d$exam_X, d$exam_subj, d$exam_grid, d$cN, d$n, d$K,
                      init, control$maxit, control$reltol)
  if (!fit$converged) {
    warning("The estimating equation did not converge in ", control$maxit,
            " iterations (relative change ", format(fit$criterion, digits = 3),
            "). Consider raising `maxit` in `pcdreg_control()`.",
            call. = FALSE)
  }
  list(beta = drop(fit$beta), Omega = fit$Omega, S = fit$S,
       scores = fit$scores, extra_vcov = list(),
       baseline = data.frame(time = d$times, mean = drop(fit$mu)),
       info = list(iterations = fit$iterations, converged = fit$converged,
                   criterion = fit$criterion))
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
