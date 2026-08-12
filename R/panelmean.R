# Shared by both fitters: pull the response and the model matrix out of an
# evaluated model frame.  Neither model has an intercept, because a constant in
# X cannot be told apart from a rescaling of the baseline function.
pcd_model_data <- function(mf) {
  mt <- attr(mf, "terms")
  y <- stats::model.response(mf)
  if (!is.PanelCount(y)) {
    stop("The left hand side of `formula` must be a `PanelCount()` object.",
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

#' Tuning parameters for the means model
#'
#' @param maxit Maximum number of Newton iterations.
#' @param reltol Convergence tolerance on the relative change in the
#'   coefficients.
#'
#' @return A list of control values.
#' @examples
#' panelmean_control(reltol = 1e-10)
#' @export
panelmean_control <- function(maxit = 100L, reltol = 1e-9) {
  stopifnot(length(maxit) == 1L, maxit >= 1,
            length(reltol) == 1L, reltol > 0)
  list(maxit = as.integer(maxit), reltol = as.numeric(reltol))
}

#' Fit the proportional means model to panel count data
#'
#' Fits
#' \deqn{E[N(t) \mid X(t)] = \mu(t) \exp(\beta' X(t))}
#' by the estimating equation of Hu, Sun and Wei (2003). This is the comparator
#' the paper behind [panelrate()] uses, and it is provided so the two model
#' families can be compared on the same data.
#'
#' @param formula A formula whose left hand side is a [PanelCount()] object.
#' @param data A data frame containing the variables in `formula`.
#' @param subset,na.action Handled as in [stats::lm()].
#' @param control A list from [panelmean_control()].
#' @param init Optional starting values for the coefficients.
#'
#' @details
#' The estimating equation is
#' \deqn{U(\beta) = \sum_i \sum_j N_i(T_{ij}) \{ X_i(T_{ij}) - \bar x(T_{ij}) \},}
#' where \eqn{\bar x(t)} averages the covariates over the subjects examined at
#' \eqn{t}, weighted by \eqn{\exp(\beta' X)}, and \eqn{N_i(T_{ij})} is the
#' cumulative count. Only the examination times enter, so no covariate values
#' between them are needed and the fit is much cheaper than [panelrate()].
#'
#' There is no likelihood here, so no log likelihood is reported and the only
#' covariance available is the sandwich \eqn{\Omega^{-1} S \Omega^{-1} / n}.
#'
#' Two cautions when comparing with [panelrate()]. The models are not
#' reparametrisations of each other when covariates vary over time, so their
#' coefficients answer different questions: \eqn{\beta} here acts on the
#' cumulative mean, and in the rate model on the instantaneous rate. And nothing
#' constrains the fitted \eqn{\hat\mu} to increase; when covariate values
#' fluctuate it may well not, which is the difficulty with the means model that
#' motivates the rate model.
#'
#' @return An object of class `"panelmean"`, inheriting from `"pcdfit"`, with
#'   components including `coefficients`, `baseline` (a data frame of the
#'   estimated mean function), `Omega` and `S`, and `vcov`.
#'
#' @references
#' Hu, X. J., Sun, J. and Wei, L. J. (2003). Regression parameter estimation
#' from panel counts. *Scandinavian Journal of Statistics* **30**(1), 25--43.
#' \doi{10.1111/1467-9469.00316}
#'
#' Sun, D., Guo, Y., Li, Y., Tu, W. and Sun, J. (2024).
#' A robust approach for regression analysis of panel count data.
#' *Bernoulli* **30**(4), 3251--3275. \doi{10.3150/23-BEJ1713}
#'
#' @examples
#' set.seed(1)
#' d <- r_panel_count(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
#' mfit <- panelmean(PanelCount(id, tstart, tstop, count) ~ x1 + x2, data = d)
#' summary(mfit)
#'
#' # The rate model on the same data. The coefficients are not comparable
#' # term by term, because the two models act on different quantities.
#' rfit <- panelrate(PanelCount(id, tstart, tstop, count) ~ x1 + x2, data = d)
#' cbind(mean = coef(mfit), rate = coef(rfit))
#'
#' @seealso [panelrate()], [PanelCount()]
#' @export
panelmean <- function(formula, data, subset, na.action,
                      control = panelmean_control(), init = NULL) {
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

  d <- prepare_panel(md$y, md$X, expand = FALSE)
  fit <- mean_fit_cpp(d$exam_X, d$exam_subj, d$exam_grid, d$cN, d$n, d$K,
                      init, control$maxit, control$reltol)
  if (!fit$converged) {
    warning("The estimating equation did not converge in ", control$maxit,
            " iterations (relative change ", format(fit$criterion, digits = 3),
            "). Consider raising `maxit` in `panelmean_control()`.",
            call. = FALSE)
  }

  beta <- drop(fit$beta)
  names(beta) <- colnames(md$X)
  Omega <- fit$Omega
  S <- fit$S
  dimnames(Omega) <- dimnames(S) <- list(names(beta), names(beta))

  vcovs <- if (p > 0) list(robust = robust_vcov(Omega, S, d$n)) else
    list(robust = Omega)
  if (!is.null(vcovs$robust)) dimnames(vcovs$robust) <- dimnames(S)

  structure(
    list(
      coefficients = beta,
      vcov = vcovs,
      Omega = Omega,
      S = S,
      scores = fit$scores,
      baseline = data.frame(time = d$times, mean = drop(fit$mu)),
      n = d$n,
      nexam = d$nexam,
      nevent = sum(d$dN),
      ngrid = d$K,
      iterations = fit$iterations,
      converged = fit$converged,
      criterion = fit$criterion,
      control = control,
      call = cl,
      terms = md$terms,
      xlevels = md$xlevels,
      contrasts = md$contrasts,
      na.action = md$na.action
    ),
    class = c("panelmean", "pcdfit")
  )
}
