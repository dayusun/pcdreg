#' Tidy a panel count model fit
#'
#' Methods for the [generics::tidy()], [generics::glance()] and
#' [generics::augment()] generics, so that a [pcdreg()] fit can be handled the
#' same way as any other model in a tidy workflow.
#'
#' @param x A fitted [pcdreg()] model.
#' @param conf.int Whether to add `conf.low` and `conf.high` columns.
#' @param conf.level Width of the interval, if one is requested.
#' @param exponentiate Whether to report \eqn{e^{\beta}} rather than
#'   \eqn{\beta}. For the rate model this is the multiplicative effect of a
#'   one unit change in the covariate on the rate, read the same way as a
#'   hazard ratio; for the means model it is the effect on the mean. Standard
#'   errors are left on the log scale, as elsewhere in the tidy ecosystem.
#' @param type Which covariance estimator to base the standard errors on. See
#'   [vcov.pcdfit()]; the default is the robust sandwich.
#' @param data,newdata Data to augment. Defaults to the data the model was
#'   fitted to. It must contain the [pcd()] response, since the fitted mean
#'   follows each subject's covariate trajectory.
#' @param ... Ignored.
#'
#' @return
#' `tidy()` gives one row per coefficient, with columns `term`, `estimate`,
#' `std.error`, `statistic` and `p.value`.
#'
#' `glance()` gives a one row summary of the fit. `logLik` is `NA` for the
#' means model, which is fitted by an estimating equation and so has none.
#'
#' `augment()` returns `data` with four columns added:
#' \describe{
#'   \item{`.linear.predictor`}{\eqn{\beta' X(t)} on each row.}
#'   \item{`.fitted`}{the predicted mean number of events by the end of the
#'     row's interval.}
#'   \item{`.observed`}{the observed cumulative count at that time, on
#'     examination rows, and `NA` on rows that only record a covariate change.}
#'   \item{`.resid`}{`.observed - .fitted`, the difference between the counts
#'     seen so far and the mean the model implies.}
#' }
#'
#' @examples
#' set.seed(1)
#' d <- sim_pcd(60, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
#' fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
#'
#' tidy(fit)
#' tidy(fit, conf.int = TRUE, exponentiate = TRUE)
#' glance(fit)
#' augment(fit)
#'
#' @name pcdreg-tidiers
NULL

#' @rdname pcdreg-tidiers
#' @export
tidy.pcdfit <- function(x, conf.int = FALSE, conf.level = 0.95,
                        exponentiate = FALSE,
                        type = c("robust", "information", "profile"), ...) {
  type <- match.arg(type)
  est <- stats::coef(x)
  if (!length(est)) {
    return(tibble::tibble(term = character(), estimate = numeric(),
                          std.error = numeric(), statistic = numeric(),
                          p.value = numeric()))
  }

  se <- sqrt(diag(stats::vcov(x, type)))
  statistic <- est / se
  out <- tibble::tibble(
    term = names(est),
    estimate = unname(est),
    std.error = unname(se),
    statistic = unname(statistic),
    p.value = 2 * stats::pnorm(-abs(unname(statistic)))
  )

  if (conf.int) {
    if (!rlang::is_double(conf.level, n = 1) || conf.level <= 0 ||
        conf.level >= 1) {
      cli::cli_abort("{.arg conf.level} must be a single number in (0, 1).",
                     class = "pcdreg_error_conf_level")
    }
    z <- stats::qnorm(1 - (1 - conf.level) / 2)
    out$conf.low <- out$estimate - z * out$std.error
    out$conf.high <- out$estimate + z * out$std.error
  }

  if (exponentiate) {
    out$estimate <- exp(out$estimate)
    if (conf.int) {
      out$conf.low <- exp(out$conf.low)
      out$conf.high <- exp(out$conf.high)
    }
  }
  out
}

#' @rdname pcdreg-tidiers
#' @export
glance.pcdfit <- function(x, ...) {
  tibble::tibble(
    model = pcd_model_label(x),
    n = x$n,
    nexam = x$nexam,
    nevent = x$nevent,
    ngrid = x$ngrid,
    logLik = if (pcd_is_mean(x)) NA_real_ else as.numeric(x$loglik),
    iterations = x$iterations,
    converged = x$converged
  )
}

#' @rdname pcdreg-tidiers
#' @export
augment.pcdfit <- function(x, data = NULL, newdata = NULL, ...) {
  at <- if (!is.null(newdata)) newdata else data
  at <- pcd_newdata(x, at)
  pcd_require_response(x, at, arg = "data")

  mf <- stats::model.frame(x$terms, at, na.action = stats::na.pass,
                           xlev = x$xlevels)
  y <- stats::model.response(mf)

  labels <- attr(y, "labels")
  ids <- if (is.null(labels)) y[, "id"] else labels[y[, "id"]]
  tstop <- y[, "tstop"]

  traj <- stats::predict(x, at, type = "mean")
  # The trajectory is a step function on the pooled grid, so the fitted mean at
  # a row's tstop is its value at the last grid time not after it.  A tstop
  # before the first grid time has fitted mean zero: both models integrate up
  # from the origin.
  fitted <- rep(NA_real_, length(tstop))
  by_subject <- split(seq_along(tstop), ids)
  for (s in names(by_subject)) {
    rows <- by_subject[[s]]
    tr <- traj[traj$id == s, , drop = FALSE]
    if (!nrow(tr)) next
    k <- findInterval(tstop[rows], tr$time)
    fitted[rows] <- ifelse(k > 0L, tr$mean[pmax(k, 1L)], 0)
  }

  # Counts are stored as zero on covariate-change rows, so a running sum over
  # all of a subject's rows is the observed cumulative count at each tstop.
  ord <- order(y[, "id"], tstop)
  observed <- rep(NA_real_, length(tstop))
  observed[ord] <- stats::ave(y[ord, "count"], y[ord, "id"], FUN = cumsum)
  observed[y[, "exam"] != 1] <- NA_real_

  out <- tibble::as_tibble(at)
  out$.linear.predictor <- pcd_linear_predictor(x, at)
  out$.fitted <- fitted
  out$.observed <- observed
  out$.resid <- observed - fitted
  out
}

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom generics augment
#' @export
generics::augment
