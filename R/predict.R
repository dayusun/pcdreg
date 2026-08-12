# The linear predictor for each row of newdata.
pcd_linear_predictor <- function(object, newdata) {
  mt <- stats::delete.response(object$terms)
  mf <- stats::model.frame(mt, newdata, na.action = stats::na.pass,
                           xlev = object$xlevels)
  X <- stats::model.matrix(mt, mf, contrasts.arg = object$contrasts)
  intercept <- match("(Intercept)", colnames(X), 0L)
  if (intercept > 0L) X <- X[, -intercept, drop = FALSE]
  drop(X %*% coef(object))
}

# Anything that follows a covariate trajectory needs the response columns as
# well as the covariates.  model.frame() would fail on its own, but with a bare
# "object 'tstart' not found" that says nothing about why the column was wanted,
# so the columns are checked here first.
pcd_require_response <- function(object, data, arg = "newdata") {
  vars <- all.vars(attr(object$terms, "variables")[[2L]])
  absent <- setdiff(vars, names(data))
  if (length(absent)) {
    cli::cli_abort(c(
      "{.arg {arg}} must contain the {.fn pcd} response.",
      "x" = "Missing column{?s}: {.field {absent}}.",
      "i" = "The predicted mean follows each subject's covariate trajectory,
             so the examination times are needed as well as the covariates."
    ), class = "pcdreg_error_no_response", call = rlang::caller_env())
  }
  invisible(data)
}

# Follow each subject's covariate trajectory across the fitted grid times that
# fall within their follow-up, and hand the per-time linear predictors to
# `combine`, which turns them into the predicted mean the model implies.
pcd_trajectory <- function(object, newdata, combine) {
  pcd_require_response(object, newdata)
  mf <- stats::model.frame(object$terms, newdata, na.action = stats::na.pass,
                           xlev = object$xlevels)
  y <- stats::model.response(mf)
  if (!is.pcd(y)) {
    cli::cli_abort(
      "{.arg newdata} must contain the {.fn pcd} variables for
       {.code type = \"mean\"}.",
      class = "pcdreg_error_no_response"
    )
  }
  X <- stats::model.matrix(object$terms, mf, contrasts.arg = object$contrasts)
  intercept <- match("(Intercept)", colnames(X), 0L)
  if (intercept > 0L) X <- X[, -intercept, drop = FALSE]

  times <- object$baseline$time
  ord <- order(y[, "id"], y[, "tstop"])
  idx <- y[ord, "id"]
  tstop <- y[ord, "tstop"]
  lp <- drop(X[ord, , drop = FALSE] %*% coef(object))
  labels <- attr(y, "labels")

  out <- lapply(split(seq_along(idx), idx), function(rows) {
    upto <- findInterval(max(tstop[rows]), times)
    if (upto == 0L) return(NULL)
    k <- seq_len(upto)
    at <- rows[findInterval(times[k], tstop[rows], left.open = TRUE) + 1L]
    data.frame(id = idx[rows[1L]], time = times[k],
               mean = combine(lp[at], k))
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  if (!is.null(labels)) out$id <- labels[out$id]
  tibble::as_tibble(out)
}

pcd_newdata <- function(object, newdata) {
  if (!missing(newdata) && !is.null(newdata)) return(newdata)
  out <- eval(object$call$data, environment(object$terms))
  if (is.null(out)) {
    cli::cli_abort(
      "{.arg newdata} is required because the original data is no longer
       available.",
      class = "pcdreg_error_no_data"
    )
  }
  out
}

#' Predictions from a panel count model
#'
#' @param object A fitted [pcdreg()] model.
#' @param newdata Data to predict for. Defaults to the data the model was fitted
#'   to. For `type = "mean"` it must also contain the variables in the
#'   [pcd()] response, since the prediction follows each subject's
#'   covariate trajectory; the counts themselves are ignored.
#' @param type `"lp"` returns the linear predictor \eqn{\beta' X(t)} for each
#'   row. `"mean"` returns the predicted mean number of events at every fitted
#'   examination time within each subject's follow-up:
#'   \eqn{\int_0^t \exp(\beta' X(s)) \, d\hat\Lambda(s)} for [pcdreg()], and
#'   \eqn{\hat\mu(t) \exp(\beta' X(t))} for [pcdreg()].
#' @param ... Ignored.
#'
#' @return For `type = "lp"`, a numeric vector with one element per row of
#'   `newdata`. For `type = "mean"`, a [tibble][tibble::tibble] with columns
#'   `id`, `time` and `mean`.
#'
#' @details
#' The rate model's predicted mean is non-decreasing in \eqn{t} by construction,
#' because it integrates a positive quantity. The means model's is not: it is a
#' fitted \eqn{\hat\mu(t)} rescaled by a covariate value that may itself move up
#' and down, so a predicted mean count that falls is possible and is the
#' practical drawback of that model with time-varying covariates.
#'
#' @examples
#' set.seed(1)
#' d <- sim_pcd(60, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
#' fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
#' head(predict(fit))
#' head(predict(fit, type = "mean"))
#' @export
predict.pcdfit <- function(object, newdata, type = c("lp", "mean"), ...) {
  type <- match.arg(type)
  newdata <- pcd_newdata(object, if (missing(newdata)) NULL else newdata)
  if (type == "lp") return(pcd_linear_predictor(object, newdata))
  if (pcd_is_mean(object)) {
    mu <- object$baseline$mean
    pcd_trajectory(object, newdata, function(lp, k) mu[k] * exp(lp))
  } else {
    jump <- object$baseline$jump
    pcd_trajectory(object, newdata, function(lp, k) cumsum(exp(lp) * jump[k]))
  }
}
