#' Predictions from a proportional rate model
#'
#' @param object A fitted [panelrate()] model.
#' @param newdata Data to predict for. Defaults to the data the model was fitted
#'   to. For `type = "mean"` it must also contain the variables in the
#'   [PanelCount()] response, since the prediction follows each subject's
#'   covariate trajectory; the counts themselves are ignored.
#' @param type `"lp"` returns the linear predictor \eqn{\beta' X(t)} for each
#'   row. `"mean"` returns the predicted mean number of events
#'   \eqn{E[N(t) \mid X] = \int_0^t \exp(\beta' X(s)) d\hat\Lambda(s)} at every
#'   fitted examination time within each subject's follow-up.
#' @param ... Ignored.
#'
#' @return For `type = "lp"`, a numeric vector with one element per row of
#'   `newdata`. For `type = "mean"`, a data frame with columns `id`, `time` and
#'   `mean`.
#'
#' @details
#' Because the rate model constrains only \eqn{\exp(\beta' X(t)) d\Lambda(t)} to
#' be positive, the predicted mean is automatically non-decreasing in \eqn{t}
#' even when covariates fluctuate, which is the practical advantage of the rate
#' model over the proportional means model.
#'
#' @examples
#' set.seed(1)
#' d <- r_panel_count(60, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
#' fit <- panelrate(PanelCount(id, tstart, tstop, count) ~ x1 + x2, data = d)
#' head(predict(fit))
#' head(predict(fit, type = "mean"))
#' @export
predict.panelrate <- function(object, newdata, type = c("lp", "mean"), ...) {
  type <- match.arg(type)
  beta <- coef(object)

  if (missing(newdata) || is.null(newdata)) {
    newdata <- eval(object$call$data, environment(object$terms))
    if (is.null(newdata)) {
      stop("`newdata` is required because the original data is no longer ",
           "available.", call. = FALSE)
    }
  }

  if (type == "lp") {
    mt <- stats::delete.response(object$terms)
    mf <- stats::model.frame(mt, newdata, na.action = stats::na.pass,
                             xlev = object$xlevels)
    X <- stats::model.matrix(mt, mf, contrasts.arg = object$contrasts)
    intercept <- match("(Intercept)", colnames(X), 0L)
    if (intercept > 0L) X <- X[, -intercept, drop = FALSE]
    return(drop(X %*% beta))
  }

  mf <- stats::model.frame(object$terms, newdata, na.action = stats::na.pass,
                           xlev = object$xlevels)
  y <- stats::model.response(mf)
  if (!is.PanelCount(y)) {
    stop("`newdata` must contain the `PanelCount()` variables for ",
         "`type = \"mean\"`.", call. = FALSE)
  }
  X <- stats::model.matrix(object$terms, mf, contrasts.arg = object$contrasts)
  intercept <- match("(Intercept)", colnames(X), 0L)
  if (intercept > 0L) X <- X[, -intercept, drop = FALSE]

  times <- object$baseline$time
  jump <- object$baseline$jump
  ord <- order(y[, "id"], y[, "tstop"])
  idx <- y[ord, "id"]
  tstop <- y[ord, "tstop"]
  lp <- drop(X[ord, , drop = FALSE] %*% beta)
  labels <- attr(y, "labels")

  rows_by_subject <- split(seq_along(idx), idx)
  out <- lapply(rows_by_subject, function(rows) {
    upto <- findInterval(max(tstop[rows]), times)
    if (upto == 0L) return(NULL)
    k <- seq_len(upto)
    at <- rows[findInterval(times[k], tstop[rows], left.open = TRUE) + 1L]
    data.frame(id = idx[rows[1L]], time = times[k],
               mean = cumsum(exp(lp[at]) * jump[k]))
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  if (!is.null(labels)) out$id <- labels[out$id]
  out
}
