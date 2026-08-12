# Both fitters return an object inheriting from "pcdfit", so the extractors
# below are written once.  Where the two models genuinely differ -- the baseline
# is a cumulative rate for one and a mean for the other, and only the rate model
# has a likelihood -- the difference is confined to these three helpers and to
# the methods that really must be separate.

pcd_is_mean <- function(object) inherits(object, "panelmean")

pcd_model_label <- function(object) {
  if (pcd_is_mean(object)) "Proportional means model for panel count data"
  else "Proportional rate model for panel count data"
}

pcd_baseline_column <- function(object) {
  if (pcd_is_mean(object)) "mean" else "cumrate"
}

#' @export
coef.pcdfit <- function(object, ...) object$coefficients

#' @export
nobs.pcdfit <- function(object, ...) object$n

#' Covariance matrix of the estimated coefficients
#'
#' @param object A fitted [panelrate()] or [panelmean()] model.
#' @param type Which estimator to return. `"robust"` is the sandwich estimator
#'   that does not rely on the Poisson assumption and is available for both
#'   models. `"information"` and `"profile"` apply to [panelrate()] only:
#'   the first is the efficient information estimator, valid under the Poisson
#'   assumption, and the second is the profile likelihood estimator, available
#'   only if the model was fitted with `profile = TRUE`.
#' @param ... Ignored.
#'
#' @return A covariance matrix.
#'
#' @examples
#' set.seed(1)
#' d <- r_panel_count(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
#' fit <- panelrate(PanelCount(id, tstart, tstop, count) ~ x1 + x2, data = d)
#' sqrt(diag(vcov(fit)))
#' sqrt(diag(vcov(fit, "information")))
#' @export
vcov.pcdfit <- function(object, type = c("robust", "information", "profile"),
                        ...) {
  type <- match.arg(type)
  out <- object$vcov[[type]]
  if (is.null(out)) {
    if (pcd_is_mean(object)) {
      stop("The means model is fitted by an estimating equation rather than a ",
           "likelihood, so only `type = \"robust\"` is available.",
           call. = FALSE)
    }
    if (type == "profile") {
      stop("The profile likelihood covariance was not computed. Refit with ",
           "`panelrate(..., profile = TRUE)`.", call. = FALSE)
    }
    stop("The ", type, " covariance is unavailable for this fit.",
         call. = FALSE)
  }
  out
}

#' @export
confint.pcdfit <- function(object, parm, level = 0.95,
                           type = c("robust", "information", "profile"), ...) {
  type <- match.arg(type)
  est <- coef(object)
  se <- sqrt(diag(vcov(object, type)))
  if (missing(parm)) {
    parm <- names(est)
  } else if (is.numeric(parm)) {
    parm <- names(est)[parm]
  }
  a <- (1 - level) / 2
  z <- stats::qnorm(1 - a)
  out <- cbind(est[parm] - z * se[parm], est[parm] + z * se[parm])
  colnames(out) <- paste(format(100 * c(a, 1 - a), trim = TRUE, digits = 3),
                         "%")
  out
}

#' @export
logLik.panelrate <- function(object, ...) {
  # The baseline contributes one free jump per grid time.
  structure(object$loglik,
            df = length(coef(object)) + object$ngrid,
            nobs = object$n, class = "logLik")
}

#' @export
logLik.panelmean <- function(object, ...) {
  stop("The means model is fitted by an estimating equation, so it has no ",
       "log likelihood. Use `panelrate()` if you need one.", call. = FALSE)
}

#' Estimated baseline function
#'
#' @param object A fitted [panelrate()] or [panelmean()] model.
#' @param ... Ignored.
#' @return A data frame with one row per examination time on the pooled grid.
#'   For [panelrate()] it gives the jump sizes and the cumulative baseline rate
#'   \eqn{\Lambda(t)}; for [panelmean()] it gives the baseline mean
#'   \eqn{\mu(t)}, which is not constrained to increase.
#' @examples
#' set.seed(1)
#' d <- r_panel_count(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
#' head(baseline(panelrate(PanelCount(id, tstart, tstop, count) ~ x1 + x2, d)))
#' head(baseline(panelmean(PanelCount(id, tstart, tstop, count) ~ x1 + x2, d)))
#' @export
baseline <- function(object, ...) UseMethod("baseline")

#' @rdname baseline
#' @export
baseline.pcdfit <- function(object, ...) object$baseline

#' @export
print.pcdfit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2L,
                  quote = FALSE)
  } else {
    cat("No coefficients (baseline only).\n")
  }
  cat("\n", x$n, " subjects, ", x$nexam, " examinations, ", x$nevent,
      " events.\n", sep = "")
  if (pcd_is_mean(x)) {
    cat("Converged in ", x$iterations, " Newton iterations.\n", sep = "")
  } else {
    cat("Log likelihood ", format(x$loglik, digits = digits), " in ",
        x$iterations, " EM iterations.\n", sep = "")
  }
  if (!x$converged) cat("Warning: the algorithm did not converge.\n")
  invisible(x)
}

#' Summarise a fitted panel count model
#'
#' @param object A fitted [panelrate()] or [panelmean()] model.
#' @param type Which covariance estimator to base the standard errors on. See
#'   [vcov.pcdfit()].
#' @param ... Ignored.
#'
#' @return An object of class `"summary.pcdfit"`, whose `coefficients`
#'   component is the usual four column table of estimates, standard errors,
#'   Wald statistics and two sided p-values.
#'
#' @examples
#' set.seed(1)
#' d <- r_panel_count(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
#' fit <- panelrate(PanelCount(id, tstart, tstop, count) ~ x1 + x2, data = d)
#' summary(fit)
#' summary(fit, "information")
#' @export
summary.pcdfit <- function(object, type = c("robust", "information",
                                            "profile"), ...) {
  type <- match.arg(type)
  est <- coef(object)
  if (length(est)) {
    se <- sqrt(diag(vcov(object, type)))
    z <- est / se
    tab <- cbind(Estimate = est, `Std. Error` = se, `z value` = z,
                 `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  } else {
    tab <- matrix(numeric(0), 0L, 4L,
                  dimnames = list(NULL, c("Estimate", "Std. Error", "z value",
                                          "Pr(>|z|)")))
  }
  structure(list(call = object$call, coefficients = tab, type = type,
                 label = pcd_model_label(object), is_mean = pcd_is_mean(object),
                 n = object$n, nexam = object$nexam, nevent = object$nevent,
                 ngrid = object$ngrid, loglik = object$loglik,
                 iterations = object$iterations, converged = object$converged),
            class = "summary.pcdfit")
}

#' @export
print.summary.pcdfit <- function(x, digits = max(3L, getOption("digits") - 3L),
                                 ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  label <- c(robust = "robust sandwich", information = "efficient information",
             profile = "profile likelihood")[x$type]
  cat(x$label, "\n")
  cat("Standard errors: ", label, "\n\n", sep = "")
  if (nrow(x$coefficients)) {
    stats::printCoefmat(x$coefficients, digits = digits, signif.stars = FALSE)
  } else {
    cat("No coefficients (baseline only).\n")
  }
  cat("\n", x$n, " subjects, ", x$nexam, " examinations, ", x$nevent,
      " events, ", x$ngrid, " distinct examination times.\n", sep = "")
  if (x$is_mean) {
    cat("Converged in ", x$iterations, " Newton iterations.\n", sep = "")
  } else {
    cat("Log likelihood ", format(x$loglik, digits = digits), " in ",
        x$iterations, " EM iterations.\n", sep = "")
  }
  if (!x$converged) cat("Warning: the algorithm did not converge.\n")
  if (x$type != "robust") {
    cat("\nNote: these standard errors assume the counts are Poisson. If the\n",
        "counts are overdispersed they will be too small; compare with\n",
        "`summary(fit, \"robust\")`.\n", sep = "")
  }
  invisible(x)
}

#' Plot the estimated baseline function
#'
#' @param x A fitted [panelrate()] or [panelmean()] model.
#' @param xlab,ylab,type Passed to [graphics::plot()], with defaults suited to a
#'   step function and to whichever baseline the model estimates.
#' @param ... Further graphical parameters.
#' @return `x`, invisibly.
#' @examples
#' set.seed(1)
#' d <- r_panel_count(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
#' plot(panelrate(PanelCount(id, tstart, tstop, count) ~ x1 + x2, d))
#' @export
plot.pcdfit <- function(x, xlab = "Time", ylab = NULL, type = "s", ...) {
  b <- x$baseline
  column <- pcd_baseline_column(x)
  if (is.null(ylab)) {
    ylab <- if (pcd_is_mean(x)) expression(hat(mu)(t)) else
      expression(hat(Lambda)(t))
  }
  graphics::plot(c(0, b$time), c(0, b[[column]]), xlab = xlab, ylab = ylab,
                 type = type, ...)
  invisible(x)
}
