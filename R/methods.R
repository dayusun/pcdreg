# Both fitters return an object inheriting from "pcdfit", so the extractors
# below are written once.  Where the two models genuinely differ -- the baseline
# is a cumulative rate for one and a mean for the other, and only the rate model
# has a likelihood -- the difference is confined to these three helpers and to
# the methods that really must be separate.

pcd_is_mean <- function(object) identical(object$model, "mean")

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
#' @param object A fitted [pcdreg()] model.
#' @param type Which estimator to return. `"robust"` is the sandwich estimator
#'   that does not rely on the Poisson assumption and is available for both
#'   models. `"information"` and `"profile"` apply to [pcdreg()] only:
#'   the first is the efficient information estimator, valid under the Poisson
#'   assumption, and the second is the profile likelihood estimator, available
#'   only if the model was fitted with `profile = TRUE`.
#' @param ... Ignored.
#'
#' @return A covariance matrix.
#'
#' @examples
#' set.seed(1)
#' d <- sim_pcd(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
#' fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
#' sqrt(diag(vcov(fit)))
#' sqrt(diag(vcov(fit, "information")))
#' @export
vcov.pcdfit <- function(object, type = c("robust", "information", "profile"),
                        ...) {
  type <- match.arg(type)
  out <- object$vcov[[type]]
  if (is.null(out)) {
    if (pcd_is_mean(object)) {
      cli::cli_abort(
        "The means model is fitted by an estimating equation rather than a
         likelihood, so only {.code type = \"robust\"} is available.",
        class = "pcdreg_error_vcov_type"
      )
    }
    if (type == "profile") {
      cli::cli_abort(c(
        "The profile likelihood covariance was not computed.",
        "i" = "Refit with {.code pcdreg(..., profile = TRUE)}."
      ), class = "pcdreg_error_no_profile")
    }
    cli::cli_abort("The {type} covariance is unavailable for this fit.",
                   class = "pcdreg_error_vcov_type")
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
logLik.pcdfit <- function(object, ...) {
  if (pcd_is_mean(object)) {
    cli::cli_abort(c(
      "The means model is fitted by an estimating equation, so it has no log
       likelihood.",
      "i" = "Refit with {.code model = \"rate\"} if you need one."
    ), class = "pcdreg_error_no_loglik")
  }
  # The baseline contributes one free jump per grid time.
  structure(object$loglik,
            df = length(coef(object)) + object$ngrid,
            nobs = object$n, class = "logLik")
}

#' Estimated baseline function
#'
#' @param object A fitted [pcdreg()] model.
#' @param ... Ignored.
#' @return A [tibble][tibble::tibble] with one row per examination time on the
#'   pooled grid.
#'   For [pcdreg()] it gives the jump sizes and the cumulative baseline rate
#'   \eqn{\Lambda(t)}; for [pcdreg()] it gives the baseline mean
#'   \eqn{\mu(t)}, which is not constrained to increase.
#' @examples
#' set.seed(1)
#' d <- sim_pcd(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
#' head(baseline(pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, d)))
#' head(baseline(pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, d, model = "mean")))
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
#' @param object A fitted [pcdreg()] model.
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
#' d <- sim_pcd(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
#' fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
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
#' @param object,x A fitted [pcdreg()] model.
#' @param ... Ignored.
#'
#' @return A [ggplot2::ggplot()] object: a step function of the cumulative
#'   baseline rate for the rate model, or of the baseline mean for the means
#'   model. `autoplot()` returns it; `plot()` draws it and returns it
#'   invisibly, so it can be used for its side effect like any other `plot()`
#'   method.
#'
#' @details
#' The rate model's \eqn{\hat\Lambda} increases by construction. The means
#' model's \eqn{\hat\mu} need not, and with covariates that fluctuate it
#' generally does not; it is drawn as estimated rather than forced upwards,
#' because that behaviour is the point of the comparison between the two models.
#'
#' @examples
#' set.seed(1)
#' d <- sim_pcd(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
#' autoplot(pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, d))
#' @export
autoplot.pcdfit <- function(object, ...) {
  b <- object$baseline
  is_mean <- pcd_is_mean(object)
  column <- pcd_baseline_column(object)
  df <- data.frame(time = c(0, b$time), value = c(0, b[[column]]))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$value)) +
    ggplot2::geom_step(colour = "#2a78d6", linewidth = 0.6) +
    ggplot2::labs(x = "Time",
                  y = if (is_mean) expression(hat(mu)(t)) else
                    expression(hat(Lambda)(t))) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
    pcd_theme()
}

#' @rdname autoplot.pcdfit
#' @export
plot.pcdfit <- function(x, ...) {
  p <- autoplot.pcdfit(x, ...)
  print(p)
  invisible(p)
}
