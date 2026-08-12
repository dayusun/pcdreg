#' Describe panel count observations
#'
#' Creates the response object used on the left hand side of a [panelrate()]
#' formula, in the same spirit as [survival::Surv()].
#'
#' @param id Subject identifier. May be numeric, character or a factor; rows
#'   sharing an identifier are treated as one subject.
#' @param time Examination time in the three argument form, or the start of the
#'   interval in the four argument form.
#' @param time2 The end of the interval in the four argument form; the event
#'   count in the three argument form.
#' @param count Number of events observed since the previous examination. Use
#'   `NA` on rows that only record a change in a covariate rather than an
#'   examination. Omitted in the three argument form.
#'
#' @details
#' There are two forms, distinguished by the number of arguments. Which one you
#' need depends on a single question: does any covariate change value during
#' follow-up?
#'
#' `pcd(id, time, count)` is for data with one row per examination, which
#' is all that time-invariant covariates need. Interval starts are filled in as
#' the previous examination time of the same subject, with the first interval
#' starting at zero.
#'
#' `pcd(id, tstart, tstop, count)` is the counting process form, needed as
#' soon as a covariate moves. Each row gives an interval `(tstart, tstop]` over
#' which the covariates in that row are constant, exactly as for a time-varying
#' [survival::coxph()] fit. The `count` column distinguishes the two kinds of
#' row:
#'
#' * a number means the row ends at an **examination**, and that many events
#'   occurred since the previous examination;
#' * `NA` means the row ends at a **covariate change** rather than an
#'   examination.
#'
#' There is no separate event indicator; the presence or absence of a count is
#' the indicator.
#'
#' @section Counts are increments, not totals:
#' `count` must be the number of events since the previous examination,
#' \eqn{\Delta N_{ij}}, not the running total \eqn{N_i(T_{ij})}. Data recorded
#' cumulatively must be differenced within subject first. This is not something
#' the package can detect for you, and getting it wrong inflates the event count
#' and biases the coefficients rather than raising an error. A quick check after
#' fitting is that `fit$nevent` matches the number of events you believe you
#' have. See `vignette("data-preparation")` for the conversion.
#'
#' @section What the data must satisfy:
#' Checked before any arithmetic, naming the offending subject when a check
#' fails:
#'
#' * follow-up starts at time zero, since the models integrate from there;
#' * intervals within a subject are contiguous, with no gaps or overlaps;
#' * no subject has two rows ending at the same time;
#' * every subject has at least one examination;
#' * counts are non-negative whole numbers.
#'
#' Row order is irrelevant: the data are sorted internally. Identifiers may be
#' numeric, character or factor, and need not be consecutive.
#'
#' A subject's follow-up ends at their last examination, because the likelihood
#' runs only over \eqn{j = 1, \ldots, J_i}. Rows lying entirely beyond it carry
#' no information and are dropped, with a message saying how many. A censoring
#' row appended after the final visit is the usual reason for seeing it.
#'
#' @return An object of class `"pcd"`: a numeric matrix with columns
#'   `id`, `tstart`, `tstop`, `count` and `exam`, carrying the original
#'   identifier labels in the `"labels"` attribute. Missing counts are stored as
#'   `count = 0` with `exam = 0`, so that the response contains no `NA` and
#'   `na.action` handling applies to the covariates alone.
#'
#' @examples
#' # One row per examination: subject 1 seen at 0.5 and 1.2 with 2 then 0
#' # events, subject 2 seen once at 0.8 with 3.
#' pcd(c(1, 1, 2), c(0.5, 1.2, 0.8), c(2, 0, 3))
#'
#' # Counting process form. Subject 1's covariates change at t = 0.9, part way
#' # through the examination interval (0.5, 1.2], so that interval needs two
#' # rows: the first ends at the change and carries NA, the second ends at the
#' # examination and carries its count.
#' pcd(
#'   id     = c(1, 1, 1, 2),
#'   time   = c(0.0, 0.5, 0.9, 0.0),
#'   time2  = c(0.5, 0.9, 1.2, 0.8),
#'   count  = c(2, NA, 0, 3)
#' )
#'
#' # Used on the left hand side of a formula, which is how you will normally
#' # meet it.
#' set.seed(1)
#' d <- r_panel_count(40)
#' fit <- panelrate(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
#' c(subjects = fit$n, examinations = fit$nexam, events = fit$nevent)
#'
#' @seealso [panelrate()], [panelmean()], and
#'   `vignette("data-preparation")` for converting data into this shape.
#' @export
pcd <- function(id, time, time2, count) {
  if (missing(id) || missing(time) || missing(time2)) {
    stop("`pcd()` needs at least `id`, `time` and a count.", call. = FALSE)
  }
  if (missing(count)) {
    tstop <- time
    count <- time2
    tstart <- NULL
  } else {
    tstart <- time
    tstop <- time2
  }

  n <- length(id)
  if (length(tstop) != n || length(count) != n ||
      (!is.null(tstart) && length(tstart) != n)) {
    stop("All arguments to `pcd()` must have the same length.",
         call. = FALSE)
  }
  if (n == 0L) stop("`pcd()` was given no observations.", call. = FALSE)

  labels <- NULL
  if (is.factor(id)) {
    labels <- levels(id)
    idx <- as.integer(id)
  } else if (is.character(id)) {
    f <- factor(id)
    labels <- levels(f)
    idx <- as.integer(f)
  } else {
    if (anyNA(id)) stop("`id` must not be missing.", call. = FALSE)
    f <- factor(id, levels = sort(unique(id)))
    labels <- as.character(levels(f))
    idx <- as.integer(f)
  }
  if (anyNA(idx)) stop("`id` must not be missing.", call. = FALSE)

  tstop <- as.numeric(tstop)
  if (anyNA(tstop) || any(!is.finite(tstop))) {
    stop("Examination times must be present and finite.", call. = FALSE)
  }

  exam <- !is.na(count)
  count <- as.numeric(count)
  count[!exam] <- 0
  if (any(count < 0) || any(count != trunc(count))) {
    stop("`count` must be a non-negative whole number at examination times.",
         call. = FALSE)
  }
  if (!any(exam)) {
    stop("No examinations found: `count` is `NA` on every row.", call. = FALSE)
  }

  if (is.null(tstart)) {
    # One row per examination: the previous examination of the same subject
    # opens the interval, and follow-up starts at zero.
    ord <- order(idx, tstop)
    tstart <- numeric(n)
    prev <- c(NA_real_, tstop[ord][-n])
    same <- c(FALSE, idx[ord][-1L] == idx[ord][-n])
    tstart[ord] <- ifelse(same, prev, 0)
  } else {
    tstart <- as.numeric(tstart)
    if (anyNA(tstart) || any(!is.finite(tstart))) {
      stop("`tstart` must be present and finite.", call. = FALSE)
    }
  }

  if (any(tstart >= tstop)) {
    stop("Every interval must have `tstart` strictly less than `tstop`; ",
         "row ", which(tstart >= tstop)[1L], " does not.", call. = FALSE)
  }

  out <- cbind(id = idx, tstart = tstart, tstop = tstop, count = count,
               exam = as.numeric(exam))
  attr(out, "labels") <- labels
  class(out) <- "pcd"
  out
}

#' @export
print.pcd <- function(x, ...) {
  print(format(x), quote = FALSE, ...)
  invisible(x)
}

#' @export
format.pcd <- function(x, ...) {
  labels <- attr(x, "labels")
  ids <- if (is.null(labels)) x[, "id"] else labels[x[, "id"]]
  count <- ifelse(x[, "exam"] == 1, format(x[, "count"]), "-")
  paste0(ids, ": (", format(x[, "tstart"]), ", ", format(x[, "tstop"]), "] ",
         count)
}

#' @export
"[.pcd" <- function(x, i, j, drop = FALSE) {
  y <- unclass(x)
  labels <- attr(x, "labels")
  attr(y, "labels") <- NULL
  if (missing(j)) {
    out <- y[i, , drop = FALSE]
    attr(out, "labels") <- labels
    class(out) <- "pcd"
    out
  } else {
    y[i, j, drop = drop]
  }
}

#' @export
length.pcd <- function(x) nrow(unclass(x))

#' @export
is.na.pcd <- function(x) {
  as.vector(rowSums(is.na(unclass(x))) > 0)
}

#' @export
as.character.pcd <- function(x, ...) format(x, ...)

#' Test for a panel count response
#'
#' @param x An object.
#' @return `TRUE` if `x` was created by [pcd()].
#' @examples
#' is.pcd(pcd(c(1, 2), c(1, 1), c(0, 2)))
#' @export
is.pcd <- function(x) inherits(x, "pcd")
