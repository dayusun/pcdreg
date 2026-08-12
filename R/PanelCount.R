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
#' There are two forms, distinguished by the number of arguments.
#'
#' `PanelCount(id, time, count)` is for data with one row per examination, which
#' is what time-invariant covariates need. Interval starts are filled in as the
#' previous examination time of the same subject, with the first interval
#' starting at zero.
#'
#' `PanelCount(id, tstart, tstop, count)` is the counting process form, needed
#' when covariates vary over time. Each row gives an interval `(tstart, tstop]`
#' over which the covariates in that row are constant. A row whose `count` is a
#' number records an examination at `tstop`, with `count` events since the
#' previous examination; a row whose `count` is `NA` records only that the
#' covariates changed at `tstop`. Within a subject the intervals must be
#' contiguous and ordered, exactly as for a time-varying [survival::coxph()] fit.
#'
#' A subject's follow-up ends at their last examination. Intervals beyond it
#' carry no information about the model and are dropped, with a message.
#'
#' @return An object of class `"PanelCount"`: a numeric matrix with columns
#'   `id`, `tstart`, `tstop`, `count` and `exam`, carrying the original
#'   identifier labels in the `"labels"` attribute. Missing counts are stored as
#'   `count = 0` with `exam = 0`, so that the response contains no `NA` and
#'   `na.action` handling applies to the covariates alone.
#'
#' @examples
#' # One row per examination.
#' PanelCount(c(1, 1, 2), c(0.5, 1.2, 0.8), c(2, 0, 3))
#'
#' # Counting process form: subject 1's covariates change at t = 0.9, between
#' # examinations at t = 0.5 and t = 1.2.
#' PanelCount(
#'   id     = c(1, 1, 1, 2),
#'   time   = c(0.0, 0.5, 0.9, 0.0),
#'   time2  = c(0.5, 0.9, 1.2, 0.8),
#'   count  = c(2, NA, 0, 3)
#' )
#'
#' @seealso [panelrate()]
#' @export
PanelCount <- function(id, time, time2, count) {
  if (missing(id) || missing(time) || missing(time2)) {
    stop("`PanelCount()` needs at least `id`, `time` and a count.", call. = FALSE)
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
    stop("All arguments to `PanelCount()` must have the same length.",
         call. = FALSE)
  }
  if (n == 0L) stop("`PanelCount()` was given no observations.", call. = FALSE)

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
  class(out) <- "PanelCount"
  out
}

#' @export
print.PanelCount <- function(x, ...) {
  print(format(x), quote = FALSE, ...)
  invisible(x)
}

#' @export
format.PanelCount <- function(x, ...) {
  labels <- attr(x, "labels")
  ids <- if (is.null(labels)) x[, "id"] else labels[x[, "id"]]
  count <- ifelse(x[, "exam"] == 1, format(x[, "count"]), "-")
  paste0(ids, ": (", format(x[, "tstart"]), ", ", format(x[, "tstop"]), "] ",
         count)
}

#' @export
"[.PanelCount" <- function(x, i, j, drop = FALSE) {
  y <- unclass(x)
  labels <- attr(x, "labels")
  attr(y, "labels") <- NULL
  if (missing(j)) {
    out <- y[i, , drop = FALSE]
    attr(out, "labels") <- labels
    class(out) <- "PanelCount"
    out
  } else {
    y[i, j, drop = drop]
  }
}

#' @export
length.PanelCount <- function(x) nrow(unclass(x))

#' @export
is.na.PanelCount <- function(x) {
  as.vector(rowSums(is.na(unclass(x))) > 0)
}

#' @export
as.character.PanelCount <- function(x, ...) format(x, ...)

#' Test for a panel count response
#'
#' @param x An object.
#' @return `TRUE` if `x` was created by [PanelCount()].
#' @examples
#' is.PanelCount(PanelCount(c(1, 2), c(1, 1), c(0, 2)))
#' @export
is.PanelCount <- function(x) inherits(x, "PanelCount")
