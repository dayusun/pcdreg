#' Simulate panel count data with a time-varying covariate
#'
#' Generates data from the design used in the simulation study of Sun et al.
#' (2024), which is a convenient source of examples and test cases.
#'
#' @param n Number of subjects.
#' @param beta Length two coefficient vector for `x1` and `x2`.
#' @param lambda Baseline rate function, vectorised over `t`.
#' @param tau Maximum follow-up time.
#' @param frailty Variance of a gamma frailty with mean one that multiplies each
#'   subject's rate. Zero, the default, gives a genuine Poisson process; a
#'   positive value produces overdispersed counts that violate the Poisson
#'   assumption while leaving the rate model intact.
#' @param exam_mean Mean of the zero-truncated Poisson distribution from which
#'   the number of examinations per subject is drawn.
#' @param lambda_max Optional upper bound for `lambda` on `[0, tau]`, used by the
#'   thinning algorithm. The default takes the maximum over a fine grid, which
#'   is reliable unless `lambda` oscillates very rapidly.
#' @param digits Number of decimal places to round examination times to, or
#'   `NA` to leave them unrounded. Rounding follows the design of the paper and
#'   of Wellner and Zhang (2007): examinations then tie across subjects, which
#'   is realistic for scheduled visits and keeps the pooled grid, and hence the
#'   cost of fitting, from growing quadratically in the sample size.
#'
#' @return A data frame in counting process form with columns `id`, `tstart`,
#'   `tstop`, `count`, `x1` and `x2`. Rows are intervals over which the
#'   covariates are constant; `count` is the number of events since the previous
#'   examination and is `NA` on rows that only record a covariate change.
#'
#' @details
#' The time-varying covariate is a single step,
#' \eqn{x_1(t) = X_{11} I(t \le U_1) + X_{12} I(t > U_1)}, with
#' \eqn{X_{11}, X_{12} \sim U(0, 1)} and \eqn{U_1 \sim U(0, \tau)}. The
#' time-invariant covariate `x2` is \eqn{U(0, 1)}. Events follow a
#' nonhomogeneous Poisson process with rate
#' \eqn{r \lambda(t) \exp(\beta' X(t))}, where \eqn{r} is one unless `frailty`
#' is positive. Examination times are the order statistics of a uniform sample
#' on \eqn{(0, \bar T)} with \eqn{\bar T \sim U(0.9\tau, \tau)}.
#'
#' @references
#' Sun, D., Guo, Y., Li, Y., Tu, W. and Sun, J. (2024).
#' A robust approach for regression analysis of panel count data.
#' *Bernoulli* **30**(4), 3251--3275. \doi{10.3150/23-BEJ1713}
#'
#' @examples
#' set.seed(42)
#' d <- r_panel_count(5, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
#' head(d)
#'
#' # Overdispersed counts: the rate model still holds, the Poisson one does not.
#' od <- r_panel_count(5, frailty = 1)
#' @export
r_panel_count <- function(n, beta = c(1, -1),
                          lambda = function(t) 8 / (1 + t), tau = 2,
                          frailty = 0, exam_mean = 4, lambda_max = NULL,
                          digits = 2) {
  stopifnot(length(n) == 1L, n >= 1, length(beta) == 2L, is.function(lambda),
            length(tau) == 1L, tau > 0, length(frailty) == 1L, frailty >= 0,
            length(exam_mean) == 1L, exam_mean > 0, length(digits) == 1L)
  if (is.null(lambda_max)) {
    lambda_max <- max(lambda(seq(0, tau, length.out = 10001L)))
  }
  if (!is.finite(lambda_max) || lambda_max <= 0) {
    stop("`lambda` must be positive and finite on [0, tau].", call. = FALSE)
  }

  subjects <- lapply(seq_len(n), function(i) {
    one_subject(i, beta, lambda, tau, frailty, exam_mean, lambda_max, digits)
  })
  out <- do.call(rbind, subjects)
  rownames(out) <- NULL
  out
}

# One subject's worth of rows.
one_subject <- function(i, beta, lambda, tau, frailty, exam_mean, lambda_max,
                        digits) {
  # Examination times.  Rounding can collapse or zero out draws, so keep
  # drawing until at least one positive examination time survives.
  end <- 0
  while (end <= 0) {
    u <- stats::runif(1, stats::ppois(0, exam_mean), 1)
    nexam <- max(1L, stats::qpois(u, exam_mean))
    tbar <- stats::runif(1, 0.9 * tau, tau)
    exams <- stats::runif(nexam, 0, tbar)
    if (!is.na(digits)) exams <- round(exams, digits)
    exams <- sort(unique(exams[exams > 0]))
    end <- if (length(exams)) exams[length(exams)] else 0
  }

  # Covariates: one step for x1, constant x2.
  x11 <- stats::runif(1)
  x12 <- stats::runif(1)
  x2 <- stats::runif(1)
  u1 <- stats::runif(1, 0, tau)

  breaks <- sort(unique(c(exams, if (u1 > 0 && u1 < end) u1)))
  tstart <- c(0, breaks[-length(breaks)])
  tstop <- breaks
  x1 <- ifelse(tstop <= u1, x11, x12)

  # Event times by thinning, with the rate bounded by lambda_max times the
  # largest relative rate over the trajectory.
  r <- if (frailty > 0) stats::rgamma(1, shape = 1 / frailty, scale = frailty)
       else 1
  x1_at <- function(t) ifelse(t <= u1, x11, x12)
  rate <- function(t) r * lambda(t) * exp(beta[1] * x1_at(t) + beta[2] * x2)
  bound <- r * lambda_max * exp(max(beta[1] * c(x11, x12)) + beta[2] * x2)
  bound <- max(bound, .Machine$double.eps)

  ncand <- stats::rpois(1, bound * end)
  events <- numeric(0)
  if (ncand > 0L) {
    cand <- stats::runif(ncand, 0, end)
    events <- cand[stats::runif(ncand) <= rate(cand) / bound]
  }

  count <- rep(NA_real_, length(tstop))
  is_exam <- tstop %in% exams
  count[is_exam] <- diff(c(0, findInterval(exams, sort(events))))

  data.frame(id = i, tstart = tstart, tstop = tstop, count = count,
             x1 = x1, x2 = x2)
}
