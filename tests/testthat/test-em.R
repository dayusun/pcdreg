make_data <- function(n = 30, seed = 20240101, ...) {
  set.seed(seed)
  r_panel_count(n, beta = c(1, -1), lambda = function(t) 8 / (1 + t), ...)
}

test_that("the compiled EM matches an independent pure R implementation", {
  # Both run the same fixed number of plain EM passes from the same start, so
  # they should agree to machine precision, iteration for iteration.  The
  # tolerance is deliberately far below anything that could hide a difference
  # in the formulas.
  d <- prep(make_data(25))
  iters <- 200L
  got <- pcdreg:::em_fit_cpp(d$X, d$subj, d$grid, d$panel, d$dN,
                                d$panelsubj, d$n, d$K, c(0, 0),
                                rep(1 / d$K, d$K), iters, 1e-15, FALSE)
  want <- ref_em(d, maxit = iters, reltol = 1e-15)

  expect_equal(got$iterations, iters)
  expect_equal(want$iterations, iters)
  expect_equal(as.vector(got$beta), as.vector(want$beta), tolerance = 1e-11)
  expect_equal(as.vector(got$lambda), as.vector(want$lambda), tolerance = 1e-11)
  expect_equal(got$loglik, ref_loglik(d, as.vector(got$beta),
                                      as.vector(got$lambda)))
})

test_that("acceleration changes the path but not the destination", {
  d <- make_data(40)
  fast <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d,
                    control = pcdreg_control(reltol = 1e-10, maxit = 5000))
  slow <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d,
                    control = pcdreg_control(reltol = 1e-10, maxit = 100000,
                                                accelerate = FALSE))
  expect_true(fast$converged)
  expect_equal(coef(fast), coef(slow), tolerance = 1e-5)
  expect_equal(fast$baseline$cumrate, slow$baseline$cumrate, tolerance = 1e-5)
  # The whole point of accelerating.
  expect_lt(fast$passes, slow$iterations)
})

test_that("the solution satisfies the score and fixed point equations", {
  d <- prep(make_data(30))
  fit <- pcdreg:::em_fit_cpp(d$X, d$subj, d$grid, d$panel, d$dN,
                                d$panelsubj, d$n, d$K, c(0, 0),
                                rep(1 / d$K, d$K), 5000L, 1e-12, TRUE)
  expect_true(fit$converged)
  beta <- as.vector(fit$beta)
  lambda <- as.vector(fit$lambda)

  Xm <- t(d$X)
  eta <- as.vector(exp(Xm %*% beta))
  S0 <- group_sum(eta, d$grid, d$K)
  S1 <- group_sum(Xm * eta, d$grid, d$K)
  xbar <- S1 / ifelse(S0 > 0, S0, 1)
  e <- ref_estep(d, eta, lambda)

  # The jump sizes are a fixed point of the M-step update (4).
  expect_equal(ifelse(S0 > 0, group_sum(e$W, d$grid, d$K) / S0, 0), lambda,
               tolerance = 1e-8)

  # The M-step score is zero at the solution.
  U <- colSums(e$W * (Xm - xbar[d$grid + 1L, , drop = FALSE]))
  expect_lt(max(abs(U)), 1e-6)

  # The fitted intensities reproduce the observed total exactly: summing the
  # jump equation over k gives sum_ij denom_ij = sum_ij Delta N_ij.
  expect_equal(sum(e$denom), sum(d$dN), tolerance = 1e-8)
})

test_that("a model with no covariates estimates the baseline alone", {
  d <- make_data(30)
  fit <- pcdreg(pcd(id, tstart, tstop, count) ~ 1, data = d)
  expect_length(coef(fit), 0L)
  expect_true(fit$converged)
  expect_true(all(fit$baseline$jump >= 0))
  expect_true(all(diff(fit$baseline$cumrate) >= 0))
  # No coefficients means no covariance to report, but no error either.
  expect_equal(dim(vcov(fit)), c(0L, 0L))
  expect_equal(nrow(summary(fit)$coefficients), 0L)
  expect_output(print(fit), "No coefficients")
})

test_that("the baseline is non-decreasing and the log likelihood is finite", {
  fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2,
                   data = make_data(40))
  expect_true(all(diff(fit$baseline$cumrate) >= 0))
  expect_true(all(fit$baseline$jump >= 0))
  expect_true(is.finite(fit$loglik))
  expect_equal(nrow(fit$baseline), fit$ngrid)
})

test_that("starting values do not change the solution", {
  d <- make_data(30)
  a <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d,
                 control = pcdreg_control(reltol = 1e-10))
  b <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d,
                 init = c(0.8, -0.7),
                 control = pcdreg_control(reltol = 1e-10))
  expect_equal(coef(a), coef(b), tolerance = 1e-5)
  expect_error(pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2,
                         data = d, init = 1), "has length 1")
})

test_that("collinear covariates raise an error rather than an approximation", {
  # Armadillo will happily return a least squares approximation for a singular
  # system, which would look like a fit but mean nothing. Both fitters ask it
  # not to, and report the problem instead.
  d <- make_data(30)
  d$copy <- d$x1
  expect_error(
    pcdreg(pcd(id, tstart, tstop, count) ~ x1 + copy, data = d),
    "singular"
  )
  expect_error(
    pcdreg(pcd(id, tstart, tstop, count) ~ x1 + copy, data = d, model = "mean"),
    "singular"
  )
})

test_that("failure to converge is reported rather than hidden", {
  expect_warning(
    pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2,
              data = make_data(30),
              control = pcdreg_control(maxit = 2L, reltol = 1e-12)),
    "did not converge"
  )
})
