make_data <- function(n = 40, seed = 20240202, ...) {
  set.seed(seed)
  sim_pcd(n, beta = c(1, -1), lambda = function(t) 8 / (1 + t), ...)
}

test_that("Omega and S match the pure R versions of the paper's formulas", {
  d <- prep(make_data(30))
  fit <- ref_em(d)
  got <- pcdreg:::covariance_cpp(d$X, d$subj, d$grid, d$panel, d$dN,
                                    d$panelsubj, d$n, d$K, fit$beta,
                                    fit$lambda)
  want <- ref_covariance(d, fit$beta, fit$lambda)

  expect_equal(got$Omega, want$Omega, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(got$S, want$S, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(got$scores, want$scores, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("Omega and S have the structure the theory requires", {
  fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2,
                   data = make_data(50),
                   control = pcdreg_control(reltol = 1e-10))

  # The nonparametric estimate sets many jump sizes to exactly zero, so an
  # examination interval can carry a total intensity that is positive but
  # subnormal.  Inverting one of those overflows, and the matching numerator has
  # underflowed to zero, so the product is NaN unless such intervals are
  # dropped.  Whether it bites depends on the last bit of the arithmetic, so it
  # showed up on one platform and not another.
  expect_true(all(is.finite(fit$Omega)))
  expect_true(all(is.finite(fit$S)))
  expect_true(all(is.finite(vcov(fit, "robust"))))
  expect_true(all(is.finite(vcov(fit, "information"))))
  expect_true(any(fit$baseline$jump == 0))

  expect_equal(fit$Omega, t(fit$Omega))
  expect_equal(fit$S, t(fit$S))
  # Omega is negative definite and S, being a sum of outer products, is
  # positive definite.
  expect_true(all(eigen(fit$Omega, only.values = TRUE)$values < 0))
  expect_true(all(eigen(fit$S, only.values = TRUE)$values > 0))

  # S is exactly the empirical covariance of the per subject score residuals,
  # and those residuals sum to zero because the score vanishes at the solution.
  expect_equal(fit$S, crossprod(fit$scores) / fit$n, ignore_attr = TRUE)
  expect_lt(max(abs(colSums(fit$scores))), 1e-5)
  expect_equal(nrow(fit$scores), fit$n)
})

test_that("the reported covariances are the formulas in Section 4", {
  fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2,
                   data = make_data(40))
  n <- fit$n
  expect_equal(vcov(fit, "robust"),
               solve(fit$Omega) %*% fit$S %*% solve(fit$Omega) / n,
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(vcov(fit, "information"), solve(fit$S) / n,
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_true(all(diag(vcov(fit, "robust")) > 0))
  expect_true(all(diag(vcov(fit, "information")) > 0))
})

test_that("the profile covariance is only available on request", {
  fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2,
                   data = make_data(30))
  expect_error(vcov(fit, "profile"), "profile = TRUE")
})

test_that("the information and profile estimators agree under Poisson data", {
  skip_on_cran()
  # Theorem 3 of the paper says both estimate the same matrix S; the simulation
  # study reports them as nearly identical replication by replication.
  fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2,
                   data = make_data(120, seed = 7), profile = TRUE,
                   control = pcdreg_control(reltol = 1e-10))
  info <- sqrt(diag(vcov(fit, "information")))
  prof <- sqrt(diag(vcov(fit, "profile")))
  expect_equal(info, prof, tolerance = 0.05)
})

test_that("overdispersion separates the robust estimator from the others", {
  skip_on_cran()
  # With a gamma frailty the rate model still holds but the counts are far from
  # Poisson, so the information estimator should understate the variability
  # while the robust one should not.
  fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2,
                   data = make_data(150, seed = 11, frailty = 1),
                   control = pcdreg_control(reltol = 1e-9))
  robust <- sqrt(diag(vcov(fit, "robust")))
  info <- sqrt(diag(vcov(fit, "information")))
  expect_true(all(robust > info))
})

test_that("confidence intervals line up with the chosen covariance", {
  fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2,
                   data = make_data(40))
  ci <- confint(fit)
  se <- sqrt(diag(vcov(fit, "robust")))
  expect_equal(ci[, 1], coef(fit) - qnorm(0.975) * se, ignore_attr = TRUE)
  expect_equal(ci[, 2], coef(fit) + qnorm(0.975) * se, ignore_attr = TRUE)

  ci_info <- confint(fit, type = "information")
  expect_false(isTRUE(all.equal(ci[, 1], ci_info[, 1])))
  expect_equal(rownames(confint(fit, parm = "x1")), "x1")
})
