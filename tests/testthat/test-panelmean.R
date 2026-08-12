make_data <- function(n = 40, seed = 20240303, ...) {
  set.seed(seed)
  r_panel_count(n, beta = c(1, -1), lambda = function(t) 8 / (1 + t), ...)
}

fit_mean <- function(d, ...) {
  panelmean(PanelCount(id, tstart, tstop, count) ~ x1 + x2, data = d,
            control = panelmean_control(reltol = 1e-12), ...)
}

test_that("the compiled means model matches an independent pure R version", {
  d <- prep(make_data(30))
  got <- pcdreg:::mean_fit_cpp(d$exam_X, d$exam_subj, d$exam_grid, d$cN,
                               d$n, d$K, c(0, 0), 100L, 1e-12)
  want <- ref_mean(d)
  expect_true(got$converged)
  expect_equal(as.vector(got$beta), as.vector(want$beta), tolerance = 1e-10)
  expect_equal(as.vector(got$mu), as.vector(want$mu), tolerance = 1e-10)
})

test_that("the estimating equation is solved at the reported estimates", {
  d <- prep(make_data(35))
  fit <- fit_mean(make_data(35))
  beta <- coef(fit)

  Xm <- t(d$exam_X)
  grid <- d$exam_grid
  eta <- as.vector(exp(Xm %*% beta))
  S0 <- group_sum(eta, grid, d$K)
  S1 <- group_sum(Xm * eta, grid, d$K)
  xbar <- S1 / ifelse(S0 > 0, S0, 1)

  U <- colSums(d$cN * (Xm - xbar[grid + 1L, , drop = FALSE]))
  expect_lt(max(abs(U)), 1e-7)

  # muhat is the ratio that makes the baseline equation hold at each grid time.
  mu <- ifelse(S0 > 0, group_sum(d$cN, grid, d$K) / S0, 0)
  expect_equal(fit$baseline$mean, mu, tolerance = 1e-8)
})

test_that("the sandwich pieces have the structure the theory requires", {
  fit <- fit_mean(make_data(60))
  expect_true(all(is.finite(fit$Omega)))
  expect_true(all(is.finite(fit$S)))
  expect_equal(fit$Omega, t(fit$Omega))
  expect_equal(fit$S, t(fit$S))
  expect_true(all(eigen(fit$Omega, only.values = TRUE)$values < 0))
  expect_true(all(eigen(fit$S, only.values = TRUE)$values > 0))

  # S is the empirical covariance of the per subject contributions, and those
  # sum to zero because both the score and the baseline equation hold.
  expect_equal(fit$S, crossprod(fit$scores) / fit$n, ignore_attr = TRUE)
  expect_lt(max(abs(colSums(fit$scores))), 1e-6)

  expect_equal(vcov(fit), solve(fit$Omega) %*% fit$S %*% solve(fit$Omega) /
                 fit$n, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("only the sandwich covariance is offered, and no log likelihood", {
  fit <- fit_mean(make_data(30))
  expect_error(vcov(fit, "information"), "estimating equation rather than a likelihood")
  expect_error(vcov(fit, "profile"), "estimating equation rather than a likelihood")
  expect_error(logLik(fit), "no log likelihood")
  expect_null(fit$loglik)
})

test_that("the shared methods work for both models", {
  d <- make_data(40)
  mfit <- fit_mean(d)
  rfit <- panelrate(PanelCount(id, tstart, tstop, count) ~ x1 + x2, data = d)

  for (fit in list(mfit, rfit)) {
    expect_s3_class(fit, "pcdfit")
    expect_named(coef(fit), c("x1", "x2"))
    expect_equal(nobs(fit), 40L)
    expect_equal(dim(vcov(fit)), c(2L, 2L))
    expect_equal(dim(confint(fit)), c(2L, 2L))
    expect_equal(nrow(summary(fit)$coefficients), 2L)
    expect_true(is.finite(sum(baseline(fit)[[2]])))
  }
  expect_output(print(summary(mfit)), "means model")
  expect_output(print(summary(rfit)), "rate model")
  expect_output(print(mfit), "Newton iterations")
  expect_output(print(rfit), "EM iterations")
})

test_that("predictions follow each subject's trajectory", {
  d <- make_data(40)
  fit <- fit_mean(d)

  lp <- predict(fit, d)
  expect_equal(lp, as.vector(as.matrix(d[, c("x1", "x2")]) %*% coef(fit)),
               ignore_attr = TRUE)

  m <- predict(fit, d, type = "mean")
  expect_named(m, c("id", "time", "mean"))
  expect_true(all(m$mean >= 0))
  # The means model prediction is muhat(t) exp(beta'X(t)) rather than an
  # integral, so unlike the rate model it need not increase.
  expect_equal(nrow(m), sum(!is.na(m$mean)))
})

test_that("results do not depend on row order or subject labels", {
  d <- make_data(35)
  base <- fit_mean(d)

  set.seed(1)
  shuffled <- d[sample(nrow(d)), ]
  rownames(shuffled) <- NULL
  expect_equal(coef(fit_mean(shuffled)), coef(base), tolerance = 1e-8)

  relabelled <- d
  relabelled$id <- paste0("s", sprintf("%03d", d$id))
  expect_equal(coef(fit_mean(relabelled)), coef(base), tolerance = 1e-8)
})

test_that("rescaling a covariate rescales its coefficient", {
  d <- make_data(35)
  base <- fit_mean(d)
  scaled <- d
  scaled$x1 <- d$x1 * 10
  got <- fit_mean(scaled)
  expect_equal(unname(coef(got)["x1"]), unname(coef(base)["x1"]) / 10,
               tolerance = 1e-7)
  expect_equal(unname(coef(got)["x2"]), unname(coef(base)["x2"]),
               tolerance = 1e-7)
})

test_that("the means model needs only the examination rows", {
  # prepare_panel(expand = FALSE) skips the grid expansion, so the examination
  # level arrays must be identical either way.
  d <- make_data(30)
  mf <- stats::model.frame(PanelCount(id, tstart, tstop, count) ~ x1 + x2, d)
  X <- stats::model.matrix(attr(mf, "terms"), mf)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  y <- stats::model.response(mf)

  full <- pcdreg:::prepare_panel(y, X, expand = TRUE)
  lean <- pcdreg:::prepare_panel(y, X, expand = FALSE)
  shared <- c("exam_X", "exam_subj", "exam_grid", "dN", "cN", "n", "K", "times")
  expect_equal(lean[shared], full[shared])
  expect_null(lean$X)

  # The cumulative counts really are cumulative within subject.
  expect_true(all(lean$cN >= lean$dN))
  expect_equal(sum(lean$dN), lean$cN[length(lean$cN)] +
                 sum(lean$cN[c(diff(lean$exam_subj) != 0, FALSE)]))
})

test_that("failure to converge is reported", {
  expect_warning(
    panelmean(PanelCount(id, tstart, tstop, count) ~ x1 + x2,
              data = make_data(30),
              control = panelmean_control(maxit = 1L, reltol = 1e-14)),
    "did not converge"
  )
})
