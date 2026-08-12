dat <- local({
  set.seed(3030)
  r_panel_count(40, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
})
fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = dat)

test_that("the extractor methods agree with the fit", {
  expect_named(coef(fit), c("x1", "x2"))
  expect_equal(nobs(fit), 40L)
  expect_equal(dim(vcov(fit)), c(2L, 2L))
  expect_equal(dimnames(vcov(fit)), list(c("x1", "x2"), c("x1", "x2")))

  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), fit$loglik)
  expect_equal(attr(ll, "df"), 2L + fit$ngrid)
  expect_equal(attr(ll, "nobs"), 40L)
})

test_that("summary reports the covariance it was asked for", {
  s <- summary(fit)
  expect_s3_class(s, "summary.pcdfit")
  expect_equal(s$coefficients[, "Std. Error"],
               sqrt(diag(vcov(fit, "robust"))), ignore_attr = TRUE)
  expect_equal(summary(fit, "information")$coefficients[, "Std. Error"],
               sqrt(diag(vcov(fit, "information"))), ignore_attr = TRUE)

  # The p-values are the usual two sided Wald ones.
  z <- s$coefficients[, "Estimate"] / s$coefficients[, "Std. Error"]
  expect_equal(s$coefficients[, "Pr(>|z|)"], 2 * pnorm(-abs(z)),
               ignore_attr = TRUE)
})

test_that("printing works and warns about the Poisson assumption", {
  expect_output(print(fit), "Coefficients")
  expect_output(print(summary(fit)), "robust sandwich")
  expect_output(print(summary(fit, "information")), "assume the counts are")
  expect_output(print(summary(fit)), "40 subjects")
})

test_that("baseline returns the estimated cumulative rate", {
  b <- baseline(fit)
  expect_named(b, c("time", "jump", "cumrate"))
  expect_equal(b$cumrate, cumsum(b$jump))
  expect_equal(nrow(b), fit$ngrid)
  expect_false(is.unsorted(b$time))
})

test_that("predict returns linear predictors and non-decreasing means", {
  lp <- predict(fit, dat)
  expect_length(lp, nrow(dat))
  expect_equal(lp, as.vector(as.matrix(dat[, c("x1", "x2")]) %*% coef(fit)),
               ignore_attr = TRUE)

  m <- predict(fit, dat, type = "mean")
  expect_named(m, c("id", "time", "mean"))
  # The rate model guarantees a non-decreasing predicted mean even though x1
  # jumps around, which the proportional means model cannot.
  by_subject <- split(m$mean, m$id)
  expect_true(all(vapply(by_subject, function(v) all(diff(v) >= 0), logical(1))))
  expect_true(all(m$mean >= 0))
})

test_that("predict rejects data without the response for type mean", {
  expect_error(predict(fit, dat[, c("x1", "x2")], type = "mean"))
})

test_that("plot draws without error", {
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  expect_silent(plot(fit))
})

test_that("factors and interactions are handled and no intercept is fitted", {
  d <- dat
  d$g <- factor(ifelse(d$id %% 2 == 0, "even", "odd"))
  f <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + g, data = d)
  expect_false("(Intercept)" %in% names(coef(f)))
  expect_named(coef(f), c("x1", "godd"))

  fi <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 * x2, data = d)
  expect_named(coef(fi), c("x1", "x2", "x1:x2"))
})

test_that("subset is honoured", {
  half <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = dat,
                    subset = id <= 20)
  expect_equal(nobs(half), 20L)
})

test_that("a non pcd response is refused", {
  expect_error(pcdreg(count ~ x1, data = dat), "pcd")
})
