dat <- local({
  set.seed(707)
  sim_pcd(50, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
})
fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = dat)
mfit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = dat,
               model = "mean")

test_that("tidy agrees with the extractors it is a view of", {
  td <- tidy(fit)
  expect_s3_class(td, "tbl_df")
  expect_named(td, c("term", "estimate", "std.error", "statistic", "p.value"))
  expect_equal(td$term, names(coef(fit)))
  expect_equal(td$estimate, unname(coef(fit)))
  expect_equal(td$std.error, unname(sqrt(diag(vcov(fit)))))

  # The same numbers summary() prints, so the two cannot drift apart.
  tab <- summary(fit)$coefficients
  expect_equal(td$statistic, unname(tab[, "z value"]))
  expect_equal(td$p.value, unname(tab[, "Pr(>|z|)"]))
})

test_that("tidy reports the covariance it was asked for", {
  robust <- tidy(fit, type = "robust")$std.error
  info <- tidy(fit, type = "information")$std.error
  expect_equal(info, unname(sqrt(diag(vcov(fit, "information")))))
  expect_false(isTRUE(all.equal(robust, info)))

  # The means model has only the sandwich, and tidy must not paper over that.
  expect_error(tidy(mfit, type = "information"),
               class = "pcdreg_error_vcov_type")
})

test_that("tidy intervals match confint and exponentiate correctly", {
  td <- tidy(fit, conf.int = TRUE)
  ci <- confint(fit)
  expect_equal(td$conf.low, unname(ci[, 1]))
  expect_equal(td$conf.high, unname(ci[, 2]))

  ex <- tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  expect_equal(ex$estimate, exp(td$estimate))
  expect_equal(ex$conf.low, exp(td$conf.low))
  expect_equal(ex$conf.high, exp(td$conf.high))
  # Standard errors stay on the log scale, as in the rest of the tidy world.
  expect_equal(ex$std.error, td$std.error)

  expect_error(tidy(fit, conf.int = TRUE, conf.level = 1.5),
               class = "pcdreg_error_conf_level")
})

test_that("glance summarises the fit and admits the means model has no logLik", {
  g <- glance(fit)
  expect_s3_class(g, "tbl_df")
  expect_equal(nrow(g), 1L)
  expect_equal(g$n, fit$n)
  expect_equal(g$nevent, fit$nevent)
  expect_equal(g$ngrid, fit$ngrid)
  expect_equal(g$logLik, as.numeric(fit$loglik))
  expect_true(g$converged)

  gm <- glance(mfit)
  expect_equal(nrow(gm), 1L)
  expect_true(is.na(gm$logLik))
  # logLik() itself still refuses, so glance is the only place this is softened.
  expect_error(logLik(mfit), class = "pcdreg_error_no_loglik")
})

test_that("augment attaches fitted means, observed counts and residuals", {
  aug <- augment(fit)
  expect_s3_class(aug, "tbl_df")
  expect_equal(nrow(aug), nrow(dat))
  expect_true(all(c(".linear.predictor", ".fitted", ".observed", ".resid") %in%
                    names(aug)))

  expect_equal(aug$.resid, aug$.observed - aug$.fitted)
  expect_equal(aug$.linear.predictor, predict(fit, dat, type = "lp"))

  # Counts are only observed at examinations.
  expect_true(all(is.na(aug$.observed[is.na(dat$count)])))
  expect_false(anyNA(aug$.observed[!is.na(dat$count)]))

  # The observed column is a running total within subject, so it ends at each
  # subject's total number of events.
  totals <- tapply(dat$count, dat$id, sum, na.rm = TRUE)
  last <- tapply(aug$.observed, aug$id, function(v) max(v, na.rm = TRUE))
  expect_equal(as.numeric(last), as.numeric(totals))

  # The rate model integrates a positive quantity, so fitted means never fall.
  by_subject <- split(aug$.fitted, aug$id)
  expect_true(all(vapply(by_subject, function(v) !is.unsorted(v), logical(1))))
})

test_that("augment matches predict at the grid times it is read from", {
  aug <- augment(fit)
  traj <- predict(fit, dat, type = "mean")
  one <- aug[aug$id == aug$id[1], ]
  tr <- traj[traj$id == one$id[1], ]
  k <- findInterval(one$tstop, tr$time)
  expect_equal(one$.fitted[k > 0], tr$mean[k[k > 0]])
})

test_that("augment needs the response, and says so", {
  expect_error(augment(fit, data = dat[, c("x1", "x2")]),
               class = "pcdreg_error_no_response")
})

test_that("the tidiers work for the means model too", {
  td <- tidy(mfit, conf.int = TRUE)
  expect_equal(td$estimate, unname(coef(mfit)))
  expect_equal(nrow(augment(mfit)), nrow(dat))
})
