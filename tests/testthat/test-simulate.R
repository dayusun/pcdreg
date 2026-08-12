test_that("the simulator returns well formed panel count data", {
  set.seed(99)
  d <- r_panel_count(30)
  expect_named(d, c("id", "tstart", "tstop", "count", "x1", "x2"))
  expect_equal(length(unique(d$id)), 30L)
  expect_true(all(d$tstart < d$tstop))
  expect_true(all(d$count[!is.na(d$count)] >= 0))
  expect_true(all(d$count[!is.na(d$count)] %% 1 == 0))

  # Every subject starts at zero and has contiguous intervals ending at an
  # examination.
  by_subject <- split(d, d$id)
  expect_true(all(vapply(by_subject, function(s) s$tstart[1] == 0, logical(1))))
  expect_true(all(vapply(by_subject,
                         function(s) all(s$tstart[-1] == s$tstop[-nrow(s)]),
                         logical(1))))
  expect_true(all(vapply(by_subject,
                         function(s) !is.na(s$count[nrow(s)]), logical(1))))
  expect_true(all(vapply(by_subject, function(s) any(!is.na(s$count)),
                         logical(1))))
})

test_that("rounding ties the examination times across subjects", {
  set.seed(7)
  tied <- r_panel_count(60, digits = 2)
  untied <- r_panel_count(60, digits = NA)
  ntied <- length(unique(tied$tstop[!is.na(tied$count)]))
  nuntied <- length(unique(untied$tstop[!is.na(untied$count)]))
  expect_lt(ntied, nuntied)
  expect_lte(ntied, 200L)
})

test_that("the time-varying covariate really varies and x2 does not", {
  set.seed(11)
  d <- r_panel_count(40)
  by_subject <- split(d, d$id)
  expect_true(any(vapply(by_subject,
                         function(s) length(unique(s$x1)) > 1, logical(1))))
  expect_true(all(vapply(by_subject,
                         function(s) length(unique(s$x2)) == 1, logical(1))))
})

test_that("a frailty produces overdispersion relative to Poisson", {
  skip_on_cran()
  set.seed(2024)
  poisson <- r_panel_count(400, frailty = 0)
  mixed <- r_panel_count(400, frailty = 1)
  total <- function(d) vapply(split(d$count, d$id),
                              function(v) sum(v, na.rm = TRUE), numeric(1))
  ratio <- function(d) { x <- total(d); stats::var(x) / mean(x) }
  expect_gt(ratio(mixed), ratio(poisson))
})

test_that("the estimator recovers the truth on simulated data", {
  skip_on_cran()
  set.seed(4321)
  d <- r_panel_count(400, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
  fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
  se <- sqrt(diag(vcov(fit)))
  # Within three standard errors of the truth.
  expect_lt(abs(coef(fit)[["x1"]] - 1), 3 * se[["x1"]])
  expect_lt(abs(coef(fit)[["x2"]] + 1), 3 * se[["x2"]])

  # And the baseline tracks Lambda(t) = 8 log(1 + t).
  truth <- 8 * log(1 + fit$baseline$time)
  expect_lt(max(abs(fit$baseline$cumrate - truth)) / max(truth), 0.2)
})

test_that("simulator arguments are validated", {
  expect_error(r_panel_count(10, beta = 1), "length\\(beta\\)")
  expect_error(r_panel_count(10, lambda = function(t) -t), "positive and finite")
})
