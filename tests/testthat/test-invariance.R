dat <- local({
  set.seed(505)
  r_panel_count(35, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
})

ctrl <- panelrate_control(reltol = 1e-10, maxit = 20000L)

fit_of <- function(d) {
  panelrate(PanelCount(id, tstart, tstop, count) ~ x1 + x2, data = d,
            control = ctrl)
}

base_fit <- fit_of(dat)

test_that("the row order of the input does not matter", {
  set.seed(1)
  shuffled <- dat[sample(nrow(dat)), ]
  rownames(shuffled) <- NULL
  got <- fit_of(shuffled)
  expect_equal(coef(got), coef(base_fit), tolerance = 1e-7)
  expect_equal(got$baseline, base_fit$baseline, tolerance = 1e-7)
  expect_equal(vcov(got), vcov(base_fit), tolerance = 1e-7)
})

test_that("relabelling the subjects does not matter", {
  relabelled <- dat
  relabelled$id <- paste0("subject-", sprintf("%03d", dat$id))
  got <- fit_of(relabelled)
  expect_equal(coef(got), coef(base_fit), tolerance = 1e-7)
  expect_equal(got$baseline, base_fit$baseline, tolerance = 1e-7)
})

test_that("shifting a covariate shifts nothing but the baseline", {
  # The rate model has no intercept, so replacing x2 by x2 + c must leave
  # exp(beta'X) dLambda unchanged: the coefficients stay put and the baseline
  # is rescaled by exp(-c * beta2).
  shifted <- dat
  shifted$x2 <- dat$x2 + 5
  got <- fit_of(shifted)
  expect_equal(coef(got), coef(base_fit), tolerance = 1e-6)
  expect_equal(got$baseline$cumrate,
               base_fit$baseline$cumrate * exp(-5 * coef(base_fit)[["x2"]]),
               tolerance = 1e-5)
})

test_that("rescaling a covariate rescales its coefficient", {
  scaled <- dat
  scaled$x1 <- dat$x1 * 10
  got <- fit_of(scaled)
  expect_equal(unname(coef(got)["x1"]), unname(coef(base_fit)["x1"]) / 10,
               tolerance = 1e-6)
  expect_equal(unname(coef(got)["x2"]), unname(coef(base_fit)["x2"]),
               tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(got)))["x1"]),
               unname(sqrt(diag(vcov(base_fit)))["x1"]) / 10, tolerance = 1e-6)
})

test_that("duplicating every subject leaves the estimates alone", {
  # Twice the data with the same empirical distribution gives the same
  # estimates, and standard errors smaller by a factor of sqrt(2).
  twin <- dat
  twin$id <- twin$id + 1000
  doubled <- rbind(dat, twin)
  got <- fit_of(doubled)
  expect_equal(coef(got), coef(base_fit), tolerance = 1e-6)
  expect_equal(sqrt(diag(vcov(got))), sqrt(diag(vcov(base_fit))) / sqrt(2),
               tolerance = 1e-4)
})
