test_that("the three argument form fills in interval starts", {
  y <- pcd(c(1, 1, 2), c(0.5, 1.2, 0.8), c(2, 0, 3))
  expect_s3_class(y, "pcd")
  expect_equal(unclass(y)[, "tstart"], c(0, 0.5, 0), ignore_attr = TRUE)
  expect_equal(unclass(y)[, "tstop"], c(0.5, 1.2, 0.8), ignore_attr = TRUE)
  expect_equal(unclass(y)[, "exam"], c(1, 1, 1), ignore_attr = TRUE)
})

test_that("the four argument form keeps the supplied intervals", {
  y <- pcd(c(1, 1, 1), c(0, 0.5, 0.9), c(0.5, 0.9, 1.2), c(2, NA, 0))
  expect_equal(unclass(y)[, "tstart"], c(0, 0.5, 0.9), ignore_attr = TRUE)
  expect_equal(unclass(y)[, "exam"], c(1, 0, 1), ignore_attr = TRUE)
  # Missing counts are stored as zero so the response never contains NA.
  expect_equal(unclass(y)[, "count"], c(2, 0, 0), ignore_attr = TRUE)
  expect_false(anyNA(unclass(y)))
  expect_false(any(is.na(y)))
})

test_that("identifiers of any type are accepted and labels are kept", {
  chr <- pcd(c("b", "a"), c(1, 1), c(0, 2))
  expect_equal(attr(chr, "labels"), c("a", "b"))
  expect_equal(unclass(chr)[, "id"], c(2, 1), ignore_attr = TRUE)

  fac <- pcd(factor(c("b", "a"), levels = c("b", "a")), c(1, 1), c(0, 2))
  expect_equal(attr(fac, "labels"), c("b", "a"))
})

test_that("malformed input is rejected with a specific message", {
  expect_error(pcd(1:2, 1:2, 1:3), "same length")
  expect_error(pcd(1:2, c(1, 1), c(-1, 2)), "non-negative whole number")
  expect_error(pcd(1:2, c(1, 1), c(0.5, 2)), "non-negative whole number")
  expect_error(pcd(1:2, c(1, 1), c(NA, NA)), "No examinations")
  expect_error(pcd(c(1, NA), c(1, 1), c(0, 2)), "must not be missing")
  # Four argument form with an interval that does not move forwards.
  expect_error(pcd(1:2, c(1, 1), c(0.5, 2), c(1, 2)),
               "strictly less than")
  expect_error(pcd(1:2, c(1, Inf), c(0, 2)), "finite")
})

# The conditions carry classes so that callers can catch a particular failure
# rather than matching on message text, which is free to change.
test_that("validation failures carry condition classes", {
  expect_error(pcd(1:2, 1:2, 1:3), class = "pcdreg_error_length")
  expect_error(pcd(1:2, c(1, 1), c(-1, 2)), class = "pcdreg_error_count")
  expect_error(pcd(1:2, c(1, 1), c(NA, NA)), class = "pcdreg_error_no_exam")
  expect_error(pcd(c(1, NA), c(1, 1), c(0, 2)),
               class = "pcdreg_error_missing_id")
  expect_error(pcd(1:2, c(1, 1), c(0.5, 2), c(1, 2)),
               class = "pcdreg_error_interval")
  expect_error(pcd(1:2, c(1, Inf), c(0, 2)), class = "pcdreg_error_time")

  bad <- data.frame(id = 1, tstart = 0.2, tstop = 1, count = 3, x = 1)
  expect_error(pcdreg(pcd(id, tstart, tstop, count) ~ x, data = bad),
               class = "pcdreg_error_origin")
})

test_that("subsetting preserves the class and the labels", {
  y <- pcd(c("a", "a", "b"), c(0.5, 1, 2), c(1, 2, 3))
  sub <- y[c(1, 3), ]
  expect_s3_class(sub, "pcd")
  expect_equal(attr(sub, "labels"), c("a", "b"))
  expect_equal(nrow(unclass(sub)), 2L)
  expect_equal(length(y), 3L)
})

test_that("printing shows intervals and marks non-examination rows", {
  y <- pcd(c(1, 1), c(0, 0.5), c(0.5, 1), c(2, NA))
  expect_match(format(y)[2], "-$")
  expect_output(print(y), "\\(0")
})

test_that("is.pcd discriminates", {
  expect_true(is.pcd(pcd(1, 1, 0)))
  expect_false(is.pcd(1:3))
})

test_that("the data plot is a ggplot with one tile per examination interval", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  d <- sim_pcd(30)
  y <- with(d, pcd(id, tstart, tstop, count))

  p <- autoplot(y)
  expect_s3_class(p, "ggplot")
  # One tile per examination, and the tiles span intervals rather than points.
  expect_equal(nrow(p$data), sum(!is.na(d$count)))
  expect_true(all(p$data$tstop > p$data$tstart))
  # Each subject's first tile starts at zero and the rest abut.
  first <- !duplicated(p$data$y)
  expect_true(all(p$data$tstart[first] == 0))

  for (by in c("followup", "events", "none")) {
    expect_s3_class(autoplot(y, order_by = by), "ggplot")
  }
  expect_error(autoplot(y, order_by = "nonsense"))

  # plot() draws the same picture as a side effect and hands the object back
  # invisibly, so it works from a script rather than only at the console.
  drawn <- tempfile(fileext = ".pdf")
  grDevices::pdf(drawn)
  on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    unlink(drawn)
  }, add = TRUE)
  res <- withVisible(plot(y))
  expect_false(res$visible)
  expect_s3_class(res$value, "ggplot")

  # It renders without warnings about missing or dropped values.
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f), add = TRUE)
  expect_silent(ggplot2::ggsave(f, p, width = 6, height = 4, dpi = 72))
  expect_gt(file.size(f), 1000)
})

test_that("large studies are thinned and identifiers survive", {
  skip_if_not_installed("ggplot2")
  set.seed(2)
  big <- with(sim_pcd(80), pcd(id, tstart, tstop, count))
  expect_message(p <- autoplot(big, max_subjects = 20), "20 of 80")
  expect_equal(length(unique(p$data$y)), 20L)

  d <- sim_pcd(20)
  d$id <- paste0("s", d$id)
  expect_s3_class(autoplot(with(d, pcd(id, tstart, tstop, count))), "ggplot")
})

test_that("the baseline plot is a ggplot for both models", {
  skip_if_not_installed("ggplot2")
  set.seed(3)
  d <- sim_pcd(40)
  for (mod in c("rate", "mean")) {
    fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d,
                  model = mod)
    p <- autoplot(fit)
    expect_s3_class(p, "ggplot")
    # Starts at the origin, one step per grid time.
    expect_equal(nrow(p$data), fit$ngrid + 1L)
    expect_equal(p$data$value[1], 0)
    expect_s3_class(plot(fit), "ggplot")
  }
})
