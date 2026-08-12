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

test_that("the data plot draws for both orderings and subsets large studies", {
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f), add = TRUE)
  set.seed(1)
  d <- r_panel_count(30)
  y <- with(d, pcd(id, tstart, tstop, count))

  for (by in c("followup", "events", "none")) {
    grDevices::png(f)
    expect_silent(plot(y, order_by = by))
    grDevices::dev.off()
    expect_gt(file.size(f), 1000)
  }

  # Large studies are thinned to keep the picture readable, and say so.
  big <- with(r_panel_count(80), pcd(id, tstart, tstop, count))
  grDevices::png(f)
  expect_message(plot(big, max_subjects = 20), "20 of 80")
  grDevices::dev.off()

  # Character identifiers survive to the axis labels.
  d$id <- paste0("s", d$id)
  grDevices::png(f)
  expect_silent(plot(with(d, pcd(id, tstart, tstop, count))))
  grDevices::dev.off()
})
