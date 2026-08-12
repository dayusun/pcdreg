test_that("the three argument form fills in interval starts", {
  y <- PanelCount(c(1, 1, 2), c(0.5, 1.2, 0.8), c(2, 0, 3))
  expect_s3_class(y, "PanelCount")
  expect_equal(unclass(y)[, "tstart"], c(0, 0.5, 0), ignore_attr = TRUE)
  expect_equal(unclass(y)[, "tstop"], c(0.5, 1.2, 0.8), ignore_attr = TRUE)
  expect_equal(unclass(y)[, "exam"], c(1, 1, 1), ignore_attr = TRUE)
})

test_that("the four argument form keeps the supplied intervals", {
  y <- PanelCount(c(1, 1, 1), c(0, 0.5, 0.9), c(0.5, 0.9, 1.2), c(2, NA, 0))
  expect_equal(unclass(y)[, "tstart"], c(0, 0.5, 0.9), ignore_attr = TRUE)
  expect_equal(unclass(y)[, "exam"], c(1, 0, 1), ignore_attr = TRUE)
  # Missing counts are stored as zero so the response never contains NA.
  expect_equal(unclass(y)[, "count"], c(2, 0, 0), ignore_attr = TRUE)
  expect_false(anyNA(unclass(y)))
  expect_false(any(is.na(y)))
})

test_that("identifiers of any type are accepted and labels are kept", {
  chr <- PanelCount(c("b", "a"), c(1, 1), c(0, 2))
  expect_equal(attr(chr, "labels"), c("a", "b"))
  expect_equal(unclass(chr)[, "id"], c(2, 1), ignore_attr = TRUE)

  fac <- PanelCount(factor(c("b", "a"), levels = c("b", "a")), c(1, 1), c(0, 2))
  expect_equal(attr(fac, "labels"), c("b", "a"))
})

test_that("malformed input is rejected with a specific message", {
  expect_error(PanelCount(1:2, 1:2, 1:3), "same length")
  expect_error(PanelCount(1:2, c(1, 1), c(-1, 2)), "non-negative whole number")
  expect_error(PanelCount(1:2, c(1, 1), c(0.5, 2)), "non-negative whole number")
  expect_error(PanelCount(1:2, c(1, 1), c(NA, NA)), "No examinations")
  expect_error(PanelCount(c(1, NA), c(1, 1), c(0, 2)), "must not be missing")
  # Four argument form with an interval that does not move forwards.
  expect_error(PanelCount(1:2, c(1, 1), c(0.5, 2), c(1, 2)),
               "strictly less than")
  expect_error(PanelCount(1:2, c(1, Inf), c(0, 2)), "finite")
})

test_that("subsetting preserves the class and the labels", {
  y <- PanelCount(c("a", "a", "b"), c(0.5, 1, 2), c(1, 2, 3))
  sub <- y[c(1, 3), ]
  expect_s3_class(sub, "PanelCount")
  expect_equal(attr(sub, "labels"), c("a", "b"))
  expect_equal(nrow(unclass(sub)), 2L)
  expect_equal(length(y), 3L)
})

test_that("printing shows intervals and marks non-examination rows", {
  y <- PanelCount(c(1, 1), c(0, 0.5), c(0.5, 1), c(2, NA))
  expect_match(format(y)[2], "-$")
  expect_output(print(y), "\\(0")
})

test_that("is.PanelCount discriminates", {
  expect_true(is.PanelCount(PanelCount(1, 1, 0)))
  expect_false(is.PanelCount(1:3))
})
