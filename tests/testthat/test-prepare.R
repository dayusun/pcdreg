# A hand-checkable data set.  Two subjects, pooled examination times
# 0.5, 1.0, 1.5.  Subject 1 is examined at 0.5 and 1.5 with a covariate change
# at 1.2; subject 2 is examined at 1.0 only.
toy <- data.frame(
  id     = c(1, 1, 1, 2),
  tstart = c(0.0, 0.5, 1.2, 0.0),
  tstop  = c(0.5, 1.2, 1.5, 1.0),
  count  = c(2, NA, 3, 1),
  x      = c(10, 10, 20, 30)
)

test_that("rows are expanded onto the pooled grid up to each subject's end", {
  d <- prep(toy, PanelCount(id, tstart, tstop, count) ~ x)

  expect_equal(d$times, c(0.5, 1.0, 1.5))
  expect_equal(d$K, 3L)
  expect_equal(d$n, 2L)

  # Subject 1 is followed to 1.5, so all three grid times apply; subject 2 is
  # followed to 1.0, so only the first two do.
  expect_equal(d$subj, c(0L, 0L, 0L, 1L, 1L))
  expect_equal(d$grid, c(0L, 1L, 2L, 0L, 1L))

  # Subject 1's panels are (0, 0.5] and (0.5, 1.5]; subject 2 has one panel.
  expect_equal(d$panel, c(0L, 1L, 1L, 2L, 2L))
  expect_equal(d$dN, c(2, 3, 1))
  expect_equal(d$panelsubj, c(0L, 0L, 1L))
})

test_that("the covariate in force at each grid time is picked up", {
  d <- prep(toy, PanelCount(id, tstart, tstop, count) ~ x)
  # Subject 1: x = 10 at t = 0.5 and t = 1.0, then 20 at t = 1.5 (the change
  # happens at 1.2).  Subject 2: x = 30 throughout.
  expect_equal(as.vector(d$X), c(10, 10, 20, 30, 30))
})

test_that("intervals beyond the last examination are dropped", {
  extra <- rbind(toy, data.frame(id = 1, tstart = 1.5, tstop = 2.0,
                                 count = NA, x = 99))
  expect_message(d <- prep(extra, PanelCount(id, tstart, tstop, count) ~ x),
                 "beyond the last examination")
  expect_false(99 %in% as.vector(d$X))
})

test_that("structural problems in the data are caught", {
  gapped <- toy
  gapped$tstart[2] <- 0.7
  expect_error(prep(gapped, PanelCount(id, tstart, tstop, count) ~ x),
               "contiguous")

  late <- toy
  late$tstart[1] <- 0.1
  expect_error(prep(late, PanelCount(id, tstart, tstop, count) ~ x),
               "start at time 0")

  none <- toy
  none$count <- c(2, NA, 3, NA)
  expect_error(prep(none, PanelCount(id, tstart, tstop, count) ~ x),
               "at least one examination")
})

test_that("the two PanelCount forms describe the same data", {
  # One row per examination, no time-varying covariates.
  simple <- data.frame(id = c(1, 1, 2), tstop = c(0.5, 1.5, 1.0),
                       count = c(2, 3, 1), x = c(10, 10, 30))
  a <- prep(simple, PanelCount(id, tstop, count) ~ x)

  explicit <- data.frame(id = c(1, 1, 2), tstart = c(0, 0.5, 0),
                         tstop = c(0.5, 1.5, 1.0), count = c(2, 3, 1),
                         x = c(10, 10, 30))
  b <- prep(explicit, PanelCount(id, tstart, tstop, count) ~ x)

  expect_equal(a[c("subj", "grid", "panel", "dN", "times")],
               b[c("subj", "grid", "panel", "dN", "times")])
  expect_equal(a$X, b$X)
})
