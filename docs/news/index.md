# Changelog

## pcdreg 0.1.0

- First release.

- [`pcdreg()`](https://www.sundayu.me/pcdreg/reference/pcdreg.md) fits
  the semiparametric proportional rate model for panel count data with
  time-varying covariates, using the EM algorithm of Sun et al. (2024).

- `pcdreg(model = "mean")` fits the proportional means model by the
  estimating equation of Hu, Sun and Wei (2003). It is the comparator
  used in the application of Sun et al. (2024), and is included so that
  the two model families can be compared on the same data. Note that the
  two are not reparametrisations of each other when covariates vary over
  time, and that the fitted mean function is not constrained to
  increase.

- Covariance estimation for the rate model follows the paper directly.
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) supports
  `type = "robust"` (the sandwich estimator, consistent without the
  Poisson assumption), `type = "information"` (the efficient information
  estimator), and `type = "profile"` (the Murphy–van der Vaart profile
  likelihood estimator, computed only when requested via
  `pcdreg(profile = TRUE)`). The means model has no likelihood, so only
  the sandwich estimator is available for it.

- Model specification uses
  [`pcd()`](https://www.sundayu.me/pcdreg/reference/pcd.md) on the left
  hand side of the formula, following the `Surv()` convention:
  `pcd(id, time, count)` for examination-time rows and
  `pcd(id, tstart, tstop, count)` for counting-process rows carrying
  time-varying covariates.

- Both fitters return an object inheriting from `"pcdfit"`, so
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`baseline()`](https://www.sundayu.me/pcdreg/reference/baseline.md),
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) work the same
  way for either.

- [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`glance()`](https://generics.r-lib.org/reference/glance.html) and
  [`augment()`](https://generics.r-lib.org/reference/augment.html)
  methods are provided for the broom generics.
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) takes the
  same `type` argument as [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  so the tidy path reports the same standard errors as
  [`summary()`](https://rdrr.io/r/base/summary.html).
  [`augment()`](https://generics.r-lib.org/reference/augment.html) adds
  the fitted mean, the observed cumulative count and their difference to
  each row.

- [`sim_pcd()`](https://www.sundayu.me/pcdreg/reference/sim_pcd.md)
  simulates from the design of the paper’s simulation study.

- [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  methods return `ggplot` objects for the data, drawn as one tile per
  examination interval shaded by its count, and for the fitted baseline.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws them
  and returns the object invisibly.

- [`baseline()`](https://www.sundayu.me/pcdreg/reference/baseline.md),
  `predict(type = "mean")` and
  [`sim_pcd()`](https://www.sundayu.me/pcdreg/reference/sim_pcd.md)
  return tibbles.

- Errors, warnings and messages use cli and carry condition classes, so
  a particular failure can be caught with
  [`tryCatch()`](https://rdrr.io/r/base/conditions.html) rather than by
  matching on message text.
