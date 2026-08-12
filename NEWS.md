# pcdreg 0.1.0

* First release.

* `pcdreg()` fits the semiparametric proportional rate model for panel count
  data with time-varying covariates, using the EM algorithm of Sun et al. (2024).

* `pcdreg(model = "mean")` fits the proportional means model by the estimating equation of
  Hu, Sun and Wei (2003). It is the comparator used in the application of Sun et
  al. (2024), and is included so that the two model families can be compared on
  the same data. Note that the two are not reparametrisations of each other when
  covariates vary over time, and that the fitted mean function is not
  constrained to increase.

* Covariance estimation for the rate model follows the paper directly. `vcov()`
  supports `type = "robust"` (the sandwich estimator, consistent without the
  Poisson assumption), `type = "information"` (the efficient information
  estimator), and `type = "profile"` (the Murphy--van der Vaart profile
  likelihood estimator, computed only when requested via
  `pcdreg(profile = TRUE)`). The means model has no likelihood, so only the
  sandwich estimator is available for it.

* Model specification uses `pcd()` on the left hand side of the formula,
  following the `Surv()` convention: `pcd(id, time, count)` for
  examination-time rows and `pcd(id, tstart, tstop, count)` for
  counting-process rows carrying time-varying covariates.

* Both fitters return an object inheriting from `"pcdfit"`, so `coef()`,
  `vcov()`, `confint()`, `summary()`, `baseline()`, `predict()` and `plot()`
  work the same way for either.
