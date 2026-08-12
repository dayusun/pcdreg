# panelrate 0.1.0

* First release.

* `panelrate()` fits the semiparametric proportional rate model for panel count
  data with time-varying covariates, using the EM algorithm of Sun et al. (2024).

* Covariance estimation follows the paper directly. `vcov()` supports
  `type = "robust"` (the sandwich estimator, consistent without the Poisson
  assumption), `type = "information"` (the efficient information estimator), and
  `type = "profile"` (the Murphy--van der Vaart profile likelihood estimator,
  computed only when requested via `panelrate(profile = TRUE)`).

* Model specification uses `PanelCount()` on the left hand side of the formula,
  following the `Surv()` convention: `PanelCount(id, time, count)` for
  examination-time rows and `PanelCount(id, tstart, tstop, count)` for
  counting-process rows carrying time-varying covariates.
