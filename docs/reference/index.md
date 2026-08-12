# Package index

## Fitting

The two models. They take the same data and the same formula, so a model
can be refitted under the other by changing the function name alone.

- [`pcdreg()`](https://www.sundayu.me/pcdreg/reference/pcdreg.md) :
  Semiparametric regression for panel count data

## Specifying the data

The response constructor, used on the left hand side of the formula. See
the data preparation article for the layouts it accepts.

- [`pcd()`](https://www.sundayu.me/pcdreg/reference/pcd.md) : Describe
  panel count observations
- [`is.pcd()`](https://www.sundayu.me/pcdreg/reference/is.pcd.md) : Test
  for a panel count response
- [`plot(`*`<pcd>`*`)`](https://www.sundayu.me/pcdreg/reference/plot.pcd.md)
  : Plot panel count data

## Inspecting a fit

Both fitters return an object inheriting from `pcdfit`, so these work
the same way for either.

- [`vcov(`*`<pcdfit>`*`)`](https://www.sundayu.me/pcdreg/reference/vcov.pcdfit.md)
  : Covariance matrix of the estimated coefficients
- [`summary(`*`<pcdfit>`*`)`](https://www.sundayu.me/pcdreg/reference/summary.pcdfit.md)
  : Summarise a fitted panel count model
- [`baseline()`](https://www.sundayu.me/pcdreg/reference/baseline.md) :
  Estimated baseline function
- [`predict(`*`<pcdfit>`*`)`](https://www.sundayu.me/pcdreg/reference/predict.pcdfit.md)
  : Predictions from a panel count model
- [`plot(`*`<pcdfit>`*`)`](https://www.sundayu.me/pcdreg/reference/plot.pcdfit.md)
  : Plot the estimated baseline function

## Tuning

- [`pcdreg_control()`](https://www.sundayu.me/pcdreg/reference/pcdreg_control.md)
  : Tuning parameters for the fitting algorithms

## Simulation

The design of the simulation study in Sun et al. (2024), useful for
examples and for checking the covariance estimators against known truth.

- [`r_panel_count()`](https://www.sundayu.me/pcdreg/reference/r_panel_count.md)
  : Simulate panel count data with a time-varying covariate

## Package

- [`pcdreg-package`](https://www.sundayu.me/pcdreg/reference/pcdreg-package.md)
  : pcdreg: Semiparametric Regression for Panel Count Data
