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
- [`autoplot(`*`<pcd>`*`)`](https://www.sundayu.me/pcdreg/reference/autoplot.pcd.md)
  [`plot(`*`<pcd>`*`)`](https://www.sundayu.me/pcdreg/reference/autoplot.pcd.md)
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
- [`autoplot(`*`<pcdfit>`*`)`](https://www.sundayu.me/pcdreg/reference/autoplot.pcdfit.md)
  [`plot(`*`<pcdfit>`*`)`](https://www.sundayu.me/pcdreg/reference/autoplot.pcdfit.md)
  : Plot the estimated baseline function

## Tidy output

Methods for the broom generics, so a fit drops into a tidy workflow: one
row per coefficient, a one row summary of the fit, and the data with
fitted means and residuals attached.

- [`tidy(`*`<pcdfit>`*`)`](https://www.sundayu.me/pcdreg/reference/pcdreg-tidiers.md)
  [`glance(`*`<pcdfit>`*`)`](https://www.sundayu.me/pcdreg/reference/pcdreg-tidiers.md)
  [`augment(`*`<pcdfit>`*`)`](https://www.sundayu.me/pcdreg/reference/pcdreg-tidiers.md)
  : Tidy a panel count model fit
- [`reexports`](https://www.sundayu.me/pcdreg/reference/reexports.md)
  [`autoplot`](https://www.sundayu.me/pcdreg/reference/reexports.md)
  [`tidy`](https://www.sundayu.me/pcdreg/reference/reexports.md)
  [`glance`](https://www.sundayu.me/pcdreg/reference/reexports.md)
  [`augment`](https://www.sundayu.me/pcdreg/reference/reexports.md) :
  Objects exported from other packages

## Tuning

- [`pcdreg_control()`](https://www.sundayu.me/pcdreg/reference/pcdreg_control.md)
  : Tuning parameters for the fitting algorithms

## Simulation

The design of the simulation study in Sun et al. (2024), useful for
examples and for checking the covariance estimators against known truth.

- [`sim_pcd()`](https://www.sundayu.me/pcdreg/reference/sim_pcd.md) :
  Simulate panel count data with a time-varying covariate

## Package

- [`pcdreg-package`](https://www.sundayu.me/pcdreg/reference/pcdreg-package.md)
  : pcdreg: Semiparametric Regression for Panel Count Data
