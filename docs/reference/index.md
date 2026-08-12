# Package index

## Fitting

The two models. They take the same data and the same formula, so a model
can be refitted under the other by changing the function name alone.

- [`panelrate()`](https://dayusun.github.io/pcdreg/reference/panelrate.md)
  : Fit the proportional rate model to panel count data
- [`panelmean()`](https://dayusun.github.io/pcdreg/reference/panelmean.md)
  : Fit the proportional means model to panel count data

## Specifying the data

The response constructor, used on the left hand side of the formula. See
the data preparation article for the layouts it accepts.

- [`pcd()`](https://dayusun.github.io/pcdreg/reference/pcd.md) :
  Describe panel count observations
- [`is.pcd()`](https://dayusun.github.io/pcdreg/reference/is.pcd.md) :
  Test for a panel count response

## Inspecting a fit

Both fitters return an object inheriting from `pcdfit`, so these work
the same way for either.

- [`vcov(`*`<pcdfit>`*`)`](https://dayusun.github.io/pcdreg/reference/vcov.pcdfit.md)
  : Covariance matrix of the estimated coefficients
- [`summary(`*`<pcdfit>`*`)`](https://dayusun.github.io/pcdreg/reference/summary.pcdfit.md)
  : Summarise a fitted panel count model
- [`baseline()`](https://dayusun.github.io/pcdreg/reference/baseline.md)
  : Estimated baseline function
- [`predict(`*`<panelrate>`*`)`](https://dayusun.github.io/pcdreg/reference/predict.panelrate.md)
  [`predict(`*`<panelmean>`*`)`](https://dayusun.github.io/pcdreg/reference/predict.panelrate.md)
  : Predictions from a panel count model
- [`plot(`*`<pcdfit>`*`)`](https://dayusun.github.io/pcdreg/reference/plot.pcdfit.md)
  : Plot the estimated baseline function

## Tuning

- [`panelrate_control()`](https://dayusun.github.io/pcdreg/reference/panelrate_control.md)
  : Tuning parameters for the EM algorithm
- [`panelmean_control()`](https://dayusun.github.io/pcdreg/reference/panelmean_control.md)
  : Tuning parameters for the means model

## Simulation

The design of the simulation study in Sun et al. (2024), useful for
examples and for checking the covariance estimators against known truth.

- [`r_panel_count()`](https://dayusun.github.io/pcdreg/reference/r_panel_count.md)
  : Simulate panel count data with a time-varying covariate

## Package

- [`pcdreg`](https://dayusun.github.io/pcdreg/reference/pcdreg-package.md)
  [`pcdreg-package`](https://dayusun.github.io/pcdreg/reference/pcdreg-package.md)
  : pcdreg: Semiparametric Regression for Panel Count Data
