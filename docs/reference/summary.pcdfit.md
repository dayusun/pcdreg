# Summarise a fitted panel count model

Summarise a fitted panel count model

## Usage

``` r
# S3 method for class 'pcdfit'
summary(object, type = c("robust", "information", "profile"), ...)
```

## Arguments

- object:

  A fitted
  [`panelrate()`](https://dayusun.github.io/pcdreg/reference/panelrate.md)
  or
  [`panelmean()`](https://dayusun.github.io/pcdreg/reference/panelmean.md)
  model.

- type:

  Which covariance estimator to base the standard errors on. See
  [`vcov.pcdfit()`](https://dayusun.github.io/pcdreg/reference/vcov.pcdfit.md).

- ...:

  Ignored.

## Value

An object of class `"summary.pcdfit"`, whose `coefficients` component is
the usual four column table of estimates, standard errors, Wald
statistics and two sided p-values.

## Examples

``` r
set.seed(1)
d <- r_panel_count(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
fit <- panelrate(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
summary(fit)
#> 
#> Call:
#> panelrate(formula = pcd(id, tstart, tstop, count) ~ x1 + x2, 
#>     data = d)
#> 
#> Proportional rate model for panel count data 
#> Standard errors: robust sandwich
#> 
#>    Estimate Std. Error z value Pr(>|z|)
#> x1   1.1811     0.1663   7.102 1.23e-12
#> x2  -1.0206     0.1655  -6.165 7.05e-10
#> 
#> 80 subjects, 320 examinations, 645 events, 156 distinct examination times.
#> Log likelihood -426.1 in 328 EM iterations.
summary(fit, "information")
#> 
#> Call:
#> panelrate(formula = pcd(id, tstart, tstop, count) ~ x1 + x2, 
#>     data = d)
#> 
#> Proportional rate model for panel count data 
#> Standard errors: efficient information
#> 
#>    Estimate Std. Error z value Pr(>|z|)
#> x1   1.1811     0.1237   9.549   <2e-16
#> x2  -1.0206     0.1213  -8.414   <2e-16
#> 
#> 80 subjects, 320 examinations, 645 events, 156 distinct examination times.
#> Log likelihood -426.1 in 328 EM iterations.
#> 
#> Note: these standard errors assume the counts are Poisson. If the
#> counts are overdispersed they will be too small; compare with
#> `summary(fit, "robust")`.
```
