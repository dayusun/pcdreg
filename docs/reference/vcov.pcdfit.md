# Covariance matrix of the estimated coefficients

Covariance matrix of the estimated coefficients

## Usage

``` r
# S3 method for class 'pcdfit'
vcov(object, type = c("robust", "information", "profile"), ...)
```

## Arguments

- object:

  A fitted
  [`panelrate()`](https://dayusun.github.io/pcdreg/reference/panelrate.md)
  or
  [`panelmean()`](https://dayusun.github.io/pcdreg/reference/panelmean.md)
  model.

- type:

  Which estimator to return. `"robust"` is the sandwich estimator that
  does not rely on the Poisson assumption and is available for both
  models. `"information"` and `"profile"` apply to
  [`panelrate()`](https://dayusun.github.io/pcdreg/reference/panelrate.md)
  only: the first is the efficient information estimator, valid under
  the Poisson assumption, and the second is the profile likelihood
  estimator, available only if the model was fitted with
  `profile = TRUE`.

- ...:

  Ignored.

## Value

A covariance matrix.

## Examples

``` r
set.seed(1)
d <- r_panel_count(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
fit <- panelrate(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
sqrt(diag(vcov(fit)))
#>        x1        x2 
#> 0.1662945 0.1655435 
sqrt(diag(vcov(fit, "information")))
#>        x1        x2 
#> 0.1236833 0.1212920 
```
