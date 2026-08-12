# Predictions from a panel count model

Predictions from a panel count model

## Usage

``` r
# S3 method for class 'pcdfit'
predict(object, newdata, type = c("lp", "mean"), ...)
```

## Arguments

- object:

  A fitted
  [`pcdreg()`](https://www.sundayu.me/pcdreg/reference/pcdreg.md) model.

- newdata:

  Data to predict for. Defaults to the data the model was fitted to. For
  `type = "mean"` it must also contain the variables in the
  [`pcd()`](https://www.sundayu.me/pcdreg/reference/pcd.md) response,
  since the prediction follows each subject's covariate trajectory; the
  counts themselves are ignored.

- type:

  `"lp"` returns the linear predictor \\\beta' X(t)\\ for each row.
  `"mean"` returns the predicted mean number of events at every fitted
  examination time within each subject's follow-up: \\\int_0^t
  \exp(\beta' X(s)) \\ d\hat\Lambda(s)\\ for
  [`pcdreg()`](https://www.sundayu.me/pcdreg/reference/pcdreg.md), and
  \\\hat\mu(t) \exp(\beta' X(t))\\ for
  [`pcdreg()`](https://www.sundayu.me/pcdreg/reference/pcdreg.md).

- ...:

  Ignored.

## Value

For `type = "lp"`, a numeric vector with one element per row of
`newdata`. For `type = "mean"`, a
[tibble](https://tibble.tidyverse.org/reference/tibble.html) with
columns `id`, `time` and `mean`.

## Details

The rate model's predicted mean is non-decreasing in \\t\\ by
construction, because it integrates a positive quantity. The means
model's is not: it is a fitted \\\hat\mu(t)\\ rescaled by a covariate
value that may itself move up and down, so a predicted mean count that
falls is possible and is the practical drawback of that model with
time-varying covariates.

## Examples

``` r
set.seed(1)
d <- sim_pcd(60, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
head(predict(fit))
#>         1         2         3         4         5         6 
#> 0.4582126 0.4582126 0.4582126 0.5172010 0.4472924 0.4472924 
head(predict(fit, type = "mean"))
#> # A tibble: 6 × 3
#>   id     time     mean
#>   <chr> <dbl>    <dbl>
#> 1 1      0.01 1.60e-45
#> 2 1      0.02 4.76e-43
#> 3 1      0.03 9.91e-41
#> 4 1      0.04 7.85e-20
#> 5 1      0.05 1.57e-19
#> 6 1      0.07 3.61e- 1
```
