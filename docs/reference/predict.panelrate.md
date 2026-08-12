# Predictions from a panel count model

Predictions from a panel count model

## Usage

``` r
# S3 method for class 'panelrate'
predict(object, newdata, type = c("lp", "mean"), ...)

# S3 method for class 'panelmean'
predict(object, newdata, type = c("lp", "mean"), ...)
```

## Arguments

- object:

  A fitted
  [`panelrate()`](https://www.sundayu.me/pcdreg/reference/panelrate.md)
  or
  [`panelmean()`](https://www.sundayu.me/pcdreg/reference/panelmean.md)
  model.

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
  [`panelrate()`](https://www.sundayu.me/pcdreg/reference/panelrate.md),
  and \\\hat\mu(t) \exp(\beta' X(t))\\ for
  [`panelmean()`](https://www.sundayu.me/pcdreg/reference/panelmean.md).

- ...:

  Ignored.

## Value

For `type = "lp"`, a numeric vector with one element per row of
`newdata`. For `type = "mean"`, a data frame with columns `id`, `time`
and `mean`.

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
d <- r_panel_count(60, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
fit <- panelrate(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
head(predict(fit))
#>         1         2         3         4         5         6 
#> 0.4582126 0.4582126 0.4582126 0.5172010 0.4472924 0.4472924 
head(predict(fit, type = "mean"))
#>   id time         mean
#> 1  1 0.01 1.595922e-45
#> 2  1 0.02 4.762185e-43
#> 3  1 0.03 9.914545e-41
#> 4  1 0.04 7.846571e-20
#> 5  1 0.05 1.569314e-19
#> 6  1 0.07 3.611803e-01
```
