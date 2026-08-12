# Estimated baseline function

Estimated baseline function

## Usage

``` r
baseline(object, ...)

# S3 method for class 'pcdfit'
baseline(object, ...)
```

## Arguments

- object:

  A fitted
  [`pcdreg()`](https://www.sundayu.me/pcdreg/reference/pcdreg.md) model.

- ...:

  Ignored.

## Value

A data frame with one row per examination time on the pooled grid. For
[`pcdreg()`](https://www.sundayu.me/pcdreg/reference/pcdreg.md) it gives
the jump sizes and the cumulative baseline rate \\\Lambda(t)\\; for
[`pcdreg()`](https://www.sundayu.me/pcdreg/reference/pcdreg.md) it gives
the baseline mean \\\mu(t)\\, which is not constrained to increase.

## Examples

``` r
set.seed(1)
d <- r_panel_count(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
head(baseline(pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, d)))
#>   time         jump   cumrate
#> 1 0.01 1.026524e-01 0.1026524
#> 2 0.02 1.598776e-81 0.1026524
#> 3 0.03 4.920883e-81 0.1026524
#> 4 0.04 0.000000e+00 0.1026524
#> 5 0.05 0.000000e+00 0.1026524
#> 6 0.06 0.000000e+00 0.1026524
head(baseline(pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, d, model = "mean")))
#>   time      mean
#> 1 0.01 0.2906572
#> 2 0.02 0.0000000
#> 3 0.03 0.0000000
#> 4 0.04 0.0000000
#> 5 0.05 0.0000000
#> 6 0.06 0.0000000
```
