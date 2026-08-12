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

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with one
row per examination time on the pooled grid. For
[`pcdreg()`](https://www.sundayu.me/pcdreg/reference/pcdreg.md) it gives
the jump sizes and the cumulative baseline rate \\\Lambda(t)\\; for
[`pcdreg()`](https://www.sundayu.me/pcdreg/reference/pcdreg.md) it gives
the baseline mean \\\mu(t)\\, which is not constrained to increase.

## Examples

``` r
set.seed(1)
d <- sim_pcd(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
head(baseline(pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, d)))
#> # A tibble: 6 × 3
#>    time     jump cumrate
#>   <dbl>    <dbl>   <dbl>
#> 1  0.01 1.03e- 1   0.103
#> 2  0.02 1.60e-81   0.103
#> 3  0.03 4.92e-81   0.103
#> 4  0.04 0          0.103
#> 5  0.05 0          0.103
#> 6  0.06 0          0.103
head(baseline(pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, d, model = "mean")))
#> # A tibble: 6 × 2
#>    time  mean
#>   <dbl> <dbl>
#> 1  0.01 0.291
#> 2  0.02 0    
#> 3  0.03 0    
#> 4  0.04 0    
#> 5  0.05 0    
#> 6  0.06 0    
```
