# Tidy a panel count model fit

Methods for the
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
and
[`generics::augment()`](https://generics.r-lib.org/reference/augment.html)
generics, so that a
[`pcdreg()`](https://www.sundayu.me/pcdreg/reference/pcdreg.md) fit can
be handled the same way as any other model in a tidy workflow.

## Usage

``` r
# S3 method for class 'pcdfit'
tidy(
  x,
  conf.int = FALSE,
  conf.level = 0.95,
  exponentiate = FALSE,
  type = c("robust", "information", "profile"),
  ...
)

# S3 method for class 'pcdfit'
glance(x, ...)

# S3 method for class 'pcdfit'
augment(x, data = NULL, newdata = NULL, ...)
```

## Arguments

- x:

  A fitted
  [`pcdreg()`](https://www.sundayu.me/pcdreg/reference/pcdreg.md) model.

- conf.int:

  Whether to add `conf.low` and `conf.high` columns.

- conf.level:

  Width of the interval, if one is requested.

- exponentiate:

  Whether to report \\e^{\beta}\\ rather than \\\beta\\. For the rate
  model this is the multiplicative effect of a one unit change in the
  covariate on the rate, read the same way as a hazard ratio; for the
  means model it is the effect on the mean. Standard errors are left on
  the log scale, as elsewhere in the tidy ecosystem.

- type:

  Which covariance estimator to base the standard errors on. See
  [`vcov.pcdfit()`](https://www.sundayu.me/pcdreg/reference/vcov.pcdfit.md);
  the default is the robust sandwich.

- ...:

  Ignored.

- data, newdata:

  Data to augment. Defaults to the data the model was fitted to. It must
  contain the [`pcd()`](https://www.sundayu.me/pcdreg/reference/pcd.md)
  response, since the fitted mean follows each subject's covariate
  trajectory.

## Value

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) gives one row
per coefficient, with columns `term`, `estimate`, `std.error`,
`statistic` and `p.value`.

[`glance()`](https://generics.r-lib.org/reference/glance.html) gives a
one row summary of the fit. `logLik` is `NA` for the means model, which
is fitted by an estimating equation and so has none.

[`augment()`](https://generics.r-lib.org/reference/augment.html) returns
`data` with four columns added:

- `.linear.predictor`:

  \\\beta' X(t)\\ on each row.

- `.fitted`:

  the predicted mean number of events by the end of the row's interval.

- `.observed`:

  the observed cumulative count at that time, on examination rows, and
  `NA` on rows that only record a covariate change.

- `.resid`:

  `.observed - .fitted`, the difference between the counts seen so far
  and the mean the model implies.

## Examples

``` r
set.seed(1)
d <- sim_pcd(60, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)

tidy(fit)
#> # A tibble: 2 × 5
#>   term  estimate std.error statistic  p.value
#>   <chr>    <dbl>     <dbl>     <dbl>    <dbl>
#> 1 x1        1.27     0.193      6.59 4.32e-11
#> 2 x2       -1.04     0.179     -5.82 5.83e- 9
tidy(fit, conf.int = TRUE, exponentiate = TRUE)
#> # A tibble: 2 × 7
#>   term  estimate std.error statistic  p.value conf.low conf.high
#>   <chr>    <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1 x1       3.58      0.193      6.59 4.32e-11    2.45      5.22 
#> 2 x2       0.354     0.179     -5.82 5.83e- 9    0.249     0.502
glance(fit)
#> # A tibble: 1 × 8
#>   model                         n nexam nevent ngrid logLik iterations converged
#>   <chr>                     <int> <int>  <dbl> <int>  <dbl>      <dbl> <lgl>    
#> 1 Proportional rate model …    60   235    487   139  -306.        184 TRUE     
augment(fit)
#> # A tibble: 283 × 10
#>       id tstart tstop count    x1    x2 .linear.predictor .fitted .observed
#>    <int>  <dbl> <dbl> <dbl> <dbl> <dbl>             <dbl>   <dbl>     <dbl>
#>  1     1   0     0.38     1 0.898 0.661             0.458    3.09         1
#>  2     1   0.38  1.07     3 0.898 0.661             0.458    8.21         4
#>  3     1   1.07  1.26    NA 0.898 0.661             0.458    9.14        NA
#>  4     1   1.26  1.7      4 0.945 0.661             0.517   11.3          8
#>  5     2   0     0.21     2 0.783 0.530             0.447    2.50         2
#>  6     2   0.21  0.8      2 0.783 0.530             0.447    6.85         4
#>  7     2   0.8   1.25     0 0.783 0.530             0.447    9.04         4
#>  8     2   1.25  1.4      1 0.783 0.530             0.447    9.78         5
#>  9     2   1.4   1.54     0 0.783 0.530             0.447   10.7          5
#> 10     2   1.54  1.58    NA 0.783 0.530             0.447   10.7         NA
#> # ℹ 273 more rows
#> # ℹ 1 more variable: .resid <dbl>
```
