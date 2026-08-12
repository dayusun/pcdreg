# Plot the estimated baseline function

Plot the estimated baseline function

## Usage

``` r
# S3 method for class 'pcdfit'
plot(x, xlab = "Time", ylab = NULL, type = "s", ...)
```

## Arguments

- x:

  A fitted
  [`panelrate()`](https://www.sundayu.me/pcdreg/reference/panelrate.md)
  or
  [`panelmean()`](https://www.sundayu.me/pcdreg/reference/panelmean.md)
  model.

- xlab, ylab, type:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html),
  with defaults suited to a step function and to whichever baseline the
  model estimates.

- ...:

  Further graphical parameters.

## Value

`x`, invisibly.

## Examples

``` r
set.seed(1)
d <- r_panel_count(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
plot(panelrate(pcd(id, tstart, tstop, count) ~ x1 + x2, d))
```
