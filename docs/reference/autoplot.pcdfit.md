# Plot the estimated baseline function

Plot the estimated baseline function

## Usage

``` r
# S3 method for class 'pcdfit'
autoplot(object, ...)

# S3 method for class 'pcdfit'
plot(x, ...)
```

## Arguments

- object, x:

  A fitted
  [`pcdreg()`](https://www.sundayu.me/pcdreg/reference/pcdreg.md) model.

- ...:

  Ignored.

## Value

A
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object: a step function of the cumulative baseline rate for the rate
model, or of the baseline mean for the means model.
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
returns it; [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
draws it and returns it invisibly, so it can be used for its side effect
like any other [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
method.

## Details

The rate model's \\\hat\Lambda\\ increases by construction. The means
model's \\\hat\mu\\ need not, and with covariates that fluctuate it
generally does not; it is drawn as estimated rather than forced upwards,
because that behaviour is the point of the comparison between the two
models.

## Examples

``` r
set.seed(1)
d <- sim_pcd(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
autoplot(pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, d))
```
