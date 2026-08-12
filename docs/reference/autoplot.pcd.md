# Plot panel count data

Draws the observation pattern as tiles: one row per subject, one tile
per examination interval spanning the time it covers, shaded by the
number of events recorded in it.

## Usage

``` r
# S3 method for class 'pcd'
autoplot(
  object,
  order_by = c("followup", "events", "none"),
  max_subjects = 60L,
  palette = NULL,
  ...
)

# S3 method for class 'pcd'
plot(x, ...)
```

## Arguments

- object, x:

  A [`pcd()`](https://www.sundayu.me/pcdreg/reference/pcd.md) object.

- order_by:

  How to order subjects up the vertical axis. `"followup"` sorts by the
  time of the last examination, which makes the censoring pattern
  legible; `"events"` sorts by the total number of events; `"none"`
  keeps the order of the identifiers.

- max_subjects:

  Draw at most this many subjects, taken as an evenly spaced sample so
  the picture stays readable for large studies. A message reports how
  many were shown.

- palette:

  Fill colours, lightest first. The default is a single-hue sequential
  ramp, which is where a magnitude belongs.

- ...:

  Ignored.

## Value

A
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object, so it can be titled, faceted or themed like any other.
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
returns it; [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
draws it and returns it invisibly, so it can be used for its side effect
like any other [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
method.

## Details

The tile is the unit of observation. A panel count is attributed to the
interval \\(T\_{i,j-1}, T\_{ij}\]\\ rather than to an instant, so
drawing the intervals rather than points at the examination times is
what the data actually says: a wide pale tile and a narrow dark one can
carry the same count while meaning very different things about the rate.

Two features decide how an analysis will go and both are visible here.
Whether examination times line up across subjects sets the size of the
pooled grid and so the cost of fitting the rate model; and how ragged
the right hand edge is shows how much of the tail of the baseline rests
on a handful of subjects.

## See also

[`autoplot.pcdfit()`](https://www.sundayu.me/pcdreg/reference/autoplot.pcdfit.md)
for the fitted baseline function.

## Examples

``` r
set.seed(1)
d <- sim_pcd(40)
autoplot(with(d, pcd(id, tstart, tstop, count)))

```
