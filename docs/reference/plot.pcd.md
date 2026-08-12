# Plot panel count data

Draws the observation pattern: one row per subject, their follow-up as a
line, and a circle at each examination with area proportional to the
number of events recorded there.

## Usage

``` r
# S3 method for class 'pcd'
plot(
  x,
  order_by = c("followup", "events", "none"),
  max_subjects = 60L,
  cex_max = 2.5,
  xlab = "Time",
  ylab = "Subject",
  col = "#2B4C7E",
  ...
)
```

## Arguments

- x:

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

- cex_max:

  Radius scaling for the largest count.

- xlab, ylab, col, ...:

  Passed to the underlying plotting calls.

## Value

`x`, invisibly.

## Details

What the picture is for: panel count data have two features that decide
which model and which estimator make sense, and both are visible here.
Whether examination times line up across subjects determines the size of
the pooled grid and so the cost of fitting the rate model, and whether
follow-up ends at widely differing times determines how much of the
baseline is estimated from few subjects.

Open circles mark examinations where no events were recorded, so a
subject with a long line of open circles is contributing information
about the baseline without contributing events.

## See also

[`plot.pcdfit()`](https://www.sundayu.me/pcdreg/reference/plot.pcdfit.md)
for the fitted baseline function.

## Examples

``` r
set.seed(1)
d <- r_panel_count(40)
plot(with(d, pcd(id, tstart, tstop, count)))

```
