# Describe panel count observations

Creates the response object used on the left hand side of a
[`pcdreg()`](https://www.sundayu.me/pcdreg/reference/pcdreg.md) formula,
in the same spirit as
[`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html).

## Usage

``` r
pcd(id, time, time2, count)
```

## Arguments

- id:

  Subject identifier. May be numeric, character or a factor; rows
  sharing an identifier are treated as one subject.

- time:

  Examination time in the three argument form, or the start of the
  interval in the four argument form.

- time2:

  The end of the interval in the four argument form; the event count in
  the three argument form.

- count:

  Number of events observed since the previous examination. Use `NA` on
  rows that only record a change in a covariate rather than an
  examination. Omitted in the three argument form.

## Value

An object of class `"pcd"`: a numeric matrix with columns `id`,
`tstart`, `tstop`, `count` and `exam`, carrying the original identifier
labels in the `"labels"` attribute. Missing counts are stored as
`count = 0` with `exam = 0`, so that the response contains no `NA` and
`na.action` handling applies to the covariates alone.

## Details

There are two forms, distinguished by the number of arguments. Which one
you need depends on a single question: does any covariate change value
during follow-up?

`pcd(id, time, count)` is for data with one row per examination, which
is all that time-invariant covariates need. Interval starts are filled
in as the previous examination time of the same subject, with the first
interval starting at zero.

`pcd(id, tstart, tstop, count)` is the counting process form, needed as
soon as a covariate moves. Each row gives an interval `(tstart, tstop]`
over which the covariates in that row are constant, exactly as for a
time-varying
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) fit.
The `count` column distinguishes the two kinds of row:

- a number means the row ends at an **examination**, and that many
  events occurred since the previous examination;

- `NA` means the row ends at a **covariate change** rather than an
  examination.

There is no separate event indicator; the presence or absence of a count
is the indicator.

## Counts are increments, not totals

`count` must be the number of events since the previous examination,
\\\Delta N\_{ij}\\, not the running total \\N_i(T\_{ij})\\. Data
recorded cumulatively must be differenced within subject first. This is
not something the package can detect for you, and getting it wrong
inflates the event count and biases the coefficients rather than raising
an error. A quick check after fitting is that `fit$nevent` matches the
number of events you believe you have. See
[`vignette("data-preparation")`](https://www.sundayu.me/pcdreg/articles/data-preparation.md)
for the conversion.

## What the data must satisfy

Checked before any arithmetic, naming the offending subject when a check
fails:

- follow-up starts at time zero, since the models integrate from there;

- intervals within a subject are contiguous, with no gaps or overlaps;

- no subject has two rows ending at the same time;

- every subject has at least one examination;

- counts are non-negative whole numbers.

Row order is irrelevant: the data are sorted internally. Identifiers may
be numeric, character or factor, and need not be consecutive.

A subject's follow-up ends at their last examination, because the
likelihood runs only over \\j = 1, \ldots, J_i\\. Rows lying entirely
beyond it carry no information and are dropped, with a message saying
how many. A censoring row appended after the final visit is the usual
reason for seeing it.

## See also

[`pcdreg()`](https://www.sundayu.me/pcdreg/reference/pcdreg.md) and
[`vignette("data-preparation")`](https://www.sundayu.me/pcdreg/articles/data-preparation.md)
for converting data into this shape.

## Examples

``` r
# One row per examination: subject 1 seen at 0.5 and 1.2 with 2 then 0
# events, subject 2 seen once at 0.8 with 3.
pcd(c(1, 1, 2), c(0.5, 1.2, 0.8), c(2, 0, 3))
#> [1] 1: (0.0, 0.5] 2 1: (0.5, 1.2] 0 2: (0.0, 0.8] 3

# Counting process form. Subject 1's covariates change at t = 0.9, part way
# through the examination interval (0.5, 1.2], so that interval needs two
# rows: the first ends at the change and carries NA, the second ends at the
# examination and carries its count.
pcd(
  id     = c(1, 1, 1, 2),
  time   = c(0.0, 0.5, 0.9, 0.0),
  time2  = c(0.5, 0.9, 1.2, 0.8),
  count  = c(2, NA, 0, 3)
)
#> [1] 1: (0.0, 0.5] 2 1: (0.5, 0.9] - 1: (0.9, 1.2] 0 2: (0.0, 0.8] 3

# Used on the left hand side of a formula, which is how you will normally
# meet it.
set.seed(1)
d <- sim_pcd(40)
fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
c(subjects = fit$n, examinations = fit$nexam, events = fit$nevent)
#>     subjects examinations       events 
#>           40          159          272 
```
