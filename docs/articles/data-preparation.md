# Preparing data and specifying the model

``` r
library(pcdreg)
```

Data preparation is most of the work in fitting these models, and most
mistakes occur here. This vignette describes the inputs required by the
two fitters, formula specification, conversion from common layouts, and
the error messages.

## What panel count data is

A subject is examined at times $T_{i1} < \ldots < T_{iJ_{i}}$. Each
examination gives the number of events since the previous examination,
but not their times. Two quantities are easy to confuse:

- $\Delta N_{ij}$, the **increment** — events in
  $(T_{i,j - 1},T_{ij}\rbrack$;
- $N_{i}\left( T_{ij} \right)$, the **cumulative** count — events up to
  $T_{ij}$.

**`pcdreg` always wants the increment.** If your data records cumulative
counts, difference them first; there is a recipe below. Passing
cumulative counts is not an error the package can detect, and it will
quietly give you the wrong answer.

## The response

[`pcd()`](https://www.sundayu.me/pcdreg/reference/pcd.md) builds the
left-hand side of the formula, as `Surv()` does for survival models.
Choose its form according to one question: **do any covariates change
value during follow-up?**

### No time-varying covariates: one row per examination

``` r
visits <- data.frame(
  id    = c(1, 1, 1, 2, 2),
  time  = c(0.5, 1.1, 1.9, 0.7, 1.4),
  count = c(2, 0, 3, 1, 4),
  age   = c(24, 24, 24, 31, 31),
  arm   = c("drug", "drug", "drug", "placebo", "placebo")
)

pcd(visits$id, visits$time, visits$count)
#> [1] 1: (0.0, 0.5] 2 1: (0.5, 1.1] 0 1: (1.1, 1.9] 3 2: (0.0, 0.7] 1
#> [5] 2: (0.7, 1.4] 4
```

Interval starts come from each subject’s previous examination, with
follow-up beginning at zero. The printed intervals show this. Subject 1
contributes $(0,0.5\rbrack$, $(0.5,1.1\rbrack$ and $(1.1,1.9\rbrack$
with 2, 0 and 3 events. Fit the model with

``` r
pcdreg(pcd(id, time, count) ~ age + arm, data = visits)
```

but not with these five rows and two subjects: `age` and `arm` are
perfectly confounded, so the information matrix is singular. Later
examples use real fits.

### Time-varying covariates: counting process rows

When a covariate changes between examinations, use one row for each
interval with constant covariates. This is the layout `coxph()` uses for
time-dependent covariates.

``` r
cp <- data.frame(
  id     = c(1, 1, 1, 2),
  tstart = c(0.0, 0.5, 0.9, 0.0),
  tstop  = c(0.5, 0.9, 1.2, 0.8),
  count  = c(2, NA, 0, 3),
  dose   = c(10, 10, 25, 25)
)

pcd(cp$id, cp$tstart, cp$tstop, cp$count)
#> [1] 1: (0.0, 0.5] 2 1: (0.5, 0.9] - 1: (0.9, 1.2] 0 2: (0.0, 0.8] 3
```

Subject 1 was examined at 0.5 with 2 events and at 1.2 with none. The
dose changed from 10 to 25 at time 0.9, within the second examination
interval.

The `count` field distinguishes the two row types:

- a **number** means “this row ends at an examination, and this many
  events happened since the previous examination”;
- **`NA`** means “this row ends at a covariate change, not an
  examination”.

There is no separate event indicator.

### Which form, in one line

Use the three argument form when every covariate is fixed at baseline.
Use the four argument form as soon as a covariate changes. When in
doubt, use the four argument form. It can express anything the three
argument form can.

## The formula

The right-hand side is an ordinary model formula.

``` r
pcd(id, tstart, tstop, count) ~ x1 + x2          # main effects
pcd(id, tstart, tstop, count) ~ x1 * x2          # interaction
pcd(id, tstart, tstop, count) ~ x1 + log(x2)     # transformation
pcd(id, tstart, tstop, count) ~ x1 + factor(site)
pcd(id, tstart, tstop, count) ~ .                # everything else
```

Factors are expanded with your contrast settings, and
[`predict()`](https://rdrr.io/r/stats/predict.html) reuses the coding
the fit used.

``` r
set.seed(1)
d <- sim_pcd(60)
d$site <- factor(c("north", "south", "east")[1 + d$id %% 3])

fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + site, data = d)
names(coef(fit))
#> [1] "x1"        "sitenorth" "sitesouth"
```

The dot expands to every column outside the response, so `id`, `tstart`,
`tstop` and `count` are excluded automatically:

``` r
names(coef(pcdreg(pcd(id, tstart, tstop, count) ~ ., data = d)))
#> [1] "x1"        "x2"        "sitenorth" "sitesouth"
```

### There is no intercept

A constant column cannot be told apart from a rescaling of the baseline
function, so both fitters drop the intercept and report the baseline
separately through
[`baseline()`](https://www.sundayu.me/pcdreg/reference/baseline.md).
Writing `~ x1 + x2 - 1` changes nothing, and adding `+ 0` is
unnecessary.

Adding a constant $c$ to a covariate leaves $\beta$ unchanged and
rescales the baseline by $e^{- c\beta}$. Centring covariates changes the
baseline, not the coefficients.

### Row order does not matter

Data are sorted internally by subject and time, so `arrange()` is
unnecessary. Subject identifiers need not be integers or consecutive.
Numbers, characters, and factors all work; `predict(type = "mean")`
returns the original labels.

## Recipes for common layouts

### Cumulative counts

This is the most common and damaging mismatch. If `count` accumulates,
difference it within subject before fitting.

``` r
# A frame whose counts accumulate: 2, 2, 5 rather than 2, 0, 3.
cumdat <- data.frame(
  id    = c(1, 1, 1, 2, 2),
  time  = c(0.5, 1.1, 1.9, 0.7, 1.4),
  count = c(2, 2, 5, 1, 5),
  age   = c(24, 24, 24, 31, 31)
)

cumdat$increment <- ave(cumdat$count, cumdat$id,
                        FUN = function(v) v - c(0, v[-length(v)]))
cumdat[, c("id", "time", "count", "increment")]
#>   id time count increment
#> 1  1  0.5     2         2
#> 2  1  1.1     2         0
#> 3  1  1.9     5         3
#> 4  2  0.7     1         1
#> 5  2  1.4     5         4
```

Then fit with `increment`. With `NA`s on covariate-change rows,
difference only the examination rows:

``` r
at <- !is.na(dat$count)
dat$count[at] <- ave(dat$count[at], dat$id[at],
                     FUN = function(v) v - c(0, v[-length(v)]))
```

### Data held in the older `panelEM` layout

This layout repeats each panel count on every row of its interval and
marks examinations with a separate indicator. A repeated count carries
no information beyond the examination row, so conversion takes one line:

``` r
new <- old
new$count[old$event == 0] <- NA
pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = new)
```

### Examinations and covariate histories in separate tables

If visits and covariate changes are in separate tables, stack them and
carry covariates forward. Examination rows retain their counts;
covariate-change rows receive `NA`.

``` r
exams <- data.frame(id = c(1, 1), time = c(0.5, 1.2), count = c(2, 0))
covar <- data.frame(id = 1, time = 0.9, dose = 25)
start <- data.frame(id = 1, time = 0, dose = 10)

both <- merge(
  data.frame(id = c(exams$id, covar$id), tstop = c(exams$time, covar$time),
             count = c(exams$count, NA)),
  rbind(start[, c("id", "time", "dose")], covar[, c("id", "time", "dose")]),
  by.x = c("id", "tstop"), by.y = c("id", "time"), all.x = TRUE
)
both <- both[order(both$id, both$tstop), ]
both$dose <- ave(both$dose, both$id,
                 FUN = function(v) { for (i in seq_along(v))
                   if (is.na(v[i]) && i > 1) v[i] <- v[i - 1]; v })
both$tstart <- ave(both$tstop, both$id,
                   FUN = function(v) c(0, v[-length(v)]))
both[, c("id", "tstart", "tstop", "count", "dose")]
#>   id tstart tstop count dose
#> 1  1    0.0   0.5     2   NA
#> 2  1    0.5   0.9    NA   25
#> 3  1    0.9   1.2     0   25
```

### One row per subject, several visit columns

Reshape wide data to long first with
[`stats::reshape()`](https://rdrr.io/r/stats/reshape.html) or
[`tidyr::pivot_longer()`](https://tidyr.tidyverse.org/reference/pivot_longer.html).
Use one row per subject and visit, then the three argument form.

## What the data must satisfy

Both fitters perform these checks before arithmetic and identify the
offending subject when a check fails.

**Follow-up starts at zero.** The models integrate from 0, so the first
interval of every subject must begin there.

``` r
bad <- data.frame(id = 1, tstart = 0.2, tstop = 1, count = 3, x = 1)
pcdreg(pcd(id, tstart, tstop, count) ~ x, data = bad)
#> Error in `prepare_panel()`:
#> ! Follow-up must start at time 0.
#> ✖ Subject "1" starts at 0.2.
#> ℹ Shift the time scale so that 0 is the start of follow-up.
```

For calendar dates, or follow-up starting at enrolment, subtract each
subject’s own origin first.

**Intervals are contiguous.** Within a subject each interval must begin
exactly where the previous one ended: no gaps, no overlaps.

``` r
bad <- data.frame(id = c(1, 1), tstart = c(0, 0.7), tstop = c(0.5, 1),
                  count = c(2, 1), x = c(1, 1))
pcdreg(pcd(id, tstart, tstop, count) ~ x, data = bad)
#> Error in `prepare_panel()`:
#> ! Intervals within a subject must be contiguous.
#> ✖ Subject "1" has a gap or overlap at time 0.7.
#> ℹ Carry the last covariate value forward rather than leaving the stretch out.
```

A gap usually indicates a covariate history that does not cover all
follow-up. Carry the last value forward instead of omitting the
interval.

**Every subject has at least one examination.** A subject with `count`
equal to `NA` on every row contributes nothing and is more likely a data
error than an intentional record.

**Counts are non-negative whole numbers.**
[`pcd()`](https://www.sundayu.me/pcdreg/reference/pcd.md) refuses other
values. This catches rates and averages passed by mistake.

``` r
pcd(1:2, c(1, 2), c(0.5, 1.2))
#> Error in `pcd()`:
#> ! `count` must be a non-negative whole number at examination times.
```

**One row per subject may end at any given time.** Duplicate rows at the
same subject and time are rejected instead of silently combined.

### What gets dropped

Follow-up ends at a subject’s last examination. The likelihood runs over
$j = 1,\ldots,J_{i}$ and has no information beyond $T_{iJ_{i}}$. Rows
wholly after that time are dropped, with a message reporting how many.

``` r
# A censoring row appended after subject 1's final visit.
last1 <- max(d$tstop[d$id == 1])
trailing <- transform(d[d$id == 1, ][1, ],
                      tstart = last1, tstop = last1 + 0.5, count = NA)

fit2 <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2,
                  data = rbind(d, trailing))
#> Dropping 1 row beyond the last examination of their subject.
```

The usual cause is a censoring row appended after the last visit. This
is harmless: the interval contains no information.

### Missing values

Missing covariates are handled by `na.action`, which by default deletes
rows — and deleting a row in the middle of a subject breaks contiguity,
so you will get the gap error rather than a silently different fit.
Resolve missing covariates before fitting, by carrying values forward or
by dropping whole subjects.

## What the package does with your rows

This processing explains the error messages.

The pooled grid $t_{1} < \ldots < t_{K}$ contains all distinct
examination times. For the **rate model**, each subject’s covariate
trajectory is evaluated at each grid time during follow-up. The
complete, contiguous trajectory is therefore required. The **means
model** uses only each subject’s examination times and is much cheaper.

``` r
fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
c(subjects = fit$n, examinations = fit$nexam,
  grid_times = fit$ngrid, events = fit$nevent)
#>     subjects examinations   grid_times       events 
#>           60          235          139          487
```

Ties across subjects reduce the grid, which drives the rate model’s
cost. Examination times recorded to the day or week tie naturally. Times
kept at full numerical precision do not, so the grid grows with sample
size. Round to the precision at which visits were scheduled.

## Checking before you fit

These quantities are a quick sanity check. `nexam` should equal the
expected number of examinations, and `nevent` the total number of
events. An unexpectedly large `nevent` usually indicates cumulative
counts.

``` r
c(rows_supplied = nrow(d),
  examination_rows = sum(!is.na(d$count)),
  events = sum(d$count, na.rm = TRUE))
#>    rows_supplied examination_rows           events 
#>              283              235              487
```

## See also

[`vignette("pcdreg")`](https://www.sundayu.me/pcdreg/articles/pcdreg.md)
for the models themselves, the choice of covariance estimator, and the
comparison between the rate and means models.
