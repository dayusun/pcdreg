# Fit the proportional rate model to panel count data

Fits \$\$E\[dN(t) \mid X(t)\] = \exp(\beta' X(t)) \\ d\Lambda(t)\$\$ by
nonparametric maximum likelihood, treating \\\Lambda\\ as a step
function with a jump at each observed examination time. Covariates may
vary over time.

## Usage

``` r
panelrate(
  formula,
  data,
  subset,
  na.action,
  control = panelrate_control(),
  profile = FALSE,
  init = NULL
)
```

## Arguments

- formula:

  A formula whose left hand side is a
  [`pcd()`](https://dayusun.github.io/pcdreg/reference/pcd.md) object,
  for example `pcd(id, tstart, tstop, count) ~ x1 + x2`. No intercept is
  fitted: it is absorbed into the baseline \\\Lambda\\.

- data:

  A data frame containing the variables in `formula`.

- subset, na.action:

  Handled as in [`stats::lm()`](https://rdrr.io/r/stats/lm.html). Note
  that missing covariate values delete rows, which can break the
  contiguity of a subject's intervals; resolve them before fitting.

- control:

  A list from
  [`panelrate_control()`](https://dayusun.github.io/pcdreg/reference/panelrate_control.md).

- profile:

  Whether to also compute the profile likelihood covariance. This is
  provided for comparison with the literature and costs one extra
  baseline-only fit per coefficient, so it is off by default.

- init:

  Optional starting values for the coefficients.

## Value

An object of class `"panelrate"`, inheriting from `"pcdfit"`, with
components including `coefficients`, `baseline` (a data frame of jump
sizes and the cumulative baseline rate), `Omega` and `S` (the two
matrices of Section 4 of the paper), `vcov` (a list of the available
covariance estimates), `loglik`, and convergence information.

## Details

The likelihood maximised is the one a nonhomogeneous Poisson process
would give. That assumption is a working device: it makes the EM
algorithm available, but the estimator remains consistent and
asymptotically normal when it fails, and the default covariance
estimator remains valid.

Three covariance estimators are available through
[`vcov.pcdfit()`](https://dayusun.github.io/pcdreg/reference/vcov.pcdfit.md).
`"robust"` is the sandwich \\\Omega^{-1} S \Omega^{-1} / n\\, which does
not rely on the Poisson assumption and is the default. `"information"`
is \\S^{-1} / n\\, valid only under that assumption but a fast
substitute for the profile likelihood when it holds. `"profile"` is the
numerical profile likelihood estimator of Murphy and van der Vaart
(2000), computed only when `profile = TRUE`. When the counts are
overdispersed relative to Poisson, the latter two understate the
standard errors, sometimes severely.

The profile estimator is included for comparison with the literature
rather than for routine use. It rests on a numerical derivative taken
with a step of order \\n^{-1/2}\\, which is coarse at small sample
sizes, and it costs one baseline-only fit per coefficient.
`"information"` estimates the same matrix far more cheaply and more
stably.

## Specifying the model

The left hand side is always a
[`pcd()`](https://dayusun.github.io/pcdreg/reference/pcd.md) object,
which carries the subject, the interval and the event count. The right
hand side is an ordinary model formula, so interactions,
transformations, factors and `.` all behave as usual:

    pcd(id, tstart, tstop, count) ~ x1 + x2
    pcd(id, tstart, tstop, count) ~ x1 * x2 + log(x3)
    pcd(id, time, count)          ~ .

The dot expands to every column that is not part of the response, so the
identifier, times and counts are excluded automatically.

No intercept is fitted. A constant column cannot be distinguished from a
rescaling of \\\Lambda\\, so it is dropped and the baseline is reported
separately by
[`baseline()`](https://dayusun.github.io/pcdreg/reference/baseline.md).
Adding `- 1` changes nothing. One consequence is that shifting a
covariate by a constant \\c\\ leaves \\\beta\\ alone and rescales the
baseline by \\e^{-c\beta}\\: centring covariates moves the baseline, not
the coefficients.

[`vignette("data-preparation")`](https://dayusun.github.io/pcdreg/articles/data-preparation.md)
covers the data layouts in full, including the conversion from
cumulative counts, which is the mistake most worth avoiding.

## References

Sun, D., Guo, Y., Li, Y., Tu, W. and Sun, J. (2024). A robust approach
for regression analysis of panel count data. *Bernoulli* **30**(4),
3251–3275. [doi:10.3150/23-BEJ1713](https://doi.org/10.3150/23-BEJ1713)

Murphy, S. A. and van der Vaart, A. W. (2000). On profile likelihood.
*Journal of the American Statistical Association* **95**, 449–465.

## See also

[`pcd()`](https://dayusun.github.io/pcdreg/reference/pcd.md),
[`panelmean()`](https://dayusun.github.io/pcdreg/reference/panelmean.md),
[`vcov.pcdfit()`](https://dayusun.github.io/pcdreg/reference/vcov.pcdfit.md),
[`r_panel_count()`](https://dayusun.github.io/pcdreg/reference/r_panel_count.md)

## Examples

``` r
set.seed(1)
d <- r_panel_count(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
fit <- panelrate(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
fit
#> 
#> Call:
#> panelrate(formula = pcd(id, tstart, tstop, count) ~ x1 + x2, 
#>     data = d)
#> 
#> Coefficients:
#>     x1      x2  
#>  1.181  -1.021  
#> 
#> 80 subjects, 320 examinations, 645 events.
#> Log likelihood -426.1 in 328 EM iterations.
summary(fit)
#> 
#> Call:
#> panelrate(formula = pcd(id, tstart, tstop, count) ~ x1 + x2, 
#>     data = d)
#> 
#> Proportional rate model for panel count data 
#> Standard errors: robust sandwich
#> 
#>    Estimate Std. Error z value Pr(>|z|)
#> x1   1.1811     0.1663   7.102 1.23e-12
#> x2  -1.0206     0.1655  -6.165 7.05e-10
#> 
#> 80 subjects, 320 examinations, 645 events, 156 distinct examination times.
#> Log likelihood -426.1 in 328 EM iterations.
```
