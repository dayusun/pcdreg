# Semiparametric regression for panel count data

Fits either of the two standard semiparametric models to panel count
data, that is, to recurrent events observed only at intermittent
examination times. Covariates may vary over time.

## Usage

``` r
pcdreg(
  formula,
  data,
  model = c("rate", "mean"),
  subset,
  na.action,
  control = pcdreg_control(),
  profile = FALSE,
  init = NULL
)
```

## Arguments

- formula:

  A formula whose left hand side is a
  [`pcd()`](https://www.sundayu.me/pcdreg/reference/pcd.md) object, for
  example `pcd(id, tstart, tstop, count) ~ x1 + x2`.

- data:

  A data frame containing the variables in `formula`.

- model:

  Which model to fit. `"rate"` is the proportional rate model and the
  default; `"mean"` is the proportional means model.

- subset, na.action:

  Handled as in [`stats::lm()`](https://rdrr.io/r/stats/lm.html). Note
  that missing covariate values delete rows, which can break the
  contiguity of a subject's intervals; resolve them before fitting.

- control:

  A list from
  [`pcdreg_control()`](https://www.sundayu.me/pcdreg/reference/pcdreg_control.md).

- profile:

  Rate model only: whether to also compute the profile likelihood
  covariance. It is provided for comparison with the literature and
  costs one extra baseline-only fit per coefficient, so it is off by
  default.

- init:

  Optional starting values for the coefficients.

## Value

An object of class `"pcdfit"`, with components including `coefficients`,
`model` (which model was fitted), `baseline` (a data frame of the
estimated baseline function), `Omega` and `S` (the two matrices behind
the sandwich), `vcov` (a list of the available covariance estimates),
and convergence information. Rate model fits also carry `loglik`.

## The two models

\$\$\textrm{rate:} \quad E\[dN(t) \mid X(t)\] = \exp(\beta' X(t)) \\
d\Lambda(t)\$\$ \$\$\textrm{mean:} \quad E\[N(t) \mid X(t)\] = \mu(t)
\exp(\beta' X(t))\$\$

The **rate model** is fitted by nonparametric maximum likelihood,
treating \\\Lambda\\ as a step function with a jump at each observed
examination time. The likelihood maximised is the one a nonhomogeneous
Poisson process would give, but that is a working device rather than an
assumption about the data: the estimator stays consistent and
asymptotically normal when it fails, and the default covariance
estimator stays valid.

The **means model** is fitted by the estimating equation of Hu, Sun and
Wei (2003). Only the examination times enter, so no covariate values
between them are needed and it is much cheaper. There is no likelihood
behind it, so it reports no log likelihood, `profile = TRUE` is
unavailable, and the only covariance is the sandwich.

The two are not reparametrisations of each other when covariates vary
over time, so their coefficients answer different questions rather than
estimating a common quantity: \\\beta\\ acts on the instantaneous rate
in one and on the cumulative mean in the other. A further practical
difference is that nothing constrains the fitted \\\hat\mu\\ to
increase, and with fluctuating covariates it generally does not, which
is the drawback of the means model that motivates the rate model.

## Specifying the model

The left hand side is always a
[`pcd()`](https://www.sundayu.me/pcdreg/reference/pcd.md) object, which
carries the subject, the interval and the event count. The right hand
side is an ordinary model formula, so interactions, transformations,
factors and `.` all behave as usual:

    pcd(id, tstart, tstop, count) ~ x1 + x2
    pcd(id, tstart, tstop, count) ~ x1 * x2 + log(x3)
    pcd(id, time, count)          ~ .

The dot expands to every column that is not part of the response, so the
identifier, times and counts are excluded automatically.

No intercept is fitted. A constant column cannot be distinguished from a
rescaling of the baseline, so it is dropped and the baseline is reported
separately by
[`baseline()`](https://www.sundayu.me/pcdreg/reference/baseline.md).
Adding `- 1` changes nothing. One consequence is that shifting a
covariate by a constant \\c\\ leaves \\\beta\\ alone and rescales the
baseline by \\e^{-c\beta}\\: centring covariates moves the baseline, not
the coefficients.

[`vignette("data-preparation")`](https://www.sundayu.me/pcdreg/articles/data-preparation.md)
covers the data layouts in full, including the conversion from
cumulative counts, which is the mistake most worth avoiding.

## Covariance estimation

Available through
[`vcov.pcdfit()`](https://www.sundayu.me/pcdreg/reference/vcov.pcdfit.md).
`"robust"` is the sandwich \\\Omega^{-1} S \Omega^{-1} / n\\, which does
not rely on the Poisson assumption, is the default, and is the only
option for the means model.

For the rate model two more are available. `"information"` is
\\S^{-1}/n\\, valid only under the Poisson assumption but a fast
substitute for the profile likelihood when it holds. `"profile"` is the
numerical profile likelihood estimator of Murphy and van der Vaart
(2000), computed only when `profile = TRUE`; it is there for comparison
with the literature rather than for routine use, since it rests on a
numerical derivative with a step of order \\n^{-1/2}\\ that is coarse at
small sample sizes.

When the counts are overdispersed relative to Poisson, which is the
usual situation, the latter two understate the standard errors,
sometimes severely.

## References

Sun, D., Guo, Y., Li, Y., Tu, W. and Sun, J. (2024). A robust approach
for regression analysis of panel count data. *Bernoulli* **30**(4),
3251–3275. [doi:10.3150/23-BEJ1713](https://doi.org/10.3150/23-BEJ1713)

Hu, X. J., Sun, J. and Wei, L. J. (2003). Regression parameter
estimation from panel counts. *Scandinavian Journal of Statistics*
**30**(1), 25–43.
[doi:10.1111/1467-9469.00316](https://doi.org/10.1111/1467-9469.00316)

Murphy, S. A. and van der Vaart, A. W. (2000). On profile likelihood.
*Journal of the American Statistical Association* **95**, 449–465.

## See also

[`pcd()`](https://www.sundayu.me/pcdreg/reference/pcd.md),
[`vcov.pcdfit()`](https://www.sundayu.me/pcdreg/reference/vcov.pcdfit.md),
[`baseline()`](https://www.sundayu.me/pcdreg/reference/baseline.md),
[`sim_pcd()`](https://www.sundayu.me/pcdreg/reference/sim_pcd.md)

## Examples

``` r
set.seed(1)
d <- sim_pcd(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))

fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
summary(fit)
#> 
#> Call:
#> pcdreg(formula = pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
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

# The same data under the means model. The coefficients are not comparable
# term by term, because the two models act on different quantities.
mfit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d,
               model = "mean")
cbind(rate = coef(fit), mean = coef(mfit))
#>         rate       mean
#> x1  1.181091  0.8753684
#> x2 -1.020581 -1.2212558
```
