# Fit the proportional means model to panel count data

Fits \$\$E\[N(t) \mid X(t)\] = \mu(t) \exp(\beta' X(t))\$\$ by the
estimating equation of Hu, Sun and Wei (2003). This is the comparator
the paper behind
[`panelrate()`](https://www.sundayu.me/pcdreg/reference/panelrate.md)
uses, and it is provided so the two model families can be compared on
the same data.

## Usage

``` r
panelmean(
  formula,
  data,
  subset,
  na.action,
  control = panelmean_control(),
  init = NULL
)
```

## Arguments

- formula:

  A formula whose left hand side is a
  [`pcd()`](https://www.sundayu.me/pcdreg/reference/pcd.md) object.

- data:

  A data frame containing the variables in `formula`.

- subset, na.action:

  Handled as in [`stats::lm()`](https://rdrr.io/r/stats/lm.html).

- control:

  A list from
  [`panelmean_control()`](https://www.sundayu.me/pcdreg/reference/panelmean_control.md).

- init:

  Optional starting values for the coefficients.

## Value

An object of class `"panelmean"`, inheriting from `"pcdfit"`, with
components including `coefficients`, `baseline` (a data frame of the
estimated mean function), `Omega` and `S`, and `vcov`.

## Details

The estimating equation is \$\$U(\beta) = \sum_i \sum_j N_i(T\_{ij}) \\
X_i(T\_{ij}) - \bar x(T\_{ij}) \\,\$\$ where \\\bar x(t)\\ averages the
covariates over the subjects examined at \\t\\, weighted by
\\\exp(\beta' X)\\, and \\N_i(T\_{ij})\\ is the cumulative count. Only
the examination times enter, so no covariate values between them are
needed and the fit is much cheaper than
[`panelrate()`](https://www.sundayu.me/pcdreg/reference/panelrate.md).

There is no likelihood here, so no log likelihood is reported and the
only covariance available is the sandwich \\\Omega^{-1} S \Omega^{-1} /
n\\.

Two cautions when comparing with
[`panelrate()`](https://www.sundayu.me/pcdreg/reference/panelrate.md).
The models are not reparametrisations of each other when covariates vary
over time, so their coefficients answer different questions: \\\beta\\
here acts on the cumulative mean, and in the rate model on the
instantaneous rate. And nothing constrains the fitted \\\hat\mu\\ to
increase; when covariate values fluctuate it may well not, which is the
difficulty with the means model that motivates the rate model.

## Specifying the model

Identical to
[`panelrate()`](https://www.sundayu.me/pcdreg/reference/panelrate.md):
the left hand side is a
[`pcd()`](https://www.sundayu.me/pcdreg/reference/pcd.md) object and the
right hand side an ordinary model formula, with no intercept. The two
fitters accept exactly the same data, so a model can be refitted under
the other by changing the function name alone. See
[`vignette("data-preparation")`](https://www.sundayu.me/pcdreg/articles/data-preparation.md)
for the data layouts.

## References

Hu, X. J., Sun, J. and Wei, L. J. (2003). Regression parameter
estimation from panel counts. *Scandinavian Journal of Statistics*
**30**(1), 25–43.
[doi:10.1111/1467-9469.00316](https://doi.org/10.1111/1467-9469.00316)

Sun, D., Guo, Y., Li, Y., Tu, W. and Sun, J. (2024). A robust approach
for regression analysis of panel count data. *Bernoulli* **30**(4),
3251–3275. [doi:10.3150/23-BEJ1713](https://doi.org/10.3150/23-BEJ1713)

## See also

[`panelrate()`](https://www.sundayu.me/pcdreg/reference/panelrate.md),
[`pcd()`](https://www.sundayu.me/pcdreg/reference/pcd.md)

## Examples

``` r
set.seed(1)
d <- r_panel_count(80, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
mfit <- panelmean(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
summary(mfit)
#> 
#> Call:
#> panelmean(formula = pcd(id, tstart, tstop, count) ~ x1 + x2, 
#>     data = d)
#> 
#> Proportional means model for panel count data 
#> Standard errors: robust sandwich
#> 
#>    Estimate Std. Error z value Pr(>|z|)
#> x1   0.8754     0.1475   5.934 2.96e-09
#> x2  -1.2213     0.1671  -7.309 2.69e-13
#> 
#> 80 subjects, 320 examinations, 645 events, 156 distinct examination times.
#> Converged in 5 Newton iterations.

# The rate model on the same data. The coefficients are not comparable
# term by term, because the two models act on different quantities.
rfit <- panelrate(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
cbind(mean = coef(mfit), rate = coef(rfit))
#>          mean      rate
#> x1  0.8753684  1.181091
#> x2 -1.2212558 -1.020581
```
