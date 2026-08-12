# Where pcdreg fits among R packages for panel count data

``` r
library(pcdreg)
```

This article describes the settings in which `pcdreg` and the available
alternatives apply.

## The landscape

Three classes of R package address recurrent event data. They are not
interchangeable.

**Exactly observed event times.** If you know *when* each event
happened, `reReg`, `survival`, and the recurrent-event literature apply.
The discussion below concerns panel count data, where only the number of
events between visits is known.

**`spef`**, by Chiou, Wang and Yan, is the established R package for
panel count data. It offers a wide range of estimators through one
function, `panelReg(formula, data, method = ...)`. These include
augmented estimating equations for an informative observation process
(`AEE`, `AEEX`), maximum pseudolikelihood and monotone spline likelihood
methods (`MPL`, `MPLs`, `MLs`), the Sun–Wei and Hu–Sun–Wei estimating
equations, and an accelerated mean model. Its response is
`PanelSurv(ID, time, count)`, one row per examination. All of these
methods target the **mean** model. At the time of writing, `spef` is
archived on CRAN — since June 2026, for an undeliverable maintainer
address rather than anything to do with the code — so it installs from
the archive rather than from
[`install.packages()`](https://rdrr.io/r/utils/install.packages.html).

**`pcdreg`** fits the rate and mean models with time-varying covariates.
Its covariance estimate does not rely on a Poisson assumption.

## What is actually different

### The rate model

This is the substantive distinction. The means model constrains the
*cumulative* count; the rate model constrains the *instantaneous* count.
With time-varying covariates, these are different models. They are not
two parametrisations of one, and has a different meaning in each.

Simulate from the rate model with a covariate that steps part way
through follow-up, and fit both:

``` r
set.seed(1)
d <- sim_pcd(300, beta = c(1, -1), lambda = function(t) 8 / (1 + t))

rate <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
mean <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d,
               model = "mean")

round(cbind(truth = c(1, -1), rate = coef(rate), mean = coef(mean)), 3)
#>    truth   rate   mean
#> x1     1  1.004  0.517
#> x2    -1 -1.074 -1.036
```

The rate model recovers the values used to generate the data. The means
model does not, and it is not meant to: it is estimating a different
quantity. If your scientific question is about the instantaneous risk —
the analogue of a hazard ratio — then only the first column of estimates
answers it, and a package offering mean models alone cannot produce it.

If the question concerns cumulative burden by time , the means model is
appropriate. `pcdreg(model = "mean")` and `spef`’s `EE.HSWc` are the
same estimator.

### Monotonicity

There is a practical consequence. Nothing constrains a fitted mean
function to increase, and with a fluctuating covariate it generally does
not:

``` r
mu <- baseline(mean)$mean
c(increasing = !is.unsorted(mu), decreases_at = sum(diff(mu) < 0))
#>   increasing decreases_at 
#>            0           92
```

A mean function that falls is not a numerical failure; it is what the
estimator returns, and it makes the fitted model awkward to interpret
and to predict from. The rate model has no such problem, because it
constrains only to be positive, so predicted means increase by
construction:

``` r
pr <- predict(rate, d, type = "mean")
all(vapply(split(pr$mean, pr$id), function(v) all(diff(v) >= 0), logical(1)))
#> [1] TRUE
```

Some `spef` methods fit with monotone I-splines, so the *baseline*
increases by construction. This resolves half the problem. The fitted
mean for a subject can still decline when their covariates decline.

### Covariates between examinations

The rate-model likelihood evaluates each subject’s covariate trajectory
at every pooled examination time during follow-up, including times at
which *that* subject was not examined. A `PanelSurv(ID, time, count)`
response has one row per examination and cannot record a covariate value
when the subject was not seen. The counting-process form of
[`pcd()`](https://www.sundayu.me/pcdreg/reference/pcd.md) can:

``` r
head(subset(d, id == 1), 4)
#> # A tibble: 4 × 6
#>      id tstart tstop count    x1    x2
#>   <int>  <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1     1   0     0.38     1 0.898 0.661
#> 2     1   0.38  1.07     3 0.898 0.661
#> 3     1   1.07  1.26    NA 0.898 0.661
#> 4     1   1.26  1.7      4 0.945 0.661
```

Row 3 has `count = NA`, recording a covariate change at that time, not
an examination. The rate model needs this information. The
one-row-per-examination layout cannot carry it.

### Variance without the Poisson assumption

The rate-model likelihood is the likelihood of a Poisson process. It is
a working device, not a claim about the data. Recurrent event data are
usually overdispersed, so standard errors that assume Poisson counts are
far too small. `pcdreg` reports a robust sandwich by default:

``` r
set.seed(2)
od <- sim_pcd(300, beta = c(1, -1), frailty = 1)  # gamma frailty
odfit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = od)

round(rbind(robust      = sqrt(diag(vcov(odfit, "robust"))),
            information = sqrt(diag(vcov(odfit, "information")))), 4)
#>                 x1     x2
#> robust      0.2017 0.2152
#> information 0.0335 0.0295
```

The paper’s simulations put coverage of nominal 95% intervals built from
the second row at 15–22%. Both estimators are available; the robust
estimator is the default.

## Seeing the data first

`pcdreg` also plots the observation pattern:

``` r
plot(with(d, pcd(id, tstart, tstop, count)), max_subjects = 40)
#> Showing 40 of 300 subjects.
```

![Panel count data plot: one row per subject, follow-up drawn as a grey
line, filled circles at examinations with area proportional to the event
count and open circles where no events were
recorded.](comparison_files/figure-html/dataplot-1.png)

Two features of this plot affect the analysis. Alignment of examination
times across subjects determines the pooled-grid size and the cost of
fitting the rate model. The raggedness of the right edge shows how much
of the baseline tail rests on a handful of subjects.

## Which to use

| If you need                                                   | Use                               |
|---------------------------------------------------------------|-----------------------------------|
| the rate model, or covariates that change between visits      | `pcdreg`                          |
| standard errors that survive overdispersion                   | `pcdreg`                          |
| the mean model with the Hu–Sun–Wei estimating equation        | either                            |
| informative observation or censoring times                    | `spef` (`AEE`, `AEEX`)            |
| a monotone spline baseline mean                               | `spef` (`MPLs`, `MLs`)            |
| an accelerated mean model, or the other mean-model estimators | `spef`                            |
| exactly observed event times                                  | not panel count data; see `reReg` |

`pcdreg` is deliberately narrow. It implements two models with
covariance estimation worked out and tested. Use `spef` when it provides
an estimator you need.

## Reference

Sun, D., Guo, Y., Li, Y., Tu, W. and Sun, J. (2024). A robust approach
for regression analysis of panel count data. *Bernoulli* **30**(4),
3251–3275. [doi:10.3150/23-BEJ1713](https://doi.org/10.3150/23-BEJ1713)
