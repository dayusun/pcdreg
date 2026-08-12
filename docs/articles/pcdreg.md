# Regression analysis of panel count data

``` r
library(pcdreg)
```

## The data and the model

Panel count data arise when a recurrent event process is followed by
periodic examinations. At each visit you learn how many events have
happened since the last visit, but not when any of them happened.
Hospitalisations counted at quarterly clinic visits, tumours counted at
scheduled sacrifices, and insurance claims tallied at renewal all have
this shape.

Write $N_{i}(t)$ for the number of events subject $i$ has experienced by
time $t$, observed at times $T_{i1} < \ldots < T_{iJ_{i}}$, and
$\Delta N_{ij}$ for the count in $(T_{i,j - 1},T_{ij}\rbrack$.
[`panelrate()`](https://www.sundayu.me/pcdreg/reference/panelrate.md)
fits the proportional rate model

$$E\left\lbrack dN(t) \mid X(t) \right\rbrack = \exp\left( \beta\prime X(t) \right)\, d\Lambda(t),$$

where $\Lambda$ is an unspecified non-decreasing baseline cumulative
rate and $X(t)$ may vary over time.

It is worth being clear about why the rate model rather than the
proportional *means* model
$E\left\lbrack N(t) \mid X(t) \right\rbrack = \mu(t)\exp\left( \beta\prime X(t) \right)$.
The means model requires $\mu$ to be non-decreasing, and when covariates
fluctuate over time that is hard to arrange and harder to interpret:
nothing stops a fitted means model from predicting a mean count that
goes down. The rate model constrains only
$\exp\left( \beta\prime X(t) \right)\, d\Lambda(t)$ to be positive, so
predicted means increase automatically, and
$\exp\left( \beta\prime X(t) \right)$ carries the same instantaneous
relative-risk reading that a hazard ratio has in survival analysis.

The means model is nonetheless the established alternative, and the
package fits it too, with
[`panelmean()`](https://www.sundayu.me/pcdreg/reference/panelmean.md).
The last section compares the two.

## Getting data into shape

The response is built by
[`pcd()`](https://www.sundayu.me/pcdreg/reference/pcd.md), which plays
the role `Surv()` plays in survival analysis and, like it, has two
forms.

When covariates do not change over time, one row per examination is
enough:

``` r
visits <- data.frame(
  id    = c(1, 1, 1, 2, 2),
  time  = c(0.5, 1.1, 1.9, 0.7, 1.4),
  count = c(2, 0, 3, 1, 4),
  dose  = c(10, 10, 10, 25, 25)
)
pcd(visits$id, visits$time, visits$count)
#> [1] 1: (0.0, 0.5] 2 1: (0.5, 1.1] 0 1: (1.1, 1.9] 3 2: (0.0, 0.7] 1
#> [5] 2: (0.7, 1.4] 4
```

Interval starts are filled in from the previous examination of the same
subject, with follow-up beginning at zero.

When covariates *do* change over time you need the counting process
form, one row per interval over which the covariates are constant. This
is the same layout `coxph()` uses for time-dependent covariates. A row
whose `count` is a number records an examination at `tstop`; a row whose
`count` is `NA` records only that a covariate changed there.

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

Subject 1 was examined at 0.5 and 1.2 with two events in the first
interval and none in the second; the dose changed at 0.9, part way
through the second interval, which is exactly the information the rate
model can use and the means model cannot.

Within a subject the intervals must be contiguous, start at zero, and
end at an examination.
[`panelrate()`](https://www.sundayu.me/pcdreg/reference/panelrate.md)
checks all of this and says which subject is at fault when it fails.

### Converting data held in the older layout

A common alternative layout repeats each panel count on every row of its
interval and marks the examinations with a separate indicator. Turning
that into the layout above takes one line, because a repeated count on a
non-examination row carries no information that the count on the
examination row does not:

``` r
new <- old
new$count[old$event == 0] <- NA
panelrate(pcd(id, tstart, tstop, count) ~ x1 + x2, data = new)
```

## Fitting

``` r
set.seed(1)
d <- r_panel_count(200, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
fit <- panelrate(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
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
#> x1  1.10111    0.10693   10.30   <2e-16
#> x2 -1.11478    0.09385  -11.88   <2e-16
#> 
#> 200 subjects, 819 examinations, 1601 events, 195 distinct examination times.
#> Log likelihood -1106 in 522 EM iterations.
```

There is no intercept: a constant in $X$ would be indistinguishable from
a rescaling of $\Lambda$, so
[`panelrate()`](https://www.sundayu.me/pcdreg/reference/panelrate.md)
drops it and reports the baseline separately.

Fitting maximises the likelihood implied by a nonhomogeneous Poisson
process, using an EM algorithm that treats the unobserved per-time
counts as missing data. That makes the awkward integral
$\int\exp\left( \beta\prime X(t) \right)d\Lambda(t)$ tractable and,
usefully, means the covariates are only ever needed at the observed
examination times rather than as complete trajectories.

The Poisson likelihood is a working device, not an assumption about the
data. The estimator stays consistent and asymptotically normal when the
counts are not Poisson at all.

## Which standard errors

That distinction matters most for the covariance. Three estimators are
available.

``` r
se <- function(type) sqrt(diag(vcov(fit, type)))
rbind(robust = se("robust"), information = se("information"))
#>                     x1         x2
#> robust      0.10693035 0.09384671
#> information 0.07937988 0.08078999
```

`"robust"`, the default, is the sandwich $\Omega^{- 1}S\Omega^{- 1}/n$
built from the two matrices

``` r
fit$Omega
#>             x1          x2
#> x1 -0.60588046  0.06284371
#> x2  0.06284371 -0.67648334
fit$S
#>            x1         x2
#> x1  0.8597886 -0.2345604
#> x2 -0.2345604  0.8300369
```

both of which are by-products of the last EM iteration, so the robust
standard errors are essentially free. It is consistent whether or not
the counts are Poisson.

`"information"` is $S^{- 1}/n$. Under the Poisson assumption $S$ is the
efficient information and this is the right answer; it is also a cheap
stand-in for the profile likelihood method, which it agrees with
closely.

`"profile"` is the numerical profile likelihood estimator of Murphy and
van der Vaart. It is included for comparison with the literature, and
because it is expensive it is computed only on request:

``` r
pfit <- panelrate(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d,
                  profile = TRUE)
rbind(information = sqrt(diag(vcov(pfit, "information"))),
      profile     = sqrt(diag(vcov(pfit, "profile"))))
#>                     x1         x2
#> information 0.07937988 0.08078999
#> profile     0.08034997 0.07755357
```

The two track each other closely, as the theory says they should, since
both estimate $S$. They are not identical: the profile version rests on
a numerical derivative with a step of order $n^{- 1/2}$, so it carries a
discretisation error that shrinks only as the sample grows. That is the
practical argument for preferring `"information"` whenever the Poisson
assumption is tenable, and it is why `"profile"` is here mainly so that
results can be compared with the literature.

### What happens when the Poisson assumption fails

Real recurrent event data are usually overdispersed: some subjects are
simply more event-prone than their covariates suggest. Simulating that
with a gamma frailty leaves the rate model intact but destroys the
Poisson structure.

``` r
set.seed(2)
od <- r_panel_count(200, beta = c(1, -1), frailty = 1)
odfit <- panelrate(pcd(id, tstart, tstop, count) ~ x1 + x2, data = od)
rbind(robust = sqrt(diag(vcov(odfit, "robust"))),
      information = sqrt(diag(vcov(odfit, "information"))))
#>                     x1         x2
#> robust      0.26146233 0.28237814
#> information 0.04003442 0.03288915
```

The information-based standard errors are several times too small. The
simulations in the paper put the resulting coverage of nominal 95%
intervals around 15–22%. This is the reason `"robust"` is the default,
and the reason to think hard before reporting anything else.

## The baseline and predicted means

[`baseline()`](https://www.sundayu.me/pcdreg/reference/baseline.md)
returns the estimated jumps and the cumulative baseline rate.

``` r
head(baseline(fit))
#>   time         jump   cumrate
#> 1 0.01 1.207759e-01 0.1207759
#> 2 0.02 0.000000e+00 0.1207759
#> 3 0.03 0.000000e+00 0.1207759
#> 4 0.04 1.355770e-46 0.1207759
#> 5 0.05 4.659539e-01 0.5867298
#> 6 0.06 0.000000e+00 0.5867298
plot(fit)
lines(fit$baseline$time, 8 * log(1 + fit$baseline$time), col = 2, lty = 2)
```

![Estimated cumulative baseline rate rising with time, with the true
curve overlaid as a dashed red line and closely tracking
it.](pcdreg_files/figure-html/baseline-1.png)

The dashed line is the truth, $\Lambda(t) = 8\log(1 + t)$.

Predicted mean counts follow each subject’s covariate trajectory,
$E\left\lbrack N(t) \mid X \right\rbrack = \int_{0}^{t}\exp\left( \beta\prime X(s) \right)\, d\widehat{\Lambda}(s)$:

``` r
pred <- predict(fit, d, type = "mean")
head(pred)
#>   id time      mean
#> 1  1 0.01 0.1554798
#> 2  1 0.02 0.1554798
#> 3  1 0.03 0.1554798
#> 4  1 0.04 0.1554798
#> 5  1 0.05 0.7553213
#> 6  1 0.06 0.7553213

one <- pred[pred$id == 1, ]
plot(one$time, one$mean, type = "s", xlab = "Time",
     ylab = "Predicted mean count", main = "Subject 1")
```

![Predicted mean count for one subject as a step function, rising
monotonically over follow-up.](pcdreg_files/figure-html/predict-1.png)

The curve is non-decreasing even though `x1` steps up or down part way
through follow-up, which is the structural advantage of the rate model.

## The means model

[`panelmean()`](https://www.sundayu.me/pcdreg/reference/panelmean.md)
fits the other classical model for these data,

$$E\left\lbrack N(t) \mid X(t) \right\rbrack = \mu(t)\exp\left( \beta\prime X(t) \right),$$

by the estimating equation of Hu, Sun and Wei (2003). It is the
comparator used in the application of the paper, and it is cheap: only
the examination times enter, so none of the grid expansion the rate
model needs is required.

``` r
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
#> x1  0.53845    0.09857   5.463 4.69e-08
#> x2 -1.04802    0.10824  -9.682  < 2e-16
#> 
#> 200 subjects, 819 examinations, 1601 events, 195 distinct examination times.
#> Converged in 4 Newton iterations.
```

There is no likelihood behind this estimator, so there is no log
likelihood to report and [`vcov()`](https://rdrr.io/r/stats/vcov.html)
offers only the sandwich.

Two things are worth care when reading the two fits side by side.

``` r
cbind(rate = coef(fit), mean = coef(mfit))
#>         rate       mean
#> x1  1.101105  0.5384467
#> x2 -1.114778 -1.0480179
```

First, these numbers are not two estimates of one quantity. With
time-varying covariates the models are not reparametrisations of each
other: $\beta$ acts on the instantaneous rate in one and on the
cumulative mean in the other, so they answer different questions and
there is no reason for them to agree. The data here were generated from
the rate model, which is why
[`panelrate()`](https://www.sundayu.me/pcdreg/reference/panelrate.md)
recovers the values used to simulate them and
[`panelmean()`](https://www.sundayu.me/pcdreg/reference/panelmean.md)
does not.

Second, nothing constrains the fitted $\widehat{\mu}$ to increase:

``` r
mu <- baseline(mfit)$mean
c(increasing = !is.unsorted(mu), decreases_at = sum(diff(mu) < 0))
#>   increasing decreases_at 
#>            0           96
```

That is the difficulty with the means model this package’s rate model is
meant to avoid, and it is visible here on ordinary simulated data rather
than only in principle. A mean function that falls is not a numerical
failure; it is what the estimator returns when covariate values move up
and down, and it is why predicted means from
[`panelmean()`](https://www.sundayu.me/pcdreg/reference/panelmean.md)
can decrease while those from
[`panelrate()`](https://www.sundayu.me/pcdreg/reference/panelrate.md)
cannot.

## Computation

The EM algorithm converges linearly and, on this problem, slowly: plain
EM can still be moving in the fifth decimal place after tens of
thousands of passes.
[`panelrate()`](https://www.sundayu.me/pcdreg/reference/panelrate.md)
therefore extrapolates the iterations using the SQUAREM scheme of
Varadhan and Roland, keeping an extrapolated point only when a
stabilising EM pass from it does at least as well on the observed data
log likelihood. The fixed point is unchanged; only the number of passes
needed to reach it falls, typically by an order of magnitude.

``` r
c(iterations = fit$iterations, EM_passes = fit$passes,
  converged = fit$converged)
#> iterations  EM_passes  converged 
#>        522       1567          1
```

Set `accelerate = FALSE` in
[`panelrate_control()`](https://www.sundayu.me/pcdreg/reference/panelrate_control.md)
to recover plain EM.

The dominant cost is the size of the pooled grid of distinct examination
times, because the covariate trajectory of every subject has to be
evaluated at every grid time within their follow-up. Examination times
that tie across subjects, which is what scheduled visits produce in
practice, keep that grid small.

## References

Sun, D., Guo, Y., Li, Y., Tu, W. and Sun, J. (2024). A robust approach
for regression analysis of panel count data. *Bernoulli* **30**(4),
3251–3275. [doi:10.3150/23-BEJ1713](https://doi.org/10.3150/23-BEJ1713)

Hu, X. J., Sun, J. and Wei, L. J. (2003). Regression parameter
estimation from panel counts. *Scandinavian Journal of Statistics*
**30**(1), 25–43.
[doi:10.1111/1467-9469.00316](https://doi.org/10.1111/1467-9469.00316)

Varadhan, R. and Roland, C. (2008). Simple and globally convergent
methods for accelerating the convergence of any EM algorithm.
*Scandinavian Journal of Statistics* **35**, 335–353.
