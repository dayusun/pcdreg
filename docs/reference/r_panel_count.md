# Simulate panel count data with a time-varying covariate

Generates data from the design used in the simulation study of Sun et
al. (2024), which is a convenient source of examples and test cases.

## Usage

``` r
r_panel_count(
  n,
  beta = c(1, -1),
  lambda = function(t) 8/(1 + t),
  tau = 2,
  frailty = 0,
  exam_mean = 4,
  lambda_max = NULL,
  digits = 2
)
```

## Arguments

- n:

  Number of subjects.

- beta:

  Length two coefficient vector for `x1` and `x2`.

- lambda:

  Baseline rate function, vectorised over `t`.

- tau:

  Maximum follow-up time.

- frailty:

  Variance of a gamma frailty with mean one that multiplies each
  subject's rate. Zero, the default, gives a genuine Poisson process; a
  positive value produces overdispersed counts that violate the Poisson
  assumption while leaving the rate model intact.

- exam_mean:

  Mean of the zero-truncated Poisson distribution from which the number
  of examinations per subject is drawn.

- lambda_max:

  Optional upper bound for `lambda` on `[0, tau]`, used by the thinning
  algorithm. The default takes the maximum over a fine grid, which is
  reliable unless `lambda` oscillates very rapidly.

- digits:

  Number of decimal places to round examination times to, or `NA` to
  leave them unrounded. Rounding follows the design of the paper and of
  Wellner and Zhang (2007): examinations then tie across subjects, which
  is realistic for scheduled visits and keeps the pooled grid, and hence
  the cost of fitting, from growing quadratically in the sample size.

## Value

A data frame in counting process form with columns `id`, `tstart`,
`tstop`, `count`, `x1` and `x2`. Rows are intervals over which the
covariates are constant; `count` is the number of events since the
previous examination and is `NA` on rows that only record a covariate
change.

## Details

The time-varying covariate is a single step, \\x_1(t) = X\_{11} I(t \le
U_1) + X\_{12} I(t \> U_1)\\, with \\X\_{11}, X\_{12} \sim U(0, 1)\\ and
\\U_1 \sim U(0, \tau)\\. The time-invariant covariate `x2` is \\U(0,
1)\\. Events follow a nonhomogeneous Poisson process with rate \\r
\lambda(t) \exp(\beta' X(t))\\, where \\r\\ is one unless `frailty` is
positive. Examination times are the order statistics of a uniform sample
on \\(0, \bar T)\\ with \\\bar T \sim U(0.9\tau, \tau)\\.

## References

Sun, D., Guo, Y., Li, Y., Tu, W. and Sun, J. (2024). A robust approach
for regression analysis of panel count data. *Bernoulli* **30**(4),
3251–3275. [doi:10.3150/23-BEJ1713](https://doi.org/10.3150/23-BEJ1713)

## Examples

``` r
set.seed(42)
d <- r_panel_count(5, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
head(d)
#>   id tstart tstop count        x1        x2
#> 1  1   0.00  0.27     2 0.7050648 0.7191123
#> 2  1   0.27  0.57     0 0.7050648 0.7191123
#> 3  1   0.57  1.03     1 0.7050648 0.7191123
#> 4  1   1.03  1.28     0 0.7050648 0.7191123
#> 5  1   1.28  1.31     0 0.7050648 0.7191123
#> 6  1   1.31  1.46     0 0.7050648 0.7191123

# Overdispersed counts: the rate model still holds, the Poisson one does not.
od <- r_panel_count(5, frailty = 1)
```
