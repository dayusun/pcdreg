# Tuning parameters for the fitting algorithms

Tuning parameters for the fitting algorithms

## Usage

``` r
pcdreg_control(
  maxit = 2000L,
  reltol = 1e-07,
  accelerate = TRUE,
  profile_maxit = 20000L,
  profile_reltol = 1e-06,
  profile_h = NULL
)
```

## Arguments

- maxit:

  Maximum number of iterations. For the rate model these are EM
  iterations; for the means model, Newton iterations of the estimating
  equation, which converges in a handful.

- reltol:

  Convergence tolerance on the largest relative change in the
  coefficients, and for the rate model in the baseline jump sizes too.

- accelerate:

  Rate model only: whether to extrapolate the EM iterations by the
  SQUAREM scheme of Varadhan and Roland (2008). This reaches the same
  fixed point in far fewer passes; set it to `FALSE` to run plain EM.

- profile_maxit, profile_reltol:

  Limits for the inner runs that maximise over the baseline with the
  coefficients held fixed, used only by the profile likelihood
  covariance. These can be looser than `reltol` without loss: the
  profile log likelihood is stationary in the baseline at its maximum,
  so an error of size \\\epsilon\\ in the inner solution moves it by
  only \\O(\epsilon^2)\\.

- profile_h:

  Step size for the numerical derivative of the profile log likelihood.
  The default, `NULL`, uses \\5 n^{-1/2}\\ as recommended by Zeng et
  al. (2016) and used in the paper.

## Value

A list of control values.

## References

Varadhan, R. and Roland, C. (2008). Simple and globally convergent
methods for accelerating the convergence of any EM algorithm.
*Scandinavian Journal of Statistics* **35**, 335–353.

## Examples

``` r
pcdreg_control(reltol = 1e-9)
#> $maxit
#> [1] 2000
#> 
#> $reltol
#> [1] 1e-09
#> 
#> $accelerate
#> [1] TRUE
#> 
#> $profile_maxit
#> [1] 20000
#> 
#> $profile_reltol
#> [1] 1e-06
#> 
#> $profile_h
#> NULL
#> 
```
