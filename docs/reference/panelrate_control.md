# Tuning parameters for the EM algorithm

Tuning parameters for the EM algorithm

## Usage

``` r
panelrate_control(
  maxit = 2000L,
  reltol = 1e-07,
  profile_maxit = 20000L,
  profile_reltol = 1e-06,
  profile_h = NULL,
  accelerate = TRUE
)
```

## Arguments

- maxit:

  Maximum number of EM iterations.

- reltol:

  Convergence tolerance. Iteration stops once the largest relative
  change in the coefficients and in the baseline jump sizes both fall
  below this value.

- profile_maxit, profile_reltol:

  Corresponding limits for the inner runs that maximise over the
  baseline with the coefficients held fixed, used only by the profile
  likelihood covariance. These can be looser than `reltol` without loss:
  the profile log likelihood is stationary in the baseline at its
  maximum, so an error of size \\\epsilon\\ in the inner solution moves
  it by only \\O(\epsilon^2)\\.

- profile_h:

  Step size for the numerical derivative of the profile log likelihood.
  The default, `NULL`, uses \\5 n^{-1/2}\\ as recommended by Zeng et
  al. (2016) and used in the paper.

- accelerate:

  Whether to extrapolate the EM iterations by the SQUAREM scheme of
  Varadhan and Roland (2008). This reaches the same fixed point in far
  fewer passes; set it to `FALSE` to run plain EM.

## Value

A list of control values.

## References

Varadhan, R. and Roland, C. (2008). Simple and globally convergent
methods for accelerating the convergence of any EM algorithm.
*Scandinavian Journal of Statistics* **35**, 335–353.

## Examples

``` r
panelrate_control(reltol = 1e-9)
#> $maxit
#> [1] 2000
#> 
#> $reltol
#> [1] 1e-09
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
#> $accelerate
#> [1] TRUE
#> 
```
