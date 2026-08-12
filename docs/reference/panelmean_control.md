# Tuning parameters for the means model

Tuning parameters for the means model

## Usage

``` r
panelmean_control(maxit = 100L, reltol = 1e-09)
```

## Arguments

- maxit:

  Maximum number of Newton iterations.

- reltol:

  Convergence tolerance on the relative change in the coefficients.

## Value

A list of control values.

## Examples

``` r
panelmean_control(reltol = 1e-10)
#> $maxit
#> [1] 100
#> 
#> $reltol
#> [1] 1e-10
#> 
```
