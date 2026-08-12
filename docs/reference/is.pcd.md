# Test for a panel count response

Test for a panel count response

## Usage

``` r
is.pcd(x)
```

## Arguments

- x:

  An object.

## Value

`TRUE` if `x` was created by
[`pcd()`](https://dayusun.github.io/pcdreg/reference/pcd.md).

## Examples

``` r
is.pcd(pcd(c(1, 2), c(1, 1), c(0, 2)))
#> [1] TRUE
```
