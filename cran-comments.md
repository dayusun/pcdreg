## Test environments

* Windows 11, R 4.5.1 (local), Rtools45
* Ubuntu 24.04, R 4.5.1 (Docker, rocker/r-ver:4.5.1)

## R CMD check results

0 errors | 0 warnings | 1 note

The note is the usual "New submission".

Two further notes appear only in the Linux container used for testing and are
properties of that environment rather than of the package: the Debian build of R
adds `-Wdate-time -Werror=format-security -Wformat` to its own compilation
flags (the package's own `src/Makevars` passes only
`$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)`, and "checking compilation flags in
Makevars" passes), and the optional V8 package needed to check math rendering in
the HTML manual is not installed there.

## Notes for the reviewer

This is a new submission. It implements two established procedures for panel
count data, both cited in the DESCRIPTION and in `inst/CITATION`:

* `pcdreg()` follows Sun, Guo, Li, Tu and Sun (2024), *Bernoulli* **30**(4),
  3251-3275, <doi:10.3150/23-BEJ1713>.
* `pcdreg(model = "mean")` follows Hu, Sun and Wei (2003), *Scandinavian
  Journal of Statistics* **30**(1), 25-43, <doi:10.1111/1467-9469.00316>.

The package compiles against RcppArmadillo and does not use OpenMP.

The examples run on data simulated by `sim_pcd()`. No dataset is shipped: the
data behind the published application is governed by a data use agreement and
is not ours to redistribute.
