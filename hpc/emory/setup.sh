#!/bin/bash
# Run once on the Emory login node.  Builds a personal library and installs
# pcdreg into it.  Compute nodes share the filesystem, so this only needs doing
# here.
#
#   bash setup.sh

set -euo pipefail

# R is picked up directly rather than through lmod: `module` is a shell
# function absent from non-interactive shells, and its MODULEPATH is not set
# when lmod's init is sourced by hand.  This mirrors the newest jobs on this
# cluster.
R_VERSION=4.5.2
export PATH=/apps/R/${R_VERSION}/bin:$PATH
export LD_LIBRARY_PATH=/apps/R/${R_VERSION}/lib64/R/lib:${LD_LIBRARY_PATH:-}
export R_LIBS_USER=${HOME}/R/pcdreg-${R_VERSION}
mkdir -p "${R_LIBS_USER}"

Rscript --vanilla -e '
  lib <- Sys.getenv("R_LIBS_USER")
  .libPaths(c(lib, .libPaths()))
  repos <- "https://cloud.r-project.org"
  need <- setdiff(c("Rcpp", "RcppArmadillo", "remotes"),
                  rownames(installed.packages()))
  if (length(need)) {
    cat("installing:", paste(need, collapse = ", "), "\n")
    install.packages(need, lib = lib, repos = repos, quiet = TRUE)
  }
  remotes::install_github("dayusun/pcdreg", lib = lib, upgrade = "never",
                          build_vignettes = FALSE, quiet = TRUE)
  library(pcdreg, lib.loc = lib)
  cat("pcdreg", as.character(packageVersion("pcdreg")), "->", lib, "\n")
  set.seed(1)
  d <- sim_pcd(30)
  fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
  cat("smoke test:", format(coef(fit), digits = 4), "\n")
'
