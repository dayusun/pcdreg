#!/bin/bash
# Run this once on a Big Red 200 login node, which has outbound network access.
# Compute nodes generally do not, so everything must be installed beforehand.
#
#   bash install.sh

set -euo pipefail

module load r

export R_LIBS_USER="${HOME}/R/pcdreg-lib"
mkdir -p "${R_LIBS_USER}"

Rscript -e '
  lib <- Sys.getenv("R_LIBS_USER")
  .libPaths(c(lib, .libPaths()))
  repos <- "https://cloud.r-project.org"
  need <- setdiff(c("Rcpp", "RcppArmadillo", "remotes"), rownames(installed.packages()))
  if (length(need)) install.packages(need, lib = lib, repos = repos)
  remotes::install_github("dayusun/pcdreg", lib = lib, upgrade = "never",
                          build_vignettes = FALSE)
  library(pcdreg, lib.loc = lib)
  cat("pcdreg", as.character(packageVersion("pcdreg")), "installed into", lib, "\n")
  set.seed(1)
  d <- r_panel_count(30)
  fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
  cat("smoke test coefficients:", format(coef(fit), digits = 4), "\n")
'
