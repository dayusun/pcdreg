#!/bin/bash
# Run once on the Emory login node.  Builds a personal library and installs
# pcdreg into it.  Compute nodes share the filesystem, so this only needs doing
# here.
#
#   bash setup.sh

set -euo pipefail

# `module` is a shell function from lmod and is not defined in a
# non-interactive shell, so make sure it exists before using it.
if ! command -v module >/dev/null 2>&1; then
  for init in /apps/lmod/lmod/init/bash /usr/share/lmod/lmod/init/bash               /etc/profile.d/modules.sh /etc/profile.d/lmod.sh; do
    [ -r "$init" ] && . "$init" && break
  done
fi

module load R/4.4.0

export R_LIBS_USER="${HOME}/R/pcdreg-4.4.0"
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
  d <- r_panel_count(30)
  fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
  cat("smoke test:", format(coef(fit), digits = 4), "\n")
'
