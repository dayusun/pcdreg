#!/bin/bash
# Combine the finished cells into Table 1.
#
#   bash collect.sh
set -euo pipefail
cd "$(dirname "$0")"

R_VERSION=4.5.2
export PATH=/apps/R/${R_VERSION}/bin:$PATH
export LD_LIBRARY_PATH=/apps/R/${R_VERSION}/lib64/R/lib:${LD_LIBRARY_PATH:-}
export R_LIBS_USER=${HOME}/R/pcdreg-${R_VERSION}

Rscript collect.R
