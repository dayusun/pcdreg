#!/bin/bash
#SBATCH --job-name=pcdreg-table1
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=06:00:00
#SBATCH --mem=32G
#SBATCH --output=table1-%j.out
#SBATCH --error=table1-%j.err
##SBATCH --account=REPLACE_WITH_YOUR_SLURM_ACCOUNT
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=dayusun@iu.edu

# Replicate Table 1 of Sun et al. (2024) on Big Red 200.
#
#   sbatch table1.slurm            # the paper's 1000 replications
#   sbatch table1.slurm 100        # a quick shakedown first
#
# Run install.sh once on the login node before submitting this.

set -euo pipefail

module load r

export R_LIBS_USER="${HOME}/R/pcdreg-lib"

# Keep the linear algebra single threaded.  The fitting code is deliberately
# free of large matrix products, so extra BLAS threads only add contention, and
# each of the Slurm-allocated cores is already running its own replication.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

NREP="${1:-1000}"

echo "host        : $(hostname)"
echo "cores       : ${SLURM_CPUS_PER_TASK}"
echo "replications: ${NREP}"
echo "started     : $(date)"

Rscript table1.R "${NREP}" "${SLURM_CPUS_PER_TASK}"

echo "finished    : $(date)"
