#!/bin/bash
#--------------------------------------------------------------------------------
#  Replicate Table 1 of Sun et al. (2024), one array task per cell.
#
#  The table's six cells (n = 100, 200, 400, each with and without the Poisson
#  assumption) are independent, so they run on six nodes at once and each uses
#  all of its node's cores for the replications inside it.
#
#    sbatch table1.sbatch          # the paper's 1000 replications
#    sbatch table1.sbatch 50       # shakedown first
#--------------------------------------------------------------------------------
#SBATCH --job-name=pcdreg-table1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time 08:00:00
#SBATCH --mail-user=dayusun@iu.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/home/dsun30/pcdreg/hpc/emory/logs/table1.%A_%a.out
#SBATCH --partition guo-mh,week-long-cpu
#SBATCH --array=1-6
#--------------------------------------------------------------------------------

echo "### Starting cell $SLURM_ARRAY_TASK_ID at: $(date) ###"

# R is picked up directly rather than through lmod: `module` is a shell
# function absent from non-interactive shells, and its MODULEPATH is not set
# when lmod's init is sourced by hand.  This mirrors the newest jobs on this
# cluster.
R_VERSION=4.5.2
export PATH=/apps/R/${R_VERSION}/bin:$PATH
export LD_LIBRARY_PATH=/apps/R/${R_VERSION}/lib64/R/lib:${LD_LIBRARY_PATH:-}
export R_LIBS_USER=${HOME}/R/pcdreg-${R_VERSION}

# One replication per core; the fitting code has no large matrix products, so
# extra BLAS threads on top would only contend.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

NREP="${1:-1000}"

cd /home/dsun30/pcdreg/hpc/emory
mkdir -p logs results

echo "### cell=$SLURM_ARRAY_TASK_ID cores=$SLURM_CPUS_PER_TASK nrep=$NREP node=$(hostname) ###"

srun Rscript table1_cell.R "$SLURM_ARRAY_TASK_ID" "$SLURM_CPUS_PER_TASK" "$NREP"

echo "### Ending at: $(date) ###"
