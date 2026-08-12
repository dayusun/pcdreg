# Reproducing Table 1 on an HPC cluster

Scripts to replicate Table 1 of Sun, Guo, Li, Tu and Sun (2024) away from a
workstation. Written for Big Red 200 at Indiana University, which runs Slurm;
the only cluster-specific parts are the `module load r` line and the `#SBATCH`
header, so adapting them elsewhere is a small edit.

These files are excluded from the package build and are not part of `pcdreg`.

## Three commands

```bash
# on a login node, once
git clone https://github.com/dayusun/pcdreg.git
cd pcdreg/hpc
bash install.sh

# then submit
sbatch table1.slurm 100      # shakedown: a few minutes, checks it all works
sbatch table1.slurm          # the real thing: 1000 replications
```

Before the first submission, open `table1.slurm` and uncomment the `--account`
line with your Slurm account, if your allocation requires one. Adjust
`--partition` if `general` is not the right queue for you.

## What it does

Six cells: sample sizes 100, 200 and 400, each with and without the Poisson
assumption, the latter a gamma frailty of variance one. Every replication fits
the rate model and computes all three covariance estimators, the profile
likelihood one included, which is the expensive part.

Results are written to `table1.rds` and `table1.csv` in the submission
directory, and the job's `.out` file carries both the table and the paper's own
numbers underneath it for comparison.

## Sizing the job

The header asks for 32 cores for 6 hours, which is generous. As a rough guide
from a 24-core workstation, one replication with the profile likelihood takes
about a second at n = 100, four at n = 200 and fifteen at n = 400, so the whole
study is on the order of 6 CPU-hours and finishes in well under an hour on 32
cores. Run the 100-replication shakedown first and scale from its timings.

Each core runs its own replication, and the fitting code contains no large
matrix products, so the batch script pins the linear algebra to a single thread.
Leaving BLAS free to spawn its own threads on top of the 32 processes only
creates contention.

## Monte Carlo error

At 1000 replications a coverage estimate carries a standard error of about 0.7
percentage points and an empirical standard deviation about 2%, which is the
resolution to hold the comparison with the paper to. Differences smaller than
that are noise, not disagreement.
