# One cell of Table 1 of Sun, Guo, Li, Tu and Sun (2024), Bernoulli 30(4).
#
# The table has six cells: sample sizes 100, 200 and 400, each with and without
# the Poisson assumption.  They are independent, so the Slurm array runs one per
# node and each uses all the cores of its node for the replications within it.
#
#   Rscript table1_cell.R <cell 1-6> <cores> [nrep]

suppressMessages({
  library(parallel)
  library(pcdreg)
})

args  <- commandArgs(trailingOnly = TRUE)
CELL  <- as.integer(args[1])
NCORE <- if (length(args) >= 2) as.integer(args[2]) else 4L
NREP  <- if (length(args) >= 3) as.integer(args[3]) else 1000L

grid  <- expand.grid(frailty = c(0, 1), n = c(100L, 200L, 400L))
stopifnot(CELL >= 1L, CELL <= nrow(grid))
N       <- grid$n[CELL]
FRAILTY <- grid$frailty[CELL]

TRUTH <- c(x1 = 1, x2 = -1)
Z     <- stats::qnorm(0.975)

cat(sprintf("cell %d | n = %d | poisson = %s | %d reps | %d cores | pcdreg %s\n",
            CELL, N, if (FRAILTY == 0) "yes" else "no", NREP, NCORE,
            as.character(utils::packageVersion("pcdreg"))))

one_rep <- function(i, n, frailty) {
  d <- pcdreg::r_panel_count(n, beta = c(1, -1),
                             lambda = function(t) 8 / (1 + t),
                             frailty = frailty)
  fit <- try(suppressWarnings(
    pcdreg::pcdreg(pcdreg::pcd(id, tstart, tstop, count) ~ x1 + x2, data = d,
                   profile = TRUE,
                   control = pcdreg::pcdreg_control(reltol = 1e-6))
  ), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, 8L))
  se <- function(ty) {
    v <- try(sqrt(diag(stats::vcov(fit, ty))), silent = TRUE)
    if (inherits(v, "try-error")) rep(NA_real_, 2L) else v
  }
  c(stats::coef(fit), se("robust"), se("information"), se("profile"))
}

cl <- makeCluster(NCORE)
on.exit(stopCluster(cl), add = TRUE)
invisible(clusterEvalQ(cl, library(pcdreg)))
# Distinct streams per cell, reproducible across reruns.
clusterSetRNGStream(cl, 20240612L + 1000L * CELL)

t0 <- Sys.time()
m <- do.call(rbind, parLapply(cl, seq_len(NREP), one_rep,
                              n = N, frailty = FRAILTY))
mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

ok <- stats::complete.cases(m)
m <- m[ok, , drop = FALSE]
est <- m[, 1:2, drop = FALSE]
covers <- function(cols) {
  100 * colMeans(abs(est - rep(TRUTH, each = nrow(est))) <=
                   Z * m[, cols, drop = FALSE])
}

out <- data.frame(
  n = N, poisson = if (FRAILTY == 0) "yes" else "no", coef = names(TRUTH),
  bias    = colMeans(est) - TRUTH,
  SE      = apply(est, 2, stats::sd),
  rob_SEE = colMeans(m[, 3:4, drop = FALSE]), rob_CP = covers(3:4),
  inf_SEE = colMeans(m[, 5:6, drop = FALSE]), inf_CP = covers(5:6),
  prf_SEE = colMeans(m[, 7:8, drop = FALSE]), prf_CP = covers(7:8),
  reps = nrow(est), minutes = round(mins, 1), row.names = NULL
)

dir.create("results", showWarnings = FALSE)
saveRDS(out, file.path("results", sprintf("cell%d.rds", CELL)))
print(format(out, digits = 3), row.names = FALSE)
cat(sprintf("done in %.1f min, %d of %d replications usable\n",
            mins, nrow(est), NREP))
