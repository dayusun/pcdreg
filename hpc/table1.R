# Replicate Table 1 of Sun, Guo, Li, Tu and Sun (2024), Bernoulli 30(4).
#
# Design of Section 5: two covariates, x1 a single step at U ~ Unif(0, tau) with
# both levels Unif(0, 1), x2 constant Unif(0, 1); beta = (1, -1);
# lambda(t) = 8 / (1 + t); tau = 2; J zero-truncated Poisson with mean 4;
# examination times the order statistics of a uniform sample on (0, Tbar) with
# Tbar ~ Unif(0.9 tau, tau).  n = 100, 200, 400, with and without the Poisson
# assumption, the latter a gamma frailty of variance one.
#
# Usage:  Rscript table1.R [nrep] [ncores]
# Defaults to the paper's 1000 replications and whatever Slurm allocated.
#
# reltol is loosened to 1e-6.  A calibration run put the cost at under 5e-6 per
# coefficient, four orders of magnitude below the Monte Carlo noise of a 1000
# replication study, while making 6000 fits with the profile likelihood
# affordable.

suppressMessages({
  library(parallel)
  library(pcdreg)
})

args  <- commandArgs(trailingOnly = TRUE)
NREP  <- if (length(args) >= 1) as.integer(args[1]) else 1000L
NCORE <- if (length(args) >= 2) as.integer(args[2]) else
  as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))

NVEC  <- c(100L, 200L, 400L)
TRUTH <- c(x1 = 1, x2 = -1)
Z     <- stats::qnorm(0.975)

cat(sprintf("pcdreg %s | %d replications | %d cores\n",
            as.character(utils::packageVersion("pcdreg")), NREP, NCORE))

one_rep <- function(i, n, frailty) {
  d <- pcdreg::sim_pcd(n, beta = c(1, -1),
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

summarise <- function(m, poisson, n) {
  m <- m[stats::complete.cases(m), , drop = FALSE]
  est <- m[, 1:2, drop = FALSE]
  covers <- function(cols) {
    100 * colMeans(abs(est - rep(TRUTH, each = nrow(est))) <=
                     Z * m[, cols, drop = FALSE])
  }
  data.frame(
    n = n, poisson = poisson, coef = names(TRUTH),
    bias    = round(colMeans(est) - TRUTH, 3),
    SE      = round(apply(est, 2, stats::sd), 3),
    rob_SEE = round(colMeans(m[, 3:4, drop = FALSE]), 3),
    rob_CP  = round(covers(3:4), 1),
    inf_SEE = round(colMeans(m[, 5:6, drop = FALSE]), 3),
    inf_CP  = round(covers(5:6), 1),
    prf_SEE = round(colMeans(m[, 7:8, drop = FALSE]), 3),
    prf_CP  = round(covers(7:8), 1),
    reps = nrow(est), row.names = NULL
  )
}

cl <- makeCluster(NCORE)
on.exit(stopCluster(cl), add = TRUE)
invisible(clusterEvalQ(cl, library(pcdreg)))
clusterSetRNGStream(cl, 20240612)

res <- list()
for (n in NVEC) {
  for (fr in c(0, 1)) {
    t0 <- Sys.time()
    m <- do.call(rbind, parLapply(cl, seq_len(NREP), one_rep,
                                  n = n, frailty = fr))
    res[[length(res) + 1L]] <- summarise(m, if (fr == 0) "yes" else "no", n)
    cat(sprintf("n = %3d, poisson = %-3s : %5.1f min\n", n,
                if (fr == 0) "yes" else "no",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

tab <- do.call(rbind, res)
rownames(tab) <- NULL

cat("\n=== Replication of Table 1 ===\n")
print(tab, row.names = FALSE)
saveRDS(tab, "table1.rds")
utils::write.csv(tab, "table1.csv", row.names = FALSE)

cat("\nPaper, Table 1, lambda(t) = 8/(1+t), 1000 replications:\n")
cat("  n=100 yes: bias -0.009/0.002  SE 0.123/0.119  rob 0.119(93.6)/0.116(94.4)  inf 0.128(95.4)/0.123(95.5)\n")
cat("  n=100 no : bias -0.008/-0.013 SE 0.365/0.387  rob 0.332(92.4)/0.367(93.2)  inf 0.048(21.9)/0.041(15.6)\n")
cat("  n=200 yes: bias -0.005/0.001  SE 0.087/0.084  rob 0.085(94.6)/0.082(94.1)  inf 0.088(95.3)/0.085(95.3)\n")
cat("  n=200 no : bias -0.007/-0.005 SE 0.251/0.273  rob 0.242(94.2)/0.266(94.6)  inf 0.032(20.1)/0.027(14.3)\n")
cat("  n=400 yes: bias -0.001/-0.001 SE 0.061/0.060  rob 0.060(94.9)/0.058(95.0)  inf 0.061(95.1)/0.059(95.3)\n")
cat("  n=400 no : bias -0.005/-0.003 SE 0.179/0.187  rob 0.174(93.6)/0.191(95.6)  inf 0.022(20.1)/0.018(15.1)\n")
