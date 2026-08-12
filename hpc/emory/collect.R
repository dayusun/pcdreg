# Combine the six cells into Table 1 and print the paper's values beside them.
#
#   Rscript collect.R

files <- sort(Sys.glob("results/cell*.rds"))
if (!length(files)) stop("no results yet")
tab <- do.call(rbind, lapply(files, readRDS))
tab <- tab[order(tab$n, tab$poisson != "yes", tab$coef), ]
rownames(tab) <- NULL

cat("=== Replication ===\n")
print(format(tab[, c("n", "poisson", "coef", "bias", "SE", "rob_SEE",
                     "rob_CP", "inf_SEE", "inf_CP", "prf_SEE", "prf_CP",
                     "reps")], digits = 3), row.names = FALSE)

# Table 1 of the paper, lambda(t) = 8/(1+t), 1000 replications.
paper <- data.frame(
  n = rep(c(100L, 200L, 400L), each = 4L),
  poisson = rep(rep(c("yes", "no"), each = 2L), 3L),
  coef = rep(c("x1", "x2"), 6L),
  bias    = c(-0.009, 0.002, -0.008, -0.013, -0.005, 0.001, -0.007, -0.005,
              -0.001, -0.001, -0.005, -0.003),
  SE      = c(0.123, 0.119, 0.365, 0.387, 0.087, 0.084, 0.251, 0.273,
              0.061, 0.060, 0.179, 0.187),
  rob_SEE = c(0.119, 0.116, 0.332, 0.367, 0.085, 0.082, 0.242, 0.266,
              0.060, 0.058, 0.174, 0.191),
  rob_CP  = c(93.6, 94.4, 92.4, 93.2, 94.6, 94.1, 94.2, 94.6,
              94.9, 95.0, 93.6, 95.6),
  inf_SEE = c(0.128, 0.123, 0.048, 0.041, 0.088, 0.085, 0.032, 0.027,
              0.061, 0.059, 0.022, 0.018),
  inf_CP  = c(95.4, 95.5, 21.9, 15.6, 95.3, 95.3, 20.1, 14.3,
              95.1, 95.3, 20.1, 15.1)
)

cat("\n=== Paper, Table 1 ===\n")
print(paper, row.names = FALSE)

key <- function(d) paste(d$n, d$poisson, d$coef)
m <- match(key(tab), key(paper))
diff <- data.frame(
  n = tab$n, poisson = tab$poisson, coef = tab$coef,
  d_SE      = round(tab$SE      - paper$SE[m], 3),
  d_rob_SEE = round(tab$rob_SEE - paper$rob_SEE[m], 3),
  d_rob_CP  = round(tab$rob_CP  - paper$rob_CP[m], 1),
  d_inf_SEE = round(tab$inf_SEE - paper$inf_SEE[m], 3),
  d_inf_CP  = round(tab$inf_CP  - paper$inf_CP[m], 1)
)
cat("\n=== Replication minus paper ===\n")
print(diff, row.names = FALSE)
cat("\nAt 1000 replications a coverage estimate carries a Monte Carlo standard\n",
    "error near 0.7 points and an empirical standard deviation about 2%, so\n",
    "differences inside that are noise rather than disagreement.\n", sep = "")

saveRDS(tab, "table1.rds")
utils::write.csv(tab, "table1.csv", row.names = FALSE)
