# A deliberately naive pure R implementation of the EM algorithm and of the two
# covariance matrices, written straight from the equations in the paper.  It
# shares no code with src/, so agreement between the two is real evidence that
# the compiled version computes what the paper says.

# Sum x within groups g (0 based, n groups), keeping empty groups as zero.
group_sum <- function(x, g, n) {
  x <- as.matrix(x)
  s <- rowsum(x, g)
  out <- matrix(0, n, ncol(x))
  out[as.integer(rownames(s)) + 1L, ] <- s
  if (ncol(x) == 1L) drop(out) else out
}

ref_estep <- function(d, eta, lambda) {
  P <- length(d$dN)
  mu <- eta * lambda[d$grid + 1L]
  denom <- group_sum(mu, d$panel, P)
  scale <- ifelse(denom > 0, d$dN / denom, 0)
  list(mu = mu, denom = denom, W = mu * scale[d$panel + 1L])
}

ref_em <- function(d, maxit = 200L, reltol = 1e-15) {
  Xm <- t(d$X)
  K <- d$K
  p <- ncol(Xm)
  beta <- numeric(p)
  lambda <- rep(1 / K, K)

  for (it in seq_len(maxit)) {
    eta <- as.vector(exp(Xm %*% beta))
    S0 <- group_sum(eta, d$grid, K)
    S1 <- group_sum(Xm * eta, d$grid, K)
    xbar <- S1 / ifelse(S0 > 0, S0, 1)
    if (p == 1L) xbar <- matrix(xbar, ncol = 1L)

    e1 <- ref_estep(d, eta, lambda)
    lambda_new <- ifelse(S0 > 0, group_sum(e1$W, d$grid, K) / S0, 0)

    e2 <- ref_estep(d, eta, lambda_new)
    centred <- Xm - xbar[d$grid + 1L, , drop = FALSE]
    U <- colSums(e2$W * centred)
    wk <- group_sum(e2$W, d$grid, K)

    H <- matrix(0, p, p)
    for (r in seq_len(nrow(Xm))) {
      k <- d$grid[r] + 1L
      if (S0[k] > 0 && wk[k] > 0) {
        H <- H + (wk[k] / S0[k]) * eta[r] * tcrossprod(Xm[r, ])
      }
    }
    for (k in seq_len(K)) {
      if (wk[k] > 0) H <- H - wk[k] * tcrossprod(xbar[k, ])
    }
    beta_new <- if (p > 0) beta + solve(H, U) else beta

    crit <- max(
      if (p > 0) max(abs(beta_new - beta)) / (max(abs(beta)) + reltol) else 0,
      max(abs(lambda_new - lambda)) / (max(abs(lambda)) + reltol)
    )
    beta <- beta_new
    lambda <- lambda_new
    if (crit < reltol) break
  }
  list(beta = beta, lambda = lambda, iterations = it)
}

ref_covariance <- function(d, beta, lambda) {
  Xm <- t(d$X)
  K <- d$K
  P <- length(d$dN)
  p <- ncol(Xm)
  eta <- as.vector(exp(Xm %*% beta))
  S0 <- group_sum(eta, d$grid, K)
  S1 <- group_sum(Xm * eta, d$grid, K)
  xbar <- S1 / ifelse(S0 > 0, S0, 1)
  if (p == 1L) xbar <- matrix(xbar, ncol = 1L)

  e <- ref_estep(d, eta, lambda)
  centred <- Xm - xbar[d$grid + 1L, , drop = FALSE]

  g <- group_sum(e$mu * centred, d$panel, P)
  if (p == 1L) g <- matrix(g, ncol = 1L)
  Omega <- matrix(0, p, p)
  for (j in seq_len(P)) {
    if (e$denom[j] > 0) Omega <- Omega - tcrossprod(g[j, ]) / e$denom[j]
  }

  scores <- group_sum((e$W - e$mu) * centred, d$subj, d$n)
  if (p == 1L) scores <- matrix(scores, ncol = 1L)

  list(Omega = Omega / d$n, S = crossprod(scores) / d$n, scores = scores)
}

# The means model of Hu, Sun and Wei (2003), again straight from the equations:
# Newton iteration on U(beta) = sum_r cN_r (X_r - xbar_k(r)), over examinations
# rather than over grid times.
ref_mean <- function(d, maxit = 100L, reltol = 1e-12) {
  Xm <- t(d$exam_X)
  K <- d$K
  p <- ncol(Xm)
  grid <- d$exam_grid
  beta <- numeric(p)

  for (it in seq_len(maxit)) {
    eta <- as.vector(exp(Xm %*% beta))
    S0 <- group_sum(eta, grid, K)
    S1 <- group_sum(Xm * eta, grid, K)
    xbar <- S1 / ifelse(S0 > 0, S0, 1)
    if (p == 1L) xbar <- matrix(xbar, ncol = 1L)

    centred <- Xm - xbar[grid + 1L, , drop = FALSE]
    U <- colSums(d$cN * centred)
    wk <- group_sum(d$cN, grid, K)

    H <- matrix(0, p, p)
    for (r in seq_len(nrow(Xm))) {
      k <- grid[r] + 1L
      if (S0[k] > 0 && wk[k] > 0) {
        H <- H + (wk[k] / S0[k]) * eta[r] * tcrossprod(Xm[r, ])
      }
    }
    for (k in seq_len(K)) {
      if (wk[k] > 0) H <- H - wk[k] * tcrossprod(xbar[k, ])
    }

    beta_new <- beta + solve(H, U)
    crit <- max(abs(beta_new - beta)) / (max(abs(beta)) + reltol)
    beta <- beta_new
    if (crit < reltol) break
  }

  eta <- as.vector(exp(Xm %*% beta))
  S0 <- group_sum(eta, grid, K)
  mu <- ifelse(S0 > 0, group_sum(d$cN, grid, K) / S0, 0)
  list(beta = beta, mu = mu, iterations = it)
}

# Observed data log likelihood, summed over subjects.
ref_loglik <- function(d, beta, lambda) {
  eta <- as.vector(exp(t(d$X) %*% beta))
  denom <- group_sum(eta * lambda[d$grid + 1L], d$panel, length(d$dN))
  sum(ifelse(d$dN > 0, d$dN * log(denom), 0) - denom - lgamma(d$dN + 1))
}

# Build the prepared arrays the way pcdreg() does, for use by the reference.
prep <- function(data, formula = pcd(id, tstart, tstop, count) ~ x1 + x2) {
  mf <- stats::model.frame(formula, data)
  X <- stats::model.matrix(attr(mf, "terms"), mf)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  pcdreg:::prepare_panel(stats::model.response(mf), X)
}
