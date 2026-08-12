#include "panelrate.h"

using namespace panelrate;

// [[Rcpp::depends(RcppArmadillo)]]

// Fit the proportional rate model by the EM algorithm of Section 2 of the
// paper.
//
// Plain EM converges linearly and, on this problem, slowly: the jump sizes can
// still be drifting in the fifth decimal place after thousands of passes.  Each
// outer iteration therefore takes two EM passes and extrapolates along them by
// the SQUAREM scheme of Varadhan and Roland (2008), keeping the extrapolated
// point only when a stabilising EM pass from it does at least as well on the
// observed data log likelihood as the two plain passes did.  The fixed point is
// unchanged; only the number of passes needed to reach it falls.  Setting
// accelerate to false recovers plain EM, one pass per iteration.
//
// [[Rcpp::export]]
Rcpp::List em_fit_cpp(const arma::mat& X, const arma::uvec& subj,
                      const arma::uvec& grid, const arma::uvec& panel,
                      const arma::vec& dN, const arma::uvec& panelsubj,
                      arma::uword n, arma::uword K, arma::vec beta,
                      arma::vec lambda, arma::uword maxit, double reltol,
                      bool accelerate) {
  const PanelData d = {X, subj, grid, panel, dN, panelsubj, n, K, dN.n_elem};
  const arma::uword p = d.p();

  Parameters cur;
  cur.beta = beta;
  cur.lambda = lambda;

  bool converged = false;
  double criterion = R_PosInf;
  arma::uword iter = 0;
  arma::uword passes = 0;

  for (iter = 0; iter < maxit; ++iter) {
    if (iter % 10 == 0) Rcpp::checkUserInterrupt();

    Parameters next;
    if (!em_step(d, cur, next)) {
      Rcpp::stop("The EM algorithm failed: the information matrix is singular "
                 "or the linear predictor overflowed. Check for collinear or "
                 "badly scaled covariates.");
    }
    ++passes;

    // Convergence is judged on the residual of one plain EM pass, which is how
    // far the current point is from being a fixed point.  Judging it on the
    // accelerated jump instead would confuse a large extrapolation with a lack
    // of convergence.
    criterion = std::max(p > 0 ? rel_change(cur.beta, next.beta, reltol) : 0.0,
                         rel_change(cur.lambda, next.lambda, reltol));

    if (accelerate) {
      Parameters second;
      if (!em_step(d, next, second)) {
        Rcpp::stop("The EM algorithm failed after one pass; check for "
                   "collinear or badly scaled covariates.");
      }
      ++passes;
      double best_ll = loglik_at(d, second);
      const arma::vec theta = pack(cur);
      const arma::vec r = pack(next) - theta;
      const arma::vec v = pack(second) - pack(next) - r;
      const double vnorm = arma::norm(v);
      Parameters best = second;

      if (vnorm > 0.0 && std::isfinite(best_ll)) {
        double alpha = -arma::norm(r) / vnorm;
        for (int back = 0; back < 10 && alpha < -1.0; ++back) {
          Parameters trial =
              unpack(theta - 2.0 * alpha * r + alpha * alpha * v, p);
          trial.lambda = arma::clamp(trial.lambda, 0.0, arma::datum::inf);
          Parameters stabilised;
          if (em_step(d, trial, stabilised)) {
            ++passes;
            const double ll = loglik_at(d, stabilised);
            if (std::isfinite(ll) && ll >= best_ll) {
              best = stabilised;
              best_ll = ll;
              break;
            }
          }
          alpha = (alpha - 1.0) / 2.0;
        }
      }
      next = best;
    }

    cur = next;

    if (criterion < reltol) {
      converged = true;
      ++iter;
      break;
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("beta") = cur.beta, Rcpp::Named("lambda") = cur.lambda,
      Rcpp::Named("loglik") = loglik_at(d, cur),
      Rcpp::Named("iterations") = iter, Rcpp::Named("passes") = passes,
      Rcpp::Named("converged") = converged,
      Rcpp::Named("criterion") = criterion);
}

// The two covariance ingredients of Section 4, evaluated at the fitted values.
//
//   Omega = -(1/n) sum_i sum_j g_ij g_ij' / denom_ij,
//   S     =  (1/n) sum_i u_i u_i',   u_i = sum_k (What_ik - mu_ik)(X_ik - xbar_k),
//
// where g_ij = sum_{t_k in (T_{i,j-1}, T_ij]} mu_ik (X_ik - xbar_k) and
// mu_ik = eta_ik lambda_k.  Note that u_i = sum_j (Delta N_ij / denom_ij - 1)
// g_ij, so u_i is subject i's contribution to the score at the solution.
//
// [[Rcpp::export]]
Rcpp::List covariance_cpp(const arma::mat& X, const arma::uvec& subj,
                          const arma::uvec& grid, const arma::uvec& panel,
                          const arma::vec& dN, const arma::uvec& panelsubj,
                          arma::uword n, arma::uword K, const arma::vec& beta,
                          const arma::vec& lambda) {
  const PanelData d = {X, subj, grid, panel, dN, panelsubj, n, K, dN.n_elem};
  const arma::uword p = d.p();

  arma::vec eta;
  if (!relative_rate(d, beta, eta)) {
    Rcpp::stop("The linear predictor overflowed while computing the "
               "covariance; check for badly scaled covariates.");
  }
  arma::vec S0;
  arma::mat S1;
  risk_sums(d, eta, S0, S1);

  arma::mat xbar(p, d.K, arma::fill::zeros);
  for (arma::uword k = 0; k < d.K; ++k) {
    if (invertible(S0(k))) xbar.col(k) = S1.col(k) / S0(k);
  }

  const arma::vec denom = panel_totals(d, eta, lambda);
  const arma::vec What = estep(d, eta, lambda, denom);

  // g_ij, accumulated over the grid times inside each examination interval.
  arma::mat g(p, d.P, arma::fill::zeros);
  arma::mat scores(d.n, p, arma::fill::zeros);
  std::vector<double> centred(p);
  for (arma::uword r = 0; r < d.M(); ++r) {
    const double mu = eta(r) * lambda(d.grid(r));
    const double* x = d.X.colptr(r);
    const double* xb = xbar.colptr(d.grid(r));
    double* gp = g.colptr(d.panel(r));
    for (arma::uword j = 0; j < p; ++j) {
      centred[j] = x[j] - xb[j];
      gp[j] += mu * centred[j];
      scores(d.subj(r), j) += (What(r) - mu) * centred[j];
    }
  }

  arma::mat Omega(p, p, arma::fill::zeros);
  double* Om = Omega.memptr();
  for (arma::uword pp = 0; pp < d.P; ++pp) {
    if (invertible(denom(pp))) {
      add_outer_upper(Om, g.colptr(pp), p, -1.0 / denom(pp));
    }
  }
  Omega = arma::symmatu(Omega) / static_cast<double>(n);

  arma::mat S(p, p, arma::fill::zeros);
  double* Sp = S.memptr();
  std::vector<double> u(p);
  for (arma::uword i = 0; i < d.n; ++i) {
    for (arma::uword j = 0; j < p; ++j) u[j] = scores(i, j);
    add_outer_upper(Sp, u.data(), p, 1.0);
  }
  S = arma::symmatu(S) / static_cast<double>(n);

  return Rcpp::List::create(Rcpp::Named("Omega") = Omega, Rcpp::Named("S") = S,
                            Rcpp::Named("scores") = scores);
}

// Profile log likelihood gradient of Murphy and van der Vaart (2000) in the
// form used by Zeng et al. (2020) and reported as "Profile" in the paper:
//
//   G[i, l] = { pl_i(betahat + h e_l) - pl_i(betahat) } / h,
//
// where pl_i(beta) = l_i(beta, Lambdahat(beta)) and Lambdahat(beta) maximises
// the likelihood over the jump sizes with beta held fixed.  The covariance
// estimate is then solve(crossprod(G)).  Each column costs one lambda-only EM
// run, warm started from the fitted lambda, rather than one run per pair of
// coefficients.
//
// [[Rcpp::export]]
Rcpp::List profile_gradient_cpp(const arma::mat& X, const arma::uvec& subj,
                                const arma::uvec& grid, const arma::uvec& panel,
                                const arma::vec& dN, const arma::uvec& panelsubj,
                                arma::uword n, arma::uword K,
                                const arma::vec& beta, const arma::vec& lambda,
                                double h, arma::uword maxit, double reltol) {
  const PanelData d = {X, subj, grid, panel, dN, panelsubj, n, K, dN.n_elem};
  const arma::uword p = d.p();

  arma::mat G(d.n, p, arma::fill::zeros);
  arma::uvec iterations(p + 1, arma::fill::zeros);
  arma::uvec converged(p + 1, arma::fill::zeros);

  // Baseline pl_i(betahat).  Re-solving for lambda costs one or two iterations
  // and keeps the baseline on the same footing as the perturbed fits.
  arma::vec eta;
  if (!relative_rate(d, beta, eta)) {
    Rcpp::stop("The linear predictor overflowed while computing the profile "
               "likelihood; check for badly scaled covariates.");
  }
  arma::vec S0;
  arma::mat S1;
  risk_sums(d, eta, S0, S1);
  arma::uword it = 0;
  bool ok = false;
  arma::vec lam = profile_lambda(d, eta, S0, lambda, maxit, reltol, it, ok);
  iterations(0) = it;
  converged(0) = ok ? 1 : 0;
  const arma::vec base = subject_loglik(d, panel_totals(d, eta, lam));

  for (arma::uword l = 0; l < p; ++l) {
    Rcpp::checkUserInterrupt();
    arma::vec beta_l = beta;
    beta_l(l) += h;
    if (!relative_rate(d, beta_l, eta)) {
      Rcpp::stop("The linear predictor overflowed while perturbing the "
                 "coefficients for the profile likelihood.");
    }
    risk_sums(d, eta, S0, S1);
    lam = profile_lambda(d, eta, S0, lambda, maxit, reltol, it, ok);
    iterations(l + 1) = it;
    converged(l + 1) = ok ? 1 : 0;
    G.col(l) = (subject_loglik(d, panel_totals(d, eta, lam)) - base) / h;
  }

  return Rcpp::List::create(Rcpp::Named("gradient") = G,
                            Rcpp::Named("h") = h,
                            Rcpp::Named("iterations") = iterations,
                            Rcpp::Named("converged") = converged);
}
