#include "pcdreg.h"

using namespace pcdreg;

// [[Rcpp::depends(RcppArmadillo)]]

// Fit the proportional means model
//
//     E[N(t) | X(t)] = mu(t) exp(beta' X(t))
//
// by the estimating equation of Hu, Sun and Wei (2003),
//
//     U(beta) = sum_i sum_j N_i(T_ij) { X_i(T_ij) - xbar(T_ij) } = 0,
//
// where xbar(t) averages the covariates over the subjects examined at t,
// weighted by exp(beta' X).  This is the shared weighted score kernel with the
// observed cumulative counts as weights, so a Newton step is
// beta + (-Udot)^{-1} U exactly as in the rate model's M-step.
//
// Unlike the rate model there is no likelihood here, so there is no log
// likelihood to report and only the sandwich covariance is available.
//
// Rows are examinations rather than grid times: row r is one (i, j) pair,
// carrying X_i(T_ij), the subject index, the index of T_ij in the pooled grid,
// and the cumulative count cN_r = N_i(T_ij).
//
// [[Rcpp::export]]
Rcpp::List mean_fit_cpp(const arma::mat& X, const arma::uvec& subj,
                        const arma::uvec& grid, const arma::vec& cN,
                        arma::uword n, arma::uword K, arma::vec beta,
                        arma::uword maxit, double reltol) {
  const arma::uword p = X.n_rows, M = X.n_cols;

  arma::vec eta, S0, U;
  arma::mat S1, negUdot;
  bool converged = false;
  double criterion = R_PosInf;
  arma::uword iter = 0;

  for (iter = 0; iter < maxit; ++iter) {
    if (iter % 20 == 0) Rcpp::checkUserInterrupt();
    if (!relative_rate(X, beta, eta)) {
      Rcpp::stop("The linear predictor overflowed; check for collinear or "
                 "badly scaled covariates.");
    }
    risk_sums(X, grid, K, eta, S0, S1);

    if (p == 0) {  // nothing to solve for; the baseline is still estimable
      converged = true;
      criterion = 0.0;
      break;
    }

    weighted_score(X, grid, K, eta, cN, S0, S1, U, negUdot);
    arma::vec step;
    if (!safe_solve(negUdot, U, step)) {
      Rcpp::stop("The estimating equation has a singular derivative; the "
                 "covariates may be collinear.");
    }
    const arma::vec beta_new = beta + step;
    criterion = rel_change(beta, beta_new, reltol);
    beta = beta_new;
    if (criterion < reltol) {
      converged = true;
      ++iter;
      break;
    }
  }

  // Final quantities at the solution.
  if (!relative_rate(X, beta, eta)) {
    Rcpp::stop("The linear predictor overflowed at the solution.");
  }
  risk_sums(X, grid, K, eta, S0, S1);
  const arma::mat xbar = covariate_means(S0, S1);
  if (p > 0) weighted_score(X, grid, K, eta, cN, S0, S1, U, negUdot);
  else negUdot.zeros(p, p);

  // Baseline mean, muhat_k = sum_r cN_r / S0_k over examinations at t_k.
  //
  // Note this is a level rather than a jump, and nothing constrains it to
  // increase.  When covariates vary over time a fitted mean function that goes
  // down is exactly the difficulty with the means model that motivates the rate
  // model, so it is reported as estimated rather than quietly forced upwards.
  arma::vec num(K, arma::fill::zeros);
  for (arma::uword r = 0; r < M; ++r) num(grid(r)) += cN(r);
  arma::vec mu(K, arma::fill::zeros);
  for (arma::uword k = 0; k < K; ++k) {
    if (invertible(S0(k))) mu(k) = num(k) / S0(k);
  }

  // Per subject contribution to the estimating function, and its empirical
  // covariance: u_i = sum_{r in i} (X_r - xbar_k) (cN_r - eta_r mu_k).
  arma::mat scores(n, p, arma::fill::zeros);
  for (arma::uword r = 0; r < M; ++r) {
    const arma::uword k = grid(r);
    const double resid = cN(r) - eta(r) * mu(k);
    const double* x = X.colptr(r);
    const double* xb = xbar.colptr(k);
    for (arma::uword j = 0; j < p; ++j) {
      scores(subj(r), j) += resid * (x[j] - xb[j]);
    }
  }

  arma::mat S(p, p, arma::fill::zeros);
  double* Sp = S.memptr();
  std::vector<double> u(p);
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = 0; j < p; ++j) u[j] = scores(i, j);
    add_outer_upper(Sp, u.data(), p, 1.0);
  }
  S = arma::symmatu(S) / static_cast<double>(n);

  // Udot itself, negative semidefinite, so that the sandwich assembled in R is
  // the same expression as for the rate model.
  const arma::mat Omega = -negUdot / static_cast<double>(n);

  return Rcpp::List::create(
      Rcpp::Named("beta") = beta, Rcpp::Named("mu") = mu,
      Rcpp::Named("Omega") = Omega, Rcpp::Named("S") = S,
      Rcpp::Named("scores") = scores, Rcpp::Named("iterations") = iter,
      Rcpp::Named("converged") = converged,
      Rcpp::Named("criterion") = criterion);
}
