#ifndef PCDREG_H
#define PCDREG_H

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// Core routines for the semiparametric proportional rate model
//
//     E[dN(t) | X(t)] = exp(beta' X(t)) dLambda(t)
//
// fitted by the EM algorithm of Sun, Guo, Li, Tu and Sun (2024),
// Bernoulli 30(4), doi:10.3150/23-BEJ1713.  Equation numbers in the comments
// refer to that paper.
//
// Lambda is a step function with jump lambda_k at the pooled distinct
// examination time t_k, k = 1, ..., K.

namespace pcdreg {

// Expanded data layout.
//
// A "row" r = 0, ..., M-1 is one (subject, grid time) pair for which the
// subject is still under observation, i.e. one (i, k) with Delta_ik = 1.
// Row r carries the covariate vector X.col(r) = X_i(t_k), the grid index
// grid(r) = k, the subject index subj(r) = i, and the index panel(r) of the
// examination interval (T_{i,j-1}, T_ij] that contains t_k.  Counts are stored
// once per examination interval: dN(panel(r)) = Delta N_ij.
struct PanelData {
  const arma::mat& X;         // p x M
  const arma::uvec& subj;     // M, values in [0, n)
  const arma::uvec& grid;     // M, values in [0, K)
  const arma::uvec& panel;    // M, values in [0, P)
  const arma::vec& dN;        // P
  const arma::uvec& panelsubj;// P, subject owning each examination interval
  arma::uword n;              // number of subjects
  arma::uword K;              // number of pooled grid times
  arma::uword P;              // number of examination intervals

  arma::uword p() const { return X.n_rows; }
  arma::uword M() const { return X.n_cols; }
};

// The parameters of the model: the coefficients and the jump sizes of Lambda.
struct Parameters {
  arma::vec beta;
  arma::vec lambda;
};

// The routines below deliberately write their inner loops out by hand instead
// of forming matrix products.  The number of covariates p is small while the
// number of rows M is large, so every such product would be a tall, thin gemm
// producing a p by p result.  A multithreaded BLAS spawns threads for each of
// those, and since they run thousands of times per fit the synchronisation cost
// swamps the arithmetic: on a machine with a threaded OpenBLAS this made fits
// roughly forty times slower, nearly all of it system time.

// eta_r = exp(beta' X_i(t_k)).  Returns false if the linear predictor overflows,
// which can happen at a trial point proposed by the accelerator.
inline bool relative_rate(const PanelData& d, const arma::vec& beta,
                          arma::vec& eta) {
  const arma::uword p = d.p(), M = d.M();
  eta.set_size(M);
  const double* b = beta.memptr();
  for (arma::uword r = 0; r < M; ++r) {
    const double* x = d.X.colptr(r);
    double lp = 0.0;
    for (arma::uword j = 0; j < p; ++j) lp += b[j] * x[j];
    eta(r) = std::exp(lp);
  }
  return eta.is_finite();
}

// Whether a denominator can be inverted without overflowing.
//
// The nonparametric maximum likelihood estimate sets many jump sizes to exactly
// zero, so an examination interval can end up with a total intensity that is
// positive but subnormal.  Inverting such a value overflows to infinity, and
// since the matching numerator has underflowed to exactly zero the product is
// NaN rather than a number.  Every quantity divided by such a denominator is
// bounded by a multiple of the denominator itself, so its limit is zero and
// dropping the term is correct as well as safe.
inline bool invertible(double x) {
  return x >= std::numeric_limits<double>::min();
}

// Accumulate c * x x' into the upper triangle of a p by p buffer.
inline void add_outer_upper(double* H, const double* x, arma::uword p,
                            double c) {
  for (arma::uword j = 0; j < p; ++j) {
    const double cxj = c * x[j];
    for (arma::uword i = 0; i <= j; ++i) H[i + j * p] += cxj * x[i];
  }
}

// S0_k = sum_i Delta_ik eta_ik and S1_k = sum_i Delta_ik eta_ik X_i(t_k),
// the unnormalised versions of S^(0) and S^(1) (the factor 1/n cancels in
// every ratio in which they are used).
inline void risk_sums(const PanelData& d, const arma::vec& eta,
                      arma::vec& S0, arma::mat& S1) {
  const arma::uword p = d.p();
  S0.zeros(d.K);
  S1.zeros(p, d.K);
  for (arma::uword r = 0; r < d.M(); ++r) {
    const arma::uword k = d.grid(r);
    S0(k) += eta(r);
    const double* x = d.X.colptr(r);
    double* s = S1.colptr(k);
    for (arma::uword j = 0; j < p; ++j) s[j] += eta(r) * x[j];
  }
}

// Total intensity of each examination interval,
//     denom_ij = sum_{T_{i,j-1} < t_k <= T_ij} eta_ik lambda_k,
// which is the Poisson mean of Delta N_ij implied by (beta, lambda).
inline arma::vec panel_totals(const PanelData& d, const arma::vec& eta,
                              const arma::vec& lambda) {
  arma::vec denom(d.P, arma::fill::zeros);
  for (arma::uword r = 0; r < d.M(); ++r) {
    denom(d.panel(r)) += eta(r) * lambda(d.grid(r));
  }
  return denom;
}

// E-step: the posterior mean of the latent Poisson count W_ik,
//
//     What_ik = eta_ik lambda_k Delta N_ij / denom_ij,
//
// which is the mean of a Binomial(Delta N_ij, eta_ik lambda_k / denom_ij).
inline arma::vec estep(const PanelData& d, const arma::vec& eta,
                       const arma::vec& lambda, const arma::vec& denom) {
  arma::vec What(d.M(), arma::fill::zeros);
  for (arma::uword r = 0; r < d.M(); ++r) {
    const arma::uword pp = d.panel(r);
    if (invertible(denom(pp)) && d.dN(pp) > 0.0) {
      What(r) = eta(r) * lambda(d.grid(r)) * d.dN(pp) / denom(pp);
    }
  }
  return What;
}

// M-step for the jump sizes, equation (4):
//
//     lambda_k = sum_i Delta_ik What_ik / sum_i Delta_ik eta_ik .
inline arma::vec mstep_lambda(const PanelData& d, const arma::vec& What,
                              const arma::vec& S0) {
  arma::vec num(d.K, arma::fill::zeros);
  for (arma::uword r = 0; r < d.M(); ++r) num(d.grid(r)) += What(r);
  arma::vec lambda(d.K, arma::fill::zeros);
  for (arma::uword k = 0; k < d.K; ++k) {
    if (invertible(S0(k))) lambda(k) = num(k) / S0(k);
  }
  return lambda;
}

// The M-step score U(beta) and its negative derivative -Udot(beta), where
//
//     U     = sum_{i,k} What_ik (X_i(t_k) - xbar_k),
//     -Udot = sum_{i,k} What_ik { S2_k / S0_k - xbar_k xbar_k' },
//
// with xbar_k = S1_k / S0_k the covariate mean at t_k weighted by the relative
// rate.  -Udot is the observed information of the M-step and is positive
// semidefinite, so the Newton step is beta + (-Udot)^{-1} U.
//
// The p x p x K array S2 is never formed: sum_k (w_k / S0_k) S2_k is
// accumulated in a single pass over rows.
inline void mstep_score(const PanelData& d, const arma::vec& eta,
                        const arma::vec& What, const arma::vec& S0,
                        const arma::mat& S1, arma::vec& U, arma::mat& negUdot) {
  const arma::uword p = d.p();
  arma::vec w(d.K, arma::fill::zeros);   // w_k = sum_i What_ik
  arma::vec xw(p, arma::fill::zeros);    // sum_{i,k} What_ik X_i(t_k)
  for (arma::uword r = 0; r < d.M(); ++r) {
    w(d.grid(r)) += What(r);
    const double* x = d.X.colptr(r);
    for (arma::uword j = 0; j < p; ++j) xw(j) += What(r) * x[j];
  }

  arma::mat xbar(p, d.K, arma::fill::zeros);
  for (arma::uword k = 0; k < d.K; ++k) {
    if (invertible(S0(k))) xbar.col(k) = S1.col(k) / S0(k);
  }

  U = xw;
  for (arma::uword k = 0; k < d.K; ++k) {
    const double* xb = xbar.colptr(k);
    for (arma::uword j = 0; j < p; ++j) U(j) -= w(k) * xb[j];
  }

  negUdot.zeros(p, p);
  double* H = negUdot.memptr();
  for (arma::uword r = 0; r < d.M(); ++r) {
    const arma::uword k = d.grid(r);
    if (invertible(S0(k)) && w(k) > 0.0) {
      add_outer_upper(H, d.X.colptr(r), p, (w(k) / S0(k)) * eta(r));
    }
  }
  for (arma::uword k = 0; k < d.K; ++k) {
    if (w(k) > 0.0) add_outer_upper(H, xbar.colptr(k), p, -w(k));
  }
  negUdot = arma::symmatu(negUdot);
}

// Observed data log likelihood contributed by each subject,
//
//     l_i = sum_j { Delta N_ij log(denom_ij) - denom_ij - log(Delta N_ij !) },
//
// the log of expression (2) evaluated at the step function Lambda.
inline arma::vec subject_loglik(const PanelData& d, const arma::vec& denom) {
  arma::vec ll(d.n, arma::fill::zeros);
  for (arma::uword pp = 0; pp < d.P; ++pp) {
    const double count = d.dN(pp);
    double term = -denom(pp) - std::lgamma(count + 1.0);
    if (count > 0.0) {
      term += (denom(pp) > 0.0) ? count * std::log(denom(pp))
                                : -std::numeric_limits<double>::infinity();
    }
    ll(d.panelsubj(pp)) += term;
  }
  return ll;
}

// Solve A x = b for a symmetric positive semidefinite A, reporting failure
// rather than returning silently wrong numbers.
inline bool safe_solve(const arma::mat& A, const arma::vec& b, arma::vec& out) {
  const arma::mat sym = arma::symmatu(A);
  if (arma::solve(out, sym, b,
                  arma::solve_opts::likely_sympd + arma::solve_opts::no_approx)) {
    return out.is_finite();
  }
  return arma::solve(out, sym, b) && out.is_finite();
}

// Relative change used as the convergence criterion, applied to both beta and
// the jump sizes.
inline double rel_change(const arma::vec& from, const arma::vec& to,
                         double eps) {
  const double scale = arma::max(arma::abs(from)) + eps;
  return arma::max(arma::abs(to - from)) / scale;
}

// One full pass of the algorithm of Section 2: E-step at the current values,
// update of the jump sizes by (4), a second E-step at the updated jump sizes,
// and one Newton step for beta.  Returns false if the pass is not numerically
// viable, which the accelerator uses to reject a trial point.
inline bool em_step(const PanelData& d, const Parameters& cur, Parameters& out) {
  arma::vec eta;
  if (!relative_rate(d, cur.beta, eta)) return false;
  arma::vec S0;
  arma::mat S1;
  risk_sums(d, eta, S0, S1);

  arma::vec denom = panel_totals(d, eta, cur.lambda);
  arma::vec What = estep(d, eta, cur.lambda, denom);
  out.lambda = mstep_lambda(d, What, S0);
  out.beta = cur.beta;

  if (d.p() > 0) {
    denom = panel_totals(d, eta, out.lambda);
    What = estep(d, eta, out.lambda, denom);
    arma::vec U;
    arma::mat negUdot;
    mstep_score(d, eta, What, S0, S1, U, negUdot);
    arma::vec step;
    if (!safe_solve(negUdot, U, step)) return false;
    out.beta = cur.beta + step;
  }
  return out.lambda.is_finite() && out.beta.is_finite();
}

// Observed data log likelihood, minus infinity if the parameters are not
// viable.  Used to guard the accelerated steps.
inline double loglik_at(const PanelData& d, const Parameters& s) {
  arma::vec eta;
  if (!relative_rate(d, s.beta, eta)) {
    return -std::numeric_limits<double>::infinity();
  }
  return arma::accu(subject_loglik(d, panel_totals(d, eta, s.lambda)));
}

inline arma::vec pack(const Parameters& s) {
  return arma::join_cols(s.beta, s.lambda);
}

inline Parameters unpack(const arma::vec& v, arma::uword p) {
  Parameters out;
  out.beta = v.head(p);
  out.lambda = v.tail(v.n_elem - p);
  return out;
}

// One pass of the baseline-only EM, holding beta (and hence eta) fixed.
inline arma::vec lambda_map(const PanelData& d, const arma::vec& eta,
                            const arma::vec& S0, const arma::vec& lambda) {
  const arma::vec denom = panel_totals(d, eta, lambda);
  const arma::vec What = estep(d, eta, lambda, denom);
  return mstep_lambda(d, What, S0);
}

// Maximise over lambda holding beta fixed, i.e. compute Lambda_hat(beta) by
// iterating the E-step and (4) alone.  This is what the profile likelihood
// needs, and it needs it accurately, because the profile log likelihoods are
// differenced and divided by a small step size.  It is accelerated the same way
// as the main algorithm; with beta fixed this really is an EM, so the observed
// data log likelihood increases at every plain pass and makes a dependable
// guard for the extrapolated ones.
inline arma::vec profile_lambda(const PanelData& d, const arma::vec& eta,
                                const arma::vec& S0, arma::vec lambda,
                                arma::uword maxit, double reltol,
                                arma::uword& iter, bool& converged) {
  converged = false;
  for (iter = 0; iter < maxit; ++iter) {
    if (iter % 20 == 0) Rcpp::checkUserInterrupt();

    const arma::vec first = lambda_map(d, eta, S0, lambda);
    const double change = rel_change(lambda, first, reltol);

    arma::vec next = lambda_map(d, eta, S0, first);
    double best_ll = arma::accu(subject_loglik(d, panel_totals(d, eta, next)));

    const arma::vec r = first - lambda;
    const arma::vec v = next - first - r;
    const double vnorm = arma::norm(v);
    if (vnorm > 0.0 && std::isfinite(best_ll)) {
      double alpha = -arma::norm(r) / vnorm;
      for (int back = 0; back < 10 && alpha < -1.0; ++back) {
        arma::vec trial = lambda - 2.0 * alpha * r + alpha * alpha * v;
        trial = arma::clamp(trial, 0.0, arma::datum::inf);
        const arma::vec stabilised = lambda_map(d, eta, S0, trial);
        const double ll =
            arma::accu(subject_loglik(d, panel_totals(d, eta, stabilised)));
        if (std::isfinite(ll) && ll >= best_ll) {
          next = stabilised;
          best_ll = ll;
          break;
        }
        alpha = (alpha - 1.0) / 2.0;
      }
    }

    lambda = next;
    if (change < reltol) {
      converged = true;
      ++iter;
      break;
    }
  }
  return lambda;
}

}  // namespace pcdreg

#endif  // PCDREG_H
