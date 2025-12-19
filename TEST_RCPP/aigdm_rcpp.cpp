// RcppArmadillo-based EM implementation
// Save as AIGDM_EM_Rcpp.cpp and compile via Rcpp::sourceCpp()

#include <RcppArmadillo.h>
#include <roptim.h>
// [[Rcpp::depends(RcppArmadillo, roptim)]]
using namespace Rcpp;
using namespace arma;

//------------------------------------------------------------------------------
// Utility: digamma, trigamma, lbeta for matrices
//------------------------------------------------------------------------------
// [[Rcpp::export]]
mat digamma_mat(const mat& X) {
  mat Y = X;
  for(uword i = 0; i < X.n_elem; ++i) {
    Y[i] = R::digamma(X[i]);
  }
  return Y;
}

// [[Rcpp::export]]
mat trigamma_mat(const mat& X) {
  mat Y = X;
  for(uword i = 0; i < X.n_elem; ++i) {
    Y[i] = R::trigamma(X[i]);
  }
  return Y;
}

mat lbeta_mat(const mat& A, const mat& B) {
  mat Y(A.n_rows, A.n_cols);
  for(uword i = 0; i < A.n_elem; ++i) {
    double a = A[i], b = B[i];
    Y[i] = lgamma(a) + lgamma(b) - lgamma(a + b);
  }
  return Y;
}

//------------------------------------------------------------------------------
// E-step: posterior
//------------------------------------------------------------------------------
// [[Rcpp::export]]
List ZIeZ_mat(const mat& pZv,
              const mat& av,
              const mat& bv,
              const mat& Y) {
  int n = Y.n_rows;
  int K = Y.n_cols - 1;
  mat Y1    = Y.cols(0, K-1);
  mat avp   = av + Y1;
  vec rowsum = sum(Y, 1);
  mat repSum = repmat(rowsum, 1, K);
  mat csum   = cumsum(Y1, 1);
  mat bvp    = bv + repSum - csum;
  
  mat logDen = lbeta_mat(av, bv);
  mat logNum = lbeta_mat(avp, bvp);
  mat betaRatio = exp(logNum - logDen);
  
  // posterior zero-inflation weights
  mat tmp     = pZv + (1.0 - pZv) % betaRatio;
  mat frac    = pZv / tmp;
  mat pv_post = zeros<mat>(n, K);
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < K; ++j) {
      if(Y1(i,j) == 0.0) {
        double denom = tmp(i,j);
        pv_post(i,j) = denom == 0.0 ? 1.0 : frac(i,j);
      }
    }
  }
  
  return List::create(
    Named("pv.post") = pv_post,
    Named("av.post") = avp,
    Named("bv.post") = bvp
  );
}

class LogitObjective : public roptim::Functor {
public:
  const vec& Del;
  const mat& W;
  LogitObjective(const vec& Del_, const mat& W_) : Del(Del_), W(W_) {}
  double operator()(const vec& par) {
    vec p = 1.0 / (1.0 + exp(-clamp(W * par, -30, 30)));
    return -sum(Del % log(p) + (1.0 - Del) % log(1.0 - p));
  }
  vec Gradient(const vec& par) {
    vec p = 1.0 / (1.0 + exp(-clamp(W * par, -30, 30)));
    return -W.t() * (Del - p);
  }
};

class AIBetaObjective : public roptim::Functor {
public:
  const vec& Del;
  const vec& A;
  const vec& B;
  const mat& Xa;
  const mat& Xb;
  AIBetaObjective(const vec& Del_, const vec& A_, const vec& B_, const mat& Xa_, const mat& Xb_) : Del(Del_), A(A_), B(B_), Xa(Xa_), Xb(Xb_) {}
  double operator()(const vec& par) {
    int da = Xa.n_cols;
    vec alpha = par.head(da);
    vec beta  = par.tail(par.n_elem - da);
    vec mu  = 1.0 / (1.0 + exp(-clamp(Xa * alpha, -30, 30)));
    vec phi = 1.0 / (1.0 + exp(-clamp(Xb * beta, -30, 30)));
    vec z = 1.0 - Del;
    vec t = (1.0/phi) - 1.0;
    vec a = max(t % mu, zeros<vec>(t.n_elem));
    vec b = max(t % (1.0 - mu), zeros<vec>(t.n_elem));
    double ll = 0.0;
    for (size_t i = 0; i < Del.n_elem; ++i) {
      ll += z(i) * (-R::lbeta(a(i), b(i)) + A(i)*(a(i)-1.0) + B(i)*(b(i)-1.0));
    }
    return -ll;
  }
  vec Gradient(const vec& par) {
    int da = Xa.n_cols;
    vec alpha = par.head(da);
    vec beta  = par.tail(par.n_elem - da);
    vec mu  = 1.0 / (1.0 + exp(-clamp(Xa * alpha, -30, 30)));
    vec phi = 1.0 / (1.0 + exp(-clamp(Xb * beta, -30, 30)));
    vec z = 1.0 - Del;
    vec t = (1.0/phi) - 1.0;
    vec a = max(t % mu, zeros<vec>(t.n_elem));
    vec b = max(t % (1.0 - mu), zeros<vec>(t.n_elem));
    vec grad(par.n_elem, fill::zeros);
    for (size_t i = 0; i < Del.n_elem; ++i) {
      double da_ = R::digamma(a(i)), db_ = R::digamma(b(i)), dab = R::digamma(a(i) + b(i));
      double zs = A(i) - B(i);
      double zd = B(i);
      double mu_i = mu(i), phi_i = phi(i);
      rowvec xa = Xa.row(i), xb = Xb.row(i);
      double dmu = mu_i * (1.0 - mu_i);
      double fac = (1.0/phi_i) - 1.0;
      vec ga = fac * dmu * (zs - (da_ - db_)) * xa.t();
      double tmpb = exp(as_scalar(xb * beta));
      vec gb = (-1.0/tmpb) * (mu_i*(zs - (da_ - db_)) + (zd - (db_ - dab))) * xb.t();
      grad.subvec(0, da-1)          -= z(i) * ga;
      grad.subvec(da, par.n_elem-1) -= z(i) * gb;
    }
    return grad;
  }
};

// [[Rcpp::export]]
List AIGDM_EM_Rcpp(const mat& Y, const mat& W, const mat& Xa, const mat& Xb, mat gamma, mat alpha, mat beta, double tol = 1e-4, int max_iter = 1000) {
  int n = Y.n_rows, K = Y.n_cols - 1, da = Xa.n_cols;
  mat DelR, AR, BR;
  mat gamma_last = gamma, alpha_last = alpha, beta_last = beta;
  for (int iter = 0; iter < max_iter; ++iter) {
    mat m1 = exp(Xa * alpha);
    mat m2 = exp(Xb * beta);
    mat m3 = exp(W * gamma);
    mat tmpMv = m1 / (1.0 + m1);
    mat tmpSv = m2 / (1.0 + m2);
    mat tmpZ  = m3 / (1.0 + m3);
    tmpZ.replace(datum::nan, 0.0);
    List post = ZIeZ_mat(tmpZ, tmpMv % (1/tmpSv - 1), (1-tmpMv) % (1/tmpSv - 1), Y);
    DelR = post["pv.post"];
    mat avp = post["av.post"], bvp = post["bv.post"], sumab = avp + bvp;
    AR.set_size(n, K); BR.set_size(n, K);
    for (int i = 0; i < n; ++i) for (int j = 0; j < K; ++j) {
      AR.at(i,j) = R::digamma(avp.at(i,j)) - R::digamma(sumab.at(i,j));
      BR.at(i,j) = R::digamma(bvp.at(i,j)) - R::digamma(sumab.at(i,j));
    }
    for (int j = 0; j < K; ++j) {
      LogitObjective logit(DelR.col(j), W);
      roptim::Roptim<LogitObjective> opt_logit(logit);
      opt_logit.set_initial(vec(gamma.col(j)));
      opt_logit.minimize();
      gamma.col(j) = opt_logit.par();
      
      AIBetaObjective aibeta(DelR.col(j), AR.col(j), BR.col(j), Xa, Xb);
      vec init = join_cols(alpha.col(j), beta.col(j));
      roptim::Roptim<AIBetaObjective> opt_beta(aibeta);
      opt_beta.set_initial(init);
      opt_beta.minimize();
      vec tmpAB = opt_beta.par();
      alpha.col(j) = tmpAB.head(da);
      beta.col(j)  = tmpAB.tail(tmpAB.n_elem - da);
    }
    double diff = norm(vectorise(gamma - gamma_last), 1) + norm(vectorise(alpha - alpha_last), 1) + norm(vectorise(beta - beta_last), 1);
    if (diff < tol) break;
    gamma_last = gamma; alpha_last = alpha; beta_last = beta;
  }
  return List::create(Named("gamma.est") = gamma, Named("alpha.est") = alpha, Named("beta.est") = beta);
}

// // ---- Logistic negative loglikelihood and score ----
// // [[Rcpp::export]]
// double LogitNegLoglik(const vec& par,
//                       const vec& Del,
//                       const mat& W) {
//   vec eta = W * par;
//   vec p   = 1.0 / (1.0 + exp(-eta));
//   vec ll  = Del % log(p) + (1.0 - Del) % log(1.0 - p);
//   uvec bad = find(p == 0.0 || p == 1.0);
//   ll.elem(bad).zeros();
//   return -sum(ll);
// }
// // [[Rcpp::export]]
// vec LogitNegScore(const vec& par,
//                   const vec& Del,
//                   const mat& W) {
//   vec p = 1.0 / (1.0 + exp(-W * par));
//   return -W.t() * (Del - p);
// }
// 
// // ---- Beta-regression negative loglikelihood and score ----
// // [[Rcpp::export]]
// double AIBetaNegLoglik(const vec& par,
//                        const vec& Del,
//                        const vec& A,
//                        const vec& B,
//                        const mat& Xa,
//                        const mat& Xb) {
//   int da = Xa.n_cols;
//   vec alpha = par.head(da);
//   vec beta  = par.tail(par.n_elem - da);
//   vec mu   = 1.0 / (1.0 + exp(-Xa * alpha));
//   vec phi  = 1.0 / (1.0 + exp(-Xb * beta));
//   vec z    = 1.0 - Del;
//   vec t    = (1.0/phi) - 1.0;
//   vec aval = max(t % mu, zeros<vec>(t.n_elem));
//   vec bval = max(t % (1.0 - mu), zeros<vec>(t.n_elem));
//   double ll = 0.0;
//   for (size_t i = 0; i < Del.n_elem; ++i) {
//     double w  = z(i);
//     double ai = aval(i), bi = bval(i);
//     double lb = R::lbeta(ai, bi);
//     ll += w * (-lb + A(i)*(ai - 1.0) + B(i)*(bi - 1.0));
//   }
//   return -ll;
// }
// // [[Rcpp::export]]
// vec AIBetaNegScore(const vec& par,
//                    const vec& Del,
//                    const vec& A,
//                    const vec& B,
//                    const mat& Xa,
//                    const mat& Xb) {
//   int da = Xa.n_cols;
//   vec alpha = par.head(da);
//   vec beta  = par.tail(par.n_elem - da);
//   vec mu   = 1.0/(1.0 + exp(-Xa * alpha));
//   vec phi  = 1.0/(1.0 + exp(-Xb * beta));
//   vec z    = 1.0 - Del;
//   vec t    = (1.0/phi) - 1.0;
//   vec aval = max(t % mu, zeros<vec>(t.n_elem));
//   vec bval = max(t % (1.0 - mu), zeros<vec>(t.n_elem));
//   vec grad(par.n_elem, fill::zeros);
//   for (size_t i = 0; i < Del.n_elem; ++i) {
//     double w   = z(i);
//     double ai  = aval(i), bi = bval(i);
//     double da_ = R::digamma(ai);
//     double db_ = R::digamma(bi);
//     double dab = R::digamma(ai + bi);
//     double zs  = A(i) - B(i);
//     double zd  = B(i);
//     double mu_i  = mu(i);
//     double phi_i = phi(i);
//     rowvec xa = Xa.row(i);
//     rowvec xb = Xb.row(i);
//     double dmu = mu_i * (1.0 - mu_i);
//     double fac = (1.0/phi_i) - 1.0;
//     vec ga = fac * dmu * (zs - (da_ - db_)) * xa.t();
//     double tmpb = exp(as_scalar(Xb.row(i) * beta));
//     vec gb = (-1.0/tmpb) * (mu_i*(zs - (da_ - db_)) + (zd - (db_ - dab))) * xb.t();
//     grad.subvec(0, da-1)          -= w * ga;
//     grad.subvec(da, par.n_elem-1) -= w * gb;
//   }
//   return grad;
// }



// // ---- Full EM driver, calling R's optim from C++ ----
// // [[Rcpp::export]]
// List AIGDM_EM_Rcpp(const mat& Y,
//                    const mat& W,
//                    const mat& Xa,
//                    const mat& Xb,
//                    mat gamma,
//                    mat alpha,
//                    mat beta,
//                    double tol = 1e-4,
//                    int max_iter = 1000) {
//   int n  = Y.n_rows;
//   int K  = Y.n_cols - 1;
//   int da = Xa.n_cols;
//   Function optim("optim");
//   mat DelR(n, K), AR(n, K), BR(n, K);
// 
//   // initialize "last" copies for convergence test
//   mat gamma_last = gamma;
//   mat alpha_last = alpha;
//   mat beta_last  = beta;
//   for (int iter = 0; iter < max_iter; ++iter) {
//     // E-step
//     mat m1 = exp(Xa * alpha);
//     mat m2 = exp(Xb * beta);
//     mat m3 = exp(W  * gamma);
//     mat tmpMv = m1 / (1.0 + m1);
//     mat tmpSv = m2 / (1.0 + m2);
//     mat tmpZ  = m3 / (1.0 + m3);
//     tmpZ.replace(datum::nan, 0.0);
//     tmpZ.elem(find_finite(tmpZ)).replace(datum::inf, 1.0);
//     List post = ZIeZ_mat(tmpZ,
//                          tmpMv % (1/tmpSv - 1),
//                          (1 - tmpMv) % (1/tmpSv - 1),
//                          Y);
//     mat DelR = post["pv.post"];
//     mat avp = post["av.post"];
//     mat bvp = post["bv.post"];
//     mat sumab = avp + bvp;
//     // elementwise digamma
//     for (int i = 0; i < n; ++i) for (int j = 0; j < K; ++j) {
//       AR(i,j) = R::digamma(avp(i,j)) - R::digamma(sumab(i,j));
//       BR(i,j) = R::digamma(bvp(i,j)) - R::digamma(sumab(i,j));
//     }
// 
//     // M-step
//     for (int j = 0; j < K; ++j) {
//       if (is_finite(gamma(0, j))) {
//         List out = optim(
//           _["par"]   = gamma.col(j),
//           _["fn"]    = Function("LogitNegLoglik"),
//           _["gr"]    = Function("LogitNegScore"),
//           _["Del"]   = DelR.col(j),
//           _["W"]     = W,
//           _["method"] = "L-BFGS-B"
//         );
//         gamma.col(j) = as<vec>(out["par"]);
//       }
//       vec init = join_cols(alpha.col(j), beta.col(j));
//       List ob = optim(
//         _["par"]   = init,
//         _["fn"]    = Function("AIBetaNegLoglik"),
//         _["gr"]    = Function("AIBetaNegScore"),
//         _["Del"]   = DelR.col(j),
//         _["A"]     = AR.col(j),
//         _["B"]     = BR.col(j),
//         _["Xa"]    = Xa,
//         _["Xb"]    = Xb,
//         _["method"] = "BFGS"
//       );
//       vec res = as<vec>(ob["par"]);
//       alpha.col(j) = res.subvec(0, da-1);
//       beta.col(j)  = res.subvec(da, res.n_elem-1);
//     }
// 
//     // convergence check against last iteration
//     double diff = norm(vectorise(gamma - gamma_last), 1)
//       + norm(vectorise(alpha - alpha_last), 1)
//       + norm(vectorise(beta  - beta_last),  1);
//       if (diff < tol) break;
//       // update last-copies
//       gamma_last = gamma;
//       alpha_last = alpha;
//       beta_last  = beta;
//   }
//   return List::create(
//     Named("gamma.est") = gamma,
//     Named("alpha.est") = alpha,
//     Named("beta.est")  = beta
//   );
// }


// ---- Intercept-only Beta-Binomial log-likelihood (BB.log.lik) ----
// [[Rcpp::export]]
List BB_log_lik(const mat& Y) {
  int n = Y.n_rows;
  // design matrix of ones
  mat X = ones<mat>(n,1);
  // initial small values and -Inf for gamma
  mat gamma0 = mat(1,1, fill::value(R_NegInf));
  mat alpha0 = mat(1,1, fill::value(0.001));
  mat beta0  = mat(1,1, fill::value(0.001));
  // run EM
  List gdm = AIGDM_EM_Rcpp(Y, X, X, X, gamma0, alpha0, beta0);
  mat alpha_est = as<mat>(gdm["alpha.est"]);
  mat beta_est  = as<mat>(gdm["beta.est"]);
  // compute mu and phi
  vec mu  = 1.0 / (1.0 + exp(-X * alpha_est));
  vec phi = 1.0 / (1.0 + exp(-X * beta_est));
  // shape parameters
  vec a = max(mu % (1.0/phi - 1.0), zeros<vec>(n));
  vec b = max((1.0 - mu) % (1.0/phi - 1.0), zeros<vec>(n));
  // original counts
  vec y = Y.col(0);
  vec m = sum(Y,1);
  // log-likelihood sum
  vec term = lgamma(m + 1.0) + lgamma(a + b) - lgamma(a) - lgamma(b)
    - lgamma(a + m + b)
    + lgamma(y + a) + lgamma(m - y + b)
    - lgamma(y + 1.0) - lgamma(m - y + 1.0);
    double ans = sum(term);
    return List::create(Named("ans") = ans,
                        Named("a.mat") = a,
                        Named("b.mat") = b);
}

// ---- Intercept-only BB and Zero-Inflated BB log-likelihood (BBnZIBB.log.lik) ----
// [[Rcpp::export]]
List BBnZIBB_log_lik(const mat& Y) {
  int n = Y.n_rows;
  mat X = ones<mat>(n,1);
  mat gamma0 = mat(1,1, fill::value(R_NegInf));
  mat alpha0 = mat(1,1, fill::value(0.001));
  mat beta0  = mat(1,1, fill::value(0.001));
  // BB fit
  List bbfit = AIGDM_EM_Rcpp(Y, X, X, X, gamma0, alpha0, beta0);
  vec mu_bb  = 1.0/(1.0+exp(-X*as<mat>(bbfit["alpha.est"])));
  vec phi_bb = 1.0/(1.0+exp(-X*as<mat>(bbfit["beta.est"])));
  vec a_bb = max(mu_bb % (1.0/phi_bb - 1.0), zeros<vec>(n));
  vec b_bb = max((1.0 - mu_bb) % (1.0/phi_bb - 1.0), zeros<vec>(n));
  vec y = Y.col(0), m = sum(Y,1);
  vec term_bb = lgamma(m + 1.0) + lgamma(a_bb + b_bb) - lgamma(a_bb) - lgamma(b_bb)
    - lgamma(a_bb + m + b_bb)
    + lgamma(y + a_bb) + lgamma(m - y + b_bb)
    - lgamma(y + 1.0) - lgamma(m - y + 1.0);
    double bb_ans = sum(term_bb);
    // ZIBB fit (positive-only)
    mat gamma1 = mat(1,1, fill::value(0.001));
    List zifit = AIGDM_EM_Rcpp(Y, X, X, X, gamma1, alpha0, beta0);
    vec mu_zi  = 1.0/(1.0+exp(-X*as<mat>(zifit["alpha.est"])));
    vec phi_zi = 1.0/(1.0+exp(-X*as<mat>(zifit["beta.est"])));
    vec a_zi = max(mu_zi % (1.0/phi_zi - 1.0), zeros<vec>(n));
    vec b_zi = max((1.0 - mu_zi) % (1.0/phi_zi - 1.0), zeros<vec>(n));
    // select positive y
    uvec idx = find(y > 0);
    vec y2 = y.elem(idx);
    vec m2 = m.elem(idx);
    vec a2 = a_zi.elem(idx);
    vec b2 = b_zi.elem(idx);
    // compute underflow-safe terms
    vec logA = lgamma(a2 + m2 + b2) + lgamma(b2);
    vec logB = lgamma(a2 + b2)     + lgamma(m2 + b2);
    // elementwise safe probability and log
    vec prob = 1.0 - exp(logB - logA);
    vec safeP = clamp(prob, DBL_MIN, datum::inf);
    vec term_zi = -log(safeP) - logA
        + lgamma(m2 + 1.0)
        + lgamma(y2 + a2) + lgamma(m2 - y2 + b2)
        + lgamma(a2 + b2)
        - lgamma(y2 + 1.0) - lgamma(m2 - y2 + 1.0) - lgamma(a2);
    double zi_ans = sum(term_zi);
    return List::create(Named("bb.ans") = bb_ans,
                        Named("zibb.ans") = zi_ans);
}


// [[Rcpp::export]]
IntegerMatrix simData_GDM(const mat& a,
                          const mat& b,
                          const vec& SeqDepth) {
  int N = a.n_rows;
  int K = a.n_cols;
  
  // Draw Z ~ Beta(a, b)
  mat Z(N, K);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < K; ++j) {
      Z(i, j) = R::rbeta(a(i, j), b(i, j));
    }
  }
  
  // Compute cumulative products of (1 - Z)
  mat P(N, K+1, fill::zeros);
  for (int i = 0; i < N; ++i) {
    double cum_prob = 1.0;
    for (int j = 0; j < K; ++j) {
      P(i, j) = Z(i, j) * cum_prob;
      cum_prob *= (1.0 - Z(i, j));
    }
    // last category
    P(i, K) = cum_prob;
    // enforce non-negativity
    for (int j = 0; j <= K; ++j) {
      if (P(i, j) < 0) P(i, j) = 0;
    }
  }
  
  // Simulate Y ~ Multinomial(SeqDepth[i], P(i,:)) rowwise
  IntegerMatrix Y(N, K+1);
  for (int i = 0; i < N; ++i) {
    // Rf_rmultinom(int n, const double *probs, int k, int *ans)
    std::vector<double> probs(K+1);
    for (int j = 0; j <= K; ++j) probs[j] = P(i, j);
    std::vector<int> counts(K+1);
    Rf_rmultinom((int)SeqDepth[i], probs.data(), K+1, counts.data());
    for (int j = 0; j <= K; ++j) {
      Y(i, j) = counts[j];
    }
  }
  
  return Y;
}
