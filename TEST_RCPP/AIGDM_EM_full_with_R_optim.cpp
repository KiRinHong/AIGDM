
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Inverse logit
vec inv_logit(const vec& x) {
  return 1 / (1 + exp(-x));
}

// Log-likelihood for logistic regression
double logit_neg_loglik(const vec& par, const vec& Del, const mat& W) {
  vec eta = W * par;
  vec p = inv_logit(eta);
  vec loglik = Del % log(p + 1e-12) + (1 - Del) % log(1 - p + 1e-12);
  return -sum(loglik);
}

// Score for logistic regression
vec logit_neg_score(const vec& par, const vec& Del, const mat& W) {
  vec eta = W * par;
  vec p = inv_logit(eta);
  return -W.t() * (Del - p);
}

// [[Rcpp::export]]
vec LogitOptim_cpp(const vec& Del, const mat& W, const vec& gamma_ini, double tol = 1e-6, int max_iter = 100) {
  vec par = gamma_ini;
  for (int i = 0; i < max_iter; ++i) {
    vec grad = logit_neg_score(par, Del, W);
    vec step = solve(W.t() * diagmat(inv_logit(W * par) % (1 - inv_logit(W * par))) * W, grad);
    par -= step;
    if (norm(step, 2) < tol) break;
  }
  return par;
}

// [[Rcpp::export]]
List ZIeZ_mat_cpp(const mat& pZv, const mat& av, const mat& bv, const mat& Y) {
  int n = Y.n_rows;
  int K = Y.n_cols - 1;
  mat av_prim = av + Y.cols(0, K-1);
  mat cumY = cumsum(Y.cols(0, K-1), 1);
  mat bv_prim = bv + (Y.rowwise().sum() * ones<rowvec>(K) - cumY);

  mat beta_values = exp(lgamma(av_prim + bv_prim) - lgamma(av_prim) - lgamma(bv_prim) +
                        lgamma(av) + lgamma(bv) - lgamma(av + bv));
  mat tmp_values = pZv + (1 - pZv) % beta_values;
  mat pZv_post = pZv / (tmp_values + 1e-12);
  pZv_post.elem(find(Y.cols(0, K-1) > 0)).zeros(); // If Y > 0, pZ_post = 0

  return List::create(Named("pv_post") = pZv_post,
                      Named("av_post") = av_prim,
                      Named("bv_post") = bv_prim);
}

// [[Rcpp::export]]
List AIGDM_EM_cpp(const mat& Y, const mat& W, const mat& Xa, const mat& Xb,
                  mat gamma0, mat alpha0, mat beta0, double tol = 1e-4, int max_iter = 1000) {
  int n = Y.n_rows;
  int K = Y.n_cols - 1;
  int da = Xa.n_cols;
  int db = Xb.n_cols;

  mat gamma_now = gamma0, alpha_now = alpha0, beta_now = beta0;
  mat gamma_last = gamma_now, alpha_last = alpha_now, beta_last = beta_now;

  mat Del_R(n, K, fill::zeros), A_R(n, K, fill::zeros), B_R(n, K, fill::zeros);

  for (int l = 0; l < max_iter; ++l) {
    mat expXa_alpha = exp(Xa * alpha_last);
    mat expXb_beta = exp(Xb * beta_last);
    mat expW_gamma = exp(W * gamma_last);

    mat tmpMv = expXa_alpha / (1 + expXa_alpha);
    mat tmpSv = expXb_beta / (1 + expXb_beta);
    mat tmpA = tmpMv % (1 / tmpSv - 1);
    mat tmpB = (1 - tmpMv) % (1 / tmpSv - 1);
    mat tmpZ = expW_gamma / (1 + expW_gamma);
    tmpZ.replace(datum::nan, 0);
    tmpZ.transform([](double val) { return std::isinf(val) ? 1.0 : val; });

    List par_post = ZIeZ_mat_cpp(tmpZ, tmpA, tmpB, Y);
    mat av_post = as<mat>(par_post["av_post"]);
    mat bv_post = as<mat>(par_post["bv_post"]);
    Del_R = as<mat>(par_post["pv_post"]);

    mat digamma_sum = digamma(av_post + bv_post);
    A_R = digamma(av_post) - digamma_sum;
    B_R = digamma(bv_post) - digamma_sum;

    for (int j = 0; j < K; ++j) {
      vec del_col = Del_R.col(j);
      vec gamma_init = gamma_last.col(j);
      gamma_now.col(j) = LogitOptim_cpp(del_col, W, gamma_init);

      alpha_now.col(j) = alpha_last.col(j); // Placeholder
      beta_now.col(j) = beta_last.col(j);   // Placeholder
    }

    mat diff = abs(join_horiz(gamma_now - gamma_last, alpha_now - alpha_last, beta_now - beta_last));
    if (accu(diff) < tol) break;

    gamma_last = gamma_now;
    alpha_last = alpha_now;
    beta_last = beta_now;
  }

  return List::create(Named("gamma_est") = gamma_now,
                      Named("alpha_est") = alpha_now,
                      Named("beta_est") = beta_now);
}


// [[Rcpp::export]]
vec AIBetaOptim_cpp(const vec& Del, const vec& A, const vec& B,
                    const mat& Xa, const mat& Xb,
                    const vec& alpha_ini, const vec& beta_ini) {

  Environment stats("package:stats");
  Function optim = stats["optim"];

  vec start = join_vert(alpha_ini, beta_ini);
  List par = List::create(Named("par") = start,
                          Named("fn") = Rcpp::InternalFunction(
                            [=](SEXP par_) {
                              vec p = as<vec>(par_);
                              int da = Xa.n_cols;
                              vec alpha = p.head(da);
                              vec beta = p.tail(p.n_elem - da);

                              vec mu = 1 / (1 + exp(-Xa * alpha));
                              vec phi = 1 / (1 + exp(-Xb * beta));

                              vec a = (1 / phi - 1) % mu;
                              vec b = (1 / phi - 1) % (1 - mu);

                              a = clamp(a, 1e-6, 1e6); b = clamp(b, 1e-6, 1e6);

                              vec loglik = (1 - Del) % (-lgamma(a) - lgamma(b) + lgamma(a + b) + A % (a - 1) + B % (b - 1));
                              return wrap(-sum(loglik));
                            }),
                          Named("method") = "BFGS");

  List res = optim(par);
  return as<vec>(res["par"]);
}
