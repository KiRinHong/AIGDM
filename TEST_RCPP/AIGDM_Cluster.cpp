
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List AIGDM_Cluster(NumericVector ID, NumericMatrix OTU, NumericMatrix X, IntegerVector X_index,
                   Nullable<NumericMatrix> Tax = R_NilValue, std::string test_type = "Omni",
                   int min_depth = 0, std::string ZI_pos = "adaptive", int n_boot = 500,
                   Nullable<IntegerVector> n_perm = R_NilValue, int n_cores = 1, double fdr_alpha = 0.05) {

    // Input validation
    if (X_index.size() == 0) {
        stop("X_index is empty.");
    }
    if (is_true(any(X_index < 1)) || is_true(any(X_index > X.ncol()))) {
        stop("X_index out of bounds.");
    }
    if (!(test_type == "Mean" || test_type == "Disp" || test_type == "Zero" || test_type == "Omni")) {
        stop("test_type should be one of {'Mean', 'Disp', 'Zero', 'Omni'}!");
    }
    if (!(ZI_pos == "no" || ZI_pos == "all" || ZI_pos == "adaptive")) {
        stop("ZI_pos should be one of {'no', 'all', 'adaptive'}!");
    }
    if (ZI_pos == "adaptive" && n_boot <= 0) {
        stop("n_boot must be > 0 when ZI_pos = 'adaptive'.");
    }

    // Remove samples with read depth less than min_depth
    LogicalVector keep_samples(OTU.nrow(), true);
    for (int i = 0; i < OTU.nrow(); ++i) {
        if (sum(OTU(i, _)) < min_depth) {
            keep_samples[i] = false;
        }
    }
    NumericVector new_ID = ID[keep_samples];
    NumericMatrix new_OTU = OTU(keep_samples, _);
    NumericMatrix new_X = X(keep_samples, _);

    // Add intercept to design matrix
    NumericMatrix X_with_intercept(new_X.nrow(), new_X.ncol() + 1);
    for (int i = 0; i < new_X.nrow(); ++i) {
        X_with_intercept(i, 0) = 1.0;
        for (int j = 0; j < new_X.ncol(); ++j) {
            X_with_intercept(i, j + 1) = new_X(i, j);
        }
    }

    return List::create(
        _["ID"] = new_ID,
        _["OTU"] = new_OTU,
        _["X"] = X_with_intercept
    );
}
