
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <random>

// [[Rcpp::export]]
Rcpp::NumericMatrix bootstrap_sampling(const Rcpp::NumericMatrix& data, int n_bootstrap) {
    int n_rows = data.nrow();
    int n_cols = data.ncol();
    Rcpp::NumericMatrix bootstrap_samples(n_bootstrap, n_cols);

    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    for (int b = 0; b < n_bootstrap; ++b) {
        std::uniform_int_distribution<> dis(0, n_rows - 1);
        for (int col = 0; col < n_cols; ++col) {
            bootstrap_samples(b, col) = data(dis(gen), col);
        }
    }
    return bootstrap_samples;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix permutation_test(const Rcpp::NumericMatrix& data, int n_permutations) {
    int n_rows = data.nrow();
    int n_cols = data.ncol();
    Rcpp::NumericMatrix permuted_data(n_permutations, n_cols);

    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    for (int p = 0; p < n_permutations; ++p) {
        std::vector<int> indices(n_rows);
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), gen);

        for (int col = 0; col < n_cols; ++col) {
            permuted_data(p, col) = data(indices[col], col);
        }
    }
    return permuted_data;
}

// [[Rcpp::export]]
Rcpp::NumericVector compute_statistics(const Rcpp::NumericMatrix& data) {
    int n_cols = data.ncol();
    Rcpp::NumericVector stats(n_cols);

    for (int col = 0; col < n_cols; ++col) {
        Rcpp::NumericVector column = data(Rcpp::_, col);
        stats[col] = Rcpp::mean(column);
    }
    return stats;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_distance_matrix(const Rcpp::NumericMatrix& data) {
    int n_rows = data.nrow();
    Rcpp::NumericMatrix distance_matrix(n_rows, n_rows);

    for (int i = 0; i < n_rows; ++i) {
        for (int j = i; j < n_rows; ++j) {
            double dist = 0.0;
            for (int k = 0; k < data.ncol(); ++k) {
                dist += std::pow(data(i, k) - data(j, k), 2);
            }
            dist = std::sqrt(dist);
            distance_matrix(i, j) = dist;
            distance_matrix(j, i) = dist;  // symmetry
        }
    }
    return distance_matrix;
}
