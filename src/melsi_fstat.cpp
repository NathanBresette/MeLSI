#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// Fused PERMANOVA pseudo-F statistic computed directly from data.
//
// Xt is the (already weighted/scaled) data transposed to features-by-samples
// (m x n) so each sample's features are contiguous in memory. group holds
// 0-based group ids of length n; k is the number of groups.
//
// Squared Euclidean distances are accumulated in a single pass over the
// n*(n-1)/2 unique sample pairs, so the full n x n distance matrix is never
// materialised. Totals use long double to limit floating-point error.
//
// Equivalent (in exact arithmetic) to the R reference:
//   D2 <- as.matrix(dist(t(Xt)))^2
//   SS_total  <- sum(D2) / (2n)
//   SS_within <- sum_g sum(D2[idx_g, idx_g]) / (2 n_g)
//   F <- (SS_between/(k-1)) / (SS_within/(n-k))
//
// [[Rcpp::export]]
double melsi_permanova_f(const NumericMatrix& Xt, const IntegerVector& group, int k) {
    const int m = Xt.nrow();
    const int n = Xt.ncol();

    std::vector<long double> within(k, 0.0L);
    std::vector<int> nsize(k, 0);
    for (int i = 0; i < n; ++i) nsize[group[i]]++;

    const double* x = &Xt[0];           // column-major: sample i starts at x + i*m
    long double total = 0.0L;

    for (int i = 0; i < n; ++i) {
        const double* xi = x + (std::size_t)i * m;
        const int gi = group[i];
        for (int j = i + 1; j < n; ++j) {
            const double* xj = x + (std::size_t)j * m;
            double s = 0.0;
            for (int d = 0; d < m; ++d) {
                const double diff = xi[d] - xj[d];
                s += diff * diff;
            }
            total += s;
            if (group[j] == gi) within[gi] += s;
        }
    }

    const long double SS_total = total / (long double)n;
    long double SS_within = 0.0L;
    for (int g = 0; g < k; ++g) {
        if (nsize[g] > 0) SS_within += within[g] / (long double)nsize[g];
    }
    const long double SS_between = SS_total - SS_within;
    const long double F = (SS_between / (long double)(k - 1)) /
                          (SS_within / (long double)(n - k));
    return (double)F;
}
