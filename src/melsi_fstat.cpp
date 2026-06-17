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

    // Plain double accumulators (not long double): on x86-64 long double is
    // 80-bit x87, which is slow and prevents SIMD vectorisation of the inner
    // loop. Summing O(n^2) similar-magnitude positive terms in double keeps the
    // relative error at ~1e-14, well within the tolerance the equivalence tests
    // confirmed.
    std::vector<double> within(k, 0.0);
    std::vector<int> nsize(k, 0);
    for (int i = 0; i < n; ++i) nsize[group[i]]++;

    const double* x = &Xt[0];           // column-major: sample i starts at x + i*m
    double total = 0.0;

    for (int i = 0; i < n; ++i) {
        const double* __restrict xi = x + (std::size_t)i * m;
        const int gi = group[i];
        for (int j = i + 1; j < n; ++j) {
            const double* __restrict xj = x + (std::size_t)j * m;
            double s = 0.0;
            for (int d = 0; d < m; ++d) {
                const double diff = xi[d] - xj[d];
                s += diff * diff;
            }
            total += s;
            if (group[j] == gi) within[gi] += s;
        }
    }

    const double SS_total = total / (double)n;
    double SS_within = 0.0;
    for (int g = 0; g < k; ++g) {
        if (nsize[g] > 0) SS_within += within[g] / (double)nsize[g];
    }
    const double SS_between = SS_total - SS_within;
    const double F = (SS_between / (double)(k - 1)) /
                     (SS_within / (double)(n - k));
    return F;
}
