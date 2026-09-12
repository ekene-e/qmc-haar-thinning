// SPDX-License-Identifier: MIT
//
// Online use with an external sample source.  Here the "source" is
// std::mt19937_64, standing in for whatever produces your i.i.d. uniform
// samples (a simulator, a file, a hardware generator).  The thinner never
// needs to know n in advance: with no `level`/`expected_n` it grows its
// resolution as points accumulate, so every prefix of the output is a good
// point set (a low-discrepancy *sequence*).

#include <cstdio>
#include <random>
#include <vector>

#include <haar/haar.hpp>

int main() {
    haar::Options opt;
    opt.dim = 3;
    opt.epsilon = 0.5;  // at most ~1.5 n samples are consumed for n points
    haar::Thinner thinner(opt);

    std::mt19937_64 source(2024);
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    std::vector<double> x(opt.dim);
    std::size_t consumed = 0;
    while (thinner.size() < 100000) {
        for (auto& c : x) c = unif(source);
        ++consumed;
        thinner.offer(x);  // true if x was kept; after a rejection the next x is always kept
    }

    const haar::Stats& s = thinner.stats();
    std::printf("%s\n", thinner.describe().c_str());
    std::printf("consumed %zu samples for %zu points (%.1f%% rejected, %llu resolution increases)\n",
                consumed, thinner.size(), 100.0 * static_cast<double>(s.rejected) / static_cast<double>(s.offered),
                static_cast<unsigned long long>(s.rebuilds));
    std::printf("largest |phi_t(H)| over all Haar functions: %lld\n",
                static_cast<long long>(thinner.max_abs_haar_discrepancy()));
    return 0;
}
