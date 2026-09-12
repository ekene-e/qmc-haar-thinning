// SPDX-License-Identifier: MIT
//
// Minimal use: draw n low-discrepancy points in [0,1)^d and integrate a
// function with them, comparing against plain Monte Carlo.

#include <cmath>
#include <cstdio>
#include <numbers>
#include <vector>

#include <haar/haar.hpp>

int main() {
    const unsigned d = 2;
    const std::size_t n = 1 << 16;

    haar::Options opt;
    opt.dim = d;
    opt.seed = 1;
    std::vector<double> pts = haar::thin(n, opt);  // n * d doubles, row-major

    // f(x) = sin(2 pi k x_1) has integral 0 and Hardy--Krause variation ~ k,
    // but smoothed-out variation only ~ sqrt(k): the beyond-HK regime.
    const double k = 32.0;
    auto f = [k](std::span<const double> x) { return std::sin(2.0 * std::numbers::pi * k * x[0]); };

    const double qmc = haar::mean({pts, d}, f);

    haar::Xoshiro256pp rng(1);
    std::vector<double> mc(n * d);
    rng.uniform_point(mc);
    const double plain = haar::mean({mc, d}, f);

    std::printf("n = %zu points in %u dimensions\n", n, d);
    std::printf("  Haar-thinning estimate: %+.3e\n", qmc);
    std::printf("  Monte Carlo estimate:   %+.3e\n", plain);
    std::printf("  star discrepancy:       thinned %.2f   Monte Carlo %.2f\n",
                *haar::star_discrepancy({pts, d}), *haar::star_discrepancy({mc, d}));
    return 0;
}
