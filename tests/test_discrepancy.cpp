// SPDX-License-Identifier: MIT
#include <algorithm>
#include <cmath>
#include <vector>

#include "check.hpp"
#include "haar/discrepancy.hpp"
#include "haar/rng.hpp"

using haar::PointsView;

namespace {

// Brute force over the grid of candidate corners: every x and y coordinate
// (plus 1), evaluating both closed and half-open boxes.
double brute_force_star_2d(const std::vector<double>& pts) {
    const std::size_t n = pts.size() / 2;
    std::vector<double> xs{1.0}, ys{1.0};
    for (std::size_t i = 0; i < n; ++i) {
        xs.push_back(pts[2 * i]);
        ys.push_back(pts[2 * i + 1]);
    }
    double best = 0.0;
    for (double x : xs) {
        for (double y : ys) {
            std::size_t closed = 0, open = 0;
            for (std::size_t i = 0; i < n; ++i) {
                const double px = pts[2 * i], py = pts[2 * i + 1];
                closed += px <= x && py <= y;
                open += px < x && py < y;
            }
            const double vol = static_cast<double>(n) * x * y;
            best = std::max({best, static_cast<double>(closed) - vol, vol - static_cast<double>(open)});
        }
    }
    return best;
}

// Numerical integral of (F_n(x) - x)^2 over [0,1] for 1-D points.
double numeric_l2_star_1d(std::vector<double> xs) {
    std::sort(xs.begin(), xs.end());
    const std::size_t n = xs.size();
    const int grid = 400000;
    double sum = 0.0;
    std::size_t k = 0;
    for (int g = 0; g < grid; ++g) {
        const double x = (g + 0.5) / grid;
        while (k < n && xs[k] <= x) ++k;
        const double diff = static_cast<double>(k) / static_cast<double>(n) - x;
        sum += diff * diff;
    }
    return static_cast<double>(n) * std::sqrt(sum / grid);
}

}  // namespace

TEST(star_1d_known_values) {
    std::vector<double> centred(16), left(16);
    for (int i = 0; i < 16; ++i) {
        centred[static_cast<std::size_t>(i)] = (2 * i + 1) / 32.0;
        left[static_cast<std::size_t>(i)] = i / 16.0;
    }
    CHECK_CLOSE(*haar::star_discrepancy({centred, 1}), 0.5, 1e-12);
    CHECK_CLOSE(*haar::star_discrepancy({left, 1}), 1.0, 1e-12);
    std::vector<double> single{0.25};
    CHECK_CLOSE(*haar::star_discrepancy({single, 1}), 0.75, 1e-12);
}

TEST(star_2d_matches_brute_force) {
    haar::Xoshiro256pp rng(3);
    for (int trial = 0; trial < 30; ++trial) {
        const std::size_t n = 1 + static_cast<std::size_t>(rng() % 40);
        std::vector<double> pts(2 * n);
        for (auto& c : pts) c = rng.uniform();
        if (trial % 3 == 0) {  // inject ties in both coordinates
            for (std::size_t i = 0; i < pts.size(); ++i) pts[i] = std::floor(pts[i] * 4) / 4;
        }
        CHECK_CLOSE(*haar::star_discrepancy({pts, 2}), brute_force_star_2d(pts), 1e-9);
    }
}

TEST(star_2d_of_a_lattice) {
    // The 4x4 centred lattice has discrepancy exactly 1 - 1/16 + ... check via brute force only.
    std::vector<double> pts;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) pts.insert(pts.end(), {(2 * i + 1) / 8.0, (2 * j + 1) / 8.0});
    CHECK_CLOSE(*haar::star_discrepancy({pts, 2}), brute_force_star_2d(pts), 1e-12);
}

TEST(unsupported_dimension_returns_nullopt) {
    std::vector<double> pts(30, 0.5);
    CHECK(!haar::star_discrepancy({pts, 3}).has_value());
}

TEST(l2_star_matches_numeric_integral_in_1d) {
    haar::Xoshiro256pp rng(5);
    std::vector<double> pts(37);
    for (auto& c : pts) c = rng.uniform();
    CHECK_CLOSE(haar::l2_star_discrepancy({pts, 1}), numeric_l2_star_1d(pts), 2e-3);
}

TEST(l2_star_of_centred_grid_is_small) {
    std::vector<double> pts(64);
    for (int i = 0; i < 64; ++i) pts[static_cast<std::size_t>(i)] = (2 * i + 1) / 128.0;
    // For the centred grid, (F_n(x) - x) is a sawtooth of amplitude 1/(2n): n * sqrt(1/(12 n^2)) = 1/sqrt(12).
    CHECK_CLOSE(haar::l2_star_discrepancy({pts, 1}), 1.0 / std::sqrt(12.0), 1e-9);
}

TEST(mean_of_constant_is_constant) {
    std::vector<double> pts(20, 0.3);
    CHECK_CLOSE(haar::mean({pts, 2}, [](std::span<const double>) { return 7.0; }), 7.0, 1e-15);
}
