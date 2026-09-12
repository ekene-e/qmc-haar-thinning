// SPDX-License-Identifier: MIT
#pragma once

/// @file discrepancy.hpp
/// Quality measures for point sets in [0,1)^d.
///
/// All quantities are *unnormalised*, matching the paper's convention
/// D_R(A) = |A ∩ R| - |A| |R|: divide by n for the usual normalised form.
///
///   star_discrepancy       exact, d = 1 (O(n log n)) and d = 2 (O(n^2));
///   l2_star_discrepancy    Warnock's closed form, any d, O(n^2 d);
///   mean                   sample mean of a function over the point set.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

namespace haar {

/// Row-major view of n points in d dimensions.
struct PointsView {
    std::span<const double> data;
    unsigned dim;

    [[nodiscard]] std::size_t size() const noexcept { return dim == 0 ? 0 : data.size() / dim; }
    [[nodiscard]] std::span<const double> operator[](std::size_t i) const noexcept {
        return data.subspan(i * dim, dim);
    }
};

namespace detail {

/// D*(A) = max over closed and half-open anchored intervals.  The extrema are
/// attained at point coordinates: (count <= x_k) - n x_k and n x_k - (count < x_k).
inline double star_discrepancy_1d(std::span<const double> xs_in) {
    std::vector<double> xs(xs_in.begin(), xs_in.end());
    std::sort(xs.begin(), xs.end());
    const double n = static_cast<double>(xs.size());
    double best = 0.0;
    for (std::size_t k = 0; k < xs.size(); ++k) {
        const bool last_of_run = k + 1 == xs.size() || xs[k + 1] != xs[k];
        const bool first_of_run = k == 0 || xs[k - 1] != xs[k];
        if (last_of_run) best = std::max(best, static_cast<double>(k + 1) - n * xs[k]);
        if (first_of_run) best = std::max(best, n * xs[k] - static_cast<double>(k));
    }
    return best;
}

/// Sweep over x-boundaries in increasing order.  For boundary xb, `open`
/// holds the y's of points with x < xb and `closed` those with x <= xb, both
/// sorted.  Closed boxes [0,xb]x[0,y] maximise count - n*vol, half-open boxes
/// [0,xb)x[0,y) maximise n*vol - count.
inline double star_discrepancy_2d(PointsView pts) {
    const std::size_t n = pts.size();
    struct P {
        double x, y;
    };
    std::vector<P> ps(n);
    for (std::size_t i = 0; i < n; ++i) ps[i] = {pts[i][0], pts[i][1]};
    std::sort(ps.begin(), ps.end(), [](const P& a, const P& b) { return a.x < b.x; });

    const double nd = static_cast<double>(n);
    double best = 0.0;
    auto evaluate = [&](double xb, const std::vector<double>& open, const std::vector<double>& closed) {
        for (std::size_t k = 0; k < closed.size(); ++k) {
            if (k + 1 < closed.size() && closed[k + 1] == closed[k]) continue;  // use count of y <= value
            best = std::max(best, static_cast<double>(k + 1) - nd * xb * closed[k]);
        }
        best = std::max(best, static_cast<double>(closed.size()) - nd * xb);
        for (std::size_t k = 0; k < open.size(); ++k) {
            if (k > 0 && open[k - 1] == open[k]) continue;  // use count of y < value
            best = std::max(best, nd * xb * open[k] - static_cast<double>(k));
        }
        best = std::max(best, nd * xb - static_cast<double>(open.size()));
    };

    std::vector<double> closed, open;
    closed.reserve(n);
    open.reserve(n);
    std::size_t i = 0;
    while (i < n) {
        const double xb = ps[i].x;
        open = closed;
        while (i < n && ps[i].x == xb) {
            closed.insert(std::upper_bound(closed.begin(), closed.end(), ps[i].y), ps[i].y);
            ++i;
        }
        evaluate(xb, open, closed);
    }
    evaluate(1.0, closed, closed);
    return best;
}

}  // namespace detail

/// Exact star discrepancy sup_R | |A ∩ R| - n |R| | over anchored boxes.
/// Implemented for d = 1 and d = 2; returns nullopt otherwise.
[[nodiscard]] inline std::optional<double> star_discrepancy(PointsView pts) {
    if (pts.dim == 1) return detail::star_discrepancy_1d(pts.data);
    if (pts.dim == 2) return detail::star_discrepancy_2d(pts);
    return std::nullopt;
}

/// Unnormalised L2 star discrepancy n * T_n(A), where (Warnock 1972)
///   T_n^2 = 3^-d - 2^(1-d)/n sum_i prod_k (1 - x_ik^2) + 1/n^2 sum_{i,j} prod_k (1 - max(x_ik, x_jk)).
[[nodiscard]] inline double l2_star_discrepancy(PointsView pts) {
    const std::size_t n = pts.size();
    const unsigned d = pts.dim;
    if (n == 0 || d == 0) return 0.0;
    double single = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double prod = 1.0;
        for (double c : pts[i]) prod *= 1.0 - c * c;
        single += prod;
    }
    double pair = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const auto xi = pts[i];
        double diag = 1.0;
        for (double c : xi) diag *= 1.0 - c;
        pair += diag;
        double row = 0.0;
        for (std::size_t j = 0; j < i; ++j) {
            const auto xj = pts[j];
            double prod = 1.0;
            for (unsigned k = 0; k < d; ++k) prod *= 1.0 - std::max(xi[k], xj[k]);
            row += prod;
        }
        pair += 2.0 * row;
    }
    const double nd = static_cast<double>(n);
    const double t2 = std::pow(3.0, -static_cast<double>(d)) - std::ldexp(1.0, 1 - static_cast<int>(d)) / nd * single +
                      pair / (nd * nd);
    return nd * std::sqrt(std::max(t2, 0.0));
}

/// Sample mean of `f` over the point set (the QMC estimate of its integral).
template <class F>
[[nodiscard]] double mean(PointsView pts, F&& f) {
    const std::size_t n = pts.size();
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) sum += f(pts[i]);
    return n == 0 ? 0.0 : sum / static_cast<double>(n);
}

}  // namespace haar
