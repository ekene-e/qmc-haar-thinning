// SPDX-License-Identifier: MIT
//
// End-to-end verification of the thinner against a deliberately naive
// reference implementation: the reference evaluates every Haar function of
// the scale set on every retained point at every step (no tables, no
// fixed-point tricks) and flips the same integer coin.  The two must agree on
// every retain/reject decision for all combinations of feedback rule, scale
// set, storage, shift and level policy.

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <map>
#include <vector>

#include "check.hpp"
#include "haar/discrepancy.hpp"
#include "haar/thinner.hpp"

namespace {

using haar::Feedback;
using haar::Options;
using haar::ScaleSet;
using haar::Shift;
using haar::Storage;
using haar::Thinner;

class Reference {
public:
    explicit Reference(const Options& opt) : opt_(opt), d_(opt.dim), coins_(opt.seed) {
        shift_.assign(d_, 0.0);
        if (opt.shift == Shift::random) {
            for (auto& s : shift_) s = coins_.uniform();
        } else if (opt.shift == Shift::fixed) {
            shift_ = opt.shift_vector;
        }
        eps_fixed_ = static_cast<std::uint64_t>(std::llround(std::ldexp(opt.epsilon, 32)));
        eps_fixed_ = std::clamp<std::uint64_t>(eps_fixed_, 1, std::uint64_t{1} << 32);
    }

    // `L` and `B` are the resolution and bound the thinner used for this step:
    // choosing them is configuration, the process given them is what we verify.
    bool offer(const std::vector<double>& x, unsigned L, std::int64_t B) {
        const unsigned R = L;
        const unsigned total = opt_.scale_set == ScaleSet::box ? d_ * L : L;
        const auto scales = all_scales(R, total);
        const std::uint64_t N = scales.size();

        std::vector<double> y = shifted(x);
        std::int64_t statistic = 0;
        for (const auto& j : scales) {
            if (std::all_of(j.begin(), j.end(), [](unsigned v) { return v == 0; })) continue;
            const auto [kx, hx] = evaluate(j, y);
            std::int64_t phi = 0;
            for (const auto& z : retained_) {
                const auto [kz, hz] = evaluate(j, shifted(z));
                if (kz == kx) phi += hz;
            }
            if (opt_.feedback == Feedback::sign) {
                statistic += ((phi > 0) - (phi < 0)) * hx;
            } else {
                statistic += phi * hx;
            }
        }

        if (forced_) {
            forced_ = false;
            retained_.push_back(x);
            return true;
        }
        const std::uint64_t coin = coins_.uniform32();
        bool accept;
        if (opt_.feedback == Feedback::sign) {
            accept = 2 * N * coin >= eps_fixed_ * static_cast<std::uint64_t>(static_cast<std::int64_t>(N) + statistic);
        } else {
            if (statistic > B) ++saturated_, statistic = B;
            if (statistic < -B) ++saturated_, statistic = -B;
            accept = 2 * static_cast<std::uint64_t>(B) * coin >= eps_fixed_ * static_cast<std::uint64_t>(B + statistic);
        }
        if (!accept) {
            forced_ = true;
            return false;
        }
        retained_.push_back(x);
        return true;
    }

    std::uint64_t saturated() const { return saturated_; }

private:
    std::vector<double> shifted(const std::vector<double>& x) const {
        std::vector<double> y(d_);
        for (unsigned i = 0; i < d_; ++i) {
            y[i] = x[i] - shift_[i];
            if (y[i] < 0) y[i] += 1.0;
            if (y[i] >= 1.0) y[i] = 0.0;
        }
        return y;
    }

    // Returns the position vector k of the unique H_{j,k} nonzero at y, and its value +-1.
    static std::pair<std::vector<std::uint64_t>, int> evaluate(const std::vector<unsigned>& j, const std::vector<double>& y) {
        std::vector<std::uint64_t> k(j.size(), 0);
        int sign = 1;
        for (std::size_t i = 0; i < j.size(); ++i) {
            if (j[i] == 0) continue;
            k[i] = static_cast<std::uint64_t>(std::floor(std::ldexp(y[i], static_cast<int>(j[i]) - 1)));
            const std::uint64_t half = static_cast<std::uint64_t>(std::floor(std::ldexp(y[i], static_cast<int>(j[i]))));
            if (half & 1) sign = -sign;
        }
        return {k, sign};
    }

    std::vector<std::vector<unsigned>> all_scales(unsigned R, unsigned total) const {
        std::vector<std::vector<unsigned>> out;
        std::vector<unsigned> j(d_, 0);
        for (;;) {
            unsigned sum = 0;
            for (unsigned v : j) sum += v;
            if (sum <= total) out.push_back(j);
            std::size_t i = 0;
            while (i < d_ && ++j[i] > R) j[i++] = 0;
            if (i == d_) break;
        }
        return out;
    }

    Options opt_;
    unsigned d_;
    haar::Xoshiro256pp coins_;
    std::vector<double> shift_;
    std::uint64_t eps_fixed_ = 0;
    std::vector<std::vector<double>> retained_;
    bool forced_ = false;
    std::uint64_t saturated_ = 0;
};

// Feeds the same external sample stream to both implementations.
void compare_with_reference(Options opt, std::size_t samples, std::uint64_t stream_seed) {
    Thinner thinner(opt);
    Reference reference(opt);
    haar::Xoshiro256pp stream(stream_seed);
    std::vector<double> x(opt.dim);
    std::size_t mismatches = 0;
    for (std::size_t s = 0; s < samples; ++s) {
        for (auto& c : x) c = stream.uniform();
        const bool a = thinner.offer(x);
        // level() and saturation_bound() still describe the configuration the
        // thinner used for this step (they only change at the next offer).
        const bool b = reference.offer(x, thinner.level(), thinner.saturation_bound());
        mismatches += a != b;
    }
    CHECK_EQ(mismatches, std::size_t{0});
    CHECK_EQ(thinner.stats().saturated, reference.saturated());
    CHECK_EQ(thinner.stats().offered, samples);
    CHECK_EQ(thinner.stats().retained + thinner.stats().rejected, samples);
}

}  // namespace

TEST(matches_reference_sign_box_dense) {
    Options o;
    o.dim = 2;
    o.feedback = Feedback::sign;
    o.scale_set = ScaleSet::box;
    o.level = 3;
    o.shift = Shift::none;
    o.storage = Storage::dense;
    o.epsilon = 0.7;
    compare_with_reference(o, 400, 1);
}

TEST(matches_reference_linear_hyperbolic_sparse_shift) {
    Options o;
    o.dim = 3;
    o.feedback = Feedback::linear;
    o.scale_set = ScaleSet::hyperbolic;
    o.level = 5;
    o.shift = Shift::random;
    o.storage = Storage::sparse;
    o.epsilon = 0.5;
    o.seed = 99;
    compare_with_reference(o, 300, 2);
}

TEST(matches_reference_linear_small_bound_saturates) {
    Options o;
    o.dim = 1;
    o.feedback = Feedback::linear;
    o.level = 6;
    o.shift = Shift::fixed;
    o.shift_vector = {0.37};
    o.saturation_bound = 3;  // deliberately tiny so clamping happens
    o.epsilon = 1.0;
    compare_with_reference(o, 500, 3);
}

TEST(matches_reference_adaptive_level_both_rules) {
    for (Feedback fb : {Feedback::sign, Feedback::linear}) {
        for (ScaleSet set : {ScaleSet::hyperbolic, ScaleSet::box}) {
            Options o;
            o.dim = 2;
            o.feedback = fb;
            o.scale_set = set;
            o.level_offset = 1;
            o.shift = Shift::random;
            o.seed = 5;
            compare_with_reference(o, 260, 4);  // several level increases and rebuilds
        }
    }
}

TEST(matches_reference_expected_n_fixes_level) {
    Options o;
    o.dim = 2;
    o.expected_n = 200;
    o.seed = 8;
    Thinner t(o);
    CHECK_EQ(t.level(), 9u);  // ceil(log2 200) + 1
    compare_with_reference(o, 250, 6);
}

TEST(never_rejects_twice_in_a_row_and_rate_below_epsilon) {
    Options o;
    o.dim = 2;
    o.epsilon = 0.3;
    o.expected_n = 20000;
    Thinner t(o);
    haar::Xoshiro256pp stream(42);
    std::vector<double> x(2);
    bool previous_rejected = false;
    while (t.size() < 20000) {
        for (auto& c : x) c = stream.uniform();
        const bool kept = t.offer(x);
        CHECK(!(previous_rejected && !kept));
        previous_rejected = !kept;
    }
    const double rate = static_cast<double>(t.stats().rejected) / static_cast<double>(t.stats().offered);
    CHECK(rate <= 0.3 * 0.5 * 1.15 + 0.01);  // rejection probability is eps/2 (1 + V/N) ~ eps/2 on average
    CHECK_EQ(t.stats().rejected + t.stats().retained, t.stats().offered);
    for (double c : t.points()) CHECK(c >= 0.0 && c < 1.0);
}

TEST(run_and_release_and_stats) {
    Options o;
    o.dim = 3;
    o.expected_n = 1000;
    Thinner t(o);
    const auto consumed = t.run(1000);
    CHECK_EQ(t.size(), std::size_t{1000});
    CHECK_EQ(consumed, t.stats().offered);
    CHECK_EQ(t.points().size(), std::size_t{3000});
    auto pts = t.release();
    CHECK_EQ(pts.size(), std::size_t{3000});
    CHECK_EQ(t.size(), std::size_t{0});
    CHECK_EQ(t.stats().offered, std::uint64_t{0});
    t.run(10);
    CHECK_EQ(t.size(), std::size_t{10});
}

TEST(rejects_bad_input) {
    Options o;
    o.dim = 2;
    Thinner t(o);
    CHECK_THROWS(t.offer({0.5}), std::invalid_argument);
    CHECK_THROWS(t.offer({0.5, 1.0}), std::invalid_argument);
    CHECK_THROWS(t.offer({-0.1, 0.2}), std::invalid_argument);
    CHECK_THROWS(t.offer({std::nan(""), 0.2}), std::invalid_argument);
    Options bad;
    bad.dim = 0;
    CHECK_THROWS(Thinner{bad}, std::invalid_argument);
    bad = Options{};
    bad.epsilon = 0.0;
    CHECK_THROWS(Thinner{bad}, std::invalid_argument);
    bad = Options{};
    bad.shift = Shift::fixed;
    bad.shift_vector = {1.5};
    CHECK_THROWS(Thinner{bad}, std::invalid_argument);
}

TEST(saturation_can_throw) {
    Options o;
    o.dim = 1;
    o.feedback = Feedback::linear;
    o.saturation_bound = 1;
    o.on_saturation = haar::OnSaturation::fail;
    o.level = 4;
    Thinner t(o);
    CHECK_THROWS(t.run(5000), haar::SaturationError);
}

TEST(haar_discrepancies_are_consistent_with_points) {
    // Sum of squared Haar discrepancies must equal what a direct recount gives.
    Options o;
    o.dim = 2;
    o.level = 4;
    o.shift = Shift::none;
    Thinner t(o);
    t.run(3000);
    std::int64_t from_table = 0;
    t.for_each_haar_discrepancy([&](std::uint64_t, std::int32_t v) { from_table += std::int64_t{v} * v; });
    // Direct recount of the level-(1,1) function: k = (0,0), value = product of signs of first bits.
    std::int64_t phi11 = 0;
    for (std::size_t i = 0; i < t.size(); ++i) {
        const auto p = t.point(i);
        phi11 += (p[0] < 0.5 ? 1 : -1) * (p[1] < 0.5 ? 1 : -1);
    }
    bool found = false;
    t.for_each_haar_discrepancy([&](std::uint64_t key, std::int32_t v) {
        // Lexicographic order over [0,4]^2 with sum <= 4: (0,0..4) occupy cells
        // 0..15 (block sizes 1,1,2,4,8), (1,0) is cell 16 and (1,1) is cell 17.
        if (key == 17) {
            found = true;
            CHECK_EQ(v, phi11);
        }
    });
    CHECK(found || phi11 == 0);
    CHECK(from_table > 0);
    CHECK(t.max_abs_haar_discrepancy() >= 1);
}

TEST(quality_beats_monte_carlo_in_1d) {
    // The star discrepancy of n i.i.d. points is the Kolmogorov statistic, with
    // mean 0.87 sqrt(n) = 222 at n = 65536.  Thinning stays polylogarithmic:
    // about 60 (linear) and 100 (sign) at this size.
    const std::size_t n = 65536;
    for (Feedback fb : {Feedback::sign, Feedback::linear}) {
        for (haar::Shift shift : {Shift::none, Shift::random}) {
            Options o;
            o.dim = 1;
            o.feedback = fb;
            o.shift = shift;
            o.seed = 2024;
            const auto pts = haar::thin(n, o);
            const double disc = *haar::star_discrepancy({pts, 1});
            CHECK(disc < (fb == Feedback::linear ? 0.4 : 0.6) * std::sqrt(static_cast<double>(n)));
        }
    }
}

TEST(quality_beats_monte_carlo_in_2d) {
    // L2-star discrepancy: E[n T_n] = sqrt(n (1/4 - 1/9)) = 0.373 sqrt(n) for
    // i.i.d. points.  At n = 16384 linear feedback lands at 0.45-0.7 of that;
    // the sign rule is only at Monte Carlo level at this size (its Haar
    // discrepancies scale like N/eps), so we merely check it is not worse than
    // a comfortable multiple.
    const std::size_t n = 16384;
    const double mc_expected = std::sqrt(static_cast<double>(n) * (0.25 - 1.0 / 9.0));
    for (Feedback fb : {Feedback::sign, Feedback::linear}) {
        Options o;
        o.dim = 2;
        o.feedback = fb;
        o.seed = 77;
        const auto pts = haar::thin(n, o);
        const double l2 = haar::l2_star_discrepancy({pts, 2});
        CHECK(l2 < (fb == Feedback::linear ? 0.85 : 1.5) * mc_expected);
    }
    // The greedy limit B = 1 (reject with probability eps whenever the field is
    // positive) is the strongest feedback and should do clearly better.
    Options greedy;
    greedy.dim = 2;
    greedy.saturation_bound = 1;
    greedy.seed = 77;
    const auto pts = haar::thin(n, greedy);
    CHECK(haar::l2_star_discrepancy({pts, 2}) < 0.6 * mc_expected);
}
