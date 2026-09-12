// SPDX-License-Identifier: MIT
#pragma once

/// @file thinner.hpp
/// Online (1 + eps)-thinning of i.i.d. uniform samples with Haar feedback.
///
/// Implements the algorithms of
///   E. Ezeunala, A. V. Jha, H. Jiang, "Quasi-Monte Carlo Beyond Hardy--Krause
///   II: (1 + eps) n Samples Suffice", building on
///   R. Dwivedi, O. N. Feldheim, O. Gurel-Gurevich, A. Ramdas, "The power of
///   online thinning in reducing discrepancy" (PTRF 2019).
///
/// Samples x_1, x_2, ... in [0,1)^d arrive one at a time.  Upon seeing a
/// sample the thinner either retains it or rejects it; after a rejection the
/// next sample is always retained.  Each candidate x is rejected with
/// probability r_t(x) = 1 + eps/2 - mu_t(x) in [0, eps], where mu_t is the
/// target density of the framework in Section 2.1 of the paper:
///
///   sign feedback   (eq. 3.1):  mu_t(x) = 1 + eps/(2N) sum_H sgn(-phi_t(H)) H(x)
///   linear feedback (eq. 4.1):  mu_t(x) = 1 - eps/(2B) sum_H phi_t(H) H(x)
///
/// where phi_t(H) = sum_{z retained} H(z) is the discrepancy of the Haar
/// function H, the sums range over all Haar functions in the configured scale
/// set, N is the number of scale vectors, and B is the saturation bound.
///
/// Everything on the per-sample path is integer arithmetic: coordinates are
/// converted to dyadic fixed point, discrepancies are int32 counters, and the
/// rejection coin is a 32-bit integer compared against an exact integer
/// threshold.  Given a seed and a sample stream, the output is therefore
/// identical on every platform.

#include <algorithm>
#include <bit>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <format>
#include <initializer_list>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "haar/detail/geometry.hpp"
#include "haar/detail/table.hpp"
#include "haar/options.hpp"
#include "haar/rng.hpp"

namespace haar {

/// Counters describing a run.
struct Stats {
    /// Samples presented through `offer` (including forced acceptances).
    std::uint64_t offered = 0;
    /// Samples retained.  Always `offered - rejected`.
    std::uint64_t retained = 0;
    /// Samples rejected.  At most about eps/2 * offered.
    std::uint64_t rejected = 0;
    /// Linear feedback only: decisions where |Phi_t(x)| exceeded B.
    std::uint64_t saturated = 0;
    /// Times the scale set was enlarged and the table rebuilt (adaptive level).
    std::uint64_t rebuilds = 0;
    /// Linear feedback only: largest |Phi_t(x)| seen at a decision.
    std::int64_t max_abs_field = 0;
};

/// Online Haar-thinning strategy.  See the file comment for the algorithm.
///
/// Typical use:
/// @code
///   haar::Options opt;  opt.dim = 2;  opt.expected_n = 100000;
///   haar::Thinner thinner(opt);
///   thinner.run(100000);                 // draws its own uniform samples
///   auto pts = thinner.points();         // row-major, 2 doubles per point
/// @endcode
/// or, feeding your own stream, `while (thinner.size() < n) thinner.offer(x);`.
class Thinner {
public:
    explicit Thinner(Options options) : options_(std::move(options)) {
        validate_options(options_);
        dim_ = options_.dim;
        rng_.reseed(options_.seed);
        // eps in 32-bit fixed point; the decision compares 32-bit coins.
        epsilon_fixed_ = static_cast<std::uint64_t>(std::llround(std::ldexp(options_.epsilon, 32)));
        epsilon_fixed_ = std::clamp<std::uint64_t>(epsilon_fixed_, 1, std::uint64_t{1} << 32);
        shift_.assign(dim_, 0.0);
        if (options_.shift == Shift::random) {
            for (auto& s : shift_) s = rng_.uniform();
        } else if (options_.shift == Shift::fixed) {
            shift_ = options_.shift_vector;
        }
        fixed_level_ = options_.level.has_value() || options_.expected_n.has_value();
        configure(options_.expected_n.value_or(1));
    }

    // ------------------------------------------------------------------ online API

    /// Presents one sample.  Returns true if it was retained.  After a
    /// rejection the next call always retains (and returns true).
    /// Throws `std::invalid_argument` if `x` has the wrong size or a
    /// coordinate outside [0, 1).
    bool offer(std::span<const double> x) {
        if (x.size() != dim_) throw std::invalid_argument("Thinner::offer: wrong dimension");
        if (stats_.retained >= kMaxPoints) throw std::length_error("Thinner: at most 2^31 points");
        ++stats_.offered;
        if (!fixed_level_ && target_level(stats_.retained + 1) != target_level_) configure(stats_.retained + 1);
        prepare_point(x, /*validate=*/true);
        // The vote pass also records, for every touched Haar function, the
        // table slot that a subsequent commit will update.
        const std::int64_t statistic = std::visit([&](auto& t) { return vote(t); }, table_);
        if (force_accept_) {
            force_accept_ = false;
        } else if (!decide(statistic)) {
            force_accept_ = true;
            ++stats_.rejected;
            return false;
        }
        std::visit([&](auto& t) { commit(t); }, table_);
        points_.insert(points_.end(), x.begin(), x.end());
        ++stats_.retained;
        return true;
    }

    bool offer(std::initializer_list<double> x) { return offer(std::span<const double>(x.begin(), x.size())); }

    /// Draws uniform samples from `gen` and offers them until `n` points are
    /// retained in total.  Returns the number of samples consumed by this call.
    template <std::uniform_random_bit_generator G>
    std::uint64_t run(std::size_t n, G& gen) {
        const std::uint64_t before = stats_.offered;
        std::vector<double> x(dim_);
        while (stats_.retained < n) {
            for (auto& c : x) c = detail::uniform01(gen);
            offer(x);
        }
        return stats_.offered - before;
    }

    /// Same as `run(n, gen)` with the internal seeded generator.
    std::uint64_t run(std::size_t n) { return run(n, rng_); }

    // ------------------------------------------------------------------ results

    [[nodiscard]] unsigned dim() const noexcept { return dim_; }

    /// Number of retained points.
    [[nodiscard]] std::size_t size() const noexcept { return static_cast<std::size_t>(stats_.retained); }

    /// All retained points in arrival order, row-major (`dim()` doubles each).
    [[nodiscard]] std::span<const double> points() const noexcept { return points_; }

    /// The i-th retained point.
    [[nodiscard]] std::span<const double> point(std::size_t i) const noexcept {
        return std::span<const double>(points_).subspan(i * dim_, dim_);
    }

    /// Moves the retained points out and resets the thinner to its initial
    /// state (same options; the random generator continues, so a fresh shift
    /// is drawn if `Shift::random`).
    [[nodiscard]] std::vector<double> release() {
        std::vector<double> out = std::move(points_);
        points_.clear();
        stats_ = Stats{};
        force_accept_ = false;
        if (options_.shift == Shift::random) {
            for (auto& s : shift_) s = rng_.uniform();
        }
        configured_ = false;  // force a fresh, empty table
        configure(options_.expected_n.value_or(1));
        return out;
    }

    [[nodiscard]] const Stats& stats() const noexcept { return stats_; }
    [[nodiscard]] const Options& options() const noexcept { return options_; }

    /// The shift s applied to the Haar system (all zeros for `Shift::none`).
    [[nodiscard]] std::span<const double> shift() const noexcept { return shift_; }

    /// Resolution parameter in effect (L or l, see `Options::level`).
    [[nodiscard]] unsigned level() const noexcept { return level_; }
    /// Resolution the automatic rule asked for before applying the work and
    /// memory caps (equal to `level()` unless capped).
    [[nodiscard]] unsigned target_level() const noexcept { return target_level_; }
    /// Number of scale vectors N (the normaliser of the sign rule).
    [[nodiscard]] std::uint64_t scale_count() const noexcept { return geometry_.scale_count; }
    /// Haar functions evaluated per sample (N - 1).
    [[nodiscard]] std::uint64_t function_count() const noexcept { return geometry_.function_count(); }
    /// Number of discrepancy counters addressable by the current scale set.
    [[nodiscard]] std::uint64_t cell_count() const noexcept { return geometry_.cell_count; }
    /// Whether counters are stored densely.
    [[nodiscard]] bool dense() const noexcept { return std::holds_alternative<detail::DenseTable>(table_); }
    /// Bytes held by the counter table.
    [[nodiscard]] std::size_t table_bytes() const noexcept {
        return std::visit([](const auto& t) { return t.bytes(); }, table_);
    }
    /// The bound B in effect (linear feedback), 0 for sign feedback.
    [[nodiscard]] std::int64_t saturation_bound() const noexcept { return bound_; }

    /// max_H |phi_t(H)| over all Haar functions in the scale set (a full
    /// table scan; for diagnostics).
    [[nodiscard]] std::int64_t max_abs_haar_discrepancy() const {
        std::int64_t best = 0;
        std::visit([&](const auto& t) {
            t.for_each_nonzero([&](std::uint64_t, std::int32_t v) {
                best = std::max<std::int64_t>(best, v < 0 ? -std::int64_t{v} : v);
            });
        }, table_);
        return best;
    }

    /// Calls fn(key, phi) for every Haar function with nonzero discrepancy.
    template <class Fn>
    void for_each_haar_discrepancy(Fn&& fn) const {
        std::visit([&](const auto& t) { t.for_each_nonzero(fn); }, table_);
    }

    /// One-line human-readable summary of the effective configuration.
    [[nodiscard]] std::string describe() const {
        std::string level_text = std::format("level={}", level_);
        if (!options_.level && level_ != target_level_) level_text += std::format(" (capped from {})", target_level_);
        if (!fixed_level_) level_text += " (adaptive)";
        return std::format(
            "dim={} eps={:.6g} feedback={} shift={} scale_set={} {} N={} cells={} storage={} ({:.1f} MiB){}",
            dim_, options_.epsilon, options_.feedback == Feedback::sign ? "sign" : "linear",
            options_.shift == Shift::none ? "none" : "on",
            options_.scale_set == ScaleSet::box ? "box" : "hyperbolic", level_text, geometry_.scale_count,
            geometry_.cell_count, dense() ? "dense" : "sparse", static_cast<double>(table_bytes()) / (1u << 20),
            options_.feedback == Feedback::linear ? std::format(" B={}", bound_) : std::string{});
    }

private:
    static constexpr std::uint64_t kMaxPoints = std::uint64_t{1} << 31;
    static constexpr std::uint64_t kMaxScaleCount = std::uint64_t{1} << 30;

    // ------------------------------------------------------------------ level policy

    /// Uncapped resolution for a run of about `n` points:
    ///   hyperbolic  L = ceil(log2 n) + offset,   box  l = ceil((log2 n + offset) / d),
    /// clamped to [1, 62].  Both make the finest cells hold about 2^-offset points.
    [[nodiscard]] unsigned target_level(std::uint64_t n) const noexcept {
        if (options_.level) return *options_.level;
        const int lg = n <= 1 ? 0 : std::bit_width(n - 1);  // ceil(log2 n)
        int level = lg + options_.level_offset;
        if (options_.scale_set == ScaleSet::box) {
            level = level <= 0 ? 0 : (level + static_cast<int>(dim_) - 1) / static_cast<int>(dim_);
        }
        return static_cast<unsigned>(std::clamp(level, 1, 62));
    }

    /// Geometry for a given level of the configured scale set.
    [[nodiscard]] detail::Geometry make_geometry(unsigned level) const {
        const bool box = options_.scale_set == ScaleSet::box;
        return detail::Geometry::make(dim_, level, box ? dim_ * level : level);
    }

    /// Closed-form number of scale vectors at `level`, as a double to avoid overflow.
    [[nodiscard]] double scale_count_estimate(unsigned level) const noexcept {
        if (options_.scale_set == ScaleSet::box) return std::pow(level + 1.0, static_cast<double>(dim_));
        double c = 1.0;  // C(level + d, d)
        for (unsigned i = 1; i <= dim_; ++i) c = c * (level + static_cast<double>(i)) / static_cast<double>(i);
        return c;
    }

    /// Largest level <= `target` whose scale set respects the work cap and
    /// whose dense table fits the memory budget.  An explicit `Options::level`
    /// is returned as is.
    [[nodiscard]] std::pair<unsigned, detail::Geometry> select_geometry(unsigned target) const {
        if (options_.level) return {*options_.level, make_geometry(*options_.level)};
        for (unsigned level = target; level >= 1; --level) {
            if (scale_count_estimate(level) > static_cast<double>(options_.max_functions_per_sample)) continue;
            try {
                detail::Geometry g = make_geometry(level);
                if (g.cell_count * sizeof(std::int32_t) <= options_.dense_budget_bytes || level == 1) {
                    return {level, std::move(g)};
                }
            } catch (const std::length_error&) {
                // too large for 64-bit keys: try a coarser level
            }
        }
        return {1u, make_geometry(1)};
    }

    /// (Re)configures the resolution for a run of about `n_hint` points:
    /// picks the level, rebuilds the table if the level changed, and refreshes
    /// the saturation bound.
    void configure(std::uint64_t n_hint) {
        target_level_ = target_level(n_hint);
        // Horizon: the number of points this configuration is meant to serve,
        // used by the automatic saturation bound.
        if (options_.expected_n) {
            horizon_ = static_cast<double>(*options_.expected_n);
        } else {
            const int bits = options_.scale_set == ScaleSet::box ? static_cast<int>(dim_ * target_level_)
                                                                 : static_cast<int>(target_level_);
            horizon_ = std::ldexp(1.0, std::clamp(bits - options_.level_offset, 0, 62));
        }
        // configure() runs only when the target changes (at most once per
        // doubling of the retained count), so building a candidate geometry
        // here is cheap.  The table is rebuilt only if the capped level moved.
        auto [level, geometry] = select_geometry(target_level_);
        if (!configured_ || level != level_) rebuild(level, std::move(geometry));
        bound_ = options_.feedback == Feedback::linear
                     ? options_.saturation_bound.value_or(automatic_bound())
                     : 0;
    }

    /// Installs a new scale set, allocates the table and re-inserts all
    /// retained points.
    void rebuild(unsigned level, detail::Geometry geometry) {
        if (geometry.scale_count > kMaxScaleCount) throw std::length_error("Thinner: too many scale vectors");
        const bool first = !configured_;
        level_ = level;
        geometry_ = std::move(geometry);
        configured_ = true;

        const unsigned per_coordinate = geometry_.per_coordinate_level;
        scale_ = std::ldexp(1.0, static_cast<int>(per_coordinate));
        const std::size_t m = static_cast<std::size_t>(geometry_.function_count());
        const std::size_t stride = per_coordinate + 1;
        cell_index_.assign(dim_ * stride, 0);
        parity_.assign(dim_ * stride, 0);
        keys_.assign(m, 0);
        signs_.assign(m, 0);
        slots_.assign(m, 0);

        const bool use_dense = options_.storage == Storage::dense ||
                               (options_.storage == Storage::automatic &&
                                geometry_.cell_count * sizeof(std::int32_t) <= options_.dense_budget_bytes);
        if (use_dense) {
            table_ = detail::DenseTable(geometry_.cell_count);
        } else {
            table_ = detail::SparseTable();
        }

        const std::size_t count = static_cast<std::size_t>(stats_.retained);
        for (std::size_t i = 0; i < count; ++i) {
            prepare_point(std::span<const double>(points_).subspan(i * dim_, dim_), /*validate=*/false);
            std::visit([&](auto& t) { (void)vote(t); commit(t); }, table_);
        }
        if (!first && count > 0) ++stats_.rebuilds;
    }

    /// Automatic saturation bound.
    ///
    /// Under linear feedback each discrepancy phi(H) is a random walk with
    /// increments of variance ||H||_2^2 and a restoring drift eps ||H||^2/(2B),
    /// so after t steps its variance is about min(B/eps, t ||H||_2^2): coarse
    /// functions have equilibrated, fine ones have not.  The field
    /// Phi(x) = sum_H phi(H) H(x) sums one such term per scale vector, hence
    /// Var Phi ~ V(B) = sum_H min(B/eps, t ||H||^2) with ||H||^2 = 1/(block size).
    /// Saturation (|Phi| > B) is rare once B >= kappa sqrt(V(B)); we take
    /// kappa = 4.5 and solve the fixed point by iteration.  The result is
    /// Theta(log^d n / eps) up to log factors, the same order as the bound of
    /// the paper with the union-bound logarithm replaced by the constant kappa.
    [[nodiscard]] std::int64_t automatic_bound() const noexcept {
        constexpr double kappa = 4.5;
        const double eps = options_.epsilon;
        // Histogram of scale vectors by log2(block size) = -log2 ||H||^2.
        std::uint64_t histogram[64] = {};
        const auto& offsets = geometry_.offsets;
        for (std::size_t s = 1; s < offsets.size(); ++s) {  // skip the constant function
            const std::uint64_t block = (s + 1 < offsets.size() ? offsets[s + 1] : geometry_.cell_count) - offsets[s];
            ++histogram[std::countr_zero(block)];
        }
        auto variance = [&](double b) {
            double v = 0.0;
            for (int bits = 0; bits < 64; ++bits) {
                if (histogram[bits]) v += static_cast<double>(histogram[bits]) * std::min(b / eps, std::ldexp(horizon_, -bits));
            }
            return v;
        };
        double b = 1.0;
        for (int iteration = 0; iteration < 200; ++iteration) {
            const double next = std::max(1.0, kappa * std::sqrt(variance(b)));
            const bool converged = std::fabs(next - b) < 0.5;
            b = next;
            if (converged) break;
        }
        return static_cast<std::int64_t>(std::min(std::ceil(b), std::ldexp(1.0, 30)));
    }

    // ------------------------------------------------------------------ hot path

    /// Shifts `x`, converts to dyadic fixed point and enumerates the keys and
    /// signs of all Haar functions that are nonzero at it.
    void prepare_point(std::span<const double> x, bool validate) {
        const unsigned R = geometry_.per_coordinate_level;
        const std::size_t stride = R + 1;
        for (unsigned i = 0; i < dim_; ++i) {
            const double xi = x[i];
            if (validate && !(xi >= 0.0 && xi < 1.0)) {
                throw std::invalid_argument("Thinner::offer: coordinates must lie in [0, 1)");
            }
            double y = xi - shift_[i];
            if (y < 0.0) y += 1.0;
            if (y >= 1.0) y = 0.0;  // guards x - s + 1 rounding up to 1
            // u = floor(y 2^R): exact, since scaling by a power of two is exact.
            const std::uint64_t u = static_cast<std::uint64_t>(y * scale_);
            std::uint64_t* k = cell_index_.data() + i * stride;
            std::uint8_t* p = parity_.data() + i * stride;
            k[0] = 0;
            p[0] = 0;
            for (unsigned j = 1; j <= R; ++j) {
                // h_{j,k} is supported on [k 2^{1-j}, (k+1) 2^{1-j}), k = floor(y 2^{j-1}),
                // and is +1 on the left half (bit R-j of u clear), -1 on the right.
                k[j] = u >> (R + 1 - j);
                p[j] = static_cast<std::uint8_t>((u >> (R - j)) & 1u);
            }
        }
        enumerate_keys();
    }

    /// Fills keys_/signs_ for the point prepared by `prepare_point`.
    void enumerate_keys() noexcept {
        const unsigned R = geometry_.per_coordinate_level;
        const std::size_t stride = R + 1;
        const std::size_t prefix_len = dim_ - 1;
        const std::uint64_t* const k_last = cell_index_.data() + prefix_len * stride;
        const std::uint8_t* const p_last = parity_.data() + prefix_len * stride;
        const std::uint64_t* offsets = geometry_.offsets.data();
        std::uint64_t* keys = keys_.data();
        std::uint8_t* signs = signs_.data();
        std::size_t out = 0;

        const std::size_t prefixes = geometry_.prefix_count();
        for (std::size_t p = 0; p < prefixes; ++p) {
            // Fold the first d-1 coordinates into the block-local cell index.
            std::uint64_t idx = 0;
            std::uint8_t parity = 0;
            const std::uint8_t* levels = geometry_.prefix_levels.data() + p * prefix_len;
            for (std::size_t i = 0; i < prefix_len; ++i) {
                const unsigned lv = levels[i];
                idx = (idx << detail::cell_shift(lv)) | cell_index_[i * stride + lv];
                parity ^= parity_[i * stride + lv];
            }
            const std::uint64_t* block = offsets + geometry_.prefix_positions[p];
            const unsigned bound = geometry_.prefix_bounds[p];
            // Prefix 0 is the all-zero prefix; its level-0 entry is the constant
            // function, which casts no vote and is skipped.
            for (unsigned j = (p == 0 ? 1u : 0u); j <= bound; ++j) {
                keys[out] = block[j] + ((idx << detail::cell_shift(j)) | k_last[j]);
                signs[out] = static_cast<std::uint8_t>(parity ^ p_last[j]);
                ++out;
            }
        }
    }

    /// Reads the discrepancies of all touched Haar functions, recording their
    /// slots.  Returns V = sum_H sgn(phi(H)) H(x) for sign feedback or
    /// Phi = sum_H phi(H) H(x) for linear feedback.
    template <class Table>
    [[nodiscard]] std::int64_t vote(Table& table) {
        const std::size_t m = keys_.size();
        table.reserve_extra(m);
        const std::uint64_t* keys = keys_.data();
        const std::uint8_t* signs = signs_.data();
        std::uint64_t* slots = slots_.data();
        for (std::size_t i = 0; i < m; ++i) table.prefetch(keys[i]);

        std::int64_t acc = 0;
        if (options_.feedback == Feedback::sign) {
            for (std::size_t i = 0; i < m; ++i) {
                const auto slot = table.locate(keys[i]);
                if constexpr (!Table::kIdentitySlots) slots[i] = slot;
                const std::int32_t phi = table.value(slot);
                const std::int32_t s = (phi > 0) - (phi < 0);
                const std::int32_t b = signs[i];  // 1 if H(x) = -1
                acc += (s ^ -b) + b;              // s * H(x), branch-free
            }
        } else {
            for (std::size_t i = 0; i < m; ++i) {
                const auto slot = table.locate(keys[i]);
                if constexpr (!Table::kIdentitySlots) slots[i] = slot;
                const std::int64_t phi = table.value(slot);
                const std::int64_t b = signs[i];
                acc += (phi ^ -b) + b;
            }
        }
        return acc;
    }

    /// Flips the rejection coin.  Returns true to retain the candidate.
    ///
    /// With E = round(eps 2^32) and a uniform 32-bit X, the candidate is
    /// rejected when X < r 2^32, i.e.
    ///   sign:    2 N X < E (N + V)      (r = eps/2 (1 + V/N))
    ///   linear:  2 B X < E (B + Phi)    (r = eps/2 (1 + Phi/B)), |Phi| <= B.
    /// All products fit in 64 bits because N, B <= 2^30.
    [[nodiscard]] bool decide(std::int64_t statistic) {
        const std::uint64_t coin = rng_.uniform32();
        if (options_.feedback == Feedback::sign) {
            const std::uint64_t n_scales = geometry_.scale_count;
            const std::uint64_t lhs = 2 * n_scales * coin;
            const std::uint64_t rhs = epsilon_fixed_ * static_cast<std::uint64_t>(static_cast<std::int64_t>(n_scales) + statistic);
            return lhs >= rhs;
        }
        std::int64_t field = statistic;
        const std::int64_t magnitude = field < 0 ? -field : field;
        stats_.max_abs_field = std::max(stats_.max_abs_field, magnitude);
        if (magnitude > bound_) {
            ++stats_.saturated;
            if (options_.on_saturation == OnSaturation::fail) {
                throw SaturationError(std::format("Haar field |Phi| = {} exceeds bound B = {}", magnitude, bound_));
            }
            field = field < 0 ? -bound_ : bound_;
        }
        const std::uint64_t lhs = 2 * static_cast<std::uint64_t>(bound_) * coin;
        const std::uint64_t rhs = epsilon_fixed_ * static_cast<std::uint64_t>(bound_ + field);
        return lhs >= rhs;
    }

    /// Adds H(x) to phi(H) for every touched Haar function.
    template <class Table>
    void commit(Table& table) noexcept {
        const std::size_t m = keys_.size();
        const std::uint64_t* keys = keys_.data();
        const std::uint8_t* signs = signs_.data();
        const std::uint64_t* slots = Table::kIdentitySlots ? keys : slots_.data();
        for (std::size_t i = 0; i < m; ++i) {
            table.add(slots[i], keys[i], 1 - 2 * static_cast<std::int32_t>(signs[i]));
        }
    }

    // ------------------------------------------------------------------ state

    Options options_;
    unsigned dim_ = 1;
    Xoshiro256pp rng_;
    std::uint64_t epsilon_fixed_ = 0;
    std::vector<double> shift_;
    bool fixed_level_ = false;
    bool configured_ = false;
    unsigned level_ = 0;
    unsigned target_level_ = 0;
    double horizon_ = 1.0;
    std::int64_t bound_ = 0;
    double scale_ = 1.0;
    detail::Geometry geometry_;
    std::variant<detail::DenseTable, detail::SparseTable> table_;

    // Per-sample scratch (sized once per rebuild, never reallocated per sample).
    std::vector<std::uint64_t> cell_index_;  // [dim][R+1]: k_i(j)
    std::vector<std::uint8_t> parity_;       // [dim][R+1]: 1 if h_{j,k}(y_i) = -1
    std::vector<std::uint64_t> keys_;        // touched cells, one per scale vector
    std::vector<std::uint8_t> signs_;        // 1 if H(x) = -1
    std::vector<std::uint64_t> slots_;       // table slots recorded by vote()

    std::vector<double> points_;
    Stats stats_;
    bool force_accept_ = false;
};

/// Convenience: retain `n` points with the given options, drawing samples
/// from the seeded internal generator.  Sets `expected_n = n` if unset.
[[nodiscard]] inline std::vector<double> thin(std::size_t n, Options options) {
    if (!options.expected_n) options.expected_n = n;
    Thinner thinner(std::move(options));
    thinner.run(n);
    return thinner.release();
}

}  // namespace haar
