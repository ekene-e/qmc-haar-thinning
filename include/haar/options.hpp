// SPDX-License-Identifier: MIT
#pragma once

/// @file options.hpp
/// Configuration of a Haar-thinning run and its validation.

#include <cstddef>
#include <cstdint>
#include <expected>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace haar {

/// Which vote each Haar function casts when a candidate sample arrives.
enum class Feedback {
    /// Haar-thinning of Dwivedi, Feldheim, Gurel-Gurevich and Ramdas (2019),
    /// eq. (3.1) of the paper: every Haar function nonzero at the candidate
    /// votes +-1 according to the *sign* of its current discrepancy.
    /// Combined with `shift = Shift::random` this is the uniformly-shifted
    /// Haar-thinning algorithm of Theorem 1.2.
    sign,
    /// Linear-feedback Haar-thinning, eq. (4.1) of the paper: each vote is
    /// weighted by the discrepancy itself, giving a restoring drift
    /// proportional to the discrepancy.  This is the algorithm of
    /// Theorem 1.3 (star discrepancy O(log^(d+1) n)).
    linear,
};

/// The set of Haar scale vectors j in Z^d whose functions cast votes.
enum class ScaleSet {
    /// All j with j_1 + ... + j_d <= level.  Dyadic boxes of volume >= 2^-level
    /// (up to the constant factor 2^d).  Storage grows like n * log^(d-1) n, so
    /// this is the practical choice; it is the resolution used by DFG+19.
    hyperbolic,
    /// All j with max_i j_i <= level: the set Pi_{<= l} used in the paper's
    /// analysis.  Storage grows like 2^(d * level), so `level` must stay small
    /// (the automatic rule keeps 2^(d * level) close to n).
    box,
};

/// Whether the Haar system is shifted by a uniform random vector s in [0,1)^d
/// (Section 3.1 of the paper).  The shift is what upgrades the guarantee from
/// Koksma--Hlawka to the smoothed-out variation bound of Theorem 1.2.
enum class Shift {
    none,
    random,
    /// Use `Options::shift_vector` verbatim.
    fixed,
};

/// Storage of the Haar-discrepancy counters.
enum class Storage {
    /// Dense when the table fits into `dense_budget_bytes`, otherwise sparse.
    automatic,
    dense,
    sparse,
};

/// What to do when the linear-feedback field |Phi_t(x)| exceeds the bound B
/// (Stopping Condition 4.2 in the paper).
enum class OnSaturation {
    /// Clamp the target density into [1 - eps/2, 1 + eps/2] and count the
    /// event in `Stats::saturated`.  The run stays a valid (1+eps)-thinning
    /// strategy; only the analysis of Section 4 assumes this never happens.
    clamp,
    /// Throw `haar::SaturationError`.
    fail,
};

/// Full configuration of a `Thinner`.  Every field has a sensible default;
/// only `dim` normally needs to be set.
struct Options {
    /// Dimension d >= 1 of the unit cube.
    unsigned dim = 1;

    /// Over-sampling parameter eps in (0, 1]: each sample is rejected with
    /// probability at most eps, so about (1 + eps) n samples produce n points.
    /// Stored internally with 32-bit fixed-point precision.
    double epsilon = 0.5;

    Feedback feedback = Feedback::linear;
    Shift shift = Shift::random;

    /// Explicit shift, used only when `shift == Shift::fixed`.  Must have
    /// `dim` entries in [0, 1).
    std::vector<double> shift_vector{};

    ScaleSet scale_set = ScaleSet::hyperbolic;

    /// Resolution parameter of the scale set (L for hyperbolic, l for box).
    /// Unset: chosen automatically from `expected_n`, or adaptively from the
    /// number of points retained so far (see `level_offset`,
    /// `max_functions_per_sample` and `dense_budget_bytes`).  An explicit
    /// value is used as given, without caps.
    std::optional<unsigned> level{};

    /// Adjusts the automatic level rule.  With `hyperbolic`,
    /// level = ceil(log2 n) + level_offset; with `box`,
    /// level = ceil((log2 n + level_offset) / d).  Larger values track finer
    /// scales at higher cost.
    int level_offset = 1;

    /// Cap on the work per sample used by the automatic level rule: the level
    /// is lowered until the number of scale vectors (= Haar functions
    /// evaluated per sample, up to one) is at most this.  Work per sample is
    /// Theta(log^d n) for the uncapped rule, which is unaffordable beyond
    /// d = 3; the cap trades resolution for speed in higher dimensions.
    std::size_t max_functions_per_sample = 2048;

    /// Number of points you intend to retain.  When set (and `level` is not),
    /// the level is fixed up front from this value; when both are unset the
    /// level grows adaptively with the retained count, producing a point
    /// *sequence* with good prefixes (footnote 7 of the paper).
    std::optional<std::uint64_t> expected_n{};

    /// The bound B of the linear-feedback rule (eq. 4.1): the density is
    /// 1 - eps/(2B) Phi_t(x), so smaller B means stronger feedback and better
    /// discrepancy, until |Phi_t(x)| starts exceeding B (see `on_saturation`).
    /// Unset: chosen automatically so that saturation is rare (see
    /// `Thinner::saturation_bound()` and the README).  Ignored by
    /// `Feedback::sign`.
    std::optional<std::int64_t> saturation_bound{};

    OnSaturation on_saturation = OnSaturation::clamp;

    /// Seed for the internal generator (rejection coins, random shift, and
    /// the samples drawn by `Thinner::run(n)`).
    std::uint64_t seed = 0x5EEDULL;

    Storage storage = Storage::automatic;

    /// Memory budget for the discrepancy counters.  The automatic level rule
    /// lowers the level until the dense table fits; an explicit `level` whose
    /// dense table exceeds the budget uses the sparse table instead
    /// (`Storage::automatic`).
    std::size_t dense_budget_bytes = std::size_t{512} << 20;
};

/// Thrown when `OnSaturation::fail` is selected and the field exceeds B.
class SaturationError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/// Checks `options` and returns a human-readable error message if invalid.
[[nodiscard]] inline std::expected<void, std::string> check_options(const Options& o) noexcept {
    if (o.dim == 0 || o.dim > 32) return std::unexpected("dim must be in [1, 32]");
    if (!(o.epsilon > 0.0) || !(o.epsilon <= 1.0)) return std::unexpected("epsilon must be in (0, 1]");
    if (o.level && (*o.level == 0 || *o.level > 62)) return std::unexpected("level must be in [1, 62]");
    if (o.level_offset < -62 || o.level_offset > 62) return std::unexpected("level_offset out of range");
    if (o.expected_n && *o.expected_n == 0) return std::unexpected("expected_n must be positive");
    if (o.expected_n && *o.expected_n > (std::uint64_t{1} << 31)) {
        return std::unexpected("expected_n must be at most 2^31 (discrepancies are 32-bit counters)");
    }
    if (o.saturation_bound && (*o.saturation_bound < 1 || *o.saturation_bound > (std::int64_t{1} << 30))) {
        return std::unexpected("saturation_bound must be in [1, 2^30]");
    }
    if (o.shift == Shift::fixed) {
        if (o.shift_vector.size() != o.dim) return std::unexpected("shift_vector must have dim entries");
        for (double s : o.shift_vector) {
            if (!(s >= 0.0 && s < 1.0)) return std::unexpected("shift_vector entries must lie in [0, 1)");
        }
    }
    return {};
}

/// Like `check_options` but throws `std::invalid_argument` on failure.
inline void validate_options(const Options& o) {
    if (auto r = check_options(o); !r) throw std::invalid_argument("haar::Options: " + r.error());
}

}  // namespace haar
