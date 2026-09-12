// SPDX-License-Identifier: MIT
#pragma once

/// @file detail/geometry.hpp
/// Enumeration of Haar scale vectors and the layout of their discrepancy
/// counters.
///
/// A d-dimensional Haar function H_{j,k} is a tensor product of 1-D Haar
/// functions h_{j_i,k_i}.  For a fixed *scale vector* j, exactly one k makes
/// H_{j,k}(x) nonzero at a given point x, so a point touches one counter per
/// scale vector.  We lay the counters out scale-major: all functions with
/// scale vector j occupy a contiguous block of 2^(sum_i max(j_i - 1, 0))
/// cells, and within the block the cell index is the mixed-radix number
/// formed by the k_i.  Coarse scales therefore live in tiny, cache-resident
/// blocks at the front of the table.
///
/// Scale vectors are enumerated lexicographically (coordinate 0 slowest).  We
/// split each vector into a *prefix* (its first d - 1 levels) and the last
/// level, which the hot loop iterates over contiguously.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace haar::detail {

/// Number of bits contributed to the cell index by a 1-D level j:
/// level 0 is the constant function (1 cell), level j >= 1 has 2^(j-1) cells.
[[nodiscard]] constexpr unsigned cell_shift(unsigned level) noexcept {
    return level == 0 ? 0u : level - 1u;
}

struct Geometry {
    unsigned dim = 0;
    /// Largest level allowed in any single coordinate (R in the code comments).
    unsigned per_coordinate_level = 0;
    /// Largest allowed sum of levels (L).  Equal to dim * R for the box set.
    unsigned total_level = 0;

    /// Number of scale vectors, including the zero vector (N_l in the paper).
    std::uint64_t scale_count = 0;
    /// Total number of counters (the dense table size).
    std::uint64_t cell_count = 0;

    /// Block offset of each scale vector, in enumeration order.
    std::vector<std::uint64_t> offsets;
    /// Concatenated prefixes, (dim - 1) bytes each.
    std::vector<std::uint8_t> prefix_levels;
    /// Largest last-coordinate level admissible for each prefix.
    std::vector<std::uint8_t> prefix_bounds;
    /// Index into `offsets` of the first scale vector of each prefix.
    std::vector<std::uint32_t> prefix_positions;

    [[nodiscard]] std::size_t prefix_count() const noexcept { return prefix_bounds.size(); }

    /// Number of Haar functions that are nonzero at any given point
    /// (all scale vectors except the zero vector).
    [[nodiscard]] std::uint64_t function_count() const noexcept { return scale_count - 1; }

    /// Builds the geometry for all j in [0, R]^dim with sum_i j_i <= L.
    /// Throws `std::length_error` when the layout does not fit in 62 bits.
    [[nodiscard]] static Geometry make(unsigned dim, unsigned per_coordinate_level, unsigned total_level) {
        if (dim == 0) throw std::invalid_argument("Geometry: dim must be positive");
        Geometry g;
        g.dim = dim;
        g.per_coordinate_level = per_coordinate_level;
        g.total_level = total_level;

        constexpr std::uint64_t kMaxCells = std::uint64_t{1} << 62;
        constexpr std::uint64_t kMaxScales = std::uint64_t{1} << 30;
        const unsigned R = per_coordinate_level;
        const unsigned L = total_level;

        std::vector<unsigned> prefix(dim - 1, 0u);

        auto emit_prefix = [&] {
            unsigned sum = 0;
            unsigned shift = 0;
            for (unsigned lv : prefix) {
                sum += lv;
                shift += cell_shift(lv);
            }
            const unsigned bound = std::min(R, L - sum);
            g.prefix_positions.push_back(static_cast<std::uint32_t>(g.offsets.size()));
            g.prefix_bounds.push_back(static_cast<std::uint8_t>(bound));
            for (unsigned lv : prefix) g.prefix_levels.push_back(static_cast<std::uint8_t>(lv));
            for (unsigned j = 0; j <= bound; ++j) {
                g.offsets.push_back(g.cell_count);
                const unsigned bits = shift + cell_shift(j);
                if (bits >= 62 || g.cell_count > kMaxCells - (std::uint64_t{1} << bits)) {
                    throw std::length_error("Geometry: level too large, table exceeds 2^62 cells");
                }
                g.cell_count += std::uint64_t{1} << bits;
            }
            if (g.offsets.size() > kMaxScales) {
                throw std::length_error("Geometry: too many scale vectors (> 2^30)");
            }
        };

        if (dim == 1) {
            emit_prefix();
        } else {
            // Odometer over the prefix, honouring both constraints.
            bool more = true;
            while (more) {
                emit_prefix();
                std::size_t i = dim - 2;
                for (;;) {
                    ++prefix[i];
                    unsigned sum = 0;
                    for (unsigned lv : prefix) sum += lv;
                    if (prefix[i] <= R && sum <= L) break;
                    prefix[i] = 0;
                    if (i == 0) {
                        more = false;
                        break;
                    }
                    --i;
                }
            }
        }
        g.scale_count = g.offsets.size();
        return g;
    }
};

}  // namespace haar::detail
