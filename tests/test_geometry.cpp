// SPDX-License-Identifier: MIT
#include <cstdint>
#include <set>
#include <vector>

#include "check.hpp"
#include "haar/detail/geometry.hpp"

using haar::detail::Geometry;
using haar::detail::cell_shift;

namespace {

// All j in [0,R]^d with sum <= L, any order.
std::vector<std::vector<unsigned>> all_scales(unsigned d, unsigned R, unsigned L) {
    std::vector<std::vector<unsigned>> out;
    std::vector<unsigned> j(d, 0);
    for (;;) {
        unsigned sum = 0;
        for (unsigned v : j) sum += v;
        if (sum <= L) out.push_back(j);
        std::size_t i = 0;
        while (i < d && ++j[i] > R) j[i++] = 0;
        if (i == d) break;
    }
    return out;
}

std::uint64_t block_size(const std::vector<unsigned>& j) {
    unsigned bits = 0;
    for (unsigned v : j) bits += cell_shift(v);
    return std::uint64_t{1} << bits;
}

std::uint64_t binomial(unsigned n, unsigned k) {
    std::uint64_t r = 1;
    for (unsigned i = 1; i <= k; ++i) r = r * (n - k + i) / i;
    return r;
}

}  // namespace

TEST(box_counts) {
    for (unsigned d = 1; d <= 4; ++d) {
        for (unsigned l = 1; l <= 5; ++l) {
            const Geometry g = Geometry::make(d, l, d * l);
            std::uint64_t scales = 1, cells = 1;
            for (unsigned i = 0; i < d; ++i) {
                scales *= l + 1;
                cells *= std::uint64_t{1} << l;
            }
            CHECK_EQ(g.scale_count, scales);
            CHECK_EQ(g.cell_count, cells);
            CHECK_EQ(g.function_count(), scales - 1);
        }
    }
}

TEST(hyperbolic_counts) {
    for (unsigned d = 1; d <= 4; ++d) {
        for (unsigned L = 1; L <= 8; ++L) {
            const Geometry g = Geometry::make(d, L, L);
            CHECK_EQ(g.scale_count, binomial(L + d, d));
            std::uint64_t cells = 0;
            for (const auto& j : all_scales(d, L, L)) cells += block_size(j);
            CHECK_EQ(g.cell_count, cells);
        }
    }
}

TEST(offsets_partition_the_table) {
    for (unsigned d = 1; d <= 3; ++d) {
        const Geometry g = Geometry::make(d, 4, 6);
        const auto scales = all_scales(d, 4, 6);
        CHECK_EQ(g.offsets.size(), scales.size());
        // Reconstruct the enumeration order from prefixes and check blocks tile [0, cell_count).
        std::uint64_t expected_offset = 0;
        std::size_t s = 0;
        for (std::size_t p = 0; p < g.prefix_count(); ++p) {
            std::vector<unsigned> j(d, 0);
            for (unsigned i = 0; i + 1 < d; ++i) j[i] = g.prefix_levels[p * (d - 1) + i];
            CHECK_EQ(g.prefix_positions[p], s);
            for (unsigned last = 0; last <= g.prefix_bounds[p]; ++last, ++s) {
                j[d - 1] = last;
                CHECK_EQ(g.offsets[s], expected_offset);
                expected_offset += block_size(j);
            }
        }
        CHECK_EQ(s, scales.size());
        CHECK_EQ(expected_offset, g.cell_count);
    }
}

TEST(enumeration_is_lexicographic_and_complete) {
    const unsigned d = 3;
    const Geometry g = Geometry::make(d, 3, 5);
    std::set<std::vector<unsigned>> seen;
    std::vector<unsigned> previous;
    for (std::size_t p = 0; p < g.prefix_count(); ++p) {
        for (unsigned last = 0; last <= g.prefix_bounds[p]; ++last) {
            std::vector<unsigned> j(d);
            for (unsigned i = 0; i + 1 < d; ++i) j[i] = g.prefix_levels[p * (d - 1) + i];
            j[d - 1] = last;
            if (!previous.empty()) CHECK(previous < j);
            previous = j;
            seen.insert(j);
        }
    }
    const auto expected = all_scales(d, 3, 5);
    CHECK_EQ(seen.size(), expected.size());
    for (const auto& j : expected) CHECK(seen.count(j) == 1);
}

TEST(rejects_oversized_layouts) {
    CHECK_THROWS(Geometry::make(2, 40, 80), std::length_error);
    CHECK_THROWS(Geometry::make(0, 1, 1), std::invalid_argument);
}
