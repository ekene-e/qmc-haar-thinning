// SPDX-License-Identifier: MIT
//
// Bit-for-bit replicability: the same options and seed must produce the same
// points, on this machine and on every other.  The golden digests below were
// produced by the reference build (GCC 13, x86-64); a change in any of them
// is a behavioural change of the algorithm and must be documented.

#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstring>
#include <span>
#include <vector>

#include "check.hpp"
#include "haar/thinner.hpp"

namespace {

std::uint64_t fnv1a(std::span<const double> values) {
    std::uint64_t h = 0xCBF29CE484222325ULL;
    for (double v : values) {
        const auto bits = std::bit_cast<std::uint64_t>(v);
        for (int i = 0; i < 8; ++i) {
            h ^= (bits >> (8 * i)) & 0xFF;
            h *= 0x100000001B3ULL;
        }
    }
    return h;
}

struct Golden {
    haar::Options options;
    std::size_t n;
    std::uint64_t digest;
};

std::vector<Golden> golden_cases() {
    std::vector<Golden> cases;
    {
        haar::Options o;
        o.dim = 2;
        o.seed = 12345;
        cases.push_back({o, 2000, 0xC1E6B5AC3E0255DBULL});
    }
    {
        haar::Options o;
        o.dim = 1;
        o.feedback = haar::Feedback::sign;
        o.shift = haar::Shift::none;
        o.epsilon = 0.25;
        o.seed = 1;
        cases.push_back({o, 3000, 0x8373065FE5D23666ULL});
    }
    {
        haar::Options o;
        o.dim = 3;
        o.scale_set = haar::ScaleSet::box;
        o.storage = haar::Storage::sparse;
        o.seed = 777;
        cases.push_back({o, 1500, 0x40439B31AAA42753ULL});
    }
    return cases;
}

}  // namespace

TEST(same_seed_same_points) {
    haar::Options o;
    o.dim = 2;
    o.seed = 42;
    const auto a = haar::thin(5000, o);
    const auto b = haar::thin(5000, o);
    CHECK(a == b);
    o.seed = 43;
    const auto c = haar::thin(5000, o);
    CHECK(a != c);
}

TEST(dense_and_sparse_storage_agree) {
    haar::Options o;
    o.dim = 2;
    o.seed = 9;
    o.storage = haar::Storage::dense;
    const auto dense = haar::thin(4000, o);
    o.storage = haar::Storage::sparse;
    const auto sparse = haar::thin(4000, o);
    CHECK(dense == sparse);
}

TEST(golden_digests) {
    for (const Golden& g : golden_cases()) {
        const auto pts = haar::thin(g.n, g.options);
        const std::uint64_t digest = fnv1a(pts);
        if (g.digest == 0) {
            std::printf("\n    digest for dim=%u n=%zu: 0x%016llXULL", g.options.dim, g.n,
                        static_cast<unsigned long long>(digest));
        } else {
            CHECK_EQ(digest, g.digest);
        }
    }
}

TEST(run_with_external_generator_equals_manual_offers) {
    // run(n, gen) is exactly "draw a point from gen, offer it" in a loop.
    haar::Options o;
    o.dim = 2;
    o.expected_n = 1000;
    o.seed = 3;

    haar::Thinner batch(o);
    haar::Xoshiro256pp g1(100);
    batch.run(1000, g1);

    haar::Thinner manual(o);
    haar::Xoshiro256pp g2(100);
    std::vector<double> x(2);
    while (manual.size() < 1000) {
        for (auto& c : x) c = haar::detail::uniform01(g2);
        manual.offer(x);
    }
    CHECK(std::ranges::equal(batch.points(), manual.points()));
    CHECK_EQ(batch.stats().offered, manual.stats().offered);
    CHECK_EQ(g1(), g2());  // both consumed the same number of draws
}
