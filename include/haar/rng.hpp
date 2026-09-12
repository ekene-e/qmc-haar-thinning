// SPDX-License-Identifier: MIT
#pragma once

/// @file rng.hpp
/// A small, fast, fully deterministic pseudo-random generator (xoshiro256++)
/// plus a portable uniform-[0,1) conversion.
///
/// The thinner uses this generator for its own randomness (rejection coins
/// and the optional random shift).  It is deliberately not `std::mt19937`:
/// xoshiro256++ is several times faster and, unlike `std::generate_canonical`,
/// our double conversion is specified bit-for-bit, so two runs with the same
/// seed produce identical output on every platform and compiler.

#include <bit>
#include <concepts>
#include <cstdint>
#include <limits>
#include <random>
#include <span>

namespace haar {

/// xoshiro256++ by Blackman and Vigna (public domain).  Period 2^256 - 1.
/// Satisfies `std::uniform_random_bit_generator`.
class Xoshiro256pp {
public:
    using result_type = std::uint64_t;

    /// Seeds all 256 bits of state from a 64-bit seed via SplitMix64.
    explicit constexpr Xoshiro256pp(std::uint64_t seed = 0) noexcept { reseed(seed); }

    constexpr void reseed(std::uint64_t seed) noexcept {
        std::uint64_t z = seed;
        for (auto& word : state_) word = splitmix64(z);
    }

    [[nodiscard]] static constexpr result_type min() noexcept { return 0; }
    [[nodiscard]] static constexpr result_type max() noexcept {
        return std::numeric_limits<result_type>::max();
    }

    constexpr result_type operator()() noexcept {
        const std::uint64_t result = std::rotl(state_[0] + state_[3], 23) + state_[0];
        const std::uint64_t t = state_[1] << 17;
        state_[2] ^= state_[0];
        state_[3] ^= state_[1];
        state_[1] ^= state_[2];
        state_[0] ^= state_[3];
        state_[2] ^= t;
        state_[3] = std::rotl(state_[3], 45);
        return result;
    }

    /// Uniform double in [0, 1) with 53 random bits.  Exact: (x >> 11) * 2^-53.
    constexpr double uniform() noexcept {
        return static_cast<double>((*this)() >> 11) * 0x1.0p-53;
    }

    /// Uniform 32-bit integer.
    constexpr std::uint32_t uniform32() noexcept {
        return static_cast<std::uint32_t>((*this)() >> 32);
    }

    /// Fills `out` with i.i.d. uniform coordinates in [0, 1).
    constexpr void uniform_point(std::span<double> out) noexcept {
        for (auto& c : out) c = uniform();
    }

    constexpr bool operator==(const Xoshiro256pp&) const noexcept = default;

private:
    static constexpr std::uint64_t splitmix64(std::uint64_t& z) noexcept {
        z += 0x9E3779B97F4A7C15ULL;
        std::uint64_t x = z;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        return x ^ (x >> 31);
    }

    std::uint64_t state_[4]{};
};

namespace detail {

/// Portable uniform [0,1) from any uniform random bit generator.
/// Uses 53 bits, drawn as one 64-bit word or two 32-bit words, so that the
/// result does not depend on `std::generate_canonical`'s implementation.
template <std::uniform_random_bit_generator G>
double uniform01(G& gen) {
    using R = typename G::result_type;
    if constexpr (G::max() == std::numeric_limits<std::uint64_t>::max() && G::min() == 0) {
        return static_cast<double>(static_cast<std::uint64_t>(gen()) >> 11) * 0x1.0p-53;
    } else if constexpr (G::max() == std::numeric_limits<std::uint32_t>::max() && G::min() == 0) {
        const std::uint64_t hi = static_cast<std::uint64_t>(gen());
        const std::uint64_t lo = static_cast<std::uint64_t>(gen());
        return static_cast<double>(((hi << 32) | lo) >> 11) * 0x1.0p-53;
    } else {
        static_assert(sizeof(R) > 0, "unreachable");
        return std::generate_canonical<double, 53>(gen);
    }
}

}  // namespace detail
}  // namespace haar
