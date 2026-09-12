// SPDX-License-Identifier: MIT
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "check.hpp"
#include "haar/detail/table.hpp"
#include "haar/rng.hpp"

using haar::detail::DenseTable;
using haar::detail::SparseTable;

TEST(dense_table_counts) {
    DenseTable t(1024);
    CHECK_EQ(t.cell_count(), std::size_t{1024});
    for (std::uint64_t k = 0; k < 1024; k += 3) t.add(t.locate(k), k, 1);
    for (std::uint64_t k = 0; k < 1024; k += 3) t.add(t.locate(k), k, -2);
    for (std::uint64_t k = 0; k < 1024; ++k) CHECK_EQ(t.value(t.locate(k)), k % 3 == 0 ? -1 : 0);
    std::size_t nonzero = 0;
    t.for_each_nonzero([&](std::uint64_t key, std::int32_t v) {
        ++nonzero;
        CHECK_EQ(key % 3, std::uint64_t{0});
        CHECK_EQ(v, -1);
    });
    CHECK_EQ(nonzero, std::size_t{342});
}

TEST(sparse_table_matches_unordered_map) {
    SparseTable t;
    std::unordered_map<std::uint64_t, std::int32_t> ref;
    haar::Xoshiro256pp rng(7);
    // Structured keys similar to real ones: block offsets plus small indices.
    for (int round = 0; round < 200000; ++round) {
        const std::uint64_t key = (rng() % 50000) * 4096 + (rng() % 64);
        const std::int32_t delta = (rng() & 1) ? 1 : -1;
        t.reserve_extra(1);
        const auto slot = t.locate(key);
        CHECK_EQ(t.value(slot), ref.count(key) ? ref[key] : 0);
        t.add(slot, key, delta);
        ref[key] += delta;
    }
    CHECK_EQ(t.entry_count(), ref.size());
    for (const auto& [key, value] : ref) CHECK_EQ(t.value(t.locate(key)), value);
    std::size_t nonzero = 0;
    t.for_each_nonzero([&](std::uint64_t key, std::int32_t v) {
        ++nonzero;
        CHECK_EQ(ref.at(key), v);
    });
    std::size_t expected_nonzero = 0;
    for (const auto& kv : ref) expected_nonzero += kv.second != 0;
    CHECK_EQ(nonzero, expected_nonzero);
}

TEST(sparse_slots_survive_reserved_insertions) {
    // Slots recorded before a batch of insertions must remain valid as long as
    // reserve_extra() was called with the batch size.
    SparseTable t;
    haar::Xoshiro256pp rng(11);
    for (int batch = 0; batch < 300; ++batch) {
        std::vector<std::uint64_t> keys(97);
        for (auto& k : keys) k = rng() >> 8;
        t.reserve_extra(keys.size());
        std::vector<std::uint64_t> slots(keys.size());
        for (std::size_t i = 0; i < keys.size(); ++i) slots[i] = t.locate(keys[i]);
        for (std::size_t i = 0; i < keys.size(); ++i) t.add(slots[i], keys[i], 1);
        for (std::size_t i = 0; i < keys.size(); ++i) CHECK(t.value(t.locate(keys[i])) >= 1);
    }
    CHECK(t.entry_count() <= 300 * 97);
    CHECK(t.entry_count() * 10 <= t.capacity() * 7);
}
