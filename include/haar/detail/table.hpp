// SPDX-License-Identifier: MIT
#pragma once

/// @file detail/table.hpp
/// Storage for the Haar discrepancies phi_t(H) = sum_{z in A_t} H(z).
///
/// Both tables expose the same tiny interface, used by the hot loop of the
/// thinner:
///
///   reserve_extra(m)  -- guarantee that m further insertions will not rehash,
///                        so slots returned by `locate` stay valid this step;
///   prefetch(key)     -- hint the cache line that `locate(key)` will touch;
///   locate(key)       -- slot holding `key`, or the empty slot it would take;
///   value(slot)       -- current counter (0 for an empty slot);
///   add(slot, key, d) -- create the entry if needed and add d.
///
/// Discrepancies are exact integers in [-n, n], so counters are `int32_t`.

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#if defined(__GNUC__) || defined(__clang__)
#define HAAR_PREFETCH(ptr) __builtin_prefetch((ptr), 0, 3)
#elif defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
#include <xmmintrin.h>
#define HAAR_PREFETCH(ptr) _mm_prefetch(reinterpret_cast<const char*>(ptr), _MM_HINT_T0)
#else
#define HAAR_PREFETCH(ptr) ((void)0)
#endif

namespace haar::detail {

/// One counter per cell; the key *is* the slot.
class DenseTable {
public:
    using Slot = std::uint64_t;
    /// `locate(key) == key`, so callers need not record slots.
    static constexpr bool kIdentitySlots = true;

    DenseTable() = default;
    explicit DenseTable(std::uint64_t cells)
        : cells_(std::make_unique<std::int32_t[]>(static_cast<std::size_t>(cells))),
          size_(static_cast<std::size_t>(cells)) {}

    void reserve_extra(std::size_t) noexcept {}

    void prefetch(std::uint64_t key) const noexcept { HAAR_PREFETCH(cells_.get() + key); }

    [[nodiscard]] Slot locate(std::uint64_t key) const noexcept { return key; }

    [[nodiscard]] std::int32_t value(Slot slot) const noexcept { return cells_[slot]; }

    void add(Slot slot, std::uint64_t, std::int32_t delta) noexcept { cells_[slot] += delta; }

    [[nodiscard]] std::size_t cell_count() const noexcept { return size_; }

    /// Memory footprint in bytes.
    [[nodiscard]] std::size_t bytes() const noexcept { return size_ * sizeof(std::int32_t); }

    /// Calls fn(key, value) for every nonzero counter.
    template <class Fn>
    void for_each_nonzero(Fn&& fn) const {
        for (std::size_t i = 0; i < size_; ++i) {
            if (cells_[i] != 0) fn(static_cast<std::uint64_t>(i), cells_[i]);
        }
    }

private:
    std::unique_ptr<std::int32_t[]> cells_;
    std::size_t size_ = 0;
};

/// Open-addressing hash table with linear probing, 16-byte entries and a
/// power-of-two capacity kept below 70% load.
class SparseTable {
public:
    using Slot = std::uint64_t;
    static constexpr bool kIdentitySlots = false;

    SparseTable() { rehash(std::size_t{1} << 12); }

    void reserve_extra(std::size_t extra) {
        if (count_ + extra > max_load_) {
            std::size_t cap = entries_.size();
            while ((count_ + extra) * 10 > cap * 7) cap *= 2;
            rehash(cap);
        }
    }

    void prefetch(std::uint64_t key) const noexcept { HAAR_PREFETCH(entries_.data() + home(key)); }

    [[nodiscard]] Slot locate(std::uint64_t key) const noexcept {
        std::size_t i = home(key);
        while (entries_[i].key != key && entries_[i].key != kEmpty) i = (i + 1) & mask_;
        return i;
    }

    [[nodiscard]] std::int32_t value(Slot slot) const noexcept { return entries_[slot].value; }

    /// `slot` may come from an earlier `locate(key)`: entries are never removed,
    /// so the probe chain of `key` still passes through it, but another key
    /// located in the same step may have taken that empty slot meanwhile.
    /// Resume probing from `slot` to find or create the entry for `key`.
    void add(Slot slot, std::uint64_t key, std::int32_t delta) noexcept {
        std::size_t i = slot;
        while (entries_[i].key != key && entries_[i].key != kEmpty) i = (i + 1) & mask_;
        Entry& e = entries_[i];
        if (e.key == kEmpty) {
            e.key = key;
            ++count_;
        }
        e.value += delta;
    }

    [[nodiscard]] std::size_t entry_count() const noexcept { return count_; }
    [[nodiscard]] std::size_t capacity() const noexcept { return entries_.size(); }
    [[nodiscard]] std::size_t bytes() const noexcept { return entries_.size() * sizeof(Entry); }

    template <class Fn>
    void for_each_nonzero(Fn&& fn) const {
        for (const Entry& e : entries_) {
            if (e.key != kEmpty && e.value != 0) fn(e.key, e.value);
        }
    }

private:
    struct Entry {
        std::uint64_t key;
        std::int32_t value;
        std::int32_t padding_;
    };
    static constexpr std::uint64_t kEmpty = ~std::uint64_t{0};

    [[nodiscard]] std::size_t home(std::uint64_t key) const noexcept {
        // SplitMix64 finaliser: keys are structured (block offset + index), so
        // a bare multiplicative hash would cluster.
        key = (key ^ (key >> 30)) * 0xBF58476D1CE4E5B9ULL;
        key = (key ^ (key >> 27)) * 0x94D049BB133111EBULL;
        key ^= key >> 31;
        return static_cast<std::size_t>(key) & mask_;
    }

    void rehash(std::size_t capacity) {
        capacity = std::bit_ceil(capacity);
        std::vector<Entry> old = std::move(entries_);
        entries_.assign(capacity, Entry{kEmpty, 0, 0});
        mask_ = capacity - 1;
        max_load_ = capacity / 10 * 7;
        count_ = 0;
        for (const Entry& e : old) {
            if (e.key != kEmpty) add(locate(e.key), e.key, e.value);
        }
    }

    std::vector<Entry> entries_;
    std::size_t mask_ = 0;
    std::size_t max_load_ = 0;
    std::size_t count_ = 0;
};

}  // namespace haar::detail
