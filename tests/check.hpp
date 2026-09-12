// SPDX-License-Identifier: MIT
#pragma once

// A minimal, dependency-free test harness.  Each test file includes this
// header, defines tests with TEST(name) { ... } and gets a main() for free.

#include <cmath>
#include <cstdio>
#include <exception>
#include <format>
#include <string>
#include <utility>
#include <vector>

namespace check {

using TestFn = void (*)();

inline std::vector<std::pair<std::string, TestFn>>& registry() {
    static std::vector<std::pair<std::string, TestFn>> tests;
    return tests;
}

struct Registrar {
    Registrar(const char* name, TestFn fn) { registry().emplace_back(name, fn); }
};

inline int failures = 0;

inline void fail(const char* file, int line, const std::string& msg) {
    ++failures;
    std::fprintf(stderr, "    FAIL %s:%d: %s\n", file, line, msg.c_str());
}

}  // namespace check

#define TEST(name)                                              \
    static void name();                                         \
    static ::check::Registrar name##_registrar(#name, &name);   \
    static void name()

#define CHECK(cond)                                                     \
    do {                                                                \
        if (!(cond)) ::check::fail(__FILE__, __LINE__, "CHECK(" #cond ")"); \
    } while (0)

#define CHECK_EQ(a, b)                                                                              \
    do {                                                                                            \
        const auto va_ = (a);                                                                       \
        const auto vb_ = (b);                                                                       \
        if (!(va_ == vb_))                                                                          \
            ::check::fail(__FILE__, __LINE__, std::format("{} == {}: got {} vs {}", #a, #b, va_, vb_)); \
    } while (0)

#define CHECK_CLOSE(a, b, tol)                                                                       \
    do {                                                                                             \
        const double va_ = static_cast<double>(a);                                                   \
        const double vb_ = static_cast<double>(b);                                                   \
        if (!(std::fabs(va_ - vb_) <= (tol)))                                                        \
            ::check::fail(__FILE__, __LINE__,                                                        \
                          std::format("{} ~ {} within {}: got {} vs {}", #a, #b, tol, va_, vb_));   \
    } while (0)

#define CHECK_THROWS(expr, Exception)                                                        \
    do {                                                                                     \
        bool thrown_ = false;                                                                \
        try {                                                                                \
            (void)(expr);                                                                    \
        } catch (const Exception&) {                                                         \
            thrown_ = true;                                                                  \
        } catch (...) {                                                                      \
        }                                                                                    \
        if (!thrown_) ::check::fail(__FILE__, __LINE__, "expected " #Exception " from " #expr); \
    } while (0)

int main() {
    int failed_tests = 0;
    for (const auto& [name, fn] : check::registry()) {
        const int before = check::failures;
        std::printf("  %-48s", name.c_str());
        std::fflush(stdout);
        try {
            fn();
        } catch (const std::exception& e) {
            check::fail("?", 0, std::string("unexpected exception: ") + e.what());
        }
        const bool ok = check::failures == before;
        failed_tests += !ok;
        std::printf("%s\n", ok ? "ok" : "FAILED");
    }
    std::printf("%zu tests, %d failed\n", check::registry().size(), failed_tests);
    return failed_tests == 0 ? 0 : 1;
}
