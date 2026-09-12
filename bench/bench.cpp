// SPDX-License-Identifier: MIT
//
// haar-bench: reproduces the quality and speed claims of the README.
//
//   haar-bench                 all sections, full size (a few minutes)
//   haar-bench --quick         smaller sizes and fewer repetitions
//   haar-bench --section disc|integ|speed|bound
//
// Output is GitHub-flavoured Markdown.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <functional>
#include <numbers>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <haar/haar.hpp>

namespace {

using Clock = std::chrono::steady_clock;

struct Method {
    std::string name;
    std::function<std::vector<double>(unsigned dim, std::size_t n, std::uint64_t seed)> generate;
};

std::vector<double> monte_carlo(unsigned dim, std::size_t n, std::uint64_t seed) {
    haar::Xoshiro256pp rng(seed);
    std::vector<double> pts(n * dim);
    rng.uniform_point(pts);
    return pts;
}

Method thinning(std::string name, haar::Feedback fb, haar::Shift shift, std::optional<std::int64_t> bound = {}) {
    return {std::move(name), [=](unsigned dim, std::size_t n, std::uint64_t seed) {
                haar::Options o;
                o.dim = dim;
                o.feedback = fb;
                o.shift = shift;
                o.saturation_bound = bound;
                o.seed = seed;
                return haar::thin(n, o);
            }};
}

std::vector<Method> methods() {
    return {
        {"Monte Carlo", monte_carlo},
        thinning("sign (DFG+19)", haar::Feedback::sign, haar::Shift::none),
        thinning("sign + shift (Thm 1.2)", haar::Feedback::sign, haar::Shift::random),
        thinning("linear (Thm 1.3)", haar::Feedback::linear, haar::Shift::none),
        thinning("linear + shift (default)", haar::Feedback::linear, haar::Shift::random),
        thinning("greedy B=1 + shift", haar::Feedback::linear, haar::Shift::random, 1),
    };
}

// ---------------------------------------------------------------- discrepancy

void section_discrepancy(bool quick) {
    const int reps = quick ? 2 : 4;
    const std::vector<std::size_t> sizes = quick ? std::vector<std::size_t>{1 << 10, 1 << 12, 1 << 14}
                                                 : std::vector<std::size_t>{1 << 10, 1 << 12, 1 << 14, 1 << 16};
    for (unsigned d : {1u, 2u}) {
        std::printf("\n### Star discrepancy, d = %u (unnormalised, mean of %d runs)\n\n", d, reps);
        std::printf("| n |");
        for (const auto& m : methods()) std::printf(" %s |", m.name.c_str());
        std::printf("\n|---|");
        for (std::size_t i = 0; i < methods().size(); ++i) std::printf("---:|");
        std::printf("\n");
        for (std::size_t n : sizes) {
            std::printf("| %zu |", n);
            for (const auto& m : methods()) {
                double sum = 0.0;
                for (int r = 0; r < reps; ++r) {
                    const auto pts = m.generate(d, n, 1000 + static_cast<std::uint64_t>(r));
                    sum += *haar::star_discrepancy({pts, d});
                }
                std::printf(" %.2f |", sum / reps);
            }
            std::printf("\n");
            std::fflush(stdout);
        }
    }

    std::printf("\n### L2-star discrepancy, d = 3 (unnormalised, mean of %d runs)\n\n", reps);
    std::printf("| n |");
    for (const auto& m : methods()) std::printf(" %s |", m.name.c_str());
    std::printf("\n|---|");
    for (std::size_t i = 0; i < methods().size(); ++i) std::printf("---:|");
    std::printf("\n");
    for (std::size_t n : sizes) {
        if (n > 20000) break;
        std::printf("| %zu |", n);
        for (const auto& m : methods()) {
            double sum = 0.0;
            for (int r = 0; r < reps; ++r) {
                const auto pts = m.generate(3, n, 2000 + static_cast<std::uint64_t>(r));
                sum += haar::l2_star_discrepancy({pts, 3});
            }
            std::printf(" %.3f |", sum / reps);
        }
        std::printf("\n");
        std::fflush(stdout);
    }
}

// ---------------------------------------------------------------- integration

struct TestFunction {
    std::string name;
    std::function<double(std::span<const double>)> f;
    std::function<double(unsigned)> exact;  // integral over [0,1)^d
};

std::vector<TestFunction> test_functions() {
    constexpr double pi = std::numbers::pi;
    return {
        {"sum_i sin(2 pi 16 x_i)",
         [](std::span<const double> x) {
             double s = 0.0;
             for (double c : x) s += std::sin(2 * pi * 16 * c);
             return s;
         },
         [](unsigned) { return 0.0; }},
        {"prod_i exp(-16 (x_i - 1/2)^2)",
         [](std::span<const double> x) {
             double p = 1.0;
             for (double c : x) p *= std::exp(-16.0 * (c - 0.5) * (c - 0.5));
             return p;
         },
         [](unsigned d) { return std::pow(std::sqrt(pi) / 4.0 * std::erf(2.0), d); }},
        {"prod_i sqrt(x_i)",
         [](std::span<const double> x) {
             double p = 1.0;
             for (double c : x) p *= std::sqrt(c);
             return p;
         },
         [](unsigned d) { return std::pow(2.0 / 3.0, d); }},
        {"1[ sum_i x_i < d/2 ]",
         [](std::span<const double> x) {
             double s = 0.0;
             for (double c : x) s += c;
             return s < static_cast<double>(x.size()) / 2.0 ? 1.0 : 0.0;
         },
         [](unsigned) { return 0.5; }},
    };
}

void section_integration(bool quick) {
    const int reps = quick ? 8 : 24;
    const std::vector<std::size_t> sizes = quick ? std::vector<std::size_t>{1 << 10, 1 << 13, 1 << 16}
                                                 : std::vector<std::size_t>{1 << 10, 1 << 12, 1 << 14, 1 << 16, 1 << 18};
    const auto ms = methods();
    for (unsigned d : {1u, 2u, 3u}) {
        for (const auto& tf : test_functions()) {
            std::printf("\n### RMS integration error, d = %u, f = %s (%d runs)\n\n", d, tf.name.c_str(), reps);
            std::printf("| n |");
            for (const auto& m : ms) std::printf(" %s |", m.name.c_str());
            std::printf("\n|---|");
            for (std::size_t i = 0; i < ms.size(); ++i) std::printf("---:|");
            std::printf("\n");
            const double exact = tf.exact(d);
            for (std::size_t n : sizes) {
                std::printf("| %zu |", n);
                for (const auto& m : ms) {
                    double sq = 0.0;
                    for (int r = 0; r < reps; ++r) {
                        const auto pts = m.generate(d, n, 5000 + static_cast<std::uint64_t>(r));
                        const double err = haar::mean({pts, d}, tf.f) - exact;
                        sq += err * err;
                    }
                    std::printf(" %.2e |", std::sqrt(sq / reps));
                }
                std::printf("\n");
                std::fflush(stdout);
            }
        }
    }
}

// ---------------------------------------------------------------- speed

void section_speed(bool quick) {
    const std::size_t n = quick ? std::size_t{200000} : std::size_t{1000000};
    std::printf("\n### Throughput, n = %zu points (linear feedback, shift on)\n\n", n);
    std::printf("| d | scale set | level | Haar functions / sample | storage | table MiB | time (s) | points / s |\n");
    std::printf("|---|---|---:|---:|---|---:|---:|---:|\n");
    for (unsigned d : {1u, 2u, 3u, 4u}) {
        for (haar::ScaleSet set : {haar::ScaleSet::hyperbolic, haar::ScaleSet::box}) {
            haar::Options o;
            o.dim = d;
            o.scale_set = set;
            o.expected_n = n;
            o.seed = 1;
            haar::Thinner t(o);
            const auto start = Clock::now();
            t.run(n);
            const double secs = std::chrono::duration<double>(Clock::now() - start).count();
            std::printf("| %u | %s | %u | %llu | %s | %.1f | %.2f | %.0f |\n", d, set == haar::ScaleSet::box ? "box" : "hyperbolic",
                        t.level(), static_cast<unsigned long long>(t.function_count()), t.dense() ? "dense" : "sparse",
                        static_cast<double>(t.table_bytes()) / (1u << 20), secs, static_cast<double>(n) / secs);
            std::fflush(stdout);
        }
    }
}

// ---------------------------------------------------------------- saturation bound calibration

void section_bound(bool quick) {
    const std::size_t n = quick ? std::size_t{65536} : std::size_t{1 << 20};
    std::printf("\n### Linear feedback: automatic bound B versus the field |Phi_t(x)|, n = %zu\n\n", n);
    std::printf("Saturation events are decisions where |Phi_t(x)| > B and the density was clamped.\n\n");
    std::printf("| d | eps | N | level | B (auto) | max abs Phi | saturations | max abs phi(H) | rejection rate |\n");
    std::printf("|---|---:|---:|---:|---:|---:|---:|---:|---:|\n");
    for (unsigned d : {1u, 2u, 3u}) {
        for (double eps : {0.1, 0.5, 1.0}) {
            haar::Options o;
            o.dim = d;
            o.epsilon = eps;
            o.expected_n = n;
            o.seed = 3;
            haar::Thinner t(o);
            t.run(n);
            const haar::Stats& s = t.stats();
            std::printf("| %u | %.1f | %llu | %u | %lld | %lld | %llu | %lld | %.4f |\n", d, eps,
                        static_cast<unsigned long long>(t.scale_count()), t.level(), static_cast<long long>(t.saturation_bound()),
                        static_cast<long long>(s.max_abs_field), static_cast<unsigned long long>(s.saturated),
                        static_cast<long long>(t.max_abs_haar_discrepancy()),
                        static_cast<double>(s.rejected) / static_cast<double>(s.offered));
            std::fflush(stdout);
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    bool quick = false;
    std::string section = "all";
    for (int i = 1; i < argc; ++i) {
        const std::string_view a = argv[i];
        if (a == "--quick") quick = true;
        else if (a == "--section" && i + 1 < argc) section = argv[++i];
        else {
            std::fprintf(stderr, "usage: haar-bench [--quick] [--section disc|integ|speed|bound]\n");
            return 2;
        }
    }
    std::printf("# haar-thinning benchmarks (%s)\n", quick ? "quick" : "full");
    if (section == "all" || section == "disc") section_discrepancy(quick);
    if (section == "all" || section == "integ") section_integration(quick);
    if (section == "all" || section == "speed") section_speed(quick);
    if (section == "all" || section == "bound") section_bound(quick);
    return 0;
}
