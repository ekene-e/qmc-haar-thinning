// SPDX-License-Identifier: MIT
//
// haar-thin: command-line front end.
//
//   haar-thin -n 100000 -d 2 > points.csv
//   haar-thin -n 1000 -d 3 --feedback sign --seed 7 --format tsv -o pts.tsv --stats
//   cat samples.txt | haar-thin --stream -d 2 -n 5000 > kept.txt
//
// Run `haar-thin --help` for all options.

#include <charconv>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <haar/haar.hpp>

namespace {

constexpr std::string_view kUsage = R"(haar-thin: online Haar-thinning of uniform samples into low-discrepancy points

Usage: haar-thin [options]

Core
  -n, --points N           number of points to retain (required unless --stream)
  -d, --dim D              dimension (default 1)
  -e, --epsilon E          over-sampling: at most (1+E)n samples used (default 0.5)
  --feedback sign|linear   voting rule (default linear); Theorem 1.2 uses sign, 1.3 linear
  --shift on|off|s1,..,sd  uniform random shift of the Haar system (default on)
  --seed S                 seed for coins, shift and generated samples (default 24301)

Resolution
  --scale-set hyperbolic|box   set of Haar scale vectors (default hyperbolic)
  --level L                fixed resolution (default: automatic from n)
  --level-offset K         adjust the automatic rule (default 1)
  --adaptive               grow the resolution with the retained count (sequence mode)
  --bound B                linear feedback saturation bound B (default automatic)
  --fail-on-saturation     abort instead of clamping when |Phi| > B

Storage
  --storage auto|dense|sparse   counter storage (default auto)
  --dense-budget-mb M      largest dense table for auto (default 256)

Input / output
  --stream                 read samples from stdin (one point per line) instead of generating them
  -i, --input FILE         like --stream but from FILE
  -o, --output FILE        write points to FILE (default stdout)
  --format csv|tsv|space|f64   text separators or raw little-endian doubles (default csv)
  --precision P            significant digits for text output (default 17, exact)
  --stats                  print run statistics as JSON to stderr
  --report                 print discrepancy measures to stderr (star: d<=2; L2-star: n<=20000)
  -q, --quiet              no informational output
  -h, --help               this message
  --version                print version
)";

[[noreturn]] void die(const std::string& msg) {
    std::fprintf(stderr, "haar-thin: %s\n", msg.c_str());
    std::exit(2);
}

template <class T>
T parse_number(std::string_view text, std::string_view flag) {
    T value{};
    const char* begin = text.data();
    const char* end = begin + text.size();
    const auto result = std::from_chars(begin, end, value);
    if (result.ec != std::errc{} || result.ptr != end) die(std::string("invalid value for ") + std::string(flag) + ": " + std::string(text));
    return value;
}

struct Cli {
    haar::Options options;
    std::optional<std::size_t> n;
    bool adaptive = false;
    bool stream = false;
    std::string input;
    std::string output;
    std::string format = "csv";
    int precision = 17;
    bool stats = false;
    bool report = false;
    bool quiet = false;
};

Cli parse(int argc, char** argv) {
    Cli cli;
    cli.options.seed = 24301;
    auto need = [&](int& i, std::string_view flag) -> std::string_view {
        if (i + 1 >= argc) die(std::string("missing value for ") + std::string(flag));
        return argv[++i];
    };
    for (int i = 1; i < argc; ++i) {
        const std::string_view a = argv[i];
        if (a == "-h" || a == "--help") {
            std::fwrite(kUsage.data(), 1, kUsage.size(), stdout);
            std::exit(0);
        } else if (a == "--version") {
            std::printf("haar-thin %s\n", HAAR_VERSION_STRING);
            std::exit(0);
        } else if (a == "-n" || a == "--points") {
            cli.n = parse_number<std::size_t>(need(i, a), a);
        } else if (a == "-d" || a == "--dim") {
            cli.options.dim = parse_number<unsigned>(need(i, a), a);
        } else if (a == "-e" || a == "--epsilon") {
            cli.options.epsilon = parse_number<double>(need(i, a), a);
        } else if (a == "--feedback") {
            const auto v = need(i, a);
            if (v == "sign") cli.options.feedback = haar::Feedback::sign;
            else if (v == "linear") cli.options.feedback = haar::Feedback::linear;
            else die("--feedback must be sign or linear");
        } else if (a == "--shift") {
            const auto v = need(i, a);
            if (v == "on") cli.options.shift = haar::Shift::random;
            else if (v == "off") cli.options.shift = haar::Shift::none;
            else {
                cli.options.shift = haar::Shift::fixed;
                std::string item;
                std::stringstream ss{std::string(v)};
                while (std::getline(ss, item, ',')) cli.options.shift_vector.push_back(parse_number<double>(item, a));
            }
        } else if (a == "--seed") {
            cli.options.seed = parse_number<std::uint64_t>(need(i, a), a);
        } else if (a == "--scale-set") {
            const auto v = need(i, a);
            if (v == "hyperbolic") cli.options.scale_set = haar::ScaleSet::hyperbolic;
            else if (v == "box") cli.options.scale_set = haar::ScaleSet::box;
            else die("--scale-set must be hyperbolic or box");
        } else if (a == "--level") {
            cli.options.level = parse_number<unsigned>(need(i, a), a);
        } else if (a == "--level-offset") {
            cli.options.level_offset = parse_number<int>(need(i, a), a);
        } else if (a == "--adaptive") {
            cli.adaptive = true;
        } else if (a == "--bound") {
            cli.options.saturation_bound = parse_number<std::int64_t>(need(i, a), a);
        } else if (a == "--fail-on-saturation") {
            cli.options.on_saturation = haar::OnSaturation::fail;
        } else if (a == "--storage") {
            const auto v = need(i, a);
            if (v == "auto") cli.options.storage = haar::Storage::automatic;
            else if (v == "dense") cli.options.storage = haar::Storage::dense;
            else if (v == "sparse") cli.options.storage = haar::Storage::sparse;
            else die("--storage must be auto, dense or sparse");
        } else if (a == "--dense-budget-mb") {
            cli.options.dense_budget_bytes = parse_number<std::size_t>(need(i, a), a) << 20;
        } else if (a == "--stream") {
            cli.stream = true;
        } else if (a == "-i" || a == "--input") {
            cli.stream = true;
            cli.input = std::string(need(i, a));
        } else if (a == "-o" || a == "--output") {
            cli.output = std::string(need(i, a));
        } else if (a == "--format") {
            cli.format = std::string(need(i, a));
            if (cli.format != "csv" && cli.format != "tsv" && cli.format != "space" && cli.format != "f64") {
                die("--format must be csv, tsv, space or f64");
            }
        } else if (a == "--precision") {
            cli.precision = parse_number<int>(need(i, a), a);
            if (cli.precision < 1 || cli.precision > 17) die("--precision must be in [1, 17]");
        } else if (a == "--stats") {
            cli.stats = true;
        } else if (a == "--report") {
            cli.report = true;
        } else if (a == "-q" || a == "--quiet") {
            cli.quiet = true;
        } else {
            die("unknown option " + std::string(a) + " (see --help)");
        }
    }
    if (!cli.stream && !cli.n) die("-n is required (or use --stream)");
    if (cli.n && !cli.adaptive && !cli.options.level) cli.options.expected_n = *cli.n;
    if (auto r = haar::check_options(cli.options); !r) die(r.error());
    return cli;
}

/// Buffered text/binary writer for points.
class Writer {
public:
    Writer(std::FILE* out, const std::string& format, int precision, unsigned dim)
        : out_(out), dim_(dim), precision_(precision), binary_(format == "f64"),
          separator_(format == "csv" ? ',' : format == "tsv" ? '\t' : ' ') {
        buffer_.reserve(1 << 16);
    }

    void write(std::span<const double> point) {
        if (binary_) {
            const char* bytes = reinterpret_cast<const char*>(point.data());
            buffer_.insert(buffer_.end(), bytes, bytes + point.size_bytes());
        } else {
            char tmp[64];
            for (unsigned i = 0; i < dim_; ++i) {
                const auto r = precision_ == 17 ? std::to_chars(tmp, tmp + sizeof tmp, point[i])
                                                : std::to_chars(tmp, tmp + sizeof tmp, point[i], std::chars_format::general, precision_);
                buffer_.insert(buffer_.end(), tmp, r.ptr);
                buffer_.push_back(i + 1 < dim_ ? separator_ : '\n');
            }
        }
        if (buffer_.size() > (1u << 16) - 1024) flush();
    }

    void flush() {
        if (!buffer_.empty()) std::fwrite(buffer_.data(), 1, buffer_.size(), out_);
        buffer_.clear();
    }

    ~Writer() { flush(); }

private:
    std::FILE* out_;
    unsigned dim_;
    int precision_;
    bool binary_;
    char separator_;
    std::vector<char> buffer_;
};

void print_stats(const haar::Thinner& t) {
    const haar::Stats& s = t.stats();
    std::fprintf(stderr,
                 "{\"config\": \"%s\", \"offered\": %llu, \"retained\": %llu, \"rejected\": %llu, "
                 "\"rejection_rate\": %.6f, \"saturated\": %llu, \"max_abs_field\": %lld, "
                 "\"saturation_bound\": %lld, \"rebuilds\": %llu, \"level\": %u, \"scale_vectors\": %llu, "
                 "\"table_bytes\": %zu, \"max_abs_haar_discrepancy\": %lld}\n",
                 t.describe().c_str(), static_cast<unsigned long long>(s.offered),
                 static_cast<unsigned long long>(s.retained), static_cast<unsigned long long>(s.rejected),
                 s.offered ? static_cast<double>(s.rejected) / static_cast<double>(s.offered) : 0.0,
                 static_cast<unsigned long long>(s.saturated), static_cast<long long>(s.max_abs_field),
                 static_cast<long long>(t.saturation_bound()), static_cast<unsigned long long>(s.rebuilds),
                 t.level(), static_cast<unsigned long long>(t.scale_count()), t.table_bytes(),
                 static_cast<long long>(t.max_abs_haar_discrepancy()));
}

void print_report(const haar::Thinner& t) {
    const haar::PointsView view{t.points(), t.dim()};
    std::fprintf(stderr, "report: n=%zu d=%u\n", view.size(), t.dim());
    if (auto star = haar::star_discrepancy(view)) {
        std::fprintf(stderr, "  star discrepancy (unnormalised): %.4f   (normalised: %.3e)\n", *star,
                     *star / static_cast<double>(view.size()));
    } else {
        std::fprintf(stderr, "  star discrepancy: exact computation available only for d <= 2\n");
    }
    if (view.size() <= 20000) {
        const double l2 = haar::l2_star_discrepancy(view);
        std::fprintf(stderr, "  L2-star discrepancy (unnormalised): %.4f   (normalised: %.3e)\n", l2,
                     l2 / static_cast<double>(view.size()));
    } else {
        std::fprintf(stderr, "  L2-star discrepancy: skipped for n > 20000 (O(n^2))\n");
    }
}

int run_generate(const Cli& cli, std::FILE* out) {
    haar::Thinner thinner(cli.options);
    if (!cli.quiet) std::fprintf(stderr, "haar-thin: %s\n", thinner.describe().c_str());
    thinner.run(*cli.n);
    Writer writer(out, cli.format, cli.precision, cli.options.dim);
    for (std::size_t i = 0; i < thinner.size(); ++i) writer.write(thinner.point(i));
    writer.flush();
    if (cli.stats) print_stats(thinner);
    if (cli.report) print_report(thinner);
    return 0;
}

int run_stream(const Cli& cli, std::FILE* out) {
    std::ifstream file;
    std::istream* in = &std::cin;
    if (!cli.input.empty()) {
        file.open(cli.input);
        if (!file) die("cannot open " + cli.input);
        in = &file;
    }
    haar::Thinner thinner(cli.options);
    if (!cli.quiet) std::fprintf(stderr, "haar-thin: %s\n", thinner.describe().c_str());
    Writer writer(out, cli.format, cli.precision, cli.options.dim);
    std::vector<double> x(cli.options.dim);
    std::string line;
    std::size_t line_number = 0;
    while ((!cli.n || thinner.size() < *cli.n) && std::getline(*in, line)) {
        ++line_number;
        std::string_view rest = line;
        unsigned parsed = 0;
        while (parsed < cli.options.dim) {
            while (!rest.empty() && (rest.front() == ' ' || rest.front() == '\t' || rest.front() == ',' || rest.front() == '\r')) rest.remove_prefix(1);
            if (rest.empty()) break;
            double v = 0.0;
            const auto r = std::from_chars(rest.data(), rest.data() + rest.size(), v);
            if (r.ec != std::errc{}) break;
            x[parsed++] = v;
            rest.remove_prefix(static_cast<std::size_t>(r.ptr - rest.data()));
        }
        if (parsed == 0 && line.find_first_not_of(" \t\r") == std::string::npos) continue;  // blank line
        if (parsed != cli.options.dim) die("line " + std::to_string(line_number) + ": expected " + std::to_string(cli.options.dim) + " coordinates");
        if (thinner.offer(x)) writer.write(x);
    }
    writer.flush();
    if (cli.n && thinner.size() < *cli.n && !cli.quiet) {
        std::fprintf(stderr, "haar-thin: input exhausted after %zu points (requested %zu)\n", thinner.size(), *cli.n);
    }
    if (cli.stats) print_stats(thinner);
    if (cli.report) print_report(thinner);
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Cli cli = parse(argc, argv);
        std::FILE* out = stdout;
        if (!cli.output.empty()) {
            out = std::fopen(cli.output.c_str(), cli.format == "f64" ? "wb" : "w");
            if (!out) die("cannot open " + cli.output + " for writing");
        }
        const int rc = cli.stream ? run_stream(cli, out) : run_generate(cli, out);
        if (out != stdout) std::fclose(out);
        return rc;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "haar-thin: error: %s\n", e.what());
        return 1;
    }
}
