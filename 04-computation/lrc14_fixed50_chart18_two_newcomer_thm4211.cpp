// Primary exact audit for THM-4211's fixed-q=50 eighteen-label chart.
// The finite branch uses literal joint walls and an upward subset-zeta sum;
// the cofinal branch uses the exact component discrepancy inequality.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

using i64 = std::int64_t;
using u32 = std::uint32_t;
using i128 = __int128_t;

namespace {

constexpr std::array<int, 18> CHART = {
    8,16,40,42,80,84,85,88,95,120,126,143,145,168,193,240,252,286};
constexpr std::array<int, 30> POOL = {
    8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
    120,126,132,143,145,168,170,176,190,193,240,252,264,286,290};
constexpr int Q = 50;
constexpr std::array<int, 9> HOSTILE_BODY_VALUES = {
    16,85,88,95,143,145,168,193,240};
constexpr std::uint64_t FNV_OFFSET = UINT64_C(0xcbf29ce484222325);
constexpr std::uint64_t FNV_PRIME = UINT64_C(0x100000001b3);

class Fnv1a64 {
  public:
    void add_u64_le(std::uint64_t word) {
        for (int shift = 0; shift < 64; shift += 8) {
            value_ ^= static_cast<std::uint8_t>((word >> shift) & UINT64_C(0xff));
            value_ *= FNV_PRIME;
        }
    }
    std::uint64_t value() const { return value_; }

  private:
    std::uint64_t value_ = FNV_OFFSET;
};

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

i64 exact_lcm(i64 left, i64 right) {
    const i64 common = std::gcd(left, right);
    const i128 value = static_cast<i128>(left / common) * right;
    require(value <= std::numeric_limits<i64>::max(), "lcm overflow");
    return static_cast<i64>(value);
}

bool safe_mid(int speed, i64 left, i64 right, i64 denominator) {
    i128 residue = static_cast<i128>(speed) * (left + right) %
                   (static_cast<i128>(2) * denominator);
    if (residue < 0) residue += static_cast<i128>(2) * denominator;
    return static_cast<i128>(7) * residue >= denominator &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * denominator;
}

struct Cell {
    i64 left;
    i64 right;
    u32 failed;
    bool q_safe;
};

std::vector<Cell> build_cells(int r, i64& denominator) {
    std::vector<int> speeds(CHART.begin(), CHART.end());
    speeds.push_back(Q);
    if (r > 0) speeds.push_back(r);
    std::sort(speeds.begin(), speeds.end());
    speeds.erase(std::unique(speeds.begin(), speeds.end()), speeds.end());
    denominator = 1;
    for (int speed : speeds) denominator = exact_lcm(denominator, 14LL * speed);
    std::vector<i64> walls{0, denominator};
    for (int speed : speeds) {
        const i64 unit = denominator / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    }
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    std::vector<Cell> cells;
    cells.reserve(walls.size() - 1);
    for (std::size_t index = 0; index + 1 < walls.size(); ++index) {
        u32 failed = 0;
        for (int vertex = 0; vertex < static_cast<int>(CHART.size()); ++vertex) {
            if (!safe_mid(CHART[vertex], walls[index], walls[index + 1], denominator)) {
                failed |= u32{1} << vertex;
            }
        }
        cells.push_back({walls[index], walls[index + 1], failed,
                         safe_mid(Q, walls[index], walls[index + 1], denominator)});
    }
    return cells;
}

std::string decimal(i128 value) {
    if (value == 0) return "0";
    bool negative = value < 0;
    if (negative) value = -value;
    std::string answer;
    while (value > 0) {
        answer.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (negative) answer.push_back('-');
    std::reverse(answer.begin(), answer.end());
    return answer;
}

std::string labels(u32 mask) {
    std::string answer;
    for (int vertex = 0; vertex < static_cast<int>(CHART.size()); ++vertex) {
        if ((mask & (u32{1} << vertex)) == 0) continue;
        if (!answer.empty()) answer += ',';
        answer += std::to_string(CHART[vertex]);
    }
    return answer;
}

u32 mask_from_values(const std::array<int, 9>& values) {
    u32 answer = 0;
    for (int value : values) {
        const auto found = std::find(CHART.begin(), CHART.end(), value);
        require(found != CHART.end(), "hostile body left chart");
        answer |= u32{1} << std::distance(CHART.begin(), found);
    }
    return answer;
}

struct Profile {
    u32 body;
    i64 mass = 0;
    int components = 0;
    i128 limit_delta = 0;
    int cutoff = -1;
    int failures = 0;
    int equalities = 0;
};

}  // namespace

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    constexpr int n = static_cast<int>(CHART.size());
    const u32 full = (u32{1} << n) - 1;
    std::vector<Profile> profiles;
    for (u32 mask = 0; mask <= full; ++mask) {
        if (std::popcount(mask) == 9) profiles.push_back(Profile{mask});
    }
    i64 base_denominator = 0;
    const auto base = build_cells(-1, base_denominator);
    int largest_cutoff = 0;
    int nonstrict = 0;
    i128 minimum_limit_delta = 0;
    u32 minimum_limit_body = 0;
    for (Profile& profile : profiles) {
        std::vector<bool> included;
        included.reserve(base.size());
        for (const Cell& cell : base) {
            const bool good = cell.q_safe && (cell.failed & profile.body) == 0;
            included.push_back(good);
            if (good) profile.mass += cell.right - cell.left;
        }
        bool previous = included.back();
        for (bool current : included) {
            if (current && !previous) ++profile.components;
            previous = current;
        }
        profile.limit_delta = static_cast<i128>(54) * profile.mass -
                              static_cast<i128>(4) * base_denominator;
        if (minimum_limit_body == 0 || profile.limit_delta < minimum_limit_delta) {
            minimum_limit_delta = profile.limit_delta;
            minimum_limit_body = profile.body;
        }
        if (profile.limit_delta <= 0) {
            ++nonstrict;
            continue;
        }
        const i128 numerator = static_cast<i128>(54) * profile.components *
                               base_denominator;
        const i128 divisor = static_cast<i128>(7) * profile.limit_delta;
        profile.cutoff = static_cast<int>((numerator + divisor - 1) / divisor);
        largest_cutoff = std::max(largest_cutoff, profile.cutoff);
    }

    int scanned_r = 0;
    int first_failure_r = -1;
    u32 first_failure_body = 0;
    i128 first_failure_delta = 0;
    std::uint64_t finite_comparisons = 0;
    i128 minimum_finite_delta = 0;
    i64 minimum_finite_denominator = 1;
    int minimum_finite_r = -1;
    u32 minimum_finite_body = 0;
    i128 hostile_r6_delta = 0;
    const u32 hostile_body = mask_from_values(HOSTILE_BODY_VALUES);
    Fnv1a64 ledger;
    for (int r = 1; r < largest_cutoff; ++r) {
        if (r == Q || std::find(POOL.begin(), POOL.end(), r) != POOL.end()) continue;
        ++scanned_r;
        i64 denominator = 0;
        const auto cells = build_cells(r, denominator);
        std::vector<i64> zeta(std::size_t{1} << n, 0);
        for (const Cell& cell : cells) {
            if (!cell.q_safe || !safe_mid(r, cell.left, cell.right, denominator)) continue;
            zeta[cell.failed] += cell.right - cell.left;
        }
        for (int bit = 0; bit < n; ++bit) {
            for (u32 mask = 0; mask <= full; ++mask) {
                if ((mask & (u32{1} << bit)) != 0) {
                    zeta[mask] += zeta[mask ^ (u32{1} << bit)];
                }
            }
        }
        for (Profile& profile : profiles) {
            if (profile.cutoff <= 0 || r >= profile.cutoff) continue;
            ++finite_comparisons;
            const i64 mass = zeta[full ^ profile.body];
            const i128 delta = static_cast<i128>(63) * mass -
                               static_cast<i128>(4) * denominator;
            if (delta < 0) {
                ++profile.failures;
                if (first_failure_r < 0) {
                    first_failure_r = r;
                    first_failure_body = profile.body;
                    first_failure_delta = delta;
                }
            }
            if (delta == 0) ++profile.equalities;
            if (minimum_finite_r < 0 ||
                delta * minimum_finite_denominator <
                    minimum_finite_delta * denominator) {
                minimum_finite_delta = delta;
                minimum_finite_denominator = denominator;
                minimum_finite_r = r;
                minimum_finite_body = profile.body;
            }
            if (r == 6 && profile.body == hostile_body) hostile_r6_delta = delta;
            ledger.add_u64_le(static_cast<std::uint64_t>(r));
            ledger.add_u64_le(profile.body);
            ledger.add_u64_le(static_cast<std::uint64_t>(mass));
            ledger.add_u64_le(static_cast<std::uint64_t>(denominator));
        }
    }

    int universal = 0;
    int total_failures = 0;
    int total_equalities = 0;
    int min_cutoff = 0;
    int max_cutoff = 0;
    for (const Profile& profile : profiles) {
        total_failures += profile.failures;
        total_equalities += profile.equalities;
        if (profile.cutoff > 0 && profile.failures == 0 && profile.equalities == 0) {
            ++universal;
        }
        if (profile.cutoff > 0) {
            if (min_cutoff == 0 || profile.cutoff < min_cutoff) min_cutoff = profile.cutoff;
            max_cutoff = std::max(max_cutoff, profile.cutoff);
        }
    }
    require(profiles.size() == 48620, "chart body universe changed");
    require(base_denominator == INT64_C(91205797082400),
            "base denominator changed");
    require(nonstrict == 0 &&
                minimum_limit_delta == static_cast<i128>(517215373867152LL) &&
                labels(minimum_limit_body) == "8,80,85,88,95,143,145,168,193",
            "strict reserve profile changed");
    require(min_cutoff == 124 && max_cutoff == 448 && scanned_r == 416,
            "finite universe or cutoff range changed");
    require(universal == 48620 && total_failures == 0 && total_equalities == 0,
            "chart contains an uncertified body");
    require(finite_comparisons > 0 && minimum_finite_delta > 0 &&
                hostile_r6_delta > 0,
            "positive finite controls changed");
    std::cout << "CHART";
    for (int label : CHART) std::cout << ' ' << label;
    std::cout << '\n';
    std::cout << "SUMMARY N " << n
              << " BODIES " << profiles.size()
              << " NONSTRICT " << nonstrict
              << " MIN_LIMIT_DELTA " << decimal(minimum_limit_delta)
              << " MIN_LIMIT_BODY " << labels(minimum_limit_body)
              << " MIN_CUTOFF " << min_cutoff
              << " MAX_CUTOFF " << max_cutoff
              << " SCANNED_R " << scanned_r
              << " FINITE_COMPARISONS " << finite_comparisons
              << " UNIVERSAL " << universal
              << " TOTAL_FAILURES " << total_failures
              << " TOTAL_EQUALITIES " << total_equalities << '\n';
    std::cout << "MIN_FINITE_NORMALIZED_SURPLUS "
              << decimal(minimum_finite_delta) << '/' << minimum_finite_denominator
              << " R " << minimum_finite_r
              << " BODY " << labels(minimum_finite_body) << '\n';
    std::cout << "DEPTH8_RESERVE_COVER_BODY_R6_DELTA "
              << decimal(hostile_r6_delta)
              << " BODY " << labels(hostile_body) << '\n';
    std::cout << "SEMANTIC_FNV1A64_LE " << std::hex << ledger.value()
              << std::dec << '\n';
    std::cout << "VERDICT ALL_48620_CHART_BODIES_ALL_OUTSIDERS\n";
    if (first_failure_r >= 0) {
        std::cout << "FIRST_FAILURE R " << first_failure_r
                  << " BODY " << labels(first_failure_body)
                  << " DELTA " << decimal(first_failure_delta) << '\n';
    }
}
