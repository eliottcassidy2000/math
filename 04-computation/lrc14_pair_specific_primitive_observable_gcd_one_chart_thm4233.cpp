// Primary exact certificate for THM-4233.
//
// Reproduction (the output is required to be identical at both optimization
// levels):
//   g++ -std=c++20 -O0 -Wall -Wextra -pedantic \
//     lrc14_pair_specific_primitive_observable_gcd_one_chart_thm4233.cpp \
//     -o /tmp/thm4233-primary-O0
//   /tmp/thm4233-primary-O0
//   g++ -std=c++20 -O3 -Wall -Wextra -pedantic \
//     lrc14_pair_specific_primitive_observable_gcd_one_chart_thm4233.cpp \
//     -o /tmp/thm4233-primary-O3
//   /tmp/thm4233-primary-O3
//
// The circle is represented on the common integer grid N=14uv.  Every wall
// of the u-comb is v(14i+-1), and every wall of the v-comb is u(14j+-1).
// Thus each open cell has an exactly decidable midpoint state.  If a is the
// indicator of G_u intersect G_v and S is its safe length in ticks, then the
// centered primitive changes across a cell of length d by
//
//                         d (N a - S) / N^2.
//
// All comparisons and reductions below use integers.  No floating-point
// operation participates in the certificate.

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

constexpr i64 MAIN_U = 3713;
constexpr i64 MAIN_V = 5149;
constexpr i64 POOL_MAX = 290;

constexpr i64 EXPECTED_GRID = INT64_C(267655318);
constexpr u64 EXPECTED_RAW_POINTS = UINT64_C(17726);
constexpr u64 EXPECTED_UNIQUE_POINTS = UINT64_C(17724);
constexpr u64 EXPECTED_SAFE_CELLS = UINT64_C(7595);
constexpr u64 EXPECTED_UNSAFE_CELLS = UINT64_C(10128);
constexpr u64 EXPECTED_COMPONENTS = UINT64_C(7595);
constexpr i64 EXPECTED_SAFE_TICKS = INT64_C(196644720);
constexpr i128 EXPECTED_PRIMITIVE_ABS =
    static_cast<i128>(INT64_C(2057535171948));
constexpr i64 EXPECTED_MIN_TICK = INT64_C(207775767);
constexpr i64 EXPECTED_MAX_TICK = INT64_C(59879551);
constexpr u64 EXPECTED_LEDGER = UINT64_C(0x4725493ac81fe903);

constexpr u64 FNV_OFFSET = UINT64_C(0xcbf29ce484222325);
constexpr u64 FNV_PRIME = UINT64_C(0x100000001b3);

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

i128 gcd128(i128 left, i128 right) {
    if (left < 0) left = -left;
    if (right < 0) right = -right;
    while (right != 0) {
        const i128 remainder = left % right;
        left = right;
        right = remainder;
    }
    return left;
}

std::string decimal(i128 value) {
    if (value == 0) return "0";
    const bool negative = value < 0;
    if (negative) value = -value;
    std::string result;
    while (value != 0) {
        result.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (negative) result.push_back('-');
    std::reverse(result.begin(), result.end());
    return result;
}

struct Fraction {
    i128 numerator = 0;
    i128 denominator = 1;
};

Fraction fraction(i128 numerator, i128 denominator) {
    require(denominator != 0, "zero fraction denominator");
    if (denominator < 0) {
        numerator = -numerator;
        denominator = -denominator;
    }
    const i128 divisor = gcd128(numerator, denominator);
    return {numerator / divisor, denominator / divisor};
}

Fraction subtract(Fraction left, Fraction right) {
    return fraction(left.numerator * right.denominator -
                        right.numerator * left.denominator,
                    left.denominator * right.denominator);
}

bool equal(Fraction left, Fraction right) {
    return left.numerator * right.denominator ==
           right.numerator * left.denominator;
}

bool less(Fraction left, Fraction right) {
    return left.numerator * right.denominator <
           right.numerator * left.denominator;
}

std::string format(Fraction value) {
    return decimal(value.numerator) + "/" + decimal(value.denominator);
}

class Fnv1a64 {
  public:
    void add_u64_le(u64 word) {
        for (int shift = 0; shift < 64; shift += 8) {
            value_ ^= static_cast<std::uint8_t>(
                (word >> shift) & UINT64_C(0xff));
            value_ *= FNV_PRIME;
        }
    }

    u64 value() const { return value_; }

  private:
    u64 value_ = FNV_OFFSET;
};

std::string hex64(u64 value) {
    std::ostringstream out;
    out << std::hex << std::setfill('0') << std::setw(16) << value;
    return out.str();
}

i64 normalized_tick(i128 value, i64 grid) {
    i128 result = value % grid;
    if (result < 0) result += grid;
    return static_cast<i64>(result);
}

// Test the midpoint (left+right)/(2N) against
// G_speed={t: ||speed*t|| >= 1/14}.  Working modulo 2N avoids fractions.
bool midpoint_safe(i64 speed, i64 grid, i64 left, i64 right) {
    const i64 period = 2 * grid;
    i128 residue = static_cast<i128>(speed) * (left + right);
    residue %= period;
    if (residue < 0) residue += period;
    return 7 * residue >= grid && 7 * residue <= 13 * grid;
}

struct Cell {
    i64 left = 0;
    i64 right = 0;
    bool safe = false;
};

struct PairStats {
    i64 u = 0;
    i64 v = 0;
    i64 grid = 0;
    u64 raw_points = 0;
    u64 unique_points = 0;
    std::vector<i64> collision_ticks;
    u64 safe_cells = 0;
    u64 unsafe_cells = 0;
    u64 components = 0;
    i64 safe_ticks = 0;
    i128 primitive_min_raw = 0;
    i64 primitive_min_tick = 0;
    i128 primitive_max_raw = 0;
    i64 primitive_max_tick = 0;
    i128 raw_denominator = 1;
    Fraction beta;
    Fraction omega;
    u64 ledger = 0;
};

PairStats exact_pair(i64 input_u, i64 input_v) {
    require(input_u > 0 && input_v > 0 && input_u != input_v,
            "pair must contain distinct positive speeds");
    i64 u = input_u;
    i64 v = input_v;
    if (u > v) std::swap(u, v);
    require(std::gcd(u, v) == 1,
            "the primitive chart expects coprime speeds");

    const i128 grid128 = static_cast<i128>(14) * u * v;
    require(grid128 <= INT64_MAX, "common grid overflows int64");
    const i64 grid = static_cast<i64>(grid128);

    std::vector<i64> points;
    points.reserve(static_cast<std::size_t>(2 * (u + v) + 2));
    points.push_back(0);
    points.push_back(grid);
    for (i64 index = 0; index < u; ++index) {
        points.push_back(normalized_tick(
            static_cast<i128>(v) * (14 * index - 1), grid));
        points.push_back(normalized_tick(
            static_cast<i128>(v) * (14 * index + 1), grid));
    }
    for (i64 index = 0; index < v; ++index) {
        points.push_back(normalized_tick(
            static_cast<i128>(u) * (14 * index - 1), grid));
        points.push_back(normalized_tick(
            static_cast<i128>(u) * (14 * index + 1), grid));
    }

    PairStats result;
    result.u = u;
    result.v = v;
    result.grid = grid;
    result.raw_points = points.size();

    std::sort(points.begin(), points.end());
    for (std::size_t index = 1; index < points.size(); ++index) {
        if (points[index] == points[index - 1]) {
            result.collision_ticks.push_back(points[index]);
        }
    }
    points.erase(std::unique(points.begin(), points.end()), points.end());
    result.unique_points = points.size();
    require(points.front() == 0 && points.back() == grid,
            "linearized circle lost a boundary");

    std::vector<Cell> cells;
    cells.reserve(points.size() - 1);
    Fnv1a64 ledger;
    for (std::size_t index = 1; index < points.size(); ++index) {
        const i64 left = points[index - 1];
        const i64 right = points[index];
        require(left < right, "zero-length unique cell");
        const bool safe = midpoint_safe(u, grid, left, right) &&
                          midpoint_safe(v, grid, left, right);
        cells.push_back({left, right, safe});
        ledger.add_u64_le(static_cast<u64>(left));
        ledger.add_u64_le(static_cast<u64>(right));
        ledger.add_u64_le(static_cast<u64>(safe));
        if (safe) {
            ++result.safe_cells;
            result.safe_ticks += right - left;
        } else {
            ++result.unsafe_cells;
        }
    }
    result.ledger = ledger.value();

    for (std::size_t index = 0; index < cells.size(); ++index) {
        const std::size_t previous =
            (index + cells.size() - 1) % cells.size();
        if (cells[index].safe && !cells[previous].safe) ++result.components;
    }

    i128 primitive = 0;
    for (const Cell& cell : cells) {
        primitive += static_cast<i128>(cell.right - cell.left) *
                     (static_cast<i128>(grid) * cell.safe -
                      result.safe_ticks);
        if (primitive < result.primitive_min_raw) {
            result.primitive_min_raw = primitive;
            result.primitive_min_tick = cell.right;
        }
        if (primitive > result.primitive_max_raw) {
            result.primitive_max_raw = primitive;
            result.primitive_max_tick = cell.right;
        }
    }
    require(primitive == 0, "centered primitive failed to close");

    result.raw_denominator = static_cast<i128>(grid) * grid;
    result.beta = fraction(result.safe_ticks, grid);
    result.omega = fraction(result.primitive_max_raw -
                                result.primitive_min_raw,
                            result.raw_denominator);
    return result;
}

void audit_main_pair(const PairStats& stats) {
    require(stats.u == MAIN_U && stats.v == MAIN_V,
            "main pair changed");
    require(stats.grid == EXPECTED_GRID, "main grid changed");
    require(stats.raw_points == EXPECTED_RAW_POINTS,
            "raw endpoint count changed");
    require(stats.unique_points == EXPECTED_UNIQUE_POINTS,
            "unique endpoint count changed");
    require(stats.raw_points - stats.unique_points == 2,
            "endpoint duplicate excess changed");
    require(stats.collision_ticks ==
                std::vector<i64>{INT64_C(95591185), INT64_C(172064133)},
            "collision ticks changed");
    require(stats.safe_cells == EXPECTED_SAFE_CELLS,
            "safe-cell count changed");
    require(stats.unsafe_cells == EXPECTED_UNSAFE_CELLS,
            "unsafe-cell count changed");
    require(stats.components == EXPECTED_COMPONENTS,
            "safe-component count changed");
    require(stats.safe_ticks == EXPECTED_SAFE_TICKS,
            "safe mass changed");
    require(equal(stats.beta, fraction(98322360, 133827659)),
            "main beta changed");
    require(equal(subtract(stats.beta, fraction(66, 91)),
                  fraction(16387038, 1739759567)),
            "main beta margin changed");
    require(stats.primitive_min_raw == -EXPECTED_PRIMITIVE_ABS &&
                stats.primitive_min_tick == EXPECTED_MIN_TICK,
            "main primitive minimum changed");
    require(stats.primitive_max_raw == EXPECTED_PRIMITIVE_ABS &&
                stats.primitive_max_tick == EXPECTED_MAX_TICK,
            "main primitive maximum changed");
    require(stats.raw_denominator ==
                static_cast<i128>(INT64_C(71639369253681124)),
            "main primitive denominator changed");
    require(equal(stats.omega,
                  fraction(277071798, INT64_C(4823550313337))),
            "main primitive oscillation changed");
    require(stats.ledger == EXPECTED_LEDGER,
            "main semantic ledger changed");
}

void audit_controls(const PairStats& hostile,
                    const PairStats& simple,
                    const PairStats& fibonacci,
                    Fraction threshold) {
    require(equal(hostile.beta, fraction(66, 91)),
            "1:13 beta control changed");
    require(equal(hostile.omega, fraction(990, 8281)),
            "1:13 oscillation control changed");
    require(equal(simple.beta, fraction(11, 14)),
            "1:2 beta control changed");
    require(equal(simple.omega, fraction(11, 98)),
            "1:2 oscillation control changed");
    require(equal(fibonacci.omega,
                  fraction(42967355, INT64_C(553336008694))),
            "Fibonacci predecessor oscillation changed");
    require(!less(fibonacci.omega, threshold) &&
                !equal(fibonacci.omega, threshold),
            "Fibonacci predecessor unexpectedly passes the transfer cutoff");
}

}  // namespace

int main() {
    try {
        const PairStats main_pair = exact_pair(MAIN_U, MAIN_V);
        const PairStats hostile = exact_pair(1, 13);
        const PairStats simple = exact_pair(1, 2);
        const PairStats fibonacci = exact_pair(2584, 4181);

        audit_main_pair(main_pair);
        const Fraction baseline = fraction(66, 91);
        const Fraction threshold = fraction(1650,
                                            static_cast<i128>(8281) * 3467);
        const Fraction beta_margin = subtract(main_pair.beta, baseline);
        const Fraction threshold_margin = subtract(threshold, main_pair.omega);
        require(beta_margin.numerator > 0,
                "main pair misses the beta baseline");
        require(threshold_margin.numerator > 0,
                "main pair misses the THM-4228 oscillation cutoff");
        require(equal(threshold, fraction(1650, 28710227)),
                "THM-4228 threshold reduction changed");
        require(equal(threshold_margin,
                      fraction(INT64_C(82934716896),
                               static_cast<i128>(INT64_C(2826229070241355051)))),
                "THM-4228 cutoff margin changed");
        audit_controls(hostile, simple, fibonacci, threshold);
        require(std::gcd(MAIN_U, MAIN_V) == 1,
                "main pair is not primitive");
        require(MAIN_U > POOL_MAX && MAIN_V > POOL_MAX,
                "main pair does not lie outside the fixed pool");

        std::cout << "THM-4233 primary integer common-grid sweep\n";
        std::cout << "pair=" << main_pair.u << ',' << main_pair.v
                  << " gcd=" << std::gcd(main_pair.u, main_pair.v)
                  << " pool_max=" << POOL_MAX
                  << " outside_pool="
                  << (main_pair.u > POOL_MAX && main_pair.v > POOL_MAX)
                  << '\n';
        std::cout << "grid=" << main_pair.grid
                  << " raw_points=" << main_pair.raw_points
                  << " unique_points=" << main_pair.unique_points
                  << " duplicate_excess="
                  << main_pair.raw_points - main_pair.unique_points
                  << " collision_ticks=";
        for (std::size_t index = 0;
             index < main_pair.collision_ticks.size(); ++index) {
            if (index != 0) std::cout << ',';
            std::cout << main_pair.collision_ticks[index];
        }
        std::cout << " safe_cells=" << main_pair.safe_cells
                  << " unsafe_cells=" << main_pair.unsafe_cells
                  << " components=" << main_pair.components << '\n';
        std::cout << "safe_ticks=" << main_pair.safe_ticks
                  << " beta=" << format(main_pair.beta)
                  << " beta_minus_66_over_91=" << format(beta_margin)
                  << '\n';
        std::cout << "primitive_min_raw="
                  << decimal(main_pair.primitive_min_raw)
                  << " at_tick=" << main_pair.primitive_min_tick
                  << " primitive_max_raw="
                  << decimal(main_pair.primitive_max_raw)
                  << " at_tick=" << main_pair.primitive_max_tick
                  << " raw_den=" << decimal(main_pair.raw_denominator)
                  << '\n';
        std::cout << "omega=" << format(main_pair.omega)
                  << " thm4228_threshold=" << format(threshold)
                  << " threshold_margin=" << format(threshold_margin)
                  << '\n';
        std::cout << "semantic_ledger_fnv1a64_le_words="
                  << hex64(main_pair.ledger) << '\n';
        std::cout << "hostile_1_13_beta=" << format(hostile.beta)
                  << " hostile_1_13_omega=" << format(hostile.omega)
                  << " fibonacci_2584_4181_omega="
                  << format(fibonacci.omega)
                  << " fibonacci_predecessor_pass="
                  << (less(fibonacci.omega, threshold) ||
                      equal(fibonacci.omega, threshold))
                  << '\n';
        std::cout << "transfer=PASS beta>=66/91,omega<=(1650/8281)/3467,"
                     "THM-4228-E3467,gcd-one,outside-pool\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM-4233 primary audit failed: "
                  << error.what() << '\n';
        return 1;
    }
}
