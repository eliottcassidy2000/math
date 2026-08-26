// Primary exact audit for THM-4214's zero-newcomer Pascal layer.
//
// The inherited THM-4191/4211 implications cover the one- and two-newcomer
// layers.  This program audits the only missing layer: all eleven-subsets of
// the displayed eighteen-label chart.  It integrates once on the full chart,
// accumulates exact failure-mask atoms, and applies a subset zeta transform.

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

using i64 = std::int64_t;
using u32 = std::uint32_t;
using i128 = __int128_t;

namespace {

constexpr std::array<int, 18> CHART = {
    8,16,40,42,80,84,85,88,95,120,126,143,145,168,193,240,252,286};
constexpr std::uint64_t FNV_OFFSET = UINT64_C(0xcbf29ce484222325);
constexpr std::uint64_t FNV_PRIME = UINT64_C(0x100000001b3);

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

std::uint64_t choose(int n, int k) {
    if (k < 0 || k > n) return 0;
    k = std::min(k, n - k);
    std::uint64_t answer = 1;
    for (int j = 1; j <= k; ++j) {
        answer = answer * static_cast<std::uint64_t>(n - k + j) /
                 static_cast<std::uint64_t>(j);
    }
    return answer;
}

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

}  // namespace

int main() {
    constexpr int n = static_cast<int>(CHART.size());
    const u32 full = (u32{1} << n) - 1;

    i64 denominator = 1;
    for (int speed : CHART) denominator = exact_lcm(denominator, 14LL * speed);

    std::vector<i64> walls{0, denominator};
    for (int speed : CHART) {
        const i64 unit = denominator / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    }
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());

    std::vector<i64> zeta(std::size_t{1} << n, 0);
    for (std::size_t index = 0; index + 1 < walls.size(); ++index) {
        u32 failed = 0;
        for (int vertex = 0; vertex < n; ++vertex) {
            if (!safe_mid(CHART[vertex], walls[index], walls[index + 1], denominator)) {
                failed |= u32{1} << vertex;
            }
        }
        zeta[failed] += walls[index + 1] - walls[index];
    }
    for (int bit = 0; bit < n; ++bit) {
        for (u32 mask = 0; mask <= full; ++mask) {
            if ((mask & (u32{1} << bit)) != 0) {
                zeta[mask] += zeta[mask ^ (u32{1} << bit)];
            }
        }
    }

    std::uint64_t bodies = 0;
    std::uint64_t failures = 0;
    std::uint64_t equalities = 0;
    i128 minimum_delta = 0;
    i64 minimum_mass = 0;
    u32 minimum_body = 0;
    std::uint64_t frontier_facing = 0;
    std::uint64_t clock_trivial = 0;
    u32 clock_trivial_body = 0;
    i128 minimum_nine_delta22 = 0;
    i64 minimum_nine_mass = 0;
    u32 minimum_nine_body = 0;
    Fnv1a64 ledger;
    for (u32 body = 0; body <= full; ++body) {
        if (std::popcount(body) != 11) continue;
        ++bodies;
        const i64 mass = zeta[full ^ body];
        const i128 delta = static_cast<i128>(63) * mass -
                           static_cast<i128>(4) * denominator;
        if (delta < 0) ++failures;
        if (delta == 0) ++equalities;
        if (minimum_body == 0 || delta < minimum_delta) {
            minimum_delta = delta;
            minimum_mass = mass;
            minimum_body = body;
        }
        bool has_multiple_six = false;
        for (int vertex = 0; vertex < n; ++vertex) {
            if ((body & (u32{1} << vertex)) != 0 && CHART[vertex] % 6 == 0) {
                has_multiple_six = true;
            }
        }
        if (has_multiple_six) {
            ++frontier_facing;
        } else {
            ++clock_trivial;
            clock_trivial_body = body;
        }
        ledger.add_u64_le(body);
        ledger.add_u64_le(static_cast<std::uint64_t>(mass));
        ledger.add_u64_le(static_cast<std::uint64_t>(delta));
    }
    for (u32 body = 0; body <= full; ++body) {
        if (std::popcount(body) != 9) continue;
        const i64 mass = zeta[full ^ body];
        const i128 delta22 = static_cast<i128>(63) * mass -
                             static_cast<i128>(22) * denominator;
        if (minimum_nine_body == 0 || delta22 < minimum_nine_delta22) {
            minimum_nine_delta22 = delta22;
            minimum_nine_mass = mass;
            minimum_nine_body = body;
        }
    }

    const std::uint64_t zero_layer = choose(18, 11);
    const std::uint64_t one_layer = 2 * choose(18, 10);
    const std::uint64_t two_layer = choose(18, 9);
    const std::uint64_t complete_layer = choose(20, 11);

    require(denominator == INT64_C(18241159416480), "denominator changed");
    require(walls.size() == 4372, "wall census changed");
    require(bodies == zero_layer, "zero-newcomer body count changed");
    require(failures == 0 && equalities == 0, "zero-newcomer layer lost strict safety");
    require(minimum_delta == static_cast<i128>(135334677938922),
            "zero-newcomer minimum delta changed");
    require(minimum_mass == INT64_C(3306338342934),
            "zero-newcomer minimum mass changed");
    require(labels(minimum_body) == "16,40,42,85,88,95,145,193,240,252,286",
            "zero-newcomer minimizer changed");
    require(zero_layer + one_layer + two_layer == complete_layer,
            "Pascal decomposition failed");
    require(frontier_facing == 31823 && clock_trivial == 1, "divisor split failed");
    require(labels(clock_trivial_body) == "8,16,40,80,85,88,95,143,145,193,286",
            "clock-trivial body changed");
    require(minimum_nine_delta22 == -static_cast<i128>(131225923051902),
            "Bonferroni hostile delta changed");
    require(minimum_nine_mass == INT64_C(4286977525566),
            "Bonferroni hostile mass changed");
    require(labels(minimum_nine_body) == "8,80,85,88,95,143,145,168,193",
            "Bonferroni hostile body changed");
    require(ledger.value() == UINT64_C(0xf594af0b446a5bf2),
            "semantic ledger changed");

    std::cout << "CHART";
    for (int speed : CHART) std::cout << ' ' << speed;
    std::cout << '\n';
    std::cout << "GEOMETRY DEN " << denominator
              << " WALLS " << walls.size()
              << " CELLS " << walls.size() - 1 << '\n';
    std::cout << "ZERO_NEWCOMER_LAYER BODIES " << bodies
              << " FAILURES " << failures
              << " EQUALITIES " << equalities
              << " MIN_DELTA " << decimal(minimum_delta)
              << " MIN_MASS " << minimum_mass
              << " MIN_BODY " << labels(minimum_body) << '\n';
    std::cout << "DIVISOR_CONTROL HAS_MULTIPLE_6 " << frontier_facing
              << " NO_MULTIPLE_6 " << clock_trivial
              << " NO_MULTIPLE_6_BODY " << labels(clock_trivial_body) << '\n';
    std::cout << "NINE_BODY_BONFERRONI MIN_DELTA_63M_MINUS_22D "
              << decimal(minimum_nine_delta22)
              << " MIN_MASS " << minimum_nine_mass
              << " MIN_BODY " << labels(minimum_nine_body) << '\n';
    std::cout << "PASCAL_LAYER ZERO " << zero_layer
              << " ONE " << one_layer
              << " TWO " << two_layer
              << " TOTAL " << zero_layer + one_layer + two_layer
              << " BINOM_20_11 " << complete_layer << '\n';
    std::cout << std::hex << "SEMANTIC_FNV1A64_LE " << ledger.value() << std::dec << '\n';
    std::cout << "VERDICT "
              << (failures == 0 && equalities == 0
                      ? "ALL_ZERO_NEWCOMER_ELEVEN_BODIES_STRICTLY_SAFE"
                      : "ZERO_NEWCOMER_LAYER_NOT_CLOSED")
              << '\n';
}
