// Independent direct-interval audit for THM-4214's zero-newcomer layer.
//
// Unlike the primary failure-mask/zeta path, this implementation enumerates
// every eleven-body and intersects its labelled safe-comb intervals by an
// event sweep.  It has no Boolean-lattice transform or shared mass table.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
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

    std::uint64_t bodies = 0;
    std::uint64_t failures = 0;
    std::uint64_t equalities = 0;
    std::uint64_t event_pairs = 0;
    i128 minimum_delta = 0;
    i64 minimum_mass = 0;
    u32 minimum_body = 0;
    Fnv1a64 ledger;

    auto direct_mass = [&](u32 body, int cardinality, std::uint64_t& pair_counter) {
        std::vector<std::pair<i64, int>> events;
        for (int vertex = 0; vertex < n; ++vertex) {
            if ((body & (u32{1} << vertex)) == 0) continue;
            const int speed = CHART[vertex];
            const i64 unit = denominator / (14LL * speed);
            for (int tooth = 0; tooth < speed; ++tooth) {
                events.emplace_back((14LL * tooth + 1) * unit, +1);
                events.emplace_back((14LL * tooth + 13) * unit, -1);
            }
        }
        pair_counter += events.size() / 2;
        std::sort(events.begin(), events.end());

        int active = 0;
        i64 previous = 0;
        i64 mass = 0;
        std::size_t index = 0;
        while (index < events.size()) {
            const i64 coordinate = events[index].first;
            if (active == cardinality) mass += coordinate - previous;
            int delta = 0;
            while (index < events.size() && events[index].first == coordinate) {
                delta += events[index].second;
                ++index;
            }
            active += delta;
            require(active >= 0 && active <= cardinality,
                    "active interval count invalid");
            previous = coordinate;
        }
        require(active == 0, "safe-comb sweep did not close");
        return mass;
    };

    for (u32 body = 0; body <= full; ++body) {
        if (std::popcount(body) != 11) continue;
        ++bodies;
        const i64 mass = direct_mass(body, 11, event_pairs);

        const i128 delta = static_cast<i128>(63) * mass -
                           static_cast<i128>(4) * denominator;
        if (delta < 0) ++failures;
        if (delta == 0) ++equalities;
        if (minimum_body == 0 || delta < minimum_delta) {
            minimum_delta = delta;
            minimum_mass = mass;
            minimum_body = body;
        }
        ledger.add_u64_le(body);
        ledger.add_u64_le(static_cast<std::uint64_t>(mass));
        ledger.add_u64_le(static_cast<std::uint64_t>(delta));
    }

    std::uint64_t nine_event_pairs = 0;
    i128 minimum_nine_delta22 = 0;
    i64 minimum_nine_mass = 0;
    u32 minimum_nine_body = 0;
    for (u32 body = 0; body <= full; ++body) {
        if (std::popcount(body) != 9) continue;
        const i64 mass = direct_mass(body, 9, nine_event_pairs);
        const i128 delta22 = static_cast<i128>(63) * mass -
                             static_cast<i128>(22) * denominator;
        if (minimum_nine_body == 0 || delta22 < minimum_nine_delta22) {
            minimum_nine_delta22 = delta22;
            minimum_nine_mass = mass;
            minimum_nine_body = body;
        }
    }

    require(denominator == INT64_C(18241159416480), "denominator changed");
    require(bodies == 31824, "body count changed");
    require(failures == 0 && equalities == 0, "zero-newcomer layer lost strict safety");
    require(minimum_delta == static_cast<i128>(135334677938922),
            "zero-newcomer minimum delta changed");
    require(minimum_mass == INT64_C(3306338342934),
            "zero-newcomer minimum mass changed");
    require(labels(minimum_body) == "16,40,42,85,88,95,145,193,240,252,286",
            "zero-newcomer minimizer changed");
    require(minimum_nine_delta22 == -static_cast<i128>(131225923051902),
            "Bonferroni hostile delta changed");
    require(minimum_nine_mass == INT64_C(4286977525566),
            "Bonferroni hostile mass changed");
    require(labels(minimum_nine_body) == "8,80,85,88,95,143,145,168,193",
            "Bonferroni hostile body changed");
    require(ledger.value() == UINT64_C(0xf594af0b446a5bf2),
            "semantic ledger changed");

    std::cout << "METHOD per-body labelled safe-comb event sweep; no zeta transform\n";
    std::cout << "GEOMETRY DEN " << denominator
              << " EVENT_PAIRS " << event_pairs << '\n';
    std::cout << "ZERO_NEWCOMER_LAYER BODIES " << bodies
              << " FAILURES " << failures
              << " EQUALITIES " << equalities
              << " MIN_DELTA " << decimal(minimum_delta)
              << " MIN_MASS " << minimum_mass
              << " MIN_BODY " << labels(minimum_body) << '\n';
    std::cout << "NINE_BODY_BONFERRONI EVENT_PAIRS " << nine_event_pairs
              << " MIN_DELTA_63M_MINUS_22D " << decimal(minimum_nine_delta22)
              << " MIN_MASS " << minimum_nine_mass
              << " MIN_BODY " << labels(minimum_nine_body) << '\n';
    std::cout << std::hex << "SEMANTIC_FNV1A64_LE " << ledger.value() << std::dec << '\n';
    std::cout << "VERDICT "
              << (failures == 0 && equalities == 0
                      ? "ALL_ZERO_NEWCOMER_ELEVEN_BODIES_STRICTLY_SAFE"
                      : "ZERO_NEWCOMER_LAYER_NOT_CLOSED")
              << '\n';
}
