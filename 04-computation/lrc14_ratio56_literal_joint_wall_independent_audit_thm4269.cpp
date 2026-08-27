#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

// Clean-room literal-wall audit of the ratio-5:6 common-deck bridge.
//
// This source does not use the primitive-prefix antiderivative from the
// primary.  For every g=59,...,128, it constructs the safe intervals of 5g
// and 6g on lcm(D,420g), intersects them, and sweeps the result against the
// exact fixed-pool wall atoms.  Candidate selection uses a streaming top-k
// priority queue instead of sorting the full repair universe.

namespace {

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

constexpr std::array<int, 30> POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
constexpr i64 D = INT64_C(18241159416480);
constexpr u64 ORDER_SEED = UINT64_C(0x4245422842334245);
constexpr std::size_t KEEP = 16384;
constexpr u64 REPAIR_UNIVERSE = UINT64_C(5852925);
constexpr u64 BODY_UNIVERSE = UINT64_C(14307150);
constexpr int G0 = 59;
constexpr int G1 = 128;
constexpr int SCALE_COUNT = G1 - G0 + 1;
constexpr std::array<std::size_t, 3> BUDGETS = {4096, 8192, 16384};
constexpr std::array<u64, 3> EXPECTED_COMMON = {1630, 3158, 6416};
constexpr std::array<u64, 3> EXPECTED_FAILURES = {116, 6, 0};
constexpr std::array<u64, 3> EXPECTED_COMMON_FNV = {
    UINT64_C(0xe0ffc2e129f24678),
    UINT64_C(0x2c493b7ae8a8a847),
    UINT64_C(0x9749384b49c2bba0),
};
constexpr std::array<u64, 3> EXPECTED_CHECKS = {
    UINT64_C(486579464), UINT64_C(486631828), UINT64_C(486639933),
};
constexpr std::array<u64, 3> EXPECTED_MAX_CHECKS = {1630, 3158, 6119};
constexpr std::array<u32, 5> EXPECTED_ADDITIONS = {
    UINT32_C(0x1320001b), UINT32_C(0x3024004b), UINT32_C(0x01608381),
    UINT32_C(0x100c0781), UINT32_C(0x0a42c081),
};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

std::string dec(i128 value) {
    if (value == 0) return "0";
    const bool negative = value < 0;
    u128 magnitude = negative ? static_cast<u128>(-(value + 1)) + 1
                              : static_cast<u128>(value);
    std::string out;
    while (magnitude != 0) {
        out.push_back(static_cast<char>('0' + magnitude % 10));
        magnitude /= 10;
    }
    if (negative) out.push_back('-');
    std::reverse(out.begin(), out.end());
    return out;
}

std::string hex64(u64 value) {
    std::ostringstream out;
    out << std::hex << std::nouppercase << std::setw(16)
        << std::setfill('0') << value;
    return out.str();
}

class Fnv64 {
  public:
    void byte(std::uint8_t value) {
        state_ ^= value;
        state_ *= UINT64_C(0x100000001b3);
    }
    void word(u64 value) {
        for (unsigned j = 0; j < 8; ++j) {
            byte(static_cast<std::uint8_t>(value >> (8 * j)));
        }
    }
    void signed128(i128 value) {
        const u128 bits = static_cast<u128>(value);
        word(static_cast<u64>(bits));
        word(static_cast<u64>(bits >> 64));
    }
    u64 value() const { return state_; }

  private:
    u64 state_ = UINT64_C(0xcbf29ce484222325);
};

u64 splitmix64(u64 x) {
    x += UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

u32 next_mask(u32 x) {
    const u32 c = x & (~x + 1u);
    const u32 r = x + c;
    if (r == 0) return 0;
    return r | (((r ^ x) >> 2) / c);
}

i64 exact_lcm(i64 a, i64 b) {
    const i64 divisor = std::gcd(a, b);
    const i128 answer = static_cast<i128>(a / divisor) * b;
    require(answer <= std::numeric_limits<i64>::max(), "lcm overflow");
    return static_cast<i64>(answer);
}

std::string labels(u32 mask) {
    std::ostringstream out;
    bool first = true;
    for (int bit = 0; bit < 30; ++bit) {
        if ((mask & (u32{1} << bit)) == 0) continue;
        if (!first) out << ',';
        first = false;
        out << POOL[bit];
    }
    return out.str();
}

struct RankedRepair {
    u64 key = 0;
    u32 mask = 0;
    bool operator<(const RankedRepair& other) const {
        return key != other.key ? key < other.key : mask < other.mask;
    }
};

std::vector<u32> select_candidates() {
    std::priority_queue<RankedRepair> best;
    u64 count = 0;
    u32 mask = (u32{1} << 8) - 1;
    const u32 limit = u32{1} << 30;
    while (mask != 0 && mask < limit) {
        ++count;
        const RankedRepair item{splitmix64(static_cast<u64>(mask) ^ ORDER_SEED),
                                mask};
        if (best.size() < KEEP) {
            best.push(item);
        } else if (item < best.top()) {
            best.pop();
            best.push(item);
        }
        const u32 following = next_mask(mask);
        if (following <= mask) break;
        mask = following;
    }
    require(count == REPAIR_UNIVERSE, "C(30,8) census mismatch");
    std::vector<RankedRepair> ordered;
    ordered.reserve(KEEP);
    while (!best.empty()) {
        ordered.push_back(best.top());
        best.pop();
    }
    std::sort(ordered.begin(), ordered.end());
    std::vector<u32> result;
    result.reserve(KEEP);
    for (const RankedRepair& item : ordered) result.push_back(item.mask);
    return result;
}

struct PoolCell {
    i64 left = 0;
    i64 right = 0;
    u32 failed = 0;
};

bool pool_speed_safe(int speed, i64 left, i64 right) {
    const i128 twice_phase =
        (static_cast<i128>(speed) * (left + right)) % (2 * D);
    return 7 * twice_phase >= D &&
           7 * twice_phase <= 13 * static_cast<i128>(D);
}

std::vector<PoolCell> pool_cells() {
    i64 denominator = 1;
    for (int speed : POOL) denominator = exact_lcm(denominator, 14LL * speed);
    require(denominator == D, "pool denominator mismatch");

    std::vector<i64> walls{0, D};
    for (int speed : POOL) {
        const i64 quantum = D / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * quantum);
            walls.push_back((14LL * tooth + 13) * quantum);
        }
    }
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    require(walls.size() == 7134, "pool wall census mismatch");

    std::vector<PoolCell> result;
    result.reserve(walls.size() - 1);
    for (std::size_t j = 0; j + 1 < walls.size(); ++j) {
        u32 failed = 0;
        for (int bit = 0; bit < 30; ++bit) {
            if (!pool_speed_safe(POOL[bit], walls[j], walls[j + 1])) {
                failed |= u32{1} << bit;
            }
        }
        result.push_back({walls[j], walls[j + 1], failed});
    }
    require(result.front().failed == (u32{1} << 30) - 1,
            "origin-neighbourhood control failed");
    require(result.back().failed == (u32{1} << 30) - 1,
            "terminal origin-neighbourhood control failed");
    return result;
}

struct Interval {
    i64 left = 0;
    i64 right = 0;
};

std::vector<Interval> speed_safe_intervals(int speed, i64 denominator) {
    const i64 wall_denominator = 14LL * speed;
    require(denominator % wall_denominator == 0, "nonintegral speed wall");
    const i64 quantum = denominator / wall_denominator;
    std::vector<Interval> result;
    result.reserve(speed);
    for (int tooth = 0; tooth < speed; ++tooth) {
        result.push_back({(14LL * tooth + 1) * quantum,
                          (14LL * tooth + 13) * quantum});
    }
    return result;
}

std::vector<Interval> literal_joint_intervals(int g, i64 denominator) {
    const auto left_set = speed_safe_intervals(5 * g, denominator);
    const auto right_set = speed_safe_intervals(6 * g, denominator);
    std::vector<Interval> result;
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < left_set.size() && j < right_set.size()) {
        const i64 left = std::max(left_set[i].left, right_set[j].left);
        const i64 right = std::min(left_set[i].right, right_set[j].right);
        if (left < right) {
            if (!result.empty() && result.back().right == left) {
                result.back().right = right;
            } else {
                result.push_back({left, right});
            }
        }
        if (left_set[i].right < right_set[j].right) {
            ++i;
        } else if (right_set[j].right < left_set[i].right) {
            ++j;
        } else {
            ++i;
            ++j;
        }
    }
    require(!result.empty(), "empty literal joint safe set");
    return result;
}

struct PrimitiveControl {
    std::vector<Interval> arcs;
    i64 safe_ticks = 0;
    i64 centered_min = 0;
    i64 centered_max = 0;
};

PrimitiveControl primitive_control() {
    PrimitiveControl result;
    result.arcs = literal_joint_intervals(1, 420);
    constexpr std::array<std::pair<i64, i64>, 10> expected = {
        std::pair<i64, i64>{6, 65}, {75, 78}, {90, 135}, {145, 162},
        {174, 205}, {215, 246}, {258, 275}, {285, 330}, {342, 345},
        {355, 414},
    };
    require(result.arcs.size() == expected.size(),
            "primitive arc count changed");
    for (std::size_t i = 0; i < expected.size(); ++i) {
        require(result.arcs[i].left == expected[i].first &&
                    result.arcs[i].right == expected[i].second,
                "primitive arc ledger changed");
        result.safe_ticks += result.arcs[i].right - result.arcs[i].left;
    }
    require(result.safe_ticks == 310, "primitive safe length changed");

    bool first = true;
    std::vector<i64> walls{0, 420};
    for (const Interval& arc : result.arcs) {
        walls.push_back(arc.left);
        walls.push_back(arc.right);
    }
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    for (i64 wall : walls) {
        i64 prefix = 0;
        for (const Interval& arc : result.arcs) {
            prefix += std::max<i64>(0, std::min(wall, arc.right) - arc.left);
        }
        const i64 centered = 420 * prefix - result.safe_ticks * wall;
        if (first) {
            first = false;
            result.centered_min = result.centered_max = centered;
        } else {
            result.centered_min = std::min(result.centered_min, centered);
            result.centered_max = std::max(result.centered_max, centered);
        }
    }
    require(result.centered_min == -4630 && result.centered_max == 4630,
            "primitive centered discrepancy changed");
    return result;
}

using WeightVector = std::array<i64, SCALE_COUNT>;

struct LiteralWeights {
    std::unordered_map<u32, WeightVector> by_failed;
    std::array<i64, SCALE_COUNT> denominator{};
    std::array<std::size_t, SCALE_COUNT> joint_components{};
};

LiteralWeights literal_weights(const std::vector<PoolCell>& cells) {
    LiteralWeights result;
    result.by_failed.reserve(cells.size());
    for (const PoolCell& cell : cells) {
        result.by_failed.try_emplace(cell.failed, WeightVector{});
    }

    for (int si = 0; si < SCALE_COUNT; ++si) {
        const int g = G0 + si;
        const i64 denominator = exact_lcm(D, 420LL * g);
        result.denominator[si] = denominator;
        const i64 pool_factor = denominator / D;
        const auto joint = literal_joint_intervals(g, denominator);
        result.joint_components[si] = joint.size();

        std::size_t cursor = 0;
        for (const PoolCell& cell : cells) {
            const i64 left = cell.left * pool_factor;
            const i64 right = cell.right * pool_factor;
            while (cursor < joint.size() && joint[cursor].right <= left) ++cursor;
            std::size_t k = cursor;
            i64 mass = 0;
            while (k < joint.size() && joint[k].left < right) {
                mass += std::max<i64>(0, std::min(right, joint[k].right) -
                                         std::max(left, joint[k].left));
                if (joint[k].right >= right) break;
                ++k;
            }
            result.by_failed[cell.failed][si] += mass;
        }

        i64 total = 0;
        for (const auto& entry : result.by_failed) total += entry.second[si];
        require(static_cast<i128>(42) * total ==
                    static_cast<i128>(31) * denominator,
                "literal 5:6 joint-density control failed");
        require(joint.size() == static_cast<std::size_t>(10 * g),
                "literal joint-component census changed");
    }
    return result;
}

struct DeckResult {
    std::array<std::vector<u32>, BUDGETS.size()> common;
    Fnv64 matrix_fnv;
    i128 least_margin = 0;
    int least_scale = 0;
    u32 least_repair = 0;
    u64 nonnegative_cells = 0;
};

DeckResult evaluate_matrix(const std::vector<u32>& candidates,
                           const LiteralWeights& weights) {
    DeckResult result;
    bool first_common = true;
    for (std::size_t index = 0; index < candidates.size(); ++index) {
        const u32 repair = candidates[index];
        std::array<i128, SCALE_COUNT> ticks{};
        u32 sub = repair;
        for (;;) {
            const auto found = weights.by_failed.find(sub);
            if (found != weights.by_failed.end()) {
                for (int si = 0; si < SCALE_COUNT; ++si) {
                    ticks[si] += found->second[si];
                }
            }
            if (sub == 0) break;
            sub = (sub - 1) & repair;
        }

        bool active = true;
        i128 local_least = 0;
        int local_scale = 0;
        bool local_first = true;
        result.matrix_fnv.word(repair);
        for (int si = 0; si < SCALE_COUNT; ++si) {
            const int g = G0 + si;
            const i128 margin = 63 * ticks[si] - 4 * weights.denominator[si];
            result.matrix_fnv.word(static_cast<u64>(g));
            result.matrix_fnv.signed128(margin);
            if (margin >= 0) ++result.nonnegative_cells;
            if (margin < 0) active = false;
            if (local_first || margin < local_least) {
                local_first = false;
                local_least = margin;
                local_scale = g;
            }
        }
        if (!active) continue;
        for (std::size_t slot = 0; slot < BUDGETS.size(); ++slot) {
            if (index < BUDGETS[slot]) result.common[slot].push_back(repair);
        }
        if (first_common || local_least < result.least_margin) {
            first_common = false;
            result.least_margin = local_least;
            result.least_scale = local_scale;
            result.least_repair = repair;
        }
    }
    return result;
}

struct BodyResult {
    u64 bodies = 0;
    u64 failures = 0;
    u64 checks = 0;
    u64 max_checks = 0;
    u32 worst_body = 0;
    u32 worst_repair = 0;
    u32 first_failure = std::numeric_limits<u32>::max();
    std::vector<u32> failed_bodies;
};

std::vector<u32> enumerate_bodies() {
    std::vector<u32> bodies;
    bodies.reserve(BODY_UNIVERSE);
    u32 body = (u32{1} << 9) - 1;
    const u32 limit = u32{1} << 30;
    while (body != 0 && body < limit) {
        bodies.push_back(body);
        const u32 following = next_mask(body);
        if (following <= body) break;
        body = following;
    }
    require(bodies.size() == BODY_UNIVERSE, "C(30,9) census mismatch");
    return bodies;
}

BodyResult scan_bodies(const std::vector<u32>& deck,
                       const std::vector<u32>& bodies) {
    constexpr unsigned lanes = 16;
    std::array<BodyResult, lanes> local{};
    std::vector<std::thread> workers;
    for (unsigned lane = 0; lane < lanes; ++lane) {
        workers.emplace_back([&, lane]() {
            BodyResult report;
            const std::size_t begin = bodies.size() * lane / lanes;
            const std::size_t end = bodies.size() * (lane + 1) / lanes;
            for (std::size_t j = begin; j < end; ++j) {
                const u32 body = bodies[j];
                u64 used = 0;
                u32 witness = 0;
                for (u32 repair : deck) {
                    ++used;
                    if ((body & repair) == 0) {
                        witness = repair;
                        break;
                    }
                }
                ++report.bodies;
                report.checks += used;
                if (witness == 0) {
                    ++report.failures;
                    report.failed_bodies.push_back(body);
                    report.first_failure = std::min(report.first_failure, body);
                }
                if (used > report.max_checks ||
                    (used == report.max_checks && body < report.worst_body)) {
                    report.max_checks = used;
                    report.worst_body = body;
                    report.worst_repair = witness;
                }
            }
            local[lane] = std::move(report);
        });
    }
    for (std::thread& worker : workers) worker.join();

    BodyResult result;
    for (const BodyResult& report : local) {
        result.bodies += report.bodies;
        result.failures += report.failures;
        result.checks += report.checks;
        result.first_failure = std::min(result.first_failure, report.first_failure);
        result.failed_bodies.insert(result.failed_bodies.end(),
                                    report.failed_bodies.begin(),
                                    report.failed_bodies.end());
        if (report.max_checks > result.max_checks ||
            (report.max_checks == result.max_checks &&
             report.worst_body < result.worst_body)) {
            result.max_checks = report.max_checks;
            result.worst_body = report.worst_body;
            result.worst_repair = report.worst_repair;
        }
    }
    std::sort(result.failed_bodies.begin(), result.failed_bodies.end());
    return result;
}

std::vector<u32> build_targeted(const std::vector<u32>& base,
                                const std::vector<u32>& full,
                                const std::vector<u32>& hostiles,
                                std::vector<u32>& additions) {
    std::vector<u32> result = base;
    for (u32 body : hostiles) {
        bool covered = false;
        for (u32 repair : result) {
            if ((repair & body) == 0) {
                covered = true;
                break;
            }
        }
        if (covered) continue;
        const auto witness = std::find_if(
            full.begin(), full.end(),
            [body](u32 repair) { return (repair & body) == 0; });
        require(witness != full.end(), "full deck lacks a targeted witness");
        result.push_back(*witness);
        additions.push_back(*witness);
    }
    return result;
}

i128 literal_margin_for(u32 repair, int si, const LiteralWeights& weights) {
    i128 ticks = 0;
    u32 sub = repair;
    for (;;) {
        const auto found = weights.by_failed.find(sub);
        if (found != weights.by_failed.end()) ticks += found->second[si];
        if (sub == 0) break;
        sub = (sub - 1) & repair;
    }
    return 63 * ticks - 4 * weights.denominator[si];
}

}  // namespace

int main() {
    try {
        const auto primitive = primitive_control();
        const auto candidates = select_candidates();
        const auto cells = pool_cells();
        const auto weights = literal_weights(cells);
        const auto deck = evaluate_matrix(candidates, weights);

        Fnv64 candidate_fnv;
        Fnv64 scale_fnv;
        for (u32 repair : candidates) candidate_fnv.word(repair);
        for (int g = G0; g <= G1; ++g) scale_fnv.word(static_cast<u64>(g));

        std::cout << "LRC14_RATIO_5_6_LITERAL_JOINT_WALL_AUDIT\n"
                  << "METHOD literal_intervals_no_primitive_prefix\n"
                  << "PRIMITIVE_N 420 SAFE_TICKS " << primitive.safe_ticks
                  << " DENSITY 31/42 ARCS";
        for (const Interval& arc : primitive.arcs) {
            std::cout << " [" << arc.left << "," << arc.right << "]";
        }
        std::cout << "\nCENTERED_MIN " << primitive.centered_min
                  << " CENTERED_MAX " << primitive.centered_max
                  << " DELTA_C "
                  << primitive.centered_max - primitive.centered_min << "\n"
                  << "POOL_WALLS " << cells.size() + 1
                  << " UNIQUE_FAILURE_MASKS " << weights.by_failed.size() << "\n"
                  << "CANDIDATES " << candidates.size()
                  << " CANDIDATE_FNV " << hex64(candidate_fnv.value()) << "\n"
                  << "SCALES " << SCALE_COUNT << " RANGE " << G0 << ".." << G1
                  << " SCALE_FNV " << hex64(scale_fnv.value())
                  << " JOINT_COMPONENTS_FIRST " << weights.joint_components.front()
                  << " JOINT_COMPONENTS_LAST " << weights.joint_components.back() << "\n"
                  << "MATRIX_CELLS " << candidates.size() * SCALE_COUNT
                  << " NONNEGATIVE " << deck.nonnegative_cells
                  << " MATRIX_FNV " << hex64(deck.matrix_fnv.value()) << "\n"
                  << "LEAST_COMMON_MARGIN " << dec(deck.least_margin)
                  << " AT_SCALE " << deck.least_scale
                  << " REPAIR {" << labels(deck.least_repair) << "}\n";

        const auto bodies = enumerate_bodies();
        std::array<BodyResult, BUDGETS.size()> scans;
        for (std::size_t slot = 0; slot < BUDGETS.size(); ++slot) {
            require(deck.common[slot].size() == EXPECTED_COMMON[slot],
                    "common budget census mismatch");
            Fnv64 common_fnv;
            for (u32 repair : deck.common[slot]) common_fnv.word(repair);
            require(common_fnv.value() == EXPECTED_COMMON_FNV[slot],
                    "common budget FNV changed");
            scans[slot] = scan_bodies(deck.common[slot], bodies);
            require(scans[slot].failures == EXPECTED_FAILURES[slot],
                    "hostile budget census mismatch");
            require(scans[slot].checks == EXPECTED_CHECKS[slot] &&
                        scans[slot].max_checks == EXPECTED_MAX_CHECKS[slot],
                    "body-scan work ledger changed");
            std::cout << "BUDGET " << BUDGETS[slot]
                      << " COMMON " << deck.common[slot].size()
                      << " COMMON_FNV " << hex64(common_fnv.value())
                      << " BODY_COUNT " << scans[slot].bodies
                      << " FAILURES " << scans[slot].failures
                      << " CHECKS " << scans[slot].checks
                      << " MAX_CHECKS " << scans[slot].max_checks
                      << " FIRST_FAILURE ";
            if (scans[slot].failures == 0) std::cout << "NONE";
            else std::cout << "{" << labels(scans[slot].first_failure) << "}";
            std::cout << " WORST_BODY {" << labels(scans[slot].worst_body) << "} ";
            if (scans[slot].failures == 0) {
                std::cout << "WITNESS {" << labels(scans[slot].worst_repair) << "}";
            } else {
                std::cout << "WITNESS NONE";
            }
            std::cout << "\n";
        }
        require(scans.back().failures == 0, "full deck misses a body");

        const auto& hostiles = scans[1].failed_bodies;
        require(hostiles.size() == 6, "targeted hostile census changed");
        std::vector<u32> additions;
        const auto targeted =
            build_targeted(deck.common[1], deck.common[2], hostiles, additions);
        const BodyResult targeted_scan = scan_bodies(targeted, bodies);
        require(targeted_scan.failures == 0, "targeted deck misses a body");
        require(additions.size() == 5, "targeted addition census changed");
        require(std::equal(additions.begin(), additions.end(),
                           EXPECTED_ADDITIONS.begin()),
                "targeted additions changed");
        Fnv64 targeted_fnv;
        for (u32 repair : targeted) targeted_fnv.word(repair);
        i128 targeted_least_margin = 0;
        int targeted_least_scale = 0;
        u32 targeted_least_repair = 0;
        bool targeted_first = true;
        for (u32 repair : targeted) {
            for (int si = 0; si < SCALE_COUNT; ++si) {
                const i128 value = literal_margin_for(repair, si, weights);
                if (targeted_first || value < targeted_least_margin) {
                    targeted_first = false;
                    targeted_least_margin = value;
                    targeted_least_scale = G0 + si;
                    targeted_least_repair = repair;
                }
            }
        }

        std::cout << "TARGETED_FROM 8192 HOSTILE_COUNT " << hostiles.size()
                  << " HOSTILE_BODIES";
        for (u32 body : hostiles) std::cout << " {" << labels(body) << "}";
        std::cout << " ADDITION_COUNT " << additions.size() << " ADDITIONS";
        for (u32 repair : additions) {
            const auto candidate_position =
                std::find(candidates.begin(), candidates.end(), repair);
            const auto common_position =
                std::find(deck.common[2].begin(), deck.common[2].end(), repair);
            require(candidate_position != candidates.end(),
                    "addition absent from candidates");
            require(common_position != deck.common[2].end(),
                    "addition absent from common deck");
            std::cout << " {" << labels(repair) << "} MASK " << hex64(repair)
                      << " CANDIDATE_RANK "
                      << (candidate_position - candidates.begin() + 1)
                      << " COMMON_RANK "
                      << (common_position - deck.common[2].begin() + 1);
        }
        std::cout << " TARGETED_COUNT " << targeted.size()
                  << " TARGETED_FNV " << hex64(targeted_fnv.value())
                  << " FAILURES " << targeted_scan.failures
                  << " CHECKS " << targeted_scan.checks
                  << " MAX_CHECKS " << targeted_scan.max_checks
                  << " LEAST_LITERAL_TARGETED_MARGIN "
                  << dec(targeted_least_margin) << " AT_SCALE "
                  << targeted_least_scale << " REPAIR {"
                  << labels(targeted_least_repair) << "}\n";

        require(candidate_fnv.value() == UINT64_C(0xadf20f0ef1cadc1f),
                "candidate FNV changed");
        require(deck.nonnegative_cells == UINT64_C(987824),
                "nonnegative matrix census changed");
        require(deck.matrix_fnv.value() == UINT64_C(0xa039beab75bfbc2c),
                "literal matrix FNV changed");
        require(targeted_fnv.value() == UINT64_C(0x4e8a91621047de5c),
                "targeted FNV changed");
        require(deck.least_margin == static_cast<i128>(3716617212LL),
                "least literal margin changed");
        require(deck.least_scale == 77, "least-margin scale changed");
        require(deck.least_repair == UINT32_C(0x25220908),
                "least-margin repair changed");
        require(targeted_scan.checks == UINT64_C(486631847) &&
                    targeted_scan.max_checks == UINT64_C(3163),
                "targeted body-scan ledger changed");
        require(targeted_least_margin == static_cast<i128>(10908182376LL) &&
                    targeted_least_scale == 64 &&
                    targeted_least_repair == UINT32_C(0x10053221),
                "targeted literal minimum changed");
        std::cout << "VERDICT PASS_LITERAL_ALL_16384x70_AND_ALL_14307150_BODIES\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << "\n";
        return 1;
    }
}
