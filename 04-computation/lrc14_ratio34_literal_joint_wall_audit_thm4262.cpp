#include <algorithm>
#include <array>
#include <atomic>
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

// Clean-room audit of the proposed LRC(14) ratio-3:4 common deck.
//
// The analytic side is deliberately NOT evaluated by a primitive-prefix
// antiderivative.  At every scale g this program constructs the literal safe
// intervals for speeds 3g and 4g, intersects those intervals, and then sweeps
// that joint wall arrangement against every atom of the 30-speed pool wall
// arrangement.  All endpoints are integer ticks on an exact common
// denominator lcm(D,168g).

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
constexpr std::size_t KEEP = 8192;
constexpr u64 REPAIR_UNIVERSE = UINT64_C(5852925);
constexpr u64 BODY_UNIVERSE = UINT64_C(14307150);
constexpr int G0 = 97;
constexpr int G1 = 192;
constexpr int SCALE_COUNT = G1 - G0 + 1;

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
    // Independent top-k selection: unlike the candidate scout, do not sort
    // all C(30,8) records.
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
    return 7 * twice_phase >= D && 7 * twice_phase <= 13 * static_cast<i128>(D);
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

    std::vector<PoolCell> cells;
    cells.reserve(walls.size() - 1);
    for (std::size_t j = 0; j + 1 < walls.size(); ++j) {
        u32 failed = 0;
        for (int bit = 0; bit < 30; ++bit) {
            if (!pool_speed_safe(POOL[bit], walls[j], walls[j + 1])) {
                failed |= u32{1} << bit;
            }
        }
        cells.push_back({walls[j], walls[j + 1], failed});
    }
    return cells;
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
    const auto a = speed_safe_intervals(3 * g, denominator);
    const auto b = speed_safe_intervals(4 * g, denominator);
    std::vector<Interval> result;
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < a.size() && j < b.size()) {
        const i64 left = std::max(a[i].left, b[j].left);
        const i64 right = std::min(a[i].right, b[j].right);
        if (left < right) {
            if (!result.empty() && result.back().right == left) {
                result.back().right = right;
            } else {
                result.push_back({left, right});
            }
        }
        if (a[i].right < b[j].right) {
            ++i;
        } else if (b[j].right < a[i].right) {
            ++j;
        } else {
            ++i;
            ++j;
        }
    }
    require(!result.empty(), "empty literal joint safe set");
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
        const i64 denominator = exact_lcm(D, 168LL * g);
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
        // The literal 3:4 joint set has six components and length 3/4 in
        // one primitive cell; scaling repeats it g times.
        require(static_cast<i128>(4) * total == static_cast<i128>(3) * denominator,
                "literal joint total-length control failed");
    }
    return result;
}

struct DeckResult {
    std::vector<u32> common;
    std::vector<u32> truncated;
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
            const auto it = weights.by_failed.find(sub);
            if (it != weights.by_failed.end()) {
                for (int si = 0; si < SCALE_COUNT; ++si) {
                    ticks[si] += it->second[si];
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
        result.common.push_back(repair);
        if (index < 4096) result.truncated.push_back(repair);
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
    const unsigned lanes = 16;
    std::array<BodyResult, lanes> local{};
    std::vector<std::thread> workers;
    for (unsigned lane = 0; lane < lanes; ++lane) {
        workers.emplace_back([&, lane]() {
            BodyResult r;
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
                ++r.bodies;
                r.checks += used;
                if (witness == 0) {
                    ++r.failures;
                    r.failed_bodies.push_back(body);
                }
                if (used > r.max_checks ||
                    (used == r.max_checks && body < r.worst_body)) {
                    r.max_checks = used;
                    r.worst_body = body;
                    r.worst_repair = witness;
                }
            }
            local[lane] = std::move(r);
        });
    }
    for (auto& worker : workers) worker.join();

    BodyResult result;
    for (const BodyResult& r : local) {
        result.bodies += r.bodies;
        result.failures += r.failures;
        result.checks += r.checks;
        result.failed_bodies.insert(result.failed_bodies.end(),
                                    r.failed_bodies.begin(), r.failed_bodies.end());
        if (r.max_checks > result.max_checks ||
            (r.max_checks == result.max_checks && r.worst_body < result.worst_body)) {
            result.max_checks = r.max_checks;
            result.worst_body = r.worst_body;
            result.worst_repair = r.worst_repair;
        }
    }
    std::sort(result.failed_bodies.begin(), result.failed_bodies.end());
    return result;
}

}  // namespace

int main() {
    try {
        const auto candidates = select_candidates();
        const auto cells = pool_cells();
        const auto weights = literal_weights(cells);
        const auto deck = evaluate_matrix(candidates, weights);

        Fnv64 candidate_fnv;
        Fnv64 deck_fnv;
        Fnv64 scale_fnv;
        for (u32 repair : candidates) candidate_fnv.word(repair);
        for (u32 repair : deck.common) deck_fnv.word(repair);
        for (int g = G0; g <= G1; ++g) scale_fnv.word(static_cast<u64>(g));

        std::cout << "LRC14_RATIO_3_4_LITERAL_JOINT_WALL_AUDIT\n"
                  << "METHOD literal_intervals_no_primitive_prefix\n"
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
                  << "COMMON_FULL " << deck.common.size()
                  << " COMMON_FULL_FNV " << hex64(deck_fnv.value()) << "\n"
                  << "LEAST_COMMON_MARGIN " << dec(deck.least_margin)
                  << " AT_SCALE " << deck.least_scale
                  << " REPAIR {" << labels(deck.least_repair) << "}\n"
                  << "HOSTILE_FIRST4096_COMMON " << deck.truncated.size() << "\n";

        require(deck.common.size() == 3751, "common deck size mismatch");
        require(deck_fnv.value() == UINT64_C(0xd708c4fa029105e1),
                "common deck FNV mismatch");

        const auto bodies = enumerate_bodies();
        const auto hostile = scan_bodies(deck.truncated, bodies);
        std::cout << "HOSTILE_BODY_COUNT " << hostile.bodies
                  << " FAILURES " << hostile.failures
                  << " CHECKS " << hostile.checks
                  << " MAX_CHECKS " << hostile.max_checks << "\n";
        for (u32 body : hostile.failed_bodies) {
            std::cout << "HOSTILE_FAILURE {" << labels(body) << "}\n";
        }
        require(hostile.failures == 2, "first-4096 hostile changed");

        const auto full = scan_bodies(deck.common, bodies);
        std::cout << "FULL_BODY_COUNT " << full.bodies
                  << " FAILURES " << full.failures
                  << " CHECKS " << full.checks
                  << " MAX_CHECKS " << full.max_checks << "\n"
                  << "WORST_BODY {" << labels(full.worst_body)
                  << "} WITNESS {" << labels(full.worst_repair) << "}\n";
        require(full.failures == 0, "common deck misses a 9-body");
        require(full.bodies == BODY_UNIVERSE, "incomplete full body scan");
        std::cout << "VERDICT PASS_LITERAL_ALL_8192x96_AND_ALL_14307150_BODIES\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << "\n";
        return 1;
    }
}
