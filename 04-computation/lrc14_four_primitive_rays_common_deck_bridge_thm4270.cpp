#include "lrc14_two_three_outsider_ray_thm4256/independent_common.hpp"

#include <set>

using namespace ray_audit;

namespace {

constexpr std::size_t CANDIDATES = 8192;

struct Ratio {
    i64 u;
    i64 v;
    std::size_t expected_common;
    u64 expected_fnv;
};

struct Primitive {
    i64 period = 0;
    i64 safe_ticks = 0;
    std::vector<std::pair<i64, i64>> arcs;
};

bool primitive_safe_midpoint(
    i64 speed, i64 left, i64 right, i64 period) {
    const i64 residue = static_cast<i64>(
        (static_cast<i128>(speed) * (left + right)) % (2 * period));
    return static_cast<i128>(7) * residue >= period &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * period;
}

Primitive build_primitive(i64 u, i64 v) {
    ensure(1 <= u && u < v && std::gcd(u, v) == 1, "ratio not primitive");
    Primitive result;
    result.period = 14 * std::lcm(u, v);
    std::vector<i64> walls = {0, result.period};
    for (i64 speed : {u, v}) {
        const i64 quantum = result.period / (14 * speed);
        for (i64 tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14 * tooth + 1) * quantum);
            walls.push_back((14 * tooth + 13) * quantum);
        }
    }
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    for (std::size_t i = 0; i + 1 < walls.size(); ++i) {
        if (!primitive_safe_midpoint(u, walls[i], walls[i + 1], result.period) ||
            !primitive_safe_midpoint(v, walls[i], walls[i + 1], result.period)) {
            continue;
        }
        if (!result.arcs.empty() && result.arcs.back().second == walls[i]) {
            result.arcs.back().second = walls[i + 1];
        } else {
            result.arcs.push_back({walls[i], walls[i + 1]});
        }
        result.safe_ticks += walls[i + 1] - walls[i];
    }
    ensure(result.safe_ticks > 0, "empty primitive safe set");
    return result;
}

i128 prefix(const Primitive& primitive, i64 remainder) {
    ensure(0 <= remainder && remainder < D, "bad remainder");
    const i128 y = static_cast<i128>(primitive.period) * remainder;
    i128 answer = 0;
    for (const auto& [a, b] : primitive.arcs) {
        const i128 left = static_cast<i128>(a) * D;
        const i128 right = static_cast<i128>(b) * D;
        if (y > left) answer += std::min(y, right) - left;
    }
    return answer;
}

i128 integral(const Primitive& primitive, i128 z) {
    ensure(z >= 0, "negative coordinate");
    const i128 whole = z / D;
    const i64 remainder = static_cast<i64>(z % D);
    return whole * primitive.safe_ticks * D + prefix(primitive, remainder);
}

i128 margin(
    const Primitive& primitive, const Geometry& geometry, i64 scale) {
    i128 numerator = 0;
    for (const auto& [left, right] : geometry.components) {
        numerator += integral(primitive, static_cast<i128>(scale) * right) -
                     integral(primitive, static_cast<i128>(scale) * left);
    }
    return 63 * numerator - 4 * primitive.period * D * scale;
}

struct OrderedRepair {
    u64 key;
    u32 mask;
};

std::vector<u32> first_candidates() {
    std::vector<OrderedRepair> universe;
    universe.reserve(REPAIR_COUNT);
    u32 mask = (u32{1} << 8) - 1;
    const u32 limit = u32{1} << 30;
    while (mask != 0 && mask < limit) {
        universe.push_back({splitmix64(static_cast<u64>(mask) ^ ORDER_SEED), mask});
        const u32 following = next_mask(mask);
        if (following <= mask) break;
        mask = following;
    }
    ensure(universe.size() == REPAIR_COUNT, "repair census mismatch");
    std::sort(universe.begin(), universe.end(),
              [](const OrderedRepair& a, const OrderedRepair& b) {
                  return a.key != b.key ? a.key < b.key : a.mask < b.mask;
              });
    std::vector<u32> result;
    result.reserve(CANDIDATES);
    for (std::size_t i = 0; i < CANDIDATES; ++i) result.push_back(universe[i].mask);
    return result;
}

std::vector<i64> bridge_scales(const Ratio& ratio) {
    const i64 first = 290 / ratio.u + 1;
    const i64 tail = (770 + ratio.v - 1) / ratio.v;
    ensure(first < tail, "empty bridge");
    std::vector<i64> result;
    for (i64 scale = first; scale < tail; ++scale) result.push_back(scale);
    return result;
}

void audit_ratio(
    const Ratio& ratio,
    const std::vector<Geometry>& geometries,
    const std::vector<u32>& bodies) {
    const Primitive primitive = build_primitive(ratio.u, ratio.v);
    const auto scales = bridge_scales(ratio);
    std::vector<u32> common;
    i128 least = 0;
    i64 least_scale = 0;
    u32 least_repair = 0;
    bool first = true;
    for (const Geometry& geometry : geometries) {
        bool active = true;
        i128 local_least = 0;
        i64 local_scale = 0;
        bool local_first = true;
        for (i64 scale : scales) {
            const i128 value = margin(primitive, geometry, scale);
            if (local_first || value < local_least) {
                local_first = false;
                local_least = value;
                local_scale = scale;
            }
            if (value < 0) {
                active = false;
                break;
            }
        }
        if (!active) continue;
        common.push_back(geometry.repair);
        if (first || local_least < least) {
            first = false;
            least = local_least;
            least_scale = local_scale;
            least_repair = geometry.repair;
        }
    }
    Fnv64 fnv;
    for (u32 repair : common) fnv.word(repair);
    const BodyScan scan = scan_bodies(common, bodies);
    std::cout << "RATIO " << ratio.u << ':' << ratio.v
              << " PERIOD " << primitive.period
              << " SAFE_TICKS " << primitive.safe_ticks
              << " ARCS " << primitive.arcs.size()
              << " RANGE " << scales.front() << ".." << scales.back()
              << " SCALES " << scales.size()
              << " COMMON " << common.size()
              << " COMMON_FNV " << hex64(fnv.value())
              << " LEAST_MARGIN " << decimal(least)
              << " AT " << least_scale
              << " REPAIR {" << labels(least_repair) << "}"
              << " BODY_FAILURES " << scan.failures
              << " CHECKS " << scan.checks
              << " MAX_CHECKS " << scan.max_checks;
    if (scan.failures != 0) {
        std::cout << " FIRST_FAILURE {" << labels(scan.first_failure) << "}";
    }
    std::cout << '\n';
    ensure(common.size() == ratio.expected_common, "common deck size mismatch");
    ensure(fnv.value() == ratio.expected_fnv, "common deck FNV mismatch");
    ensure(least > 0, "common deck has a nonpositive bridge margin");
    ensure(scan.failures == 0, "common deck misses a nine-body");
    ensure(scan.bodies == BODY_COUNT, "incomplete nine-body scan");
}

}  // namespace

int main() {
    try {
        const std::vector<Ratio> ratios = {
            {3, 5, 4178, UINT64_C(0x7a9034824db27f98)},
            {7, 8, 4011, UINT64_C(0x7ecafda9de695f6d)},
            {8, 9, 3046, UINT64_C(0x35a89d2659e92442)},
            {11, 12, 4073, UINT64_C(0x5c48ba1405aa1453)}};
        const auto cells = build_pool_arrangement();
        const auto candidates = first_candidates();
        std::vector<Geometry> geometries;
        geometries.reserve(candidates.size());
        for (u32 repair : candidates) geometries.push_back(direct_geometry(repair, cells));
        const auto bodies = all_bodies();
        Fnv64 candidate_fnv;
        for (u32 repair : candidates) candidate_fnv.word(repair);
        ensure(candidate_fnv.value() == UINT64_C(0x60148ca1fc61dbcb),
               "candidate deck FNV mismatch");
        std::cout << "LRC14_FOUR_PRIMITIVE_RAYS_COMMON_DECK_BRIDGE_THM4270\n"
                  << "CANDIDATES " << candidates.size()
                  << " CANDIDATE_FNV " << hex64(candidate_fnv.value())
                  << " BODIES " << bodies.size() << '\n';
        for (const Ratio& ratio : ratios) audit_ratio(ratio, geometries, bodies);
        std::cout << "SCOPE fixed_pool_nine_bodies strict_above_pool finite_bridge_only\n"
                  << "VERDICT PASS_PREFIX_INTEGRAL_FOUR_RAYS_ALL_8192_MATRICES_"
                     "AND_ALL_14307150_BODIES_EACH\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
