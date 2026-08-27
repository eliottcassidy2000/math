#include "lrc14_two_three_outsider_ray_thm4256/independent_common.hpp"

#include <set>

using namespace ray_audit;

namespace {

constexpr i64 U = 3;
constexpr i64 V = 4;
constexpr i64 N34 = 14 * U * V;
constexpr std::size_t CANDIDATES = 8192;

std::vector<i64> residual_scales() {
    std::vector<i64> result;
    for (i64 scale = 97; scale <= 166; ++scale) result.push_back(scale);
    result.push_back(172);
    result.push_back(180);
    ensure(result.size() == 72, "residual scale census mismatch");
    return result;
}
std::vector<i64> full_scales() {
    std::vector<i64> result;
    for (i64 scale = 97; scale <= 192; ++scale) result.push_back(scale);
    ensure(result.size() == 96, "full scale census mismatch");
    return result;
}

struct Primitive {
    std::vector<std::pair<i64, i64>> arcs;
    i64 safe_ticks = 0;
};

bool primitive_safe_midpoint(i64 speed, i64 left, i64 right) {
    const i64 residue = (speed * (left + right)) % (2 * N34);
    return 7 * residue >= N34 && 7 * residue <= 13 * N34;
}

Primitive build_primitive() {
    std::vector<i64> walls = {0, N34};
    for (i64 speed : {U, V}) {
        const i64 factor = N34 / (14 * speed);
        for (i64 tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14 * tooth + 1) * factor);
            walls.push_back((14 * tooth + 13) * factor);
        }
    }
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());

    Primitive result;
    for (std::size_t i = 0; i + 1 < walls.size(); ++i) {
        if (!primitive_safe_midpoint(U, walls[i], walls[i + 1]) ||
            !primitive_safe_midpoint(V, walls[i], walls[i + 1])) {
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

i128 primitive_prefix(const Primitive& primitive, i64 remainder) {
    ensure(0 <= remainder && remainder < D, "bad primitive remainder");
    const i128 y = static_cast<i128>(N34) * remainder;
    i128 answer = 0;
    for (const auto& [a, b] : primitive.arcs) {
        const i128 left = static_cast<i128>(a) * D;
        const i128 right = static_cast<i128>(b) * D;
        if (y > left) answer += std::min(y, right) - left;
    }
    return answer;
}

i128 primitive_integral(const Primitive& primitive, i128 z) {
    ensure(z >= 0, "negative primitive coordinate");
    const i128 whole = z / D;
    const i64 remainder = static_cast<i64>(z % D);
    return whole * primitive.safe_ticks * D +
           primitive_prefix(primitive, remainder);
}

i128 margin(const Primitive& primitive, const Geometry& geometry, i64 scale) {
    i128 mass_numerator = 0;
    for (const auto& [left, right] : geometry.components) {
        mass_numerator +=
            primitive_integral(primitive, static_cast<i128>(scale) * right) -
            primitive_integral(primitive, static_cast<i128>(scale) * left);
    }
    return 63 * mass_numerator - 4 * N34 * D * scale;
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
        universe.push_back(
            {splitmix64(static_cast<u64>(mask) ^ ORDER_SEED), mask});
        const u32 next = next_mask(mask);
        if (next <= mask) break;
        mask = next;
    }
    ensure(universe.size() == REPAIR_COUNT, "repair census mismatch");
    std::sort(universe.begin(), universe.end(),
              [](const OrderedRepair& a, const OrderedRepair& b) {
                  return a.key != b.key ? a.key < b.key : a.mask < b.mask;
              });
    std::vector<u32> result;
    for (std::size_t i = 0; i < CANDIDATES; ++i) {
        result.push_back(universe[i].mask);
    }
    return result;
}

}  // namespace

int main() {
    try {
        const Primitive primitive = build_primitive();
        const auto cells = build_pool_arrangement();
        const auto repairs = first_candidates();
        std::vector<Geometry> geometry;
        geometry.reserve(repairs.size());
        for (u32 repair : repairs) {
            geometry.push_back(direct_geometry(repair, cells));
        }

        std::vector<u32> common_residual;
        std::vector<u32> common_full;
        std::vector<u32> truncated_full;
        const auto residual = residual_scales();
        const auto full = full_scales();
        i128 least_margin = 0;
        i64 least_scale = 0;
        u32 least_repair = 0;
        bool first = true;
        for (std::size_t i = 0; i < geometry.size(); ++i) {
            bool active_residual = true;
            bool active_full = true;
            i128 local_least = 0;
            i64 local_scale = 0;
            bool local_first = true;
            for (i64 scale : full) {
                const i128 value = margin(primitive, geometry[i], scale);
                if (local_first || value < local_least) {
                    local_first = false;
                    local_least = value;
                    local_scale = scale;
                }
                if (value < 0) active_full = false;
            }
            for (i64 scale : residual) {
                if (margin(primitive, geometry[i], scale) < 0) {
                    active_residual = false;
                    break;
                }
            }
            if (active_residual) common_residual.push_back(geometry[i].repair);
            if (!active_full) continue;
            common_full.push_back(geometry[i].repair);
            if (i < 4096) truncated_full.push_back(geometry[i].repair);
            if (first || local_least < least_margin) {
                first = false;
                least_margin = local_least;
                least_scale = local_scale;
                least_repair = geometry[i].repair;
            }
        }

        Fnv64 candidate_fnv;
        Fnv64 residual_fnv;
        Fnv64 full_fnv;
        Fnv64 scale_fnv;
        for (u32 repair : repairs) candidate_fnv.word(repair);
        for (u32 repair : common_residual) residual_fnv.word(repair);
        for (u32 repair : common_full) full_fnv.word(repair);
        for (i64 scale : full) scale_fnv.word(static_cast<u64>(scale));

        std::cout << "LRC14_RATIO_3_4_COMMON_DECK_SCOUT\n"
                  << "PRIMITIVE_N " << N34 << " SAFE_TICKS "
                  << primitive.safe_ticks << " ARCS";
        for (const auto& [left, right] : primitive.arcs) {
            std::cout << " [" << left << "," << right << "]";
        }
        std::cout << "\nCANDIDATES " << repairs.size() << " CANDIDATE_FNV "
                  << hex64(candidate_fnv.value()) << "\nCOMMON_RESIDUAL "
                  << common_residual.size() << " COMMON_RESIDUAL_FNV "
                  << hex64(residual_fnv.value()) << " RESIDUAL_SCALES "
                  << residual.size() << "\nCOMMON_FULL " << common_full.size()
                  << " COMMON_FULL_FNV " << hex64(full_fnv.value())
                  << " FULL_SCALES " << full.size() << " RANGE 97..192"
                  << " SCALE_FNV " << hex64(scale_fnv.value()) << "\n";
        std::cout << "LEAST_COMMON_MARGIN " << decimal(least_margin)
                  << " AT_SCALE " << least_scale << " REPAIR {"
                  << labels(least_repair) << "}\n";

        const auto bodies = all_bodies();
        const BodyScan hostile = scan_bodies(truncated_full, bodies);
        std::cout << "HOSTILE_FIRST4096_COMMON " << truncated_full.size()
                  << " BODY_COUNT " << hostile.bodies << " FAILURES "
                  << hostile.failures << " FIRST_FAILURE {"
                  << labels(hostile.first_failure) << "}\n";
        ensure(hostile.failures == 2, "candidate-truncation hostile changed");

        const BodyScan scan = scan_bodies(common_full, bodies);
        std::cout << "FULL_BODY_COUNT " << scan.bodies << " FAILURES "
                  << scan.failures << " CHECKS " << scan.checks
                  << " MAX_CHECKS " << scan.max_checks << "\n";
        ensure(scan.failures == 0, "full common deck misses a body");
        std::cout << "WORST_BODY {" << labels(scan.worst_body)
                  << "} WITNESS {" << labels(scan.worst_repair)
                  << "}\nVERDICT COMMON_DECK_CLOSES_EVERY_BODY_ALL_SCALES_97_192\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << "\n";
        return 1;
    }
}
