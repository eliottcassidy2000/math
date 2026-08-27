#include "lrc14_two_three_outsider_ray_thm4256/independent_common.hpp"

#include <array>

using namespace ray_audit;

namespace {

constexpr i64 U = 4;
constexpr i64 V = 5;
constexpr i64 N45 = 14 * U * V;
constexpr std::size_t CANDIDATES = 65536;
constexpr std::array<std::size_t, 5> BUDGETS = {
    4096, 8192, 16384, 32768, 65536,
};
constexpr std::array<u64, 5> EXPECTED_COMMON = {
    1245, 2417, 4942, 9805, 19763,
};
constexpr std::array<u64, 5> EXPECTED_FAILURES = {
    1377, 219, 30, 2, 0,
};

std::vector<i64> full_scales() {
    std::vector<i64> result;
    for (i64 scale = 73; scale <= 153; ++scale) result.push_back(scale);
    ensure(result.size() == 81, "full scale census mismatch");
    return result;
}

struct Primitive {
    std::vector<std::pair<i64, i64>> arcs;
    i64 safe_ticks = 0;
};

bool primitive_safe_midpoint(i64 speed, i64 left, i64 right) {
    const i64 residue = (speed * (left + right)) % (2 * N45);
    return 7 * residue >= N45 && 7 * residue <= 13 * N45;
}

Primitive build_primitive() {
    std::vector<i64> walls = {0, N45};
    for (i64 speed : {U, V}) {
        const i64 factor = N45 / (14 * speed);
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
    ensure(result.safe_ticks == 208, "primitive safe length changed");
    ensure(result.arcs.size() == 8, "primitive arc count changed");
    return result;
}

i128 primitive_prefix(const Primitive& primitive, i64 remainder) {
    ensure(0 <= remainder && remainder < D, "bad primitive remainder");
    const i128 y = static_cast<i128>(N45) * remainder;
    i128 answer = 0;
    for (const auto& [left_tick, right_tick] : primitive.arcs) {
        const i128 left = static_cast<i128>(left_tick) * D;
        const i128 right = static_cast<i128>(right_tick) * D;
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
    return 63 * mass_numerator - 4 * N45 * D * scale;
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
              [](const OrderedRepair& left, const OrderedRepair& right) {
                  return left.key != right.key ? left.key < right.key
                                               : left.mask < right.mask;
              });
    std::vector<u32> result;
    result.reserve(CANDIDATES);
    for (std::size_t i = 0; i < CANDIDATES; ++i) {
        result.push_back(universe[i].mask);
    }
    return result;
}

std::vector<u32> uncovered_bodies(const std::vector<u32>& deck,
                                  const std::vector<u32>& bodies) {
    const unsigned available = std::thread::hardware_concurrency();
    const unsigned threads = std::max(1u, std::min(16u, available ? available : 1u));
    std::vector<std::vector<u32>> local(threads);
    std::vector<std::thread> workers;
    for (unsigned lane = 0; lane < threads; ++lane) {
        workers.emplace_back([&, lane]() {
            const std::size_t begin = bodies.size() * lane / threads;
            const std::size_t end = bodies.size() * (lane + 1) / threads;
            for (std::size_t index = begin; index < end; ++index) {
                const u32 body = bodies[index];
                bool covered = false;
                for (u32 repair : deck) {
                    if ((repair & body) == 0) {
                        covered = true;
                        break;
                    }
                }
                if (!covered) local[lane].push_back(body);
            }
        });
    }
    for (std::thread& worker : workers) worker.join();
    std::vector<u32> result;
    for (const auto& lane : local) result.insert(result.end(), lane.begin(), lane.end());
    std::sort(result.begin(), result.end());
    return result;
}

}  // namespace

int main() {
    try {
        const Primitive primitive = build_primitive();
        const auto scales = full_scales();
        const auto cells = build_pool_arrangement();
        const auto repairs = first_candidates();

        std::vector<Geometry> geometry;
        geometry.reserve(repairs.size());
        for (u32 repair : repairs) {
            geometry.push_back(direct_geometry(repair, cells));
        }

        std::array<std::vector<u32>, BUDGETS.size()> common;
        i128 least_margin = 0;
        i64 least_scale = 0;
        u32 least_repair = 0;
        bool first_common = true;
        for (std::size_t index = 0; index < geometry.size(); ++index) {
            bool active = true;
            i128 local_least = 0;
            i64 local_scale = 0;
            bool local_first = true;
            for (i64 scale : scales) {
                const i128 value = margin(primitive, geometry[index], scale);
                if (local_first || value < local_least) {
                    local_first = false;
                    local_least = value;
                    local_scale = scale;
                }
                if (value < 0) active = false;
            }
            if (!active) continue;
            for (std::size_t slot = 0; slot < BUDGETS.size(); ++slot) {
                if (index < BUDGETS[slot]) {
                    common[slot].push_back(geometry[index].repair);
                }
            }
            if (first_common || local_least < least_margin) {
                first_common = false;
                least_margin = local_least;
                least_scale = local_scale;
                least_repair = geometry[index].repair;
            }
        }

        Fnv64 candidate_fnv;
        Fnv64 common_fnv;
        Fnv64 scale_fnv;
        for (u32 repair : repairs) candidate_fnv.word(repair);
        for (u32 repair : common.back()) common_fnv.word(repair);
        for (i64 scale : scales) scale_fnv.word(static_cast<u64>(scale));

        std::cout << "LRC14_RATIO_4_5_COMMON_DECK_BRIDGE\n"
                  << "PRIMITIVE_N " << N45 << " SAFE_TICKS "
                  << primitive.safe_ticks << " ARCS";
        for (const auto& [left, right] : primitive.arcs) {
            std::cout << " [" << left << "," << right << "]";
        }
        std::cout << "\nCANDIDATES " << repairs.size() << " CANDIDATE_FNV "
                  << hex64(candidate_fnv.value()) << "\nCOMMON_FULL "
                  << common.back().size() << " COMMON_FULL_FNV "
                  << hex64(common_fnv.value()) << " FULL_SCALES "
                  << scales.size() << " RANGE 73..153 SCALE_FNV "
                  << hex64(scale_fnv.value()) << "\nLEAST_COMMON_MARGIN "
                  << decimal(least_margin) << " AT_SCALE " << least_scale
                  << " REPAIR {" << labels(least_repair) << "}\n";

        const auto bodies = all_bodies();
        for (std::size_t slot = 0; slot < BUDGETS.size(); ++slot) {
            ensure(common[slot].size() == EXPECTED_COMMON[slot],
                   "common-deck budget census changed");
            const BodyScan scan = scan_bodies(common[slot], bodies);
            ensure(scan.failures == EXPECTED_FAILURES[slot],
                   "candidate-budget hostile changed");
            std::cout << "BUDGET " << BUDGETS[slot] << " COMMON "
                      << common[slot].size() << " BODY_COUNT " << scan.bodies
                      << " FAILURES " << scan.failures << " CHECKS "
                      << scan.checks << " MAX_CHECKS " << scan.max_checks
                      << " FIRST_FAILURE ";
            if (scan.failures == 0) {
                std::cout << "NONE";
            } else {
                std::cout << "{" << labels(scan.first_failure) << "}";
            }
            std::cout << " WORST_BODY {" << labels(scan.worst_body) << "} ";
            if (scan.failures == 0) {
                std::cout << "WITNESS {" << labels(scan.worst_repair) << "}";
            } else {
                std::cout << "WITNESS NONE";
            }
            std::cout << "\n";
        }

        const std::vector<u32> hostile_bodies = uncovered_bodies(common[3], bodies);
        ensure(hostile_bodies.size() == 2, "32K hostile-body census changed");
        std::vector<u32> targeted = common[3];
        std::vector<u32> targeted_additions;
        for (u32 body : hostile_bodies) {
            bool now_covered = false;
            for (u32 repair : targeted) {
                if ((repair & body) == 0) {
                    now_covered = true;
                    break;
                }
            }
            if (now_covered) continue;
            const auto witness = std::find_if(
                common.back().begin(), common.back().end(),
                [body](u32 repair) { return (repair & body) == 0; });
            ensure(witness != common.back().end(), "full deck lacks targeted repair");
            targeted.push_back(*witness);
            targeted_additions.push_back(*witness);
        }
        const BodyScan targeted_scan = scan_bodies(targeted, bodies);
        ensure(targeted_scan.failures == 0, "targeted deck misses a body");
        Fnv64 targeted_fnv;
        for (u32 repair : targeted) targeted_fnv.word(repair);
        std::cout << "TARGETED_FROM_32768 HOSTILE_BODIES";
        for (u32 body : hostile_bodies) std::cout << " {" << labels(body) << "}";
        std::cout << " ADDITIONS";
        for (u32 repair : targeted_additions) {
            const auto candidate_position = std::find(repairs.begin(), repairs.end(), repair);
            const auto common_position =
                std::find(common.back().begin(), common.back().end(), repair);
            ensure(candidate_position != repairs.end(), "target repair missing from candidates");
            ensure(common_position != common.back().end(), "target repair missing from full common deck");
            std::cout << " {" << labels(repair) << "} MASK "
                      << hex64(repair) << " CANDIDATE_RANK "
                      << (candidate_position - repairs.begin() + 1)
                      << " COMMON_RANK "
                      << (common_position - common.back().begin() + 1);
        }
        std::cout << " TARGETED_COUNT " << targeted.size() << " TARGETED_FNV "
                  << hex64(targeted_fnv.value()) << " FAILURES "
                  << targeted_scan.failures << " CHECKS " << targeted_scan.checks
                  << " MAX_CHECKS " << targeted_scan.max_checks << "\n";
        std::cout << "VERDICT COMMON_DECK_CLOSES_EVERY_BODY_ALL_SCALES_73_153\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << "\n";
        return 1;
    }
}
