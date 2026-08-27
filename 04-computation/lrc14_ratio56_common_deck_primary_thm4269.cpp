#include "independent_common.hpp"

#include <array>
#include <fstream>

using namespace ray_audit;

namespace {

constexpr i64 U = 5;
constexpr i64 V = 6;
constexpr i64 N56 = 14 * U * V;
constexpr std::size_t CANDIDATES = 16384;
constexpr std::array<std::size_t, 3> BUDGETS = {
    4096, 8192, 16384,
};
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

std::vector<i64> bridge_scales() {
    std::vector<i64> result;
    // Strictly above the fixed pool means 5g>290, hence g>=59.
    // THM-4231 takes over once 6g>=770, hence g>=129.
    for (i64 scale = 59; scale <= 128; ++scale) result.push_back(scale);
    ensure(result.size() == 70, "bridge scale census mismatch");
    return result;
}

struct Primitive {
    std::vector<std::pair<i64, i64>> arcs;
    i64 safe_ticks = 0;
    i64 centered_min = 0;
    i64 centered_max = 0;
};

bool primitive_safe_midpoint(i64 speed, i64 left, i64 right) {
    const i64 residue = (speed * (left + right)) % (2 * N56);
    return 7 * residue >= N56 && 7 * residue <= 13 * N56;
}

i64 prefix_ticks(const std::vector<std::pair<i64, i64>>& arcs, i64 point) {
    i64 answer = 0;
    for (const auto& [left, right] : arcs) {
        answer += std::max<i64>(0, std::min(point, right) - left);
    }
    return answer;
}

Primitive build_primitive() {
    std::vector<i64> walls = {0, N56};
    for (i64 speed : {U, V}) {
        const i64 factor = N56 / (14 * speed);
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
    ensure(result.safe_ticks == 310, "primitive safe length changed");

    bool first = true;
    for (i64 wall : walls) {
        const i64 centered =
            N56 * prefix_ticks(result.arcs, wall) - result.safe_ticks * wall;
        if (first) {
            first = false;
            result.centered_min = result.centered_max = centered;
        } else {
            result.centered_min = std::min(result.centered_min, centered);
            result.centered_max = std::max(result.centered_max, centered);
        }
    }
    ensure(result.centered_max - result.centered_min == 9260,
           "primitive centered range changed");
    return result;
}

i128 primitive_prefix(const Primitive& primitive, i64 remainder) {
    ensure(0 <= remainder && remainder < D, "bad primitive remainder");
    const i128 y = static_cast<i128>(N56) * remainder;
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
    return 63 * mass_numerator - 4 * N56 * D * scale;
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

std::vector<u32> failed_bodies(const std::vector<u32>& deck,
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

void write_mask_ledger(const std::string& path, const std::vector<u32>& deck) {
    std::ofstream out(path);
    ensure(static_cast<bool>(out), "cannot open mask ledger");
    for (u32 repair : deck) {
        out << hex64(repair) << "," << labels(repair) << "\n";
    }
    ensure(static_cast<bool>(out), "mask ledger write failed");
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Primitive primitive = build_primitive();
        const auto scales = bridge_scales();
        const auto cells = build_pool_arrangement();
        const auto repairs = first_candidates();

        std::vector<Geometry> geometry;
        geometry.reserve(repairs.size());
        for (u32 repair : repairs) geometry.push_back(direct_geometry(repair, cells));

        std::array<std::vector<u32>, BUDGETS.size()> common;
        i128 least_margin = 0;
        i64 least_scale = 0;
        u32 least_repair = 0;
        bool first_common = true;
        i128 least_literal_margin = 0;
        i64 least_literal_scale = 0;
        u32 least_literal_repair = 0;
        bool first_literal_common = true;
        u64 nonnegative_cells = 0;
        Fnv64 matrix_fnv;
        Fnv64 literal_matrix_fnv;
        for (std::size_t index = 0; index < geometry.size(); ++index) {
            bool active = true;
            i128 local_least = 0;
            i64 local_scale = 0;
            bool local_first = true;
            i128 local_literal_least = 0;
            i64 local_literal_scale = 0;
            bool local_literal_first = true;
            matrix_fnv.word(geometry[index].repair);
            literal_matrix_fnv.word(geometry[index].repair);
            for (i64 scale : scales) {
                const i128 value = margin(primitive, geometry[index], scale);
                const i64 denominator = exact_lcm(D, N56 * scale);
                const i128 conversion =
                    static_cast<i128>(N56) * D * scale / denominator;
                ensure(value % conversion == 0,
                       "literal-margin conversion is nonintegral");
                const i128 literal_value = value / conversion;
                matrix_fnv.word(static_cast<u64>(scale));
                matrix_fnv.word(static_cast<u64>(value));
                matrix_fnv.word(static_cast<u64>(static_cast<__uint128_t>(value) >> 64));
                literal_matrix_fnv.word(static_cast<u64>(scale));
                literal_matrix_fnv.word(static_cast<u64>(literal_value));
                literal_matrix_fnv.word(static_cast<u64>(
                    static_cast<__uint128_t>(literal_value) >> 64));
                if (value >= 0) ++nonnegative_cells;
                if (local_first || value < local_least) {
                    local_first = false;
                    local_least = value;
                    local_scale = scale;
                }
                if (local_literal_first || literal_value < local_literal_least) {
                    local_literal_first = false;
                    local_literal_least = literal_value;
                    local_literal_scale = scale;
                }
                if (value < 0) active = false;
            }
            if (!active) continue;
            for (std::size_t slot = 0; slot < BUDGETS.size(); ++slot) {
                if (index < BUDGETS[slot]) common[slot].push_back(geometry[index].repair);
            }
            if (first_common || local_least < least_margin) {
                first_common = false;
                least_margin = local_least;
                least_scale = local_scale;
                least_repair = geometry[index].repair;
            }
            if (first_literal_common ||
                local_literal_least < least_literal_margin) {
                first_literal_common = false;
                least_literal_margin = local_literal_least;
                least_literal_scale = local_literal_scale;
                least_literal_repair = geometry[index].repair;
            }
        }

        Fnv64 candidate_fnv;
        Fnv64 common_fnv;
        Fnv64 scale_fnv;
        for (u32 repair : repairs) candidate_fnv.word(repair);
        for (u32 repair : common.back()) common_fnv.word(repair);
        for (i64 scale : scales) scale_fnv.word(static_cast<u64>(scale));
        ensure(candidate_fnv.value() == UINT64_C(0xadf20f0ef1cadc1f),
               "candidate FNV changed");
        ensure(common_fnv.value() == UINT64_C(0x9749384b49c2bba0),
               "full common FNV changed");
        ensure(scale_fnv.value() == UINT64_C(0x4146f5007ba813de),
               "scale FNV changed");
        ensure(nonnegative_cells == UINT64_C(987824),
               "nonnegative matrix census changed");
        ensure(matrix_fnv.value() == UINT64_C(0xee5331e52b0387bc),
               "primary matrix FNV changed");
        ensure(literal_matrix_fnv.value() == UINT64_C(0xa039beab75bfbc2c),
               "literal-coordinate matrix FNV changed");
        ensure(least_margin == static_cast<i128>(36651492783360LL) &&
                   least_scale == 64 && least_repair == UINT32_C(0x10053221),
               "least primary margin changed");
        ensure(least_literal_margin == static_cast<i128>(3716617212LL) &&
                   least_literal_scale == 77 &&
                   least_literal_repair == UINT32_C(0x25220908),
               "least literal common margin changed");

        std::cout << "LRC14_RATIO_5_6_COMMON_DECK_PRIMARY\n"
                  << "PRIMITIVE_N " << N56 << " SAFE_TICKS "
                  << primitive.safe_ticks << " DENSITY 31/42 ARCS";
        for (const auto& [left, right] : primitive.arcs) {
            std::cout << " [" << left << "," << right << "]";
        }
        std::cout << "\nCENTERED_MIN " << primitive.centered_min
                  << " CENTERED_MAX " << primitive.centered_max
                  << " DELTA_C " << primitive.centered_max - primitive.centered_min
                  << "\nCANDIDATES " << repairs.size() << " CANDIDATE_FNV "
                  << hex64(candidate_fnv.value()) << "\nCOMMON_FULL "
                  << common.back().size() << " COMMON_FULL_FNV "
                  << hex64(common_fnv.value()) << " BRIDGE_SCALES "
                  << scales.size() << " RANGE 59..128 SCALE_FNV "
                  << hex64(scale_fnv.value()) << "\nMATRIX_CELLS "
                  << repairs.size() * scales.size() << " NONNEGATIVE "
                  << nonnegative_cells << " MATRIX_FNV "
                  << hex64(matrix_fnv.value()) << " LITERAL_MATRIX_FNV "
                  << hex64(literal_matrix_fnv.value())
                  << "\nLEAST_COMMON_MARGIN "
                  << decimal(least_margin) << " AT_SCALE " << least_scale
                  << " REPAIR {" << labels(least_repair)
                  << "}\nLEAST_LITERAL_COMMON_MARGIN "
                  << decimal(least_literal_margin) << " AT_SCALE "
                  << least_literal_scale << " REPAIR {"
                  << labels(least_literal_repair) << "}\n";

        const auto bodies = all_bodies();
        std::array<BodyScan, BUDGETS.size()> scans;
        for (std::size_t slot = 0; slot < BUDGETS.size(); ++slot) {
            ensure(common[slot].size() == EXPECTED_COMMON[slot],
                   "common-deck budget census changed");
            Fnv64 budget_fnv;
            for (u32 repair : common[slot]) budget_fnv.word(repair);
            ensure(budget_fnv.value() == EXPECTED_COMMON_FNV[slot],
                   "common-deck budget FNV changed");
            scans[slot] = scan_bodies(common[slot], bodies);
            ensure(scans[slot].failures == EXPECTED_FAILURES[slot],
                   "hostile-body budget census changed");
            ensure(scans[slot].checks == EXPECTED_CHECKS[slot] &&
                       scans[slot].max_checks == EXPECTED_MAX_CHECKS[slot],
                   "body-scan work ledger changed");
            std::cout << "BUDGET " << BUDGETS[slot] << " COMMON "
                      << common[slot].size() << " COMMON_FNV "
                      << hex64(budget_fnv.value()) << " BODY_COUNT "
                      << scans[slot].bodies
                      << " FAILURES " << scans[slot].failures << " CHECKS "
                      << scans[slot].checks << " MAX_CHECKS " << scans[slot].max_checks
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
        ensure(scans.back().failures == 0, "full common deck misses a body");

        constexpr std::size_t base_slot = 1;
        const std::vector<u32> hostiles = failed_bodies(common[base_slot], bodies);
        ensure(hostiles.size() == scans[base_slot].failures,
               "hostile-body list mismatch");
        std::vector<u32> targeted = common[base_slot];
        std::vector<u32> targeted_additions;
        for (u32 body : hostiles) {
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
        ensure(targeted.size() == 3163 && targeted_additions.size() == 5 &&
                   std::equal(targeted_additions.begin(), targeted_additions.end(),
                              EXPECTED_ADDITIONS.begin()),
               "targeted deck construction changed");
        Fnv64 targeted_fnv;
        for (u32 repair : targeted) targeted_fnv.word(repair);
        ensure(targeted_fnv.value() == UINT64_C(0x4e8a91621047de5c),
               "targeted deck FNV changed");
        i128 targeted_least_margin = 0;
        i64 targeted_least_scale = 0;
        u32 targeted_least_repair = 0;
        bool targeted_first = true;
        i128 targeted_literal_least_margin = 0;
        i64 targeted_literal_least_scale = 0;
        u32 targeted_literal_least_repair = 0;
        bool targeted_literal_first = true;
        for (u32 repair : targeted) {
            const auto position = std::find(repairs.begin(), repairs.end(), repair);
            ensure(position != repairs.end(), "targeted repair absent from candidates");
            const Geometry& item = geometry[position - repairs.begin()];
            for (i64 scale : scales) {
                const i128 value = margin(primitive, item, scale);
                const i64 denominator = exact_lcm(D, N56 * scale);
                const i128 conversion =
                    static_cast<i128>(N56) * D * scale / denominator;
                ensure(value % conversion == 0,
                       "targeted literal-margin conversion is nonintegral");
                const i128 literal_value = value / conversion;
                if (targeted_first || value < targeted_least_margin) {
                    targeted_first = false;
                    targeted_least_margin = value;
                    targeted_least_scale = scale;
                    targeted_least_repair = repair;
                }
                if (targeted_literal_first ||
                    literal_value < targeted_literal_least_margin) {
                    targeted_literal_first = false;
                    targeted_literal_least_margin = literal_value;
                    targeted_literal_least_scale = scale;
                    targeted_literal_least_repair = repair;
                }
            }
        }
        std::cout << "TARGETED_FROM " << BUDGETS[base_slot]
                  << " HOSTILE_COUNT " << hostiles.size() << " HOSTILE_BODIES";
        for (u32 body : hostiles) std::cout << " {" << labels(body) << "}";
        std::cout << " ADDITION_COUNT " << targeted_additions.size() << " ADDITIONS";
        for (u32 repair : targeted_additions) {
            const auto candidate_position = std::find(repairs.begin(), repairs.end(), repair);
            const auto common_position =
                std::find(common.back().begin(), common.back().end(), repair);
            ensure(candidate_position != repairs.end(), "target repair missing from candidates");
            ensure(common_position != common.back().end(), "target repair missing from common deck");
            std::cout << " {" << labels(repair) << "} MASK " << hex64(repair)
                      << " CANDIDATE_RANK " << (candidate_position - repairs.begin() + 1)
                      << " COMMON_RANK " << (common_position - common.back().begin() + 1);
        }
        std::cout << " TARGETED_COUNT " << targeted.size() << " TARGETED_FNV "
                  << hex64(targeted_fnv.value()) << " FAILURES "
                  << targeted_scan.failures << " CHECKS " << targeted_scan.checks
                  << " MAX_CHECKS " << targeted_scan.max_checks
                  << " LEAST_TARGETED_MARGIN " << decimal(targeted_least_margin)
                  << " AT_SCALE " << targeted_least_scale << " REPAIR {"
                  << labels(targeted_least_repair)
                  << "} LEAST_LITERAL_TARGETED_MARGIN "
                  << decimal(targeted_literal_least_margin) << " AT_SCALE "
                  << targeted_literal_least_scale << " REPAIR {"
                  << labels(targeted_literal_least_repair) << "}\n";

        ensure(targeted_scan.checks == UINT64_C(486631847) &&
                   targeted_scan.max_checks == UINT64_C(3163),
               "targeted body-scan ledger changed");
        ensure(targeted_least_margin == static_cast<i128>(36651492783360LL) &&
                   targeted_least_scale == 64 &&
                   targeted_least_repair == UINT32_C(0x10053221),
               "targeted primary minimum changed");
        ensure(targeted_literal_least_margin ==
                   static_cast<i128>(10908182376LL) &&
                   targeted_literal_least_scale == 64 &&
                   targeted_literal_least_repair == UINT32_C(0x10053221),
               "targeted literal minimum changed");

        if (argc == 2) write_mask_ledger(argv[1], targeted);
        std::cout << "VERDICT COMMON_DECK_CLOSES_EVERY_BODY_ALL_SCALES_59_128\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << "\n";
        return 1;
    }
}
