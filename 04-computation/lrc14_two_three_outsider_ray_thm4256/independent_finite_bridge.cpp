#include "independent_common.hpp"

namespace {
using namespace ray_audit;

constexpr i64 FIRST_G = 146;
constexpr i64 LAST_G = 1015;
constexpr std::size_t CANDIDATES = 4096;
constexpr u32 ROBUST_1015_BODY =
    (u32{1} << 13) | (u32{1} << 14) | (u32{1} << 21) |
    (u32{1} << 24) | (u32{1} << 25) | (u32{1} << 26) |
    (u32{1} << 27) | (u32{1} << 28) | (u32{1} << 29);

struct OrderedRepair {
    u64 key = 0;
    u32 mask = 0;
};

std::vector<u32> globally_first_repairs() {
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
    std::sort(universe.begin(), universe.end(), [](const OrderedRepair& a,
                                                   const OrderedRepair& b) {
        return a.key != b.key ? a.key < b.key : a.mask < b.mask;
    });
    std::vector<u32> answer;
    answer.reserve(CANDIDATES);
    for (std::size_t i = 0; i < CANDIDATES; ++i) answer.push_back(universe[i].mask);
    return answer;
}

struct Witness {
    u64 checks = 0;
    u32 repair = 0;
};

Witness first_disjoint(u32 body, const std::vector<u32>& deck) {
    Witness answer;
    for (u32 repair : deck) {
        ++answer.checks;
        if ((repair & body) == 0) {
            answer.repair = repair;
            return answer;
        }
    }
    return answer;
}

void smoke_mode(const std::vector<Geometry>& geometry, u64 candidate_fnv) {
    constexpr std::array<i64, 6> scales = {1, 145, 146, 256, 1015, 1016};
    Fnv64 ledger;
    for (i64 g : scales) {
        for (std::size_t index : {std::size_t(0), std::size_t(17),
                                  std::size_t(873), geometry.size() - 1}) {
            const i128 mass = exact_mass_numerator(geometry[index], g);
            const i128 margin = static_cast<i128>(63) * mass
                              - static_cast<i128>(4) * N * D * g;
            ledger.word(static_cast<u64>(g));
            ledger.word(geometry[index].repair);
            ledger.word(static_cast<u64>(mass));
            ledger.word(static_cast<u64>(mass >> 64));
            ledger.word(static_cast<u64>(margin));
            ledger.word(static_cast<u64>(margin >> 64));
        }
    }
    std::cout << "LRC_23_RAY_INDEPENDENT_FINITE_SMOKE_V1\n"
              << "CANDIDATES " << geometry.size() << " CANDIDATE_FNV "
              << hex64(candidate_fnv) << " SAMPLE_LEDGER "
              << hex64(ledger.value()) << "\n"
              << "VERDICT EXACT_DIRECT_ARC_INTEGRATOR_SMOKE_PASS\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        using namespace ray_audit;
        const bool smoke = argc == 2 && std::string(argv[1]) == "--smoke";
        ensure(argc == 1 || smoke, "usage: independent_finite_bridge [--smoke]");
        check_primitive_constants();
        const std::vector<Cell> cells = build_pool_arrangement();
        const std::vector<u32> repairs = globally_first_repairs();
        Fnv64 candidate_ledger;
        std::vector<Geometry> geometry;
        geometry.reserve(repairs.size());
        for (u32 repair : repairs) {
            candidate_ledger.word(repair);
            geometry.push_back(direct_geometry(repair, cells));
        }
        ensure(candidate_ledger.value() == UINT64_C(0x225a0c689159c41e),
               "candidate order ledger mismatch");
        if (smoke) {
            smoke_mode(geometry, candidate_ledger.value());
            return 0;
        }

        const std::vector<u32> bodies = all_bodies();
        Fnv64 transcript;
        u64 total_bodies = 0;
        u64 total_checks = 0;
        u64 minimum_active = std::numeric_limits<u64>::max();
        i64 minimum_active_g = 0;
        u64 maximum_checks = 0;
        i64 maximum_checks_g = 0;
        u32 maximum_checks_body = 0;
        Witness sharpness_witness;
        i128 sharpness_margin = 0;

        std::cout << "LRC_23_RAY_INDEPENDENT_FINITE_BRIDGE_V1 FIRST_G "
                  << FIRST_G << " LAST_G " << LAST_G << " CANDIDATES "
                  << geometry.size() << " CANDIDATE_FNV "
                  << hex64(candidate_ledger.value()) << std::endl;

        for (i64 g = FIRST_G; g <= LAST_G; ++g) {
            std::vector<u32> active;
            active.reserve(geometry.size());
            Fnv64 active_ledger;
            for (const Geometry& item : geometry) {
                if (exact_margin(item, g) >= 0) {
                    active.push_back(item.repair);
                    active_ledger.word(item.repair);
                }
            }

            const BodyScan scan = scan_bodies(active, bodies);
            transcript.word(static_cast<u64>(g));
            transcript.word(active.size());
            transcript.word(active_ledger.value());
            transcript.word(scan.failures);
            transcript.word(scan.checks);
            transcript.word(scan.max_checks);
            transcript.word(scan.worst_body);
            total_bodies += scan.bodies;
            total_checks += scan.checks;
            if (active.size() < minimum_active) {
                minimum_active = active.size();
                minimum_active_g = g;
            }
            if (scan.max_checks > maximum_checks) {
                maximum_checks = scan.max_checks;
                maximum_checks_g = g;
                maximum_checks_body = scan.worst_body;
            }
            ensure(scan.failures == 0, "uncovered body in finite bridge");

            if (g == 1015) {
                sharpness_witness = first_disjoint(ROBUST_1015_BODY, active);
                ensure(sharpness_witness.repair != 0,
                       "robust-boundary body lacks exact witness");
                const auto found = std::find(repairs.begin(), repairs.end(),
                                             sharpness_witness.repair);
                ensure(found != repairs.end(), "sharpness witness left candidate deck");
                const std::size_t index = static_cast<std::size_t>(found - repairs.begin());
                sharpness_margin = exact_margin(geometry[index], g);
                ensure(sharpness_margin >= 0, "sharpness witness is inactive");
            }

            if (g == FIRST_G || g == LAST_G || g % 100 == 0) {
                std::cout << "CHECKPOINT G " << g << " ACTIVE " << active.size()
                          << " ACTIVE_FNV " << hex64(active_ledger.value())
                          << " MAX_CHECKS " << scan.max_checks << std::endl;
            }
        }

        ensure(total_bodies == UINT64_C(12447220500),
               "body-scale census mismatch");
        ensure(transcript.value() == UINT64_C(0x266a67dd758f4e75),
               "finite transcript ledger mismatch");
        std::cout << "SCALES " << (LAST_G - FIRST_G + 1)
                  << " TOTAL_BODIES " << total_bodies
                  << " TOTAL_CHECKS " << total_checks << '\n';
        std::cout << "MIN_ACTIVE " << minimum_active << " AT_G "
                  << minimum_active_g << " MAX_CHECKS " << maximum_checks
                  << " AT_G " << maximum_checks_g << " BODY {"
                  << labels(maximum_checks_body) << "}\n";
        std::cout << "ROBUST_1015_BODY {" << labels(ROBUST_1015_BODY)
                  << "} EXACT_WITNESS_CHECK " << sharpness_witness.checks
                  << " REPAIR {" << labels(sharpness_witness.repair)
                  << "} EXACT_MARGIN_NUM " << decimal(sharpness_margin) << '\n';
        std::cout << "TRANSCRIPT_FNV " << hex64(transcript.value()) << '\n';
        std::cout << "VERDICT EVERY_BODY_CLOSED_AT_EVERY_INTEGER_SCALE_IN_BRIDGE\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
