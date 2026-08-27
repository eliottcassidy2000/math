// Exact endpoint-cocycle audit for the two THM-4261 hostile
// bodies at (542,732).  This imports the maintained proof checker as a
// library, but changes the question: exhaust the rank-eight common complement
// of the two missed bodies and classify minimum carrier augmentation size.

#define CASCADE_PREFIX_LIBRARY_ONLY
#include "cascade_prefix_union_exact_verifier.cpp"
#undef CASCADE_PREFIX_LIBRARY_ONLY

#include <set>

namespace {

constexpr u64 EXPECTED_UNION_FNV_LOCAL = UINT64_C(0xce4e76ec11df057c);
constexpr u64 EXPECTED_COMMON_ACTIVE_FNV = UINT64_C(0x829ae906d6b54c9a);
constexpr u32 EXPECTED_BODY_A = UINT32_C(0x151a5400);
constexpr u32 EXPECTED_BODY_B = UINT32_C(0x153a4400);
constexpr u32 EXPECTED_CANONICAL = UINT32_C(0x1aa89);

std::vector<u32> inherited_prefix_union_local(
    const std::filesystem::path& replay) {
    std::vector<u32> answer;
    std::set<u32> seen;
    for (const PairKey key : BAND) {
        const FrozenPrefix frozen = read_prefix(
            replay / (std::to_string(key.q) + "_" +
                      std::to_string(key.r) + ".out"), key);
        for (u32 repair : frozen.masks) {
            if (seen.insert(repair).second) answer.push_back(repair);
        }
    }
    require(answer.size() == 4675, "inherited prefix union changed");
    FnvLocal ledger;
    for (u32 repair : answer) ledger.add(repair);
    require(ledger.state == EXPECTED_UNION_FNV_LOCAL,
            "inherited prefix-union FNV changed");
    return answer;
}

struct ActiveDeckLocal {
    std::vector<u32> masks;
};

ActiveDeckLocal active_union_deck_local(const std::vector<u32>& candidates,
                                        const std::vector<Cell>& pool_cells,
                                        int q, int r) {
    const i64 g = std::gcd(q, r);
    const PrimitivePair primitive = build_primitive(q / g, r / g);
    const AtomData atoms = build_cocycle_atoms(pool_cells, primitive, g);
    const i128 denominator = static_cast<i128>(primitive.grid) * g * COMMON;
    ActiveDeckLocal answer;
    for (u32 repair : candidates) {
        const i128 mass = direct_atom_mass(atoms, repair);
        const i128 margin = static_cast<i128>(63) * mass -
                            static_cast<i128>(4) * denominator;
        if (margin < 0) continue;
        require(margin > 0, "zero active margin entered inherited carrier");
        answer.masks.push_back(repair);
    }
    return answer;
}

std::string labels_local(u32 mask) {
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

std::vector<u32> exact_failures(const std::vector<u32>& active,
                                u64& checks) {
    std::vector<u32> failures;
    checks = 0;
    u32 body = (u32{1} << 9) - 1;
    const u32 limit = u32{1} << 30;
    u64 bodies = 0;
    while (body != 0 && body < limit) {
        bool covered = false;
        for (u32 repair : active) {
            ++checks;
            if ((body & repair) == 0) {
                covered = true;
                break;
            }
        }
        if (!covered) failures.push_back(body);
        ++bodies;
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    require(bodies == EXPECTED_BODIES, "rank-nine body universe changed");
    return failures;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "usage: min-atom-primary THM4254_PACKET_ROOT");
        init_choose8_local();
        const std::filesystem::path replay =
            std::filesystem::path(argv[1]) / "replay_band";
        const std::vector<u32> candidates = inherited_prefix_union_local(replay);
        const std::vector<Cell> pool_cells = build_pool_cells();

        constexpr int q = 542;
        constexpr int r = 732;
        const ActiveDeckLocal inherited =
            active_union_deck_local(candidates, pool_cells, q, r);
        require(inherited.masks.size() == 3227,
                "hostile inherited active count changed");

        u64 inherited_checks = 0;
        const std::vector<u32> failures =
            exact_failures(inherited.masks, inherited_checks);
        require(failures.size() == 2, "hostile failure count changed");
        require(inherited_checks == UINT64_C(474614827),
                "hostile body-check count changed");

        const u32 body_a = failures[0];
        const u32 body_b = failures[1];
        require(body_a == EXPECTED_BODY_A, "first hostile body changed");
        require(body_b == EXPECTED_BODY_B, "second hostile body changed");
        const u32 body_union = body_a | body_b;
        const u32 body_intersection = body_a & body_b;
        const u32 pool_mask = (u32{1} << 30) - 1;
        const u32 common_complement = pool_mask & ~body_union;
        const int complement_size = std::popcount(common_complement);

        const i64 g = std::gcd(q, r);
        const PrimitivePair primitive = build_primitive(q / g, r / g);
        const AtomData atoms = build_cocycle_atoms(pool_cells, primitive, g);
        const i128 denominator = static_cast<i128>(primitive.grid) * g * COMMON;

        u64 common_candidates = 0;
        u64 common_active = 0;
        u64 common_zero_margin = 0;
        u64 common_carrier_overlap = 0;
        u32 canonical = 0;
        i128 canonical_mass = 0;
        i128 minimum_margin = 0;
        i128 maximum_margin = 0;
        u32 minimum_margin_mask = 0;
        u32 maximum_margin_mask = 0;
        FnvLocal active_mask_ledger;

        u32 repair = (u32{1} << 8) - 1;
        const u32 limit = u32{1} << 30;
        while (repair != 0 && repair < limit) {
            if ((repair & body_union) == 0) {
                ++common_candidates;
                const i128 mass = direct_atom_mass(atoms, repair);
                const i128 margin = static_cast<i128>(63) * mass -
                                    static_cast<i128>(4) * denominator;
                if (margin >= 0) {
                    ++common_active;
                    if (std::find(candidates.begin(), candidates.end(), repair) !=
                        candidates.end()) {
                        ++common_carrier_overlap;
                    }
                    active_mask_ledger.add(repair);
                    if (margin == 0) ++common_zero_margin;
                    if (canonical == 0) {
                        canonical = repair;
                        canonical_mass = mass;
                        minimum_margin = maximum_margin = margin;
                        minimum_margin_mask = maximum_margin_mask = repair;
                    }
                    if (margin < minimum_margin ||
                        (margin == minimum_margin && repair < minimum_margin_mask)) {
                        minimum_margin = margin;
                        minimum_margin_mask = repair;
                    }
                    if (margin > maximum_margin ||
                        (margin == maximum_margin && repair < maximum_margin_mask)) {
                        maximum_margin = margin;
                        maximum_margin_mask = repair;
                    }
                }
            }
            const u32 next = next_combination(repair);
            if (next <= repair) break;
            repair = next;
        }

        require(complement_size == 20, "common complement size changed");
        require(common_candidates == choose8_local[complement_size][8],
                "common-complement enumeration changed");
        require(common_candidates == UINT64_C(125970),
                "common-complement candidate count changed");
        require(common_active == UINT64_C(2172),
                "common-active repair count changed");
        require(common_zero_margin == 0, "common-active equality appeared");
        require(active_mask_ledger.state == EXPECTED_COMMON_ACTIVE_FNV,
                "common-active repair FNV changed");
        require(canonical == EXPECTED_CANONICAL,
                "canonical minimum repair changed");
        require(common_carrier_overlap == 0,
                "common active repair unexpectedly overlaps inherited carrier");
        require(std::find(candidates.begin(), candidates.end(), canonical) ==
                    candidates.end(),
                "canonical new atom unexpectedly belongs to inherited carrier");
        require((canonical & body_a) == 0 && (canonical & body_b) == 0,
                "canonical atom does not cover both failures");

        std::vector<u32> augmented = inherited.masks;
        augmented.push_back(canonical);
        const ScanResult augmented_scan = scan_bodies(augmented);
        require(augmented_scan.bodies == EXPECTED_BODIES &&
                    augmented_scan.failures == 0,
                "one-atom augmentation did not close every body");
        require(augmented_scan.checks == UINT64_C(474614829) &&
                    augmented_scan.max_checks == UINT64_C(3228),
                "augmented body-scan work ledger changed");

        const Reduced canonical_fraction = reduced(canonical_mass, denominator);
        require(canonical_fraction.num ==
                    static_cast<i128>(2395416707526053LL) &&
                    canonical_fraction.den ==
                    static_cast<i128>(37693075789228860LL),
                "canonical exact mass changed");
        require(minimum_margin == static_cast<i128>(592365863849760LL) &&
                    minimum_margin_mask == UINT32_C(0x08c12a40),
                "minimum common-active margin changed");
        require(maximum_margin ==
                    static_cast<i128>(1988586177999240556LL) * 10 + 8 &&
                    maximum_margin_mask == UINT32_C(0x08412a81),
                "maximum common-active margin changed");

        std::cout << "LRC14_ENDPOINT_542_732_MIN_ATOM_PRIMARY_THM4273\n";
        std::cout << "PAIR " << q << ',' << r
                  << " PRIMITIVE " << primitive.u << ':' << primitive.v
                  << " SCALE " << g
                  << " GRID " << primitive.grid
                  << " POOL_CELLS " << pool_cells.size() << '\n';
        std::cout << "POOL_SIZE 30 REPAIR_RANK 8 BODY_RANK 9 REPAIRS "
                  << EXPECTED_REPAIRS << " BODIES " << EXPECTED_BODIES << '\n';
        std::cout << "INHERITED_CARRIER " << candidates.size()
                  << " ACTIVE " << inherited.masks.size()
                  << " FAILURES " << failures.size()
                  << " CHECKS " << inherited_checks << '\n';
        std::cout << "BODY_A_HEX " << mask_hex(body_a)
                  << " LABELS {" << labels_local(body_a) << "}\n";
        std::cout << "BODY_B_HEX " << mask_hex(body_b)
                  << " LABELS {" << labels_local(body_b) << "}\n";
        std::cout << "INTERSECTION_SIZE " << std::popcount(body_intersection)
                  << " LABELS {" << labels_local(body_intersection) << "}"
                  << " UNION_SIZE " << std::popcount(body_union)
                  << " LABELS {" << labels_local(body_union) << "}\n";
        std::cout << "COMMON_COMPLEMENT_SIZE " << complement_size
                  << " LABELS {" << labels_local(common_complement) << "}"
                  << " RANK8_CANDIDATES " << common_candidates << '\n';
        std::cout << "COMMON_ACTIVE " << common_active
                  << " ZERO_MARGIN " << common_zero_margin
                  << " INHERITED_CARRIER_OVERLAP " << common_carrier_overlap
                  << " ACTIVE_MASK_FNV " << std::hex << active_mask_ledger.state
                  << std::dec << '\n';
        std::cout << "CANONICAL_MASK " << mask_hex(canonical)
                  << " COLEX_RANK0 " << colex_rank8_local(canonical)
                  << " LABELS {" << labels_local(canonical) << "}"
                  << " MASS " << show(canonical_fraction)
                  << " MARGIN63_NUM "
                  << decimal(static_cast<i128>(63) * canonical_mass -
                             static_cast<i128>(4) * denominator) << '\n';
        std::cout << "MIN_MARGIN " << decimal(minimum_margin)
                  << " MASK " << mask_hex(minimum_margin_mask)
                  << " LABELS {" << labels_local(minimum_margin_mask) << "}\n";
        std::cout << "MAX_MARGIN " << decimal(maximum_margin)
                  << " MASK " << mask_hex(maximum_margin_mask)
                  << " LABELS {" << labels_local(maximum_margin_mask) << "}\n";
        std::cout << "AUGMENTED_ACTIVE " << augmented.size()
                  << " BODIES " << augmented_scan.bodies
                  << " FAILURES " << augmented_scan.failures
                  << " CHECKS " << augmented_scan.checks
                  << " MAX_CHECKS " << augmented_scan.max_checks << '\n';
        std::cout << "MINIMUM_AUGMENTATION_CARDINALITY 1\n";
        std::cout << "VERDICT PASS EXACT_TWO_FAILURE_CLASSIFICATION; "
                     "ONE_COMMON_ACTIVE_ATOM_IS_CARDINALITY_MINIMAL\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "PRIMARY_ERROR " << error.what() << '\n';
        return 1;
    }
}
