// Exact signed-endpoint-cocycle verifier for the THM-4261 semantic band.
//
// The 4,675 candidate repairs are the set-theoretic union of the 59 short
// THM-4254 prefixes.  For each target pair supplied on stdin, this checker
// recomputes every candidate mass from signed endpoint atoms, retains only
// pair-active candidates, and exhausts all C(30,9) labelled bodies.  It does
// not run the discovery superset transform or assume that a candidate active
// for one pair is active for another.

#define CASCADE_PREFIX_LIBRARY_ONLY
#include "cascade_prefix_union_exact_verifier.cpp"
#undef CASCADE_PREFIX_LIBRARY_ONLY

#include <iostream>

namespace {

constexpr u64 EXPECTED_TARGET_FNV = UINT64_C(0xe923d1494185b820);
constexpr u64 EXPECTED_UNION_FNV = UINT64_C(0xce4e76ec11df057c);
constexpr u64 EXPECTED_PAIR_SUMMARY_FNV = UINT64_C(0xc2657c090f345f51);
constexpr u64 EXPECTED_ACTIVE_INCIDENCE_FNV = UINT64_C(0x247b8bb8de3112c0);
constexpr u64 EXPECTED_ACTIVE_INCIDENCES = UINT64_C(1344725);
constexpr u64 EXPECTED_BODY_CASES = UINT64_C(4249223550);
constexpr u64 EXPECTED_BODY_CHECKS = UINT64_C(125275100068);

std::vector<u32> inherited_prefix_union(const std::filesystem::path& replay) {
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
    require(ledger.state == EXPECTED_UNION_FNV,
            "inherited prefix-union FNV changed");
    return answer;
}

std::vector<PairKey> read_targets() {
    std::vector<PairKey> targets;
    int q = 0;
    int r = 0;
    while (std::cin >> q >> r) {
        require(q > 0 && q < r && 733 <= r && r <= 754,
                "target leaves semantic endpoint band");
        if (!targets.empty()) {
            const PairKey prior = targets.back();
            require(prior.q < q || (prior.q == q && prior.r < r),
                    "targets are not strictly lexicographic");
        }
        targets.push_back({q, r});
    }
    require(std::cin.eof(), "malformed target stream");
    require(targets.size() == 297, "semantic target count changed");
    FnvLocal ledger;
    for (const PairKey key : targets) {
        ledger.add(static_cast<u64>(key.q));
        ledger.add(static_cast<u64>(key.r));
    }
    require(ledger.state == EXPECTED_TARGET_FNV,
            "semantic target FNV changed");
    return targets;
}

struct ActiveDeck {
    std::vector<u32> masks;
    i128 minimum_margin = 0;
    u32 minimum_repair = 0;
};

ActiveDeck active_union_deck(const std::vector<u32>& candidates,
                             const std::vector<Cell>& pool_cells,
                             int q, int r) {
    const i64 g = std::gcd(q, r);
    const PrimitivePair primitive = build_primitive(q / g, r / g);
    const AtomData atoms = build_cocycle_atoms(pool_cells, primitive, g);
    const i128 denominator = static_cast<i128>(primitive.grid) * g * COMMON;
    ActiveDeck answer;
    for (u32 repair : candidates) {
        const i128 mass = direct_atom_mass(atoms, repair);
        const i128 margin = static_cast<i128>(63) * mass -
                            static_cast<i128>(4) * denominator;
        if (margin < 0) continue;
        require(margin > 0, "zero active margin entered the certificate");
        if (answer.masks.empty() || margin < answer.minimum_margin) {
            answer.minimum_margin = margin;
            answer.minimum_repair = repair;
        }
        answer.masks.push_back(repair);
    }
    return answer;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2,
                "usage: exact-verifier THM4254_PACKET_ROOT < TARGET_PAIRS");
        const std::filesystem::path replay =
            std::filesystem::path(argv[1]) / "replay_band";
        const std::vector<u32> candidates = inherited_prefix_union(replay);
        const std::vector<PairKey> targets = read_targets();
        const std::vector<Cell> pool_cells = build_pool_cells();
        require(pool_cells.size() == 7133, "fixed-pool cell count changed");

        FnvLocal pair_summary_ledger;
        FnvLocal active_incidence_ledger;
        u64 total_active = 0;
        u64 total_bodies = 0;
        u64 total_checks = 0;
        std::size_t minimum_active = candidates.size();
        std::size_t maximum_active = 0;
        PairKey minimum_active_pair{};
        PairKey maximum_active_pair{};

        std::cout << "LRC14_PREFIX_UNION_SEMANTIC_BAND_EXACT_V1\n";
        std::cout << "TARGETS 297 ENDPOINTS 733..754 FNV " << std::hex
                  << EXPECTED_TARGET_FNV << " PREFIX_UNION 4675 FNV "
                  << EXPECTED_UNION_FNV << std::dec << '\n';

        for (const PairKey key : targets) {
            const ActiveDeck active =
                active_union_deck(candidates, pool_cells, key.q, key.r);
            require(!active.masks.empty(), "pair has no active union repair");
            const ScanResult scan = scan_bodies(active.masks);
            require(scan.bodies == EXPECTED_BODIES && scan.failures == 0,
                    "semantic-band body closure failed");

            pair_summary_ledger.add(static_cast<u64>(key.q));
            pair_summary_ledger.add(static_cast<u64>(key.r));
            pair_summary_ledger.add(active.masks.size());
            active_incidence_ledger.add(static_cast<u64>(key.q));
            active_incidence_ledger.add(static_cast<u64>(key.r));
            active_incidence_ledger.add(active.masks.size());
            for (u32 repair : active.masks) active_incidence_ledger.add(repair);
            total_active += active.masks.size();
            total_bodies += scan.bodies;
            total_checks += scan.checks;
            if (active.masks.size() < minimum_active) {
                minimum_active = active.masks.size();
                minimum_active_pair = key;
            }
            if (active.masks.size() > maximum_active) {
                maximum_active = active.masks.size();
                maximum_active_pair = key;
            }

            std::cout << "EDGE " << key.q << ',' << key.r
                      << " ACTIVE_UNION " << active.masks.size()
                      << " MIN_MARGIN_NUM " << decimal(active.minimum_margin)
                      << " MIN_REPAIR " << mask_hex(active.minimum_repair)
                      << " BODIES " << scan.bodies
                      << " FAILURES " << scan.failures
                      << " CHECKS " << scan.checks
                      << " MAX_CHECKS " << scan.max_checks << '\n';
        }

        require(total_active == EXPECTED_ACTIVE_INCIDENCES,
                "active-incidence count changed");
        require(total_bodies == EXPECTED_BODY_CASES,
                "body-case count changed");
        require(total_checks == EXPECTED_BODY_CHECKS,
                "body-check count changed");
        require(pair_summary_ledger.state == EXPECTED_PAIR_SUMMARY_FNV,
                "pair-summary FNV changed");
        require(active_incidence_ledger.state == EXPECTED_ACTIVE_INCIDENCE_FNV,
                "active-incidence FNV changed");

        const ActiveDeck hostile =
            active_union_deck(candidates, pool_cells, 542, 732);
        const ScanResult hostile_scan = scan_bodies(hostile.masks);
        require(hostile.masks.size() == 3227 && hostile_scan.failures == 2 &&
                    hostile_scan.checks == UINT64_C(474614827) &&
                    hostile_scan.max_checks == hostile.masks.size() &&
                    hostile_scan.first_failure == UINT32_C(0x151a5400),
                "endpoint-732 hostile boundary changed");

        std::cout << "SUMMARY PAIRS " << targets.size()
                  << " ACTIVE_INCIDENCES " << total_active
                  << " MIN_ACTIVE " << minimum_active << " AT "
                  << minimum_active_pair.q << ',' << minimum_active_pair.r
                  << " MAX_ACTIVE " << maximum_active << " AT "
                  << maximum_active_pair.q << ',' << maximum_active_pair.r
                  << " BODY_CASES " << total_bodies
                  << " BODY_CHECKS " << total_checks
                  << " PAIR_SUMMARY_FNV " << std::hex
                  << pair_summary_ledger.state
                  << " ACTIVE_INCIDENCE_FNV "
                  << active_incidence_ledger.state << std::dec << '\n';
        std::cout << "HOSTILE 542,732 ACTIVE_UNION " << hostile.masks.size()
                  << " FAILURES " << hostile_scan.failures
                  << " FIRST_FAILURE " << mask_hex(hostile_scan.first_failure)
                  << " CHECKS " << hostile_scan.checks
                  << " MAX_CHECKS " << hostile_scan.max_checks << '\n';
        std::cout << "VERDICT PASS ALL_297_EVERY_BODY; "
                     "ENDPOINT_732_UNION_NOT_UNIFORM\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
