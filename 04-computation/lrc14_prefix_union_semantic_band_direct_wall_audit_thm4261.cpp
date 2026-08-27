// Clean-room literal-wall audit for THM-4261.
//
// This program imports only the established THM-4254 transcript parser and
// its direct joint-wall machinery.  It reconstructs the 4,675-mask carrier,
// builds fresh walls for every target pair, aggregates literal safe-cell
// widths by their failed-pool mask, thresholds each candidate directly, and
// recursively enumerates every labelled nine-body.  It does not use the
// primitive endpoint cocycle or signed endpoint atoms of the primary checker.

#define main thm4254_direct_wall_audit_main
#include "lrc14_endpoint_cascade_direct_wall_body_audit.cpp"
#undef main

#include <set>

namespace {

constexpr u64 TARGET_FNV_4261 = UINT64_C(0xe923d1494185b820);
constexpr u64 UNION_FNV_4261 = UINT64_C(0xce4e76ec11df057c);
constexpr u64 PAIR_SUMMARY_FNV_4261 = UINT64_C(0xc2657c090f345f51);
constexpr u64 ACTIVE_INCIDENCE_FNV_4261 = UINT64_C(0x247b8bb8de3112c0);
constexpr u64 ACTIVE_INCIDENCES_4261 = UINT64_C(1344725);
constexpr u64 BODY_CASES_4261 = UINT64_C(4249223550);
constexpr u64 BODY_CHECKS_4261 = UINT64_C(125275100068);

std::vector<u32> union_from_transcripts(
    const std::map<std::pair<int, int>, Parsed>& parsed) {
    std::vector<u32> answer;
    std::set<u32> seen;
    for (const auto& key : EXPECTED_PAIRS) {
        const auto found = parsed.find(key);
        require(found != parsed.end(), "inherited transcript disappeared");
        for (u32 repair : found->second.deck) {
            if (seen.insert(repair).second) answer.push_back(repair);
        }
    }
    require(answer.size() == 4675, "inherited prefix union changed");
    Fnv ledger;
    for (u32 repair : answer) ledger.add(repair);
    require(ledger.state == UNION_FNV_4261,
            "inherited prefix-union FNV changed");
    return answer;
}

std::vector<std::pair<int, int>> read_targets_4261() {
    std::vector<std::pair<int, int>> targets;
    int q = 0;
    int r = 0;
    while (std::cin >> q >> r) {
        require(q > 0 && q < r && 733 <= r && r <= 754,
                "target leaves semantic endpoint band");
        if (!targets.empty()) {
            require(targets.back() < std::make_pair(q, r),
                    "targets are not strictly lexicographic");
        }
        targets.emplace_back(q, r);
    }
    require(std::cin.eof(), "malformed target stream");
    require(targets.size() == 297, "semantic target count changed");
    Fnv ledger;
    for (const auto& [q_value, r_value] : targets) {
        ledger.add(static_cast<u64>(q_value));
        ledger.add(static_cast<u64>(r_value));
    }
    require(ledger.state == TARGET_FNV_4261,
            "semantic target FNV changed");
    return targets;
}

struct LiteralDeck {
    std::vector<u32> masks;
    i128 minimum_margin = 0;
    u32 minimum_repair = 0;
    std::size_t atom_count = 0;
    i64 pair_ticks = 0;
    i64 grid = 0;
    std::size_t joint_cells = 0;
};

LiteralDeck literal_active_deck(const std::vector<u32>& candidates,
                                int q, int r) {
    const Geometry geometry = build_joint_geometry(q, r);
    std::map<u32, i64> atoms;
    i64 atom_total = 0;
    for (const Cell& cell : geometry.cells) {
        if (!cell.pair_safe) continue;
        atoms[cell.failed_pool] += cell.width;
        atom_total += cell.width;
    }
    require(atom_total == geometry.pair_ticks,
            "literal atom aggregation lost pair-safe mass");

    LiteralDeck answer;
    answer.atom_count = atoms.size();
    answer.pair_ticks = geometry.pair_ticks;
    answer.grid = geometry.grid;
    answer.joint_cells = geometry.cells.size();
    for (u32 repair : candidates) {
        i64 mass = 0;
        for (const auto& [failure_mask, width] : atoms) {
            if ((failure_mask & ~repair) == 0) mass += width;
        }
        const i128 margin = static_cast<i128>(63) * mass -
                            static_cast<i128>(4) * geometry.grid;
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
                "usage: direct-wall-audit THM4254_REPLAY_DIR < TARGET_PAIRS");
        const auto parsed = parse_probe_directory(argv[1]);
        const std::vector<u32> candidates = union_from_transcripts(parsed);
        const std::vector<std::pair<int, int>> targets = read_targets_4261();

        std::vector<u32> bodies;
        bodies.reserve(EXPECTED_BODIES);
        enumerate_bodies(0, 9, 0, bodies);
        require(bodies.size() == EXPECTED_BODIES,
                "recursive body universe changed");

        Fnv pair_summary_ledger;
        Fnv active_incidence_ledger;
        u64 total_active = 0;
        u64 total_body_cases = 0;
        u64 total_body_checks = 0;
        u64 total_joint_cells = 0;
        u64 total_literal_atoms = 0;

        std::cout << "LRC14_PREFIX_UNION_SEMANTIC_BAND_DIRECT_WALL_V1\n";
        std::cout << "TARGETS 297 ENDPOINTS 733..754 FNV " << std::hex
                  << TARGET_FNV_4261 << " PREFIX_UNION 4675 FNV "
                  << UNION_FNV_4261 << std::dec << '\n';

        for (const auto& [q, r] : targets) {
            const LiteralDeck active = literal_active_deck(candidates, q, r);
            require(!active.masks.empty(), "pair has no active union repair");
            const BodyAudit body = audit_bodies(bodies, active.masks);
            require(body.bodies == EXPECTED_BODIES && body.failures == 0,
                    "literal-wall semantic-band closure failed");

            pair_summary_ledger.add(static_cast<u64>(q));
            pair_summary_ledger.add(static_cast<u64>(r));
            pair_summary_ledger.add(active.masks.size());
            active_incidence_ledger.add(static_cast<u64>(q));
            active_incidence_ledger.add(static_cast<u64>(r));
            active_incidence_ledger.add(active.masks.size());
            for (u32 repair : active.masks) active_incidence_ledger.add(repair);
            total_active += active.masks.size();
            total_body_cases += body.bodies;
            total_body_checks += body.checks;
            total_joint_cells += active.joint_cells;
            total_literal_atoms += active.atom_count;

            std::cout << "EDGE " << q << ',' << r
                      << " GRID " << active.grid
                      << " JOINT_CELLS " << active.joint_cells
                      << " LITERAL_ATOMS " << active.atom_count
                      << " PAIR_TICKS " << active.pair_ticks
                      << " ACTIVE_UNION " << active.masks.size()
                      << " MIN_MARGIN_NUM " << decimal(active.minimum_margin)
                      << " MIN_REPAIR {" << labels(active.minimum_repair) << "}"
                      << " BODIES " << body.bodies
                      << " FAILURES " << body.failures
                      << " CHECKS " << body.checks
                      << " MAX_CHECKS " << body.max_checks << '\n';
        }

        require(total_active == ACTIVE_INCIDENCES_4261,
                "active-incidence count disagrees with primary");
        require(total_body_cases == BODY_CASES_4261,
                "body-case count disagrees with primary");
        require(total_body_checks == BODY_CHECKS_4261,
                "body-check count disagrees with primary");
        require(pair_summary_ledger.state == PAIR_SUMMARY_FNV_4261,
                "pair-summary FNV disagrees with primary");
        require(active_incidence_ledger.state == ACTIVE_INCIDENCE_FNV_4261,
                "active-incidence FNV disagrees with primary");

        const LiteralDeck hostile = literal_active_deck(candidates, 542, 732);
        const BodyAudit hostile_body = audit_bodies(bodies, hostile.masks);
        require(hostile.masks.size() == 3227 && hostile_body.failures == 2 &&
                    hostile_body.checks == UINT64_C(474614827) &&
                    hostile_body.max_checks == hostile.masks.size(),
                "literal-wall endpoint-732 hostile changed");

        std::cout << "SUMMARY PAIRS " << targets.size()
                  << " JOINT_CELLS " << total_joint_cells
                  << " LITERAL_ATOMS " << total_literal_atoms
                  << " ACTIVE_INCIDENCES " << total_active
                  << " BODY_CASES " << total_body_cases
                  << " BODY_CHECKS " << total_body_checks
                  << " PAIR_SUMMARY_FNV " << std::hex
                  << pair_summary_ledger.state
                  << " ACTIVE_INCIDENCE_FNV "
                  << active_incidence_ledger.state << std::dec << '\n';
        std::cout << "HOSTILE 542,732 ACTIVE_UNION " << hostile.masks.size()
                  << " FAILURES " << hostile_body.failures
                  << " CHECKS " << hostile_body.checks
                  << " MAX_CHECKS " << hostile_body.max_checks << '\n';
        std::cout << "VERDICT PASS FRESH_JOINT_WALLS ALL_297_EVERY_BODY; "
                     "ENDPOINT_732_UNION_NOT_UNIFORM\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "DIRECT_WALL_ERROR " << error.what() << '\n';
        return 1;
    }
}
