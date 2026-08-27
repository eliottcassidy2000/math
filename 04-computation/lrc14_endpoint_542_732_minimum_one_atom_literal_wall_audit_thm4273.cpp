// Clean-room literal joint-wall audit of THM-4273's minimum one-atom
// augmentation.

#define main thm4254_direct_wall_maintained_main
#include "lrc14_endpoint_cascade_direct_wall_body_audit.cpp"
#undef main

#include <set>

namespace {

constexpr u64 EXPECTED_UNION_FNV_LOCAL = UINT64_C(0xce4e76ec11df057c);
constexpr u64 EXPECTED_COMMON_ACTIVE_FNV = UINT64_C(0x829ae906d6b54c9a);
constexpr u32 EXPECTED_CANONICAL = UINT32_C(0x1aa89);

std::vector<u32> union_from_transcripts_local(
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
    require(ledger.state == EXPECTED_UNION_FNV_LOCAL,
            "inherited prefix-union FNV changed");
    return answer;
}

struct LiteralAtoms {
    Geometry geometry;
    std::map<u32, i64> mass;
    i64 total = 0;
};

LiteralAtoms build_literal_atoms(int q, int r) {
    LiteralAtoms answer;
    answer.geometry = build_joint_geometry(q, r);
    for (const Cell& cell : answer.geometry.cells) {
        if (!cell.pair_safe) continue;
        answer.mass[cell.failed_pool] += cell.width;
        answer.total += cell.width;
    }
    require(answer.total == answer.geometry.pair_ticks,
            "literal atom aggregation lost pair-safe mass");
    return answer;
}

i64 literal_mass(const LiteralAtoms& atoms, u32 repair) {
    i64 mass = 0;
    for (const auto& [failure_mask, width] : atoms.mass) {
        if ((failure_mask & ~repair) == 0) mass += width;
    }
    return mass;
}

std::vector<u32> literal_active_deck_local(const std::vector<u32>& candidates,
                                           const LiteralAtoms& atoms) {
    std::vector<u32> active;
    for (u32 repair : candidates) {
        const i64 mass = literal_mass(atoms, repair);
        const i128 margin = static_cast<i128>(63) * mass -
                            static_cast<i128>(4) * atoms.geometry.grid;
        if (margin < 0) continue;
        require(margin > 0, "zero literal margin entered active carrier");
        active.push_back(repair);
    }
    return active;
}

std::vector<u32> collect_failures(const std::vector<u32>& bodies,
                                  const std::vector<u32>& active,
                                  u64& checks) {
    std::vector<u32> failures;
    checks = 0;
    for (u32 body : bodies) {
        bool covered = false;
        for (u32 repair : active) {
            ++checks;
            if ((body & repair) == 0) {
                covered = true;
                break;
            }
        }
        if (!covered) failures.push_back(body);
    }
    return failures;
}

void enumerate_allowed_repairs(const std::vector<int>& allowed, int index,
                               int need, u32 mask,
                               std::vector<u32>& repairs) {
    if (need == 0) {
        repairs.push_back(mask);
        return;
    }
    if (static_cast<int>(allowed.size()) - index < need) return;
    for (int cursor = index;
         cursor <= static_cast<int>(allowed.size()) - need; ++cursor) {
        enumerate_allowed_repairs(allowed, cursor + 1, need - 1,
                                  mask | (u32{1} << allowed[cursor]), repairs);
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "usage: min-atom-direct-wall THM4254_REPLAY_DIR");
        const auto parsed = parse_probe_directory(argv[1]);
        const std::vector<u32> carrier = union_from_transcripts_local(parsed);

        constexpr int q = 542;
        constexpr int r = 732;
        const LiteralAtoms atoms = build_literal_atoms(q, r);
        require(atoms.geometry.grid == INT64_C(301544606313830880) &&
                    atoms.geometry.cells.size() == 9653 &&
                    atoms.mass.size() == 2545,
                "literal joint-wall ledger changed");
        const std::vector<u32> active =
            literal_active_deck_local(carrier, atoms);
        require(active.size() == 3227, "literal active carrier changed");

        std::vector<u32> bodies;
        bodies.reserve(EXPECTED_BODIES);
        enumerate_bodies(0, 9, 0, bodies);
        require(bodies.size() == EXPECTED_BODIES,
                "recursive rank-nine universe changed");
        u64 inherited_checks = 0;
        std::vector<u32> failures =
            collect_failures(bodies, active, inherited_checks);
        require(failures.size() == 2 &&
                    inherited_checks == UINT64_C(474614827),
                "literal hostile failure boundary changed");
        std::sort(failures.begin(), failures.end());
        require(failures[0] == UINT32_C(0x151a5400) &&
                    failures[1] == UINT32_C(0x153a4400),
                "literal hostile bodies disagree with primary");

        const u32 body_union = failures[0] | failures[1];
        std::vector<int> allowed;
        for (int bit = 0; bit < 30; ++bit) {
            if ((body_union & (u32{1} << bit)) == 0) allowed.push_back(bit);
        }
        require(allowed.size() == 20, "common-complement size changed");
        std::vector<u32> common_repairs;
        common_repairs.reserve(125970);
        enumerate_allowed_repairs(allowed, 0, 8, 0, common_repairs);
        require(common_repairs.size() == 125970,
                "literal common-complement rank-eight count changed");

        std::vector<u32> common_active;
        i128 minimum_margin = 0;
        i128 maximum_margin = 0;
        u64 zero_margin = 0;
        for (u32 repair : common_repairs) {
            const i64 mass = literal_mass(atoms, repair);
            const i128 margin = static_cast<i128>(63) * mass -
                                static_cast<i128>(4) * atoms.geometry.grid;
            if (margin < 0) continue;
            if (margin == 0) ++zero_margin;
            if (common_active.empty()) minimum_margin = maximum_margin = margin;
            minimum_margin = std::min(minimum_margin, margin);
            maximum_margin = std::max(maximum_margin, margin);
            common_active.push_back(repair);
        }
        std::sort(common_active.begin(), common_active.end());
        require(common_active.size() == 2172,
                "literal common-active count disagrees with primary");
        require(zero_margin == 0, "literal activation equality appeared");
        require(common_active.front() == EXPECTED_CANONICAL,
                "literal canonical minimum disagrees with primary");
        Fnv active_ledger;
        for (u32 repair : common_active) active_ledger.add(repair);
        require(active_ledger.state == EXPECTED_COMMON_ACTIVE_FNV,
                "literal common-active FNV disagrees with primary");
        std::set<u32> carrier_set(carrier.begin(), carrier.end());
        const u64 carrier_overlap = static_cast<u64>(std::count_if(
            common_active.begin(), common_active.end(),
            [&](u32 repair) { return carrier_set.contains(repair); }));
        require(carrier_overlap == 0,
                "literal common-active set overlaps inherited carrier");

        const i64 canonical_mass = literal_mass(atoms, EXPECTED_CANONICAL);
        const i128 canonical_margin =
            static_cast<i128>(63) * canonical_mass -
            static_cast<i128>(4) * atoms.geometry.grid;
        require(canonical_margin > 0, "literal canonical atom is inactive");
        require(fraction(canonical_mass, atoms.geometry.grid) ==
                    "2395416707526053/37693075789228860",
                "literal canonical exact mass changed");

        std::vector<u32> augmented = active;
        augmented.push_back(EXPECTED_CANONICAL);
        const BodyAudit augmented_body = audit_bodies(bodies, augmented);
        require(augmented_body.failures == 0 &&
                    augmented_body.bodies == EXPECTED_BODIES &&
                    augmented_body.checks == UINT64_C(474614829) &&
                    augmented_body.max_checks == UINT64_C(3228),
                "literal one-atom closure changed");

        std::cout << "LRC14_ENDPOINT_542_732_MIN_ATOM_LITERAL_WALL_THM4273\n";
        std::cout << "PAIR " << q << ',' << r
                  << " GRID " << atoms.geometry.grid
                  << " JOINT_CELLS " << atoms.geometry.cells.size()
                  << " LITERAL_ATOMS " << atoms.mass.size()
                  << " PAIR_TICKS " << atoms.geometry.pair_ticks << '\n';
        std::cout << "INHERITED_CARRIER " << carrier.size()
                  << " ACTIVE " << active.size()
                  << " BODIES " << bodies.size()
                  << " FAILURES " << failures.size()
                  << " CHECKS " << inherited_checks << '\n';
        std::cout << "BODY_A {" << labels(failures[0]) << "} HEX "
                  << std::hex << failures[0] << std::dec << '\n';
        std::cout << "BODY_B {" << labels(failures[1]) << "} HEX "
                  << std::hex << failures[1] << std::dec << '\n';
        std::cout << "COMMON_COMPLEMENT_SIZE " << allowed.size()
                  << " RANK8_CANDIDATES " << common_repairs.size()
                  << " ACTIVE " << common_active.size()
                  << " ZERO_MARGIN " << zero_margin
                  << " INHERITED_CARRIER_OVERLAP " << carrier_overlap
                  << " ACTIVE_MASK_FNV " << std::hex << active_ledger.state
                  << std::dec << '\n';
        std::cout << "CANONICAL_MASK " << std::hex << EXPECTED_CANONICAL
                  << std::dec << " LABELS {" << labels(EXPECTED_CANONICAL) << "}"
                  << " MASS " << fraction(canonical_mass, atoms.geometry.grid)
                  << " MARGIN63_NUM " << decimal(canonical_margin) << '\n';
        std::cout << "MARGIN_RANGE [" << decimal(minimum_margin) << ','
                  << decimal(maximum_margin) << "]\n";
        std::cout << "AUGMENTED_ACTIVE " << augmented.size()
                  << " BODIES " << augmented_body.bodies
                  << " FAILURES " << augmented_body.failures
                  << " CHECKS " << augmented_body.checks
                  << " MAX_CHECKS " << augmented_body.max_checks << '\n';
        std::cout << "VERDICT PASS INDEPENDENT_LITERAL_WALL_ONE_ATOM_MINIMUM\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "DIRECT_WALL_ERROR " << error.what() << '\n';
        return 1;
    }
}
