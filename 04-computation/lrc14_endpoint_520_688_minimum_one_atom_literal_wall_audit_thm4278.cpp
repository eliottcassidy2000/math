// Clean-room literal-wall audit for THM-4278 at (520,688). This does not use
// endpoint cocycle atoms or the primary rank-eight zeta transform.

#define main thm4254_direct_wall_maintained_main
#include "lrc14_endpoint_cascade_direct_wall_body_audit.cpp"
#undef main

#include <set>

namespace {

constexpr u64 EXPECTED_BASE_FNV = UINT64_C(0xce4e76ec11df057c);
constexpr u64 EXPECTED_FINAL_FNV = UINT64_C(0xe08b227730f6793c);
constexpr u64 EXPECTED_COMMON_FNV = UINT64_C(0xed1bfbaf6eaa06a3);
constexpr u32 EXPECTED_CANONICAL = UINT32_C(0x00048ec1);

std::vector<u32> base_union(
    const std::map<std::pair<int, int>, Parsed>& parsed) {
    std::vector<u32> answer;
    std::set<u32> seen;
    for (const auto& key : EXPECTED_PAIRS) {
        const auto found = parsed.find(key);
        require(found != parsed.end(), "base transcript disappeared");
        for (u32 repair : found->second.deck) {
            if (seen.insert(repair).second) answer.push_back(repair);
        }
    }
    Fnv ledger;
    for (u32 repair : answer) ledger.add(repair);
    require(answer.size() == 4675 && ledger.state == EXPECTED_BASE_FNV,
            "base carrier changed");
    return answer;
}

struct PrefixLiteral {
    int q = 0;
    int r = 0;
    u64 stated_size = 0;
    u64 stated_fnv = 0;
    std::vector<u32> masks;
};

u64 parse_hex_literal(const std::string& word) {
    std::size_t used = 0;
    const u64 value = std::stoull(word, &used, 16);
    require(used == word.size(), "bad hexadecimal token");
    return value;
}

PrefixLiteral read_prefix_literal(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open source prefix");
    PrefixLiteral answer;
    bool saw_pair = false;
    bool saw_prefix = false;
    bool saw_masks = false;
    bool saw_verdict = false;
    std::string line;
    while (std::getline(input, line)) {
        if (line.rfind("PAIR ", 0) == 0) {
            std::istringstream row(line.substr(5));
            char comma = 0;
            row >> answer.q >> comma >> answer.r;
            require(static_cast<bool>(row) && comma == ',', "bad pair row");
            saw_pair = true;
        } else if (line.rfind("PREFIX_CERT ", 0) == 0) {
            std::istringstream row(line);
            std::string tag, size_tag, fnv_tag, fnv_word;
            row >> tag >> size_tag >> answer.stated_size >> fnv_tag >> fnv_word;
            require(static_cast<bool>(row) && tag == "PREFIX_CERT" &&
                        size_tag == "SIZE" && fnv_tag == "FNV",
                    "bad prefix declaration");
            answer.stated_fnv = parse_hex_literal(fnv_word);
            saw_prefix = true;
        } else if (line.rfind("PREFIX_MASKS_HEX ", 0) == 0) {
            const std::string rest = line.substr(17);
            std::size_t begin = 0;
            while (true) {
                const std::size_t comma = rest.find(',', begin);
                const std::string word = rest.substr(
                    begin, comma == std::string::npos ? std::string::npos
                                                      : comma - begin);
                require(!word.empty(), "empty prefix mask");
                const u64 wide = parse_hex_literal(word);
                require(wide < (UINT64_C(1) << 30), "mask leaves pool");
                answer.masks.push_back(static_cast<u32>(wide));
                if (comma == std::string::npos) break;
                begin = comma + 1;
            }
            saw_masks = true;
        } else if (line.rfind("VERDICT EVERY_BODY_CLOSED", 0) == 0) {
            saw_verdict = true;
        }
    }
    require(saw_pair && saw_prefix && saw_masks && saw_verdict,
            "incomplete source prefix");
    require(answer.masks.size() == answer.stated_size,
            "prefix size changed");
    Fnv ledger;
    std::set<u32> distinct;
    for (u32 mask : answer.masks) {
        require(std::popcount(mask) == 8, "source repair rank changed");
        require(distinct.insert(mask).second, "source prefix repeats a mask");
        ledger.add(mask);
    }
    require(ledger.state == answer.stated_fnv, "source prefix FNV changed");
    return answer;
}

std::vector<u32> final_carrier(
    const std::map<std::pair<int, int>, Parsed>& parsed,
    const std::filesystem::path& result_root) {
    std::vector<u32> answer = base_union(parsed);
    std::set<u32> seen(answer.begin(), answer.end());
    struct Source { const char* name; int q; int r; std::size_t size; };
    constexpr std::array<Source, 4> sources = {{
        {"full_discovery_416_704_O3.semantic.out", 416, 704, 4733},
        {"full_discovery_416_700_O3.semantic.out", 416, 700, 7703},
        {"full_discovery_520_700_O3.semantic.out", 520, 700, 7986},
        {"full_discovery_384_694_O3.semantic.out", 384, 694, 8319},
    }};
    for (const Source& source : sources) {
        const PrefixLiteral prefix = read_prefix_literal(result_root / source.name);
        require(prefix.q == source.q && prefix.r == source.r,
                "source pair changed");
        for (u32 mask : prefix.masks) {
            if (seen.insert(mask).second) answer.push_back(mask);
        }
        require(answer.size() == source.size, "carrier union size changed");
    }
    Fnv ledger;
    for (u32 mask : answer) ledger.add(mask);
    require(ledger.state == EXPECTED_FINAL_FNV, "final carrier FNV changed");
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

std::vector<u32> literal_active(const std::vector<u32>& candidates,
                                const LiteralAtoms& atoms) {
    std::vector<u32> active;
    for (u32 repair : candidates) {
        const i64 mass = literal_mass(atoms, repair);
        const i128 margin = static_cast<i128>(63) * mass -
                            static_cast<i128>(4) * atoms.geometry.grid;
        if (margin < 0) continue;
        require(margin > 0, "literal activation equality appeared");
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

void enumerate_allowed(const std::vector<int>& allowed, int index, int need,
                       u32 mask, std::vector<u32>& answer) {
    if (need == 0) {
        answer.push_back(mask);
        return;
    }
    for (int cursor = index;
         cursor <= static_cast<int>(allowed.size()) - need; ++cursor) {
        enumerate_allowed(allowed, cursor + 1, need - 1,
                          mask | (u32{1} << allowed[cursor]), answer);
    }
}

}  // namespace

#ifndef TOP_LITERAL_LIBRARY_ONLY
int main(int argc, char** argv) {
    try {
        require(argc == 3,
                "usage: literal-audit THM4254_REPLAY_DIR THM4266_RESULTS");
        const auto parsed = parse_probe_directory(argv[1]);
        const std::vector<u32> carrier = final_carrier(parsed, argv[2]);

        constexpr int q = 520;
        constexpr int r = 688;
        const LiteralAtoms atoms = build_literal_atoms(q, r);
        const std::vector<u32> active = literal_active(carrier, atoms);
        require(active.size() == 5934, "literal active carrier changed");

        std::vector<u32> bodies;
        bodies.reserve(EXPECTED_BODIES);
        enumerate_bodies(0, 9, 0, bodies);
        require(bodies.size() == EXPECTED_BODIES, "body universe changed");
        u64 inherited_checks = 0;
        std::vector<u32> failures =
            collect_failures(bodies, active, inherited_checks);
        std::sort(failures.begin(), failures.end());
        require(failures == std::vector<u32>({UINT32_C(0x07187008),
                                               UINT32_C(0x27503008)}) &&
                    inherited_checks == UINT64_C(486545452),
                "literal failure ledger disagrees with primary");

        const u32 body_union = failures[0] | failures[1];
        std::vector<int> allowed;
        for (int bit = 0; bit < 30; ++bit) {
            if ((body_union & (u32{1} << bit)) == 0) allowed.push_back(bit);
        }
        require(allowed.size() == 19, "common complement size changed");
        std::vector<u32> candidates;
        candidates.reserve(75582);
        enumerate_allowed(allowed, 0, 8, 0, candidates);
        require(candidates.size() == 75582, "common candidate count changed");

        std::vector<u32> common_active;
        u64 equalities = 0;
        i128 minimum_margin = 0;
        i128 maximum_margin = 0;
        for (u32 repair : candidates) {
            const i64 mass = literal_mass(atoms, repair);
            const i128 margin = static_cast<i128>(63) * mass -
                                static_cast<i128>(4) * atoms.geometry.grid;
            if (margin < 0) continue;
            if (margin == 0) ++equalities;
            if (common_active.empty()) minimum_margin = maximum_margin = margin;
            minimum_margin = std::min(minimum_margin, margin);
            maximum_margin = std::max(maximum_margin, margin);
            common_active.push_back(repair);
        }
        std::sort(common_active.begin(), common_active.end());
        require(common_active.size() == 72 && equalities == 0 &&
                    common_active.front() == EXPECTED_CANONICAL,
                "literal common-active classification disagrees with primary");
        Fnv common_fnv;
        for (u32 repair : common_active) common_fnv.add(repair);
        require(common_fnv.state == EXPECTED_COMMON_FNV,
                "literal common-active FNV disagrees with primary");
        std::set<u32> carrier_set(carrier.begin(), carrier.end());
        const u64 overlap = static_cast<u64>(std::count_if(
            common_active.begin(), common_active.end(),
            [&](u32 repair) { return carrier_set.contains(repair); }));
        require(overlap == 0, "common repair already belongs to carrier");

        const i64 canonical_mass = literal_mass(atoms, EXPECTED_CANONICAL);
        const i128 canonical_margin =
            static_cast<i128>(63) * canonical_mass -
            static_cast<i128>(4) * atoms.geometry.grid;
        require(canonical_margin > 0 &&
                    fraction(canonical_mass, atoms.geometry.grid) ==
                        "1559831620541/24511557965895",
                "literal canonical mass disagrees with primary");
        std::vector<u32> augmented = active;
        augmented.push_back(EXPECTED_CANONICAL);
        const BodyAudit body_audit = audit_bodies(bodies, augmented);
        require(body_audit.failures == 0 &&
                    body_audit.bodies == EXPECTED_BODIES &&
                    body_audit.checks == UINT64_C(486545454) &&
                    body_audit.max_checks == UINT64_C(5935),
                "literal augmented body closure changed");

        std::cout << "LRC14_ENDPOINT_520_688_MINIMUM_ONE_ATOM_LITERAL_WALL_AUDIT_THM4278\n";
        std::cout << "PAIR " << q << ',' << r << " GRID "
                  << atoms.geometry.grid << " JOINT_CELLS "
                  << atoms.geometry.cells.size() << " LITERAL_ATOMS "
                  << atoms.mass.size() << " PAIR_TICKS "
                  << atoms.geometry.pair_ticks << '\n';
        std::cout << "CARRIER " << carrier.size() << " ACTIVE "
                  << active.size() << " BODIES " << bodies.size()
                  << " FAILURES " << failures.size() << " CHECKS "
                  << inherited_checks << '\n';
        for (std::size_t i = 0; i < failures.size(); ++i) {
            std::cout << "FAILURE_" << i << " MASK " << std::hex
                      << failures[i] << std::dec << " LABELS {"
                      << labels(failures[i]) << "}\n";
        }
        std::cout << "COMMON_COMPLEMENT_SIZE " << allowed.size()
                  << " CANDIDATES " << candidates.size() << " ACTIVE "
                  << common_active.size() << " EQUALITIES " << equalities
                  << " OVERLAP " << overlap << " ACTIVE_FNV " << std::hex
                  << common_fnv.state << std::dec << '\n';
        std::cout << "CANONICAL MASK " << std::hex << EXPECTED_CANONICAL
                  << std::dec << " LABELS {" << labels(EXPECTED_CANONICAL)
                  << "} MASS " << fraction(canonical_mass, atoms.geometry.grid)
                  << " MARGIN63_NUM " << decimal(canonical_margin) << '\n';
        std::cout << "MARGIN_RANGE [" << decimal(minimum_margin) << ','
                  << decimal(maximum_margin) << "]\n";
        std::cout << "AUGMENTED_ACTIVE " << augmented.size()
                  << " BODIES " << body_audit.bodies << " FAILURES "
                  << body_audit.failures << " CHECKS " << body_audit.checks
                  << " MAX_CHECKS " << body_audit.max_checks << '\n';
        std::cout << "COMMON_ACTIVE_MASKS_HEX";
        for (u32 mask : common_active) std::cout << ' ' << std::hex << mask;
        std::cout << std::dec << '\n';
        std::cout << "STATUS FINITE_EXACT INDEPENDENT_LITERAL_WALL PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "LITERAL_ERROR " << error.what() << '\n';
        return 1;
    }
}
#endif
