// Exact proof checker for the frozen 59-edge endpoint-cocycle cascade.
//
// The discovery transcripts were produced by cascade_pair_exhaustive_primary.cpp,
// which evaluates all C(30,8) repairs by a labelled superset transform.  This
// checker deliberately does not repeat that transform.  It reads only the
// emitted short prefixes, recomputes each listed repair mass directly from the
// signed endpoint-cocycle atoms, and exhausts all C(30,9) bodies against each
// prefix.  Thus discovery and proof checking have materially different inner
// loops.

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#pragma clang diagnostic ignored "-Wunneeded-internal-declaration"
#endif
#define CASCADE_LIBRARY_ONLY
#include "cascade_pair_exhaustive_primary.cpp"
#undef CASCADE_LIBRARY_ONLY
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>

namespace {

struct PairKey { int q; int r; };

constexpr PairKey BAND[] = {
    {616,755},{616,756},{616,757},{616,758},{616,759},{616,760},{616,761},
    {616,762},{616,763},{616,764},{616,765},{616,766},{616,767},{616,768},
    {698,755},{698,757},
    {704,755},{704,757},{704,758},{704,759},{704,761},{704,762},{704,763},
    {704,764},{704,765},
    {721,755},{721,757},{721,758},{721,759},{721,761},{721,762},{721,763},
    {721,764},{721,765},{721,766},{721,767},{721,768},
    {726,755},{726,757},{726,758},{726,761},
    {732,755},{732,757},{732,761},{732,762},{732,763},
    {744,762},{744,763},{744,765},{744,766},{744,768},
    {750,762},{750,763},{750,765},{750,766},{750,768},
    {765,766},{765,768},{766,768}
};

struct FrozenPrefix {
    int q = 0;
    int r = 0;
    std::size_t stated_size = 0;
    u64 stated_fnv = 0;
    std::vector<u32> masks;
};

u64 parse_hex64(const std::string& word) {
    std::size_t used = 0;
    const u64 value = std::stoull(word, &used, 16);
    require(used == word.size(), "trailing characters in hexadecimal word");
    return value;
}

FrozenPrefix read_prefix(const std::filesystem::path& path, PairKey expected) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open frozen primary transcript");
    FrozenPrefix out;
    out.q = expected.q;
    out.r = expected.r;
    bool saw_pair = false;
    bool saw_header = false;
    bool saw_masks = false;
    bool saw_verdict = false;
    std::string line;
    while (std::getline(input, line)) {
        if (line.rfind("PAIR ", 0) == 0) {
            std::istringstream row(line.substr(5));
            int q = 0, r = 0;
            char comma = 0;
            row >> q >> comma >> r;
            require(q == expected.q && comma == ',' && r == expected.r,
                    "transcript pair mismatch");
            saw_pair = true;
        } else if (line.rfind("PREFIX_CERT ", 0) == 0) {
            std::istringstream row(line);
            std::string tag, size_tag, fnv_tag, fnv_word;
            row >> tag >> size_tag >> out.stated_size >> fnv_tag >> fnv_word;
            require(tag == "PREFIX_CERT" && size_tag == "SIZE" &&
                    fnv_tag == "FNV", "malformed prefix header");
            out.stated_fnv = parse_hex64(fnv_word);
            saw_header = true;
        } else if (line.rfind("PREFIX_MASKS_HEX ", 0) == 0) {
            std::string rest = line.substr(std::string("PREFIX_MASKS_HEX ").size());
            std::size_t begin = 0;
            while (begin <= rest.size()) {
                const std::size_t comma = rest.find(',', begin);
                const std::string word = rest.substr(
                    begin, comma == std::string::npos ? std::string::npos : comma - begin);
                require(!word.empty(), "empty mask in frozen prefix");
                const u64 wide = parse_hex64(word);
                require(wide < (UINT64_C(1) << 30), "prefix mask leaves pool universe");
                out.masks.push_back(static_cast<u32>(wide));
                if (comma == std::string::npos) break;
                begin = comma + 1;
            }
            saw_masks = true;
        } else if (line == "VERDICT EVERY_BODY_CLOSED" ||
                   line.rfind("VERDICT EVERY_BODY_CLOSED ", 0) == 0) {
            saw_verdict = true;
        }
    }
    require(saw_pair && saw_header && saw_masks && saw_verdict,
            "incomplete frozen primary transcript");
    require(out.masks.size() == out.stated_size, "prefix size mismatch");
    std::set<u32> distinct;
    FnvLocal ledger;
    for (u32 mask : out.masks) {
        require(std::popcount(mask) == 8, "prefix repair is not an eight-set");
        require(distinct.insert(mask).second, "duplicate prefix repair");
        ledger.add(mask);
    }
    require(ledger.state == out.stated_fnv, "prefix FNV mismatch");
    return out;
}

i128 direct_atom_mass(const AtomData& atoms, u32 repair) {
    i128 mass = 0;
    for (const auto& [failure_mask, value] : atoms.mass) {
        if ((failure_mask & ~repair) == 0) mass += value;
    }
    return mass;
}

std::string mask_hex(u32 mask) {
    std::ostringstream out;
    out << std::hex << mask;
    return out.str();
}

} // namespace

#ifndef CASCADE_PREFIX_LIBRARY_ONLY
int main(int argc, char** argv) {
    try {
        require(argc == 2, "usage: verifier PACKET_ROOT");
        const std::filesystem::path root(argv[1]);
        const std::filesystem::path replay = root / "replay_band";
        const std::vector<Cell> pool_cells = build_pool_cells();
        require(pool_cells.size() == 7133, "pool cell count changed");
        require(std::size(BAND) == 59, "frozen band cardinality changed");

        std::vector<u32> prefix_union;
        std::set<u32> union_seen;
        FnvLocal pair_summary_ledger;
        FnvLocal pair_labelled_incidence_ledger;
        u64 total_prefix_incidences = 0;
        u64 total_body_cases = 0;
        u64 total_body_checks = 0;
        u64 total_atom_tests = 0;

        std::cout << "ENDPOINT_CASCADE_PREFIX_UNION_EXACT_V1\n";
        for (const PairKey key : BAND) {
            const auto path = replay /
                (std::to_string(key.q) + "_" + std::to_string(key.r) + ".out");
            const FrozenPrefix frozen = read_prefix(path, key);
            const i64 g = std::gcd(key.q, key.r);
            const PrimitivePair primitive = build_primitive(key.q / g, key.r / g);
            const AtomData atoms = build_cocycle_atoms(pool_cells, primitive, g);
            const i128 denominator = static_cast<i128>(primitive.grid) * g * COMMON;
            const Reduced beta = reduced(primitive.safe_ticks, primitive.grid);
            const Reduced omega = reduced(primitive.max_raw - primitive.min_raw,
                                          static_cast<i128>(primitive.grid) * primitive.grid);
            const Reduced omega_over_g = reduced(omega.num, omega.den * g);
            require(static_cast<i128>(91) * beta.num >=
                        static_cast<i128>(66) * beta.den,
                    "beta gate unexpectedly fails in frozen band");
            require(omega_over_g.num * THRESHOLD_DEN >
                        THRESHOLD_NUM * omega_over_g.den,
                    "pair unexpectedly enters scalar gate");

            i128 minimum_margin = 0;
            u32 minimum_repair = 0;
            i128 last_mass = 0;
            for (std::size_t i = 0; i < frozen.masks.size(); ++i) {
                const u32 repair = frozen.masks[i];
                const i128 mass = direct_atom_mass(atoms, repair);
                const i128 margin = static_cast<i128>(63) * mass -
                                    static_cast<i128>(4) * denominator;
                require(margin >= 0, "frozen prefix contains inactive repair");
                if (i == 0 || margin < minimum_margin) {
                    minimum_margin = margin;
                    minimum_repair = repair;
                }
                if (i + 1 == frozen.masks.size()) last_mass = mass;
                total_atom_tests += atoms.mass.size();
                if (union_seen.insert(repair).second) prefix_union.push_back(repair);
            }

            const ComponentReport last_component = component_report(
                pool_cells, primitive, g, frozen.masks.back());
            require(last_component.exact_mass_num == last_mass,
                    "last-repair component sum disagrees with direct atom mass");

            const ScanResult scan = scan_bodies(frozen.masks);
            require(scan.bodies == EXPECTED_BODIES && scan.failures == 0,
                    "frozen prefix does not cover every labelled body");
            require(scan.max_checks == frozen.masks.size(),
                    "prefix lacks an order-minimality witness");
            require((scan.worst_body & frozen.masks.back()) == 0,
                    "minimality witness is not disjoint from final repair");
            for (std::size_t i = 0; i + 1 < frozen.masks.size(); ++i) {
                require((scan.worst_body & frozen.masks[i]) != 0,
                        "minimality witness misses an earlier repair");
            }

            pair_summary_ledger.add(static_cast<u64>(key.q));
            pair_summary_ledger.add(static_cast<u64>(key.r));
            pair_summary_ledger.add(static_cast<u64>(frozen.masks.size()));
            pair_summary_ledger.add(frozen.stated_fnv);
            pair_labelled_incidence_ledger.add(static_cast<u64>(key.q));
            pair_labelled_incidence_ledger.add(static_cast<u64>(key.r));
            pair_labelled_incidence_ledger.add(
                static_cast<u64>(frozen.masks.size()));
            for (u32 repair : frozen.masks) {
                pair_labelled_incidence_ledger.add(repair);
            }
            total_prefix_incidences += frozen.masks.size();
            total_body_cases += scan.bodies;
            total_body_checks += scan.checks;

            std::cout << "EDGE " << key.q << ',' << key.r
                      << " G " << g
                      << " N " << primitive.grid
                      << " BETA " << show(beta)
                      << " OMEGA_OVER_G " << show(omega_over_g)
                      << " SCALAR_GATE 0"
                      << " ATOMS " << atoms.mass.size()
                      << " PREFIX " << frozen.masks.size()
                      << " PREFIX_FNV " << std::hex << frozen.stated_fnv << std::dec
                      << " MIN_MARGIN_NUM " << decimal(minimum_margin)
                      << " MIN_MARGIN_REPAIR " << mask_hex(minimum_repair)
                      << " DEN " << decimal(denominator)
                      << " BODIES " << scan.bodies
                      << " CHECKS " << scan.checks
                      << " MAX_CHECKS " << scan.max_checks
                      << " MINIMALITY_BODY " << mask_hex(scan.worst_body)
                      << " LAST_REPAIR " << mask_hex(frozen.masks.back())
                      << " DIRECT_CELL_AUDIT PASS COMPONENT_AUDIT PASS\n";
        }

        FnvLocal union_ledger;
        for (u32 mask : prefix_union) union_ledger.add(mask);
        std::cout << "SUMMARY EDGES " << std::size(BAND)
                  << " REPAIRS_FULL_PER_EDGE " << EXPECTED_REPAIRS
                  << " BODIES_PER_EDGE " << EXPECTED_BODIES
                  << " PREFIX_INCIDENCES " << total_prefix_incidences
                  << " PREFIX_UNION " << prefix_union.size()
                  << " PREFIX_UNION_FNV " << std::hex << union_ledger.state << std::dec
                  << " PAIR_SUMMARY_FNV " << std::hex
                  << pair_summary_ledger.state
                  << " PAIR_LABELLED_PREFIX_INCIDENCE_FNV "
                  << pair_labelled_incidence_ledger.state << std::dec
                  << " DIRECT_ATOM_TESTS " << total_atom_tests
                  << " BODY_CASES " << total_body_cases
                  << " BODY_CHECKS " << total_body_checks << '\n';
        std::cout << "VERDICT PASS ALL_59_EVERY_BODY PREFIX_UNION_NO_ZETA\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
#endif
