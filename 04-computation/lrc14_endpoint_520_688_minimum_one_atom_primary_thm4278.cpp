// Exact primary audit for THM-4278: minimum rank-eight repair augmentation of
// THM-4266's frozen 8,319-mask carrier at the former top pair (520,688).

#define CARRIER_CEGAR_LIBRARY_ONLY
#include "lrc14_three_round_learned_carrier_thm4266/carrier_cegar_descent.cpp"
#undef CARRIER_CEGAR_LIBRARY_ONLY

#include <iomanip>

namespace {

std::vector<u32> reconstruct_round3_carrier(
    const std::filesystem::path& packet,
    const std::filesystem::path& result_root) {
    std::vector<u32> carrier;
    std::set<u32> seen;
    for (const EdgeCegar expected : FROZEN_BAND_CEGAR) {
        const PrefixCegar prefix = read_prefix_cegar(
            packet / "replay_band" /
            (std::to_string(expected.q) + "_" +
             std::to_string(expected.r) + ".out"));
        require(prefix.q == expected.q && prefix.r == expected.r,
                "base prefix pair changed");
        for (u32 mask : prefix.masks) {
            if (seen.insert(mask).second) carrier.push_back(mask);
        }
    }
    FnvLocal ledger;
    for (u32 mask : carrier) ledger.add(mask);
    require(carrier.size() == 4675 &&
                ledger.state == UINT64_C(0xce4e76ec11df057c),
            "base carrier changed");

    struct Source {
        const char* name;
        int q;
        int r;
        std::size_t final_size;
        u64 final_fnv;
    };
    constexpr std::array<Source, 4> sources = {{
        {"full_discovery_416_704_O3.semantic.out", 416, 704, 4733,
         UINT64_C(0xa7b046289655c733)},
        {"full_discovery_416_700_O3.semantic.out", 416, 700, 7703, 0},
        {"full_discovery_520_700_O3.semantic.out", 520, 700, 7986,
         UINT64_C(0xbaef1d2f49444638)},
        {"full_discovery_384_694_O3.semantic.out", 384, 694, 8319,
         UINT64_C(0xe08b227730f6793c)},
    }};
    for (const Source& source : sources) {
        const PrefixCegar prefix = read_prefix_cegar(result_root / source.name);
        require(prefix.q == source.q && prefix.r == source.r,
                "source prefix pair changed");
        for (u32 mask : prefix.masks) {
            if (seen.insert(mask).second) carrier.push_back(mask);
        }
        require(carrier.size() == source.final_size,
                "source prefix union size changed");
        if (source.final_fnv != 0) {
            FnvLocal current;
            for (u32 mask : carrier) current.add(mask);
            require(current.state == source.final_fnv,
                    "source prefix union FNV changed");
        }
    }
    return carrier;
}

std::vector<u32> list_failures(const std::vector<u32>& active,
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
    require(bodies == EXPECTED_BODIES, "body universe changed");
    return failures;
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

std::string hex_local(u32 mask) {
    std::ostringstream out;
    out << std::hex << mask;
    return out.str();
}

struct DisjointLedger {
    u64 candidates = 0;
    u64 active = 0;
    u64 equalities = 0;
    u32 least_active = 0;
    i128 least_mass = 0;
    i128 least_margin = 0;
    u32 minimum_margin_mask = 0;
    i128 minimum_margin = 0;
    u32 maximum_margin_mask = 0;
    i128 maximum_margin = 0;
    FnvLocal active_fnv;
    std::vector<u32> active_masks;
};

}  // namespace

#ifndef TOP_PAIR_MIN_AUG_LIBRARY_ONLY
int main(int argc, char** argv) {
    try {
        require(argc == 3,
                "usage: top-pair-min-augmentation PACKET_ROOT THM4266_RESULTS");
        init_choose8_local();
        const std::filesystem::path packet(argv[1]);
        const std::filesystem::path result_root(argv[2]);
        const std::vector<u32> carrier =
            reconstruct_round3_carrier(packet, result_root);
        const std::vector<Cell> pool_cells = build_pool_cells();
        require(pool_cells.size() == 7133, "pool cells changed");

        constexpr int q = 520;
        constexpr int r = 688;
        const i64 g = std::gcd(q, r);
        const PrimitivePair primitive = build_primitive(q / g, r / g);
        const AtomData atoms = build_cocycle_atoms(pool_cells, primitive, g);
        std::vector<i128> masses(EXPECTED_REPAIRS, 0);
        u64 zeta_operations = 0;
        for (const auto& [mask, value] : atoms.mass) {
            add_supersets_pair(mask, 8 - std::popcount(mask), 0, 0,
                               value, masses, zeta_operations);
        }
        require(zeta_operations == UINT64_C(152170690),
                "zeta operation count changed");
        const i128 denominator =
            static_cast<i128>(primitive.grid) * g * COMMON;

        std::vector<u32> active_carrier;
        for (u32 mask : carrier) {
            const i128 margin = static_cast<i128>(63) *
                                    masses[colex_rank8_local(mask)] -
                                static_cast<i128>(4) * denominator;
            if (margin >= 0) active_carrier.push_back(mask);
        }
        require(active_carrier.size() == 5934,
                "top-pair active carrier changed");
        u64 inherited_checks = 0;
        const std::vector<u32> failures =
            list_failures(active_carrier, inherited_checks);
        require(failures.size() == 2, "top-pair failure count changed");

        const u32 failure_union = failures[0] | failures[1];
        const u32 failure_intersection = failures[0] & failures[1];
        const u32 pool_mask = (u32{1} << 30) - 1;
        const u32 common_complement = pool_mask & ~failure_union;

        std::array<DisjointLedger, 3> ledgers;
        u32 repair = (u32{1} << 8) - 1;
        const u32 limit = u32{1} << 30;
        u64 rank = 0;
        u64 full_active = 0;
        FnvLocal full_active_fnv;
        while (repair != 0 && repair < limit) {
            const i128 mass = masses[rank];
            const i128 margin = static_cast<i128>(63) * mass -
                                static_cast<i128>(4) * denominator;
            const bool active = margin >= 0;
            if (active) {
                ++full_active;
                full_active_fnv.add(repair);
            }
            const std::array<u32, 3> forbidden = {
                failures[0], failures[1], failure_union};
            for (std::size_t i = 0; i < forbidden.size(); ++i) {
                if ((repair & forbidden[i]) != 0) continue;
                DisjointLedger& ledger = ledgers[i];
                ++ledger.candidates;
                if (!active) continue;
                ++ledger.active;
                ledger.active_fnv.add(repair);
                ledger.active_masks.push_back(repair);
                if (margin == 0) ++ledger.equalities;
                if (ledger.least_active == 0) {
                    ledger.least_active = repair;
                    ledger.least_mass = mass;
                    ledger.least_margin = margin;
                    ledger.minimum_margin_mask = repair;
                    ledger.minimum_margin = margin;
                    ledger.maximum_margin_mask = repair;
                    ledger.maximum_margin = margin;
                }
                if (margin < ledger.minimum_margin ||
                    (margin == ledger.minimum_margin &&
                     repair < ledger.minimum_margin_mask)) {
                    ledger.minimum_margin = margin;
                    ledger.minimum_margin_mask = repair;
                }
                if (margin > ledger.maximum_margin ||
                    (margin == ledger.maximum_margin &&
                     repair < ledger.maximum_margin_mask)) {
                    ledger.maximum_margin = margin;
                    ledger.maximum_margin_mask = repair;
                }
            }
            ++rank;
            const u32 next = next_combination(repair);
            if (next <= repair) break;
            repair = next;
        }
        require(rank == EXPECTED_REPAIRS, "repair enumeration incomplete");
        require(full_active == UINT64_C(2504029),
                "full active count changed");

        const DisjointLedger& common = ledgers[2];
        require(common.active > 0, "no one-atom common repair exists");
        require(std::find(carrier.begin(), carrier.end(), common.least_active) ==
                    carrier.end(),
                "least common repair unexpectedly in inherited carrier");
        std::vector<u32> augmented = active_carrier;
        augmented.push_back(common.least_active);
        u64 augmented_checks = 0;
        const std::vector<u32> augmented_failures =
            list_failures(augmented, augmented_checks);
        require(augmented_failures.empty(),
                "least one-atom augmentation did not close every body");

        FnvLocal carrier_fnv;
        for (u32 mask : carrier) carrier_fnv.add(mask);
        FnvLocal failure_fnv;
        for (u32 body : failures) failure_fnv.add(body);
        const Reduced least_fraction =
            reduced(common.least_mass, denominator);

        std::cout << "LRC14_ENDPOINT_520_688_MINIMUM_ONE_ATOM_PRIMARY_THM4278\n";
        std::cout << "PAIR " << q << ',' << r << " PRIMITIVE "
                  << primitive.u << ':' << primitive.v << " SCALE " << g
                  << " GRID " << primitive.grid << " POOL_CELLS "
                  << pool_cells.size() << '\n';
        std::cout << "CARRIER " << carrier.size() << " CARRIER_FNV "
                  << std::hex << carrier_fnv.state << std::dec << " ACTIVE "
                  << active_carrier.size() << " FAILURES " << failures.size()
                  << " CHECKS " << inherited_checks << " FAILURE_FNV "
                  << std::hex << failure_fnv.state << std::dec << '\n';
        for (std::size_t i = 0; i < failures.size(); ++i) {
            std::cout << "FAILURE_" << i << " MASK "
                      << hex_local(failures[i]) << " LABELS {"
                      << labels_local(failures[i]) << "}\n";
        }
        std::cout << "INTERSECTION SIZE " << std::popcount(failure_intersection)
                  << " MASK " << hex_local(failure_intersection)
                  << " LABELS {" << labels_local(failure_intersection) << "}"
                  << " UNION SIZE " << std::popcount(failure_union)
                  << " MASK " << hex_local(failure_union) << " LABELS {"
                  << labels_local(failure_union) << "}\n";
        std::cout << "COMMON_COMPLEMENT SIZE "
                  << std::popcount(common_complement) << " MASK "
                  << hex_local(common_complement) << " LABELS {"
                  << labels_local(common_complement) << "}\n";
        std::cout << "FULL_ACTIVE " << full_active << " FULL_ACTIVE_FNV "
                  << std::hex << full_active_fnv.state << std::dec
                  << " EQUALITIES 0 ZETA_OPS " << zeta_operations << '\n';
        const std::array<const char*, 3> names = {
            "FAILURE_0", "FAILURE_1", "COMMON"};
        for (std::size_t i = 0; i < ledgers.size(); ++i) {
            const DisjointLedger& x = ledgers[i];
            std::cout << "DISJOINT_" << names[i]
                      << " CANDIDATES " << x.candidates
                      << " ACTIVE " << x.active << " EQUALITIES "
                      << x.equalities << " ACTIVE_FNV " << std::hex
                      << x.active_fnv.state << std::dec << " LEAST "
                      << hex_local(x.least_active) << " LABELS {"
                      << labels_local(x.least_active) << "} MIN_MARGIN "
                      << decimal(x.minimum_margin) << " MIN_MASK "
                      << hex_local(x.minimum_margin_mask) << " MAX_MARGIN "
                      << decimal(x.maximum_margin) << " MAX_MASK "
                      << hex_local(x.maximum_margin_mask) << '\n';
        }
        std::cout << "CANONICAL_ONE_ATOM MASK "
                  << hex_local(common.least_active) << " COLEX_RANK0 "
                  << colex_rank8_local(common.least_active) << " LABELS {"
                  << labels_local(common.least_active) << "} MASS "
                  << show(least_fraction) << " MARGIN63_NUM "
                  << decimal(common.least_margin) << '\n';
        std::cout << "COMMON_ACTIVE_MASKS_HEX";
        for (u32 mask : common.active_masks) {
            std::cout << ' ' << std::hex << mask;
        }
        std::cout << std::dec << '\n';
        std::cout << "AUGMENTATION_NUMBER 1 ZERO_FAILS " << failures.size()
                  << " ONE_FAILS " << augmented_failures.size()
                  << " AUGMENTED_ACTIVE " << augmented.size()
                  << " CHECKS " << augmented_checks << '\n';
        std::cout << "STATUS FINITE_EXACT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
#endif
