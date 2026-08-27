// Hostile boundary probe for the frozen 4,675-repair prefix union.
// This tests resistance of the compressed certificate, not of the full
// 5,852,925-repair deck.  A failure is therefore not a counterexample to the
// endpoint-cocycle method; it is a precise stopping witness for this packet's
// reusable union.

#define CASCADE_PREFIX_LIBRARY_ONLY
#include "cascade_prefix_union_exact_verifier.cpp"
#undef CASCADE_PREFIX_LIBRARY_ONLY

int main(int argc, char** argv) {
    try {
        require(argc == 4, "usage: boundary_probe PACKET_ROOT q r");
        const std::filesystem::path root(argv[1]);
        int q = std::stoi(argv[2]);
        int r = std::stoi(argv[3]);
        if (q > r) std::swap(q, r);
        require(q > 0 && q < r, "pair ordering invalid");

        std::vector<u32> prefix_union;
        std::set<u32> seen;
        for (const PairKey key : BAND) {
            const auto path = root / "replay_band" /
                (std::to_string(key.q) + "_" + std::to_string(key.r) + ".out");
            const FrozenPrefix frozen = read_prefix(path, key);
            for (u32 repair : frozen.masks) {
                if (seen.insert(repair).second) prefix_union.push_back(repair);
            }
        }
        require(prefix_union.size() == 4675, "prefix union changed");

        const std::vector<Cell> pool_cells = build_pool_cells();
        const i64 g = std::gcd(q, r);
        const PrimitivePair primitive = build_primitive(q / g, r / g);
        const AtomData atoms = build_cocycle_atoms(pool_cells, primitive, g);
        const i128 denominator = static_cast<i128>(primitive.grid) * g * COMMON;
        std::vector<u32> active;
        active.reserve(prefix_union.size());
        std::vector<i128> margins;
        margins.reserve(prefix_union.size());
        for (u32 repair : prefix_union) {
            const i128 mass = direct_atom_mass(atoms, repair);
            const i128 margin = static_cast<i128>(63) * mass -
                                static_cast<i128>(4) * denominator;
            margins.push_back(margin);
            if (margin >= 0) active.push_back(repair);
        }
        const ScanResult scan = scan_bodies(active);

        i128 best_disjoint_margin = 0;
        u32 best_disjoint_repair = 0;
        if (scan.failures != 0) {
            bool first = true;
            for (std::size_t i = 0; i < prefix_union.size(); ++i) {
                if ((prefix_union[i] & scan.first_failure) != 0) continue;
                if (first || margins[i] > best_disjoint_margin) {
                    first = false;
                    best_disjoint_margin = margins[i];
                    best_disjoint_repair = prefix_union[i];
                }
            }
            require(best_disjoint_repair != 0 && best_disjoint_margin < 0,
                    "failure witness has an active disjoint union repair");
        }

        std::cout << "CASCADE_UNION_BOUNDARY_V1 PAIR " << q << ',' << r
                  << " UNION " << prefix_union.size()
                  << " ACTIVE_UNION " << active.size()
                  << " BODIES " << scan.bodies
                  << " FAILURES " << scan.failures
                  << " CHECKS " << scan.checks
                  << " MAX_CHECKS " << scan.max_checks
                  << " FIRST_FAILURE " << mask_hex(
                         scan.failures ? scan.first_failure : u32{0})
                  << " BEST_DISJOINT_REPAIR " << mask_hex(best_disjoint_repair)
                  << " BEST_DISJOINT_MARGIN_NUM " << decimal(best_disjoint_margin)
                  << " DEN " << decimal(denominator)
                  << " VERDICT " << (scan.failures ? "RESISTANT" : "CLOSED")
                  << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
