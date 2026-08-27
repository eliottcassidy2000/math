// Second and third CEGAR rounds: enrich the 4,733-mask carrier with masks
// novel in the order-minimal full prefixes for (416,700) and (520,700), then
// optionally enrich the resulting 7,986-mask carrier from (384,694).

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#endif
#define CARRIER_CEGAR_LIBRARY_ONLY
#include "carrier_cegar_descent.cpp"
#undef CARRIER_CEGAR_LIBRARY_ONLY
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

int main(int argc, char** argv) {
    try {
        require(argc == 8 || argc == 9,
                "usage: round2 PACKET_ROOT PREFIX_416_704 PREFIX_416_700 "
                "PREFIX_520_700 RESIDUAL_CSV START_ENDPOINT MIN_ENDPOINT "
                "[PREFIX_384_694]");
        const std::filesystem::path packet(argv[1]);
        const PrefixCegar prefix_704 = read_prefix_cegar(argv[2]);
        const PrefixCegar prefix_416 = read_prefix_cegar(argv[3]);
        const PrefixCegar prefix_520 = read_prefix_cegar(argv[4]);
        require(prefix_704.q == 416 && prefix_704.r == 704 &&
                    prefix_704.masks.size() == 2608 &&
                    prefix_704.stated_fnv == UINT64_C(0x18ff663a123e684e),
                "round-one prefix changed");
        require(prefix_416.q == 416 && prefix_416.r == 700 &&
                    prefix_416.masks.size() == 5894 &&
                    prefix_416.stated_fnv == UINT64_C(0x701c0233c8f8abeb),
                "416,700 prefix changed");
        require(prefix_520.q == 520 && prefix_520.r == 700 &&
                    prefix_520.masks.size() == 4557 &&
                    prefix_520.stated_fnv == UINT64_C(0xd8466558bcd8ef9d),
                "520,700 prefix changed");
        const int start_endpoint = std::stoi(argv[6]);
        const int min_endpoint = std::stoi(argv[7]);
        require(start_endpoint >= min_endpoint, "bad endpoint range");

        std::vector<u32> base;
        std::set<u32> base_seen;
        for (const EdgeCegar& expected : FROZEN_BAND_CEGAR) {
            const PrefixCegar prefix = read_prefix_cegar(
                packet / "replay_band" /
                (std::to_string(expected.q) + "_" +
                 std::to_string(expected.r) + ".out"));
            require(prefix.q == expected.q && prefix.r == expected.r,
                    "band transcript mismatch");
            for (u32 mask : prefix.masks) {
                if (base_seen.insert(mask).second) base.push_back(mask);
            }
        }
        FnvLocal base_ledger;
        for (u32 mask : base) base_ledger.add(mask);
        require(base.size() == 4675 &&
                    base_ledger.state == UINT64_C(0xce4e76ec11df057c),
                "base carrier changed");

        std::vector<u32> round1 = base;
        std::set<u32> round1_seen = base_seen;
        std::vector<u32> novel_704;
        for (u32 mask : prefix_704.masks) {
            if (round1_seen.insert(mask).second) {
                novel_704.push_back(mask);
                round1.push_back(mask);
            }
        }
        FnvLocal round1_ledger;
        for (u32 mask : round1) round1_ledger.add(mask);
        require(novel_704.size() == 58 && round1.size() == 4733 &&
                    round1_ledger.state == UINT64_C(0xa7b046289655c733),
                "round-one carrier changed");

        std::vector<u32> novel_416;
        std::vector<u32> novel_520;
        std::set<u32> novel_416_set;
        std::set<u32> novel_520_set;
        for (u32 mask : prefix_416.masks) {
            if (!round1_seen.count(mask)) {
                novel_416.push_back(mask);
                novel_416_set.insert(mask);
            }
        }
        for (u32 mask : prefix_520.masks) {
            if (!round1_seen.count(mask)) {
                novel_520.push_back(mask);
                novel_520_set.insert(mask);
            }
        }
        require(novel_416.size() == 2970 && novel_520.size() == 1194,
                "round-two novel counts changed");
        std::vector<u32> overlap;
        for (u32 mask : novel_416) {
            if (novel_520_set.count(mask)) overlap.push_back(mask);
        }
        require(overlap.size() == 911, "round-two overlap changed");

        FnvLocal novel_416_ledger;
        FnvLocal novel_520_ledger;
        FnvLocal overlap_ledger;
        for (u32 mask : novel_416) novel_416_ledger.add(mask);
        for (u32 mask : novel_520) novel_520_ledger.add(mask);
        for (u32 mask : overlap) overlap_ledger.add(mask);

        std::vector<u32> round2 = round1;
        std::set<u32> round2_seen = round1_seen;
        for (u32 mask : novel_416) {
            require(round2_seen.insert(mask).second,
                    "416-novel prefix repeats itself");
            round2.push_back(mask);
        }
        u64 appended_520_only = 0;
        for (u32 mask : novel_520) {
            if (round2_seen.insert(mask).second) {
                round2.push_back(mask);
                ++appended_520_only;
            }
        }
        require(appended_520_only == 283 && round2.size() == 7986,
                "round-two union count changed");
        FnvLocal round2_ledger;
        for (u32 mask : round2) round2_ledger.add(mask);

        const std::vector<Cell> pool_cells = build_pool_cells();
        require(pool_cells.size() == 7133, "pool-cell count changed");
        const ScanResult prefix_416_scan = scan_bodies(prefix_416.masks);
        const ScanResult prefix_520_scan = scan_bodies(prefix_520.masks);
        require(prefix_416_scan.failures == 0 &&
                    prefix_416_scan.max_checks == prefix_416.masks.size() &&
                    prefix_520_scan.failures == 0 &&
                    prefix_520_scan.max_checks == prefix_520.masks.size(),
                "round-two source prefix lost order minimality");
        const CarrierScanCegar old_416 = scan_carrier_cegar(
            pool_cells, 416, 700, round1);
        const CarrierScanCegar old_520 = scan_carrier_cegar(
            pool_cells, 520, 700, round1);
        require(old_416.body.failures == 38 &&
                    old_416.body.first_failure == UINT32_C(0x051a7400) &&
                    old_520.body.failures == 1 &&
                    old_520.body.first_failure == UINT32_C(0x2d183001),
                "round-two hostiles changed");
        const CarrierScanCegar new_416 = scan_carrier_cegar(
            pool_cells, 416, 700, round2);
        const CarrierScanCegar new_520 = scan_carrier_cegar(
            pool_cells, 520, 700, round2);
        require(new_416.body.failures == 0 && new_520.body.failures == 0,
                "round-two union does not close both hostiles");

        FnvLocal atom_component_ledger;
        for (std::size_t index = 4733; index < round2.size(); ++index) {
            const u32 mask = round2[index];
            for (const EdgeCegar edge :
                 std::array<EdgeCegar, 2>{{{416,700},{520,700}}}) {
                const i64 g = std::gcd(edge.q, edge.r);
                const PrimitivePair primitive = build_primitive(
                    edge.q / g, edge.r / g);
                const AtomData atoms = build_cocycle_atoms(
                    pool_cells, primitive, g);
                const i128 mass = atom_mass_cegar(atoms, mask);
                const ComponentReport component = component_report(
                    pool_cells, primitive, g, mask);
                require(component.exact_mass_num == mass,
                        "round-two atom/component mismatch");
                atom_component_ledger.add(mask);
                atom_component_ledger.add(static_cast<u64>(edge.q));
                atom_component_ledger.add(static_cast<u64>(edge.r));
                fnv_add_i128_cegar(atom_component_ledger, mass);
            }
        }

        std::vector<u32> descent_carrier = round2;
        bool has_round3 = argc == 9;
        u64 round3_novel_count = 0;
        u64 round3_prefix_size = 0;
        u64 round3_prefix_fnv = 0;
        u64 round3_carrier_fnv = 0;
        u64 round3_atom_component_fnv = 0;
        u64 round3_old_failures = 0;
        u64 round3_new_failures = 0;
        u64 round3_new_active = 0;
        if (has_round3) {
            const PrefixCegar prefix_384 = read_prefix_cegar(argv[8]);
            require(prefix_384.q == 384 && prefix_384.r == 694 &&
                        prefix_384.masks.size() == 5307 &&
                        prefix_384.stated_fnv ==
                            UINT64_C(0x7b9a46c60e8514f8),
                    "384,694 prefix changed");
            const ScanResult prefix_384_scan = scan_bodies(prefix_384.masks);
            require(prefix_384_scan.failures == 0 &&
                        prefix_384_scan.max_checks == prefix_384.masks.size(),
                    "384,694 prefix lost order minimality");
            std::vector<u32> novel_384;
            for (u32 mask : prefix_384.masks) {
                if (round2_seen.insert(mask).second) {
                    novel_384.push_back(mask);
                    descent_carrier.push_back(mask);
                }
            }
            require(novel_384.size() == 333 &&
                        descent_carrier.size() == 8319,
                    "round-three novel count changed");
            FnvLocal novel_384_ledger;
            for (u32 mask : novel_384) novel_384_ledger.add(mask);
            require(novel_384_ledger.state ==
                        UINT64_C(0x1d1422eaa770226d),
                    "round-three novel FNV changed");
            FnvLocal round3_ledger;
            for (u32 mask : descent_carrier) round3_ledger.add(mask);
            FnvLocal round3_component_ledger;
            const i64 g384 = std::gcd(384, 694);
            const PrimitivePair primitive384 = build_primitive(
                384 / g384, 694 / g384);
            const AtomData atoms384 = build_cocycle_atoms(
                pool_cells, primitive384, g384);
            for (u32 mask : novel_384) {
                const i128 mass = atom_mass_cegar(atoms384, mask);
                const ComponentReport component = component_report(
                    pool_cells, primitive384, g384, mask);
                require(component.exact_mass_num == mass,
                        "round-three atom/component mismatch");
                round3_component_ledger.add(mask);
                fnv_add_i128_cegar(round3_component_ledger, mass);
            }
            const CarrierScanCegar old_384 = scan_carrier_cegar(
                pool_cells, 384, 694, round2);
            const CarrierScanCegar new_384 = scan_carrier_cegar(
                pool_cells, 384, 694, descent_carrier);
            require(old_384.body.failures == 1 &&
                        old_384.body.first_failure == UINT32_C(0x0d186401) &&
                        new_384.body.failures == 0,
                    "round-three seed closure changed");
            round3_novel_count = novel_384.size();
            round3_prefix_size = prefix_384.masks.size();
            round3_prefix_fnv = prefix_384.stated_fnv;
            round3_carrier_fnv = round3_ledger.state;
            round3_atom_component_fnv = round3_component_ledger.state;
            round3_old_failures = old_384.body.failures;
            round3_new_failures = new_384.body.failures;
            round3_new_active = new_384.active;
        }

        std::cout << (has_round3 ? "CARRIER_CEGAR_ROUND23_V1\n"
                                 : "CARRIER_CEGAR_ROUND2_V1\n");
        std::cout << "ROUND1 UNION " << round1.size() << " FNV " << std::hex
                  << round1_ledger.state << std::dec << "\n";
        std::cout << "PREFIX_416_700 SIZE " << prefix_416.masks.size()
                  << " FNV " << std::hex << prefix_416.stated_fnv << std::dec
                  << " OVERLAP_ROUND1 "
                  << prefix_416.masks.size() - novel_416.size()
                  << " NOVEL " << novel_416.size() << " NOVEL_FNV "
                  << std::hex << novel_416_ledger.state << std::dec << "\n";
        std::cout << "PREFIX_520_700 SIZE " << prefix_520.masks.size()
                  << " FNV " << std::hex << prefix_520.stated_fnv << std::dec
                  << " OVERLAP_ROUND1 "
                  << prefix_520.masks.size() - novel_520.size()
                  << " NOVEL " << novel_520.size() << " NOVEL_FNV "
                  << std::hex << novel_520_ledger.state << std::dec << "\n";
        std::cout << "NOVEL_OVERLAP " << overlap.size() << " FNV "
                  << std::hex << overlap_ledger.state << std::dec
                  << " ONLY_416 " << novel_416.size() - overlap.size()
                  << " ONLY_520 " << novel_520.size() - overlap.size()
                  << " UNION_NOVEL " << round2.size() - round1.size()
                  << "\n";
        std::cout << "ROUND2 UNION " << round2.size() << " FNV " << std::hex
                  << round2_ledger.state << " ATOM_COMPONENT_LEDGER "
                  << atom_component_ledger.state << std::dec << "\n";
        std::cout << "HOSTILE_416_700 OLD_FAILURES "
                  << old_416.body.failures << " NEW_FAILURES "
                  << new_416.body.failures << " NEW_ACTIVE " << new_416.active
                  << "\n";
        std::cout << "HOSTILE_520_700 OLD_FAILURES "
                  << old_520.body.failures << " NEW_FAILURES "
                  << new_520.body.failures << " NEW_ACTIVE " << new_520.active
                  << "\n";
        if (has_round3) {
            std::cout << "ROUND3 PREFIX_384_694 SIZE "
                      << round3_prefix_size << " FNV " << std::hex
                      << round3_prefix_fnv << std::dec
                      << " OVERLAP_ROUND2 "
                      << round3_prefix_size - round3_novel_count
                      << " NOVEL " << round3_novel_count
                      << " UNION " << descent_carrier.size()
                      << " UNION_FNV " << std::hex << round3_carrier_fnv
                      << " ATOM_COMPONENT_LEDGER "
                      << round3_atom_component_fnv << std::dec << "\n";
            std::cout << "HOSTILE_384_694 OLD_FAILURES "
                      << round3_old_failures << " NEW_FAILURES "
                      << round3_new_failures << " NEW_ACTIVE "
                      << round3_new_active << "\n";
        }

        const std::vector<EdgeCegar> all_edges = read_edges_cegar(argv[5]);
        int current_layer = std::numeric_limits<int>::max();
        u64 layer_rows = 0;
        u64 layer_resistant = 0;
        u64 layer_active_sum = 0;
        u64 layer_checks = 0;
        u64 layer_max_checks = 0;
        FnvLocal layer_ledger;
        std::vector<std::string> resistant_lines;
        bool found_layer = false;

        auto flush_layer = [&]() {
            if (current_layer == std::numeric_limits<int>::max()) return false;
            std::cout << "LAYER " << current_layer << " ROWS " << layer_rows
                      << " RESISTANT " << layer_resistant
                      << " ACTIVE_SUM " << layer_active_sum
                      << " CHECKS " << layer_checks
                      << " MAX_CHECKS " << layer_max_checks
                      << " ROW_LEDGER " << std::hex << layer_ledger.state
                      << std::dec << "\n";
            for (const std::string& line : resistant_lines) {
                std::cout << line << '\n';
            }
            return layer_resistant != 0;
        };

        for (const EdgeCegar& edge : all_edges) {
            if (edge.r > start_endpoint || edge.r < min_endpoint) continue;
            if (current_layer != edge.r) {
                if (flush_layer()) {
                    found_layer = true;
                    break;
                }
                current_layer = edge.r;
                layer_rows = layer_resistant = layer_active_sum = 0;
                layer_checks = layer_max_checks = 0;
                layer_ledger = FnvLocal{};
                resistant_lines.clear();
            }
            const CarrierScanCegar scan = scan_carrier_cegar(
                pool_cells, edge.q, edge.r, descent_carrier);
            ++layer_rows;
            if (scan.body.failures != 0) ++layer_resistant;
            layer_active_sum += scan.active;
            layer_checks += scan.body.checks;
            layer_max_checks = std::max(layer_max_checks,
                                        scan.body.max_checks);
            layer_ledger.add(static_cast<u64>(edge.q));
            layer_ledger.add(static_cast<u64>(edge.r));
            layer_ledger.add(scan.active);
            layer_ledger.add(scan.body.failures);
            layer_ledger.add(scan.body.checks);
            layer_ledger.add(scan.body.max_checks);
            layer_ledger.add(scan.body.first_failure);
            layer_ledger.add(scan.best_disjoint_repair);
            fnv_add_i128_cegar(layer_ledger, scan.best_disjoint_margin);
            fnv_add_i128_cegar(layer_ledger, scan.denominator);
            if (scan.body.failures != 0) {
                std::ostringstream row;
                row << "RESISTANT PAIR " << edge.q << ',' << edge.r
                    << " ACTIVE " << scan.active
                    << " FAILURES " << scan.body.failures
                    << " FIRST_FAILURE " << std::hex
                    << scan.body.first_failure << std::dec
                    << " FIRST_FAILURE_LABELS {"
                    << labels(scan.body.first_failure) << "}"
                    << " BEST_DISJOINT_REPAIR " << std::hex
                    << scan.best_disjoint_repair << std::dec
                    << " BEST_DISJOINT_LABELS {"
                    << labels(scan.best_disjoint_repair) << "}"
                    << " BEST_DISJOINT_MARGIN_NUM "
                    << decimal(scan.best_disjoint_margin)
                    << " DEN " << decimal(scan.denominator);
                resistant_lines.push_back(row.str());
            }
        }
        if (!found_layer) found_layer = flush_layer();
        std::cout << "STOP FIRST_RESISTANT_LAYER "
                  << (found_layer ? current_layer : -1)
                  << " VERDICT " << (found_layer ? "BOUNDARY_FOUND"
                                                   : "RANGE_CLOSED")
                  << "\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "CARRIER_CEGAR_ROUND2_ERROR " << error.what() << '\n';
        return 1;
    }
}
