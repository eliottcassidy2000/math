// Exact primary endpoint scan after appending the order-frozen two-mask
// minimum repair of the sole endpoint-669 hostile row to the THM-4276 carrier
// plus the 421-mask rectangle-common/r670 joint deck.

#define ENDPOINT669_SCAN_LIBRARY_ONLY
#include "scan_augmented_carrier_endpoint669.cpp"
#undef ENDPOINT669_SCAN_LIBRARY_ONLY

namespace {

#ifndef TARGET_ENDPOINT
#define TARGET_ENDPOINT 668
#endif
#ifndef INCLUDE_ENDPOINT667_MINIMUM
#define INCLUDE_ENDPOINT667_MINIMUM 0
#endif
#ifndef INCLUDE_ENDPOINT665_MINIMUM
#define INCLUDE_ENDPOINT665_MINIMUM 0
#endif
#ifndef INCLUDE_ENDPOINT664_MINIMUM
#define INCLUDE_ENDPOINT664_MINIMUM 0
#endif
#ifndef QUIET_ENDPOINT_PROGRESS
#define QUIET_ENDPOINT_PROGRESS 0
#endif
constexpr int SCAN_ENDPOINT = TARGET_ENDPOINT;

constexpr std::array<u32, 2> ENDPOINT669_MINIMUM = {
    UINT32_C(0x003884c8), UINT32_C(0x00c4c124)};
constexpr std::array<u32, 2> ENDPOINT667_MINIMUM = {
    UINT32_C(0x02a05206), UINT32_C(0x10e05240)};
constexpr std::array<u32, 1> ENDPOINT665_MINIMUM = {
    UINT32_C(0x0004409f)};
constexpr std::array<u32, 1> ENDPOINT664_MINIMUM = {
    UINT32_C(0x00c0c125)};

void add_scan_row(FnvLocal& ledger, const EdgeCegar& edge,
                  const CarrierScanCegar& scan) {
    ledger.add(edge.q);
    ledger.add(edge.r);
    ledger.add(scan.active);
    ledger.add(scan.body.failures);
    ledger.add(scan.body.checks);
    ledger.add(scan.body.max_checks);
    ledger.add(scan.body.first_failure);
    ledger.add(scan.best_disjoint_repair);
    fnv_add_i128_cegar(ledger, scan.best_disjoint_margin);
    fnv_add_i128_cegar(ledger, scan.denominator);
}

std::string hostile_row(const EdgeCegar& edge,
                        const CarrierScanCegar& scan) {
    std::ostringstream row;
    row << "RESISTANT PAIR " << edge.q << ',' << edge.r
        << " ACTIVE " << scan.active << " FAILURES "
        << scan.body.failures << " FIRST_FAILURE " << std::hex
        << scan.body.first_failure << std::dec << " FIRST_FAILURE_LABELS {"
        << labels(scan.body.first_failure) << "} BEST_DISJOINT_REPAIR "
        << std::hex << scan.best_disjoint_repair << std::dec
        << " BEST_DISJOINT_LABELS {" << labels(scan.best_disjoint_repair)
        << "} BEST_DISJOINT_MARGIN_NUM "
        << decimal(scan.best_disjoint_margin) << " DEN "
        << decimal(scan.denominator);
    return row.str();
}

}  // namespace

#ifndef ENDPOINT_DESCENT_SCAN_LIBRARY_ONLY
int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 9,
                "usage: endpoint668 REPLAY_ROOT PREFIX_416_704 PREFIX_416_700 PREFIX_520_700 PREFIX_384_694 PREFIX_520_688 RESIDUAL_CSV JOINT_DECK");
        const PrefixCegar prefix704 = read_prefix_cegar(argv[2]);
        const PrefixCegar prefix416 = read_prefix_cegar(argv[3]);
        const PrefixCegar prefix520 = read_prefix_cegar(argv[4]);
        const PrefixCegar prefix384 = read_prefix_cegar(argv[5]);
        const PrefixCegar prefix688 = read_prefix_cegar(argv[6]);
        std::vector<u32> carrier = build_thm4276(
            argv[1], prefix704, prefix416, prefix520, prefix384, prefix688);
        const std::vector<u32> joint = read_joint(argv[8]);
        std::set<u32> seen(carrier.begin(), carrier.end());
        append_unseen(joint, seen, carrier);
        require(carrier.size() == 8945 &&
                    mask_fnv(carrier) == UINT64_C(0x3212efa05dd18c00),
                "THM-4276 plus joint421 carrier changed");
        std::vector<u32> minimum;
        for (u32 mask : ENDPOINT669_MINIMUM) {
            require(std::popcount(mask) == 8, "minimum repair arity changed");
            require(seen.insert(mask).second,
                    "minimum repair overlaps prior carrier");
            carrier.push_back(mask);
            minimum.push_back(mask);
        }
        require(carrier.size() == 8947,
                "final augmented carrier count changed");

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool cell count changed");
        const CarrierScanCegar repaired669 =
            scan_carrier_cegar(cells, 256, 669, carrier);
        require(repaired669.body.failures == 0,
                "two-mask certificate does not repair 256,669");
        CarrierScanCegar repaired667;
        std::vector<u32> minimum667;
        if constexpr (INCLUDE_ENDPOINT667_MINIMUM != 0) {
            for (u32 mask : ENDPOINT667_MINIMUM) {
                require(std::popcount(mask) == 8,
                        "endpoint667 minimum repair arity changed");
                require(seen.insert(mask).second,
                        "endpoint667 minimum repair overlaps prior carrier");
                carrier.push_back(mask);
                minimum667.push_back(mask);
            }
            require(carrier.size() == 8949,
                    "endpoint667 augmented carrier count changed");
            repaired667 = scan_carrier_cegar(cells, 256, 667, carrier);
            require(repaired667.body.failures == 0,
                    "two-mask certificate does not repair 256,667");
        }
        CarrierScanCegar repaired665;
        std::vector<u32> minimum665;
        if constexpr (INCLUDE_ENDPOINT665_MINIMUM != 0) {
            require(INCLUDE_ENDPOINT667_MINIMUM != 0,
                    "endpoint665 repair requires endpoint667 carrier");
            for (u32 mask : ENDPOINT665_MINIMUM) {
                require(std::popcount(mask) == 8,
                        "endpoint665 minimum repair arity changed");
                require(seen.insert(mask).second,
                        "endpoint665 minimum repair overlaps prior carrier");
                carrier.push_back(mask);
                minimum665.push_back(mask);
            }
            require(carrier.size() == 8950,
                    "endpoint665 augmented carrier count changed");
            repaired665 = scan_carrier_cegar(cells, 520, 665, carrier);
            require(repaired665.body.failures == 0,
                    "one-mask certificate does not repair 520,665");
        }
        CarrierScanCegar repaired664;
        std::vector<u32> minimum664;
        if constexpr (INCLUDE_ENDPOINT664_MINIMUM != 0) {
            require(INCLUDE_ENDPOINT665_MINIMUM != 0,
                    "endpoint664 repair requires endpoint665 carrier");
            for (u32 mask : ENDPOINT664_MINIMUM) {
                require(std::popcount(mask) == 8,
                        "endpoint664 minimum repair arity changed");
                require(seen.insert(mask).second,
                        "endpoint664 minimum repair overlaps prior carrier");
                carrier.push_back(mask);
                minimum664.push_back(mask);
            }
            require(carrier.size() == 8951,
                    "endpoint664 augmented carrier count changed");
            repaired664 = scan_carrier_cegar(cells, 256, 664, carrier);
            require(repaired664.body.failures == 0,
                    "one-mask certificate does not repair 256,664");
        }

        const std::vector<EdgeCegar> edges = read_edges_cegar(argv[7]);
        u64 rows = 0;
        u64 resistant = 0;
        u64 active_sum = 0;
        u64 checks = 0;
        u64 max_checks = 0;
        FnvLocal ledger;
        std::vector<std::string> hostile;
        std::cout << "AUGMENTED_ENDPOINT_DESCENT_V1\n";
        std::cout << "BASE_CARRIER 8524 FNV 5ddb84a44f5d2ad7\n";
        std::cout << "JOINT_DECK " << joint.size() << " FNV " << std::hex
                  << mask_fnv(joint) << std::dec << " OVERLAP 0\n";
        std::cout << "ENDPOINT669_MINIMUM 2 MASKS " << std::hex
                  << ENDPOINT669_MINIMUM[0] << ',' << ENDPOINT669_MINIMUM[1]
                  << " FNV " << mask_fnv(minimum) << std::dec
                  << " OVERLAP 0\n";
        std::cout << "AUGMENTED_CARRIER " << carrier.size() << " FNV "
                  << std::hex << mask_fnv(carrier) << std::dec << '\n';
        std::cout << "REPAIRED PAIR 256,669 ACTIVE " << repaired669.active
                  << " FAILURES " << repaired669.body.failures
                  << " CHECKS " << repaired669.body.checks
                  << " MAX_CHECKS " << repaired669.body.max_checks << '\n';
        if constexpr (INCLUDE_ENDPOINT667_MINIMUM != 0) {
            std::cout << "ENDPOINT667_MINIMUM 2 MASKS " << std::hex
                      << ENDPOINT667_MINIMUM[0] << ','
                      << ENDPOINT667_MINIMUM[1] << " FNV "
                      << mask_fnv(minimum667) << std::dec << " OVERLAP 0\n";
            std::cout << "REPAIRED PAIR 256,667 ACTIVE "
                      << repaired667.active << " FAILURES "
                      << repaired667.body.failures << " CHECKS "
                      << repaired667.body.checks << " MAX_CHECKS "
                      << repaired667.body.max_checks << '\n';
        }
        if constexpr (INCLUDE_ENDPOINT665_MINIMUM != 0) {
            std::cout << "ENDPOINT665_MINIMUM 1 MASK " << std::hex
                      << ENDPOINT665_MINIMUM[0] << " FNV "
                      << mask_fnv(minimum665) << std::dec << " OVERLAP 0\n";
            std::cout << "REPAIRED PAIR 520,665 ACTIVE "
                      << repaired665.active << " FAILURES "
                      << repaired665.body.failures << " CHECKS "
                      << repaired665.body.checks << " MAX_CHECKS "
                      << repaired665.body.max_checks << '\n';
        }
        if constexpr (INCLUDE_ENDPOINT664_MINIMUM != 0) {
            std::cout << "ENDPOINT664_MINIMUM 1 MASK " << std::hex
                      << ENDPOINT664_MINIMUM[0] << " FNV "
                      << mask_fnv(minimum664) << std::dec << " OVERLAP 0\n";
            std::cout << "REPAIRED PAIR 256,664 ACTIVE "
                      << repaired664.active << " FAILURES "
                      << repaired664.body.failures << " CHECKS "
                      << repaired664.body.checks << " MAX_CHECKS "
                      << repaired664.body.max_checks << '\n';
        }
        for (const EdgeCegar& edge : edges) {
            if (edge.r != SCAN_ENDPOINT) continue;
            const CarrierScanCegar scan =
                scan_carrier_cegar(cells, edge.q, edge.r, carrier);
            ++rows;
            resistant += scan.body.failures != 0;
            active_sum += scan.active;
            checks += scan.body.checks;
            max_checks = std::max(max_checks, scan.body.max_checks);
            add_scan_row(ledger, edge, scan);
            if (scan.body.failures != 0)
                hostile.push_back(hostile_row(edge, scan));
            if constexpr (QUIET_ENDPOINT_PROGRESS == 0) {
                if (rows % 16 == 0)
                std::cerr << "PROGRESS ROWS " << rows << " RESISTANT "
                          << resistant << '\n';
            }
        }
        const u64 expected_rows = SCAN_ENDPOINT == 668 ? 170 :
                                  SCAN_ENDPOINT == 667 ? 176 :
                                  SCAN_ENDPOINT == 666 ? 179 :
                                  SCAN_ENDPOINT == 665 ? 183 :
                                  SCAN_ENDPOINT == 664 ? 190 : 0;
        require(expected_rows != 0 && rows == expected_rows,
                "endpoint row count changed");
        const u64 expected_resistant =
            SCAN_ENDPOINT == 667 || SCAN_ENDPOINT == 665 ? 1 :
            SCAN_ENDPOINT == 664 && INCLUDE_ENDPOINT664_MINIMUM == 0 ? 1 : 0;
        const u64 expected_active_sum =
            SCAN_ENDPOINT == 668 ? UINT64_C(1502443) :
            SCAN_ENDPOINT == 667 ? UINT64_C(1536773) :
            SCAN_ENDPOINT == 666 ? UINT64_C(1583997) :
            SCAN_ENDPOINT == 665 ? UINT64_C(1565443) :
            INCLUDE_ENDPOINT664_MINIMUM != 0 ? UINT64_C(1660533) :
                                               UINT64_C(1660345);
        const u64 expected_checks =
            SCAN_ENDPOINT == 668 ? UINT64_C(71368851159) :
            SCAN_ENDPOINT == 667 ? UINT64_C(74521089350) :
            SCAN_ENDPOINT == 666 ? UINT64_C(75099164361) :
            SCAN_ENDPOINT == 665 ? UINT64_C(78347641341) :
            INCLUDE_ENDPOINT664_MINIMUM != 0 ? UINT64_C(80180889036) :
                                               UINT64_C(80180889034);
        const u64 expected_max_checks =
            SCAN_ENDPOINT == 668 ? UINT64_C(6511) :
            SCAN_ENDPOINT == 667 ? UINT64_C(4199) :
            SCAN_ENDPOINT == 666 ? UINT64_C(5344) :
            SCAN_ENDPOINT == 665 ? UINT64_C(5907) : UINT64_C(6307);
        const u64 expected_ledger =
            SCAN_ENDPOINT == 668 ? UINT64_C(0x05b5c447f482c7d7) :
            SCAN_ENDPOINT == 667 ? UINT64_C(0x94237267aec84f12) :
            SCAN_ENDPOINT == 666 ? UINT64_C(0x6b834428d86c0231) :
            SCAN_ENDPOINT == 665 ? UINT64_C(0x22942b499df769d0) :
            INCLUDE_ENDPOINT664_MINIMUM != 0 ?
                UINT64_C(0x2973e6c6d40d5af8) :
                UINT64_C(0x172cd995dd68cdd6);
        require(resistant == expected_resistant &&
                    active_sum == expected_active_sum &&
                    checks == expected_checks &&
                    max_checks == expected_max_checks &&
                    ledger.state == expected_ledger,
                "frozen endpoint-layer ledger changed");
        std::cout << "LAYER " << SCAN_ENDPOINT << " ROWS " << rows
                  << " RESISTANT " << resistant
                  << " ACTIVE_SUM " << active_sum << " CHECKS " << checks
                  << " MAX_CHECKS " << max_checks << " ROW_LEDGER "
                  << std::hex << ledger.state << std::dec << '\n';
        for (const std::string& row : hostile) std::cout << row << '\n';
        std::cout << "VERDICT "
                  << (resistant == 0 ? "RANGE_CLOSED" : "BOUNDARY_FOUND")
                  << " EXACT_PRIMARY_ENDPOINT_DESCENT\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "AUGMENTED_ENDPOINT_DESCENT_ERROR " << error.what() << '\n';
        return 1;
    }
}
#endif
