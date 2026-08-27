// Exact primary endpoint-atom scan of the THM-4276 carrier augmented by the
// frozen 421-mask rectangle-common/r670 joint deck, restricted to the entire
// current endpoint-669 layer of the post-THM-4271 semantic residual.

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#endif
#define CARRIER_CEGAR_LIBRARY_ONLY
#include "04-computation/lrc14_three_round_learned_carrier_thm4266/carrier_cegar_descent.cpp"
#undef CARRIER_CEGAR_LIBRARY_ONLY
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

namespace {

constexpr std::size_t EXPECTED_JOINT = 421;
constexpr u64 EXPECTED_JOINT_FNV = UINT64_C(0x20d63dd42fe8150e);

void append_unseen(const std::vector<u32>& source, std::set<u32>& seen,
                   std::vector<u32>& target) {
    for (u32 mask : source) {
        require(seen.insert(mask).second, "unexpected carrier/deck overlap");
        target.push_back(mask);
    }
}

u64 mask_fnv(const std::vector<u32>& masks) {
    FnvLocal ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

std::vector<u32> read_joint(const std::filesystem::path& path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open joint deck");
    std::vector<u32> masks;
    std::set<u32> seen;
    std::string word;
    while (in >> word) {
        const u64 wide = std::stoull(word, nullptr, 16);
        require(wide < (UINT64_C(1) << 30), "joint mask leaves pool");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8, "joint mask not rank eight");
        require(seen.insert(mask).second, "duplicate joint mask");
        masks.push_back(mask);
    }
    require(masks.size() == EXPECTED_JOINT &&
                mask_fnv(masks) == EXPECTED_JOINT_FNV,
            "joint-deck ledger changed");
    return masks;
}

void append_allow_seen(const std::vector<u32>& source, std::set<u32>& seen,
                       std::vector<u32>& target) {
    for (u32 mask : source) {
        if (seen.insert(mask).second) target.push_back(mask);
    }
}

std::vector<u32> build_thm4276(
    const std::filesystem::path& replay_root,
    const PrefixCegar& prefix704, const PrefixCegar& prefix416,
    const PrefixCegar& prefix520, const PrefixCegar& prefix384,
    const PrefixCegar& prefix688) {
    require(prefix704.q == 416 && prefix704.r == 704 &&
                prefix704.masks.size() == 2608 &&
                prefix704.stated_fnv == UINT64_C(0x18ff663a123e684e),
            "416,704 prefix changed");
    require(prefix416.q == 416 && prefix416.r == 700 &&
                prefix416.masks.size() == 5894 &&
                prefix416.stated_fnv == UINT64_C(0x701c0233c8f8abeb),
            "416,700 prefix changed");
    require(prefix520.q == 520 && prefix520.r == 700 &&
                prefix520.masks.size() == 4557 &&
                prefix520.stated_fnv == UINT64_C(0xd8466558bcd8ef9d),
            "520,700 prefix changed");
    require(prefix384.q == 384 && prefix384.r == 694 &&
                prefix384.masks.size() == 5307 &&
                prefix384.stated_fnv == UINT64_C(0x7b9a46c60e8514f8),
            "384,694 prefix changed");
    require(prefix688.q == 520 && prefix688.r == 688 &&
                prefix688.masks.size() == 5398 &&
                prefix688.stated_fnv == UINT64_C(0x6ab471da88c8e1d1),
            "520,688 prefix changed");

    std::vector<u32> carrier;
    std::set<u32> seen;
    for (const EdgeCegar expected : FROZEN_BAND_CEGAR) {
        const PrefixCegar prefix = read_prefix_cegar(
            replay_root / "replay_band" /
            (std::to_string(expected.q) + "_" + std::to_string(expected.r) +
             ".out"));
        require(prefix.q == expected.q && prefix.r == expected.r,
                "base transcript changed");
        append_allow_seen(prefix.masks, seen, carrier);
    }
    require(carrier.size() == 4675 &&
                mask_fnv(carrier) == UINT64_C(0xce4e76ec11df057c),
            "base carrier changed");
    append_allow_seen(prefix704.masks, seen, carrier);
    require(carrier.size() == 4733 &&
                mask_fnv(carrier) == UINT64_C(0xa7b046289655c733),
            "round-one carrier changed");
    const std::set<u32> round1_seen = seen;
    std::vector<u32> novel416;
    std::vector<u32> novel520;
    for (u32 mask : prefix416.masks) {
        if (!round1_seen.count(mask)) novel416.push_back(mask);
    }
    for (u32 mask : prefix520.masks) {
        if (!round1_seen.count(mask)) novel520.push_back(mask);
    }
    append_allow_seen(novel416, seen, carrier);
    append_allow_seen(novel520, seen, carrier);
    require(carrier.size() == 7986 &&
                mask_fnv(carrier) == UINT64_C(0xbaef1d2f49444638),
            "round-two carrier changed");
    append_allow_seen(prefix384.masks, seen, carrier);
    require(carrier.size() == 8319 &&
                mask_fnv(carrier) == UINT64_C(0xe08b227730f6793c),
            "round-three carrier changed");
    append_allow_seen(prefix688.masks, seen, carrier);
    require(carrier.size() == 8518 &&
                mask_fnv(carrier) == UINT64_C(0x1603e3fe970f8428),
            "THM-4271 carrier changed");
    const std::vector<u32> promoted_six = {
        UINT32_C(0x00289285), UINT32_C(0x0260812c),
        UINT32_C(0x18689040), UINT32_C(0x20c0c124),
        UINT32_C(0x302c1006), UINT32_C(0x30580888)};
    append_unseen(promoted_six, seen, carrier);
    require(carrier.size() == 8524 &&
                mask_fnv(carrier) == UINT64_C(0x5ddb84a44f5d2ad7),
            "THM-4276 carrier changed");
    return carrier;
}

}  // namespace

#ifndef ENDPOINT669_SCAN_LIBRARY_ONLY
int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 9,
                "usage: endpoint669 REPLAY_ROOT PREFIX_416_704 PREFIX_416_700 PREFIX_520_700 PREFIX_384_694 PREFIX_520_688 RESIDUAL_CSV JOINT_DECK");
        const PrefixCegar prefix704 = read_prefix_cegar(argv[2]);
        const PrefixCegar prefix416 = read_prefix_cegar(argv[3]);
        const PrefixCegar prefix520 = read_prefix_cegar(argv[4]);
        const PrefixCegar prefix384 = read_prefix_cegar(argv[5]);
        const PrefixCegar prefix688 = read_prefix_cegar(argv[6]);
        std::vector<u32> carrier = build_thm4276(
            argv[1], prefix704, prefix416, prefix520, prefix384, prefix688);
        const std::vector<u32> joint = read_joint(argv[8]);
        std::set<u32> seen(carrier.begin(), carrier.end());
        const std::size_t before = seen.size();
        append_unseen(joint, seen, carrier);
        require(seen.size() == before + joint.size() && carrier.size() == 8945,
                "augmented carrier count changed");

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool cell count changed");
        const std::vector<EdgeCegar> edges = read_edges_cegar(argv[7]);
        u64 rows = 0;
        u64 resistant = 0;
        u64 active_sum = 0;
        u64 checks = 0;
        u64 max_checks = 0;
        FnvLocal ledger;
        std::vector<std::string> hostile;
        std::cout << "AUGMENTED_ENDPOINT669_V1\n";
        std::cout << "BASE_CARRIER 8524 FNV 5ddb84a44f5d2ad7\n";
        std::cout << "JOINT_DECK " << joint.size() << " FNV " << std::hex
                  << mask_fnv(joint) << std::dec << " OVERLAP 0\n";
        std::cout << "AUGMENTED_CARRIER " << carrier.size() << " FNV "
                  << std::hex << mask_fnv(carrier) << std::dec << '\n';
        for (const EdgeCegar& edge : edges) {
            if (edge.r != 669) continue;
            const CarrierScanCegar scan =
                scan_carrier_cegar(cells, edge.q, edge.r, carrier);
            ++rows;
            resistant += scan.body.failures != 0;
            active_sum += scan.active;
            checks += scan.body.checks;
            max_checks = std::max(max_checks, scan.body.max_checks);
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
            if (scan.body.failures != 0) {
                std::ostringstream row;
                row << "RESISTANT PAIR " << edge.q << ',' << edge.r
                    << " ACTIVE " << scan.active << " FAILURES "
                    << scan.body.failures << " FIRST_FAILURE " << std::hex
                    << scan.body.first_failure << std::dec
                    << " BEST_DISJOINT_REPAIR " << std::hex
                    << scan.best_disjoint_repair << std::dec
                    << " BEST_DISJOINT_MARGIN_NUM "
                    << decimal(scan.best_disjoint_margin) << " DEN "
                    << decimal(scan.denominator);
                hostile.push_back(row.str());
            }
            if (rows % 16 == 0) {
                std::cerr << "PROGRESS ROWS " << rows << " RESISTANT "
                          << resistant << '\n';
            }
        }
        require(rows == 168, "endpoint-669 row count changed");
        require(resistant == 1 && active_sum == UINT64_C(1455317) &&
                    checks == UINT64_C(71524673001) &&
                    max_checks == UINT64_C(3359) &&
                    ledger.state == UINT64_C(0x1b12eb2fb263db7e),
                "frozen endpoint-669 layer ledger changed");
        std::cout << "LAYER 669 ROWS " << rows << " RESISTANT " << resistant
                  << " ACTIVE_SUM " << active_sum << " CHECKS " << checks
                  << " MAX_CHECKS " << max_checks << " ROW_LEDGER "
                  << std::hex << ledger.state << std::dec << '\n';
        for (const std::string& row : hostile) std::cout << row << '\n';
        std::cout << "VERDICT " << (resistant == 0 ? "RANGE_CLOSED" : "BOUNDARY_FOUND")
                  << " EXACT_PRIMARY_ENDPOINT669\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "AUGMENTED_ENDPOINT669_ERROR " << error.what() << '\n';
        return 1;
    }
}
#endif
