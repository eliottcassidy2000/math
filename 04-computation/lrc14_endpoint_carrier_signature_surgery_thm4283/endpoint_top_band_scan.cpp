// Exact top-band carrier scan for the reserved THM-4283 namespace.
//
// This wrapper reuses THM-4282's source-pinned joint-exposure implementation
// but replaces its frozen 90-row main program.  It verifies the complete
// post-THM-4282 residual, scans endpoint layers in descending order, compares
// the inherited 8,996-mask carrier with the nested one-mask augmentation, and
// stops after the first layer on which the nested carrier has a failure.

#include "carrier_scan_support.cpp"

#include <algorithm>
#include <tuple>

namespace {

constexpr std::size_t kResidualCount = 23373;
constexpr u64 kResidualFnv = UINT64_C(0xc6ab0ae49ee32273);
constexpr u32 kNestedRepair = UINT32_C(0x014c9084);
constexpr u64 kNestedCarrierFnv = UINT64_C(0x8e1860a25d0fcf87);

std::vector<Pair> read_residual_band(const std::filesystem::path& path,
                                     int lower_endpoint) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open final residual");
    std::vector<Pair> all;
    std::vector<Pair> band;
    FnvLocal all_ledger;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos &&
                    line.find(',', comma + 1) == std::string::npos,
                "malformed residual row");
        std::size_t used_q = 0;
        std::size_t used_r = 0;
        Pair pair;
        pair.q = std::stoi(line.substr(0, comma), &used_q);
        pair.r = std::stoi(line.substr(comma + 1), &used_r);
        require(used_q == comma && used_r == line.size() - comma - 1 &&
                    pair.q > 0 && pair.q < pair.r,
                "invalid residual pair");
        if (!all.empty())
            require(std::tie(all.back().q, all.back().r) <
                        std::tie(pair.q, pair.r),
                    "residual order changed");
        all.push_back(pair);
        all_ledger.add(pair.q);
        all_ledger.add(pair.r);
        if (pair.r >= lower_endpoint) band.push_back(pair);
    }
    require(all.size() == kResidualCount && all_ledger.state == kResidualFnv,
            "post-THM4282 residual identity changed");
    std::sort(band.begin(), band.end(), [](const Pair& left, const Pair& right) {
        if (left.r != right.r) return left.r > right.r;
        return left.q < right.q;
    });
    require(!band.empty(), "selected band is empty");
    return band;
}

u64 pair_list_fnv(const std::vector<Pair>& pairs) {
    FnvLocal ledger;
    for (const Pair pair : pairs) {
        ledger.add(pair.q);
        ledger.add(pair.r);
    }
    return ledger.state;
}

std::vector<u32> nested_failures(const PairAudit& base, bool repair_active) {
    std::vector<u32> failures;
    for (u32 body : base.failure_bodies)
        if (!repair_active || (body & kNestedRepair) != 0)
            failures.push_back(body);
    return failures;
}

u64 body_fnv(const std::vector<u32>& bodies) {
    FnvLocal ledger;
    for (u32 body : bodies) ledger.add(body);
    return ledger.state;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 8,
                "usage: endpoint-top-band JOINT421 BASE8951 ADDITIONS45 "
                "FINAL23373 LOWER_ENDPOINT BASE_FAILURES NESTED_FAILURES");
        init_choose8_local();
        const int lower_endpoint = std::stoi(argv[5]);
        require(lower_endpoint > 0, "lower endpoint must be positive");
        const std::vector<u32> joint =
            read_masks(argv[1], kJointCount, kJointFnv);
        std::vector<u32> base =
            read_masks(argv[2], kCarrierCount, kCarrierFnv);
        const std::vector<u32> additions = read_additions(argv[3]);
        std::set<u32> base_set(base.begin(), base.end());
        for (u32 mask : additions) {
            require(base_set.insert(mask).second,
                    "carrier addition overlaps inherited carrier");
            base.push_back(mask);
        }
        require(base.size() == kCarrierCount + kAdditionCount &&
                    mask_fnv(base) == kAugmentedCarrierFnv,
                "final8996 identity changed");
        require(std::popcount(kNestedRepair) == 8 &&
                    !base_set.contains(kNestedRepair),
                "nested repair rank or novelty changed");
        FnvLocal nested_carrier_ledger;
        for (u32 mask : base) nested_carrier_ledger.add(mask);
        nested_carrier_ledger.add(kNestedRepair);
        require(nested_carrier_ledger.state == kNestedCarrierFnv,
                "nested8997 identity changed");
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        for (u32 mask : joint)
            require(base_set.contains(mask), "joint mask absent from carrier");
        const std::vector<Pair> band =
            read_residual_band(argv[4], lower_endpoint);
        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool-cell count changed");

        std::ofstream base_output(argv[6]);
        std::ofstream nested_output(argv[7]);
        require(static_cast<bool>(base_output) &&
                    static_cast<bool>(nested_output),
                "cannot create failure ledgers");
        base_output << "q,r,body_hex\n";
        nested_output << "q,r,body_hex\n";

        std::cout << "THM4283_ENDPOINT_TOP_BAND_SCAN_V1\n"
                  << "RESIDUAL " << kResidualCount << " FNV " << std::hex
                  << kResidualFnv << std::dec << " REQUESTED_LOWER "
                  << lower_endpoint << " SELECTED " << band.size()
                  << " PAIR_FNV " << std::hex << pair_list_fnv(band)
                  << std::dec << " BASE8996_FNV " << std::hex
                  << kAugmentedCarrierFnv << " NESTED8997_FNV "
                  << kNestedCarrierFnv << std::dec << '\n';

        FnvLocal pair_ledger;
        u64 total_base_failures = 0;
        u64 total_nested_failures = 0;
        bool stopped = false;
        std::size_t index = 0;
        while (index < band.size() && !stopped) {
            const int endpoint = band[index].r;
            u64 layer_base_failures = 0;
            u64 layer_nested_failures = 0;
            std::size_t layer_rows = 0;
            FnvLocal layer_ledger;
            while (index < band.size() && band[index].r == endpoint) {
                const Pair pair = band[index++];
                const PairAudit audit =
                    audit_pair(cells, joint, base, joint_set, pair);
                const ActiveUniverse active =
                    build_active_universe(cells, pair.q, pair.r);
                const bool repair_active =
                    active.active[colex_rank8_local(kNestedRepair)] != 0;
                const std::vector<u32> nested =
                    nested_failures(audit, repair_active);
                for (u32 body : audit.failure_bodies)
                    base_output << pair.q << ',' << pair.r << ',' << std::hex
                                << std::setw(8) << std::setfill('0') << body
                                << std::dec << std::setfill(' ') << '\n';
                for (u32 body : nested)
                    nested_output << pair.q << ',' << pair.r << ','
                                  << std::hex << std::setw(8)
                                  << std::setfill('0') << body << std::dec
                                  << std::setfill(' ') << '\n';
                layer_base_failures += audit.failures;
                layer_nested_failures += nested.size();
                total_base_failures += audit.failures;
                total_nested_failures += nested.size();
                ++layer_rows;
                layer_ledger.add(pair.q);
                layer_ledger.add(pair.r);
                pair_ledger.add(pair.q);
                pair_ledger.add(pair.r);
                pair_ledger.add(audit.exposed);
                pair_ledger.add(audit.exposed_fnv);
                pair_ledger.add(audit.failures);
                pair_ledger.add(audit.failure_fnv);
                pair_ledger.add(repair_active);
                pair_ledger.add(nested.size());
                pair_ledger.add(body_fnv(nested));
                std::cout << "PAIR " << pair.q << ',' << pair.r
                          << " EXPOSED " << audit.exposed << " EXPOSED_FNV "
                          << std::hex << audit.exposed_fnv << std::dec
                          << " BASE_FAILURES " << audit.failures
                          << " BASE_FAILURE_FNV " << std::hex
                          << audit.failure_fnv << std::dec
                          << " REPAIR_ACTIVE " << repair_active
                          << " NESTED_FAILURES " << nested.size()
                          << " NESTED_FAILURE_FNV " << std::hex
                          << body_fnv(nested) << std::dec;
                if (audit.failures != 0 && audit.failures < 64)
                    std::cout << " RESPONSE_CLASSES "
                              << audit.complete_response_classes
                              << " FULL_RESPONDERS "
                              << audit.full_response_masks << " LEAST_FULL "
                              << std::hex << std::setw(8)
                              << std::setfill('0')
                              << audit.least_full_response << std::dec
                              << std::setfill(' ') << " EXACT_MINIMUM "
                              << audit.exact_minimum_replacements;
                std::cout << '\n';
            }
            std::cout << "LAYER " << endpoint << " ROWS " << layer_rows
                      << " FNV " << std::hex << layer_ledger.state
                      << std::dec << " BASE_FAILURES "
                      << layer_base_failures << " NESTED_FAILURES "
                      << layer_nested_failures << '\n';
            stopped = layer_nested_failures != 0;
        }
        require(base_output.good() && nested_output.good(),
                "failure-ledger write failed");
        std::cout << "PAIR_LEDGER_FNV " << std::hex << pair_ledger.state
                  << std::dec << " TOTAL_BASE_FAILURES "
                  << total_base_failures << " TOTAL_NESTED_FAILURES "
                  << total_nested_failures << " STOPPED_AT_FAILURE "
                  << stopped << '\n'
                  << "VERDICT PASS EXACT_DESCENDING_SCAN\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4283_TOP_BAND_ERROR " << error.what() << '\n';
        return 1;
    }
}
