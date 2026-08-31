// Exact audit: append the exact-minimum 14-mask response witness to the
// THM-4287 carrier and descend the live 22,647-row residual by whole layers.
#include "04-computation/lrc14_endpoint_carrier_signature_surgery_thm4283/carrier_scan_support.cpp"

namespace {

constexpr std::size_t kLiveCount = 22647;
constexpr u64 kLiveFnv = UINT64_C(0xdf5374d4aca67677);
constexpr std::size_t kAppendCount = 14;
constexpr u32 kPriorRepair = UINT32_C(0x014c9084);
constexpr std::size_t kBoundaryWitnessCount = 9;
constexpr u64 kBoundaryWitnessFnv = UINT64_C(0x02b936529030e4bc);
constexpr u64 kRepairedCarrierFnv = UINT64_C(0xfdc1c57ae4dc1bb6);

u64 local_mask_fnv(const std::vector<u32>& masks) {
    FnvLocal ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

u64 pair_list_fnv(const std::vector<Pair>& pairs) {
    FnvLocal ledger;
    for (Pair pair : pairs) {
        ledger.add(pair.q);
        ledger.add(pair.r);
    }
    return ledger.state;
}

std::vector<u32> read_append14(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open append14 mask ledger");
    std::vector<u32> masks;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        const u64 wide = std::stoull(token, nullptr, 16);
        require(wide < (UINT64_C(1) << 30), "append mask outside pool");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "append mask rank/distinctness changed");
        masks.push_back(mask);
    }
    require(input.eof() && masks.size() == kAppendCount,
            "append14 mask count changed");
    return masks;
}

std::vector<Pair> read_live_band(const std::filesystem::path& path,
                                 int lower_endpoint) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open live residual");
    std::vector<Pair> all;
    std::vector<Pair> band;
    FnvLocal ledger;
    std::string line;
    while (std::getline(input, line)) {
        require(!line.empty(), "blank live residual row");
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos &&
                    line.find(',', comma + 1) == std::string::npos,
                "malformed live residual row");
        Pair pair{std::stoi(line.substr(0, comma)),
                  std::stoi(line.substr(comma + 1))};
        require(pair.q > 0 && pair.q < pair.r &&
                    (all.empty() ||
                     std::tie(all.back().q, all.back().r) <
                         std::tie(pair.q, pair.r)),
                "invalid live residual order");
        all.push_back(pair);
        ledger.add(pair.q);
        ledger.add(pair.r);
        if (pair.r >= lower_endpoint) band.push_back(pair);
    }
    require(input.eof() && all.size() == kLiveCount &&
                ledger.state == kLiveFnv,
            "live THM-4287 residual identity changed");
    std::sort(band.begin(), band.end(), [](Pair a, Pair b) {
        if (a.r != b.r) return a.r > b.r;
        return a.q < b.q;
    });
    require(!band.empty() && band.front().q == 100 && band.front().r == 636,
            "live residual top changed");
    return band;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 9,
                "usage: append14-scan JOINT BASE ADDITIONS45 WITNESS9 "
                "LIVE22647 APPEND14 LOWER FAILURES.csv");
        init_choose8_local();
        const auto joint = read_masks(argv[1], kJointCount, kJointFnv);
        auto carrier = read_masks(argv[2], kCarrierCount, kCarrierFnv);
        const auto additions = read_additions(argv[3]);
        const auto witness9 =
            read_masks(argv[4], kBoundaryWitnessCount, kBoundaryWitnessFnv);
        std::set<u32> distinct(carrier.begin(), carrier.end());
        for (u32 mask : additions) {
            require(distinct.insert(mask).second, "addition45 overlap");
            carrier.push_back(mask);
        }
        require(carrier.size() == 8996 &&
                    local_mask_fnv(carrier) == kAugmentedCarrierFnv,
                "base8996 changed");
        require(distinct.insert(kPriorRepair).second, "prior repair overlap");
        carrier.push_back(kPriorRepair);
        for (u32 mask : witness9) {
            require(distinct.insert(mask).second, "witness9 overlap");
            carrier.push_back(mask);
        }
        require(carrier.size() == 9006 &&
                    local_mask_fnv(carrier) == kRepairedCarrierFnv,
                "repaired9006 changed");
        const auto append14 = read_append14(argv[6]);
        for (u32 mask : append14) {
            require(distinct.insert(mask).second,
                    "append14 overlaps repaired carrier");
            carrier.push_back(mask);
        }
        const u64 carrier_fnv = local_mask_fnv(carrier);
        const u64 append_fnv = local_mask_fnv(append14);
        const int lower = std::stoi(argv[7]);
        require(lower > 0 && lower <= 636, "bad lower endpoint");
        const auto band = read_live_band(argv[5], lower);
        const auto cells = build_pool_cells();
        require(cells.size() == 7133, "pool cells changed");
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        std::ofstream failures(argv[8]);
        require(static_cast<bool>(failures), "cannot create failure ledger");
        failures << "q,r,body_hex\n";

        std::cout << "THM4295_APPEND14_LIVE_DESCENDING_V1\n"
                  << "LIVE_RESIDUAL " << kLiveCount << " FNV " << std::hex
                  << kLiveFnv << std::dec << " SELECTED " << band.size()
                  << " PAIR_FNV " << std::hex << pair_list_fnv(band)
                  << std::dec << " LOWER " << lower << '\n'
                  << "BASE9006_FNV " << std::hex << kRepairedCarrierFnv
                  << " APPEND14_FNV " << append_fnv << " CARRIER9020_FNV "
                  << carrier_fnv << std::dec << '\n';

        FnvLocal pair_ledger;
        u64 total_failures = 0;
        std::size_t completed_rows = 0;
        bool stopped = false;
        std::size_t index = 0;
        while (index < band.size() && !stopped) {
            const int endpoint = band[index].r;
            u64 layer_failures = 0;
            std::size_t layer_rows = 0;
            FnvLocal layer_ledger;
            while (index < band.size() && band[index].r == endpoint) {
                const Pair pair = band[index++];
                const PairAudit audit = audit_pair(
                    cells, joint, carrier, joint_set, pair);
                for (u32 body : audit.failure_bodies) {
                    failures << pair.q << ',' << pair.r << ',' << std::hex
                             << std::setw(8) << std::setfill('0') << body
                             << std::dec << std::setfill(' ') << '\n';
                }
                layer_failures += audit.failures;
                total_failures += audit.failures;
                ++layer_rows;
                ++completed_rows;
                layer_ledger.add(pair.q);
                layer_ledger.add(pair.r);
                pair_ledger.add(pair.q);
                pair_ledger.add(pair.r);
                pair_ledger.add(audit.active_carrier);
                pair_ledger.add(audit.active_carrier_fnv);
                pair_ledger.add(audit.exposed);
                pair_ledger.add(audit.exposed_fnv);
                pair_ledger.add(audit.failures);
                pair_ledger.add(audit.failure_fnv);
                pair_ledger.add(audit.hit_incidences);
                pair_ledger.add(audit.minimum_hits);
                pair_ledger.add(audit.maximum_hits);
                pair_ledger.add(audit.response_fnv);
                std::cout << "PAIR " << pair.q << ',' << pair.r
                          << " ACTIVE_CARRIER " << audit.active_carrier
                          << " ACTIVE_FNV " << std::hex
                          << audit.active_carrier_fnv << std::dec
                          << " ACTIVE_JOINT " << audit.active_joint
                          << " ACTIVE_NONJOINT " << audit.active_nonjoint
                          << " INACTIVE_JOINT " << audit.inactive_joint
                          << " EXPOSED " << audit.exposed << " EXPOSED_FNV "
                          << std::hex << audit.exposed_fnv << std::dec
                          << " FAILURES " << audit.failures << " FAILURE_FNV "
                          << std::hex << audit.failure_fnv << std::dec
                          << " HIT_INCIDENCES " << audit.hit_incidences
                          << " HIT_RANGE " << audit.minimum_hits << ".."
                          << audit.maximum_hits << " RESPONSE_FNV " << std::hex
                          << audit.response_fnv << std::dec << '\n';
            }
            std::cout << "LAYER " << endpoint << " ROWS " << layer_rows
                      << " FNV " << std::hex << layer_ledger.state << std::dec
                      << " FAILURES " << layer_failures << '\n';
            stopped = layer_failures != 0;
        }
        require(failures.good(), "failure-ledger write failed");
        std::cout << "PAIR_LEDGER_FNV " << std::hex << pair_ledger.state
                  << std::dec << " COMPLETED_ROWS " << completed_rows
                  << " TOTAL_FAILURES " << total_failures
                  << " STOPPED_AT_FAILURE " << stopped << '\n'
                  << "SCOPE FIXED_POOL_LIVE_RESIDUAL_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS EXACT_APPEND14_DESCENDING_AUDIT\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4295_APPEND14_DESCENDING_ERROR " << error.what() << '\n';
        return 1;
    }
}
