// Scratch independent all-body replay for the 9,019-mask r>=617 exchanged
// carrier.  The sign audit only selects deletions; this program re-enumerates
// all C(30,9)=14,307,150 bodies at every residual pair in the prefix.

#define ENDPOINT617_EXCHANGE_MAIN endpoint617_exchange_hidden_main
#include "endpoint617_size_preserving_exchange_scan.cpp"
#undef ENDPOINT617_EXCHANGE_MAIN

#include <limits>

namespace {

constexpr u64 kBodyCount617 = UINT64_C(14307150);
constexpr u64 kSelectedDeleteFnv617 = UINT64_C(0x9e07ae6ed26a77df);
constexpr u64 kExchangedCarrierFnv617 = UINT64_C(0x34282dbb418a1b86);

struct PairAudit617 {
    u64 active = 0;
    u64 active_fnv = 0;
    u64 active_joint = 0;
    u64 active_nonjoint = 0;
    u64 exposed = 0;
    u64 exposed_fnv = 0;
    u64 hit_incidences = 0;
    u64 minimum_hits = std::numeric_limits<u64>::max();
    u64 maximum_hits = 0;
    u64 failures = 0;
    u64 failure_fnv = 0;
    std::vector<u32> failure_bodies;
};

std::vector<u32> read_mixed617(const std::filesystem::path& path,
                               std::size_t expected_count,
                               u64 expected_fnv) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open mixed carrier ledger");
    std::vector<u32> masks;
    std::set<u32> distinct;
    Fnv ledger;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_mask_agent(token);
        require((std::popcount(mask) == 8 || std::popcount(mask) == 9) &&
                    distinct.insert(mask).second,
                "mixed carrier rank/distinctness changed");
        masks.push_back(mask);
        ledger.add(mask);
    }
    require(masks.size() == expected_count && ledger.state == expected_fnv,
            "mixed carrier identity changed");
    return masks;
}

PairAudit617 audit_pair617(
    PairAgent pair, const std::vector<u32>& joint,
    const std::vector<u32>& carrier,
    const std::unordered_set<u32>& joint_set) {
    const Geometry geometry = build_geometry(pair.q, pair.r);
    std::vector<u32> active_joint;
    std::vector<u32> active_nonjoint;
    PairAudit617 audit;
    Fnv active_ledger;
    for (u32 mask : joint)
        if (margin(geometry, mask).ticks >= 0) active_joint.push_back(mask);
    for (u32 mask : carrier) {
        if (margin(geometry, mask).ticks < 0) continue;
        ++audit.active;
        active_ledger.add(mask);
        if (!joint_set.contains(mask)) active_nonjoint.push_back(mask);
    }
    audit.active_fnv = active_ledger.state;
    audit.active_joint = active_joint.size();
    audit.active_nonjoint = active_nonjoint.size();
    require(audit.active == audit.active_joint + audit.active_nonjoint,
            "active joint partition changed");

    Fnv exposed_ledger;
    Fnv failure_ledger;
    u64 bodies = 0;
    const u32 limit = u32{1} << 30;
    for (u32 body = (u32{1} << 9) - 1; body < limit;
         body = next_combination(body)) {
        ++bodies;
        bool joint_hit = false;
        for (u32 mask : active_joint)
            if ((mask & body) == 0) {
                joint_hit = true;
                break;
            }
        if (joint_hit) continue;
        ++audit.exposed;
        exposed_ledger.add(body);
        u64 hits = 0;
        for (u32 mask : active_nonjoint)
            if ((mask & body) == 0) ++hits;
        audit.hit_incidences += hits;
        if (hits == 0) {
            ++audit.failures;
            audit.failure_bodies.push_back(body);
            failure_ledger.add(body);
        } else {
            audit.minimum_hits = std::min(audit.minimum_hits, hits);
            audit.maximum_hits = std::max(audit.maximum_hits, hits);
        }
    }
    require(bodies == kBodyCount617, "body universe changed");
    audit.exposed_fnv = exposed_ledger.state;
    audit.failure_fnv = failure_ledger.state;
    if (audit.exposed == 0) audit.minimum_hits = 0;
    return audit;
}

}  // namespace

#ifndef ENDPOINT617_RAW_VERIFY_MAIN
#define ENDPOINT617_RAW_VERIFY_MAIN main
#endif

int ENDPOINT617_RAW_VERIFY_MAIN(int argc, char** argv) {
    try {
        require(argc == 10,
                "usage: r617_raw_verify JOINT BASE8951 ADD45 SUFFIX9 "
                "RESIDUAL REPAIRS67 DELETE60 EXCHANGED9019 FAILURES");
        const std::vector<u32> joint = read_masks_agent(
            argv[1], kJointCountAgent, kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> repairs = read_repairs617(argv[6]);
        const std::vector<u32> selected = read_masks_agent(
            argv[7], 60, kSelectedDeleteFnv617, 8);
        const std::set<u32> selected_set(selected.begin(), selected.end());
        for (u32 mask : selected)
            require(!joint_set.contains(mask),
                    "selected deletion entered joint deck");

        std::vector<u32> reconstructed =
            build_mixed_carrier(argv[2], argv[3], argv[4]);
        std::set<u32> distinct(reconstructed.begin(), reconstructed.end());
        for (u32 repair : repairs) {
            require(distinct.insert(repair).second, "repair overlap");
            reconstructed.push_back(repair);
        }
        require(reconstructed.size() == 9079 &&
                    masks_fnv_agent(reconstructed) ==
                        kAugmentedCarrierFnv617,
                "augmented carrier changed");
        std::vector<u32> exchanged;
        for (u32 mask : reconstructed)
            if (!selected_set.contains(mask)) exchanged.push_back(mask);
        require(exchanged.size() == 9019 &&
                    masks_fnv_agent(exchanged) == kExchangedCarrierFnv617,
                "reconstructed exchanged carrier changed");
        const std::vector<u32> serialized =
            read_mixed617(argv[8], 9019, kExchangedCarrierFnv617);
        require(serialized == exchanged,
                "serialized/reconstructed carrier order differs");
        for (u32 mask : joint)
            require(std::find(exchanged.begin(), exchanged.end(), mask) !=
                        exchanged.end(),
                    "joint mask absent from exchanged carrier");

        const std::vector<PairAgent> band = read_band_agent(argv[5], 617);
        require(band.size() == 183, "endpoint-617 prefix changed");
        std::ofstream failures_out(argv[9]);
        require(static_cast<bool>(failures_out),
                "cannot create failure ledger");
        failures_out << "q,r,body_hex\n";

        const u64 rank8_count = std::count_if(
            exchanged.begin(), exchanged.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        const u64 rank9_count = exchanged.size() - rank8_count;
        std::cout << "ENDPOINT617_EXCHANGED_CARRIER_RAW_VERIFY_V1\n"
                  << "CARRIER_SIZE " << exchanged.size() << " FNV "
                  << std::hex << masks_fnv_agent(exchanged) << std::dec
                  << " RANK8 " << rank8_count << " RANK9 " << rank9_count
                  << " LOWER 617 SELECTED " << band.size()
                  << " BODY_UNIVERSE_PER_PAIR " << kBodyCount617 << '\n';
        u64 total_failures = 0;
        u64 total_exposed = 0;
        u64 total_hit_incidences = 0;
        Fnv ledger;
        std::size_t index = 0;
        while (index < band.size()) {
            const int endpoint = band[index].r;
            u64 layer_rows = 0;
            u64 layer_failures = 0;
            Fnv layer_ledger;
            while (index < band.size() && band[index].r == endpoint) {
                const PairAgent pair = band[index++];
                const PairAudit617 audit =
                    audit_pair617(pair, joint, exchanged, joint_set);
                for (u32 body : audit.failure_bodies)
                    failures_out << pair.q << ',' << pair.r << ','
                                 << hex8(body) << '\n';
                ++layer_rows;
                layer_failures += audit.failures;
                total_failures += audit.failures;
                total_exposed += audit.exposed;
                total_hit_incidences += audit.hit_incidences;
                layer_ledger.add(pair.q);
                layer_ledger.add(pair.r);
                ledger.add(pair.q);
                ledger.add(pair.r);
                ledger.add(audit.active);
                ledger.add(audit.active_fnv);
                ledger.add(audit.exposed);
                ledger.add(audit.exposed_fnv);
                ledger.add(audit.failures);
                ledger.add(audit.failure_fnv);
                std::cout << "PAIR " << pair.q << ',' << pair.r
                          << " ACTIVE " << audit.active << " ACTIVE_FNV "
                          << std::hex << audit.active_fnv << std::dec
                          << " ACTIVE_JOINT " << audit.active_joint
                          << " ACTIVE_NONJOINT " << audit.active_nonjoint
                          << " EXPOSED " << audit.exposed << " HIT_RANGE "
                          << audit.minimum_hits << ".." << audit.maximum_hits
                          << " HIT_INCIDENCES " << audit.hit_incidences
                          << " FAILURES " << audit.failures
                          << " FAILURE_FNV " << std::hex << audit.failure_fnv
                          << std::dec << '\n';
            }
            std::cout << "LAYER " << endpoint << " ROWS " << layer_rows
                      << " FNV " << std::hex << layer_ledger.state << std::dec
                      << " FAILURES " << layer_failures << '\n';
        }
        require(failures_out.good(), "failure ledger write failed");
        require(total_failures == 0,
                "exchanged carrier has a prefix failure");
        std::cout << "COMPLETED_ROWS " << band.size()
                  << " TOTAL_BODY_TESTS " << kBodyCount617 * band.size()
                  << " TOTAL_EXPOSED " << total_exposed
                  << " TOTAL_HIT_INCIDENCES " << total_hit_incidences
                  << " TOTAL_FAILURES " << total_failures << " LEDGER_FNV "
                  << std::hex << ledger.state << std::dec << '\n'
                  << "JOINT_PARTITION 421_OF_421_RETAINED; FAST_EXPOSED_"
                     "BODY_REPLAY_DEPENDS_ON_INTACT_JOINT_DECK\n"
                  << "SCOPE FINITE_EXACT_FIXED_POOL_COMPLETE_R_GE_617_"
                     "PREFIX_NO_BELOW617_CLAIM_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS EXACT_9019_MASK_EXCHANGED_CARRIER_"
                     "CLOSES_PREFIX\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT617_RAW_VERIFY_ERROR " << error.what() << '\n';
        return 1;
    }
}
