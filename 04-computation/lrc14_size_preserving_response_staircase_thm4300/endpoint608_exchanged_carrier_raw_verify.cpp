// Independent all-body replay of the 9,019-mask exchange on all 272 current
// residual rows with endpoint at least 608.

#define ENDPOINT617_RAW_VERIFY_MAIN endpoint617_raw_hidden_for_r608
#include "endpoint617_exchanged_carrier_raw_verify.cpp"
#undef ENDPOINT617_RAW_VERIFY_MAIN

namespace {

#ifndef ENDPOINT608_RAW_REPAIR_COUNT
#define ENDPOINT608_RAW_REPAIR_COUNT 72
#endif
#ifndef ENDPOINT608_RAW_REPAIR_FNV
#define ENDPOINT608_RAW_REPAIR_FNV UINT64_C(0x0217ad40e4b94b2f)
#endif
#ifndef ENDPOINT608_RAW_DELETE_COUNT
#define ENDPOINT608_RAW_DELETE_COUNT 65
#endif
#ifndef ENDPOINT608_RAW_DELETE_FNV
#define ENDPOINT608_RAW_DELETE_FNV UINT64_C(0x6445c29cb1bf9526)
#endif
#ifndef ENDPOINT608_RAW_AUGMENTED_SIZE
#define ENDPOINT608_RAW_AUGMENTED_SIZE 9084
#endif
#ifndef ENDPOINT608_RAW_AUGMENTED_FNV
#define ENDPOINT608_RAW_AUGMENTED_FNV UINT64_C(0x92fe8eeb7a6a69e8)
#endif
#ifndef ENDPOINT608_RAW_EXCHANGED_FNV
#define ENDPOINT608_RAW_EXCHANGED_FNV UINT64_C(0x2b05a51d92ec8983)
#endif

constexpr u64 kRepair72Fnv608 = ENDPOINT608_RAW_REPAIR_FNV;
constexpr u64 kDelete65Fnv608 = ENDPOINT608_RAW_DELETE_FNV;
constexpr u64 kAugmentedFnv608 = ENDPOINT608_RAW_AUGMENTED_FNV;
constexpr u64 kExchangedFnv608 = ENDPOINT608_RAW_EXCHANGED_FNV;

}  // namespace

#ifndef ENDPOINT608_RAW_VERIFY_MAIN
#define ENDPOINT608_RAW_VERIFY_MAIN main
#endif
#ifndef ENDPOINT608_RAW_READ_LOWER
#define ENDPOINT608_RAW_READ_LOWER 608
#endif
#ifndef ENDPOINT608_RAW_MIN_ENDPOINT
#define ENDPOINT608_RAW_MIN_ENDPOINT 608
#endif
#ifndef ENDPOINT608_RAW_MAX_ENDPOINT
#define ENDPOINT608_RAW_MAX_ENDPOINT 1000000
#endif
#ifndef ENDPOINT608_RAW_PAIR_COUNT
#define ENDPOINT608_RAW_PAIR_COUNT 272
#endif
#ifndef ENDPOINT608_RAW_HEADER
#define ENDPOINT608_RAW_HEADER "ENDPOINT608_EXCHANGED_CARRIER_RAW_VERIFY_V1"
#endif
#ifndef ENDPOINT608_RAW_RANGE_TEXT
#define ENDPOINT608_RAW_RANGE_TEXT " LOWER 608 SELECTED "
#endif
#ifndef ENDPOINT608_RAW_SCOPE
#define ENDPOINT608_RAW_SCOPE \
    "SCOPE FINITE_EXACT_FIXED_POOL_COMPLETE_R_GE_608_PREFIX_NO_BELOW608_" \
    "CLAIM_NO_PHYSICAL_ENTRY_NO_LRC14"
#endif
#ifndef ENDPOINT608_RAW_INITIAL_LEDGER
#define ENDPOINT608_RAW_INITIAL_LEDGER kEmptyFnv
#endif

int ENDPOINT608_RAW_VERIFY_MAIN(int argc, char** argv) {
    try {
        require(argc == 10,
                "usage: r608_raw_verify JOINT BASE8951 ADD45 SUFFIX9 "
                "RESIDUAL REPAIRS72 DELETE65 EXCHANGED9019 FAILURES");
        const std::vector<u32> joint = read_masks_agent(
            argv[1], kJointCountAgent, kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> repairs =
            read_mixed617(argv[6], ENDPOINT608_RAW_REPAIR_COUNT,
                          kRepair72Fnv608);
        const std::vector<u32> selected = read_masks_agent(
            argv[7], ENDPOINT608_RAW_DELETE_COUNT, kDelete65Fnv608, 8);
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
        require(reconstructed.size() == ENDPOINT608_RAW_AUGMENTED_SIZE &&
                    masks_fnv_agent(reconstructed) == kAugmentedFnv608,
                "augmented carrier changed");
        std::vector<u32> exchanged;
        for (u32 mask : reconstructed)
            if (!selected_set.contains(mask)) exchanged.push_back(mask);
        require(exchanged.size() == 9019 &&
                    masks_fnv_agent(exchanged) == kExchangedFnv608,
                "reconstructed exchanged carrier changed");
        const std::vector<u32> serialized =
            read_mixed617(argv[8], 9019, kExchangedFnv608);
        require(serialized == exchanged,
                "serialized/reconstructed carrier order differs");
        for (u32 mask : joint)
            require(std::find(exchanged.begin(), exchanged.end(), mask) !=
                        exchanged.end(),
                    "joint mask absent from exchanged carrier");

        std::vector<PairAgent> band =
            read_band_agent(argv[5], ENDPOINT608_RAW_READ_LOWER);
        std::erase_if(band, [](PairAgent pair) {
            return pair.r < ENDPOINT608_RAW_MIN_ENDPOINT ||
                   pair.r > ENDPOINT608_RAW_MAX_ENDPOINT;
        });
        require(band.size() == ENDPOINT608_RAW_PAIR_COUNT,
                "selected raw-replay band changed");
        std::ofstream failures_out(argv[9]);
        require(static_cast<bool>(failures_out),
                "cannot create failure ledger");
        failures_out << "q,r,body_hex\n";
        const u64 rank8_count = std::count_if(
            exchanged.begin(), exchanged.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        const u64 rank9_count = exchanged.size() - rank8_count;
        std::cout << ENDPOINT608_RAW_HEADER << '\n'
                  << "CARRIER_SIZE " << exchanged.size() << " FNV "
                  << std::hex << masks_fnv_agent(exchanged) << std::dec
                  << " RANK8 " << rank8_count << " RANK9 " << rank9_count
                  << ENDPOINT608_RAW_RANGE_TEXT << band.size()
                  << " BODY_UNIVERSE_PER_PAIR " << kBodyCount617 << '\n';

        u64 total_failures = 0;
        u64 total_exposed = 0;
        u64 total_hit_incidences = 0;
        Fnv ledger;
        ledger.state = ENDPOINT608_RAW_INITIAL_LEDGER;
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
                  << ENDPOINT608_RAW_SCOPE << '\n'
                  << "VERDICT PASS EXACT_9019_MASK_EXCHANGED_CARRIER_"
                     "CLOSES_PREFIX\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT608_RAW_VERIFY_ERROR " << error.what() << '\n';
        return 1;
    }
}
