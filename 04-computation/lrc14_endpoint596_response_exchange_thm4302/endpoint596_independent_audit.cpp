// Structurally independent audit for THM-4302.  Responders are generated as
// subsets of the 21-label complements of the 24 obligations rather than by
// a scan of the full mask universes.  The program then independently replays
// all C(30,9) labelled bodies before and after the frozen exchange.

#define ENDPOINT617_RAW_VERIFY_MAIN thm4302_raw_hidden_main
#include "../lrc14_size_preserving_response_staircase_thm4300/endpoint617_exchanged_carrier_raw_verify.cpp"
#undef ENDPOINT617_RAW_VERIFY_MAIN

#include <fstream>
#include <unordered_map>

namespace {

constexpr u64 kRepairFnv = UINT64_C(0x64ce5f9d1ec8c4c2);
constexpr u64 kFailureFnv = UINT64_C(0x3dbd5b39673070ff);
constexpr u64 kAdditionFnv = UINT64_C(0xdc0eebaebf688c65);
constexpr u64 kDeleteFnv = UINT64_C(0x9240b264ab65aa62);
constexpr u64 kInactiveFnv = UINT64_C(0xfa143e58f59119f8);
constexpr u64 kExchangedFnv = UINT64_C(0x892fef44a9e6b37e);

u32 response_independent(u32 mask, const std::vector<u32>& failures) {
    u32 reply = 0;
    for (std::size_t i = 0; i < failures.size(); ++i)
        if ((mask & failures[i]) == 0) reply |= u32{1} << i;
    return reply;
}

template <class Consumer>
void enumerate_complement_masks(u32 body, unsigned rank, Consumer consume) {
    std::array<unsigned, 30> free{};
    unsigned count = 0;
    for (unsigned bit = 0; bit < 30; ++bit)
        if ((body & (u32{1} << bit)) == 0) free[count++] = bit;
    require(count == 21, "body complement size changed");
    const u32 limit = u32{1} << count;
    for (u32 local = (u32{1} << rank) - 1; local < limit;
         local = next_combination(local)) {
        u32 mask = 0;
        for (unsigned index = 0; index < count; ++index)
            if (local & (u32{1} << index)) mask |= u32{1} << free[index];
        consume(mask);
    }
}

u64 sorted_fnv(std::vector<u32> values) {
    std::sort(values.begin(), values.end());
    Fnv ledger;
    for (u32 value : values) ledger.add(value);
    return ledger.state;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 9,
                "usage: independent JOINT BASE8951 ADD45 SUFFIX9 RESIDUAL "
                "REPAIRS76 ADDITIONS4 DELETE73");
        const auto joint = read_masks_agent(
            argv[1], kJointCountAgent, kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const auto repairs = read_mixed617(argv[6], 76, kRepairFnv);
        const auto additions = read_mixed617(argv[7], 4, kAdditionFnv);
        const auto deletions = read_mixed617(argv[8], 73, kDeleteFnv);

        std::vector<u32> augmented =
            build_mixed_carrier(argv[2], argv[3], argv[4]);
        std::set<u32> distinct(augmented.begin(), augmented.end());
        for (u32 repair : repairs) {
            require(distinct.insert(repair).second, "repair overlap");
            augmented.push_back(repair);
        }
        require(augmented.size() == 9088 &&
                    masks_fnv_agent(augmented) ==
                        UINT64_C(0x55e8588798885ae5),
                "augmented carrier changed");
        const std::unordered_set<u32> augmented_set(
            augmented.begin(), augmented.end());

        // Derive the complete failure universe from a fresh labelled-body
        // replay.  No maintained failure CSV is an input to this path.
        const PairAgent pair{210, 596};
        const PairAudit617 before =
            audit_pair617(pair, joint, augmented, joint_set);
        require(before.failures == 24 && before.failure_fnv == kFailureFnv,
                "independent raw failure replay changed");
        const std::vector<u32> failures = before.failure_bodies;

        const Geometry local = build_geometry(210, 596);
        std::unordered_set<u32> responders8;
        std::unordered_set<u32> responders9;
        for (u32 body : failures) {
            enumerate_complement_masks(body, 8, [&](u32 mask) {
                if (margin(local, mask).ticks >= 0) responders8.insert(mask);
            });
            enumerate_complement_masks(body, 9, [&](u32 mask) {
                if (margin(local, mask).ticks >= 0) responders9.insert(mask);
            });
        }
        require(responders8.size() == 9090 && responders9.size() == 138019 &&
                    sorted_fnv(std::vector<u32>(responders8.begin(),
                                                responders8.end())) ==
                        UINT64_C(0x2dddbe3405491cdd) &&
                    sorted_fnv(std::vector<u32>(responders9.begin(),
                                                responders9.end())) ==
                        UINT64_C(0x2202e93c739926df),
                "complement-generated responder universe changed");
        for (u32 mask : responders8)
            require(!augmented_set.contains(mask),
                    "existing rank-eight responder contradicts failure");
        for (u32 mask : responders9)
            require(!augmented_set.contains(mask),
                    "existing rank-nine responder contradicts failure");

        std::unordered_map<u32, std::array<u64, 2>> atlas;
        for (u32 mask : responders8)
            ++atlas[response_independent(mask, failures)][0];
        for (u32 mask : responders9)
            ++atlas[response_independent(mask, failures)][1];
        require(atlas.size() == 718, "response-type count changed");
        u64 rank8_types = 0;
        std::array<u64, 3> mixed_hist{};
        std::array<u64, 3> rank8_hist{};
        constexpr u32 mixed_half = 0x00288146;
        constexpr u32 rank8_half = 0x00d88143;
        constexpr u32 rank8_double = 0x00000004;
        for (const auto& [reply, counts] : atlas) {
            const unsigned mixed_load = std::popcount(reply & mixed_half);
            require(mixed_load <= 2, "mixed dual overload");
            ++mixed_hist[mixed_load];
            if (counts[0] == 0) continue;
            ++rank8_types;
            const unsigned rank8_load =
                std::popcount(reply & rank8_half) +
                2 * std::popcount(reply & rank8_double);
            require(rank8_load <= 2, "rank-eight dual overload");
            ++rank8_hist[rank8_load];
        }
        require(mixed_hist == std::array<u64, 3>{255, 355, 108} &&
                    rank8_types == 220 &&
                    rank8_hist == std::array<u64, 3>{78, 104, 38},
                "independent dual certificates changed");
        u32 addition_cover = 0;
        for (u32 mask : additions) {
            require(std::popcount(mask) == 9 && responders9.contains(mask),
                    "frozen addition is not a responder");
            addition_cover |= response_independent(mask, failures);
        }
        require(addition_cover == 0x00ffffff,
                "four frozen additions fail cover");

        const auto band = read_band_agent(argv[5], 596);
        require(band.size() == 363, "prefix changed");
        std::vector<Geometry> geometries;
        geometries.reserve(band.size());
        for (PairAgent pair : band)
            geometries.push_back(build_geometry(pair.q, pair.r));
        std::vector<u32> inactive;
        for (u32 mask : augmented) {
            bool ever_active = false;
            for (const Geometry& geometry : geometries)
                ever_active |= margin(geometry, mask).ticks >= 0;
            if (!ever_active) inactive.push_back(mask);
        }
        std::sort(inactive.begin(), inactive.end());
        require(inactive.size() == 75 && sorted_fnv(inactive) == kInactiveFnv,
                "independent inactive pool changed");
        require(std::equal(deletions.begin(), deletions.end(), inactive.begin()),
                "frozen deletion is not inactive prefix");

        const std::set<u32> deletion_set(deletions.begin(), deletions.end());
        std::vector<u32> exchanged;
        for (u32 mask : augmented)
            if (!deletion_set.contains(mask)) exchanged.push_back(mask);
        for (u32 mask : additions) exchanged.push_back(mask);
        require(exchanged.size() == 9019 &&
                    masks_fnv_agent(exchanged) == kExchangedFnv,
                "independent exchange identity changed");
        for (u32 mask : joint)
            require(std::find(exchanged.begin(), exchanged.end(), mask) !=
                        exchanged.end(),
                    "exchange deleted a joint coordinate");

        std::vector<u32> augmented_plus = augmented;
        augmented_plus.insert(augmented_plus.end(), additions.begin(),
                              additions.end());
        const PairAudit617 after_additions =
            audit_pair617(pair, joint, augmented_plus, joint_set);
        const PairAudit617 after_exchange =
            audit_pair617(pair, joint, exchanged, joint_set);
        require(after_additions.failures == 0 &&
                    after_additions.minimum_hits == 1 &&
                    after_additions.maximum_hits == 254 &&
                    after_exchange.failures == 0 &&
                    after_exchange.minimum_hits == 1 &&
                    after_exchange.maximum_hits == 254,
                "independent repaired raw replay changed");

        Fnv deletion_sign_ledger;
        u64 equality_cells = 0;
        for (u32 mask : deletions) {
            require(!joint_set.contains(mask), "deleted joint coordinate");
            for (const Geometry& geometry : geometries) {
                const Margin exact = margin(geometry, mask);
                require(exact.ticks < 0, "deleted mask active in prefix");
                equality_cells += exact.ticks == 0;
                deletion_sign_ledger.add(mask);
                deletion_sign_ledger.add(false);
            }
        }
        require(equality_cells == 0 && deletion_sign_ledger.state ==
                    UINT64_C(0xee7857ffb11111b2),
                "deletion sign ledger changed");

        std::cout << "THM4302_ENDPOINT596_INDEPENDENT_AUDIT_V1\n"
                  << "RAW_REPLAY BODY_UNIVERSE 14307150 "
                     "NO_FAILURE_CSV_INPUT YES\n"
                  << "ENUMERATION complement_subsets_of_24_bodies\n"
                  << "RESPONDERS_RANK8 " << responders8.size() << " FNV "
                  << std::hex << UINT64_C(0x2dddbe3405491cdd) << std::dec
                  << " RESPONDERS_RANK9 " << responders9.size() << " FNV "
                  << std::hex << UINT64_C(0x2202e93c739926df) << std::dec
                  << '\n'
                  << "RESPONSE_TYPES " << atlas.size()
                  << " MIXED_DUAL_VALUE 7/2 RANK8_DUAL_VALUE 11/2\n"
                  << "BEFORE ACTIVE " << before.active << " JOINT "
                  << before.active_joint << " NONJOINT "
                  << before.active_nonjoint << " EXPOSED " << before.exposed
                  << " FAILURES " << before.failures << " FNV " << std::hex
                  << before.failure_fnv << std::dec << '\n'
                  << "AFTER_ADDITIONS ACTIVE " << after_additions.active
                  << " JOINT " << after_additions.active_joint << " NONJOINT "
                  << after_additions.active_nonjoint << " EXPOSED "
                  << after_additions.exposed << " HIT_RANGE "
                  << after_additions.minimum_hits << ".."
                  << after_additions.maximum_hits << " FAILURES "
                  << after_additions.failures << '\n'
                  << "AFTER_EXCHANGE ACTIVE " << after_exchange.active
                  << " JOINT " << after_exchange.active_joint << " NONJOINT "
                  << after_exchange.active_nonjoint << " EXPOSED "
                  << after_exchange.exposed << " HIT_RANGE "
                  << after_exchange.minimum_hits << ".."
                  << after_exchange.maximum_hits << " FAILURES "
                  << after_exchange.failures
                  << '\n'
                  << "INACTIVE 75 FNV " << std::hex << kInactiveFnv
                  << " DELETE 73 FNV " << kDeleteFnv << " ADD 4 FNV "
                  << kAdditionFnv << std::dec << '\n'
                  << "DELETE_SIGN_CELLS " << deletions.size() * band.size()
                  << " SIGN_FNV " << std::hex << deletion_sign_ledger.state
                  << std::dec << " EQUALITIES " << equality_cells << '\n'
                  << "EXCHANGED 9019 FNV " << std::hex << kExchangedFnv
                  << std::dec << '\n'
                  << "SCOPE INDEPENDENT_FINITE_EXACT_CARRIER_AUDIT_"
                     "NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4302_INDEPENDENT_AUDIT_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}
