// Exact sign-only exchange audit for THM-4300.  The inherited literal-wall
// implementation supplies the fixed-pool geometry and the 9,019-mask carrier.

#define ENDPOINT626_EXCHANGE_MAIN thm4296_endpoint626_unused_main
#include "endpoint_exchange_primitives.cpp"
#undef ENDPOINT626_EXCHANGE_MAIN

#include <fstream>

namespace {

constexpr std::size_t kRepairCount4300 = 76;
constexpr u64 kRepairFnv4300 = UINT64_C(0x64ce5f9d1ec8c4c2);
constexpr std::size_t kDeleteCount4300 = 69;
constexpr u64 kDeleteFnv4300 = UINT64_C(0xce08bd348c1ac6c7);
constexpr u64 kAugmentedFnv4300 = UINT64_C(0x55e8588798885ae5);
constexpr u64 kExchangedFnv4300 = UINT64_C(0xe0fcd06628e1aa37);
constexpr u64 kOriginalFnv4300 = UINT64_C(0xd7f0e06e154e78c2);
constexpr u64 kInactiveFnv4300 = UINT64_C(0x70cf014b73f82b0e);
constexpr u64 kSignFnv4300 = UINT64_C(0x31605c3bc8dc5311);

constexpr std::array<u32, 7> kInheritedRepairs4300 = {
    0x0010e125, 0x002ac4c0, 0x3882a082, 0x0041c325,
    0x08c28e40, 0x02008327, 0x0006e281};

std::vector<u32> read_mixed4300(const std::filesystem::path& path,
                                std::size_t expected_count,
                                u64 expected_fnv) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open mixed-mask ledger");
    std::vector<u32> masks;
    std::set<u32> distinct;
    Fnv ledger;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_mask_agent(token);
        require((std::popcount(mask) == 8 || std::popcount(mask) == 9) &&
                    distinct.insert(mask).second,
                "mixed-mask rank/distinctness changed");
        masks.push_back(mask);
        ledger.add(mask);
    }
    require(input.eof() && masks.size() == expected_count &&
                ledger.state == expected_fnv,
            "mixed-mask identity changed");
    return masks;
}

void write_masks4300(const std::filesystem::path& path,
                     const std::vector<u32>& masks) {
    std::ofstream output(path, std::ios::binary);
    require(static_cast<bool>(output), "cannot create mask ledger");
    for (u32 mask : masks) output << hex8(mask) << '\n';
    require(output.good(), "mask-ledger write failed");
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 10,
                "usage: audit JOINT BASE8951 ADD45 SUFFIX9 RESIDUAL "
                "REPAIRS76 INACTIVE_OUT DELETE_OUT EXCHANGED_OUT");
        const std::vector<u32> joint = read_masks_agent(
            argv[1], kJointCountAgent, kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> repairs =
            read_mixed4300(argv[6], kRepairCount4300, kRepairFnv4300);
        for (std::size_t i = 0; i < kInheritedRepairs4300.size(); ++i)
            require(repairs[i] == kInheritedRepairs4300[i],
                    "inherited repair prefix changed");

        std::vector<u32> original =
            build_mixed_carrier(argv[2], argv[3], argv[4]);
        std::set<u32> distinct(original.begin(), original.end());
        for (u32 mask : kInheritedRepairs4300) {
            require(distinct.insert(mask).second,
                    "inherited repair overlaps carrier");
            original.push_back(mask);
        }
        require(original.size() == 9019 &&
                    masks_fnv_agent(original) == kOriginalFnv4300,
                "inherited carrier identity changed");

        std::vector<u32> augmented = original;
        for (std::size_t i = kInheritedRepairs4300.size();
             i < repairs.size(); ++i) {
            require(distinct.insert(repairs[i]).second,
                    "new response overlaps carrier");
            augmented.push_back(repairs[i]);
        }
        require(augmented.size() == 9088 &&
                    masks_fnv_agent(augmented) == kAugmentedFnv4300,
                "augmented carrier identity changed");

        const std::vector<PairAgent> band = read_band_agent(argv[5], 597);
        require(band.size() == 354,
                "complete endpoint-at-least-597 prefix changed");
        std::vector<unsigned char> ever_active(original.size(), 0);
        u64 equalities = 0;
        Fnv sign_ledger;
        for (PairAgent pair : band) {
            const Geometry geometry = build_geometry(pair.q, pair.r);
            for (std::size_t i = 0; i < original.size(); ++i) {
                const Margin exact = margin(geometry, original[i]);
                const bool active = exact.ticks >= 0;
                equalities += exact.ticks == 0;
                ever_active[i] |= active;
                sign_ledger.add(original[i]);
                sign_ledger.add(active);
            }
        }
        require(equalities == 0 && sign_ledger.state == kSignFnv4300,
                "literal sign ledger changed");

        std::vector<u32> inactive;
        for (std::size_t i = 0; i < original.size(); ++i)
            if (!ever_active[i]) inactive.push_back(original[i]);
        std::sort(inactive.begin(), inactive.end());
        Fnv inactive_ledger;
        for (u32 mask : inactive) {
            require(!joint_set.contains(mask),
                    "common-inactive mask entered retained joint deck");
            require(std::popcount(mask) == 8,
                    "common-inactive mask rank changed");
            inactive_ledger.add(mask);
        }
        require(inactive.size() == 76 &&
                    inactive_ledger.state == kInactiveFnv4300,
                "common-inactive pool changed");
        write_masks4300(argv[7], inactive);

        const std::vector<u32> selected(inactive.begin(),
                                        inactive.begin() + kDeleteCount4300);
        Fnv selected_ledger;
        for (u32 mask : selected) selected_ledger.add(mask);
        require(selected_ledger.state == kDeleteFnv4300,
                "lexicographic deletion ledger changed");
        write_masks4300(argv[8], selected);

        const std::set<u32> selected_set(selected.begin(), selected.end());
        std::vector<u32> exchanged;
        for (u32 mask : augmented)
            if (!selected_set.contains(mask)) exchanged.push_back(mask);
        require(exchanged.size() == 9019 &&
                    masks_fnv_agent(exchanged) == kExchangedFnv4300,
                "size-preserving carrier identity changed");
        for (u32 mask : joint)
            require(std::find(exchanged.begin(), exchanged.end(), mask) !=
                        exchanged.end(),
                    "exchange deleted a joint-deck coordinate");
        write_masks4300(argv[9], exchanged);

        const u64 rank8 = std::count_if(
            exchanged.begin(), exchanged.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        std::cout << "THM4300_ENDPOINT597_EXCHANGE_AUDIT_V1\n"
                  << "ORIGINAL 9019 FNV " << std::hex << kOriginalFnv4300
                  << std::dec << " PREFIX_ROWS " << band.size()
                  << " SIGN_CELLS " << original.size() * band.size()
                  << " SIGN_FNV " << std::hex << sign_ledger.state
                  << std::dec << " EQUALITIES " << equalities << '\n'
                  << "REPAIRS " << repairs.size() << " FNV " << std::hex
                  << kRepairFnv4300 << " AUGMENTED 9088 FNV "
                  << kAugmentedFnv4300 << std::dec << '\n'
                  << "COMMON_INACTIVE " << inactive.size() << " FNV "
                  << std::hex << inactive_ledger.state << std::dec
                  << " NONJOINT " << inactive.size() << " RANK8 "
                  << inactive.size() << " RANK9 0\n"
                  << "SELECTED_DELETE " << selected.size() << " FNV "
                  << std::hex << selected_ledger.state << std::dec << '\n'
                  << "EXCHANGED " << exchanged.size() << " FNV "
                  << std::hex << masks_fnv_agent(exchanged) << std::dec
                  << " RANK8 " << rank8 << " RANK9 "
                  << exchanged.size() - rank8 << '\n'
                  << "TRANSFER deleted masks are inactive on every prefix "
                     "row, so deleting them preserves each active subfamily "
                     "and every body-cover consequence of the augmented "
                     "carrier\n"
                  << "SCOPE FINITE_EXACT_FIXED_POOL_COMPLETE_R_GE_597_"
                     "PREFIX_NO_BELOW597_CLAIM_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS EXACT_SIZE_PRESERVING_EXCHANGE\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4300_EXCHANGE_AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}
