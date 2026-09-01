// Scratch exact sign audit for a size-preserving exchange on the complete
// THM-4296 residual prefix with endpoint at least 617.  This deliberately
// consumes only the sign of each literal margin at each pair: raw margins
// are never ordered across distinct pairs (MISTAKE-532).

#define ENDPOINT626_EXCHANGE_MAIN endpoint626_exchange_scan_hidden_main
#include "endpoint_exchange_primitives.cpp"
#undef ENDPOINT626_EXCHANGE_MAIN

#include <fstream>

namespace {

constexpr std::size_t kRepairCount617 = 67;
constexpr std::size_t kDeleteCount617 = 60;
constexpr u64 kAugmentedCarrierFnv617 = UINT64_C(0xbfc37d5b0cbb6744);

std::vector<u32> read_repairs617(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open pinned repair ledger");
    std::vector<u32> repairs;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_mask_agent(token);
        require((std::popcount(mask) == 8 || std::popcount(mask) == 9) &&
                    distinct.insert(mask).second,
                "repair rank/distinctness changed");
        repairs.push_back(mask);
    }
    require(repairs.size() == kRepairCount617,
            "pinned repair count changed");
    return repairs;
}

void write_masks617(const std::filesystem::path& path,
                    const std::vector<u32>& masks) {
    std::ofstream output(path);
    require(static_cast<bool>(output), "cannot create mask ledger");
    for (u32 mask : masks) output << hex8(mask) << '\n';
    require(output.good(), "mask ledger write failed");
}

}  // namespace

#ifndef ENDPOINT617_EXCHANGE_MAIN
#define ENDPOINT617_EXCHANGE_MAIN main
#endif

int ENDPOINT617_EXCHANGE_MAIN(int argc, char** argv) {
    try {
        require(argc == 9,
                "usage: r617_exchange_scan JOINT BASE8951 ADD45 SUFFIX9 "
                "RESIDUAL REPAIRS67 DELETE60_OUT EXCHANGED9019_OUT");
        const std::vector<u32> joint = read_masks_agent(
            argv[1], kJointCountAgent, kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> repairs = read_repairs617(argv[6]);
        Fnv repair_ledger;
        for (u32 mask : repairs) repair_ledger.add(mask);

        std::vector<u32> original =
            build_mixed_carrier(argv[2], argv[3], argv[4]);
        std::set<u32> distinct(original.begin(), original.end());
        for (std::size_t index = 0; index < kPost632Repairs.size(); ++index) {
            require(repairs[index] == kPost632Repairs[index],
                    "inherited seven-mask prefix changed");
            require(distinct.insert(repairs[index]).second,
                    "inherited repair overlap");
            original.push_back(repairs[index]);
        }
        require(original.size() == 9019 &&
                    masks_fnv_agent(original) ==
                        UINT64_C(0xd7f0e06e154e78c2),
                "original 9019-mask carrier identity changed");
        for (u32 mask : joint)
            require(std::find(original.begin(), original.end(), mask) !=
                        original.end(),
                    "joint mask absent from original carrier");

        std::vector<u32> augmented = original;
        for (std::size_t index = kPost632Repairs.size();
             index < repairs.size(); ++index) {
            require(distinct.insert(repairs[index]).second,
                    "new repair overlap");
            augmented.push_back(repairs[index]);
        }
        require(augmented.size() == 9079 &&
                    masks_fnv_agent(augmented) == kAugmentedCarrierFnv617,
                "augmented 9079-mask carrier identity changed");

        const std::vector<PairAgent> band = read_band_agent(argv[5], 617);
        require(band.size() == 183, "endpoint-617 prefix changed");
        std::vector<unsigned char> ever_active(original.size(), 0);
        std::vector<u64> active_rows(original.size(), 0);
        u64 equality_cells = 0;
        Fnv sign_ledger;
        for (PairAgent pair : band) {
            const Geometry geometry = build_geometry(pair.q, pair.r);
            for (std::size_t index = 0; index < original.size(); ++index) {
                const Margin value = margin(geometry, original[index]);
                const bool active = value.ticks >= 0;
                equality_cells += value.ticks == 0;
                ever_active[index] |= active;
                active_rows[index] += active;
                sign_ledger.add(original[index]);
                sign_ledger.add(active);
            }
        }

        std::vector<u32> inactive_all;
        std::vector<u32> inactive_nonjoint;
        std::vector<u64> row_histogram(band.size() + 1, 0);
        for (std::size_t index = 0; index < original.size(); ++index) {
            require(active_rows[index] <= band.size(),
                    "active-row count escaped");
            ++row_histogram[active_rows[index]];
            if (ever_active[index]) continue;
            inactive_all.push_back(original[index]);
            if (!joint_set.contains(original[index]))
                inactive_nonjoint.push_back(original[index]);
        }
        std::sort(inactive_all.begin(), inactive_all.end());
        std::sort(inactive_nonjoint.begin(), inactive_nonjoint.end());
        require(inactive_nonjoint.size() >= kDeleteCount617,
                "fewer than sixty original nonjoint masks are inactive");
        std::vector<u32> selected(inactive_nonjoint.begin(),
                                  inactive_nonjoint.begin() +
                                      kDeleteCount617);
        const std::set<u32> selected_set(selected.begin(), selected.end());
        for (u32 mask : selected)
            require(!joint_set.contains(mask),
                    "selected deletion entered joint deck");
        std::vector<u32> exchanged;
        for (u32 mask : augmented)
            if (!selected_set.contains(mask)) exchanged.push_back(mask);
        require(exchanged.size() == 9019,
                "size-preserving exchange changed size");
        for (u32 mask : joint)
            require(std::find(exchanged.begin(), exchanged.end(), mask) !=
                        exchanged.end(),
                    "exchange deleted a joint mask");
        write_masks617(argv[7], selected);
        write_masks617(argv[8], exchanged);

        Fnv inactive_all_ledger;
        for (u32 mask : inactive_all) inactive_all_ledger.add(mask);
        Fnv inactive_nonjoint_ledger;
        for (u32 mask : inactive_nonjoint)
            inactive_nonjoint_ledger.add(mask);
        Fnv selected_ledger;
        for (u32 mask : selected) selected_ledger.add(mask);
        const u64 inactive_rank8 = std::count_if(
            inactive_nonjoint.begin(), inactive_nonjoint.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        const u64 inactive_rank9 = inactive_nonjoint.size() - inactive_rank8;
        const u64 exchanged_rank8 = std::count_if(
            exchanged.begin(), exchanged.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        const u64 exchanged_rank9 = exchanged.size() - exchanged_rank8;

        std::cout << "ENDPOINT617_SIZE_PRESERVING_EXCHANGE_SCAN_V1\n"
                  << "REPAIRS " << repairs.size() << " FNV " << std::hex
                  << repair_ledger.state << std::dec << " AUGMENTED "
                  << augmented.size() << " FNV " << std::hex
                  << masks_fnv_agent(augmented) << std::dec << '\n'
                  << "ORIGINAL " << original.size() << " FNV " << std::hex
                  << masks_fnv_agent(original) << std::dec << " PAIRS "
                  << band.size() << " SIGN_CELLS "
                  << original.size() * band.size() << " SIGN_FNV " << std::hex
                  << sign_ledger.state << std::dec << " EQUALITIES "
                  << equality_cells << '\n'
                  << "INACTIVE_ALL " << inactive_all.size() << " FNV "
                  << std::hex << inactive_all_ledger.state << std::dec
                  << " INACTIVE_NONJOINT " << inactive_nonjoint.size()
                  << " FNV " << std::hex << inactive_nonjoint_ledger.state
                  << std::dec << " RANK8 " << inactive_rank8 << " RANK9 "
                  << inactive_rank9 << '\n'
                  << "ACTIVE_ROW_HISTOGRAM";
        for (std::size_t rows = 0; rows < row_histogram.size(); ++rows)
            if (row_histogram[rows] != 0)
                std::cout << ' ' << rows << ':' << row_histogram[rows];
        std::cout << '\n' << "SELECTED_DELETE";
        for (u32 mask : selected) std::cout << ' ' << hex8(mask);
        std::cout << " FNV " << std::hex << selected_ledger.state << std::dec
                  << '\n'
                  << "RESULT_SIZE " << exchanged.size() << " RESULT_FNV "
                  << std::hex << masks_fnv_agent(exchanged) << std::dec
                  << " RANK8 " << exchanged_rank8 << " RANK9 "
                  << exchanged_rank9 << '\n'
                  << "JOINT_PARTITION 421_OF_421_RETAINED; DELETIONS_"
                     "NONJOINT_ONLY; FAST_EXPOSED_BODY_REPLAY_REQUIRES_"
                     "THIS_INTACT_JOINT_DECK\n"
                  << "SCOPE EXACT_LITERAL_SIGN_SCAN_COMPLETE_R_GE_617_"
                     "PREFIX_NO_CROSS_GRID_MARGIN_ORDER_NO_BELOW617_CLAIM_"
                     "NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS EXACT_60_MASK_SIZE_PRESERVING_EXCHANGE_"
                     "POOL\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT617_EXCHANGE_SCAN_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}
