// Activity support census for THM-4302's C596 carrier on the exact
// 391-row carrier-preservation prefix K596 union the residual r=595 layer.

#define ENDPOINT617_RAW_VERIFY_MAIN endpoint595_activity_hidden_main
#include "04-computation/lrc14_size_preserving_response_staircase_thm4300/endpoint617_exchanged_carrier_raw_verify.cpp"
#undef ENDPOINT617_RAW_VERIFY_MAIN

#include <fstream>

namespace {
constexpr u64 kRepairFnv = UINT64_C(0x64ce5f9d1ec8c4c2);
constexpr u64 kAdditionFnv = UINT64_C(0xdc0eebaebf688c65);
constexpr u64 kDeleteFnv = UINT64_C(0x9240b264ab65aa62);
constexpr u64 kAugmentedFnv = UINT64_C(0x55e8588798885ae5);
constexpr u64 kCarrierFnv = UINT64_C(0x892fef44a9e6b37e);
constexpr u64 kOldUnionFnv = UINT64_C(0x11414a33ab91fef6);

std::set<std::pair<int, int>> read_pair_set_activity(
    const std::filesystem::path& path, std::size_t expected, u64 expected_fnv) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open typed row set");
    std::set<std::pair<int, int>> rows;
    Fnv ledger;
    std::string line;
    while (std::getline(input, line)) {
        const auto comma = line.find(',');
        require(comma != std::string::npos, "malformed typed row");
        const int q = std::stoi(line.substr(0, comma));
        const int r = std::stoi(line.substr(comma + 1));
        require(q > 0 && q < r && rows.insert({q, r}).second,
                "typed row invalid/duplicate");
        ledger.add(q);
        ledger.add(r);
    }
    require(input.eof() && rows.size() == expected && ledger.state == expected_fnv,
            "typed row identity changed");
    return rows;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 11,
                "usage: activity JOINT BASE8951 ADD45 SUFFIX9 RESIDUAL "
                "OLD_UNION REPAIRS76 ADDITIONS4 DELETE73 OUTPUT_CSV");
        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const auto repairs = read_mixed617(argv[7], 76, kRepairFnv);
        const auto additions = read_mixed617(argv[8], 4, kAdditionFnv);
        const auto deletions = read_mixed617(argv[9], 73, kDeleteFnv);
        std::vector<u32> augmented = build_mixed_carrier(argv[2], argv[3], argv[4]);
        std::set<u32> distinct(augmented.begin(), augmented.end());
        for (u32 mask : repairs) {
            require(distinct.insert(mask).second, "repair overlap");
            augmented.push_back(mask);
        }
        require(augmented.size() == 9088 && masks_fnv_agent(augmented) == kAugmentedFnv,
                "augmented carrier changed");
        const std::set<u32> deleted(deletions.begin(), deletions.end());
        std::vector<u32> carrier;
        for (u32 mask : augmented)
            if (!deleted.contains(mask)) carrier.push_back(mask);
        for (u32 mask : additions) carrier.push_back(mask);
        require(carrier.size() == 9019 && masks_fnv_agent(carrier) == kCarrierFnv,
                "C596 carrier changed");

        const auto old_union = read_pair_set_activity(argv[6], 1624, kOldUnionFnv);
        const auto raw_band = read_band_agent(argv[5], 595);
        require(raw_band.size() == 394, "raw band changed");
        std::vector<PairAgent> band;
        for (PairAgent pair : raw_band)
            if (pair.r >= 596 || !old_union.contains({pair.q, pair.r}))
                band.push_back(pair);
        require(band.size() == 391, "391-row prefix changed");
        std::vector<Geometry> geometries;
        for (PairAgent pair : band) geometries.push_back(build_geometry(pair.q, pair.r));

        std::ofstream output(argv[10]);
        require(static_cast<bool>(output), "cannot create activity output");
        output << "mask_hex,rank,joint,active_rows,first_q,first_r,row_support_fnv\n";
        std::array<u64, 392> histogram{};
        Fnv full_ledger;
        for (u32 mask : carrier) {
            u64 count = 0;
            PairAgent first{};
            Fnv support;
            for (std::size_t index = 0; index < band.size(); ++index) {
                const bool active = margin(geometries[index], mask).ticks >= 0;
                if (!active) continue;
                if (count == 0) first = band[index];
                ++count;
                support.add(band[index].q);
                support.add(band[index].r);
            }
            ++histogram[count];
            full_ledger.add(mask);
            full_ledger.add(count);
            full_ledger.add(support.state);
            output << hex8(mask) << ',' << std::popcount(mask) << ','
                   << joint_set.contains(mask) << ',' << count << ',';
            if (count) output << first.q << ',' << first.r;
            else output << ',';
            output << ',' << std::hex << support.state << std::dec << '\n';
        }
        require(output.good(), "activity output write failed");
        std::cout << "LRC14_ENDPOINT595_C596_ACTIVITY_SUPPORT_V1\n"
                  << "PREFIX_ROWS 391 CARRIER 9019 SIGN_CELLS "
                  << carrier.size() * band.size() << '\n'
                  << "SUPPORT_HIST";
        for (std::size_t count = 0; count < histogram.size(); ++count)
            if (histogram[count]) std::cout << ' ' << count << ':' << histogram[count];
        std::cout << "\nLEDGER_FNV " << std::hex << full_ledger.state << std::dec
                  << "\nVERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ACTIVITY_CANDIDATES_ERROR " << error.what() << '\n';
        return 1;
    }
}
