// Serialization adapter for the already-canonical THM-4318 carrier.
// This file deliberately does no endpoint-588 auditing: the independent
// audit consumes the resulting checked mask ledger as an inherited input.

#define ENDPOINT590_QUOTIENT_MAIN endpoint590_quotient_hidden_for_export
#include "04-computation/lrc14_endpoint590_exact_nine_response_exchange_thm4318/endpoint590_protected_deletion_quotient.cpp"
#undef ENDPOINT590_QUOTIENT_MAIN

namespace {
std::vector<u32> read_delete9_export(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open delete-nine ledger");
    std::vector<u32> result;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_mask_agent(token);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "delete-nine rank/distinctness changed");
        result.push_back(mask);
    }
    require(input.eof() && result.size() == 9,
            "delete-nine ledger size changed");
    return result;
}
}

int main(int argc, char** argv) {
    try {
        require(argc == 19,
                "usage: export JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE OLD_UNION "
                "REPAIRS76 ADD4 DELETE73 ADD10 FINAL_DELETE ADD593 DELETE593 "
                "COVER43 DELETE43 COVER9 DELETE9 OUTPUT");
        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> before =
            reconstruct_final_carrier_quotient(argv, joint);
        const std::unordered_set<u32> before_set(before.begin(), before.end());
        const auto additions = read_cover_quotient(argv[16]);
        const auto deletions = read_delete9_export(argv[17]);
        const std::set<u32> deletion_set(deletions.begin(), deletions.end());
        std::vector<u32> after;
        for (u32 mask : before) {
            if (deletion_set.contains(mask)) {
                require(!joint_set.contains(mask), "deleted a joint mask");
            } else {
                after.push_back(mask);
            }
        }
        require(after.size() + deletions.size() == before.size(),
                "delete-nine mask absent from inherited carrier");
        for (u32 mask : additions) {
            require(!before_set.contains(mask) && !joint_set.contains(mask),
                    "addition overlaps inherited carrier or joint deck");
            after.push_back(mask);
        }
        const auto rank8 = std::count_if(after.begin(), after.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        require(after.size() == 3925 && rank8 == 3809 &&
                    masks_fnv_agent(after) == UINT64_C(0xeeae5518d84ccac5) &&
                    std::set<u32>(after.begin(), after.end()).size() == after.size(),
                "THM-4318 carrier identity changed");
        std::ofstream output(argv[18], std::ios::binary);
        require(static_cast<bool>(output), "cannot create carrier ledger");
        for (u32 mask : after) output << hex8(mask) << '\n';
        require(output.good(), "carrier ledger write failed");
        std::cout << "INHERITED_THM4318_CARRIER_EXPORT PASS\n"
                  << "MASKS " << after.size() << " RANK8 " << rank8
                  << " RANK9 " << after.size() - rank8 << " FNV "
                  << std::hex << masks_fnv_agent(after) << std::dec << '\n'
                  << "JOINT_RETAINED " << joint.size() << '\n'
                  << "SCOPE SERIALIZATION_ONLY_NO_ENDPOINT588_CLAIM\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "CARRIER_EXPORT_ERROR " << error.what() << '\n';
        return 1;
    }
}
