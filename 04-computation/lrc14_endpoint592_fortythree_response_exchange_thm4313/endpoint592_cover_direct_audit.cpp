// Direct arithmetic replay of an explicit endpoint-592 response cover.

#define ENDPOINT592_RESPONSE_POOL_SCOUT_MAIN endpoint592_pool_hidden_main
#include "endpoint592_response_multiplicity.cpp"
#undef ENDPOINT592_RESPONSE_POOL_SCOUT_MAIN

namespace {
constexpr u64 kCover43Fnv592 = UINT64_C(0xca3cb80f471f2e7e);

std::vector<u32> read_cover592(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open endpoint592 cover");
    std::string line;
    require(std::getline(input, line) &&
                line == "mask_hex,rank,total_weight,q96,q100,q105,q192,q210,q256,q416",
            "endpoint592 cover header changed");
    std::vector<u32> masks;
    std::set<u32> distinct;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos, "malformed cover row");
        const u32 mask = parse_mask_agent(line.substr(0, comma));
        require((std::popcount(mask) == 8 || std::popcount(mask) == 9) &&
                    distinct.insert(mask).second,
                "cover mask rank/distinctness changed");
        masks.push_back(mask);
    }
    require(input.eof() && masks.size() == 43 &&
                masks_fnv_agent(masks) == kCover43Fnv592,
            "cover identity changed");
    return masks;
}
}

int main(int argc, char** argv) {
    try {
        require(argc == 3, "usage: audit FAILURES2468 COVER43_CSV");
        const FailureUniverse592 failures = read_failures592(argv[1]);
        const std::vector<u32> cover = read_cover592(argv[2]);
        std::array<Geometry, 7> geometry;
        for (std::size_t row = 0; row < 7; ++row)
            geometry[row] = build_geometry(kQ592[row], 592);
        Bits592 covered{};
        u64 incidences = 0;
        std::array<u64, 7> row_incidences{};
        for (u32 mask : cover) {
            const Candidate592 candidate = response592(mask, geometry, failures);
            require(candidate.weight != 0, "cover contains null response");
            incidences += candidate.weight;
            for (std::size_t row = 0; row < 7; ++row)
                row_incidences[row] += row_weight592(candidate, row, failures);
            for (std::size_t word = 0; word < kTotalWords592; ++word)
                covered[word] |= candidate.response[word];
        }
        u64 missing = 0;
        for (std::size_t word = 0; word < kTotalWords592; ++word)
            missing += std::popcount(failures.valid[word] & ~covered[word]);
        require(missing == 0, "explicit cover missed frozen obligations");
        const u64 rank8 = std::count_if(
            cover.begin(), cover.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        require(rank8 == 2, "cover rank profile changed");
        std::cout << "LRC14_ENDPOINT592_COVER43_DIRECT_AUDIT_V1\n"
                  << "COVER 43 FNV " << std::hex << masks_fnv_agent(cover)
                  << std::dec << " RANK8 " << rank8 << " RANK9 "
                  << cover.size() - rank8 << '\n'
                  << "OBLIGATIONS " << kTotalFailures592 << " FNV "
                  << std::hex << kFailureFnv592 << std::dec << " INCIDENCES "
                  << incidences << " ROW_INCIDENCES";
        for (u64 value : row_incidences) std::cout << ' ' << value;
        std::cout << " MISSING " << missing
                  << "\nSCOPE DIRECT_FIXED_FAILURE_COVER_REPLAY_"
                     "NO_MINIMUM_CLAIM_NO_FULL467_EXCHANGE_NO_PHYSICAL_ENTRY_"
                     "NO_LRC14\nVERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT592_COVER43_ERROR " << error.what() << '\n';
        return 1;
    }
}
