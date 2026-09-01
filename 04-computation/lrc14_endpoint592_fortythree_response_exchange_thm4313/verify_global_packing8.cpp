// Certificate-only verifier for an eight-obligation pair-tagged response
// packing in the frozen endpoint-592 failure universe.  It scans each rank-8
// and rank-9 mask exactly once and contains no search heuristic.

#define main endpoint592_packing8_hidden_main
#include "04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/r632_detached_hostile_survivor.cpp"
#undef main

namespace {

constexpr u64 kFailureFnv = UINT64_C(0x2209b8d6760280cc);
constexpr std::array<u32, 7> kQ256Bodies = {
    0x0530e401, 0x06167400, 0x07246088, 0x0f141208,
    0x13116408, 0x1430240b, 0x16624401};
constexpr u32 kQ105Body = UINT32_C(0x06d81088);
constexpr std::array<std::array<u64, 2>, 7> kQ256Multiplicity = {{
    {229, 4269}, {93, 3029}, {155, 5296}, {700, 12247},
    {294, 8328}, {229, 9018}, {864, 10920}}};
constexpr std::array<u64, 2> kQ105Multiplicity = {819, 13251};

struct Tagged {
    int q = 0;
    int r = 0;
    u32 body = 0;
    auto operator<=>(const Tagged&) const = default;
};

std::set<Tagged> read_universe(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open failure ledger");
    std::string line;
    require(std::getline(input, line) && line == "q,r,body_hex",
            "failure header changed");
    std::set<Tagged> universe;
    Fnv ledger;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        int q = 0, r = 0;
        std::string token;
        fields >> q >> r >> token;
        require(static_cast<bool>(fields) && r == 592,
                "malformed failure row");
        const u32 body = static_cast<u32>(std::stoul(token, nullptr, 16));
        require(std::popcount(body) == 9 &&
                    universe.insert({q, r, body}).second,
                "body rank/distinctness changed");
        ledger.add(q); ledger.add(r); ledger.add(body);
    }
    require(input.eof() && universe.size() == 2468 &&
                ledger.state == kFailureFnv,
            "failure universe identity changed");
    return universe;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "usage: verify-global-packing8 FAILURES2468");
        const std::set<Tagged> universe = read_universe(argv[1]);
        std::set<Tagged> certificate;
        certificate.insert({105, 592, kQ105Body});
        for (u32 body : kQ256Bodies) certificate.insert({256, 592, body});
        require(certificate.size() == 8, "certificate distinctness changed");
        for (const Tagged& item : certificate)
            require(universe.contains(item), "certificate body left universe");

        const Geometry geometry105 = build_geometry(105, 592);
        const Geometry geometry256 = build_geometry(256, 592);
        std::array<std::array<u64, 2>, 7> q256_responses{};
        std::array<u64, 2> q105_responses{};
        std::array<u64, 2> active256{};
        std::array<u64, 2> q256_multi_response{};
        std::array<u64, 2> cross_co_response{};
        std::array<std::array<u64, 7>, 7> q256_pair_counts{};
        std::array<u64, 7> cross_pair_counts{};
        std::array<u64, 2> all_masks{};
        const u32 limit = u32{1} << 30;
        for (unsigned rank : {8U, 9U}) {
            const std::size_t slot = rank - 8;
            for (u32 mask = (u32{1} << rank) - 1; mask < limit;
                 mask = next_combination(mask)) {
                ++all_masks[slot];
                const bool is_active256 = margin(geometry256, mask).ticks >= 0;
                active256[slot] += is_active256;
                unsigned response256 = 0;
                if (is_active256)
                    for (unsigned index = 0; index < kQ256Bodies.size(); ++index)
                        if ((mask & kQ256Bodies[index]) == 0) {
                            response256 |= 1U << index;
                            ++q256_responses[index][slot];
                        }
                const bool response105 = (mask & kQ105Body) == 0 &&
                    margin(geometry105, mask).ticks >= 0;
                q105_responses[slot] += response105;
                if (std::popcount(response256) >= 2) {
                    ++q256_multi_response[slot];
                    for (unsigned i = 0; i < 7; ++i)
                        for (unsigned j = i + 1; j < 7; ++j)
                            if ((response256 >> i & 1U) &&
                                (response256 >> j & 1U))
                                ++q256_pair_counts[i][j];
                }
                if (response105 && response256 != 0) {
                    ++cross_co_response[slot];
                    for (unsigned i = 0; i < 7; ++i)
                        if (response256 >> i & 1U)
                            ++cross_pair_counts[i];
                }
            }
        }
        require(all_masks == std::array<u64, 2>{5852925, 14307150},
                "rank universes changed");
        require(q256_responses == kQ256Multiplicity &&
                    q105_responses == kQ105Multiplicity,
                "certificate responder multiplicity changed");
        require(q256_multi_response == std::array<u64, 2>{0, 0} &&
                    cross_co_response == std::array<u64, 2>{0, 0},
                "certificate has a co-response");
        for (const auto& row : q256_pair_counts)
            for (u64 count : row) require(count == 0, "q256 pair count nonzero");
        for (u64 count : cross_pair_counts)
            require(count == 0, "cross pair count nonzero");

        Fnv certificate_ledger;
        for (const Tagged& item : certificate) {
            certificate_ledger.add(item.q);
            certificate_ledger.add(item.r);
            certificate_ledger.add(item.body);
        }
        std::cout << "LRC14_ENDPOINT592_GLOBAL_PACKING8_VERIFY_V1\n"
                  << "UNIVERSE 2468 FNV " << std::hex << kFailureFnv
                  << std::dec << " CERTIFICATE 8 FNV " << std::hex
                  << certificate_ledger.state << std::dec << '\n'
                  << "GEOMETRY105 CELLS " << geometry105.cells << " CLASSES "
                  << geometry105.classes.size() << " GEOMETRY256 CELLS "
                  << geometry256.cells << " CLASSES "
                  << geometry256.classes.size() << '\n';
        for (unsigned rank : {8U, 9U}) {
            const std::size_t slot = rank - 8;
            std::cout << "RANK " << rank << " ALL " << all_masks[slot]
                      << " ACTIVE256 " << active256[slot]
                      << " Q105_RESPONSES "
                      << q105_responses[slot] << " Q256_RESPONSES";
            for (unsigned index = 0; index < 7; ++index)
                std::cout << ' ' << q256_responses[index][slot];
            std::cout << " Q256_MULTI " << q256_multi_response[slot]
                      << " CROSS_MULTI " << cross_co_response[slot] << '\n';
        }
        std::cout << "OBLIGATIONS 105:" << hex8(kQ105Body);
        for (u32 body : kQ256Bodies) std::cout << " 256:" << hex8(body);
        std::cout << "\nPAIR_AUDIT Q256_PAIRS 21 CROSS_PAIRS 7 TOTAL 28 "
                     "CO_RESPONSES 0\n"
                  << "LOWER_BOUND 8 MECHANISM pair-tagged-response-packing\n"
                  << "SCOPE FIXED_FAILURE_UNIVERSE_COMPLETE_RANK8_RANK9_"
                     "RESPONSE_COVER_LOWER_BOUND_ONLY_NO_MINIMUM_NO_LRC14\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "GLOBAL_PACKING8_ERROR " << error.what() << '\n';
        return 1;
    }
}
