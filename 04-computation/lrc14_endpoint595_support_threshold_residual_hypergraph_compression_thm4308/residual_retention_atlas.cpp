// Generic deleted-mask response atlas for a raw failure ledger (up to 128
// pair-tagged obligations).  Used to reinsert a minimal subset after an
// aggressive support-threshold deletion.

#define ENDPOINT617_RAW_VERIFY_MAIN generic_retention_hidden_main
#include "04-computation/lrc14_size_preserving_response_staircase_thm4300/endpoint617_exchanged_carrier_raw_verify.cpp"
#undef ENDPOINT617_RAW_VERIFY_MAIN

#include <fstream>

namespace {
struct ObligationG { PairAgent pair; u32 body = 0; };
struct ResponseG {
    std::array<u64, 2> word{};
    bool operator<(const ResponseG& other) const { return word < other.word; }
};
struct TypeG { u64 count = 0; u32 least = UINT32_MAX; };

std::vector<u32> read_masks_generic(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open deletion ledger");
    std::vector<u32> masks;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_mask_agent(token);
        require((std::popcount(mask) == 8 || std::popcount(mask) == 9) &&
                    distinct.insert(mask).second,
                "deletion mask invalid/duplicate");
        masks.push_back(mask);
    }
    require(input.eof() && !masks.empty(), "empty deletion ledger");
    return masks;
}

std::vector<ObligationG> read_failures_generic(const std::filesystem::path& path,
                                                u64& ledger_state) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open failure ledger");
    std::string line;
    require(std::getline(input, line) && line == "q,r,body_hex",
            "failure header changed");
    std::vector<ObligationG> obligations;
    std::set<std::tuple<int, int, u32>> distinct;
    Fnv ledger;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        int q, r; std::string token;
        fields >> q >> r >> token;
        const u32 body = parse_mask_agent(token);
        require(fields && q > 0 && q < r && std::popcount(body) == 9 &&
                    distinct.insert({q, r, body}).second,
                "failure invalid/duplicate");
        obligations.push_back({PairAgent{q, r}, body});
        ledger.add(q); ledger.add(r); ledger.add(body);
    }
    require(input.eof() && !obligations.empty() && obligations.size() <= 128,
            "failure count escaped 1..128");
    ledger_state = ledger.state;
    return obligations;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 4, "usage: generic_retention DELETE FAILURES ATLAS_CSV");
        const auto deleted = read_masks_generic(argv[1]);
        Fnv deletion_ledger;
        for (u32 mask : deleted) deletion_ledger.add(mask);
        u64 failure_fnv = 0;
        const auto obligations = read_failures_generic(argv[2], failure_fnv);
        std::vector<Geometry> geometry;
        for (const auto& obligation : obligations)
            geometry.push_back(build_geometry(obligation.pair.q, obligation.pair.r));

        std::map<ResponseG, TypeG> types;
        Fnv responder_ledger;
        u64 responders = 0;
        for (u32 mask : deleted) {
            ResponseG response;
            for (std::size_t index = 0; index < obligations.size(); ++index)
                if ((mask & obligations[index].body) == 0 &&
                    margin(geometry[index], mask).ticks >= 0)
                    response.word[index / 64] |= UINT64_C(1) << (index % 64);
            if (response.word[0] == 0 && response.word[1] == 0) continue;
            ++responders;
            responder_ledger.add(mask);
            auto& type = types[response];
            ++type.count;
            type.least = std::min(type.least, mask);
        }
        std::ofstream output(argv[3]);
        require(static_cast<bool>(output), "cannot create retention atlas");
        output << "w0,w1,count,least_mask\n";
        Fnv type_ledger;
        for (const auto& [response, type] : types) {
            type_ledger.add(response.word[0]); type_ledger.add(response.word[1]);
            type_ledger.add(type.count); type_ledger.add(type.least);
            output << std::hex << std::setfill('0') << std::setw(16)
                   << response.word[0] << ',' << std::setw(16) << response.word[1]
                   << std::dec << ',' << type.count << ',' << hex8(type.least) << '\n';
        }
        require(output.good(), "retention atlas write failed");
        std::cout << "LRC14_GENERIC_RETENTION_ATLAS_V1\nDELETED "
                  << deleted.size() << " FNV " << std::hex
                  << deletion_ledger.state << std::dec << " OBLIGATIONS "
                  << obligations.size() << " FNV " << std::hex << failure_fnv
                  << std::dec << '\n'
                  << "RESPONDERS " << responders << " FNV " << std::hex
                  << responder_ledger.state << std::dec << " TYPES "
                  << types.size() << " TYPE_FNV " << std::hex
                  << type_ledger.state << std::dec << "\nVERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "GENERIC_RETENTION_ERROR " << error.what() << '\n';
        return 1;
    }
}
