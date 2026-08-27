#define RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
#include "04-computation/lrc14_joint421_global_common_carrier_thm4281/response_pattern_atlas.cpp"
#undef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY

#include <fstream>
#include <iomanip>
#include <set>

namespace {
std::vector<u32> read_rebuilt(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open rebuilt deck");
    std::vector<u32> deck;
    std::set<u32> seen;
    std::string token;
    while (input >> token) {
        std::size_t used = 0;
        const u64 wide = std::stoull(token, &used, 16);
        require(used == token.size() && wide < (UINT64_C(1) << 30),
                "invalid rebuilt mask token");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && seen.insert(mask).second,
                "bad/duplicate rebuilt mask");
        deck.push_back(mask);
    }
    require(!deck.empty() && deck.size() <= 1000, "bad rebuilt deck size");
    return deck;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 4, "usage: rebuilt-pair-activity DECK Q R");
        init_choose8_local();
        const std::vector<u32> deck = read_rebuilt(argv[1]);
        const int q = std::stoi(argv[2]);
        const int r = std::stoi(argv[3]);
        require(q > 0 && q < r, "invalid target pair");
        const std::vector<Cell> cells = build_pool_cells();
        const ActiveUniverse active = build_active_universe(cells, q, r);
        FnvLocal deck_ledger;
        u64 deck_active = 0;
        std::vector<std::pair<std::size_t, u32>> inactive;
        for (std::size_t index = 0; index < deck.size(); ++index) {
            deck_ledger.add(deck[index]);
            const u64 rank = colex_rank8_local(deck[index]);
            if (active.active[rank]) ++deck_active;
            else inactive.push_back({index, deck[index]});
        }
        std::cout << "THM4281_REBUILT_PAIR_ACTIVITY_V1\nPAIR " << q << ',' << r
                  << " RANK8_ACTIVE " << active.count << " ACTIVE_FNV "
                  << std::hex << active.fnv << std::dec
                  << " ZETA_OPERATIONS " << active.zeta_operations << '\n'
                  << "DECK " << deck.size() << " FNV " << std::hex
                  << deck_ledger.state << std::dec << " ACTIVE " << deck_active
                  << " INACTIVE " << inactive.size() << '\n';
        for (const auto& [index, mask] : inactive)
            std::cout << "INACTIVE_INDEX " << index << " MASK " << std::hex
                      << std::setw(8) << std::setfill('0') << mask << std::dec
                      << std::setfill(' ') << '\n';
        const std::size_t replacement_begin = deck.size() >= 6
            ? deck.size() - 6 : deck.size();
        for (std::size_t index = replacement_begin; index < deck.size(); ++index)
            std::cout << "REPLACEMENT_INDEX " << index << " MASK " << std::hex
                      << std::setw(8) << std::setfill('0') << deck[index]
                      << std::dec << std::setfill(' ') << " ACTIVE "
                      << static_cast<unsigned>(
                             active.active[colex_rank8_local(deck[index])])
                      << '\n';
        std::cout << "VERDICT " << (inactive.empty() ? "PASS ALL_ACTIVE" :
                                                        "FAIL INACTIVE_PRESENT")
                  << '\n';
        return inactive.empty() ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "REBUILT_PAIR_ACTIVITY_ERROR " << error.what() << '\n';
        return 1;
    }
}
