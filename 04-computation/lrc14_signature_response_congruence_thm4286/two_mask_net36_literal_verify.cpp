#define TWO_MASK_LITERAL_VERIFY_LIBRARY_ONLY
#include "two_mask_literal_verify.cpp"
#undef TWO_MASK_LITERAL_VERIFY_LIBRARY_ONLY

#include <set>

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 7,
                "usage: net36-literal DECK SIGS PRIOR REP1 REP2 REP3");
        const std::vector<u32> deck = read_deck(argv[1]);
        const std::vector<SigRow> rows = read_signatures(argv[2]);
        std::ifstream residual_in(argv[3]);
        require(static_cast<bool>(residual_in), "cannot open prior residual");
        std::set<std::pair<int, int>> prior_residual;
        FnvLocal residual_ledger;
        std::pair<int, int> previous{-1, -1};
        std::string line;
        while (std::getline(residual_in, line)) {
            const auto f = split_csv(line);
            require(f.size() == 2, "bad prior residual row");
            const std::pair<int, int> pair{std::stoi(f[0]), std::stoi(f[1])};
            require(previous < pair, "prior residual not unique and ordered");
            prior_residual.insert(pair);
            residual_ledger.add(pair.first); residual_ledger.add(pair.second);
            previous = pair;
        }
        require(prior_residual.size() == 23373 &&
                    residual_ledger.state == UINT64_C(0xc6ab0ae49ee32273),
                "prior residual count changed");

        constexpr std::array<std::size_t, 3> deleted = {107, 318, 374};
        std::array<u64, 7> allowed{};
        for (std::size_t index : deleted)
            allowed[index / 64] |= UINT64_C(1) << (index % 64);
        std::vector<u32> replacements;
        std::set<u32> replacement_set;
        for (int i = 4; i < 7; ++i) {
            const u32 mask =
                static_cast<u32>(std::stoul(argv[i], nullptr, 16));
            require(std::popcount(mask) == 8 && mask < (u32{1} << 30) &&
                        std::find(deck.begin(), deck.end(), mask) == deck.end() &&
                        replacement_set.insert(mask).second,
                    "invalid replacement");
            replacements.push_back(mask);
        }
        require(replacements == std::vector<u32>{UINT32_C(0x32043014),
                                                  UINT32_C(0x20807016),
                                                  UINT32_C(0x128c8012)},
                "replacement identity changed");

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool cells changed");
        std::vector<i64> walls(cells.size() + 1);
        walls.front() = cells.front().left;
        for (std::size_t i = 0; i < cells.size(); ++i) {
            require(cells[i].left == walls[i], "pool cells not contiguous");
            walls[i + 1] = cells[i].right;
        }
        require(walls.front() == 0 && walls.back() == COMMON, "walls changed");
        std::vector<IndexedRepair> geometry;
        for (u32 mask : deck) geometry.push_back(index_repair(mask, cells));
        std::vector<IndexedRepair> replacement_geometry;
        for (u32 mask : replacements)
            replacement_geometry.push_back(index_repair(mask, cells));

        u64 candidates = 0, common = 0, net = 0, equalities = 0;
        FnvLocal candidate_ledger, common_ledger, net_ledger;
        FnvLocal literal_row_ledger;
        i128 minimum_margin = 0, minimum_grid = 1;
        std::pair<int, int> minimum_pair{};
        u32 minimum_mask = 0;
        bool minimum_set = false;
        std::cout << "THM4286_TWO_MASK_NET36_LITERAL_V1\n";
        for (const SigRow& row : rows) {
            bool subset = true;
            for (std::size_t word = 0; word < row.words.size(); ++word)
                subset &= (row.words[word] & ~allowed[word]) == 0;
            if (!subset) continue;
            ++candidates;
            candidate_ledger.add(row.q); candidate_ledger.add(row.r);
            const LiteralPair pair = build_literal_pair(row.q, row.r);
            std::vector<i128> prefix(walls.size());
            for (std::size_t w = 0; w < walls.size(); ++w)
                prefix[w] = prefix_at(pair,
                    static_cast<i128>(walls[w]) * pair.pool_scale);

            std::array<u64, 7> literal_signature{};
            for (std::size_t i = 0; i < geometry.size(); ++i) {
                const i128 margin = repair_margin(geometry[i], prefix,
                                                   pair.grid);
                equalities += margin == 0;
                if (margin < 0)
                    literal_signature[i / 64] |= UINT64_C(1) << (i % 64);
            }
            require(literal_signature == row.words,
                    "literal old signature disagrees with frozen CSV");

            bool all_active = true;
            std::vector<i128> margins;
            for (const IndexedRepair& replacement : replacement_geometry) {
                const i128 margin = repair_margin(replacement, prefix,
                                                   pair.grid);
                margins.push_back(margin);
                equalities += margin == 0;
                all_active &= margin > 0;
            }
            literal_row_ledger.add(row.q); literal_row_ledger.add(row.r);
            for (u64 word : literal_signature) literal_row_ledger.add(word);
            for (std::size_t i = 0; i < replacements.size(); ++i) {
                literal_row_ledger.add(replacements[i]);
                add_i128(literal_row_ledger, margins[i]);
            }
            add_i128(literal_row_ledger, pair.grid);
            if (!all_active) continue;
            ++common;
            common_ledger.add(row.q); common_ledger.add(row.r);
            if (!prior_residual.contains({row.q, row.r})) continue;
            ++net;
            net_ledger.add(row.q); net_ledger.add(row.r);
            std::cout << "NET " << row.q << ',' << row.r << " MARGINS";
            for (i128 margin : margins)
                std::cout << ' ' << decimal(margin);
            std::cout << " DEN " << decimal(pair.grid) << '\n';
            for (std::size_t i = 0; i < replacements.size(); ++i) {
                if (!minimum_set || fraction_less_exact(
                                        margins[i], pair.grid,
                                        minimum_margin, minimum_grid)) {
                    minimum_set = true;
                    minimum_margin = margins[i];
                    minimum_grid = pair.grid;
                    minimum_pair = {row.q, row.r};
                    minimum_mask = replacements[i];
                }
            }
        }
        require(candidates == 112 && common == 106 && net == 36 &&
                    equalities == 0 &&
                    candidate_ledger.state == UINT64_C(0xd2fcb3bf714762ef) &&
                    common_ledger.state == UINT64_C(0x2dab535d17359dfb) &&
                    net_ledger.state == UINT64_C(0x8da60395e47e11a3) &&
                    literal_row_ledger.state == UINT64_C(0xdccb7a781070d099) &&
                    minimum_pair == std::pair<int, int>{332, 496} &&
                    minimum_mask == UINT32_C(0x128c8012) &&
                    minimum_margin == static_cast<i128>(344114666009478) &&
                    minimum_grid == static_cast<i128>(46934503178603040),
                "literal census/ledger changed");

        std::vector<u32> rebuilt;
        for (std::size_t i = 0; i < deck.size(); ++i)
            if (std::find(deleted.begin(), deleted.end(), i) == deleted.end())
                rebuilt.push_back(deck[i]);
        rebuilt.insert(rebuilt.end(), replacements.begin(), replacements.end());
        require(rebuilt.size() == 421 &&
                    std::set<u32>(rebuilt.begin(), rebuilt.end()).size() == 421 &&
                    mask_fnv(rebuilt) == UINT64_C(0xc9ac86709cda10df),
                "rebuilt deck identity changed");
        u64 checks = 0, max_checks = 0, failures = 0, body_row_fnv = 0;
        const u64 cover_fnv = scan_body_cover(rebuilt, checks, max_checks,
                                               failures, body_row_fnv);
        require(failures == 0 && checks == UINT64_C(405320191) &&
                    max_checks == 421 &&
                    cover_fnv == UINT64_C(0xe9584c43ebc26a0b) &&
                    body_row_fnv == UINT64_C(0x2984a661ec30ab5e),
                "rebuilt global body scan changed");

        std::cout << "SUMMARY CANDIDATES " << candidates
                  << " CANDIDATE_FNV " << std::hex << candidate_ledger.state
                  << std::dec << " COMMON " << common << " COMMON_FNV "
                  << std::hex << common_ledger.state << std::dec
                  << " NET " << net << " NET_FNV " << std::hex
                  << net_ledger.state << " LITERAL_ROW_FNV "
                  << literal_row_ledger.state << std::dec
                  << " EQUALITIES " << equalities << '\n'
                  << "MIN_NET_GAP " << minimum_pair.first << ','
                  << minimum_pair.second << " MASK " << std::hex
                  << minimum_mask << std::dec << " NUM "
                  << decimal(minimum_margin) << " DEN "
                  << decimal(minimum_grid) << '\n'
                  << "REBUILT421 FNV " << std::hex << mask_fnv(rebuilt)
                  << std::dec << " BODY_CHECKS " << checks
                  << " MAX_CHECKS " << max_checks << " FAILURES " << failures
                  << " COVER_FNV " << std::hex << cover_fnv
                  << " BODY_ROW_FNV " << body_row_fnv << std::dec << '\n'
                  << "VERDICT PASS DETACHED_LITERAL_NET36_AND_GLOBAL_BODY\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "NET36_LITERAL_ERROR " << error.what() << '\n';
        return 1;
    }
}
