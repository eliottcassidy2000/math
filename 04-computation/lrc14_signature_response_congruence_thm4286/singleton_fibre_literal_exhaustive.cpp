#define SINGLETON_LITERAL_LIBRARY_ONLY
#include "04-computation/lrc14_signature_response_congruence_thm4286/singleton_fibre_literal_verify.cpp"
#undef SINGLETON_LITERAL_LIBRARY_ONLY

#include <map>

namespace {

struct FibreLiteral {
    int q = 0;
    int r = 0;
    i128 grid = 0;
    std::vector<i128> prefix;
};

std::vector<u32> private_bodies(const std::vector<u32>& retained,
                                u64& checks, u64& body_count) {
    std::vector<u32> failures;
    checks = body_count = 0;
    u32 body = (u32{1} << 9) - 1;
    while (body < (u32{1} << 30)) {
        ++body_count;
        bool covered = false;
        for (u32 mask : retained) {
            ++checks;
            if ((body & mask) == 0) {
                covered = true;
                break;
            }
        }
        if (!covered) failures.push_back(body);
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    require(body_count == UINT64_C(14307150), "body enumeration changed");
    return failures;
}

u64 mask_list_fnv(const std::vector<u32>& masks) {
    FnvLocal ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 6,
                "usage: literal-exhaustive DECK SIGNATURES INDEX Q R");
        const std::vector<u32> deck = read_deck(argv[1]);
        const std::vector<SigRow> rows = read_signatures(argv[2]);
        const std::size_t index = std::stoul(argv[3]);
        const int target_q = std::stoi(argv[4]);
        const int target_r = std::stoi(argv[5]);
        require(index < deck.size(), "bad index");

        const auto wanted = singleton_signature(index);
        std::vector<std::pair<int, int>> fibre;
        std::size_t target_position = rows.size();
        for (const SigRow& row : rows) {
            if (row.words != wanted) continue;
            if (row.q == target_q && row.r == target_r)
                target_position = fibre.size();
            fibre.push_back({row.q, row.r});
        }
        require(target_position < fibre.size() && fibre.size() <= 64,
                "target missing from singleton fibre");

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool cells changed");
        std::vector<i64> walls(cells.size() + 1);
        walls.front() = cells.front().left;
        for (std::size_t i = 0; i < cells.size(); ++i) {
            require(cells[i].left == walls[i], "pool cells not contiguous");
            walls[i + 1] = cells[i].right;
        }
        require(walls.front() == 0 && walls.back() == COMMON, "walls changed");

        std::vector<FibreLiteral> literal;
        literal.reserve(fibre.size());
        FnvLocal fibre_ledger;
        for (const auto [q, r] : fibre) {
            fibre_ledger.add(q);
            fibre_ledger.add(r);
            const LiteralPair pair = build_literal_pair(q, r);
            FibreLiteral item;
            item.q = q;
            item.r = r;
            item.grid = pair.grid;
            item.prefix.resize(walls.size());
            for (std::size_t w = 0; w < walls.size(); ++w)
                item.prefix[w] = prefix_at(
                    pair, static_cast<i128>(walls[w]) * pair.pool_scale);
            literal.push_back(std::move(item));
        }

        std::vector<u32> retained;
        retained.reserve(deck.size() - 1);
        for (std::size_t i = 0; i < deck.size(); ++i)
            if (i != index) retained.push_back(deck[i]);
        require(retained.size() == 420, "retained deck size changed");

        u64 body_checks = 0;
        u64 body_count = 0;
        const std::vector<u32> obligations =
            private_bodies(retained, body_checks, body_count);
        FnvLocal obligation_ledger;
        u32 obligation_union = 0;
        for (u32 body : obligations) {
            obligation_ledger.add(body);
            obligation_union |= body;
        }

        std::vector<u32> responders;
        std::vector<u32> fibre_wide;
        std::map<u64, u64> activity_classes;
        std::vector<std::array<u64, 8>> row_profiles(fibre.size());
        FnvLocal margin_ledger;
        u64 active_incidences = 0;
        u64 equalities = 0;
        u64 mask_count = 0;
        u32 candidate = (u32{1} << 8) - 1;
        while (candidate < (u32{1} << 30)) {
            ++mask_count;
            bool covers = true;
            for (u32 body : obligations) {
                if ((candidate & body) != 0) {
                    covers = false;
                    break;
                }
            }
            if (covers) {
                const IndexedRepair repair = index_repair(candidate, cells);
                const std::size_t responder_index = responders.size();
                responders.push_back(candidate);
                u64 pattern = 0;
                for (std::size_t i = 0; i < literal.size(); ++i) {
                    const i128 m = repair_margin(
                        repair, literal[i].prefix, literal[i].grid);
                    equalities += m == 0;
                    if (m > 0) {
                        pattern |= UINT64_C(1) << i;
                        row_profiles[i][responder_index / 64] |=
                            UINT64_C(1) << (responder_index % 64);
                        ++active_incidences;
                    }
                    margin_ledger.add(candidate);
                    margin_ledger.add(literal[i].q);
                    margin_ledger.add(literal[i].r);
                    add_i128(margin_ledger, m);
                    add_i128(margin_ledger, literal[i].grid);
                }
                ++activity_classes[pattern];
                const u64 full = (UINT64_C(1) << literal.size()) - 1;
                if (pattern == full) fibre_wide.push_back(candidate);
            }
            const u32 next = next_combination(candidate);
            if (next <= candidate) break;
            candidate = next;
        }

        require(mask_count == UINT64_C(5852925), "mask enumeration changed");
        require(fibre.size() == 36 &&
                    fibre_ledger.state == UINT64_C(0x3d92ab45b46a72c0),
                "singleton fibre ledger changed");
        require(obligations.size() == 8 &&
                    obligation_ledger.state == UINT64_C(0xe3fd616a4aa21839) &&
                    obligation_union == UINT32_C(0x33dfdc06),
                "private body ledger changed");
        require(responders.size() == 495 &&
                    mask_list_fnv(responders) == UINT64_C(0x3c76be2ab12086ed),
                "one-mask responder census changed");
        require(activity_classes.size() == 83 && active_incidences == 1193 &&
                    equalities == 0 &&
                    margin_ledger.state == UINT64_C(0x35e5d65bc81ae0e4),
                "literal activity census changed");
        require(fibre_wide == std::vector<u32>{UINT32_C(0x042022c9),
                                                UINT32_C(0x0c202289)},
                "fibre-wide responder census changed");

        std::map<std::array<u64, 8>, std::vector<std::pair<int, int>>>
            profile_classes;
        FnvLocal profile_ledger;
        for (std::size_t i = 0; i < fibre.size(); ++i) {
            profile_classes[row_profiles[i]].push_back(fibre[i]);
            profile_ledger.add(fibre[i].first);
            profile_ledger.add(fibre[i].second);
            for (u64 word : row_profiles[i]) profile_ledger.add(word);
        }
        std::vector<std::vector<std::pair<int, int>>> collisions;
        for (const auto& [profile, pairs] : profile_classes)
            if (pairs.size() > 1) collisions.push_back(pairs);
        require(profile_classes.size() == 35 &&
                    profile_ledger.state == UINT64_C(0x553dca3fbf7bea86) &&
                    collisions.size() == 1 &&
                    collisions.front() ==
                        std::vector<std::pair<int, int>>{{301, 366},
                                                         {366, 547}},
                "row-response refinement changed");

        std::cout << "THM4286_SINGLETON_FIBRE_LITERAL_EXHAUSTIVE_V1\n"
                  << "TARGET " << target_q << ',' << target_r << " INDEX "
                  << index << " OLD_MASK " << std::hex << deck[index]
                  << std::dec << '\n'
                  << "FIBRE " << fibre.size() << " FNV " << std::hex
                  << fibre_ledger.state << std::dec << '\n'
                  << "BODY_SCAN BODIES " << body_count << " CHECKS "
                  << body_checks << " OBLIGATIONS " << obligations.size()
                  << " FNV " << std::hex << obligation_ledger.state
                  << " UNION " << obligation_union << std::dec << '\n'
                  << "MASK_SCAN " << mask_count << " ONE_MASK_RESPONDERS "
                  << responders.size() << " FNV " << std::hex
                  << mask_list_fnv(responders) << std::dec << '\n'
                  << "LITERAL_ACTIVITY_TESTS "
                  << responders.size() * fibre.size() << " ACTIVE_INCIDENCES "
                  << active_incidences << " EQUALITIES " << equalities
                  << " CLASSES " << activity_classes.size() << " ROW_FNV "
                  << std::hex << margin_ledger.state << std::dec << '\n'
                  << "FIBRE_WIDE " << fibre_wide.size() << " FNV "
                  << std::hex << mask_list_fnv(fibre_wide) << " MASKS";
        for (u32 mask : fibre_wide)
            std::cout << ' ' << std::setw(8) << std::setfill('0') << mask;
        std::cout << std::dec << "\nROW_RESPONSE_PROFILES "
                  << profile_classes.size() << " FNV " << std::hex
                  << profile_ledger.state << std::dec << " COLLISION "
                  << collisions.front()[0].first << ','
                  << collisions.front()[0].second << ';'
                  << collisions.front()[1].first << ','
                  << collisions.front()[1].second
                  << "\nVERDICT PASS DETACHED_LITERAL_EXHAUSTIVE_CENSUS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "LITERAL_EXHAUSTIVE_ERROR " << error.what() << '\n';
        return 1;
    }
}
