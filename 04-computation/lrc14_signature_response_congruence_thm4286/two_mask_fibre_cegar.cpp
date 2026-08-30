#define RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
#include "04-computation/lrc14_joint421_global_common_carrier_thm4281/response_pattern_atlas.cpp"
#undef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY

#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>

namespace {

struct SigRow {
    int q = 0;
    int r = 0;
    std::array<u64, 7> words{};
};

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> out;
    std::istringstream in(line);
    std::string field;
    while (std::getline(in, field, ',')) out.push_back(field);
    return out;
}

std::vector<SigRow> read_signatures(const std::string& path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open signatures");
    std::string line;
    require(static_cast<bool>(std::getline(in, line)) &&
                line == "q,r,inactive_count,w0,w1,w2,w3,w4,w5,w6",
            "signature header changed");
    std::vector<SigRow> rows;
    while (std::getline(in, line)) {
        const auto f = split_csv(line);
        require(f.size() == 10, "bad signature row");
        SigRow row;
        row.q = std::stoi(f[0]);
        row.r = std::stoi(f[1]);
        unsigned count = 0;
        for (int i = 0; i < 7; ++i) {
            row.words[i] = std::stoull(f[3 + i], nullptr, 16);
            count += std::popcount(row.words[i]);
        }
        require(count == static_cast<unsigned>(std::stoul(f[2])),
                "signature count changed");
        rows.push_back(row);
    }
    require(rows.size() == 24223, "signature universe changed");
    return rows;
}

u64 mask_fnv(const std::vector<u32>& masks) {
    FnvLocal ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

std::vector<u32> read_deck(const std::string& path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open deck");
    std::vector<u32> deck;
    std::set<u32> seen;
    std::string token;
    while (in >> token) {
        const u32 mask = static_cast<u32>(std::stoul(token, nullptr, 16));
        require(std::popcount(mask) == 8 && mask < (u32{1} << 30) &&
                    seen.insert(mask).second,
                "bad deck mask");
        deck.push_back(mask);
    }
    require(deck.size() == 421 &&
                mask_fnv(deck) == UINT64_C(0x20d63dd42fe8150e),
            "deck identity changed");
    return deck;
}

std::array<u64, 7> signature_for(const std::vector<std::size_t>& indices) {
    std::array<u64, 7> out{};
    for (std::size_t index : indices) {
        require(index < 421, "bad deleted index");
        out[index / 64] |= UINT64_C(1) << (index % 64);
    }
    return out;
}

std::vector<u32> exposed_bodies(const std::vector<u32>& deck,
                                const std::vector<std::size_t>& deleted,
                                u64& candidate_bodies, u64& retained_checks) {
    std::array<unsigned char, 421> is_deleted{};
    for (std::size_t index : deleted) is_deleted[index] = 1;
    std::vector<u32> out;
    candidate_bodies = retained_checks = 0;
    u64 bodies = 0;
    u32 body = (u32{1} << 9) - 1;
    while (body < (u32{1} << 30)) {
        ++bodies;
        bool deleted_witness = false;
        for (std::size_t index : deleted)
            deleted_witness |= (body & deck[index]) == 0;
        if (deleted_witness) {
            ++candidate_bodies;
            bool retained_witness = false;
            for (std::size_t i = 0; i < deck.size(); ++i) {
                if (is_deleted[i]) continue;
                ++retained_checks;
                if ((body & deck[i]) == 0) {
                    retained_witness = true;
                    break;
                }
            }
            if (!retained_witness) out.push_back(body);
        }
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    require(bodies == UINT64_C(14307150), "body universe changed");
    return out;
}

struct CoverResult {
    int minimum = -1;
    std::vector<u32> masks;
};

CoverResult exact_cover(const std::map<u64, std::pair<u64, u32>>& classes,
                        u64 full) {
    std::vector<std::pair<u64, u32>> maximal;
    for (const auto& [pattern, data] : classes) {
        if (pattern == 0) continue;
        bool dominated = false;
        for (const auto& [other, other_data] : classes) {
            if (pattern != other && (pattern | other) == other) {
                dominated = true;
                break;
            }
        }
        if (!dominated) maximal.push_back({pattern, data.second});
    }
    CoverResult out;
    std::unordered_map<u64, std::pair<int, std::pair<u64, u32>>> state;
    state[0] = {0, {0, 0}};
    std::vector<u64> frontier{0};
    for (int depth = 0; depth <= 20 && !frontier.empty(); ++depth) {
        if (state.contains(full)) break;
        std::vector<u64> next_frontier;
        for (u64 covered : frontier) {
            for (const auto& [pattern, mask] : maximal) {
                const u64 next = covered | pattern;
                if (next == covered || state.contains(next)) continue;
                state[next] = {depth + 1, {covered, mask}};
                next_frontier.push_back(next);
            }
        }
        frontier = std::move(next_frontier);
    }
    if (!state.contains(full)) return out;
    out.minimum = state[full].first;
    for (u64 s = full; s != 0;) {
        const auto prev = state[s].second;
        out.masks.push_back(prev.second);
        s = prev.first;
    }
    std::reverse(out.masks.begin(), out.masks.end());
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 8,
                "usage: two-mask-cegar DECK SIGS INDEX1 INDEX2 Q R MODE");
        init_choose8_local();
        const std::vector<u32> deck = read_deck(argv[1]);
        const std::vector<SigRow> rows = read_signatures(argv[2]);
        std::vector<std::size_t> deleted{std::stoul(argv[3]),
                                         std::stoul(argv[4])};
        std::sort(deleted.begin(), deleted.end());
        require(deleted[0] != deleted[1], "duplicate deleted index");
        const int target_q = std::stoi(argv[5]);
        const int target_r = std::stoi(argv[6]);
        const std::string mode = argv[7];
        require(mode == "target" || mode == "fibre", "bad mode");
        require(deleted == std::vector<std::size_t>{107, 374} &&
                    target_q == 512 && target_r == 644 && mode == "fibre",
                "target surgery identity changed");
        const auto wanted = signature_for(deleted);
        std::vector<std::pair<int, int>> fibre;
        std::size_t target_position = 0;
        bool target_found = false;
        for (const SigRow& row : rows) {
            if (row.words != wanted) continue;
            if (row.q == target_q && row.r == target_r) {
                target_found = true;
                target_position = fibre.size();
            }
            fibre.push_back({row.q, row.r});
        }
        require(target_found && !fibre.empty(), "target missing from fibre");

        u64 candidate_bodies = 0;
        u64 retained_checks = 0;
        const std::vector<u32> obligations = exposed_bodies(
            deck, deleted, candidate_bodies, retained_checks);
        require(!obligations.empty() && obligations.size() <= 63,
                "unsupported obligation count");
        FnvLocal obligation_ledger;
        u32 obligation_union = 0;
        for (u32 body : obligations) {
            obligation_ledger.add(body);
            obligation_union |= body;
        }
        const u64 full = (UINT64_C(1) << obligations.size()) - 1;

        const std::vector<Cell> cells = build_pool_cells();
        std::vector<ActiveUniverse> active;
        for (const auto [q, r] : fibre) {
            active.push_back(build_active_universe(cells, q, r));
            std::cout << "ACTIVE " << q << ',' << r << " COUNT "
                      << active.back().count << " FNV " << std::hex
                      << active.back().fnv << std::dec << '\n';
        }

        std::vector<u64> responses(EXPECTED_REPAIRS, 0);
        u64 response_incidences = 0;
        for (std::size_t o = 0; o < obligations.size(); ++o) {
            enumerate_disjoint_repairs(obligations[o], [&](u32, u64 rank) {
                responses[rank] |= UINT64_C(1) << o;
                ++response_incidences;
            });
        }

        std::map<u64, std::pair<u64, u32>> target_classes;
        std::map<u64, std::pair<u64, u32>> fibre_classes;
        std::map<std::pair<u64, u64>, std::pair<u64, u32>> joint_classes;
        u64 target_union = 0;
        u64 fibre_union = 0;
        u64 nonempty = 0;
        u64 fibrewide = 0;
        for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
            const u64 response = responses[rank];
            if (response == 0) continue;
            ++nonempty;
            u64 activity_word = 0;
            for (std::size_t i = 0; i < fibre.size(); ++i)
                if (active[i].active[rank]) activity_word |= UINT64_C(1) << i;
            const u32 mask = unrank_colex8(rank);
            joint_classes[{response, activity_word}].first++;
            auto& joint_least = joint_classes[{response, activity_word}].second;
            if (joint_least == 0 || mask < joint_least) joint_least = mask;
            if (active[target_position].active[rank]) {
                target_union |= response;
                auto& data = target_classes[response];
                ++data.first;
                if (data.second == 0 || mask < data.second) data.second = mask;
            }
            const u64 all_rows = (UINT64_C(1) << fibre.size()) - 1;
            if (activity_word == all_rows) {
                ++fibrewide;
                fibre_union |= response;
                auto& data = fibre_classes[response];
                ++data.first;
                if (data.second == 0 || mask < data.second) data.second = mask;
            }
        }

        const CoverResult target_cover = exact_cover(target_classes, full);
        const CoverResult fibre_cover = exact_cover(fibre_classes, full);
        FnvLocal fibre_ledger;
        for (const auto [q, r] : fibre) {
            fibre_ledger.add(q); fibre_ledger.add(r);
        }
        FnvLocal class_ledger;
        for (const auto& [key, data] : joint_classes) {
            class_ledger.add(key.first); class_ledger.add(key.second);
            class_ledger.add(data.first); class_ledger.add(data.second);
        }
        const std::array<u64, 5> expected_active_counts = {
            2837507, 2534838, 3011668, 2509053, 2271751};
        const std::array<u64, 5> expected_active_fnvs = {
            UINT64_C(0xd2a05625a4d1d334),
            UINT64_C(0xc0fb529b905944cd),
            UINT64_C(0x0162b1e23d99cdf2),
            UINT64_C(0xbce628bad1959acb),
            UINT64_C(0x9f37456f9594fb36)};
        for (std::size_t i = 0; i < active.size(); ++i) {
            require(active[i].count == expected_active_counts[i] &&
                        active[i].fnv == expected_active_fnvs[i],
                    "active-universe identity changed");
        }
        require(fibre.size() == 5 &&
                    fibre_ledger.state == UINT64_C(0xd1822d3aaeb858d1) &&
                    deck[107] == UINT32_C(0x128c8900) &&
                    deck[374] == UINT32_C(0x20027016) &&
                    candidate_bodies == UINT64_C(992838) &&
                    retained_checks == UINT64_C(26377164) &&
                    obligations.size() == 20 &&
                    obligation_ledger.state == UINT64_C(0x534bd49ebf6aa688) &&
                    obligation_union == UINT32_C(0x3ffff7ef) &&
                    response_incidences == UINT64_C(4069800) &&
                    nonempty == UINT64_C(2139042) &&
                    joint_classes.size() == 12087 &&
                    class_ledger.state == UINT64_C(0x5b4176259e5d38a7) &&
                    target_classes.size() == 956 &&
                    std::popcount(target_union) == 20 &&
                    target_cover.minimum == 3 &&
                    target_cover.masks ==
                        std::vector<u32>{UINT32_C(0x32043014),
                                         UINT32_C(0x00827016),
                                         UINT32_C(0x12848122)} &&
                    fibrewide == UINT64_C(483491) &&
                    fibre_classes.size() == 915 &&
                    std::popcount(fibre_union) == 20 &&
                    fibre_cover.minimum == 3 &&
                    fibre_cover.masks ==
                        std::vector<u32>{UINT32_C(0x32043014),
                                         UINT32_C(0x20807016),
                                         UINT32_C(0x128c8012)},
                "two-mask response audit changed");
        std::cout << "THM4286_TWO_MASK_FIBRE_CEGAR_V1\n"
                  << "TARGET " << target_q << ',' << target_r << " DELETED "
                  << deleted[0] << ':' << std::hex << deck[deleted[0]] << ' '
                  << std::dec << deleted[1] << ':' << std::hex
                  << deck[deleted[1]] << std::dec << '\n'
                  << "FIBRE " << fibre.size() << " FNV " << std::hex
                  << fibre_ledger.state << std::dec << " ROWS";
        for (const auto [q, r] : fibre) std::cout << ' ' << q << ',' << r;
        std::cout << '\n'
                  << "BODY_CANDIDATES " << candidate_bodies
                  << " RETAINED_CHECKS " << retained_checks
                  << " OBLIGATIONS " << obligations.size() << " FNV "
                  << std::hex << obligation_ledger.state << " UNION "
                  << obligation_union << std::dec << " UNION_SIZE "
                  << std::popcount(obligation_union) << '\n'
                  << "RESPONSE_INCIDENCES " << response_incidences
                  << " NONEMPTY_MASKS " << nonempty
                  << " JOINT_CLASSES " << joint_classes.size() << " FNV "
                  << std::hex << class_ledger.state << std::dec << '\n'
                  << "TARGET_RESPONSE_CLASSES " << target_classes.size()
                  << " UNION_POPCOUNT " << std::popcount(target_union)
                  << " MINIMUM " << target_cover.minimum << " MASKS";
        for (u32 mask : target_cover.masks)
            std::cout << ' ' << std::hex << mask << std::dec;
        std::cout << '\n'
                  << "FIBREWIDE_MASKS " << fibrewide
                  << " RESPONSE_CLASSES " << fibre_classes.size()
                  << " UNION_POPCOUNT " << std::popcount(fibre_union)
                  << " MINIMUM " << fibre_cover.minimum << " MASKS";
        for (u32 mask : fibre_cover.masks)
            std::cout << ' ' << std::hex << mask << std::dec;
        std::cout << "\nVERDICT PASS COMPLETE_TWO_MASK_RESPONSE_QUOTIENT\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "TWO_MASK_CEGAR_ERROR " << error.what() << '\n';
        return 1;
    }
}
