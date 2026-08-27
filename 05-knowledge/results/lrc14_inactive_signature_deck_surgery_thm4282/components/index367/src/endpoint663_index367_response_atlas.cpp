#define RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
#include "04-computation/lrc14_joint421_global_common_carrier_thm4281/response_pattern_atlas.cpp"
#undef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY

#include <fstream>
#include <iomanip>
#include <set>

namespace {

constexpr std::size_t TARGET_INDEX = 367;
constexpr u32 TARGET_MASK = UINT32_C(0x02188125);
constexpr u64 EXPECTED_DECK_FNV = UINT64_C(0x20d63dd42fe8150e);
constexpr std::array<u32, 3> PRIVATE_BODIES = {
    UINT32_C(0x05646408), UINT32_C(0x25065480), UINT32_C(0x35a24080)};
constexpr u32 FULL_RESPONSE = (u32{1} << PRIVATE_BODIES.size()) - 1;

u64 mask_fnv(const std::vector<u32>& masks) {
    FnvLocal ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

std::vector<u32> read_deck(const std::string& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open deck");
    std::vector<u32> deck;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        std::size_t used = 0;
        const u64 wide = std::stoull(token, &used, 16);
        require(used == token.size() && wide < (UINT64_C(1) << 30),
                "invalid deck token");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "invalid deck mask");
        deck.push_back(mask);
    }
    require(deck.size() == 421 && mask_fnv(deck) == EXPECTED_DECK_FNV &&
                deck[TARGET_INDEX] == TARGET_MASK,
            "deck identity changed");
    return deck;
}

u64 scan_body_cover(const std::vector<u32>& deck, u64& checks,
                    u64& maximum_checks, u64& failures, u64& row_fnv) {
    FnvLocal covered_ledger;
    FnvLocal row_ledger;
    checks = 0;
    maximum_checks = 0;
    failures = 0;
    u64 bodies = 0;
    u32 body = (u32{1} << 9) - 1;
    while (body < (u32{1} << 30)) {
        ++bodies;
        u64 tested = 0;
        u32 witness = 0;
        for (u32 repair : deck) {
            ++tested;
            if ((body & repair) == 0) {
                witness = repair;
                break;
            }
        }
        checks += tested;
        maximum_checks = std::max(maximum_checks, tested);
        failures += witness == 0;
        row_ledger.add(body);
        row_ledger.add(tested);
        row_ledger.add(witness);
        if (witness != 0) {
            covered_ledger.add(body);
            covered_ledger.add(witness);
        }
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    require(bodies == UINT64_C(14307150), "body enumeration changed");
    row_fnv = row_ledger.state;
    return covered_ledger.state;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 2, "usage: endpoint663-atlas JOINT_DECK");
        init_choose8_local();
        const std::vector<u32> deck = read_deck(argv[1]);
        FnvLocal failure_ledger;
        u32 union_mask = 0;
        for (u32 body : PRIVATE_BODIES) {
            require(std::popcount(body) == 9 && (body & TARGET_MASK) == 0,
                    "private body identity changed");
            for (std::size_t index = 0; index < deck.size(); ++index) {
                if (index == TARGET_INDEX) continue;
                require((body & deck[index]) != 0,
                        "private body has retained-deck witness");
            }
            failure_ledger.add(body);
            union_mask |= body;
        }
        require(failure_ledger.state == UINT64_C(0x22bd212c1ffec6e2),
                "private-body ledger changed");

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool cells changed");
        const ActiveUniverse active = build_active_universe(cells, 366, 663);
        const u64 target_rank = colex_rank8_local(TARGET_MASK);
        require(target_rank < EXPECTED_REPAIRS && !active.active[target_rank],
                "deleted target unexpectedly active");
        u64 deck_inactive = 0;
        u64 deck_active = 0;
        FnvLocal inactive_ledger;
        for (std::size_t index = 0; index < deck.size(); ++index) {
            const bool is_active = active.active[colex_rank8_local(deck[index])];
            deck_active += is_active;
            deck_inactive += !is_active;
            if (!is_active) {
                inactive_ledger.add(index);
                inactive_ledger.add(deck[index]);
            }
        }
        require(deck_active == 420 && deck_inactive == 1,
                "pair singleton-inactive signature changed");

        std::vector<u32> response(EXPECTED_REPAIRS, 0);
        u64 active_incidences = 0;
        for (std::size_t obligation = 0; obligation < PRIVATE_BODIES.size();
             ++obligation) {
            enumerate_disjoint_repairs(
                PRIVATE_BODIES[obligation], [&](u32, u64 rank) {
                    if (!active.active[rank]) return;
                    response[rank] |= u32{1} << obligation;
                    ++active_incidences;
                });
        }

        std::array<u64, 1 << PRIVATE_BODIES.size()> multiplicity{};
        std::array<u32, 1 << PRIVATE_BODIES.size()> least{};
        FnvLocal nonempty_ledger;
        u64 nonempty = 0;
        std::vector<u32> full_response_masks;
        FnvLocal full_response_ledger;
        for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
            const u32 pattern = response[rank];
            ++multiplicity[pattern];
            const u32 repair = unrank_colex8(rank);
            if (least[pattern] == 0 || repair < least[pattern])
                least[pattern] = repair;
            if (pattern != 0) {
                ++nonempty;
                nonempty_ledger.add(repair);
                nonempty_ledger.add(pattern);
            }
            if (pattern == FULL_RESPONSE) {
                full_response_masks.push_back(repair);
                full_response_ledger.add(repair);
            }
        }
        require(std::accumulate(multiplicity.begin(), multiplicity.end(),
                                u64{0}) == EXPECTED_REPAIRS,
                "response multiplicity changed");
        FnvLocal class_ledger;
        std::vector<u32> present;
        for (u32 pattern = 1; pattern <= FULL_RESPONSE; ++pattern) {
            if (multiplicity[pattern] == 0) continue;
            present.push_back(pattern);
            class_ledger.add(pattern);
            class_ledger.add(multiplicity[pattern]);
            class_ledger.add(least[pattern]);
        }
        std::vector<u32> maximal;
        FnvLocal maximal_ledger;
        for (u32 pattern : present) {
            bool dominated = false;
            for (u32 other : present) {
                if (pattern != other && (pattern | other) == other) {
                    dominated = true;
                    break;
                }
            }
            if (!dominated) {
                maximal.push_back(pattern);
                maximal_ledger.add(pattern);
                maximal_ledger.add(multiplicity[pattern]);
                maximal_ledger.add(least[pattern]);
            }
        }

        constexpr int INF = 1000;
        std::array<int, 1 << PRIVATE_BODIES.size()> distance{};
        std::array<u32, 1 << PRIVATE_BODIES.size()> previous{};
        std::array<u32, 1 << PRIVATE_BODIES.size()> previous_pattern{};
        distance.fill(INF);
        distance[0] = 0;
        for (u32 covered = 0; covered <= FULL_RESPONSE; ++covered) {
            if (distance[covered] == INF) continue;
            for (u32 pattern : maximal) {
                const u32 next = covered | pattern;
                if (distance[next] > distance[covered] + 1) {
                    distance[next] = distance[covered] + 1;
                    previous[next] = covered;
                    previous_pattern[next] = pattern;
                }
            }
        }
        require(distance[FULL_RESPONSE] < INF,
                "private-body response has no active cover");
        std::vector<u32> cover_patterns;
        for (u32 state = FULL_RESPONSE; state != 0; state = previous[state]) {
            require(previous_pattern[state] != 0, "cover predecessor missing");
            cover_patterns.push_back(previous_pattern[state]);
        }
        std::vector<u32> cover_masks;
        for (u32 pattern : cover_patterns) cover_masks.push_back(least[pattern]);
        std::sort(cover_masks.begin(), cover_masks.end());
        require(cover_masks.size() ==
                    static_cast<std::size_t>(distance[FULL_RESPONSE]),
                "cover recovery changed");

        std::vector<u32> rebuilt = deck;
        rebuilt.erase(rebuilt.begin() + TARGET_INDEX);
        rebuilt.insert(rebuilt.end(), cover_masks.begin(), cover_masks.end());
        require(std::set<u32>(rebuilt.begin(), rebuilt.end()).size() ==
                    rebuilt.size(),
                "rebuilt deck has duplicate repair");
        u64 body_checks = 0;
        u64 maximum_checks = 0;
        u64 body_failures = 0;
        u64 body_row_fnv = 0;
        const u64 body_cover_fnv = scan_body_cover(
            rebuilt, body_checks, maximum_checks, body_failures, body_row_fnv);
        require(body_failures == 0, "rebuilt deck does not cover every body");

        FnvLocal cover_ledger;
        for (u32 mask : cover_masks) cover_ledger.add(mask);
        require(inactive_ledger.state == UINT64_C(0xeb3f16a1b127b399) &&
                    union_mask == UINT32_C(0x35e67488) &&
                    std::popcount(union_mask) == 15 &&
                    active.count == UINT64_C(3095104) &&
                    active.fnv == UINT64_C(0xb6646ef482f50d9c) &&
                    active.zeta_operations == UINT64_C(152170690) &&
                    active_incidences == UINT64_C(77172) &&
                    nonempty == UINT64_C(73605) &&
                    multiplicity == std::array<u64, 8>{
                        UINT64_C(5779320), UINT64_C(18559),
                        UINT64_C(19210), UINT64_C(719),
                        UINT64_C(32330), UINT64_C(266),
                        UINT64_C(2460), UINT64_C(61)} &&
                    nonempty_ledger.state == UINT64_C(0x41818be95d3b963d) &&
                    present == std::vector<u32>{1,2,3,4,5,6,7} &&
                    class_ledger.state == UINT64_C(0xeb621f2bb7fc9ee3) &&
                    maximal == std::vector<u32>{7} &&
                    maximal_ledger.state == UINT64_C(0xda850b41adc17035) &&
                    full_response_masks.size() == 61 &&
                    full_response_ledger.state ==
                        UINT64_C(0x1b2bd26ec685ee47) &&
                    distance[FULL_RESPONSE] == 1 &&
                    cover_masks == std::vector<u32>{UINT32_C(0x02108325)} &&
                    cover_ledger.state == UINT64_C(0x44fa6227992e30bb) &&
                    rebuilt.size() == 421 &&
                    mask_fnv(rebuilt) == UINT64_C(0x5072d6c29f0d8a08) &&
                    body_checks == UINT64_C(405169849) &&
                    maximum_checks == UINT64_C(421) && body_failures == 0 &&
                    body_cover_fnv == UINT64_C(0x1ac5ced486c204b9) &&
                    body_row_fnv == UINT64_C(0x98f96181bcfbeff8),
                "frozen endpoint663 index367 response atlas changed");
        std::cout << "ENDPOINT663_INDEX367_RESPONSE_ATLAS_V1\n"
                  << "PAIR 366,663 TARGET_INDEX " << TARGET_INDEX
                  << " TARGET_MASK " << std::hex << TARGET_MASK << std::dec
                  << " DECK_ACTIVE " << deck_active << " DECK_INACTIVE "
                  << deck_inactive << " INACTIVE_FNV " << std::hex
                  << inactive_ledger.state << std::dec << '\n'
                  << "PRIVATE_BODIES " << PRIVATE_BODIES.size()
                  << " FAILURE_FNV " << std::hex << failure_ledger.state
                  << " UNION " << union_mask << std::dec
                  << " UNION_SIZE " << std::popcount(union_mask) << '\n'
                  << "REPAIRS " << EXPECTED_REPAIRS << " ACTIVE "
                  << active.count << " ACTIVE_FNV " << std::hex << active.fnv
                  << std::dec << " ZETA_OPERATIONS "
                  << active.zeta_operations << '\n'
                  << "COMPLEMENT_CHECKS "
                  << PRIVATE_BODIES.size() * DISJOINT_REPAIRS_PER_BODY
                  << " ACTIVE_INCIDENCES " << active_incidences
                  << " NONEMPTY_REPAIRS " << nonempty
                  << " EMPTY_REPAIRS " << multiplicity[0]
                  << " RESPONSE_FNV " << std::hex << nonempty_ledger.state
                  << std::dec << '\n'
                  << "NONEMPTY_CLASSES " << present.size()
                  << " CLASS_FNV " << std::hex << class_ledger.state
                  << std::dec << " MAXIMAL_CLASSES " << maximal.size()
                  << " MAXIMAL_FNV " << std::hex << maximal_ledger.state
                  << std::dec << '\n';
        for (u32 pattern : present) {
            std::cout << "PATTERN " << std::hex << pattern << " LEAST "
                      << least[pattern] << std::dec << " MULTIPLICITY "
                      << multiplicity[pattern] << " COVER "
                      << std::popcount(pattern) << '\n';
        }
        std::cout << "FULL_RESPONSE_CLASS " << full_response_masks.size()
                  << " FNV " << std::hex << full_response_ledger.state
                  << std::dec << " MASKS_HEX ";
        for (std::size_t index = 0; index < full_response_masks.size(); ++index) {
            if (index) std::cout << ',';
            std::cout << std::hex << std::setw(8) << std::setfill('0')
                      << full_response_masks[index] << std::dec
                      << std::setfill(' ');
        }
        std::cout << '\n';
        std::cout << "MINIMUM_REPLACEMENTS " << distance[FULL_RESPONSE]
                  << " COVER_MASKS ";
        for (std::size_t index = 0; index < cover_masks.size(); ++index) {
            if (index) std::cout << ',';
            std::cout << std::hex << cover_masks[index] << std::dec;
        }
        std::cout << " COVER_FNV " << std::hex << cover_ledger.state
                  << std::dec << '\n'
                  << "REBUILT_DECK " << rebuilt.size() << " FNV " << std::hex
                  << mask_fnv(rebuilt) << std::dec << " BODY_CHECKS "
                  << body_checks << " MAX_CHECKS " << maximum_checks
                  << " FAILURES " << body_failures << " BODY_COVER_FNV "
                  << std::hex << body_cover_fnv << " BODY_ROW_FNV "
                  << body_row_fnv << std::dec << '\n'
                  << "VERDICT PASS EXACT_COMPLETE_ACTIVE_RESPONSE_QUOTIENT_AND_REBUILT_BODY_COVER\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT663_INDEX367_ATLAS_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}
