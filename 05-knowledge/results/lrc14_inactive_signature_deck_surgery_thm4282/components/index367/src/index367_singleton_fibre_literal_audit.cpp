#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#pragma clang diagnostic ignored "-Wunneeded-internal-declaration"
#endif
#define CARRIER_CEGAR_LIBRARY_ONLY
#include "04-computation/lrc14_three_round_learned_carrier_thm4266/carrier_cegar_descent.cpp"
#undef CARRIER_CEGAR_LIBRARY_ONLY
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#include <fstream>
#include <iomanip>
#include <set>

namespace {

constexpr std::size_t TARGET_INDEX = 367;
constexpr u32 TARGET_MASK = UINT32_C(0x02188125);
constexpr u64 EXPECTED_DECK_FNV = UINT64_C(0x20d63dd42fe8150e);
constexpr std::array<u32, 61> FULL_RESPONSE_CANDIDATES = {{
    0x02108325,0x02188207,0x02188213,0x02188a03,0x02198203,0x02198a01,
    0x0a108125,0x0a108207,0x0a108213,0x0a108321,0x0a108807,0x0a108813,
    0x0a108903,0x0a108921,0x0a108a03,0x0a108a05,0x0a108a11,0x0a108a21,
    0x0a108a41,0x0a108b01,0x0a118203,0x0a118205,0x0a118211,0x0a118221,
    0x0a118241,0x0a118301,0x0a118803,0x0a118a01,0x0a188007,0x0a188013,
    0x0a188016,0x0a188023,0x0a188043,0x0a188103,0x0a188121,0x0a188203,
    0x0a188205,0x0a188206,0x0a188211,0x0a188212,0x0a188221,0x0a188241,
    0x0a188301,0x0a188803,0x0a188805,0x0a188806,0x0a188811,0x0a188812,
    0x0a188821,0x0a188841,0x0a188842,0x0a188901,0x0a188a01,0x0a188a02,
    0x0a198003,0x0a198101,0x0a198201,0x0a198202,0x0a198801,0x0a198802,
    0x0a198a00
}};
constexpr std::array<std::pair<int, int>, 26> FIBRE = {{
    {17,112}, {31,112}, {45,437}, {62,112}, {67,109}, {75,112},
    {75,113}, {90,332}, {109,122}, {109,127}, {109,151}, {109,302},
    {109,381}, {109,462}, {109,516}, {112,127}, {112,144}, {112,151},
    {112,274}, {112,302}, {127,428}, {127,462}, {127,518}, {227,370},
    {308,381}, {366,663}
}};

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
    for (u32 candidate : FULL_RESPONSE_CANDIDATES)
        require(std::popcount(candidate) == 8 && !distinct.contains(candidate),
                "replacement candidate invalid/already present");
    return deck;
}

struct IndexedRepairLiteral {
    u32 mask = 0;
    std::vector<std::pair<std::uint16_t, std::uint16_t>> components;
};

IndexedRepairLiteral index_repair_literal(u32 mask,
                                          const std::vector<Cell>& cells) {
    IndexedRepairLiteral out;
    out.mask = mask;
    std::size_t index = 0;
    while (index < cells.size()) {
        if ((cells[index].failed_pool & ~mask) != 0) {
            ++index;
            continue;
        }
        const std::size_t begin = index;
        do {
            ++index;
        } while (index < cells.size() &&
                 (cells[index].failed_pool & ~mask) == 0);
        require(begin < UINT16_MAX && index < UINT16_MAX,
                "pool wall index exceeds u16");
        out.components.push_back({static_cast<std::uint16_t>(begin),
                                  static_cast<std::uint16_t>(index)});
    }
    require(!out.components.empty(), "empty repair geometry");
    return out;
}

i128 literal_lcm(i128 left, i128 right) {
    require(left > 0 && right > 0, "nonpositive literal lcm input");
    return left / gcd_i128(left, right) * right;
}

struct LiteralArc {
    i128 left = 0;
    i128 right = 0;
    i128 prefix = 0;
};

struct LiteralPair {
    i128 grid = 0;
    i128 pool_scale = 0;
    std::vector<LiteralArc> arcs;
};

std::vector<std::pair<i128, i128>> speed_intervals_literal(int speed,
                                                            i128 grid) {
    const i128 denominator = static_cast<i128>(14) * speed;
    require(speed > 0 && grid % denominator == 0,
            "literal grid not speed-divisible");
    const i128 unit = grid / denominator;
    std::vector<std::pair<i128, i128>> arcs;
    arcs.reserve(speed);
    for (int tooth = 0; tooth < speed; ++tooth) {
        arcs.push_back({static_cast<i128>(14 * tooth + 1) * unit,
                        static_cast<i128>(14 * tooth + 13) * unit});
    }
    return arcs;
}

LiteralPair build_literal_pair(int q, int r) {
    require(q > 0 && q < r, "invalid literal pair");
    LiteralPair out;
    out.grid = literal_lcm(literal_lcm(static_cast<i128>(COMMON),
                                      static_cast<i128>(14) * q),
                           static_cast<i128>(14) * r);
    require(out.grid % COMMON == 0, "literal grid lost pool scale");
    out.pool_scale = out.grid / COMMON;
    const auto q_arcs = speed_intervals_literal(q, out.grid);
    const auto r_arcs = speed_intervals_literal(r, out.grid);
    std::size_t i = 0;
    std::size_t j = 0;
    i128 running = 0;
    while (i < q_arcs.size() && j < r_arcs.size()) {
        const i128 left = std::max(q_arcs[i].first, r_arcs[j].first);
        const i128 right = std::min(q_arcs[i].second, r_arcs[j].second);
        if (left < right) {
            if (!out.arcs.empty() && out.arcs.back().right == left) {
                out.arcs.back().right = right;
                running += right - left;
            } else {
                out.arcs.push_back({left, right, running});
                running += right - left;
            }
        }
        if (q_arcs[i].second < r_arcs[j].second) {
            ++i;
        } else if (r_arcs[j].second < q_arcs[i].second) {
            ++j;
        } else {
            ++i;
            ++j;
        }
    }
    require(!out.arcs.empty() && running > 0,
            "literal pair safe intersection empty");
    return out;
}

i128 literal_prefix(const LiteralPair& pair, i128 tick) {
    require(tick >= 0 && tick <= pair.grid,
            "literal prefix outside circle");
    std::size_t low = 0;
    std::size_t high = pair.arcs.size();
    while (low < high) {
        const std::size_t middle = low + (high - low) / 2;
        if (pair.arcs[middle].right <= tick) low = middle + 1;
        else high = middle;
    }
    if (low == pair.arcs.size()) {
        const LiteralArc& last = pair.arcs.back();
        return last.prefix + last.right - last.left;
    }
    const LiteralArc& arc = pair.arcs[low];
    return arc.prefix +
           (tick > arc.left ? std::min(tick, arc.right) - arc.left : 0);
}

i128 repair_margin(const IndexedRepairLiteral& repair,
                   const std::vector<i128>& prefix, i128 grid) {
    i128 mass = 0;
    for (const auto [left, right] : repair.components)
        mass += prefix[right] - prefix[left];
    return static_cast<i128>(63) * mass - static_cast<i128>(4) * grid;
}

bool rational_less_nonnegative(i128 a, i128 b, i128 c, i128 d) {
    require(a >= 0 && c >= 0 && b > 0 && d > 0,
            "invalid rational comparison");
    bool reversed = false;
    while (true) {
        const i128 aq = a / b;
        const i128 cq = c / d;
        if (aq != cq) return reversed ? aq > cq : aq < cq;
        a %= b;
        c %= d;
        if (a == 0 || c == 0) {
            if (a == 0 && c == 0) return false;
            const bool raw_less = a == 0;
            return reversed ? !raw_less : raw_less;
        }
        std::swap(a, b);
        std::swap(c, d);
        reversed = !reversed;
    }
}

void add_i128(FnvLocal& ledger, i128 value) {
    const __uint128_t bits = static_cast<__uint128_t>(value);
    ledger.add(static_cast<u64>(bits));
    ledger.add(static_cast<u64>(bits >> 64));
}

std::string pool_labels(u32 mask) {
    std::ostringstream out;
    bool first = true;
    for (std::size_t bit = 0; bit < POOL.size(); ++bit) {
        if ((mask & (u32{1} << bit)) == 0) continue;
        if (!first) out << ',';
        first = false;
        out << POOL[bit];
    }
    return out.str();
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
        require(argc == 2, "usage: fibre-literal-audit JOINT_DECK");
        const std::vector<u32> deck = read_deck(argv[1]);
        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool cells changed");
        std::vector<i64> walls(cells.size() + 1);
        walls.front() = cells.front().left;
        for (std::size_t index = 0; index < cells.size(); ++index) {
            require(cells[index].left == walls[index],
                    "pool cells not contiguous");
            walls[index + 1] = cells[index].right;
        }
        require(walls.front() == 0 && walls.back() == COMMON,
                "pool walls changed");
        std::vector<IndexedRepairLiteral> geometry;
        geometry.reserve(deck.size());
        for (u32 mask : deck)
            geometry.push_back(index_repair_literal(mask, cells));
        std::vector<IndexedRepairLiteral> candidate_geometry;
        candidate_geometry.reserve(FULL_RESPONSE_CANDIDATES.size());
        for (u32 candidate : FULL_RESPONSE_CANDIDATES)
            candidate_geometry.push_back(index_repair_literal(candidate, cells));

        FnvLocal fibre_ledger;
        std::array<i128, FIBRE.size()> grids{};
        std::array<i128, FIBRE.size()> target_margins{};
        std::array<std::array<i128, FULL_RESPONSE_CANDIDATES.size()>,
                   FIBRE.size()> candidate_margins{};
        std::array<u64, FULL_RESPONSE_CANDIDATES.size()> active_rows{};
        std::array<bool, FULL_RESPONSE_CANDIDATES.size()> all_active{};
        all_active.fill(true);
        u64 equalities = 0;
        std::cout << "INDEX367_FIBRE_LITERAL_AUDIT_V1\n";
        for (std::size_t row = 0; row < FIBRE.size(); ++row) {
            const auto [q, r] = FIBRE[row];
            fibre_ledger.add(q);
            fibre_ledger.add(r);
            const LiteralPair pair = build_literal_pair(q, r);
            grids[row] = pair.grid;
            std::vector<i128> prefix(walls.size());
            for (std::size_t wall = 0; wall < walls.size(); ++wall)
                prefix[wall] = literal_prefix(
                    pair, static_cast<i128>(walls[wall]) * pair.pool_scale);
            u64 inactive = 0;
            std::size_t inactive_index = deck.size();
            u32 inactive_mask = 0;
            for (std::size_t index = 0; index < geometry.size(); ++index) {
                const i128 margin = repair_margin(geometry[index], prefix,
                                                  pair.grid);
                equalities += margin == 0;
                if (margin < 0) {
                    ++inactive;
                    inactive_index = index;
                    inactive_mask = geometry[index].mask;
                }
                if (index == TARGET_INDEX) target_margins[row] = margin;
            }
            require(inactive == 1 && inactive_index == TARGET_INDEX &&
                        inactive_mask == TARGET_MASK &&
                        target_margins[row] < 0,
                    "literal singleton fibre changed");
            for (std::size_t candidate = 0;
                 candidate < candidate_geometry.size(); ++candidate) {
                const i128 margin = repair_margin(
                    candidate_geometry[candidate], prefix, pair.grid);
                candidate_margins[row][candidate] = margin;
                equalities += margin == 0;
                active_rows[candidate] += margin >= 0;
                all_active[candidate] = all_active[candidate] && margin >= 0;
            }
        }
        require(equalities == 0, "literal fibre equality changed");

        std::vector<std::size_t> survivors;
        FnvLocal survivor_ledger;
        FnvLocal candidate_activity_ledger;
        for (std::size_t candidate = 0;
             candidate < FULL_RESPONSE_CANDIDATES.size(); ++candidate) {
            candidate_activity_ledger.add(FULL_RESPONSE_CANDIDATES[candidate]);
            candidate_activity_ledger.add(active_rows[candidate]);
            candidate_activity_ledger.add(all_active[candidate]);
            if (!all_active[candidate]) continue;
            survivors.push_back(candidate);
            survivor_ledger.add(FULL_RESPONSE_CANDIDATES[candidate]);
        }
        require(!survivors.empty(),
                "no full-response replacement survives singleton fibre");
        const std::size_t chosen_index = *std::min_element(
            survivors.begin(), survivors.end(), [](std::size_t left,
                                                    std::size_t right) {
                return FULL_RESPONSE_CANDIDATES[left] <
                       FULL_RESPONSE_CANDIDATES[right];
            });
        const u32 replacement_mask = FULL_RESPONSE_CANDIDATES[chosen_index];

        FnvLocal row_ledger;
        bool minimum_set = false;
        std::pair<int, int> minimum_pair;
        i128 minimum_margin = 0;
        i128 minimum_grid = 1;
        for (std::size_t row = 0; row < FIBRE.size(); ++row) {
            const auto [q, r] = FIBRE[row];
            const i128 margin = candidate_margins[row][chosen_index];
            require(margin > 0, "chosen replacement not strictly active");
            if (!minimum_set || rational_less_nonnegative(
                    margin, grids[row], minimum_margin, minimum_grid)) {
                minimum_set = true;
                minimum_pair = {q, r};
                minimum_margin = margin;
                minimum_grid = grids[row];
            }
            row_ledger.add(q);
            row_ledger.add(r);
            row_ledger.add(TARGET_INDEX);
            row_ledger.add(TARGET_MASK);
            add_i128(row_ledger, target_margins[row]);
            row_ledger.add(replacement_mask);
            add_i128(row_ledger, margin);
            add_i128(row_ledger, grids[row]);
            const i128 target_gcd = gcd_i128(-target_margins[row], grids[row]);
            const i128 replacement_gcd = gcd_i128(margin, grids[row]);
            std::cout << "ROW " << q << ',' << r
                      << " INACTIVE_INDEX " << TARGET_INDEX
                      << " INACTIVE_MASK " << std::hex << TARGET_MASK
                      << std::dec << " TARGET_MARGIN_NUM "
                      << decimal(target_margins[row]) << " REPLACEMENT "
                      << std::hex << replacement_mask << std::dec
                      << " REPLACEMENT_MARGIN_NUM " << decimal(margin)
                      << " DEN " << decimal(grids[row])
                      << " TARGET_NORMALIZED "
                      << decimal(target_margins[row] / target_gcd) << '/'
                      << decimal(grids[row] / target_gcd)
                      << " REPLACEMENT_NORMALIZED "
                      << decimal(margin / replacement_gcd) << '/'
                      << decimal(grids[row] / replacement_gcd) << '\n';
        }
        require(minimum_set, "replacement gap minimum missing");

        std::vector<u32> rebuilt = deck;
        rebuilt.erase(rebuilt.begin() + TARGET_INDEX);
        rebuilt.push_back(replacement_mask);
        require(rebuilt.size() == 421 &&
                    std::set<u32>(rebuilt.begin(), rebuilt.end()).size() == 421,
                "rebuilt deck identity changed");
        u64 body_checks = 0;
        u64 maximum_checks = 0;
        u64 body_failures = 0;
        u64 body_row_fnv = 0;
        const u64 body_cover_fnv = scan_body_cover(
            rebuilt, body_checks, maximum_checks, body_failures, body_row_fnv);
        require(body_failures == 0, "literal rebuilt deck fails body cover");
        require(survivors.size() == 2 &&
                    FULL_RESPONSE_CANDIDATES[survivors[0]] ==
                        UINT32_C(0x0a188803) &&
                    FULL_RESPONSE_CANDIDATES[survivors[1]] ==
                        UINT32_C(0x0a188a01) &&
                    survivor_ledger.state == UINT64_C(0x9ecd12e99ccfb0a1) &&
                    candidate_activity_ledger.state ==
                        UINT64_C(0x4da2ad8141f0fae6) &&
                    replacement_mask == UINT32_C(0x0a188803) &&
                    fibre_ledger.state == UINT64_C(0x9758d0de94c51f75) &&
                    row_ledger.state == UINT64_C(0x73f3ca7c94282e97) &&
                    equalities == 0 &&
                    minimum_pair == std::pair{227, 370} &&
                    minimum_margin == static_cast<i128>(2213184690321570LL) &&
                    minimum_grid == static_cast<i128>(153207497939015520LL) &&
                    mask_fnv(rebuilt) == UINT64_C(0x87b42cf8a2069177) &&
                    body_checks == UINT64_C(405169849) &&
                    maximum_checks == UINT64_C(421) && body_failures == 0 &&
                    body_cover_fnv == UINT64_C(0xe067558cb36fd4e6) &&
                    body_row_fnv == UINT64_C(0xdf7a749839a21c83),
                "frozen detached singleton-fibre audit changed");
        std::cout << "FULL_RESPONSE_CANDIDATES "
                  << FULL_RESPONSE_CANDIDATES.size() << " SURVIVORS "
                  << survivors.size() << " SURVIVOR_FNV " << std::hex
                  << survivor_ledger.state << std::dec << '\n';
        for (std::size_t candidate = 0;
             candidate < FULL_RESPONSE_CANDIDATES.size(); ++candidate) {
            std::cout << "CANDIDATE " << std::hex
                      << FULL_RESPONSE_CANDIDATES[candidate] << std::dec
                      << " ACTIVE_ROWS " << active_rows[candidate]
                      << " ALL_ACTIVE "
                      << (all_active[candidate] ? "YES" : "NO") << '\n';
        }
        const std::size_t top_row = FIBRE.size() - 1;
        const i128 local_margin = candidate_margins[top_row][0];
        const i128 local_gcd = gcd_i128(local_margin, grids[top_row]);
        const i128 top_target_gcd =
            gcd_i128(-target_margins[top_row], grids[top_row]);
        require(FIBRE[top_row] == std::pair{366, 663} && local_margin > 0,
                "top local exchange changed");
        std::cout << "CANDIDATE_ACTIVITY_FNV " << std::hex
                  << candidate_activity_ledger.state << std::dec << '\n'
                  << "TOP_LOCAL_ONE_SWAP PAIR 366,663 REPLACEMENT "
                  << std::hex << FULL_RESPONSE_CANDIDATES[0] << std::dec
                  << " REMOVE_BITS 19 REMOVE_LABELS 145"
                  << " ADD_BITS 9 ADD_LABELS 63"
                  << " TARGET_MARGIN_NUM " << decimal(target_margins[top_row])
                  << " REPLACEMENT_MARGIN_NUM " << decimal(local_margin)
                  << " DEN " << decimal(grids[top_row])
                  << " TARGET_NORMALIZED "
                  << decimal(target_margins[top_row] / top_target_gcd) << '/'
                  << decimal(grids[top_row] / top_target_gcd)
                  << " REPLACEMENT_NORMALIZED "
                  << decimal(local_margin / local_gcd) << '/'
                  << decimal(grids[top_row] / local_gcd) << '\n'
                  << "CHOSEN_REPLACEMENT " << std::hex << replacement_mask
                  << " REMOVE_MASK " << (TARGET_MASK & ~replacement_mask)
                  << " ADD_MASK " << (replacement_mask & ~TARGET_MASK)
                  << std::dec << " REMOVE_LABELS "
                  << pool_labels(TARGET_MASK & ~replacement_mask)
                  << " ADD_LABELS "
                  << pool_labels(replacement_mask & ~TARGET_MASK) << '\n'
                  << "FIBRE_ROWS " << FIBRE.size() << " FIBRE_FNV "
                  << std::hex << fibre_ledger.state << " ROW_FNV "
                  << row_ledger.state << std::dec << " EQUALITIES "
                  << equalities << '\n'
                  << "MIN_REPLACEMENT_GAP PAIR " << minimum_pair.first << ','
                  << minimum_pair.second << " MARGIN_NUM "
                  << decimal(minimum_margin) << " DEN "
                  << decimal(minimum_grid) << '\n'
                  << "REBUILT_DECK " << rebuilt.size() << " FNV " << std::hex
                  << mask_fnv(rebuilt) << std::dec << " BODY_CHECKS "
                  << body_checks << " MAX_CHECKS " << maximum_checks
                  << " FAILURES " << body_failures << " BODY_COVER_FNV "
                  << std::hex << body_cover_fnv << " BODY_ROW_FNV "
                  << body_row_fnv << std::dec << '\n'
                  << "VERDICT PASS DETACHED_LITERAL_SINGLETON_FIBRE_AND_BODY_COVER\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "INDEX367_FIBRE_LITERAL_ERROR " << error.what() << '\n';
        return 1;
    }
}
