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
                "signature popcount changed");
        rows.push_back(row);
    }
    require(rows.size() == 24223, "signature universe changed");
    return rows;
}

std::array<u64, 7> pair_signature(std::size_t a, std::size_t b) {
    std::array<u64, 7> out{};
    out[a / 64] |= UINT64_C(1) << (a % 64);
    out[b / 64] |= UINT64_C(1) << (b % 64);
    return out;
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
                "invalid deck mask");
        deck.push_back(mask);
    }
    require(deck.size() == 421 &&
                mask_fnv(deck) == UINT64_C(0x20d63dd42fe8150e),
            "deck identity changed");
    return deck;
}

struct IndexedRepair {
    u32 mask = 0;
    std::vector<std::pair<std::uint16_t, std::uint16_t>> components;
};

IndexedRepair index_repair(u32 mask, const std::vector<Cell>& cells) {
    IndexedRepair out;
    out.mask = mask;
    std::size_t i = 0;
    while (i < cells.size()) {
        if ((cells[i].failed_pool & ~mask) != 0) {
            ++i;
            continue;
        }
        const std::size_t begin = i;
        do {
            ++i;
        } while (i < cells.size() && (cells[i].failed_pool & ~mask) == 0);
        require(begin < UINT16_MAX && i < UINT16_MAX, "component overflow");
        out.components.push_back({static_cast<std::uint16_t>(begin),
                                  static_cast<std::uint16_t>(i)});
    }
    require(!out.components.empty(), "empty repair geometry");
    return out;
}

i128 lcm128(i128 a, i128 b) { return a / gcd_i128(a, b) * b; }

struct Arc {
    i128 left = 0;
    i128 right = 0;
    i128 prefix = 0;
};

struct LiteralPair {
    i128 grid = 0;
    i128 pool_scale = 0;
    std::vector<Arc> arcs;
};

std::vector<std::pair<i128, i128>> speed_intervals(int speed, i128 grid) {
    const i128 den = static_cast<i128>(14) * speed;
    require(speed > 0 && grid % den == 0, "bad speed grid");
    const i128 unit = grid / den;
    std::vector<std::pair<i128, i128>> out;
    for (int tooth = 0; tooth < speed; ++tooth)
        out.push_back({static_cast<i128>(14 * tooth + 1) * unit,
                       static_cast<i128>(14 * tooth + 13) * unit});
    return out;
}

LiteralPair build_literal_pair(int q, int r) {
    require(q > 0 && q < r, "bad literal pair");
    LiteralPair out;
    out.grid = lcm128(lcm128(static_cast<i128>(COMMON),
                            static_cast<i128>(14) * q),
                      static_cast<i128>(14) * r);
    out.pool_scale = out.grid / COMMON;
    const auto qa = speed_intervals(q, out.grid);
    const auto ra = speed_intervals(r, out.grid);
    std::size_t i = 0, j = 0;
    i128 running = 0;
    while (i < qa.size() && j < ra.size()) {
        const i128 left = std::max(qa[i].first, ra[j].first);
        const i128 right = std::min(qa[i].second, ra[j].second);
        if (left < right) {
            if (!out.arcs.empty() && out.arcs.back().right == left) {
                out.arcs.back().right = right;
                running += right - left;
            } else {
                out.arcs.push_back({left, right, running});
                running += right - left;
            }
        }
        if (qa[i].second < ra[j].second) ++i;
        else if (ra[j].second < qa[i].second) ++j;
        else { ++i; ++j; }
    }
    require(!out.arcs.empty() && running > 0, "empty literal pair");
    return out;
}

i128 prefix_at(const LiteralPair& pair, i128 tick) {
    std::size_t lo = 0, hi = pair.arcs.size();
    while (lo < hi) {
        const std::size_t mid = lo + (hi - lo) / 2;
        if (pair.arcs[mid].right <= tick) lo = mid + 1;
        else hi = mid;
    }
    if (lo == pair.arcs.size()) {
        const Arc& a = pair.arcs.back();
        return a.prefix + a.right - a.left;
    }
    const Arc& a = pair.arcs[lo];
    return a.prefix +
           (tick > a.left ? std::min(tick, a.right) - a.left : 0);
}

i128 repair_margin(const IndexedRepair& repair,
                   const std::vector<i128>& prefix, i128 grid) {
    i128 mass = 0;
    for (const auto [left, right] : repair.components)
        mass += prefix[right] - prefix[left];
    return static_cast<i128>(63) * mass - static_cast<i128>(4) * grid;
}

void add_i128(FnvLocal& ledger, i128 value) {
    const __uint128_t bits = static_cast<__uint128_t>(value);
    ledger.add(static_cast<u64>(bits));
    ledger.add(static_cast<u64>(bits >> 64));
}

bool fraction_less_exact(i128 a, i128 b, i128 c, i128 d) {
    require(a >= 0 && b > 0 && c >= 0 && d > 0,
            "fraction comparison domain error");
    bool reversed = false;
    while (true) {
        const i128 qa = a / b;
        const i128 qc = c / d;
        if (qa != qc) return reversed ? qa > qc : qa < qc;
        a %= b;
        c %= d;
        if (a == 0 || c == 0) {
            if (a == 0 && c == 0) return false;
            const bool current_less = a == 0;
            return reversed ? !current_less : current_less;
        }
        std::swap(a, b);
        std::swap(c, d);
        reversed = !reversed;
    }
}

u64 scan_body_cover(const std::vector<u32>& deck, u64& checks,
                    u64& max_checks, u64& failures, u64& row_fnv) {
    FnvLocal cover, rows;
    checks = max_checks = failures = 0;
    u64 bodies = 0;
    u32 body = (u32{1} << 9) - 1;
    while (body < (u32{1} << 30)) {
        ++bodies;
        u64 tested = 0;
        u32 witness = 0;
        for (u32 mask : deck) {
            ++tested;
            if ((body & mask) == 0) { witness = mask; break; }
        }
        checks += tested;
        max_checks = std::max(max_checks, tested);
        failures += witness == 0;
        rows.add(body); rows.add(tested); rows.add(witness);
        if (witness) { cover.add(body); cover.add(witness); }
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    require(bodies == UINT64_C(14307150), "body enumeration changed");
    row_fnv = rows.state;
    return cover.state;
}

struct ShrinkAudit {
    u64 bodies = 0;
    u64 retained_checks = 0;
    u64 base_obligations = 0;
    u64 base_uncovered = 0;
    u64 unique_obligations = 0;
    u64 unique_uncovered = 0;
    std::vector<u64> blockers;
};

ShrinkAudit audit_extra_deletion(const std::vector<u32>& deck,
                                 std::size_t deleted1,
                                 std::size_t deleted2,
                                 const std::vector<u32>& replacements) {
    ShrinkAudit out;
    out.blockers.assign(deck.size(), 0);
    u32 body = (u32{1} << 9) - 1;
    while (body < (u32{1} << 30)) {
        ++out.bodies;
        std::size_t unique = deck.size();
        unsigned retained_witnesses = 0;
        for (std::size_t i = 0; i < deck.size(); ++i) {
            if (i == deleted1 || i == deleted2) continue;
            ++out.retained_checks;
            if ((body & deck[i]) != 0) continue;
            unique = i;
            if (++retained_witnesses == 2) break;
        }
        if (retained_witnesses < 2) {
            bool appended_cover = false;
            for (u32 mask : replacements)
                appended_cover |= (body & mask) == 0;
            if (retained_witnesses == 0) {
                ++out.base_obligations;
                out.base_uncovered += !appended_cover;
            } else {
                ++out.unique_obligations;
                if (!appended_cover) {
                    ++out.unique_uncovered;
                    ++out.blockers[unique];
                }
            }
        }
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    return out;
}

}  // namespace

#ifndef TWO_MASK_LITERAL_VERIFY_LIBRARY_ONLY
int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 10,
                "usage: verify DECK SIGS INDEX1 INDEX2 Q R REP1 REP2 REP3");
        const std::vector<u32> deck = read_deck(argv[1]);
        const std::vector<SigRow> rows = read_signatures(argv[2]);
        std::size_t deleted1 = std::stoul(argv[3]);
        std::size_t deleted2 = std::stoul(argv[4]);
        require(deleted1 < deleted2 && deleted2 < deck.size(), "bad indices");
        const int target_q = std::stoi(argv[5]);
        const int target_r = std::stoi(argv[6]);
        require(deleted1 == 107 && deleted2 == 374 &&
                    target_q == 512 && target_r == 644,
                "target surgery identity changed");
        std::vector<u32> replacements;
        for (int i = 7; i < 10; ++i)
            replacements.push_back(
                static_cast<u32>(std::stoul(argv[i], nullptr, 16)));
        std::set<u32> replacement_set;
        for (u32 mask : replacements)
            require(std::popcount(mask) == 8 && mask < (u32{1} << 30) &&
                        std::find(deck.begin(), deck.end(), mask) == deck.end() &&
                        replacement_set.insert(mask).second,
                    "invalid replacement");
        require(replacements == std::vector<u32>{UINT32_C(0x32043014),
                                                  UINT32_C(0x20807016),
                                                  UINT32_C(0x128c8012)},
                "replacement identity changed");

        const auto wanted = pair_signature(deleted1, deleted2);
        std::vector<std::pair<int, int>> fibre;
        bool target_found = false;
        for (const SigRow& row : rows) {
            if (row.words != wanted) continue;
            fibre.push_back({row.q, row.r});
            target_found |= row.q == target_q && row.r == target_r;
        }
        require(target_found && !fibre.empty(), "target missing from fibre");

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

        FnvLocal fibre_ledger, row_ledger;
        u64 equalities = 0;
        i128 minimum_margin = 0, minimum_grid = 1;
        std::pair<int, int> minimum_pair{};
        u32 minimum_mask = 0;
        bool minimum_set = false;
        std::cout << "THM4286_TWO_MASK_LITERAL_V1\n";
        for (const auto [q, r] : fibre) {
            fibre_ledger.add(q); fibre_ledger.add(r);
            const LiteralPair pair = build_literal_pair(q, r);
            std::vector<i128> prefix(walls.size());
            for (std::size_t w = 0; w < walls.size(); ++w)
                prefix[w] = prefix_at(pair,
                    static_cast<i128>(walls[w]) * pair.pool_scale);
            std::vector<std::size_t> inactive;
            for (std::size_t i = 0; i < geometry.size(); ++i) {
                const i128 margin = repair_margin(geometry[i], prefix, pair.grid);
                equalities += margin == 0;
                if (margin < 0) inactive.push_back(i);
            }
            require(inactive == std::vector<std::size_t>({deleted1, deleted2}),
                    "literal inactive signature mismatch");
            std::cout << "ROW " << q << ',' << r;
            for (std::size_t j = 0; j < replacements.size(); ++j) {
                const i128 margin = repair_margin(replacement_geometry[j],
                                                   prefix, pair.grid);
                equalities += margin == 0;
                require(margin > 0, "replacement not strictly active");
                if (!minimum_set || fraction_less_exact(
                                        margin, pair.grid,
                                        minimum_margin, minimum_grid)) {
                    minimum_set = true;
                    minimum_margin = margin;
                    minimum_grid = pair.grid;
                    minimum_pair = {q, r};
                    minimum_mask = replacements[j];
                }
                row_ledger.add(q); row_ledger.add(r);
                row_ledger.add(replacements[j]); add_i128(row_ledger, margin);
                add_i128(row_ledger, pair.grid);
                std::cout << " M" << j + 1 << ' ' << decimal(margin)
                          << '/' << decimal(pair.grid);
            }
            std::cout << '\n';
        }
        require(equalities == 0 && fibre.size() == 5 &&
                    fibre_ledger.state == UINT64_C(0xd1822d3aaeb858d1) &&
                    row_ledger.state == UINT64_C(0xbb1f2c71b7183247) &&
                    minimum_pair == std::pair<int, int>{512, 644} &&
                    minimum_mask == UINT32_C(0x32043014) &&
                    minimum_margin == static_cast<i128>(308658430702440) &&
                    minimum_grid == static_cast<i128>(13425493330529280),
                "literal fibre summary changed");

        std::vector<u32> rebuilt;
        for (std::size_t i = 0; i < deck.size(); ++i)
            if (i != deleted1 && i != deleted2) rebuilt.push_back(deck[i]);
        rebuilt.insert(rebuilt.end(), replacements.begin(), replacements.end());
        require(rebuilt.size() == 422 &&
                    std::set<u32>(rebuilt.begin(), rebuilt.end()).size() == 422,
                "rebuilt deck invalid");
        u64 checks = 0, max_checks = 0, failures = 0, body_row_fnv = 0;
        const u64 cover_fnv = scan_body_cover(rebuilt, checks, max_checks,
                                               failures, body_row_fnv);
        require(failures == 0 && checks == UINT64_C(405321845) &&
                    max_checks == 422 &&
                    mask_fnv(rebuilt) == UINT64_C(0x35d476c57331f14f) &&
                    cover_fnv == UINT64_C(0x6d3884ea6b84094c) &&
                    body_row_fnv == UINT64_C(0x7138f103d586741f),
                "rebuilt deck body scan changed");

        const ShrinkAudit shrink = audit_extra_deletion(
            deck, deleted1, deleted2, replacements);
        require(shrink.bodies == UINT64_C(14307150) &&
                    shrink.retained_checks == UINT64_C(839366930) &&
                    shrink.base_obligations == 20 &&
                    shrink.base_uncovered == 0 &&
                    shrink.unique_obligations == 3763 &&
                    shrink.unique_uncovered == 3393,
                "shrink audit base failure");
        std::vector<std::size_t> redundant;
        for (std::size_t i = 0; i < deck.size(); ++i)
            if (i != deleted1 && i != deleted2 && shrink.blockers[i] == 0)
                redundant.push_back(i);
        require(redundant == std::vector<std::size_t>{318},
                "extra-deletion identity changed");
        const std::size_t extra = redundant.front();
        std::vector<u32> shrunk;
        for (std::size_t i = 0; i < deck.size(); ++i)
            if (i != deleted1 && i != deleted2 && i != extra)
                shrunk.push_back(deck[i]);
        shrunk.insert(shrunk.end(), replacements.begin(), replacements.end());
        require(shrunk.size() == 421 &&
                    std::set<u32>(shrunk.begin(), shrunk.end()).size() == 421,
                "shrunk deck invalid");
        u64 shrink_checks = 0, shrink_max = 0, shrink_failures = 0;
        u64 shrink_row_fnv = 0;
        const u64 shrink_cover_fnv = scan_body_cover(
            shrunk, shrink_checks, shrink_max, shrink_failures, shrink_row_fnv);
        require(shrink_failures == 0, "shrunk deck body failure");

        FnvLocal blocker_ledger;
        for (std::size_t i = 0; i < shrink.blockers.size(); ++i) {
            blocker_ledger.add(i); blocker_ledger.add(shrink.blockers[i]);
        }
        require(blocker_ledger.state == UINT64_C(0x7abe3d9521893a79) &&
                    mask_fnv(shrunk) == UINT64_C(0xc9ac86709cda10df) &&
                    shrink_checks == UINT64_C(405320191) &&
                    shrink_max == 421 && shrink_failures == 0 &&
                    shrink_cover_fnv == UINT64_C(0xe9584c43ebc26a0b) &&
                    shrink_row_fnv == UINT64_C(0x2984a661ec30ab5e),
                "shrunk deck audit changed");
        std::cout << "SUMMARY TARGET " << target_q << ',' << target_r
                  << " DELETED " << deleted1 << ':' << std::hex
                  << deck[deleted1] << ' ' << std::dec << deleted2 << ':'
                  << std::hex << deck[deleted2] << " REPLACEMENTS";
        for (u32 mask : replacements) std::cout << ' ' << mask;
        std::cout << " FIBRE " << std::dec << fibre.size()
                  << " FIBRE_FNV " << std::hex << fibre_ledger.state
                  << " ROW_FNV " << row_ledger.state << std::dec
                  << " EQUALITIES " << equalities << '\n'
                  << "MIN_GAP_PAIR " << minimum_pair.first << ','
                  << minimum_pair.second << " MASK " << std::hex
                  << minimum_mask << std::dec << " NUM "
                  << decimal(minimum_margin) << " DEN "
                  << decimal(minimum_grid) << '\n'
                  << "REBUILT422_FNV " << std::hex << mask_fnv(rebuilt)
                  << std::dec << " BODY_CHECKS " << checks
                  << " MAX_CHECKS " << max_checks << " FAILURES " << failures
                  << " COVER_FNV " << std::hex << cover_fnv
                  << " BODY_ROW_FNV " << body_row_fnv << std::dec << '\n'
                  << "SHRINK_AUDIT RETAINED_CHECKS " << shrink.retained_checks
                  << " BASE_OBLIGATIONS " << shrink.base_obligations
                  << " BASE_UNCOVERED " << shrink.base_uncovered
                  << " UNIQUE_OBLIGATIONS " << shrink.unique_obligations
                  << " UNIQUE_UNCOVERED " << shrink.unique_uncovered
                  << " BLOCKER_FNV " << std::hex << blocker_ledger.state
                  << std::dec << " REDUNDANT " << redundant.size()
                  << " INDICES";
        for (std::size_t i : redundant) std::cout << ' ' << i;
        std::cout << '\n'
                  << "SHRUNK421 EXTRA_DELETED " << extra << ':' << std::hex
                  << deck[extra] << " FNV " << mask_fnv(shrunk) << std::dec
                  << " BODY_CHECKS " << shrink_checks
                  << " MAX_CHECKS " << shrink_max
                  << " FAILURES " << shrink_failures << " COVER_FNV "
                  << std::hex << shrink_cover_fnv << " BODY_ROW_FNV "
                  << shrink_row_fnv << std::dec << '\n'
                  << "VERDICT PASS DETACHED_LITERAL_FIBRE_BODY_AND_SHRINK\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "TWO_MASK_LITERAL_ERROR " << error.what() << '\n';
        return 1;
    }
}
#endif
