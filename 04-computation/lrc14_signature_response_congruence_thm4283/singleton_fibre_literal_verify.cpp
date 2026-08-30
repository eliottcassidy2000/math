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

std::array<u64, 7> singleton_signature(std::size_t index) {
    std::array<u64, 7> out{};
    out[index / 64] = UINT64_C(1) << (index % 64);
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

i128 lcm128(i128 a, i128 b) {
    return a / gcd_i128(a, b) * b;
}

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
    for (int tooth = 0; tooth < speed; ++tooth) {
        out.push_back({static_cast<i128>(14 * tooth + 1) * unit,
                       static_cast<i128>(14 * tooth + 13) * unit});
    }
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
    std::size_t i = 0;
    std::size_t j = 0;
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
    std::size_t lo = 0;
    std::size_t hi = pair.arcs.size();
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

u64 scan_body_cover(const std::vector<u32>& deck, u64& checks,
                    u64& max_checks, u64& failures, u64& row_fnv) {
    FnvLocal cover;
    FnvLocal rows;
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

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 7,
                "usage: literal-verify DECK SIGNATURES INDEX Q R REPLACEMENT");
        const std::vector<u32> deck = read_deck(argv[1]);
        const std::vector<SigRow> rows = read_signatures(argv[2]);
        const std::size_t index = std::stoul(argv[3]);
        const int target_q = std::stoi(argv[4]);
        const int target_r = std::stoi(argv[5]);
        const u32 replacement = static_cast<u32>(std::stoul(argv[6], nullptr, 16));
        require(index < deck.size() && std::popcount(replacement) == 8 &&
                    replacement < (u32{1} << 30) &&
                    std::find(deck.begin(), deck.end(), replacement) == deck.end(),
                "invalid replacement");
        const auto wanted = singleton_signature(index);
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
        const IndexedRepair replacement_geometry = index_repair(replacement, cells);

        FnvLocal fibre_ledger;
        FnvLocal row_ledger;
        u64 equalities = 0;
        i128 minimum_margin = 0;
        i128 minimum_grid = 1;
        std::pair<int, int> minimum_pair{};
        bool minimum_set = false;
        std::cout << "THM4283_SINGLETON_FIBRE_LITERAL_V1\n";
        for (const auto [q, r] : fibre) {
            fibre_ledger.add(q); fibre_ledger.add(r);
            const LiteralPair pair = build_literal_pair(q, r);
            std::vector<i128> prefix(walls.size());
            for (std::size_t w = 0; w < walls.size(); ++w)
                prefix[w] = prefix_at(pair,
                    static_cast<i128>(walls[w]) * pair.pool_scale);
            u64 inactive = 0;
            std::size_t inactive_index = deck.size();
            i128 old_margin = 0;
            for (std::size_t i = 0; i < geometry.size(); ++i) {
                const i128 m = repair_margin(geometry[i], prefix, pair.grid);
                equalities += m == 0;
                if (m < 0) { ++inactive; inactive_index = i; }
                if (i == index) old_margin = m;
            }
            const i128 new_margin = repair_margin(replacement_geometry, prefix,
                                                   pair.grid);
            equalities += new_margin == 0;
            require(inactive == 1 && inactive_index == index &&
                        old_margin < 0 && new_margin > 0,
                    "literal fibre/replacement check failed");
            if (!minimum_set || new_margin * minimum_grid <
                                minimum_margin * pair.grid) {
                minimum_set = true;
                minimum_margin = new_margin;
                minimum_grid = pair.grid;
                minimum_pair = {q, r};
            }
            row_ledger.add(q); row_ledger.add(r); row_ledger.add(index);
            row_ledger.add(deck[index]); add_i128(row_ledger, old_margin);
            row_ledger.add(replacement); add_i128(row_ledger, new_margin);
            add_i128(row_ledger, pair.grid);
            std::cout << "ROW " << q << ',' << r << " OLD_MARGIN "
                      << decimal(old_margin) << " NEW_MARGIN "
                      << decimal(new_margin) << " DEN " << decimal(pair.grid)
                      << '\n';
        }
        require(equalities == 0, "unexpected equality");

        std::vector<u32> rebuilt;
        rebuilt.reserve(421);
        for (std::size_t i = 0; i < deck.size(); ++i)
            if (i != index) rebuilt.push_back(deck[i]);
        rebuilt.push_back(replacement);
        require(rebuilt.size() == 421 &&
                    std::set<u32>(rebuilt.begin(), rebuilt.end()).size() == 421,
                "rebuilt deck invalid");
        u64 checks = 0, max_checks = 0, failures = 0, body_row_fnv = 0;
        const u64 cover_fnv = scan_body_cover(rebuilt, checks, max_checks,
                                               failures, body_row_fnv);
        require(failures == 0, "rebuilt deck body failure");
        std::cout << "SUMMARY TARGET " << target_q << ',' << target_r
                  << " INDEX " << index << " OLD_MASK " << std::hex
                  << deck[index] << " REPLACEMENT " << replacement << std::dec
                  << " FIBRE " << fibre.size() << " FIBRE_FNV " << std::hex
                  << fibre_ledger.state << " ROW_FNV " << row_ledger.state
                  << std::dec << " EQUALITIES " << equalities << '\n'
                  << "MIN_GAP_PAIR " << minimum_pair.first << ','
                  << minimum_pair.second << " NUM " << decimal(minimum_margin)
                  << " DEN " << decimal(minimum_grid) << '\n'
                  << "REBUILT_DECK_FNV " << std::hex << mask_fnv(rebuilt)
                  << std::dec << " BODY_CHECKS " << checks
                  << " MAX_CHECKS " << max_checks << " FAILURES " << failures
                  << " COVER_FNV " << std::hex << cover_fnv
                  << " BODY_ROW_FNV " << body_row_fnv << std::dec << '\n'
                  << "VERDICT PASS DETACHED_LITERAL_FIBRE_AND_BODY_COVER\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "LITERAL_VERIFY_ERROR " << error.what() << '\n';
        return 1;
    }
}
