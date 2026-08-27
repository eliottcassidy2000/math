// Detached literal-grid scout for the global common-activity set of the
// frozen 421-mask joint deck.  Pair safety is built from the two actual
// speeds on L=lcm(D,14q,14r); no primitive endpoint cocycle or atom mass is
// used.  Fixed-pool repair components are inherited only as interval geometry.

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#pragma clang diagnostic ignored "-Wunneeded-internal-declaration"
#endif
#define CARRIER_CEGAR_LIBRARY_ONLY
#include "lrc14_three_round_learned_carrier_thm4266/carrier_cegar_descent.cpp"
#undef CARRIER_CEGAR_LIBRARY_ONLY
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#include <atomic>
#include <fstream>
#include <map>
#include <set>

namespace {

constexpr std::size_t EXPECTED_POST4271 = 174904;
constexpr u64 EXPECTED_POST4271_FNV = UINT64_C(0xb3db855040bcf19e);
constexpr std::size_t EXPECTED_POST4277 = 172322;
constexpr u64 EXPECTED_POST4277_FNV = UINT64_C(0x30b2a7e597ac548c);
constexpr std::size_t EXPECTED_DECK = 420;
constexpr u64 EXPECTED_DECK_FNV = UINT64_C(0xe72c3c8b50ec6c6e);

struct LitEdge {
    int q = 0;
    int r = 0;
    auto operator<=>(const LitEdge&) const = default;
};

void add_edge(FnvLocal& ledger, const LitEdge& edge) {
    ledger.add(edge.q);
    ledger.add(edge.r);
}

u64 edge_fnv(const std::vector<LitEdge>& edges) {
    FnvLocal ledger;
    for (const LitEdge& edge : edges) add_edge(ledger, edge);
    return ledger.state;
}

std::vector<LitEdge> read_post4277(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open post-THM-4271 residual");
    std::vector<LitEdge> post4271;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::istringstream row(line);
        LitEdge edge;
        char comma = 0;
        char trailing = 0;
        row >> edge.q >> comma >> edge.r;
        require(static_cast<bool>(row) && comma == ',' && edge.q > 0 &&
                    edge.q < edge.r && !(row >> trailing),
                "malformed residual edge");
        post4271.push_back(edge);
    }
    require(std::is_sorted(post4271.begin(), post4271.end()) &&
                std::adjacent_find(post4271.begin(), post4271.end()) ==
                    post4271.end() &&
                post4271.size() == EXPECTED_POST4271 &&
                edge_fnv(post4271) == EXPECTED_POST4271_FNV,
            "post-THM-4271 identity changed");
    std::vector<LitEdge> post4277;
    post4277.reserve(EXPECTED_POST4277);
    for (const LitEdge& edge : post4271) {
        const bool retained4276 =
            edge.r < 670 ||
            (edge.r == 670 && (edge.q == 256 || edge.q == 384));
        if (!retained4276) continue;
        const bool rectangle = 450 <= edge.q && edge.q <= 499 &&
                               600 <= edge.r && edge.r <= 650;
        if (!rectangle) post4277.push_back(edge);
    }
    require(post4277.size() == EXPECTED_POST4277 &&
                edge_fnv(post4277) == EXPECTED_POST4277_FNV,
            "post-THM-4277 identity changed");
    return post4277;
}

u64 mask_fnv(const std::vector<u32>& masks) {
    FnvLocal ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

std::vector<u32> read_deck(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open joint deck");
    std::vector<u32> deck;
    std::set<u32> seen;
    std::string word;
    while (input >> word) {
        std::size_t used = 0;
        const u64 wide = std::stoull(word, &used, 16);
        require(used == word.size() && wide < (UINT64_C(1) << 30),
                "invalid joint mask token");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && seen.insert(mask).second,
                "joint mask arity/distinctness changed");
        deck.push_back(mask);
    }
    require(deck.size() == EXPECTED_DECK &&
                mask_fnv(deck) == EXPECTED_DECK_FNV,
            "joint-deck identity changed");
    return deck;
}

struct IndexedRepair {
    u32 mask = 0;
    std::vector<std::pair<std::uint16_t, std::uint16_t>> components;
};

IndexedRepair index_repair(u32 mask, const std::vector<Cell>& cells) {
    IndexedRepair out;
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
        out.components.push_back(
            {static_cast<std::uint16_t>(begin),
             static_cast<std::uint16_t>(index)});
    }
    require(!out.components.empty(), "empty repair geometry");
    return out;
}

i128 lcm128(i128 left, i128 right) {
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

std::vector<std::pair<i128, i128>> speed_intervals(int speed, i128 grid) {
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
    out.grid = lcm128(lcm128(static_cast<i128>(COMMON),
                            static_cast<i128>(14) * q),
                      static_cast<i128>(14) * r);
    require(out.grid % COMMON == 0, "literal grid lost pool scale");
    out.pool_scale = out.grid / COMMON;
    const auto q_arcs = speed_intervals(q, out.grid);
    const auto r_arcs = speed_intervals(r, out.grid);
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
            "literal prefix outside full circle");
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

bool rational_less_nonnegative(i128 a, i128 b, i128 c, i128 d) {
    require(a >= 0 && c >= 0 && b > 0 && d > 0,
            "invalid rational gap comparison");
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

struct PairAudit {
    u64 tests = 0;
    u64 equalities = 0;
    bool common = false;
    u32 first_inactive = 0;
    i128 first_margin = 0;
    u32 weakest_mask = 0;
    i128 weakest_margin = 0;
    i128 grid = 0;
};

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 4,
                "usage: global-literal POST4271 JOINT_DECK COMMON_OUT");
        const std::vector<LitEdge> edges = read_post4277(argv[1]);
        const std::vector<u32> deck = read_deck(argv[2]);
        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool-cell count changed");
        std::vector<i64> walls(cells.size() + 1);
        walls.front() = cells.front().left;
        for (std::size_t index = 0; index < cells.size(); ++index) {
            require(cells[index].left == walls[index],
                    "pool cells not contiguous");
            walls[index + 1] = cells[index].right;
        }
        require(walls.front() == 0 && walls.back() == COMMON,
                "pool wall endpoints changed");
        std::vector<IndexedRepair> geometry;
        geometry.reserve(deck.size());
        for (u32 mask : deck) geometry.push_back(index_repair(mask, cells));

        std::vector<PairAudit> audits(edges.size());
        std::atomic<std::size_t> next{0};
        const unsigned hardware = std::thread::hardware_concurrency();
        const unsigned threads =
            std::max(1u, std::min(8u, hardware == 0 ? 1u : hardware));
        std::vector<std::thread> workers;
        for (unsigned lane = 0; lane < threads; ++lane) {
            workers.emplace_back([&]() {
                std::vector<i128> prefix(walls.size());
                for (;;) {
                    const std::size_t row = next.fetch_add(1);
                    if (row >= edges.size()) break;
                    const LiteralPair pair =
                        build_literal_pair(edges[row].q, edges[row].r);
                    for (std::size_t wall = 0; wall < walls.size(); ++wall)
                        prefix[wall] = literal_prefix(
                            pair, static_cast<i128>(walls[wall]) *
                                      pair.pool_scale);
                    PairAudit audit;
                    audit.common = true;
                    audit.grid = pair.grid;
                    bool weakest_set = false;
                    for (const IndexedRepair& repair : geometry) {
                        i128 mass = 0;
                        for (const auto [left, right] : repair.components)
                            mass += prefix[right] - prefix[left];
                        const i128 margin = static_cast<i128>(63) * mass -
                                            static_cast<i128>(4) * pair.grid;
                        ++audit.tests;
                        audit.equalities += margin == 0;
                        if (margin < 0) {
                            audit.common = false;
                            audit.first_inactive = repair.mask;
                            audit.first_margin = margin;
                            break;
                        }
                        if (!weakest_set || margin < audit.weakest_margin ||
                            (margin == audit.weakest_margin &&
                             repair.mask < audit.weakest_mask)) {
                            weakest_set = true;
                            audit.weakest_mask = repair.mask;
                            audit.weakest_margin = margin;
                        }
                    }
                    require(!audit.common ||
                                (weakest_set && audit.tests == deck.size() &&
                                 audit.weakest_margin > 0),
                            "literal common-pair audit incomplete");
                    audits[row] = audit;
                }
            });
        }
        for (auto& worker : workers) worker.join();

        std::ofstream common_out(argv[3]);
        require(static_cast<bool>(common_out), "cannot create common CSV");
        std::vector<LitEdge> common;
        u64 tests = 0;
        u64 equalities = 0;
        bool global_set = false;
        std::size_t global_row = 0;
        for (std::size_t row = 0; row < edges.size(); ++row) {
            tests += audits[row].tests;
            equalities += audits[row].equalities;
            if (!audits[row].common) continue;
            common.push_back(edges[row]);
            common_out << edges[row].q << ',' << edges[row].r << '\n';
            if (!global_set || rational_less_nonnegative(
                    audits[row].weakest_margin, audits[row].grid,
                    audits[global_row].weakest_margin,
                    audits[global_row].grid) ||
                (!rational_less_nonnegative(
                     audits[global_row].weakest_margin,
                     audits[global_row].grid,
                     audits[row].weakest_margin, audits[row].grid) &&
                 std::tie(edges[row].q, edges[row].r,
                          audits[row].weakest_mask) <
                     std::tie(edges[global_row].q, edges[global_row].r,
                              audits[global_row].weakest_mask))) {
                global_set = true;
                global_row = row;
            }
        }
        require(common_out.good() && global_set,
                "failed writing literal common CSV");

        std::cout << "JOINT421_GLOBAL_LITERAL_ACTIVITY_SCOUT_V1\n"
                  << "POST_THM4277 " << edges.size() << " FNV " << std::hex
                  << edge_fnv(edges) << std::dec << '\n'
                  << "JOINT_DECK " << deck.size() << " FNV " << std::hex
                  << mask_fnv(deck) << std::dec << '\n'
                  << "SCOUT COMMON_ACTIVE " << common.size() << " COMMON_FNV "
                  << std::hex << edge_fnv(common) << std::dec
                  << " REPAIR_TESTS " << tests << " EQUALITIES "
                  << equalities << '\n'
                  << "MIN_COMMON_GAP PAIR " << edges[global_row].q << ','
                  << edges[global_row].r << " MASK " << std::hex
                  << audits[global_row].weakest_mask << std::dec
                  << " MARGIN_NUM "
                  << decimal(audits[global_row].weakest_margin)
                  << " DEN " << decimal(audits[global_row].grid) << '\n'
                  << "VERDICT PASS FINITE_EXACT_DETACHED_LITERAL_SCOUT_ONLY\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "JOINT421_GLOBAL_LITERAL_ERROR " << error.what() << '\n';
        return 1;
    }
}
