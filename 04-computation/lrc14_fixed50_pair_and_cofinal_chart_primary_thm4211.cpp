// Primary fixed-prefix pair and cofinal-chart audit for THM-4211.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

namespace {

constexpr std::array<int, 30> P = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
constexpr i64 D = INT64_C(18241159416480);
constexpr int Q1 = 50;

std::array<std::array<u64, 10>, 31> choose_table{};

void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

i64 lcm_exact(i64 a, i64 b) {
    const i64 g = std::gcd(a, b);
    const i128 x = static_cast<i128>(a / g) * b;
    require(x <= std::numeric_limits<i64>::max(), "lcm overflow");
    return static_cast<i64>(x);
}

u32 next_combination(u32 mask) {
    const u32 low = mask & (~mask + 1u);
    const u32 ripple = mask + low;
    if (ripple == 0) return 0;
    return ripple | (((mask ^ ripple) >> 2) / low);
}

void init_choose() {
    for (int n = 0; n <= 30; ++n) {
        choose_table[n][0] = 1;
        for (int k = 1; k <= 9; ++k) {
            choose_table[n][k] = (n == 0 ? 0 : choose_table[n - 1][k] +
                                  choose_table[n - 1][k - 1]);
        }
    }
    require(choose_table[30][8] == UINT64_C(5852925), "C(30,8) changed");
    require(choose_table[30][9] == UINT64_C(14307150), "C(30,9) changed");
}

u64 colex_rank_k(u32 mask, int target_weight) {
    u64 rank = 0;
    int ordinal = 1;
    for (int bit = 0; bit < 30; ++bit) {
        if ((mask & (u32{1} << bit)) == 0) continue;
        rank += choose_table[bit][ordinal];
        ++ordinal;
    }
    require(ordinal == target_weight + 1 &&
                rank < choose_table[30][target_weight], "bad colex mask");
    return rank;
}

u64 colex_rank(u32 mask) { return colex_rank_k(mask, 8); }

bool safe_mid(int speed, i64 left, i64 right, i64 denom) {
    const i128 period = static_cast<i128>(2) * denom;
    i128 residue = static_cast<i128>(speed) * (left + right) % period;
    if (residue < 0) residue += period;
    return static_cast<i128>(7) * residue >= denom &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * denom;
}

struct Cell {
    i64 left;
    i64 right;
    u32 failed;
    bool q1_safe;
};

struct Base {
    i64 denom;
    std::vector<Cell> all_cells;
    std::vector<Cell> q1_safe_cells;
    u64 wall_cells;
};

Base build_base() {
    Base base;
    base.denom = lcm_exact(D, 14LL * Q1);
    std::vector<i64> walls = {0, base.denom};
    auto add_walls = [&](int speed) {
        require(base.denom % (14LL * speed) == 0, "base wall not integral");
        const i64 unit = base.denom / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int p : P) add_walls(p);
    add_walls(Q1);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    base.wall_cells = walls.size() - 1;
    for (std::size_t j = 0; j + 1 < walls.size(); ++j) {
        const bool q1_safe = safe_mid(Q1, walls[j], walls[j + 1], base.denom);
        u32 failed = 0;
        for (int v = 0; v < 30; ++v) {
            if (!safe_mid(P[v], walls[j], walls[j + 1], base.denom)) {
                failed |= u32{1} << v;
            }
        }
        const Cell cell{walls[j], walls[j + 1], failed, q1_safe};
        base.all_cells.push_back(cell);
        if (q1_safe) base.q1_safe_cells.push_back(cell);
    }
    require(base.denom == INT64_C(91205797082400), "D50 changed");
    require(base.wall_cells > 7133 && base.wall_cells < 7300,
            "P+50 wall-cell count out of range");
    return base;
}

// F_q(t) is 14q times the Lebesgue length of G_q intersect [0,t/denom],
// represented in denom-ticks. Thus F_q(denom)=12q*denom exactly.
i64 safe_prefix(i64 tick, int q, i64 denom) {
    const i128 product = static_cast<i128>(q) * tick;
    const i64 whole = static_cast<i64>(product / denom);
    const i64 rem = static_cast<i64>(product % denom);
    const i128 scaled = static_cast<i128>(14) * rem;
    i128 partial = 0;
    if (scaled <= denom) partial = 0;
    else if (scaled >= static_cast<i128>(13) * denom) partial =
        static_cast<i128>(12) * denom;
    else partial = scaled - denom;
    const i128 answer = static_cast<i128>(12) * whole * denom + partial;
    require(answer <= std::numeric_limits<i64>::max(), "safe-prefix overflow");
    return static_cast<i64>(answer);
}

struct AtomData {
    std::unordered_map<u32, i64> mass;
    std::array<u64, 9> unique_by_weight{};
    std::array<u64, 9> cells_by_weight{};
};

AtomData build_atoms(const Base& base, int q2) {
    require(q2 > 0 && q2 != Q1, "invalid second newcomer");
    require(std::find(P.begin(), P.end(), q2) == P.end(), "q2 lies in pool");
    AtomData result;
    result.mass.reserve(8192);
    for (const Cell& cell : base.q1_safe_cells) {
        const i64 contribution = safe_prefix(cell.right, q2, base.denom) -
                                 safe_prefix(cell.left, q2, base.denom);
        require(contribution >= 0, "negative q2 contribution");
        const int w = std::popcount(cell.failed);
        if (w <= 8 && contribution != 0) {
            result.mass[cell.failed] += contribution;
            ++result.cells_by_weight[w];
        }
    }
    for (const auto& [mask, mass] : result.mass) {
        require(mass > 0, "nonpositive atom mass");
        ++result.unique_by_weight[std::popcount(mask)];
    }
    return result;
}

void add_supersets_rec(u32 atom, int need, int start, u32 extra, i64 value,
                       std::vector<i64>& zeta, u64& operations) {
    if (need == 0) {
        zeta[colex_rank(atom | extra)] += value;
        ++operations;
        return;
    }
    for (int bit = start; bit <= 30 - need; ++bit) {
        const u32 b = u32{1} << bit;
        if ((atom & b) != 0) continue;
        add_supersets_rec(atom, need - 1, bit + 1, extra | b, value,
                          zeta, operations);
    }
}

struct Layer {
    std::vector<u32> edges;
    u64 equalities = 0;
    u64 zeta_operations = 0;
};

Layer build_e8(const Base& base, const AtomData& atoms, int q2) {
    std::vector<i64> masses(choose_table[30][8], 0);
    u64 operations = 0;
    for (const auto& [atom, mass] : atoms.mass) {
        const int w = std::popcount(atom);
        add_supersets_rec(atom, 8 - w, 0, 0, mass, masses, operations);
    }
    Layer layer;
    layer.zeta_operations = operations;
    u32 deletion = (u32{1} << 8) - 1;
    const u32 limit = u32{1} << 30;
    u64 index = 0;
    while (deletion != 0 && deletion < limit) {
        require(colex_rank(deletion) == index, "colex enumeration mismatch");
        const i128 delta = static_cast<i128>(9) * masses[index] -
                           static_cast<i128>(8) * q2 * base.denom;
        if (delta == 0) ++layer.equalities;
        if (delta >= 0) layer.edges.push_back(deletion);
        ++index;
        const u32 next = next_combination(deletion);
        if (next <= deletion) break;
        deletion = next;
    }
    require(index == choose_table[30][8], "rank-8 universe incomplete");
    return layer;
}

u64 mix64(u64 value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

struct Scan {
    u64 bodies = 0;
    u64 checks = 0;
    u64 max_checks = 0;
    u32 closest = 0;
    u32 missed = 0;
    u32 first_cover = 0;
    u64 covers = 0;
};

Scan scan_nine(std::vector<u32> edges, bool stop_at_first_cover) {
    std::sort(edges.begin(), edges.end(), [](u32 a, u32 b) {
        const u64 x = mix64(a), y = mix64(b);
        return x == y ? a < b : x < y;
    });
    Scan scan;
    u32 body = (u32{1} << 9) - 1;
    const u32 limit = u32{1} << 30;
    while (body < limit) {
        u32 missed = 0;
        u64 checked = 0;
        for (u32 edge : edges) {
            ++checked;
            if ((body & edge) == 0) {
                missed = edge;
                break;
            }
        }
        ++scan.bodies;
        scan.checks += checked;
        if (checked > scan.max_checks) {
            scan.max_checks = checked;
            scan.closest = body;
            scan.missed = missed;
        }
        if (missed == 0) {
            if (scan.first_cover == 0) scan.first_cover = body;
            ++scan.covers;
            if (stop_at_first_cover) break;
        }
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    return scan;
}

u32 mask_of(std::initializer_list<int> values) {
    u32 mask = 0;
    for (int x : values) {
        auto it = std::find(P.begin(), P.end(), x);
        require(it != P.end(), "label absent from P");
        mask |= u32{1} << std::distance(P.begin(), it);
    }
    return mask;
}

std::string labels(u32 mask) {
    std::string out;
    for (int v = 0; v < 30; ++v) {
        if ((mask & (u32{1} << v)) == 0) continue;
        if (!out.empty()) out += ',';
        out += std::to_string(P[v]);
    }
    return out;
}

std::string decimal(i128 x) {
    if (x == 0) return "0";
    const bool neg = x < 0;
    if (neg) x = -x;
    std::string out;
    while (x != 0) {
        out.push_back(static_cast<char>('0' + x % 10));
        x /= 10;
    }
    if (neg) out.push_back('-');
    std::reverse(out.begin(), out.end());
    return out;
}

i128 body_delta(const Base& base, int q2, u32 body) {
    i128 mass = 0;
    for (const Cell& cell : base.q1_safe_cells) {
        if ((cell.failed & body) != 0) continue;
        mass += safe_prefix(cell.right, q2, base.denom) -
                safe_prefix(cell.left, q2, base.denom);
    }
    return static_cast<i128>(63) * mass -
           static_cast<i128>(56) * q2 * base.denom;
}

// Frozen theorem output is semantic only; wall-clock timing is deliberately
// excluded so independent optimization modes can byte-match.

void run_one(const Base& base, int q2, bool all_covers) {
    const AtomData atoms = build_atoms(base, q2);
    const Layer layer = build_e8(base, atoms, q2);
    const Scan scan = scan_nine(layer.edges, !all_covers);
    std::cout << "QPAIR 50," << q2
              << " BASE_DEN " << base.denom
              << " BASE_CELLS " << base.wall_cells
              << " Q1_SAFE_CELLS " << base.q1_safe_cells.size()
              << " ATOMS " << atoms.mass.size()
              << " ZETA_OPS " << layer.zeta_operations
              << " E8 " << layer.edges.size()
              << " EQ " << layer.equalities
              << " BODIES_SCANNED " << scan.bodies
              << " CHECKS " << scan.checks
              << " MAX " << scan.max_checks
              << " FIRST_COVER " << labels(scan.first_cover)
              << " COVERS " << scan.covers
              << " CLOSEST " << labels(scan.closest)
              << " MISSED " << labels(scan.missed)
              << '\n';
    if (q2 == 49) {
        require(layer.edges.size() == UINT64_C(1536023), "49 E8 edge mismatch");
        require(layer.equalities == 0, "49 equality mismatch");
        if (all_covers) {
            require(scan.bodies == UINT64_C(14307150) && scan.covers == 0,
                    "49 no-cover mismatch");
        }
    }
    if (q2 == 6) {
        require(layer.edges.size() == UINT64_C(13497), "6 E8 edge mismatch");
        require(layer.equalities == 0, "6 equality mismatch");
        const u32 hostile = mask_of({8,63,80,84,85,88,120,143,145});
        require(std::all_of(layer.edges.begin(), layer.edges.end(),
                            [&](u32 e) { return (e & hostile) != 0; }),
                "explicit q6 hostile ceased to cover E8");
        const i128 scaled = body_delta(base, q2, hostile);
        // The joint-wall scout uses denominator 5D; this prefix route uses
        // 14*q2*(5D), so its printed numerator is larger by exactly 84.
        require(scaled == static_cast<i128>(601044495065784LL) * 84,
                "q6 hostile body surplus mismatch");
        std::cout << "Q6_HOSTILE " << labels(hostile)
                  << " PREFIX_SCALED_BODY_DELTA " << decimal(scaled)
                  << " JOINT_WALL_DELTA " << decimal(scaled / 84) << '\n';
    }
}

struct LimitAtom {
    i64 length = 0;
    i64 safe_cells = 0;
    i64 safe_adjacencies = 0;
};

void add_limit_supersets_rec(u32 atom, int need, int start, u32 extra,
                             const LimitAtom& value,
                             std::vector<i64>& lengths,
                             std::vector<i64>& cell_counts,
                             std::vector<i64>& adjacency_counts,
                             u64& operations) {
    if (need == 0) {
        const u64 rank = colex_rank(atom | extra);
        lengths[rank] += value.length;
        cell_counts[rank] += value.safe_cells;
        adjacency_counts[rank] += value.safe_adjacencies;
        ++operations;
        return;
    }
    for (int bit = start; bit <= 30 - need; ++bit) {
        const u32 b = u32{1} << bit;
        if ((atom & b) != 0) continue;
        add_limit_supersets_rec(atom, need - 1, bit + 1, extra | b, value,
                                lengths, cell_counts, adjacency_counts, operations);
    }
}

struct ChartEdge {
    u32 mask = 0;
    u64 cutoff = 0;
    i64 strict_surplus = 0;
    i64 components = 0;
};

void run_limit_chart(const Base& base) {
    std::unordered_map<u32, LimitAtom> atoms;
    atoms.reserve(8192);
    for (const Cell& cell : base.all_cells) {
        if (!cell.q1_safe || std::popcount(cell.failed) > 8) continue;
        atoms[cell.failed].length += cell.right - cell.left;
        ++atoms[cell.failed].safe_cells;
    }
    for (std::size_t j = 0; j < base.all_cells.size(); ++j) {
        const Cell& previous = base.all_cells[(j + base.all_cells.size() - 1) %
                                              base.all_cells.size()];
        const Cell& current = base.all_cells[j];
        if (!previous.q1_safe || !current.q1_safe) continue;
        const u32 joined = previous.failed | current.failed;
        if (std::popcount(joined) <= 8) ++atoms[joined].safe_adjacencies;
    }
    std::vector<i64> lengths(choose_table[30][8], 0);
    std::vector<i64> cell_counts(choose_table[30][8], 0);
    std::vector<i64> adjacency_counts(choose_table[30][8], 0);
    u64 zeta_operations = 0;
    for (const auto& [atom, value] : atoms) {
        add_limit_supersets_rec(atom, 8 - std::popcount(atom), 0, 0, value,
                                lengths, cell_counts, adjacency_counts,
                                zeta_operations);
    }
    std::vector<ChartEdge> edges;
    edges.reserve(5000000);
    u64 equalities = 0;
    i64 minimum_surplus = 0;
    u32 minimum_surplus_edge = 0;
    u32 deletion = (u32{1} << 8) - 1;
    const u32 limit = u32{1} << 30;
    u64 rank = 0;
    while (deletion != 0 && deletion < limit) {
        const i128 surplus128 = static_cast<i128>(27) * lengths[rank] -
                                static_cast<i128>(2) * base.denom;
        if (surplus128 == 0) ++equalities;
        if (surplus128 > 0) {
            require(surplus128 <= std::numeric_limits<i64>::max(),
                    "limit surplus overflow");
            const i64 surplus = static_cast<i64>(surplus128);
            const i64 components = cell_counts[rank] - adjacency_counts[rank];
            require(components > 0, "nonpositive component count");
            const i128 numerator = static_cast<i128>(27) * components * base.denom;
            const i128 denominator = static_cast<i128>(7) * surplus;
            const i128 cutoff128 = (numerator + denominator - 1) / denominator;
            require(cutoff128 > 0 && cutoff128 <= std::numeric_limits<u64>::max(),
                    "chart cutoff overflow");
            edges.push_back({deletion, static_cast<u64>(cutoff128), surplus,
                             components});
            if (minimum_surplus_edge == 0 || surplus < minimum_surplus) {
                minimum_surplus = surplus;
                minimum_surplus_edge = deletion;
            }
        }
        ++rank;
        const u32 next = next_combination(deletion);
        if (next <= deletion) break;
        deletion = next;
    }
    require(rank == choose_table[30][8], "limit rank-8 universe incomplete");
    std::sort(edges.begin(), edges.end(), [](const ChartEdge& a, const ChartEdge& b) {
        if (a.cutoff != b.cutoff) return a.cutoff < b.cutoff;
        const u64 ah = mix64(a.mask), bh = mix64(b.mask);
        return ah == bh ? a.mask < b.mask : ah < bh;
    });
    u64 bodies = 0;
    u64 checks = 0;
    u64 maximum_checks = 0;
    u64 chart_number = 0;
    u32 chart_body = 0;
    ChartEdge chart_edge{};
    u32 limit_cover = 0;
    i64 limit_cover_max_surplus = std::numeric_limits<i64>::min();
    u32 limit_cover_best_deletion = 0;
    i64 limit_cover_body_mass = 0;
    u32 body = (u32{1} << 9) - 1;
    while (body < limit) {
        u64 checked = 0;
        const ChartEdge* first = nullptr;
        for (const ChartEdge& edge : edges) {
            ++checked;
            if ((body & edge.mask) == 0) {
                first = &edge;
                break;
            }
        }
        if (first == nullptr) {
            limit_cover = body;
            u32 candidate = (u32{1} << 8) - 1;
            while (candidate != 0 && candidate < limit) {
                if ((candidate & body) == 0) {
                    const i64 surplus = static_cast<i64>(
                        static_cast<i128>(27) * lengths[colex_rank(candidate)] -
                        static_cast<i128>(2) * base.denom);
                    if (limit_cover_best_deletion == 0 ||
                        surplus > limit_cover_max_surplus) {
                        limit_cover_max_surplus = surplus;
                        limit_cover_best_deletion = candidate;
                    }
                }
                const u32 next_candidate = next_combination(candidate);
                if (next_candidate <= candidate) break;
                candidate = next_candidate;
            }
            for (const Cell& cell : base.q1_safe_cells) {
                if ((cell.failed & body) == 0) {
                    limit_cover_body_mass += cell.right - cell.left;
                }
            }
            break;
        }
        ++bodies;
        checks += checked;
        maximum_checks = std::max(maximum_checks, checked);
        if (first->cutoff > chart_number) {
            chart_number = first->cutoff;
            chart_body = body;
            chart_edge = *first;
        }
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    if (limit_cover != 0) {
        const i128 body_limit_delta =
            static_cast<i128>(54) * limit_cover_body_mass -
            static_cast<i128>(4) * base.denom;
        std::cout << "FIXED50_LIMIT_OBSTRUCTION BASE_DEN " << base.denom
                  << " BASE_CELLS " << base.all_cells.size()
                  << " ATOMS " << atoms.size()
                  << " ZETA_OPS " << zeta_operations
                  << " STRICT_E8 " << edges.size()
                  << " EQUALITIES " << equalities
                  << " MIN_SURPLUS " << minimum_surplus
                  << " MIN_SURPLUS_EDGE " << labels(minimum_surplus_edge)
                  << " BODIES_BEFORE_FIRST_COVER " << bodies
                  << " FIRST_COVER " << labels(limit_cover)
                  << " BEST_DISJOINT_DELETION " << labels(limit_cover_best_deletion)
                  << " BEST_DISJOINT_SURPLUS " << limit_cover_max_surplus
                  << " COVER_BODY_BASE_MASS " << limit_cover_body_mass
                  << " COVER_BODY_LIMIT_DELTA " << decimal(body_limit_delta) << '\n';
        return;
    }
    require(bodies == UINT64_C(14307150), "chart nine-body universe incomplete");
    const u64 below_count = std::count_if(edges.begin(), edges.end(),
        [&](const ChartEdge& edge) { return edge.cutoff < chart_number; });
    const u64 at_count = std::count_if(edges.begin(), edges.end(),
        [&](const ChartEdge& edge) { return edge.cutoff <= chart_number; });
    require(std::all_of(edges.begin(), edges.begin() + below_count,
                        [&](const ChartEdge& edge) {
                            return (edge.mask & chart_body) != 0;
                        }), "chart witness does not cover below-chart deck");
    std::cout << "FIXED50_LIMIT_CHART BASE_DEN " << base.denom
              << " BASE_CELLS " << base.all_cells.size()
              << " ATOMS " << atoms.size()
              << " ZETA_OPS " << zeta_operations
              << " STRICT_E8 " << edges.size()
              << " EQUALITIES " << equalities
              << " MIN_SURPLUS " << minimum_surplus
              << " MIN_SURPLUS_EDGE " << labels(minimum_surplus_edge)
              << " BODIES " << bodies
              << " CHECKS " << checks
              << " MAX_CHECKS " << maximum_checks
              << " CHART_NUMBER " << chart_number
              << " BELOW_EDGES " << below_count
              << " AT_EDGES " << at_count
              << " WITNESS_BODY " << labels(chart_body)
              << " FIRST_DISJOINT_AT_CHART " << labels(chart_edge.mask)
              << " EDGE_SURPLUS " << chart_edge.strict_surplus
              << " EDGE_COMPONENTS " << chart_edge.components << '\n';
}

void add_e9_supersets_rec(u32 atom, int need, int start, u32 extra,
                          const LimitAtom& value,
                          std::vector<i64>& lengths,
                          std::vector<std::uint16_t>& cell_counts,
                          std::vector<std::uint16_t>& adjacency_counts,
                          u64& operations) {
    if (need == 0) {
        const u64 rank = colex_rank_k(atom | extra, 9);
        lengths[rank] += value.length;
        cell_counts[rank] = static_cast<std::uint16_t>(
            cell_counts[rank] + value.safe_cells);
        adjacency_counts[rank] = static_cast<std::uint16_t>(
            adjacency_counts[rank] + value.safe_adjacencies);
        ++operations;
        return;
    }
    for (int bit = start; bit <= 30 - need; ++bit) {
        const u32 b = u32{1} << bit;
        if ((atom & b) != 0) continue;
        add_e9_supersets_rec(atom, need - 1, bit + 1, extra | b, value,
                             lengths, cell_counts, adjacency_counts, operations);
    }
}

struct ChartEdge9 {
    u64 cutoff = 0;
    u32 mask = 0;
    i64 surplus = 0;
    std::uint16_t components = 0;
};

void run_limit_chart_e9(const Base& base) {
    std::unordered_map<u32, LimitAtom> atoms;
    atoms.reserve(8192);
    std::array<u64, 10> support_by_weight{};
    for (const Cell& cell : base.all_cells) {
        if (!cell.q1_safe || std::popcount(cell.failed) > 9) continue;
        atoms[cell.failed].length += cell.right - cell.left;
        ++atoms[cell.failed].safe_cells;
    }
    for (std::size_t j = 0; j < base.all_cells.size(); ++j) {
        const Cell& previous = base.all_cells[(j + base.all_cells.size() - 1) %
                                              base.all_cells.size()];
        const Cell& current = base.all_cells[j];
        if (!previous.q1_safe || !current.q1_safe) continue;
        const u32 joined = previous.failed | current.failed;
        if (std::popcount(joined) <= 9) ++atoms[joined].safe_adjacencies;
    }
    for (const auto& [mask, value] : atoms) {
        (void)value;
        ++support_by_weight[std::popcount(mask)];
    }
    std::cout << "E9_SUPPORT ATOMS " << atoms.size() << " BY_WEIGHT";
    for (int w = 0; w <= 9; ++w) std::cout << ' ' << w << ':' << support_by_weight[w];
    std::cout << std::endl;

    const u64 universe = choose_table[30][9];
    std::vector<i64> lengths(universe, 0);
    std::vector<std::uint16_t> cell_counts(universe, 0);
    std::vector<std::uint16_t> adjacency_counts(universe, 0);
    u64 zeta_operations = 0;
    for (const auto& [atom, value] : atoms) {
        add_e9_supersets_rec(atom, 9 - std::popcount(atom), 0, 0, value,
                             lengths, cell_counts, adjacency_counts,
                             zeta_operations);
    }
    std::cout << "E9_ZETA_DONE OPS " << zeta_operations << std::endl;

    std::vector<ChartEdge9> edges;
    edges.reserve(10000000);
    u64 equalities = 0;
    i64 minimum_surplus = 0;
    u32 minimum_surplus_edge = 0;
    u32 deletion = (u32{1} << 9) - 1;
    const u32 limit = u32{1} << 30;
    u64 rank = 0;
    while (deletion != 0 && deletion < limit) {
        const i128 surplus128 = static_cast<i128>(27) * lengths[rank] -
                                static_cast<i128>(2) * base.denom;
        if (surplus128 == 0) ++equalities;
        if (surplus128 > 0) {
            require(surplus128 <= std::numeric_limits<i64>::max(),
                    "E9 surplus overflow");
            const i64 surplus = static_cast<i64>(surplus128);
            const i64 components64 = static_cast<i64>(cell_counts[rank]) -
                                     static_cast<i64>(adjacency_counts[rank]);
            require(components64 > 0 && components64 <= 65535,
                    "E9 component count invalid");
            const i128 numerator = static_cast<i128>(27) * components64 * base.denom;
            const i128 denominator = static_cast<i128>(7) * surplus;
            const i128 cutoff128 = (numerator + denominator - 1) / denominator;
            require(cutoff128 > 0 && cutoff128 <= std::numeric_limits<u64>::max(),
                    "E9 cutoff overflow");
            edges.push_back({static_cast<u64>(cutoff128), deletion, surplus,
                             static_cast<std::uint16_t>(components64)});
            if (minimum_surplus_edge == 0 || surplus < minimum_surplus) {
                minimum_surplus = surplus;
                minimum_surplus_edge = deletion;
            }
        }
        ++rank;
        const u32 next = next_combination(deletion);
        if (next <= deletion) break;
        deletion = next;
    }
    require(rank == universe, "E9 deletion universe incomplete");
    std::sort(edges.begin(), edges.end(), [](const ChartEdge9& a, const ChartEdge9& b) {
        if (a.cutoff != b.cutoff) return a.cutoff < b.cutoff;
        const u64 ah = mix64(a.mask), bh = mix64(b.mask);
        return ah == bh ? a.mask < b.mask : ah < bh;
    });
    require(!edges.empty(), "E9 strict deck empty");
    std::cout << "E9_EDGES_SORTED COUNT " << edges.size()
              << " EQUALITIES " << equalities
              << " MIN_ACTIVATION " << edges.front().cutoff
              << " MAX_ACTIVATION " << edges.back().cutoff << std::endl;

    struct ThresholdResult {
        bool no_cover = false;
        u32 first_cover = 0;
        u64 active_edges = 0;
        u64 bodies = 0;
        u64 checks = 0;
        u64 maximum_checks = 0;
    };
    auto threshold_test = [&](u64 cutoff) {
        ThresholdResult result;
        std::vector<u32> active;
        const auto end = std::upper_bound(edges.begin(), edges.end(), cutoff,
            [](u64 value, const ChartEdge9& edge) { return value < edge.cutoff; });
        active.reserve(std::distance(edges.begin(), end));
        for (auto it = edges.begin(); it != end; ++it) active.push_back(it->mask);
        std::sort(active.begin(), active.end(), [](u32 a, u32 b) {
            const u64 ah = mix64(a), bh = mix64(b);
            return ah == bh ? a < b : ah < bh;
        });
        result.active_edges = active.size();
        u32 body = (u32{1} << 9) - 1;
        while (body < limit) {
            u64 checked = 0;
            u32 missed = 0;
            for (u32 edge : active) {
                ++checked;
                if ((body & edge) == 0) {
                    missed = edge;
                    break;
                }
            }
            ++result.bodies;
            result.checks += checked;
            result.maximum_checks = std::max(result.maximum_checks, checked);
            if (missed == 0) {
                result.first_cover = body;
                return result;
            }
            const u32 next = next_combination(body);
            if (next <= body) break;
            body = next;
        }
        require(result.bodies == universe, "threshold body universe incomplete");
        result.no_cover = true;
        return result;
    };

    std::vector<u64> activation_values;
    activation_values.reserve(edges.size());
    for (const ChartEdge9& edge : edges) {
        if (activation_values.empty() || activation_values.back() != edge.cutoff) {
            activation_values.push_back(edge.cutoff);
        }
    }
    ThresholdResult high_result = threshold_test(activation_values.back());
    if (!high_result.no_cover) {
        std::cout << "FIXED50_E9_LIMIT_COVER COVER "
                  << labels(high_result.first_cover)
                  << " STRICT_EDGES " << edges.size() << '\n';
        return;
    }
    std::size_t low_index = 0;
    std::size_t high_index = activation_values.size() - 1;
    ThresholdResult low_result = threshold_test(activation_values[low_index]);
    if (low_result.no_cover) {
        high_index = low_index;
        high_result = low_result;
    } else {
        while (high_index - low_index > 1) {
            const std::size_t middle = low_index + (high_index - low_index) / 2;
            ThresholdResult middle_result = threshold_test(activation_values[middle]);
            std::cout << "E9_THRESHOLD_PROBE Q " << activation_values[middle]
                      << " ACTIVE " << middle_result.active_edges
                      << " NO_COVER " << middle_result.no_cover
                      << " FIRST_COVER " << labels(middle_result.first_cover)
                      << " BODIES " << middle_result.bodies
                      << " CHECKS " << middle_result.checks << std::endl;
            if (middle_result.no_cover) {
                high_index = middle;
                high_result = middle_result;
            } else {
                low_index = middle;
                low_result = middle_result;
            }
        }
    }
    const u64 chart_number = activation_values[high_index];
    const u64 previous_activation =
        (high_index == 0 ? chart_number - 1 : activation_values[high_index - 1]);
    if (high_index != low_index || low_result.no_cover) {
        low_result = threshold_test(chart_number - 1);
    }
    require(high_result.no_cover && !low_result.no_cover &&
                low_result.first_cover != 0,
            "E9 chart bracketing failed");
    const u32 chart_body = low_result.first_cover;
    const auto below = std::lower_bound(edges.begin(), edges.end(), chart_number,
        [](const ChartEdge9& edge, u64 cutoff) { return edge.cutoff < cutoff; });
    const auto at = std::upper_bound(edges.begin(), edges.end(), chart_number,
        [](u64 cutoff, const ChartEdge9& edge) { return cutoff < edge.cutoff; });
    require(std::all_of(edges.begin(), below,
                        [&](const ChartEdge9& edge) {
                            return (edge.mask & chart_body) != 0;
                        }), "E9 sharp chart witness invalid");
    auto chart_it = std::find_if(below, at, [&](const ChartEdge9& edge) {
        return (edge.mask & chart_body) == 0;
    });
    require(chart_it != at, "E9 chart witness has no activating edge");
    const ChartEdge9 chart_edge = *chart_it;
    std::cout << "FIXED50_E9_LIMIT_CHART BASE_DEN " << base.denom
              << " BASE_CELLS " << base.all_cells.size()
              << " STRICT_E9 " << edges.size()
              << " EQUALITIES " << equalities
              << " MIN_SURPLUS " << minimum_surplus
              << " MIN_SURPLUS_EDGE " << labels(minimum_surplus_edge)
              << " MIN_ACTIVATION " << edges.front().cutoff
              << " MAX_ACTIVATION " << edges.back().cutoff
              << " PREVIOUS_ACTIVATION " << previous_activation
              << " BODIES " << high_result.bodies
              << " CHECKS " << high_result.checks
              << " MAX_CHECKS " << high_result.maximum_checks
              << " CHART_NUMBER " << chart_number
              << " BELOW_EDGES " << std::distance(edges.begin(), below)
              << " AT_EDGES " << std::distance(edges.begin(), at)
              << " WITNESS_BODY " << labels(chart_body)
              << " FIRST_DISJOINT_AT_CHART " << labels(chart_edge.mask)
              << " EDGE_SURPLUS " << chart_edge.surplus
              << " EDGE_COMPONENTS " << chart_edge.components << '\n';
}

} // namespace

int main(int argc, char** argv) {
    init_choose();
    const Base base = build_base();
    require(argc >= 2, "usage: fixed50 q2 [q2 ...] [--all-covers]");
    bool all_covers = false;
    for (int j = 1; j < argc; ++j) {
        if (std::string(argv[j]) == "--all-covers") all_covers = true;
    }
    for (int j = 1; j < argc; ++j) {
        if (std::string(argv[j]) == "--all-covers") continue;
        if (std::string(argv[j]) == "limit-chart") {
            run_limit_chart(base);
            continue;
        }
        if (std::string(argv[j]) == "limit-chart-e9") {
            run_limit_chart_e9(base);
            continue;
        }
        run_one(base, std::stoi(argv[j]), all_covers);
    }
}
