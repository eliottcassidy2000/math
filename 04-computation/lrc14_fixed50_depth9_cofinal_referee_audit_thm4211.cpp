#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Independent literal-wall referee for the fixed-q=50 E9 reserve chart.
// This source does not include or depend on the primary implementation.

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

namespace {

constexpr std::array<int, 30> P = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
constexpr int Q0 = 50;
constexpr int TEST_Q_MINUS = 5681;
constexpr int TEST_Q = 5682;

std::array<std::array<u32, 11>, 31> CHOOSE{};
std::array<u32, 1u << 15> LOW_RANK{};
std::array<std::array<u32, 1u << 15>, 11> HIGH_PART{};

void need(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

i64 lcm_exact(i64 a, i64 b) {
    const i64 g = std::gcd(a, b);
    const i128 v = static_cast<i128>(a / g) * b;
    need(v <= std::numeric_limits<i64>::max(), "lcm overflow");
    return static_cast<i64>(v);
}

std::string dec(i128 x) {
    if (x == 0) return "0";
    bool neg = x < 0;
    if (neg) x = -x;
    std::string s;
    while (x) {
        s.push_back(static_cast<char>('0' + x % 10));
        x /= 10;
    }
    if (neg) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}

std::string labels(u32 mask) {
    std::ostringstream out;
    bool first = true;
    for (int i = 0; i < 30; ++i) {
        if (!(mask & (u32{1} << i))) continue;
        if (!first) out << ',';
        first = false;
        out << P[i];
    }
    return out.str();
}

u32 mask_of(std::initializer_list<int> values) {
    u32 mask = 0;
    for (int value : values) {
        auto it = std::find(P.begin(), P.end(), value);
        need(it != P.end(), "label outside P");
        mask |= u32{1} << static_cast<unsigned>(it - P.begin());
    }
    return mask;
}

u32 next_combination(u32 mask) {
    const u32 low = mask & (~mask + 1u);
    const u32 ripple = mask + low;
    if (ripple == 0) return 0;
    return ripple | (((mask ^ ripple) >> 2) / low);
}

void init_choose_and_rank() {
    CHOOSE[0][0] = 1;
    for (int n = 1; n <= 30; ++n) {
        CHOOSE[n][0] = 1;
        for (int k = 1; k <= 10; ++k) {
            CHOOSE[n][k] = CHOOSE[n - 1][k] + CHOOSE[n - 1][k - 1];
        }
    }
    for (u32 mask = 0; mask < (u32{1} << 15); ++mask) {
        u32 rank = 0;
        int ordinal = 0;
        for (int pos = 0; pos < 15; ++pos) {
            if (mask & (u32{1} << pos)) {
                ++ordinal;
                rank += CHOOSE[pos][ordinal];
            }
        }
        LOW_RANK[mask] = rank;
        for (int lower_count = 0; lower_count <= 10; ++lower_count) {
            u32 high = 0;
            int high_ordinal = 0;
            for (int pos = 0; pos < 15; ++pos) {
                if (mask & (u32{1} << pos)) {
                    ++high_ordinal;
                    const int ordinal2 = lower_count + high_ordinal;
                    if (ordinal2 <= 10) high += CHOOSE[15 + pos][ordinal2];
                }
            }
            HIGH_PART[lower_count][mask] = high;
        }
    }
}

u32 colex_rank(u32 mask) {
    const u32 low = mask & ((u32{1} << 15) - 1);
    const u32 high = mask >> 15;
    const int count_low = std::popcount(low);
    return LOW_RANK[low] + HIGH_PART[count_low][high];
}

struct Cell {
    i64 left;
    i64 right;
    bool qsafe;
    u32 failed;
};

bool safe_midpoint(int speed, i64 left, i64 right, i64 D) {
    const i64 residue = static_cast<i64>(
        (static_cast<i128>(speed) * (left + right)) %
        (2 * static_cast<i128>(D)));
    return static_cast<i128>(7) * residue >= D &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * D;
}

struct AtomDatum {
    i64 length = 0;
    int component_delta = 0;
};

struct Geometry {
    i64 denominator;
    std::vector<Cell> cells;
    std::map<u32, AtomDatum> atoms;
};

Geometry literal_geometry(int max_arity) {
    i64 D = 1;
    for (int speed : P) D = lcm_exact(D, 14LL * speed);
    D = lcm_exact(D, 14LL * Q0);

    std::vector<i64> walls = {0, D};
    auto add_walls = [&](int speed) {
        const i64 unit = D / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int speed : P) add_walls(speed);
    add_walls(Q0);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());

    std::vector<Cell> cells;
    cells.reserve(walls.size() - 1);
    for (std::size_t j = 0; j + 1 < walls.size(); ++j) {
        u32 failed = 0;
        for (int i = 0; i < 30; ++i) {
            if (!safe_midpoint(P[i], walls[j], walls[j + 1], D)) {
                failed |= u32{1} << i;
            }
        }
        cells.push_back({walls[j], walls[j + 1],
                         safe_midpoint(Q0, walls[j], walls[j + 1], D), failed});
    }

    std::map<u32, AtomDatum> atoms;
    for (const Cell& cell : cells) {
        if (cell.qsafe && std::popcount(cell.failed) <= max_arity) {
            atoms[cell.failed].length += cell.right - cell.left;
            ++atoms[cell.failed].component_delta;
        }
    }
    for (std::size_t j = 0; j < cells.size(); ++j) {
        const Cell& here = cells[j];
        const Cell& before = cells[(j + cells.size() - 1) % cells.size()];
        if (!here.qsafe || !before.qsafe) continue;
        const u32 joined = here.failed | before.failed;
        if (std::popcount(joined) <= max_arity) {
            --atoms[joined].component_delta;
        }
    }
    for (auto it = atoms.begin(); it != atoms.end();) {
        if (it->second.length == 0 && it->second.component_delta == 0) {
            it = atoms.erase(it);
        } else {
            ++it;
        }
    }
    return {D, std::move(cells), std::move(atoms)};
}

struct LayerArrays {
    int arity;
    std::vector<i64> mass;
    std::vector<int> components;
    u64 updates = 0;
};

struct SupersetAdder {
    int arity;
    u32 base;
    i64 mass_value;
    int component_value;
    std::array<int, 30> available{};
    int available_count = 0;
    LayerArrays* layer = nullptr;

    void visit(int start, int left, u32 selected) {
        if (left == 0) {
            const u32 index = colex_rank(base | selected);
            layer->mass[index] += mass_value;
            layer->components[index] += component_value;
            ++layer->updates;
            return;
        }
        const int final_start = available_count - left;
        for (int j = start; j <= final_start; ++j) {
            visit(j + 1, left - 1, selected | (u32{1} << available[j]));
        }
    }

    void run() {
        const int base_count = std::popcount(base);
        const int choose_count = arity - base_count;
        if (choose_count < 0) return;
        for (int pos = 0; pos < 30; ++pos) {
            if (!(base & (u32{1} << pos))) available[available_count++] = pos;
        }
        visit(0, choose_count, 0);
    }
};

LayerArrays build_arrays(const Geometry& geometry, int arity) {
    const u32 count = CHOOSE[30][arity];
    LayerArrays layer{arity, std::vector<i64>(count, 0),
                      std::vector<int>(count, 0), 0};
    for (const auto& [mask, datum] : geometry.atoms) {
        if (std::popcount(mask) > arity) continue;
        SupersetAdder adder;
        adder.arity = arity;
        adder.base = mask;
        adder.mass_value = datum.length;
        adder.component_value = datum.component_delta;
        adder.layer = &layer;
        adder.run();
    }
    return layer;
}

i128 ceil_div(i128 numerator, i128 denominator) {
    need(numerator >= 0 && denominator > 0, "bad ceil operands");
    return (numerator + denominator - 1) / denominator;
}

struct Deck {
    int arity;
    u64 bodies = 0;
    u64 strict_edges = 0;
    u64 equalities = 0;
    u64 active_minus_count = 0;
    u64 active_count = 0;
    u64 activating_count = 0;
    i128 minimum_cutoff = 0;
    i128 maximum_cutoff = 0;
    std::vector<u32> active_minus;
    std::vector<u32> active;
};

Deck classify(const Geometry& geometry, const LayerArrays& layer,
              bool retain_edges) {
    Deck deck;
    deck.arity = layer.arity;
    deck.bodies = layer.mass.size();
    if (retain_edges) {
        deck.active_minus.reserve(1000000);
        deck.active.reserve(1000000);
    }
    const u32 limit = u32{1} << 30;
    u32 mask = (u32{1} << layer.arity) - 1;
    u32 index = 0;
    while (mask && mask < limit) {
        need(colex_rank(mask) == index, "colex indexing failure");
        const i128 surplus = static_cast<i128>(27) * layer.mass[index] -
                             static_cast<i128>(2) * geometry.denominator;
        if (surplus == 0) ++deck.equalities;
        if (surplus > 0) {
            ++deck.strict_edges;
            const int c = layer.components[index];
            need(c > 0, "strict edge has no components");
            const i128 cutoff = ceil_div(
                static_cast<i128>(54) * c * geometry.denominator,
                static_cast<i128>(7) *
                    (static_cast<i128>(54) * layer.mass[index] -
                     static_cast<i128>(4) * geometry.denominator));
            if (deck.minimum_cutoff == 0 || cutoff < deck.minimum_cutoff) {
                deck.minimum_cutoff = cutoff;
            }
            if (cutoff > deck.maximum_cutoff) deck.maximum_cutoff = cutoff;
            if (cutoff <= TEST_Q_MINUS) {
                ++deck.active_minus_count;
                if (retain_edges) deck.active_minus.push_back(mask);
            }
            if (cutoff <= TEST_Q) {
                ++deck.active_count;
                if (retain_edges) deck.active.push_back(mask);
            }
            if (cutoff == TEST_Q) ++deck.activating_count;
        }
        ++index;
        const u32 next = next_combination(mask);
        if (next <= mask) break;
        mask = next;
    }
    need(index == layer.mass.size(), "combination count mismatch");
    return deck;
}

u64 splitmix(u64 x) {
    x += UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

void independent_order(std::vector<u32>& edges) {
    std::sort(edges.begin(), edges.end(), [](u32 a, u32 b) {
        const u64 ha = splitmix(static_cast<u64>(a) ^ UINT64_C(0x50e95682));
        const u64 hb = splitmix(static_cast<u64>(b) ^ UINT64_C(0x50e95682));
        return ha < hb || (ha == hb && a < b);
    });
}

bool is_cover(u32 body, const std::vector<u32>& edges, u64* checks = nullptr,
              u32* missed = nullptr) {
    u64 local = 0;
    for (u32 edge : edges) {
        ++local;
        if ((body & edge) == 0) {
            if (checks) *checks = local;
            if (missed) *missed = edge;
            return false;
        }
    }
    if (checks) *checks = local;
    if (missed) *missed = 0;
    return true;
}

struct Exhaustion {
    u64 bodies = 0;
    u64 covers = 0;
    u64 checks = 0;
    u64 maximum_checks = 0;
    u32 closest_body = 0;
    u32 closest_missed = 0;
    u32 first_cover = 0;
};

Exhaustion exhaust_nine(const std::vector<u32>& ordered_edges,
                        bool stop_at_first_cover) {
    Exhaustion out;
    const u32 limit = u32{1} << 30;
    u32 body = (u32{1} << 9) - 1;
    while (body && body < limit) {
        ++out.bodies;
        u64 checked = 0;
        u32 missed = 0;
        const bool cover = is_cover(body, ordered_edges, &checked, &missed);
        out.checks += checked;
        if (checked > out.maximum_checks) {
            out.maximum_checks = checked;
            out.closest_body = body;
            out.closest_missed = missed;
        }
        if (cover) {
            ++out.covers;
            if (!out.first_cover) out.first_cover = body;
            if (stop_at_first_cover) break;
        }
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    return out;
}

}  // namespace

int main() {
    init_choose_and_rank();
    const Geometry geometry = literal_geometry(9);
    need(geometry.denominator == INT64_C(91205797082400), "D50 changed");
    need(geometry.cells.size() == 7213, "joint atom count changed");

    const u32 claimed_body = mask_of({16,85,88,95,143,168,193,240,290});
    const u32 claimed_edge = mask_of({40,63,80,120,126,145,176,190,252});
    need(std::popcount(claimed_body) == 9 && std::popcount(claimed_edge) == 9,
         "claimed masks wrong arity");
    need((claimed_body & claimed_edge) == 0, "claimed masks not disjoint");

    LayerArrays e9_arrays = build_arrays(geometry, 9);
    Deck e9 = classify(geometry, e9_arrays, true);
    const u32 edge_index = colex_rank(claimed_edge);
    const i64 claimed_mass = e9_arrays.mass[edge_index];
    const int claimed_components = e9_arrays.components[edge_index];
    const i128 claimed_surplus = static_cast<i128>(27) * claimed_mass -
                                   static_cast<i128>(2) * geometry.denominator;
    const i128 claimed_cutoff = ceil_div(
        static_cast<i128>(54) * claimed_components * geometry.denominator,
        static_cast<i128>(7) *
            (static_cast<i128>(54) * claimed_mass -
             static_cast<i128>(4) * geometry.denominator));
    need(claimed_surplus == INT64_C(10402250150088),
         "claimed surplus mismatch");
    need(claimed_components == 168, "claimed component count mismatch");
    need(claimed_cutoff == TEST_Q, "claimed cutoff mismatch");
    need(std::find(e9.active_minus.begin(), e9.active_minus.end(), claimed_edge) ==
             e9.active_minus.end(),
         "claimed edge activates too early");
    need(std::find(e9.active.begin(), e9.active.end(), claimed_edge) !=
             e9.active.end(),
         "claimed edge absent at Q");

    need(is_cover(claimed_body, e9.active_minus),
         "claimed Q-1 body is not a cover");
    need(!is_cover(claimed_body, e9.active),
         "claimed Q body remained a cover");
    independent_order(e9.active);
    const Exhaustion e9_exhaust = exhaust_nine(e9.active, false);
    need(e9_exhaust.bodies == CHOOSE[30][9], "E9 body universe incomplete");
    need(e9_exhaust.covers == 0, "E9 Q deck has a nine-cover");

    LayerArrays e8_arrays = build_arrays(geometry, 8);
    Deck e8 = classify(geometry, e8_arrays, true);
    independent_order(e8.active);
    const Exhaustion e8_search = exhaust_nine(e8.active, true);

    std::cout << "AUDIT FIXED_Q50_E9_LITERAL_WALL_REFEREE\n";
    std::cout << "P_SIZE 30 Q0 50 D50 " << geometry.denominator
              << " WALLS " << geometry.cells.size() + 1
              << " OPEN_ATOMS " << geometry.cells.size()
              << " COMPRESSED_MASKS " << geometry.atoms.size() << "\n";
    std::cout << "E9 UNIVERSE " << e9.bodies
              << " SUPERSET_UPDATES " << e9_arrays.updates
              << " STRICT_EDGES " << e9.strict_edges
              << " EQUALITIES " << e9.equalities
              << " MIN_CUTOFF " << dec(e9.minimum_cutoff)
              << " MAX_CUTOFF " << dec(e9.maximum_cutoff) << "\n";
    std::cout << "E9 Q " << TEST_Q_MINUS
              << " ACTIVE_EDGES " << e9.active_minus_count
              << " COVER " << labels(claimed_body) << "\n";
    std::cout << "E9 Q " << TEST_Q
              << " ACTIVE_EDGES " << e9.active_count
              << " ACTIVATING_AT_Q " << e9.activating_count
              << " BODIES " << e9_exhaust.bodies
              << " COVERS " << e9_exhaust.covers
              << " EDGE_CHECKS " << e9_exhaust.checks
              << " MAX_CHECKS " << e9_exhaust.maximum_checks
              << " CLOSEST_BODY " << labels(e9_exhaust.closest_body)
              << " MISSED " << labels(e9_exhaust.closest_missed) << "\n";
    std::cout << "ACTIVATING_EDGE " << labels(claimed_edge)
              << " MASS " << claimed_mass
              << " SURPLUS_27M_MINUS_2D " << dec(claimed_surplus)
              << " COMPONENTS " << claimed_components
              << " CUTOFF " << dec(claimed_cutoff)
              << " DISJOINT_FROM_COVER YES\n";
    std::cout << "E8 UNIVERSE " << e8.bodies
              << " SUPERSET_UPDATES " << e8_arrays.updates
              << " STRICT_EDGES " << e8.strict_edges
              << " EQUALITIES " << e8.equalities
              << " Q " << TEST_Q
              << " ACTIVE_EDGES " << e8.active_count
              << " SEARCHED_BODIES " << e8_search.bodies
              << " COVER_FOUND " << (e8_search.covers ? "YES" : "NO");
    if (e8_search.covers) {
        std::cout << " COVER " << labels(e8_search.first_cover);
    }
    std::cout << "\n";
    std::cout << "VERDICT E9_Q5681_COVER_Q5682_NO_COVER_ACCEPT\n";
    return 0;
}
