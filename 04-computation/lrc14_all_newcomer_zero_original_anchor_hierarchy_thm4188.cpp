#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

// Primary exact audit for THM-4188.  This path keeps the fixed-pool wall
// arrangement and integrates each newcomer comb by an exact safe-prefix
// formula.  The independent audit instead builds the full joint arrangement.

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

namespace {

constexpr std::array<int, 30> POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
constexpr i64 COMMON = INT64_C(18241159416480);
constexpr int FINITE_LIMIT = 2586;
constexpr std::array<int, 15> EXPECTED_Q5_FAILURES = {
    6, 24, 25, 48, 70, 96, 105, 110, 128, 140, 192, 206, 210, 256, 366};
constexpr std::array<int, 19> EXPECTED_Q6_FAILURES = {
    6, 22, 24, 25, 48, 70, 72, 96, 105, 110, 128, 130, 140, 192,
    206, 210, 256, 260, 366};
constexpr std::array<int, 23> EXPECTED_Q7_FAILURES = {
    6, 22, 24, 25, 48, 70, 72, 96, 100, 105, 110, 128, 130, 140,
    186, 192, 206, 210, 220, 256, 260, 294, 366};
constexpr u64 FNV1A64_OFFSET = UINT64_C(0xcbf29ce484222325);
constexpr u64 FNV1A64_PRIME = UINT64_C(0x100000001b3);
constexpr u64 EXPECTED_LEDGER = UINT64_C(0x188e752eabb43991);

class Fnv1a64 {
  public:
    void add_u64_le(u64 word) {
        for (int shift = 0; shift < 64; shift += 8) {
            value_ ^= static_cast<std::uint8_t>((word >> shift) & UINT64_C(0xff));
            value_ *= FNV1A64_PRIME;
        }
    }
    u64 value() const { return value_; }

  private:
    u64 value_ = FNV1A64_OFFSET;
};

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

i64 exact_lcm(i64 left, i64 right) {
    const i64 divisor = std::gcd(left, right);
    const i128 value = static_cast<i128>(left / divisor) * right;
    require(value <= std::numeric_limits<i64>::max(), "LCM overflow");
    return static_cast<i64>(value);
}

u32 next_combination(u32 mask) {
    const u32 low = mask & (~mask + 1u);
    const u32 ripple = mask + low;
    if (ripple == 0) return 0;
    return ripple | (((mask ^ ripple) >> 2) / low);
}

bool midpoint_safe(int speed, i64 left, i64 right) {
    const i64 period = 2 * COMMON;
    const i64 residue = static_cast<i64>(
        (static_cast<i128>(speed) * (left + right)) % period);
    return static_cast<i128>(7) * residue >= COMMON &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * COMMON;
}

struct Cell {
    i64 left;
    i64 right;
    u32 failed_pool;
};

struct FixedSafeStats {
    std::unordered_map<u32, i64> length;
    std::unordered_map<u32, i64> cell_count;
    std::unordered_map<u32, i64> adjacency_union_count;
};

std::vector<Cell> build_pool_cells() {
    i64 check = 1;
    for (int speed : POOL) check = exact_lcm(check, 14LL * speed);
    require(check == COMMON, "fixed pool LCM changed");

    std::vector<i64> walls = {0, COMMON};
    for (int speed : POOL) {
        const i64 unit = COMMON / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    }
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    require(walls.size() == 7134, "fixed pool wall count changed");

    std::vector<Cell> cells;
    cells.reserve(walls.size() - 1);
    for (std::size_t index = 0; index + 1 < walls.size(); ++index) {
        u32 failed = 0;
        for (int vertex = 0; vertex < 30; ++vertex) {
            if (!midpoint_safe(POOL[vertex], walls[index], walls[index + 1])) {
                failed |= u32{1} << vertex;
            }
        }
        cells.push_back({walls[index], walls[index + 1], failed});
    }
    return cells;
}

FixedSafeStats build_fixed_safe_stats(const std::vector<Cell>& cells,
                                      int max_arity) {
    FixedSafeStats stats;
    stats.length.reserve(4096);
    stats.cell_count.reserve(4096);
    stats.adjacency_union_count.reserve(8192);
    for (const Cell& cell : cells) {
        if (std::popcount(cell.failed_pool) <= max_arity) {
            stats.length[cell.failed_pool] += cell.right - cell.left;
            ++stats.cell_count[cell.failed_pool];
        }
    }
    for (std::size_t index = 0; index < cells.size(); ++index) {
        const u32 joined = cells[index].failed_pool |
                           cells[(index + 1) % cells.size()].failed_pool;
        if (std::popcount(joined) <= max_arity) {
            ++stats.adjacency_union_count[joined];
        }
    }
    return stats;
}

i64 safe_prefix(i64 tick, int q) {
    const i128 product = static_cast<i128>(q) * tick;
    const i64 whole = static_cast<i64>(product / COMMON);
    const i64 rem = static_cast<i64>(product % COMMON);
    const i64 scaled = 14 * rem;
    i64 partial = 0;
    if (scaled <= COMMON) {
        partial = 0;
    } else if (scaled >= 13 * COMMON) {
        partial = 12 * COMMON;
    } else {
        partial = scaled - COMMON;
    }
    const i128 answer = static_cast<i128>(12) * whole * COMMON + partial;
    require(answer <= std::numeric_limits<i64>::max(), "prefix overflow");
    return static_cast<i64>(answer);
}

using AtomMass = std::unordered_map<u32, i64>;

struct AtomSupport {
    std::vector<u32> masks;
    std::unordered_map<u32, u32> index;
};

AtomSupport build_atom_support(const std::vector<Cell>& cells, int max_arity) {
    AtomSupport support;
    for (const Cell& cell : cells) {
        if (std::popcount(cell.failed_pool) <= max_arity) {
            support.masks.push_back(cell.failed_pool);
        }
    }
    std::sort(support.masks.begin(), support.masks.end());
    support.masks.erase(std::unique(support.masks.begin(), support.masks.end()),
                        support.masks.end());
    support.index.reserve(2 * support.masks.size());
    for (u32 index = 0; index < support.masks.size(); ++index) {
        support.index.emplace(support.masks[index], index);
    }
    return support;
}

AtomMass build_atom_mass(const std::vector<Cell>& cells, int q, int max_arity) {
    AtomMass mass;
    mass.reserve(4096);
    i64 previous = safe_prefix(0, q);
    i128 total = 0;
    for (const Cell& cell : cells) {
        const i64 current = safe_prefix(cell.right, q);
        const i64 contribution = current - previous;
        require(contribution >= 0, "negative safe contribution");
        if (std::popcount(cell.failed_pool) <= max_arity) {
            mass[cell.failed_pool] += contribution;
        }
        total += contribution;
        previous = current;
    }
    require(total == static_cast<i128>(12) * q * COMMON,
            "q-safe cell mass did not sum to 6/7");
    return mass;
}

std::vector<i64> build_dense_atom_mass(const std::vector<Cell>& cells,
                                       const AtomSupport& support,
                                       int q) {
    std::vector<i64> mass(support.masks.size(), 0);
    i64 previous = safe_prefix(0, q);
    i128 total = 0;
    for (const Cell& cell : cells) {
        const i64 current = safe_prefix(cell.right, q);
        const i64 contribution = current - previous;
        require(contribution >= 0, "negative safe contribution");
        const auto found = support.index.find(cell.failed_pool);
        if (found != support.index.end()) mass[found->second] += contribution;
        total += contribution;
        previous = current;
    }
    require(total == static_cast<i128>(12) * q * COMMON,
            "q-safe dense cell mass did not sum to 6/7");
    return mass;
}

i64 deletion_mass(u32 deletion, const AtomMass& atom_mass) {
    i64 mass = 0;
    u32 subset = deletion;
    while (true) {
        const auto found = atom_mass.find(subset);
        if (found != atom_mass.end()) mass += found->second;
        if (subset == 0) break;
        subset = (subset - 1) & deletion;
    }
    return mass;
}

i128 delta(u32 deletion, const AtomMass& atom_mass, int q) {
    const i64 mass = deletion_mass(deletion, atom_mass);
    return static_cast<i128>(9) * mass - static_cast<i128>(8) * q * COMMON;
}

std::string labels(u32 mask) {
    std::string text;
    for (int vertex = 0; vertex < 30; ++vertex) {
        if ((mask & (u32{1} << vertex)) == 0) continue;
        if (!text.empty()) text += ',';
        text += std::to_string(POOL[vertex]);
    }
    return text;
}

std::string decimal(i128 value) {
    if (value == 0) return "0";
    bool negative = value < 0;
    if (negative) value = -value;
    std::string text;
    while (value != 0) {
        text.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (negative) text.push_back('-');
    std::reverse(text.begin(), text.end());
    return text;
}

std::string hex64(u64 value) {
    std::ostringstream out;
    out << std::hex << std::nouppercase << std::setw(16) << std::setfill('0')
        << value;
    return out.str();
}

template <std::size_t N>
bool expected_contains(const std::array<int, N>& values, int q) {
    return std::find(values.begin(), values.end(), q) != values.end();
}

std::string integer_list(const std::vector<int>& values) {
    std::string text;
    for (int value : values) {
        if (!text.empty()) text += ',';
        text += std::to_string(value);
    }
    return text;
}

struct EdgeLayer {
    int arity;
    std::vector<u32> edges;
    u64 equalities = 0;
};

struct CompressedLayer {
    int arity;
    std::vector<u32> edges;
    std::vector<u32> offsets;
    std::vector<u32> atom_indices;
};

CompressedLayer compress_layer(const EdgeLayer& layer,
                               const AtomSupport& support) {
    CompressedLayer compressed;
    compressed.arity = layer.arity;
    compressed.edges = layer.edges;
    compressed.offsets.reserve(layer.edges.size() + 1);
    compressed.offsets.push_back(0);
    for (u32 edge : layer.edges) {
        u32 subset = edge;
        while (true) {
            const auto found = support.index.find(subset);
            if (found != support.index.end()) {
                compressed.atom_indices.push_back(found->second);
            }
            if (subset == 0) break;
            subset = (subset - 1) & edge;
        }
        compressed.offsets.push_back(
            static_cast<u32>(compressed.atom_indices.size()));
    }
    return compressed;
}

EdgeLayer build_layer(const AtomMass& atom_mass, int q, int arity) {
    EdgeLayer layer{arity, {}, 0};
    const u32 limit = u32{1} << 30;
    u32 deletion = (u32{1} << arity) - 1;
    while (deletion != 0 && deletion < limit) {
        const i128 value = delta(deletion, atom_mass, q);
        if (value == 0) ++layer.equalities;
        if (value >= 0) layer.edges.push_back(deletion);
        const u32 next = next_combination(deletion);
        if (next <= deletion) break;
        deletion = next;
    }
    return layer;
}

struct InclusionResult {
    u64 missing = 0;
    u64 equalities = 0;
    std::vector<u32> missing_edges;
    u32 first = 0;
    i128 first_delta = 0;
    i128 worst_delta = 0;
    u32 worst = 0;
};

InclusionResult test_inclusion(const CompressedLayer& base,
                               const std::vector<i64>& atom_mass,
                               int q) {
    InclusionResult result;
    bool have_worst = false;
    for (std::size_t edge_index = 0; edge_index < base.edges.size(); ++edge_index) {
        i64 mass = 0;
        for (u32 position = base.offsets[edge_index];
             position < base.offsets[edge_index + 1]; ++position) {
            mass += atom_mass[base.atom_indices[position]];
        }
        const i128 value = static_cast<i128>(9) * mass -
                           static_cast<i128>(8) * q * COMMON;
        const u32 edge = base.edges[edge_index];
        if (!have_worst || value < result.worst_delta) {
            result.worst_delta = value;
            result.worst = edge;
            have_worst = true;
        }
        if (value == 0) ++result.equalities;
        if (value < 0) {
            if (result.missing == 0) {
                result.first = edge;
                result.first_delta = value;
            }
            ++result.missing;
            result.missing_edges.push_back(edge);
        }
    }
    return result;
}

bool good(u32 mask) {
    int common_gcd = 0;
    for (int vertex = 0; vertex < 30; ++vertex) {
        if ((mask & (u32{1} << vertex)) != 0) {
            common_gcd = std::gcd(common_gcd, POOL[vertex]);
        }
    }
    if (common_gcd != 1) return false;
    for (int modulus = 2; modulus <= 14; ++modulus) {
        bool covered = false;
        for (int vertex = 0; vertex < 30; ++vertex) {
            if ((mask & (u32{1} << vertex)) != 0 &&
                POOL[vertex] % modulus == 0) {
                covered = true;
                break;
            }
        }
        if (!covered) return false;
    }
    return true;
}

using AnchorLayers = std::array<std::vector<u32>, 7>;

AnchorLayers enumerate_anchor_layers() {
    std::vector<int> ground;
    for (int vertex = 0; vertex < 30; ++vertex) {
        if (POOL[vertex] != 120 && POOL[vertex] != 126 && POOL[vertex] != 143) {
            ground.push_back(vertex);
        }
    }
    AnchorLayers layers;
    std::vector<u32> previous;
    for (int size = 1; size <= 6; ++size) {
        const u32 limit = u32{1} << ground.size();
        u32 local = (u32{1} << size) - 1;
        while (local != 0 && local < limit) {
            u32 mask = 0;
            u32 bits = local;
            while (bits != 0) {
                const u32 bit = bits & (~bits + 1u);
                mask |= u32{1} << ground[std::countr_zero(bit)];
                bits ^= bit;
            }
            if (good(mask)) {
                bool minimal = true;
                for (u32 old : previous) {
                    if ((mask & old) == old) {
                        minimal = false;
                        break;
                    }
                }
                if (minimal) layers[size].push_back(mask);
            }
            const u32 next = next_combination(local);
            if (next <= local) break;
            local = next;
        }
        previous.insert(previous.end(), layers[size].begin(), layers[size].end());
    }
    require(layers[4].size() == 32, "M4 anchor count changed");
    require(layers[5].size() == 297, "M5 anchor count changed");
    require(layers[6].size() == 24, "M6 anchor count changed");
    return layers;
}

struct RestrictedFailure {
    u64 affected = 0;
    u32 first_anchor = 0;
    u32 first_edge = 0;
};

RestrictedFailure restricted_failures(const std::vector<u32>& anchors,
                                      const std::vector<u32>& missing_edges) {
    RestrictedFailure result;
    for (u32 anchor : anchors) {
        for (u32 edge : missing_edges) {
            if ((anchor & edge) == 0) {
                if (result.affected == 0) {
                    result.first_anchor = anchor;
                    result.first_edge = edge;
                }
                ++result.affected;
                break;
            }
        }
    }
    return result;
}

template <class Callback>
void for_each_outside_anchor(u32 anchor, int cardinality, Callback&& callback) {
    std::vector<int> vertices;
    for (int vertex = 0; vertex < 30; ++vertex) {
        if ((anchor & (u32{1} << vertex)) == 0) vertices.push_back(vertex);
    }
    if (cardinality == 0) {
        callback(u32{0});
        return;
    }
    const u32 limit = u32{1} << vertices.size();
    u32 local = (u32{1} << cardinality) - 1;
    while (local != 0 && local < limit) {
        u32 global = 0;
        u32 bits = local;
        while (bits != 0) {
            const u32 bit = bits & (~bits + 1u);
            global |= u32{1} << vertices[std::countr_zero(bit)];
            bits ^= bit;
        }
        callback(global);
        const u32 next = next_combination(local);
        if (next <= local) break;
        local = next;
    }
}

bool restricted_cover(u32 chosen, u32 anchor, const std::vector<u32>& edges) {
    for (u32 edge : edges) {
        if ((edge & anchor) == 0 && (edge & chosen) == 0) return false;
    }
    return true;
}

std::vector<u32> complete_restricted_covers(u32 anchor,
                                            const std::vector<u32>& edges,
                                            int budget) {
    std::vector<u32> covers;
    for (int size = 0; size <= budget; ++size) {
        for_each_outside_anchor(anchor, size, [&](u32 candidate) {
            if (restricted_cover(candidate, anchor, edges)) {
                covers.push_back(candidate);
            }
        });
    }
    return covers;
}

std::vector<u32> filter_restricted_covers(const std::vector<u32>& candidates,
                                          u32 anchor,
                                          const std::vector<u32>& edges) {
    std::vector<u32> survivors;
    for (u32 candidate : candidates) {
        if (restricted_cover(candidate, anchor, edges)) {
            survivors.push_back(candidate);
        }
    }
    return survivors;
}

struct StaircaseAudit {
    u64 shallow_covers = 0;
    u64 middle_covers = 0;
    u64 final_covers = 0;
    u64 final_rows = 0;
    u32 first_final_anchor = 0;
    u32 first_final_cover = 0;
};

StaircaseAudit audit_staircase(const std::vector<u32>& anchors,
                               const std::vector<u32>& shallow_edges,
                               const std::vector<u32>* middle_edges,
                               const std::vector<u32>& final_edges,
                               int budget) {
    StaircaseAudit audit;
    for (u32 anchor : anchors) {
        std::vector<u32> covers =
            complete_restricted_covers(anchor, shallow_edges, budget);
        audit.shallow_covers += covers.size();
        if (middle_edges != nullptr) {
            covers = filter_restricted_covers(covers, anchor, *middle_edges);
        }
        audit.middle_covers += covers.size();
        covers = filter_restricted_covers(covers, anchor, final_edges);
        audit.final_covers += covers.size();
        if (!covers.empty()) {
            ++audit.final_rows;
            if (audit.first_final_anchor == 0) {
                audit.first_final_anchor = anchor;
                audit.first_final_cover = covers.front();
            }
        }
    }
    return audit;
}

struct LimitResult {
    u64 nonpositive = 0;
    i128 minimum_surplus_numerator = 0;
    u32 minimum_surplus_edge = 0;
    i128 maximum_bound_ceiling = 0;
    u32 maximum_bound_edge = 0;
    i64 maximum_bound_components = 0;
    i64 maximum_bound_mass = 0;
};

i64 subset_sum(u32 deletion, const std::unordered_map<u32, i64>& values) {
    i64 result = 0;
    u32 subset = deletion;
    while (true) {
        const auto found = values.find(subset);
        if (found != values.end()) result += found->second;
        if (subset == 0) break;
        subset = (subset - 1) & deletion;
    }
    return result;
}

LimitResult analyze_limit(const EdgeLayer& layer,
                          const FixedSafeStats& stats) {
    LimitResult result;
    bool first = true;
    for (u32 edge : layer.edges) {
        const i64 mass = subset_sum(edge, stats.length);
        const i64 safe_cells = subset_sum(edge, stats.cell_count);
        const i64 safe_adjacencies = subset_sum(edge, stats.adjacency_union_count);
        const i64 components = safe_cells - safe_adjacencies;
        require(mass > 0 && components > 0, "invalid fixed safe-set geometry");
        const i128 surplus_numerator =
            static_cast<i128>(27) * mass - static_cast<i128>(2) * COMMON;
        if (first || surplus_numerator < result.minimum_surplus_numerator) {
            result.minimum_surplus_numerator = surplus_numerator;
            result.minimum_surplus_edge = edge;
        }
        if (surplus_numerator <= 0) {
            ++result.nonpositive;
        } else {
            // From |mu(U cap G_q)-6 mu(U)/7| <= 6c/(49q),
            // q >= 27 c L / (7(27m-2L)) is sufficient.
            const i128 numerator = static_cast<i128>(27) * components * COMMON;
            const i128 denominator = static_cast<i128>(7) * surplus_numerator;
            const i128 ceiling = (numerator + denominator - 1) / denominator;
            if (ceiling > result.maximum_bound_ceiling) {
                result.maximum_bound_ceiling = ceiling;
                result.maximum_bound_edge = edge;
                result.maximum_bound_components = components;
                result.maximum_bound_mass = mass;
            }
        }
        first = false;
    }
    return result;
}

}  // namespace

int main() {
    try {
        const std::vector<Cell> cells = build_pool_cells();
        const AnchorLayers anchors = enumerate_anchor_layers();
        const FixedSafeStats fixed_stats = build_fixed_safe_stats(cells, 7);
        const AtomSupport support = build_atom_support(cells, 7);
        const AtomMass mass50 = build_atom_mass(cells, 50, 7);
        const std::array<EdgeLayer, 3> raw_base = {
            build_layer(mass50, 50, 5),
            build_layer(mass50, 50, 6),
            build_layer(mass50, 50, 7)};
        require(raw_base[0].edges.size() == 3017, "q50 E5 count changed");
        require(raw_base[1].edges.size() == 85324, "q50 E6 count changed");
        require(raw_base[2].edges.size() == 821737, "q50 E7 count changed");
        require(raw_base[0].equalities == 0 && raw_base[1].equalities == 0 &&
                    raw_base[2].equalities == 0,
                "q50 threshold equality appeared");
        const std::array<CompressedLayer, 3> base = {
            compress_layer(raw_base[0], support),
            compress_layer(raw_base[1], support),
            compress_layer(raw_base[2], support)};

        std::array<LimitResult, 3> limits;
        const std::array<i128, 3> expected_ceilings = {2297, 2443, 2587};
        for (int index = 0; index < 3; ++index) {
            limits[index] = analyze_limit(raw_base[index], fixed_stats);
            require(limits[index].nonpositive == 0,
                    "q50 repair has nonpositive strict-limit surplus");
            require(limits[index].maximum_bound_ceiling == expected_ceilings[index],
                    "discrepancy ceiling changed");
        }

        std::vector<int> failures5;
        std::vector<int> failures6;
        std::vector<int> failures7;
        std::array<u64, 3> inclusion_equalities = {0, 0, 0};
        u64 newcomer_count = 0;
        u64 residual_q5_equalities = 0;
        u64 residual_q6_equalities = 0;
        Fnv1a64 ledger;

        std::cout << "AUDIT all_newcomer_zero_original_anchor_hierarchy_pool_wall\n";
        std::cout << "POOL_CELLS " << cells.size()
                  << " COMMON " << COMMON
                  << " ATOM_SUPPORT " << support.masks.size()
                  << " ANCHORS_M4_M5_M6 " << anchors[4].size() << ','
                  << anchors[5].size() << ',' << anchors[6].size() << '\n';
        std::cout << "Q50_LAYERS E5 " << raw_base[0].edges.size()
                  << " E6 " << raw_base[1].edges.size()
                  << " E7 " << raw_base[2].edges.size()
                  << " EQUALITIES 0,0,0\n";
        for (int index = 0; index < 3; ++index) {
            const LimitResult& limit = limits[index];
            std::cout << "LIMIT D" << raw_base[index].arity
                      << " NONPOSITIVE " << limit.nonpositive
                      << " MIN_SURPLUS_NUM "
                      << decimal(limit.minimum_surplus_numerator)
                      << " MIN_SURPLUS_EDGE " << labels(limit.minimum_surplus_edge)
                      << " MAX_BOUND_CEIL "
                      << decimal(limit.maximum_bound_ceiling)
                      << " MAX_BOUND_EDGE " << labels(limit.maximum_bound_edge)
                      << " MAX_BOUND_COMPONENTS " << limit.maximum_bound_components
                      << " MAX_BOUND_MASS_NUM " << limit.maximum_bound_mass
                      << '\n';
        }

        for (int q = 1; q <= FINITE_LIMIT; ++q) {
            if (std::find(POOL.begin(), POOL.end(), q) != POOL.end()) continue;
            ++newcomer_count;
            const std::vector<i64> mass = build_dense_atom_mass(cells, support, q);
            std::array<InclusionResult, 3> inclusion;
            for (int index = 0; index < 3; ++index) {
                inclusion[index] = test_inclusion(base[index], mass, q);
                inclusion_equalities[index] += inclusion[index].equalities;
                ledger.add_u64_le(static_cast<u64>(q));
                ledger.add_u64_le(static_cast<u64>(raw_base[index].arity));
                ledger.add_u64_le(inclusion[index].missing);
                ledger.add_u64_le(inclusion[index].equalities);
                require(inclusion[index].worst_delta >=
                            std::numeric_limits<i64>::min() &&
                            inclusion[index].worst_delta <=
                            std::numeric_limits<i64>::max(),
                        "finite inclusion delta exceeds signed 64-bit");
                ledger.add_u64_le(static_cast<u64>(
                    static_cast<i64>(inclusion[index].worst_delta)));
            }
            if (inclusion[0].missing != 0) failures5.push_back(q);
            if (inclusion[1].missing != 0) failures6.push_back(q);
            if (inclusion[2].missing != 0) failures7.push_back(q);

            if (!expected_contains(EXPECTED_Q7_FAILURES, q)) continue;
            const RestrictedFailure affected4 =
                restricted_failures(anchors[4], inclusion[2].missing_edges);
            const RestrictedFailure affected5 =
                restricted_failures(anchors[5], inclusion[2].missing_edges);
            const RestrictedFailure affected6 =
                restricted_failures(anchors[6], inclusion[1].missing_edges);
            const AtomMass sparse_mass = build_atom_mass(cells, q, 7);
            const EdgeLayer q5 = build_layer(sparse_mass, q, 5);
            const EdgeLayer q6 = build_layer(sparse_mass, q, 6);
            residual_q5_equalities += q5.equalities;
            residual_q6_equalities += q6.equalities;
            const StaircaseAudit m4 = audit_staircase(
                anchors[4], q6.edges, nullptr, q6.edges, 6);
            const StaircaseAudit m5 = audit_staircase(
                anchors[5], q5.edges, nullptr, q5.edges, 5);
            const StaircaseAudit m6 = audit_staircase(
                anchors[6], q5.edges, nullptr, q5.edges, 4);
            require(m4.shallow_covers == 0 && m5.shallow_covers == 0 &&
                        m6.shallow_covers == 0,
                    "resonance retained a forbidden shallow cover");
            ledger.add_u64_le(q5.edges.size());
            ledger.add_u64_le(q6.edges.size());
            ledger.add_u64_le(q5.equalities);
            ledger.add_u64_le(q6.equalities);
            ledger.add_u64_le(affected4.affected);
            ledger.add_u64_le(affected5.affected);
            ledger.add_u64_le(affected6.affected);
            std::cout << "RESONANCE q " << q
                      << " LOST_E5_E6_E7 " << inclusion[0].missing << ','
                      << inclusion[1].missing << ',' << inclusion[2].missing
                      << " AFFECTED_M4_M5_M6 " << affected4.affected << ','
                      << affected5.affected << ',' << affected6.affected
                      << " NATIVE_E5_E6 " << q5.edges.size() << ','
                      << q6.edges.size()
                      << " NATIVE_EQUALITIES " << q5.equalities << ','
                      << q6.equalities
                      << " COVERS_M4E6_M5E5_M6E5 0,0,0\n";
        }

        require(newcomer_count == 2556, "finite newcomer count changed");
        require(failures5 == std::vector<int>(EXPECTED_Q5_FAILURES.begin(),
                                              EXPECTED_Q5_FAILURES.end()),
                "E5 inclusion failure set changed");
        require(failures6 == std::vector<int>(EXPECTED_Q6_FAILURES.begin(),
                                              EXPECTED_Q6_FAILURES.end()),
                "E6 inclusion failure set changed");
        require(failures7 == std::vector<int>(EXPECTED_Q7_FAILURES.begin(),
                                              EXPECTED_Q7_FAILURES.end()),
                "E7 inclusion failure set changed");

        require(ledger.value() == EXPECTED_LEDGER,
                "finite semantic ledger changed");

        std::cout << "FINITE_UNIVERSE Q_LE " << FINITE_LIMIT
                  << " NEWCOMERS " << newcomer_count << '\n';
        std::cout << "FAILURE_Q5 " << failures5.size() << ' '
                  << integer_list(failures5) << '\n';
        std::cout << "FAILURE_Q6 " << failures6.size() << ' '
                  << integer_list(failures6) << '\n';
        std::cout << "FAILURE_Q7 " << failures7.size() << ' '
                  << integer_list(failures7) << '\n';
        std::cout << "INCLUSION_EQUALITIES E5_E6_E7 "
                  << inclusion_equalities[0] << ',' << inclusion_equalities[1]
                  << ',' << inclusion_equalities[2] << '\n';
        std::cout << "RESONANCE_NATIVE_EQUALITIES E5_E6 "
                  << residual_q5_equalities << ',' << residual_q6_equalities
                  << '\n';
        std::cout << "FINITE_LEDGER_FNV1A64_LE " << hex64(ledger.value()) << '\n';
        std::cout << "COFINAL_INCLUSION q>=2587 BY_COMPONENT_DISCREPANCY\n";
        std::cout << "BODY_CONSEQUENCE ALL_Q_OUTSIDE_P GOOD_ZERO_ORIGINAL_BODIES_PER_Q "
                  << 1491665 << '\n';
        std::cout << "VERDICT PASS\n";
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
