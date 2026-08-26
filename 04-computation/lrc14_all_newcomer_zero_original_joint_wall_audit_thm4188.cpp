#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Independent exact audit for THM-4188.
//
// Unlike the primary pool-cell/prefix integration, this program explicitly
// forms the common refinement of every fixed-pool wall and every newcomer-q
// wall.  It measures each joint atom in the joint integer denominator.  The
// residual cover calculation is also independent: it is a bounded recursive
// transversal search, with an exhaustive candidate control on q=6.

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

namespace {

constexpr std::array<int, 30> POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
constexpr i64 POOL_DENOMINATOR = INT64_C(18241159416480);
constexpr int FINITE_LIMIT = 2586;

constexpr std::array<int, 15> EXPECTED_Q5 = {
    6, 24, 25, 48, 70, 96, 105, 110, 128, 140, 192, 206, 210, 256, 366};
constexpr std::array<int, 19> EXPECTED_Q6 = {
    6, 22, 24, 25, 48, 70, 72, 96, 105, 110, 128, 130, 140, 192,
    206, 210, 256, 260, 366};
constexpr std::array<int, 23> EXPECTED_Q7 = {
    6, 22, 24, 25, 48, 70, 72, 96, 100, 105, 110, 128, 130, 140,
    186, 192, 206, 210, 220, 256, 260, 294, 366};

struct NativeExpected {
    int q;
    u64 e5;
    u64 e6;
};

constexpr std::array<NativeExpected, 23> EXPECTED_NATIVE = {{
    {6, 50160, 389365}, {22, 54396, 421289},
    {24, 51188, 390537}, {25, 34684, 320668},
    {48, 77327, 475101}, {70, 73996, 469924},
    {72, 59802, 418692}, {96, 21457, 247722},
    {100, 26410, 284863}, {105, 37139, 325808},
    {110, 44484, 365850}, {128, 31275, 301247},
    {130, 105407, 539067}, {140, 59351, 416495},
    {186, 61671, 432488}, {192, 43980, 364769},
    {206, 65683, 451190}, {210, 26533, 261505},
    {220, 77964, 479384}, {256, 25765, 273944},
    {260, 61262, 422570}, {294, 68909, 454315},
    {366, 57331, 422399}
}};

constexpr u64 FNV_OFFSET = UINT64_C(0xcbf29ce484222325);
constexpr u64 FNV_PRIME = UINT64_C(0x100000001b3);
constexpr u64 EXPECTED_LEDGER = UINT64_C(0x6016229d0317ff06);

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

i64 checked_lcm(i64 a, i64 b) {
    const i64 g = std::gcd(a, b);
    const i128 result = static_cast<i128>(a / g) * b;
    require(result <= std::numeric_limits<i64>::max(), "LCM overflow");
    return static_cast<i64>(result);
}

bool safe_at_midpoint(int speed, i64 left, i64 right, i64 denominator) {
    const i128 twice_denominator = static_cast<i128>(2) * denominator;
    const i128 raw = static_cast<i128>(speed) * (left + right);
    const i64 residue = static_cast<i64>(raw % twice_denominator);
    return static_cast<i128>(7) * residue >= denominator &&
           static_cast<i128>(7) * residue <=
               static_cast<i128>(13) * denominator;
}

std::string decimal(i128 value) {
    if (value == 0) return "0";
    bool negative = value < 0;
    if (negative) value = -value;
    std::string result;
    while (value != 0) {
        result.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (negative) result.push_back('-');
    std::reverse(result.begin(), result.end());
    return result;
}

std::string hex64(u64 value) {
    std::ostringstream out;
    out << std::hex << std::nouppercase << std::setfill('0') << std::setw(16)
        << value;
    return out.str();
}

class Ledger {
  public:
    void add(u64 word) {
        for (int shift = 0; shift < 64; shift += 8) {
            value_ ^= static_cast<std::uint8_t>((word >> shift) & 0xffu);
            value_ *= FNV_PRIME;
        }
    }
    u64 value() const { return value_; }

  private:
    u64 value_ = FNV_OFFSET;
};

std::string labels(u32 mask) {
    std::string result;
    for (int index = 0; index < 30; ++index) {
        if ((mask & (u32{1} << index)) == 0) continue;
        if (!result.empty()) result += ',';
        result += std::to_string(POOL[index]);
    }
    return result;
}

u32 label_mask(std::initializer_list<int> values) {
    u32 result = 0;
    for (int value : values) {
        const auto found = std::find(POOL.begin(), POOL.end(), value);
        require(found != POOL.end(), "hostile-control label absent from pool");
        result |= u32{1} << std::distance(POOL.begin(), found);
    }
    return result;
}

bool is_pool_label(int q) {
    return std::find(POOL.begin(), POOL.end(), q) != POOL.end();
}

template <class Callback>
void choose_recursive(int n, int need, int next, u32 mask, Callback& callback) {
    if (need == 0) {
        callback(mask);
        return;
    }
    for (int vertex = next; vertex <= n - need; ++vertex) {
        choose_recursive(n, need - 1, vertex + 1,
                         mask | (u32{1} << vertex), callback);
    }
}

template <class Callback>
void for_each_k_subset(int n, int k, Callback&& callback) {
    auto local = std::forward<Callback>(callback);
    choose_recursive(n, k, 0, 0, local);
}

struct FixedCell {
    i64 left;
    i64 right;
    u32 failed;
};

struct FixedGeometry {
    std::vector<i64> walls;
    std::vector<FixedCell> cells;
};

FixedGeometry build_fixed_geometry() {
    i64 lcm = 1;
    for (int speed : POOL) lcm = checked_lcm(lcm, 14LL * speed);
    require(lcm == POOL_DENOMINATOR, "fixed-pool denominator changed");

    FixedGeometry geometry;
    geometry.walls = {0, POOL_DENOMINATOR};
    for (int speed : POOL) {
        const i64 unit = POOL_DENOMINATOR / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            geometry.walls.push_back((14LL * tooth + 1) * unit);
            geometry.walls.push_back((14LL * tooth + 13) * unit);
        }
    }
    std::sort(geometry.walls.begin(), geometry.walls.end());
    geometry.walls.erase(
        std::unique(geometry.walls.begin(), geometry.walls.end()),
        geometry.walls.end());
    require(geometry.walls.size() == 7134, "fixed wall count changed");

    geometry.cells.reserve(geometry.walls.size() - 1);
    for (std::size_t i = 0; i + 1 < geometry.walls.size(); ++i) {
        u32 failed = 0;
        for (int vertex = 0; vertex < 30; ++vertex) {
            if (!safe_at_midpoint(POOL[vertex], geometry.walls[i],
                                  geometry.walls[i + 1], POOL_DENOMINATOR)) {
                failed |= u32{1} << vertex;
            }
        }
        geometry.cells.push_back(
            {geometry.walls[i], geometry.walls[i + 1], failed});
    }
    return geometry;
}

struct Support {
    std::vector<u32> masks;
    std::unordered_map<u32, u32> position;
};

Support build_support(const FixedGeometry& fixed) {
    Support support;
    for (const FixedCell& cell : fixed.cells) {
        if (std::popcount(cell.failed) <= 7) support.masks.push_back(cell.failed);
    }
    std::sort(support.masks.begin(), support.masks.end());
    support.masks.erase(
        std::unique(support.masks.begin(), support.masks.end()),
        support.masks.end());
    support.position.reserve(2 * support.masks.size());
    for (u32 i = 0; i < support.masks.size(); ++i) {
        support.position.emplace(support.masks[i], i);
    }
    require(support.masks.size() == 2535, "low-arity atom support changed");
    return support;
}

std::vector<i64> joint_walls(const FixedGeometry& fixed, int q,
                             i64 denominator) {
    require(denominator % POOL_DENOMINATOR == 0,
            "joint denominator does not refine fixed pool");
    const i64 scale = denominator / POOL_DENOMINATOR;
    std::vector<i64> walls;
    walls.reserve(fixed.walls.size() + 2 * static_cast<std::size_t>(q));
    for (i64 wall : fixed.walls) walls.push_back(wall * scale);
    const i64 unit = denominator / (14LL * q);
    for (int tooth = 0; tooth < q; ++tooth) {
        walls.push_back((14LL * tooth + 1) * unit);
        walls.push_back((14LL * tooth + 13) * unit);
    }
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    return walls;
}

struct JointMass {
    i64 denominator = 0;
    std::vector<i64> atom;
    u64 joint_cells = 0;
};

JointMass build_joint_mass(const FixedGeometry& fixed, const Support& support,
                           int q) {
    JointMass result;
    result.denominator = checked_lcm(POOL_DENOMINATOR, 14LL * q);
    result.atom.assign(support.masks.size(), 0);
    const i64 scale = result.denominator / POOL_DENOMINATOR;
    const std::vector<i64> walls = joint_walls(fixed, q, result.denominator);
    result.joint_cells = walls.size() - 1;

    std::size_t fixed_index = 0;
    i128 safe_total = 0;
    for (std::size_t i = 0; i + 1 < walls.size(); ++i) {
        const i64 left = walls[i];
        const i64 right = walls[i + 1];
        while (fixed_index + 1 < fixed.walls.size() &&
               fixed.walls[fixed_index + 1] * scale <= left) {
            ++fixed_index;
        }
        require(fixed_index < fixed.cells.size(),
                "joint sweep escaped fixed cells");
        require(right <= fixed.walls[fixed_index + 1] * scale,
                "joint atom crossed a fixed-pool wall");
        if (!safe_at_midpoint(q, left, right, result.denominator)) continue;
        const i64 length = right - left;
        safe_total += length;
        const auto found = support.position.find(fixed.cells[fixed_index].failed);
        if (found != support.position.end()) result.atom[found->second] += length;
    }
    require(static_cast<i128>(7) * safe_total ==
                static_cast<i128>(6) * result.denominator,
            "joint q-safe mass changed");
    return result;
}

i64 subset_mass(u32 deletion, const JointMass& mass, const Support& support) {
    i64 result = 0;
    u32 subset = deletion;
    while (true) {
        const auto found = support.position.find(subset);
        if (found != support.position.end()) result += mass.atom[found->second];
        if (subset == 0) break;
        subset = (subset - 1) & deletion;
    }
    return result;
}

i128 threshold_delta(u32 deletion, const JointMass& mass,
                     const Support& support) {
    return static_cast<i128>(63) * subset_mass(deletion, mass, support) -
           static_cast<i128>(4) * mass.denominator;
}

// This hostile path deliberately ignores fixed-cell failure masks.  It tests
// every undeleted pool runner afresh on the explicit joint partition.
i64 direct_joint_mass(const FixedGeometry& fixed, int q, u32 deletion) {
    const i64 denominator = checked_lcm(POOL_DENOMINATOR, 14LL * q);
    const std::vector<i64> walls = joint_walls(fixed, q, denominator);
    i64 total = 0;
    for (std::size_t i = 0; i + 1 < walls.size(); ++i) {
        if (!safe_at_midpoint(q, walls[i], walls[i + 1], denominator)) continue;
        bool safe = true;
        for (int vertex = 0; vertex < 30; ++vertex) {
            if ((deletion & (u32{1} << vertex)) != 0) continue;
            if (!safe_at_midpoint(POOL[vertex], walls[i], walls[i + 1],
                                  denominator)) {
                safe = false;
                break;
            }
        }
        if (safe) total += walls[i + 1] - walls[i];
    }
    return total;
}

struct Layer {
    int arity = 0;
    std::vector<u32> edges;
    u64 equalities = 0;
};

Layer build_layer(const JointMass& mass, const Support& support, int arity) {
    Layer result;
    result.arity = arity;
    for_each_k_subset(30, arity, [&](u32 deletion) {
        const i128 delta = threshold_delta(deletion, mass, support);
        if (delta == 0) ++result.equalities;
        if (delta >= 0) result.edges.push_back(deletion);
    });
    return result;
}

struct IncidenceLayer {
    int arity = 0;
    std::vector<u32> edges;
    std::vector<u32> offsets;
    std::vector<u32> atom_positions;
};

IncidenceLayer build_incidence(const Layer& layer, const Support& support) {
    IncidenceLayer result;
    result.arity = layer.arity;
    result.edges = layer.edges;
    result.offsets.reserve(layer.edges.size() + 1);
    result.offsets.push_back(0);
    for (u32 edge : layer.edges) {
        u32 subset = edge;
        while (true) {
            const auto found = support.position.find(subset);
            if (found != support.position.end()) {
                result.atom_positions.push_back(found->second);
            }
            if (subset == 0) break;
            subset = (subset - 1) & edge;
        }
        result.offsets.push_back(
            static_cast<u32>(result.atom_positions.size()));
    }
    return result;
}

struct Inclusion {
    u64 missing = 0;
    u64 equalities = 0;
    std::vector<u32> missing_edges;
};

Inclusion test_inclusion(const IncidenceLayer& base, const JointMass& mass) {
    Inclusion result;
    for (std::size_t i = 0; i < base.edges.size(); ++i) {
        i64 measured = 0;
        const u32 begin = base.offsets[i];
        const u32 end = base.offsets[i + 1];
        for (u32 j = begin; j < end; ++j) {
            measured += mass.atom[base.atom_positions[j]];
        }
        const i128 delta = static_cast<i128>(63) * measured -
                           static_cast<i128>(4) * mass.denominator;
        if (delta == 0) ++result.equalities;
        if (delta < 0) {
            ++result.missing;
            result.missing_edges.push_back(base.edges[i]);
        }
    }
    return result;
}

bool good_anchor(u32 mask) {
    int divisor = 0;
    for (int vertex = 0; vertex < 30; ++vertex) {
        if ((mask & (u32{1} << vertex)) != 0) {
            divisor = std::gcd(divisor, POOL[vertex]);
        }
    }
    if (divisor != 1) return false;
    for (int modulus = 2; modulus <= 14; ++modulus) {
        bool hit = false;
        for (int vertex = 0; vertex < 30; ++vertex) {
            if ((mask & (u32{1} << vertex)) != 0 &&
                POOL[vertex] % modulus == 0) {
                hit = true;
                break;
            }
        }
        if (!hit) return false;
    }
    return true;
}

using AnchorLayers = std::array<std::vector<u32>, 7>;

AnchorLayers build_anchor_layers() {
    u32 allowed = (u32{1} << 30) - 1;
    for (int label : {120, 126, 143}) {
        allowed &= ~label_mask({label});
    }
    std::vector<int> relabel;
    for (int vertex = 0; vertex < 30; ++vertex) {
        if ((allowed & (u32{1} << vertex)) != 0) relabel.push_back(vertex);
    }

    AnchorLayers layers;
    std::vector<u32> earlier;
    for (int size = 1; size <= 6; ++size) {
        for_each_k_subset(static_cast<int>(relabel.size()), size, [&](u32 local) {
            u32 global = 0;
            for (int bit = 0; bit < static_cast<int>(relabel.size()); ++bit) {
                if ((local & (u32{1} << bit)) != 0) {
                    global |= u32{1} << relabel[bit];
                }
            }
            if (!good_anchor(global)) return;
            for (u32 old : earlier) {
                if ((global & old) == old) return;
            }
            layers[size].push_back(global);
        });
        earlier.insert(earlier.end(), layers[size].begin(), layers[size].end());
    }
    require(layers[4].size() == 32, "M4 count changed");
    require(layers[5].size() == 297, "M5 count changed");
    require(layers[6].size() == 24, "M6 count changed");
    return layers;
}

struct LimitAudit {
    u64 nonpositive = 0;
    i128 minimum_surplus = 0;
    u32 minimum_edge = 0;
    i128 maximum_ceiling = 0;
    u32 maximum_edge = 0;
    i64 maximum_components = 0;
    i64 maximum_mass = 0;
};

LimitAudit audit_limit(const FixedGeometry& fixed, const Layer& layer) {
    LimitAudit result;
    bool first = true;
    std::vector<unsigned char> safe(fixed.cells.size());
    for (u32 edge : layer.edges) {
        i64 mass = 0;
        u64 safe_cells = 0;
        for (std::size_t i = 0; i < fixed.cells.size(); ++i) {
            safe[i] = static_cast<unsigned char>((fixed.cells[i].failed & ~edge) == 0);
            if (safe[i] != 0) {
                mass += fixed.cells[i].right - fixed.cells[i].left;
                ++safe_cells;
            }
        }
        i64 components = 0;
        for (std::size_t i = 0; i < safe.size(); ++i) {
            const std::size_t previous = (i + safe.size() - 1) % safe.size();
            if (safe[i] != 0 && safe[previous] == 0) ++components;
        }
        if (safe_cells == safe.size()) components = 1;
        require(mass > 0 && components > 0, "invalid limiting safe geometry");

        const i128 surplus = static_cast<i128>(27) * mass -
                             static_cast<i128>(2) * POOL_DENOMINATOR;
        if (first || surplus < result.minimum_surplus) {
            result.minimum_surplus = surplus;
            result.minimum_edge = edge;
        }
        if (surplus <= 0) {
            ++result.nonpositive;
        } else {
            const i128 numerator = static_cast<i128>(27) * components *
                                   POOL_DENOMINATOR;
            const i128 denominator = static_cast<i128>(7) * surplus;
            const i128 ceiling = (numerator + denominator - 1) / denominator;
            if (ceiling > result.maximum_ceiling) {
                result.maximum_ceiling = ceiling;
                result.maximum_edge = edge;
                result.maximum_components = components;
                result.maximum_mass = mass;
            }
        }
        first = false;
    }
    return result;
}

u64 restricted_affected(const std::vector<u32>& anchors,
                        const std::vector<u32>& missing_edges) {
    u64 result = 0;
    for (u32 anchor : anchors) {
        if (std::any_of(missing_edges.begin(), missing_edges.end(),
                        [anchor](u32 edge) { return (anchor & edge) == 0; })) {
            ++result;
        }
    }
    return result;
}

std::vector<u32> restricted_edges(u32 anchor, const Layer& layer) {
    std::vector<u32> result;
    for (u32 edge : layer.edges) {
        if ((edge & anchor) == 0) result.push_back(edge);
    }
    return result;
}

bool cover_search(const std::vector<u32>& edges, u32 chosen, int remaining,
                  std::unordered_set<u32>& dead, u64& nodes) {
    ++nodes;
    u32 uncovered = 0;
    for (u32 edge : edges) {
        if ((edge & chosen) == 0) {
            uncovered = edge;
            break;
        }
    }
    if (uncovered == 0) return true;
    if (remaining == 0) return false;
    if (dead.find(chosen) != dead.end()) return false;

    u32 options = uncovered;
    while (options != 0) {
        const u32 bit = options & (~options + 1u);
        if (cover_search(edges, chosen | bit, remaining - 1, dead, nodes)) {
            return true;
        }
        options ^= bit;
    }
    dead.insert(chosen);
    return false;
}

struct CoverAudit {
    u64 rows_with_cover = 0;
    u64 search_nodes = 0;
};

CoverAudit audit_covers(const std::vector<u32>& anchors, const Layer& layer,
                        int budget) {
    CoverAudit result;
    for (u32 anchor : anchors) {
        const std::vector<u32> edges = restricted_edges(anchor, layer);
        require(!edges.empty(), "restricted repair layer unexpectedly empty");
        std::unordered_set<u32> dead;
        if (cover_search(edges, 0, budget, dead, result.search_nodes)) {
            ++result.rows_with_cover;
        }
    }
    return result;
}

template <class Callback>
void choose_from_vertices(const std::vector<int>& vertices, int need,
                          int next, u32 mask, Callback& callback) {
    if (need == 0) {
        callback(mask);
        return;
    }
    for (int i = next; i <= static_cast<int>(vertices.size()) - need; ++i) {
        choose_from_vertices(vertices, need - 1, i + 1,
                             mask | (u32{1} << vertices[i]), callback);
    }
}

bool exhaustive_exact_cover(u32 anchor, const Layer& layer, int budget) {
    std::vector<int> vertices;
    for (int vertex = 0; vertex < 30; ++vertex) {
        if ((anchor & (u32{1} << vertex)) == 0) vertices.push_back(vertex);
    }
    const std::vector<u32> edges = restricted_edges(anchor, layer);
    bool found = false;
    auto callback = [&](u32 candidate) {
        if (found) return;
        bool cover = true;
        for (u32 edge : edges) {
            if ((edge & candidate) == 0) {
                cover = false;
                break;
            }
        }
        if (cover) found = true;
    };
    choose_from_vertices(vertices, budget, 0, 0, callback);
    return found;
}

void cover_solver_controls() {
    const std::vector<u32> positive = {
        (u32{1} << 0) | (u32{1} << 1),
        (u32{1} << 1) | (u32{1} << 2)};
    const std::vector<u32> negative = {
        (u32{1} << 0) | (u32{1} << 1),
        (u32{1} << 2) | (u32{1} << 3)};
    u64 nodes = 0;
    std::unordered_set<u32> dead;
    require(cover_search(positive, 0, 1, dead, nodes),
            "positive cover control failed");
    nodes = 0;
    dead.clear();
    require(!cover_search(negative, 0, 1, dead, nodes),
            "negative cover control failed");
}

struct FiniteRecord {
    bool newcomer = false;
    std::array<Inclusion, 3> inclusion;
};

std::string integer_list(const std::vector<int>& values) {
    std::string result;
    for (int value : values) {
        if (!result.empty()) result += ',';
        result += std::to_string(value);
    }
    return result;
}

}  // namespace

int main() {
    try {
        const FixedGeometry fixed = build_fixed_geometry();
        const Support support = build_support(fixed);
        const AnchorLayers anchors = build_anchor_layers();
        cover_solver_controls();

        const JointMass mass50 = build_joint_mass(fixed, support, 50);
        const std::array<Layer, 3> base = {
            build_layer(mass50, support, 5),
            build_layer(mass50, support, 6),
            build_layer(mass50, support, 7)};
        require(base[0].edges.size() == 3017, "q50 E5 count changed");
        require(base[1].edges.size() == 85324, "q50 E6 count changed");
        require(base[2].edges.size() == 821737, "q50 E7 count changed");
        require(base[0].equalities == 0 && base[1].equalities == 0 &&
                    base[2].equalities == 0,
                "q50 threshold equality appeared");

        const std::array<IncidenceLayer, 3> incidence = {
            build_incidence(base[0], support),
            build_incidence(base[1], support),
            build_incidence(base[2], support)};

        std::array<LimitAudit, 3> limit;
        constexpr std::array<i128, 3> EXPECTED_CEIL = {2297, 2443, 2587};
        for (int i = 0; i < 3; ++i) {
            limit[i] = audit_limit(fixed, base[i]);
            require(limit[i].nonpositive == 0,
                    "nonpositive limiting surplus appeared");
            require(limit[i].maximum_ceiling == EXPECTED_CEIL[i],
                    "component-discrepancy ceiling changed");
        }

        const std::array<std::pair<int, u32>, 4> direct_controls = {{
            {50, base[0].edges.front()},
            {6, label_mask({40, 80, 145, 168, 290})},
            {366, label_mask({8, 15, 145, 193, 290})},
            {367, label_mask({8, 95, 145, 168, 240})}
        }};
        for (const auto& [q, edge] : direct_controls) {
            const JointMass mass = build_joint_mass(fixed, support, q);
            require(subset_mass(edge, mass, support) ==
                        direct_joint_mass(fixed, q, edge),
                    "direct joint-wall hostile control disagreed");
        }

        std::vector<FiniteRecord> records(FINITE_LIMIT + 1);
        std::atomic<int> next_q{1};
        std::exception_ptr thread_error;
        std::mutex error_mutex;
        const unsigned detected = std::thread::hardware_concurrency();
        const unsigned worker_count = std::min(8u, std::max(1u, detected));
        std::vector<std::thread> workers;
        workers.reserve(worker_count);
        for (unsigned worker = 0; worker < worker_count; ++worker) {
            workers.emplace_back([&]() {
                try {
                    while (true) {
                        const int q = next_q.fetch_add(1);
                        if (q > FINITE_LIMIT) break;
                        if (is_pool_label(q)) continue;
                        JointMass mass = build_joint_mass(fixed, support, q);
                        FiniteRecord record;
                        record.newcomer = true;
                        for (int layer = 0; layer < 3; ++layer) {
                            record.inclusion[layer] =
                                test_inclusion(incidence[layer], mass);
                        }
                        records[q] = std::move(record);
                    }
                } catch (...) {
                    std::lock_guard<std::mutex> lock(error_mutex);
                    if (thread_error == nullptr) thread_error = std::current_exception();
                }
            });
        }
        for (std::thread& worker : workers) worker.join();
        if (thread_error != nullptr) std::rethrow_exception(thread_error);

        std::vector<int> failures5;
        std::vector<int> failures6;
        std::vector<int> failures7;
        std::array<u64, 3> inclusion_equalities = {0, 0, 0};
        u64 newcomer_count = 0;
        for (int q = 1; q <= FINITE_LIMIT; ++q) {
            if (!records[q].newcomer) continue;
            ++newcomer_count;
            for (int layer = 0; layer < 3; ++layer) {
                inclusion_equalities[layer] += records[q].inclusion[layer].equalities;
            }
            if (records[q].inclusion[0].missing != 0) failures5.push_back(q);
            if (records[q].inclusion[1].missing != 0) failures6.push_back(q);
            if (records[q].inclusion[2].missing != 0) failures7.push_back(q);
        }
        require(newcomer_count == 2556, "finite newcomer count changed");
        require(failures5 == std::vector<int>(EXPECTED_Q5.begin(), EXPECTED_Q5.end()),
                "independent E5 failure set changed");
        require(failures6 == std::vector<int>(EXPECTED_Q6.begin(), EXPECTED_Q6.end()),
                "independent E6 failure set changed");
        require(failures7 == std::vector<int>(EXPECTED_Q7.begin(), EXPECTED_Q7.end()),
                "independent E7 failure set changed");
        require(inclusion_equalities == std::array<u64, 3>{0, 0, 0},
                "finite inclusion equality appeared");

        std::cout << "AUDIT all_newcomer_zero_original_joint_wall_independent\n";
        std::cout << "FIXED_POOL_CELLS " << fixed.cells.size()
                  << " POOL_DENOMINATOR " << POOL_DENOMINATOR
                  << " LOW_ARITY_SUPPORT " << support.masks.size()
                  << " ANCHORS_M4_M5_M6 " << anchors[4].size() << ','
                  << anchors[5].size() << ',' << anchors[6].size() << '\n';
        std::cout << "METHOD EXPLICIT_JOINT_WALLS PARALLEL_FINITE_SCAN "
                     "RECURSIVE_TRANSVERSALS\n";
        std::cout << "Q50_LAYERS E5 " << base[0].edges.size()
                  << " E6 " << base[1].edges.size()
                  << " E7 " << base[2].edges.size()
                  << " EQUALITIES 0,0,0\n";
        for (int i = 0; i < 3; ++i) {
            std::cout << "LIMIT D" << base[i].arity
                      << " NONPOSITIVE " << limit[i].nonpositive
                      << " MIN_SURPLUS_NUM " << decimal(limit[i].minimum_surplus)
                      << " MIN_SURPLUS_EDGE " << labels(limit[i].minimum_edge)
                      << " MAX_BOUND_CEIL " << decimal(limit[i].maximum_ceiling)
                      << " MAX_BOUND_EDGE " << labels(limit[i].maximum_edge)
                      << " MAX_BOUND_COMPONENTS " << limit[i].maximum_components
                      << " MAX_BOUND_MASS_NUM " << limit[i].maximum_mass << '\n';
        }
        std::cout << "DIRECT_MASS_CONTROLS q=50,6,366,367 PASS\n";

        Ledger ledger;
        u64 native_equalities5 = 0;
        u64 native_equalities6 = 0;
        for (std::size_t row = 0; row < EXPECTED_NATIVE.size(); ++row) {
            const int q = EXPECTED_NATIVE[row].q;
            const JointMass mass = build_joint_mass(fixed, support, q);
            const Layer native5 = build_layer(mass, support, 5);
            const Layer native6 = build_layer(mass, support, 6);
            require(native5.edges.size() == EXPECTED_NATIVE[row].e5,
                    "native E5 count changed");
            require(native6.edges.size() == EXPECTED_NATIVE[row].e6,
                    "native E6 count changed");
            native_equalities5 += native5.equalities;
            native_equalities6 += native6.equalities;

            const CoverAudit m4 = audit_covers(anchors[4], native6, 6);
            const CoverAudit m5 = audit_covers(anchors[5], native5, 5);
            const CoverAudit m6 = audit_covers(anchors[6], native5, 4);
            require(m4.rows_with_cover == 0 && m5.rows_with_cover == 0 &&
                        m6.rows_with_cover == 0,
                    "resonance retained a shallow native cover");

            if (q == 6) {
                require(!exhaustive_exact_cover(anchors[4].front(), native6, 6) &&
                            !exhaustive_exact_cover(anchors[5].front(), native5, 5) &&
                            !exhaustive_exact_cover(anchors[6].front(), native5, 4),
                        "q6 exhaustive cover cross-check disagreed");
            }

            const u64 affected4 = restricted_affected(
                anchors[4], records[q].inclusion[2].missing_edges);
            const u64 affected5 = restricted_affected(
                anchors[5], records[q].inclusion[2].missing_edges);
            const u64 affected6 = restricted_affected(
                anchors[6], records[q].inclusion[1].missing_edges);
            ledger.add(static_cast<u64>(q));
            ledger.add(records[q].inclusion[0].missing);
            ledger.add(records[q].inclusion[1].missing);
            ledger.add(records[q].inclusion[2].missing);
            ledger.add(native5.edges.size());
            ledger.add(native6.edges.size());
            ledger.add(affected4);
            ledger.add(affected5);
            ledger.add(affected6);
            ledger.add(m4.search_nodes);
            ledger.add(m5.search_nodes);
            ledger.add(m6.search_nodes);

            std::cout << "RESONANCE q " << q
                      << " LOST_E5_E6_E7 " << records[q].inclusion[0].missing
                      << ',' << records[q].inclusion[1].missing << ','
                      << records[q].inclusion[2].missing
                      << " AFFECTED_M4_M5_M6 " << affected4 << ',' << affected5
                      << ',' << affected6
                      << " NATIVE_E5_E6 " << native5.edges.size() << ','
                      << native6.edges.size()
                      << " NATIVE_EQUALITIES " << native5.equalities << ','
                      << native6.equalities
                      << " COVER_ROWS_M4_M5_M6 0,0,0"
                      << " SEARCH_NODES " << m4.search_nodes << ','
                      << m5.search_nodes << ',' << m6.search_nodes << '\n';
        }
        require(native_equalities5 == 0 && native_equalities6 == 0,
                "native threshold equality appeared");

        require(ledger.value() == EXPECTED_LEDGER,
                "independent semantic ledger changed");

        std::cout << "FINITE_UNIVERSE Q_LE " << FINITE_LIMIT
                  << " NEWCOMERS " << newcomer_count << '\n';
        std::cout << "FAILURE_Q5 " << failures5.size() << ' '
                  << integer_list(failures5) << '\n';
        std::cout << "FAILURE_Q6 " << failures6.size() << ' '
                  << integer_list(failures6) << '\n';
        std::cout << "FAILURE_Q7 " << failures7.size() << ' '
                  << integer_list(failures7) << '\n';
        std::cout << "INCLUSION_EQUALITIES E5_E6_E7 0,0,0\n";
        std::cout << "NATIVE_EQUALITIES E5_E6 0,0\n";
        std::cout << "Q6_EXHAUSTIVE_COVER_CONTROLS M4_M5_M6 0,0,0\n";
        std::cout << "INDEPENDENT_LEDGER_FNV1A64_LE " << hex64(ledger.value())
                  << '\n';
        std::cout << "COFINAL_INCLUSION q>=2587 BY_COMPONENT_DISCREPANCY\n";
        std::cout << "VERDICT PASS\n";
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
