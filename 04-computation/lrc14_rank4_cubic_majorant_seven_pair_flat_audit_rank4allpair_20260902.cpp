// Cleanroom seven-pair audit for the all-pair rank-four cubic majorant.
//
// This program imports no primary census code.  It constructs literal walls,
// uses an unfurled two-sided midpoint predicate, and flat-enumerates every one
// of C(30,8)=5,852,925 bodies for seven fixed universe indices.  It contains
// no optimizer pruning.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

constexpr std::array<int, 30> kPool = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
constexpr unsigned kBodySize = 8;
constexpr unsigned kTripleCube = 30 * 30 * 30;

struct Control {
    std::size_t universe_index;
    int q;
    int r;
    i128 expected_reward;
    u32 expected_body;
};

constexpr std::array<Control, 7> kControls = {{
    {0, 1, 2, i128{73888684844670LL}, UINT32_C(0x00bc2401)},
    {20910, 50, 70, i128{357589704156176LL}, UINT32_C(0x0138c402)},
    {60000, 147, 466, i128{111153401410931214LL}, UINT32_C(0x21740600)},
    {120000, 321, 612, i128{7111517651650272LL}, UINT32_C(0x10b82401)},
    {163336, 509, 640, i128{274811872704948582LL}, UINT32_C(0x10386401)},
    {181103, 721, 746, i128{2573965188596712114LL}, UINT32_C(0x21b02401)},
    {181193, 766, 768, i128{418404562977997824LL}, UINT32_C(0x05a82401)},
}};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}
void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

std::string decimal(i128 value) {
    if (value == 0) return "0";
    const bool negative = value < 0;
    u128 magnitude = static_cast<u128>(value);
    if (negative) magnitude = u128{0} - magnitude;
    std::string result;
    while (magnitude != 0) {
        result.push_back(static_cast<char>('0' + magnitude % 10));
        magnitude /= 10;
    }
    if (negative) result.push_back('-');
    std::reverse(result.begin(), result.end());
    return result;
}

std::string hex8(u32 mask) {
    std::ostringstream out;
    out << std::hex << std::setw(8) << std::setfill('0') << mask;
    return out.str();
}

i128 gcd128(i128 a, i128 b) {
    while (b != 0) {
        const i128 remainder = a % b;
        a = b;
        b = remainder;
    }
    return a;
}
i128 lcm128(i128 a, i128 b) { return a / gcd128(a, b) * b; }

unsigned triple_index(unsigned a, unsigned b, unsigned c) {
    return (a * 30 + b) * 30 + c;
}

bool midpoint_safe_unfurled(int speed, i128 grid, i128 left, i128 right) {
    i128 residue = i128(speed) * (left + right) % (2 * grid);
    if (residue < 0) residue += 2 * grid;
    return grid <= 7 * residue && 7 * residue <= 13 * grid;
}

struct Cell {
    u32 failures = 0;
    i128 width = 0;
};

struct Graph {
    i128 grid = 1;
    i128 retained = 0;
    std::array<i128, 5> rank_mass{};
    std::array<i128, 30> degree{};
    std::array<std::array<i128, 30>, 30> codegree{};
    std::vector<i128> triple = std::vector<i128>(kTripleCube, 0);
    std::vector<Cell> cells;
};

Graph literal_graph(int q, int r) {
    Graph graph;
    for (int speed : kPool) graph.grid = lcm128(graph.grid, 14 * i128(speed));
    graph.grid = lcm128(graph.grid, 14 * i128(q));
    graph.grid = lcm128(graph.grid, 14 * i128(r));

    std::vector<i128> walls = {0, graph.grid};
    auto append_walls = [&](int speed) {
        const i128 unit = graph.grid / (14 * i128(speed));
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14 * i128(tooth) + 1) * unit);
            walls.push_back((14 * i128(tooth) + 13) * unit);
        }
    };
    for (int speed : kPool) append_walls(speed);
    append_walls(q);
    append_walls(r);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());

    for (std::size_t cursor = 1; cursor < walls.size(); ++cursor) {
        const i128 left = walls[cursor - 1];
        const i128 right = walls[cursor];
        if (!midpoint_safe_unfurled(q, graph.grid, left, right) ||
            !midpoint_safe_unfurled(r, graph.grid, left, right))
            continue;
        u32 failures = 0;
        for (unsigned bit = 0; bit < 30; ++bit)
            if (!midpoint_safe_unfurled(
                    kPool[bit], graph.grid, left, right))
                failures |= u32{1} << bit;
        const unsigned rank = std::popcount(failures);
        if (rank > 4) continue;
        const i128 width = right - left;
        graph.retained += width;
        graph.rank_mass[rank] += width;
        graph.cells.push_back({failures, width});
        std::array<unsigned, 4> vertices{};
        unsigned count = 0;
        u32 copy = failures;
        while (copy != 0) {
            vertices[count++] = std::countr_zero(copy);
            copy &= copy - 1;
        }
        require(count == rank, "literal rank extraction changed");
        for (unsigned i = 0; i < count; ++i)
            graph.degree[vertices[i]] += width;
        for (unsigned i = 0; i < count; ++i)
            for (unsigned j = i + 1; j < count; ++j) {
                graph.codegree[vertices[i]][vertices[j]] += width;
                graph.codegree[vertices[j]][vertices[i]] += width;
            }
        for (unsigned i = 0; i < count; ++i)
            for (unsigned j = i + 1; j < count; ++j)
                for (unsigned k = j + 1; k < count; ++k) {
                    const unsigned a = vertices[i];
                    const unsigned b = vertices[j];
                    const unsigned c = vertices[k];
                    const std::array<std::array<unsigned, 3>, 6> perms = {
                        std::array<unsigned, 3>{a, b, c}, {a, c, b},
                        {b, a, c}, {b, c, a}, {c, a, b}, {c, b, a}};
                    for (const auto& p : perms)
                        graph.triple[triple_index(p[0], p[1], p[2])] += width;
                }
    }
    require(graph.retained == std::accumulate(graph.rank_mass.begin(),
                                              graph.rank_mass.end(), i128{0}),
            "literal retained partition changed");
    return graph;
}

i128 direct_reward(const Graph& graph, u32 body) {
    constexpr std::array<int, 5> g = {0, 14, 19, 18, 14};
    i128 reward = 0;
    for (const Cell& cell : graph.cells) {
        const unsigned hits = std::popcount(cell.failures & body);
        reward += i128(g[hits]) * cell.width;
    }
    return reward;
}

struct FlatEnumeration {
    const Graph& graph;
    i128 best_reward = -1;
    u32 best_body = 0;
    u64 nodes = 0;
    u64 leaves = 0;

    i128 marginal(unsigned vertex,
                  const std::array<unsigned, kBodySize>& selected,
                  unsigned count) const {
        i128 gain = 14 * graph.degree[vertex];
        for (unsigned i = 0; i < count; ++i)
            gain -= 9 * graph.codegree[vertex][selected[i]];
        for (unsigned i = 0; i < count; ++i)
            for (unsigned j = i + 1; j < count; ++j)
                gain += 3 * graph.triple[triple_index(
                                vertex, selected[i], selected[j])];
        return gain;
    }

    void visit(unsigned start, unsigned count,
               std::array<unsigned, kBodySize>& selected,
               i128 reward, u32 body) {
        ++nodes;
        if (count == kBodySize) {
            ++leaves;
            if (reward > best_reward ||
                (reward == best_reward && body < best_body)) {
                best_reward = reward;
                best_body = body;
            }
            return;
        }
        const unsigned need = kBodySize - count;
        for (unsigned vertex = start; vertex <= 30 - need; ++vertex) {
            const i128 gain = marginal(vertex, selected, count);
            selected[count] = vertex;
            visit(vertex + 1, count + 1, selected, reward + gain,
                  body | (u32{1} << vertex));
        }
    }
};

std::vector<std::pair<int, int>> read_universe(const std::string& path) {
    std::ifstream input(path);
    require(bool(input), "cannot open universe");
    std::vector<std::pair<int, int>> result;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos, "malformed universe row");
        result.emplace_back(std::stoi(line.substr(0, comma)),
                            std::stoi(line.substr(comma + 1)));
    }
    require(result.size() == 181194, "universe cardinality changed");
    return result;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "usage: seven_pair_flat_audit PAIR_UNIVERSE");
        const auto universe = read_universe(argv[1]);
        std::cout << "LRC14_RANK4_CUBIC_MAJORANT_SEVEN_PAIR_FLAT_AUDIT_V1\n"
                  << "METHOD LITERAL_WALLS_UNFURLED_TWO_SIDED_MIDPOINT_"
                     "FLAT_C30_8_NO_PRUNING\n";
        for (const Control& control : kControls) {
            require(universe[control.universe_index] ==
                        std::pair{control.q, control.r},
                    "control universe index changed");
            const Graph graph = literal_graph(control.q, control.r);
            FlatEnumeration flat{graph};
            std::array<unsigned, kBodySize> selected{};
            flat.visit(0, 0, selected, 0, 0);
            require(flat.leaves == UINT64_C(5852925),
                    "flat body cardinality changed");
            require(flat.nodes == UINT64_C(7888725),
                    "flat recursion node count changed");
            require(flat.best_reward == control.expected_reward,
                    "flat optimum reward changed");
            require(flat.best_body == control.expected_body,
                    "flat least maximizing body changed");
            require(direct_reward(graph, flat.best_body) == flat.best_reward,
                    "flat tensor/direct reward mismatch");
            const i128 lower14 = 14 * graph.retained - flat.best_reward;
            const i128 ticks = 81 * lower14 - 98 * graph.grid;
            require(ticks > 0, "flat control misses 7/81 target");
            std::cout << "INDEX " << control.universe_index << " PAIR "
                      << control.q << ',' << control.r << " GRID "
                      << decimal(graph.grid) << " RETAINED "
                      << decimal(graph.retained) << " REWARD "
                      << decimal(flat.best_reward) << " LEAST_BODY "
                      << hex8(flat.best_body) << " LOWER14 "
                      << decimal(lower14) << " TICKS " << decimal(ticks)
                      << " LEAVES " << flat.leaves << " STATUS MATCH\n";
        }
        std::cout << "SCOPE SEVEN_STRATIFIED_PAIR_CONTROLS_NOT_GLOBAL_CENSUS\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "RANK4_SEVEN_PAIR_AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}
