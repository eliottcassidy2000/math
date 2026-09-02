// Independent hostile audit: literal midpoint walls plus flat C(30,9)
// enumeration.  This does not import the event sweep or branch-and-bound.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
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
constexpr std::array<int, 30> pool = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
constexpr unsigned body_size = 9;

[[noreturn]] void die(const std::string& message) {
    throw std::runtime_error(message);
}
void check(bool value, const std::string& message) {
    if (!value) die(message);
}
std::string dec(i128 value) {
    if (value == 0) return "0";
    const bool neg = value < 0;
    u128 mag = static_cast<u128>(value);
    if (neg) mag = u128{0} - mag;
    std::string out;
    while (mag) {
        out.push_back(char('0' + mag % 10));
        mag /= 10;
    }
    if (neg) out.push_back('-');
    std::reverse(out.begin(), out.end());
    return out;
}
std::string hex8(u32 mask) {
    std::ostringstream out;
    out << std::hex << std::setw(8) << std::setfill('0') << mask;
    return out.str();
}
i128 gcd128(i128 a, i128 b) {
    while (b) {
        const i128 r = a % b;
        a = b;
        b = r;
    }
    return a;
}
i128 lcm128(i128 a, i128 b) { return a / gcd128(a, b) * b; }
bool safe(int speed, i128 grid, i128 left, i128 right) {
    i128 residue = i128(speed) * (left + right) % (2 * grid);
    if (residue < 0) residue += 2 * grid;
    return 7 * residue >= grid && 7 * residue <= 13 * grid;
}
struct Graph {
    i128 grid = 1, zero = 0, total = 0;
    std::array<i128, 30> one{}, degree{};
    std::array<std::array<i128, 30>, 30> two{}, codegree{};
    std::array<std::array<std::array<i128, 30>, 30>, 30> three{};
    std::array<u64, 4> cells{};
};

Graph build() {
    constexpr int q = 50, r = 70;
    Graph graph;
    for (int speed : pool) graph.grid = lcm128(graph.grid, 14 * i128(speed));
    graph.grid = lcm128(graph.grid, 14 * i128(q));
    graph.grid = lcm128(graph.grid, 14 * i128(r));
    std::vector<i128> walls = {0, graph.grid};
    auto add = [&](int speed) {
        const i128 unit = graph.grid / (14 * i128(speed));
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14 * i128(tooth) + 1) * unit);
            walls.push_back((14 * i128(tooth) + 13) * unit);
        }
    };
    for (int speed : pool) add(speed);
    add(q);
    add(r);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i128 left = walls[index - 1], right = walls[index];
        if (!safe(q, graph.grid, left, right) ||
            !safe(r, graph.grid, left, right))
            continue;
        u32 failure = 0;
        for (unsigned bit = 0; bit < 30; ++bit)
            if (!safe(pool[bit], graph.grid, left, right))
                failure |= u32{1} << bit;
        const unsigned rank = std::popcount(failure);
        if (rank > 3) continue;
        const i128 width = right - left;
        ++graph.cells[rank];
        if (rank == 0) {
            graph.zero += width;
        } else if (rank == 1) {
            graph.one[std::countr_zero(failure)] += width;
        } else if (rank == 2) {
            const unsigned a = std::countr_zero(failure);
            const unsigned b = std::countr_zero(failure & (failure - 1));
            graph.two[a][b] += width;
            graph.two[b][a] += width;
        } else {
            u32 copy = failure;
            const unsigned a = std::countr_zero(copy);
            copy &= copy - 1;
            const unsigned b = std::countr_zero(copy);
            copy &= copy - 1;
            const unsigned c = std::countr_zero(copy);
            for (const auto permutation : std::array<std::array<unsigned, 3>, 6>{
                     std::array<unsigned, 3>{a, b, c}, {a, c, b}, {b, a, c},
                     {b, c, a}, {c, a, b}, {c, b, a}})
                graph.three[permutation[0]][permutation[1]][permutation[2]] +=
                    width;
        }
    }
    graph.total = graph.zero;
    for (unsigned a = 0; a < 30; ++a) {
        graph.total += graph.one[a];
        graph.degree[a] += graph.one[a];
    }
    for (unsigned a = 0; a < 30; ++a)
        for (unsigned b = a + 1; b < 30; ++b) {
            const i128 weight = graph.two[a][b];
            graph.total += weight;
            graph.degree[a] += weight;
            graph.degree[b] += weight;
            graph.codegree[a][b] += weight;
            graph.codegree[b][a] += weight;
        }
    for (unsigned a = 0; a < 30; ++a)
        for (unsigned b = a + 1; b < 30; ++b)
            for (unsigned c = b + 1; c < 30; ++c) {
                const i128 weight = graph.three[a][b][c];
                graph.total += weight;
                graph.degree[a] += weight;
                graph.degree[b] += weight;
                graph.degree[c] += weight;
                for (const auto pair : std::array<std::pair<unsigned, unsigned>, 6>{
                         std::pair{a, b}, {b, a}, {a, c}, {c, a}, {b, c}, {c, b}})
                    graph.codegree[pair.first][pair.second] += weight;
            }
    return graph;
}

struct Flat {
    const Graph& graph;
    i128 best_reward = -1;
    u32 best_body = 0;
    u64 leaves = 0, nodes = 0;

    i128 marginal(unsigned vertex,
                  const std::array<unsigned, body_size>& selected,
                  unsigned count) const {
        i128 gain = graph.degree[vertex];
        for (unsigned i = 0; i < count; ++i)
            gain -= graph.codegree[vertex][selected[i]];
        for (unsigned i = 0; i < count; ++i)
            for (unsigned j = i + 1; j < count; ++j)
                gain += graph.three[vertex][selected[i]][selected[j]];
        return gain;
    }

    void visit(unsigned start, unsigned count,
               std::array<unsigned, body_size>& selected,
               i128 reward, u32 body) {
        ++nodes;
        if (count == body_size) {
            ++leaves;
            if (reward > best_reward ||
                (reward == best_reward && body < best_body)) {
                best_reward = reward;
                best_body = body;
            }
            return;
        }
        const unsigned need = body_size - count;
        for (unsigned vertex = start; vertex <= 30 - need; ++vertex) {
            const i128 gain = marginal(vertex, selected, count);
            selected[count] = vertex;
            visit(vertex + 1, count + 1, selected, reward + gain,
                  body | (u32{1} << vertex));
        }
    }
};

i128 direct_survivor(const Graph& graph, u32 body) {
    i128 mass = graph.zero;
    for (unsigned a = 0; a < 30; ++a)
        if (!(body & (u32{1} << a))) mass += graph.one[a];
    for (unsigned a = 0; a < 30; ++a)
        if (!(body & (u32{1} << a)))
            for (unsigned b = a + 1; b < 30; ++b)
                if (!(body & (u32{1} << b))) mass += graph.two[a][b];
    for (unsigned a = 0; a < 30; ++a)
        if (!(body & (u32{1} << a)))
            for (unsigned b = a + 1; b < 30; ++b)
                if (!(body & (u32{1} << b)))
                    for (unsigned c = b + 1; c < 30; ++c)
                        if (!(body & (u32{1} << c)))
                            mass += graph.three[a][b][c];
    return mass;
}
}  // namespace

int main() {
    try {
        Graph graph = build();
        Flat flat{graph};
        std::array<unsigned, body_size> selected{};
        flat.visit(0, 0, selected, 0, 0);
        const i128 minimum = graph.total - flat.best_reward;
        const i128 direct = direct_survivor(graph, flat.best_body);
        const i128 ticks = 27 * minimum - 2 * graph.grid;
        check(flat.leaves == UINT64_C(14307150), "body count changed");
        check(direct == minimum, "direct survivor disagrees");
        check(graph.grid == i128(91205797082400LL), "grid changed");
        check(graph.total == i128(22084503304648LL), "rank3 total changed");
        check(minimum == i128(7799459654598LL), "minimum changed");
        check(ticks == i128(28173816509346LL), "ticks changed");
        check(flat.best_body == UINT32_C(0x011cd402), "least body changed");
        std::cout << "INDEPENDENT_PAIR5070_FLAT_AUDIT_V1\n"
                  << "GRID " << dec(graph.grid) << " CELLS " << graph.cells[0]
                  << ',' << graph.cells[1] << ',' << graph.cells[2] << ','
                  << graph.cells[3] << " RANK3_TOTAL " << dec(graph.total)
                  << '\n'
                  << "BODIES " << flat.leaves << " NODES " << flat.nodes
                  << " MIN_MASS " << dec(minimum) << " TICKS " << dec(ticks)
                  << " LEAST_BODY " << hex8(flat.best_body) << '\n'
                  << "METHOD LITERAL_MIDPOINT_WALLS_PLUS_UNPRUNED_COMBINATION_ENUMERATION\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "INDEPENDENT_FLAT_ERROR " << error.what() << '\n';
        return 1;
    }
}

