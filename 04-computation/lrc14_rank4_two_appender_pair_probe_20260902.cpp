// Exact rank-four retained-mass probe for the two THM-4333 control pairs.
//
// For every eight-label body in the fixed thirty-label pool, this program
// computes the mass retained from pair-safe wall cells of pool-failure rank
// at most four.  It flat-enumerates all C(30,8) bodies, with no optimizer
// pruning, and tests the two-appender target L_4/D > 7/81.  The complete
// safe mass is reported separately on the rank-four minimizer; it is not
// used in the retained-truncation test.

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

constexpr std::array<int, 30> kPool = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
constexpr unsigned kBodySize = 8;
constexpr unsigned kN = 30;
constexpr unsigned kCube3 = kN * kN * kN;
constexpr unsigned kCube4 = kCube3 * kN;

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

std::string body_labels(u32 mask) {
    std::ostringstream out;
    bool first = true;
    for (unsigned bit = 0; bit < kN; ++bit) {
        if (((mask >> bit) & 1U) == 0) continue;
        if (!first) out << ',';
        first = false;
        out << kPool[bit];
    }
    return out.str();
}

i128 gcd128(i128 a, i128 b);

std::string reduced_fraction(i128 numerator, i128 denominator) {
    require(denominator > 0, "fraction denominator is not positive");
    const i128 divisor = gcd128(numerator < 0 ? -numerator : numerator,
                                denominator);
    return decimal(numerator / divisor) + "/" + decimal(denominator / divisor);
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

unsigned index3(unsigned a, unsigned b, unsigned c) {
    return (a * kN + b) * kN + c;
}
unsigned index4(unsigned a, unsigned b, unsigned c, unsigned d) {
    return ((a * kN + b) * kN + c) * kN + d;
}

// Evaluate at the midpoint of an open wall cell without fractions.  The
// residue is speed*(left+right) modulo 2D, so distance >= 1/14 becomes
// D <= 7*residue <= 13D after folding to the [0,D] half-circle.
bool midpoint_safe(int speed, i128 grid, i128 left, i128 right) {
    i128 residue = i128(speed) * (left + right) % (2 * grid);
    if (residue < 0) residue += 2 * grid;
    if (residue > grid) residue = 2 * grid - residue;
    return 7 * residue >= grid;
}

struct Cell {
    u32 failures = 0;
    i128 width = 0;
};

struct Graph4 {
    i128 grid = 1;
    i128 retained = 0;
    i128 pair_safe_mass = 0;
    i128 rank0 = 0;
    std::array<i128, kN> degree{};
    std::array<std::array<i128, kN>, kN> codegree{};
    std::vector<i128> triple_degree = std::vector<i128>(kCube3, 0);
    std::vector<i128> rank4 = std::vector<i128>(kCube4, 0);
    std::array<u64, kN + 1> cell_count{};
    std::array<i128, kN + 1> rank_mass{};
    std::vector<Cell> cells;
};

i128 quadruple_weight(const Graph4& graph, unsigned a, unsigned b,
                     unsigned c, unsigned d) {
    std::array<unsigned, 4> vertices = {a, b, c, d};
    std::sort(vertices.begin(), vertices.end());
    return graph.rank4[index4(vertices[0], vertices[1], vertices[2],
                              vertices[3])];
}

i128 triple_degree(const Graph4& graph, unsigned a, unsigned b, unsigned c) {
    return graph.triple_degree[index3(a, b, c)];
}

Graph4 build_graph(int q, int r) {
    Graph4 graph;
    for (int speed : kPool)
        graph.grid = lcm128(graph.grid, 14 * i128(speed));
    graph.grid = lcm128(graph.grid, 14 * i128(q));
    graph.grid = lcm128(graph.grid, 14 * i128(r));

    std::vector<i128> walls = {0, graph.grid};
    auto add_walls = [&](int speed) {
        const i128 unit = graph.grid / (14 * i128(speed));
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14 * i128(tooth) + 1) * unit);
            walls.push_back((14 * i128(tooth) + 13) * unit);
        }
    };
    for (int speed : kPool) add_walls(speed);
    add_walls(q);
    add_walls(r);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());

    for (std::size_t wall = 1; wall < walls.size(); ++wall) {
        const i128 left = walls[wall - 1];
        const i128 right = walls[wall];
        if (!midpoint_safe(q, graph.grid, left, right) ||
            !midpoint_safe(r, graph.grid, left, right))
            continue;
        u32 failures = 0;
        for (unsigned bit = 0; bit < kN; ++bit) {
            if (!midpoint_safe(kPool[bit], graph.grid, left, right))
                failures |= u32{1} << bit;
        }
        const i128 width = right - left;
        const unsigned rank = std::popcount(failures);
        ++graph.cell_count[rank];
        graph.rank_mass[rank] += width;
        graph.pair_safe_mass += width;
        graph.cells.push_back({failures, width});
        if (rank > 4) continue;

        graph.retained += width;
        if (rank == 0) {
            graph.rank0 += width;
            continue;
        }
        std::array<unsigned, 4> vertices{};
        unsigned count = 0;
        u32 copy = failures;
        while (copy != 0) {
            vertices[count++] = std::countr_zero(copy);
            copy &= copy - 1;
        }
        require(count == rank, "failure-rank extraction changed");
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
                    const std::array<std::array<unsigned, 3>, 6> permutations = {
                        std::array<unsigned, 3>{vertices[i], vertices[j], vertices[k]},
                        {vertices[i], vertices[k], vertices[j]},
                        {vertices[j], vertices[i], vertices[k]},
                        {vertices[j], vertices[k], vertices[i]},
                        {vertices[k], vertices[i], vertices[j]},
                        {vertices[k], vertices[j], vertices[i]}};
                    for (const auto& p : permutations)
                        graph.triple_degree[index3(p[0], p[1], p[2])] += width;
                }
        if (rank == 4) {
            graph.rank4[index4(vertices[0], vertices[1], vertices[2],
                                vertices[3])] += width;
        }
    }
    require(graph.retained == graph.rank_mass[0] + graph.rank_mass[1] +
                                  graph.rank_mass[2] + graph.rank_mass[3] +
                                  graph.rank_mass[4],
            "retained-rank mass identity failed");
    i128 all_rank_mass = 0;
    for (i128 mass : graph.rank_mass) all_rank_mass += mass;
    require(all_rank_mass == graph.pair_safe_mass,
            "pair-safe rank partition failed");
    return graph;
}

struct FlatMinimum {
    const Graph4& graph;
    i128 best_reward = -1;
    u32 best_body = 0;
    u64 leaves = 0;
    u64 nodes = 0;

    i128 marginal(unsigned vertex,
                  const std::array<unsigned, kBodySize>& selected,
                  unsigned count) const {
        i128 gain = graph.degree[vertex];
        for (unsigned i = 0; i < count; ++i)
            gain -= graph.codegree[vertex][selected[i]];
        for (unsigned i = 0; i < count; ++i)
            for (unsigned j = i + 1; j < count; ++j)
                gain += triple_degree(graph, vertex, selected[i], selected[j]);
        for (unsigned i = 0; i < count; ++i)
            for (unsigned j = i + 1; j < count; ++j)
                for (unsigned k = j + 1; k < count; ++k)
                    gain -= quadruple_weight(graph, vertex, selected[i],
                                             selected[j], selected[k]);
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
        for (unsigned vertex = start; vertex <= kN - need; ++vertex) {
            const i128 gain = marginal(vertex, selected, count);
            selected[count] = vertex;
            visit(vertex + 1, count + 1, selected, reward + gain,
                  body | (u32{1} << vertex));
        }
    }
};

// Independent exact optimizer on the literal wall-cell coverage objective.
// Unlike FlatMinimum, this path does not use inclusion-exclusion tensors:
// a marginal is the direct weight of still-uncovered cells incident to v.
struct DirectCoverageBnB {
    const Graph4& graph;
    std::array<std::vector<unsigned>, kN> incident{};
    std::array<unsigned, kN> order{};
    i128 best_reward = -1;
    u32 best_body = 0;
    u64 nodes = 0;
    u64 prunes = 0;

    explicit DirectCoverageBnB(const Graph4& source) : graph(source) {
        for (unsigned cell = 0; cell < graph.cells.size(); ++cell) {
            const u32 failures = graph.cells[cell].failures;
            const unsigned rank = std::popcount(failures);
            if (rank == 0 || rank > 4) continue;
            for (unsigned vertex = 0; vertex < kN; ++vertex)
                if ((failures >> vertex) & 1U)
                    incident[vertex].push_back(cell);
        }
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](unsigned a, unsigned b) {
            if (graph.degree[a] != graph.degree[b])
                return graph.degree[a] > graph.degree[b];
            return a < b;
        });
        greedy_seed();
    }

    i128 marginal(unsigned vertex, u32 selected) const {
        i128 gain = 0;
        for (unsigned cell : incident[vertex])
            if ((graph.cells[cell].failures & selected) == 0)
                gain += graph.cells[cell].width;
        return gain;
    }

    void consider(u32 body, i128 reward) {
        require(std::popcount(body) == static_cast<int>(kBodySize),
                "direct optimizer body rank changed");
        if (reward > best_reward ||
            (reward == best_reward && body < best_body)) {
            best_reward = reward;
            best_body = body;
        }
    }

    void greedy_seed() {
        u32 selected = 0;
        i128 reward = 0;
        for (unsigned step = 0; step < kBodySize; ++step) {
            unsigned choice = kN;
            i128 gain = -1;
            for (unsigned vertex = 0; vertex < kN; ++vertex) {
                if ((selected >> vertex) & 1U) continue;
                const i128 candidate = marginal(vertex, selected);
                if (candidate > gain ||
                    (candidate == gain && vertex < choice)) {
                    choice = vertex;
                    gain = candidate;
                }
            }
            require(choice < kN, "direct greedy exhausted vertices");
            selected |= u32{1} << choice;
            reward += gain;
        }
        consider(selected, reward);
    }

    void descend(unsigned position, unsigned need, u32 selected, i128 reward) {
        ++nodes;
        if (need == 0) {
            consider(selected, reward);
            return;
        }
        if (kN - position < need) return;

        std::array<i128, kN> future{};
        unsigned future_count = 0;
        for (unsigned index = position; index < kN; ++index)
            future[future_count++] = marginal(order[index], selected);
        std::sort(future.begin(), future.begin() + future_count,
                  std::greater<>());
        i128 upper = reward;
        for (unsigned i = 0; i < need; ++i) upper += future[i];
        if (upper < best_reward) {
            ++prunes;
            return;
        }

        const unsigned vertex = order[position];
        const i128 gain = marginal(vertex, selected);
        descend(position + 1, need - 1, selected | (u32{1} << vertex),
                reward + gain);
        descend(position + 1, need, selected, reward);
    }

    void solve() { descend(0, kBodySize, 0, 0); }
};

i128 direct_mass(const Graph4& graph, u32 body, unsigned maximum_rank) {
    i128 total = 0;
    for (const Cell& cell : graph.cells) {
        if (static_cast<unsigned>(std::popcount(cell.failures)) <= maximum_rank &&
            (cell.failures & body) == 0)
            total += cell.width;
    }
    return total;
}

void run_pair(int q, int r) {
    Graph4 graph = build_graph(q, r);
    const i128 inherited_rank3 = graph.rank_mass[0] + graph.rank_mass[1] +
                                 graph.rank_mass[2] + graph.rank_mass[3];
    const i128 expected_rank3 =
        q == 50 && r == 70
            ? i128(22084503304648LL)
            : i128(21053923224812998LL);
    require(inherited_rank3 == expected_rank3,
            "rank-three prefix disagrees with THM-4333 packet");
    FlatMinimum flat{graph};
    std::array<unsigned, kBodySize> selected{};
    flat.visit(0, 0, selected, 0, 0);
    require(flat.leaves == UINT64_C(5852925), "body count changed");
    DirectCoverageBnB direct_optimizer(graph);
    direct_optimizer.solve();
    require(direct_optimizer.best_reward == flat.best_reward,
            "direct optimizer reward disagrees with flat inclusion-exclusion");
    require(direct_optimizer.best_body == flat.best_body,
            "direct optimizer least body disagrees with flat enumeration");
    const i128 minimum = graph.retained - flat.best_reward;
    const i128 direct_rank4 = direct_mass(graph, flat.best_body, 4);
    const i128 candidate_full = direct_mass(graph, flat.best_body, kN);
    require(direct_rank4 == minimum, "direct rank-four mass disagrees");
    require(candidate_full >= minimum,
            "complete safe mass below retained truncation");
    const i128 target_ticks = 81 * minimum - 7 * graph.grid;
    const i128 candidate_full_ticks = 81 * candidate_full - 7 * graph.grid;
    std::array<i128, kN> sorted_degrees = graph.degree;
    std::sort(sorted_degrees.begin(), sorted_degrees.end(), std::greater<>());
    const i128 top_eight_degree_sum = std::accumulate(
        sorted_degrees.begin(), sorted_degrees.begin() + kBodySize, i128{0});
    const i128 degree_bound = graph.retained - top_eight_degree_sum;
    const i128 degree_bound_ticks = 81 * degree_bound - 7 * graph.grid;

    std::cout << "PAIR " << q << ',' << r << " GRID " << decimal(graph.grid)
              << " PAIR_SAFE_MASS " << decimal(graph.pair_safe_mass)
              << " RETAINED_RANK4_MASS " << decimal(graph.retained) << '\n';
    std::cout << "RANK_CELLS";
    for (unsigned rank = 0; rank <= kN; ++rank)
        if (graph.cell_count[rank] != 0)
            std::cout << ' ' << rank << ':' << graph.cell_count[rank];
    std::cout << '\n';
    std::cout << "RANK_MASSES_0_TO_4";
    for (unsigned rank = 0; rank <= 4; ++rank)
        std::cout << ' ' << rank << ':' << decimal(graph.rank_mass[rank]);
    std::cout << '\n';
    std::cout << "THM4333_RANK3_PREFIX " << decimal(inherited_rank3)
              << " STATUS MATCH\n";
    std::cout << "BODIES " << flat.leaves << " NODES " << flat.nodes
              << " MIN_RETAINED_L4 " << decimal(minimum)
              << " MIN_RETAINED_RATIO "
              << reduced_fraction(minimum, graph.grid)
              << " TARGET_TICKS_81L4_MINUS_7D " << decimal(target_ticks)
              << " SURPLUS_OVER_7_81 "
              << reduced_fraction(target_ticks, 81 * graph.grid)
              << " LEAST_BODY " << hex8(flat.best_body)
              << " LABELS " << body_labels(flat.best_body) << '\n';
    std::cout << "ROOT_DEGREE_BOUND " << decimal(degree_bound)
              << " DEGREE_BOUND_TARGET_TICKS " << decimal(degree_bound_ticks)
              << " DEGREE_SCREEN_STATUS "
              << (degree_bound_ticks > 0 ? "STRICT_PASS" : "INCONCLUSIVE")
              << '\n';
    std::cout << "INDEPENDENT_DIRECT_COVERAGE_BNB NODES "
              << direct_optimizer.nodes << " PRUNES " << direct_optimizer.prunes
              << " OPTIMUM_REWARD " << decimal(direct_optimizer.best_reward)
              << " LEAST_BODY " << hex8(direct_optimizer.best_body)
              << " STATUS MATCH\n";
    std::cout << "SAME_BODY_COMPLETE_SAFE_MASS " << decimal(candidate_full)
              << " COMPLETE_TARGET_TICKS " << decimal(candidate_full_ticks)
              << " OMITTED_SAFE_MASS " << decimal(candidate_full - minimum)
              << '\n';
    std::cout << "TARGET_STATUS " << (target_ticks > 0 ? "STRICT_PASS" : "FAIL")
              << " FULL_MASS_NOTE SAME_RANK4_MINIMIZER_NOT_GLOBAL_FULL_MINIMUM\n";
}
}  // namespace

int main() {
    try {
        std::cout << "LRC14_RANK4_TWO_APPENDER_PAIR_PROBE_V1\n"
                  << "TARGET FOR_ALL_K_IN_C(P,8): 81*L4(K)-7*D>0\n"
                  << "IDENTITY (6/7)*(7/81)=2/27; (6/7)*(2/27)=4/63\n"
                  << "METHOD LITERAL_MIDPOINT_WALLS_PLUS_UNPRUNED_C(30,8)_ENUMERATION\n";
        run_pair(50, 70);
        run_pair(509, 640);
        std::cout << "SCOPE TWO_THM4333_RESIDUAL_PAIRS_ONLY_RETAINED_RANK4_NOT_FULL_MASS_MINIMUM_NO_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "RANK4_PAIR_PROBE_ERROR " << error.what() << '\n';
        return 1;
    }
}
