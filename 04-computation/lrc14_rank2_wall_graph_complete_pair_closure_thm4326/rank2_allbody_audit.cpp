// Scratch exact all-C(30,9) optimizer for the rank-at-most-two wall skeleton.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
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
using i64 = i128;

constexpr std::array<int, 30> kPool = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
constexpr u64 kBodyCount = UINT64_C(14307150);

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
    std::string answer;
    while (magnitude != 0) {
        answer.push_back(static_cast<char>('0' + magnitude % 10));
        magnitude /= 10;
    }
    if (negative) answer.push_back('-');
    std::reverse(answer.begin(), answer.end());
    return answer;
}

std::string hex8(u32 value) {
    std::ostringstream output;
    output << std::hex << std::setw(8) << std::setfill('0') << value;
    return output.str();
}

 i64 wide_gcd(i64 left, i64 right) {
    while (right != 0) {
        const i64 remainder = left % right;
        left = right;
        right = remainder;
    }
    return left;
}

i64 checked_lcm(i64 left, i64 right) {
    const i64 divisor = wide_gcd(left, right);
    const i128 result = static_cast<i128>(left / divisor) * right;
    require(result > 0, "LCM overflow");
    return result;
}

i64 base_grid() {
    i64 grid = 1;
    for (int speed : kPool) grid = checked_lcm(grid, 14LL * speed);
    require(grid == INT64_C(18241159416480), "base grid changed");
    return grid;
}

bool safe_midpoint(int speed, i64 grid, i64 left, i64 right) {
    i128 residue = static_cast<i128>(speed) *
                   (static_cast<i128>(left) + right);
    residue %= static_cast<i128>(2) * grid;
    if (residue < 0) residue += static_cast<i128>(2) * grid;
    return static_cast<i128>(7) * residue >= grid &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * grid;
}

u32 next_combination(u32 value) {
    const u32 smallest = value & (u32{0} - value);
    const u32 ripple = value + smallest;
    return ripple | (((value ^ ripple) >> 2) / smallest);
}

struct Graph {
    i64 grid = 0;
    i64 zero = 0;
    std::array<i64, 30> singleton{};
    std::array<std::array<i64, 30>, 30> edge{};
    std::array<i64, 30> incident{};
    i64 total = 0;
    u64 rank0_cells = 0;
    u64 rank1_cells = 0;
    u64 rank2_cells = 0;
};

Graph build_graph(int q, int r) {
    Graph graph;
    graph.grid = checked_lcm(base_grid(), 14LL * q);
    graph.grid = checked_lcm(graph.grid, 14LL * r);
    std::vector<i64> walls = {0, graph.grid};
    auto add_walls = [&](int speed) {
        const i64 unit = graph.grid / (14LL * speed);
        require(unit * 14LL * speed == graph.grid, "nonintegral unit");
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int speed : kPool) add_walls(speed);
    add_walls(q);
    add_walls(r);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1], right = walls[index];
        if (!safe_midpoint(q, graph.grid, left, right) ||
            !safe_midpoint(r, graph.grid, left, right))
            continue;
        u32 failure = 0;
        for (unsigned bit = 0; bit < kPool.size(); ++bit)
            if (!safe_midpoint(kPool[bit], graph.grid, left, right))
                failure |= u32{1} << bit;
        const unsigned rank = std::popcount(failure);
        const i64 width = right - left;
        if (rank == 0) {
            graph.zero += width;
            ++graph.rank0_cells;
        } else if (rank == 1) {
            const unsigned vertex = std::countr_zero(failure);
            graph.singleton[vertex] += width;
            ++graph.rank1_cells;
        } else if (rank == 2) {
            const unsigned left_vertex = std::countr_zero(failure);
            const unsigned right_vertex = std::countr_zero(
                failure & (failure - 1));
            graph.edge[left_vertex][right_vertex] += width;
            graph.edge[right_vertex][left_vertex] += width;
            ++graph.rank2_cells;
        }
    }
    graph.total = graph.zero;
    for (unsigned vertex = 0; vertex < 30; ++vertex) {
        graph.total += graph.singleton[vertex];
        graph.incident[vertex] += graph.singleton[vertex];
    }
    for (unsigned left = 0; left < 30; ++left)
        for (unsigned right = left + 1; right < 30; ++right) {
            graph.total += graph.edge[left][right];
            graph.incident[left] += graph.edge[left][right];
            graph.incident[right] += graph.edge[left][right];
        }
    return graph;
}

struct Audit {
    u64 bodies = 0;
    u64 positive = 0;
    u64 zero = 0;
    i128 minimum = 0;
    i128 maximum = 0;
    u32 minimum_body = 0;
    u32 maximum_body = 0;
};

i64 degree_lower_bound(const Graph& graph) {
    std::array<i64, 30> ordered = graph.incident;
    std::sort(ordered.begin(), ordered.end(), std::greater<i64>());
    i64 lower = graph.total;
    for (unsigned index = 0; index < 9; ++index) lower -= ordered[index];
    return lower;
}

Audit audit(const Graph& graph) {
    Audit audit;
    const u32 limit = u32{1} << 30;
    for (u32 body = (u32{1} << 9) - 1; body < limit;
         body = next_combination(body)) {
        i64 mass = graph.total;
        std::array<unsigned, 9> vertices{};
        unsigned count = 0;
        u32 copy = body;
        while (copy != 0) {
            const unsigned vertex = std::countr_zero(copy);
            vertices[count++] = vertex;
            mass -= graph.incident[vertex];
            copy &= copy - 1;
        }
        require(count == 9, "body rank changed");
        for (unsigned left = 0; left < 9; ++left)
            for (unsigned right = left + 1; right < 9; ++right)
                mass += graph.edge[vertices[left]][vertices[right]];
        const i128 ticks = static_cast<i128>(63) * mass -
                           static_cast<i128>(4) * graph.grid;
        if (audit.bodies == 0 || ticks < audit.minimum ||
            (ticks == audit.minimum && body < audit.minimum_body)) {
            audit.minimum = ticks;
            audit.minimum_body = body;
        }
        if (audit.bodies == 0 || ticks > audit.maximum ||
            (ticks == audit.maximum && body < audit.maximum_body)) {
            audit.maximum = ticks;
            audit.maximum_body = body;
        }
        audit.positive += ticks > 0;
        audit.zero += ticks == 0;
        ++audit.bodies;
    }
    require(audit.bodies == kBodyCount, "body universe changed");
    return audit;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc >= 2, "usage: rank2_allbody q,r...");
        std::vector<std::string> pair_tokens;
        if (argc == 3 && std::string(argv[1]) == "--degree-csv") {
            std::ifstream input(argv[2]);
            require(static_cast<bool>(input), "cannot open degree CSV");
            std::string line;
            require(std::getline(input, line) &&
                        line == "q,r,grid,rank2_total,degree_bound_mass,degree_bound_ticks,positive,top9_hex",
                    "degree CSV header changed");
            while (std::getline(input, line)) {
                if (line.empty()) continue;
                std::vector<std::string> fields;
                std::size_t start = 0;
                while (true) {
                    const std::size_t comma = line.find(',', start);
                    fields.push_back(line.substr(start, comma - start));
                    if (comma == std::string::npos) break;
                    start = comma + 1;
                }
                require(fields.size() == 8, "malformed degree CSV row");
                if (fields[6] == "0")
                    pair_tokens.push_back(fields[0] + "," + fields[1]);
                else
                    require(fields[6] == "1", "invalid degree positive flag");
            }
            require(input.eof() && !pair_tokens.empty(),
                    "no coarse candidates in degree CSV");
        } else {
            for (int index = 1; index < argc; ++index)
                pair_tokens.push_back(argv[index]);
        }
        std::cout << "LRC14_RANK2_ALLBODY_GRAPH_AUDIT_V1\n";
        for (const std::string& token : pair_tokens) {
            const std::size_t comma = token.find(',');
            require(comma != std::string::npos, "pair token missing comma");
            const int q = std::stoi(token.substr(0, comma));
            const int r = std::stoi(token.substr(comma + 1));
            require(q > 0 && q < r, "pair ordering failed");
            const Graph graph = build_graph(q, r);
            const Audit result = audit(graph);
            const i64 coarse_mass = degree_lower_bound(graph);
            const i128 coarse_ticks = static_cast<i128>(63) * coarse_mass -
                                      static_cast<i128>(4) * graph.grid;
            std::cout << "PAIR " << q << ',' << r << " GRID " << decimal(graph.grid)
                      << " CELL_COUNTS " << graph.rank0_cells << ','
                      << graph.rank1_cells << ',' << graph.rank2_cells
                      << " RANK2_TOTAL " << decimal(graph.total) << " BODIES "
                      << result.bodies << " POSITIVE " << result.positive
                      << " ZERO " << result.zero << " MIN_TICKS "
                      << decimal(result.minimum) << " MIN_BODY "
                      << hex8(result.minimum_body) << " MAX_TICKS "
                      << decimal(result.maximum) << " MAX_BODY "
                      << hex8(result.maximum_body) << " DEGREE_BOUND_TICKS "
                      << decimal(coarse_ticks) << '\n';
        }
        std::cout << "GRAPH_IDENTITY L2(B)=W0+SUM_I_NOTIN_B_WI+SUM_I_J_NOTIN_B_WIJ\n"
                  << "SCOPE FINITE_EXACT_ALL_RANK9_BODIES_LISTED_PAIRS_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "RANK2_ALLBODY_AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}
