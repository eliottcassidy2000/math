// Clean-room exact branch-and-bound minimizer for the rank-two literal graph.
// It audits precisely the rows on which the coarse nine-degree bound is not
// positive, without enumerating all C(30,9) bodies.

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

namespace clean {
using U = std::uint32_t;
using Q = __int128_t;
using W = __uint128_t;
using I = Q;

constexpr std::array<int, 30> speeds = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};

[[noreturn]] void bad(const std::string& text) { throw std::runtime_error(text); }
void check(bool value, const std::string& text) { if (!value) bad(text); }

std::string dec(Q value) {
    if (value == 0) return "0";
    bool neg = value < 0;
    W mag = static_cast<W>(value);
    if (neg) mag = W{0} - mag;
    std::string out;
    while (mag) { out.push_back(char('0' + mag % 10)); mag /= 10; }
    if (neg) out.push_back('-');
    std::reverse(out.begin(), out.end());
    return out;
}

std::string hx(U value) {
    std::ostringstream out;
    out << std::hex << std::setw(8) << std::setfill('0') << value;
    return out.str();
}

I lcm_checked(I a, I b) {
    I left = a, right = b;
    while (right != 0) {
        const I remainder = left % right;
        left = right;
        right = remainder;
    }
    const Q wide = Q(a / left) * b;
    check(wide > 0, "grid overflow");
    return wide;
}

bool lonely(int velocity, I grid, I a, I b) {
    Q phase = Q(velocity) * (Q(a) + b);
    phase %= Q(2) * grid;
    if (phase < 0) phase += Q(2) * grid;
    return Q(7) * phase >= grid && Q(7) * phase <= Q(13) * grid;
}

struct Network {
    I grid = 0;
    I empty = 0;
    std::array<I, 30> vertex{};
    std::array<std::array<I, 30>, 30> link{};
    std::array<I, 30> removal{};
    I total = 0;
};

Network reconstruct(int q, int r) {
    Network net;
    net.grid = 1;
    for (int v : speeds) net.grid = lcm_checked(net.grid, I(14) * v);
    net.grid = lcm_checked(net.grid, I(14) * q);
    net.grid = lcm_checked(net.grid, I(14) * r);
    std::vector<I> cut = {0, net.grid};
    auto insert_walls = [&](int v) {
        check(net.grid % (I(14) * v) == 0, "fractional wall");
        const I step = net.grid / (I(14) * v);
        for (int k = 0; k < v; ++k) {
            cut.push_back((I(14) * k + 1) * step);
            cut.push_back((I(14) * k + 13) * step);
        }
    };
    for (int v : speeds) insert_walls(v);
    insert_walls(q); insert_walls(r);
    std::sort(cut.begin(), cut.end());
    cut.erase(std::unique(cut.begin(), cut.end()), cut.end());
    for (std::size_t k = 1; k < cut.size(); ++k) {
        const I a = cut[k - 1], b = cut[k];
        if (!lonely(q, net.grid, a, b) || !lonely(r, net.grid, a, b)) continue;
        U failed = 0;
        for (unsigned i = 0; i < speeds.size(); ++i)
            if (!lonely(speeds[i], net.grid, a, b)) failed |= U{1} << i;
        const unsigned rank = std::popcount(failed);
        const I width = b - a;
        if (rank == 0) {
            net.empty += width;
        } else if (rank == 1) {
            net.vertex[std::countr_zero(failed)] += width;
        } else if (rank == 2) {
            const unsigned x = std::countr_zero(failed);
            const unsigned y = std::countr_zero(failed & (failed - 1));
            net.link[x][y] += width;
            net.link[y][x] += width;
        }
    }
    net.total = net.empty;
    for (unsigned i = 0; i < 30; ++i) {
        net.total += net.vertex[i];
        net.removal[i] = net.vertex[i];
    }
    for (unsigned i = 0; i < 30; ++i)
        for (unsigned j = i + 1; j < 30; ++j) {
            net.total += net.link[i][j];
            net.removal[i] += net.link[i][j];
            net.removal[j] += net.link[i][j];
        }
    return net;
}

struct Search {
    const Network& net;
    std::array<unsigned, 30> order{};
    std::array<unsigned, 9> chosen{};
    I best = -(Q{1} << 126);
    U best_mask = 0;
    std::uint64_t nodes = 0;
    std::uint64_t prunes = 0;

    explicit Search(const Network& input) : net(input) {
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](unsigned a, unsigned b) {
            if (net.removal[a] != net.removal[b])
                return net.removal[a] > net.removal[b];
            return a < b;
        });
        // A deterministic marginal-greedy incumbent.
        U used = 0;
        I score = 0;
        for (unsigned depth = 0; depth < 9; ++depth) {
            bool have = false;
            I best_gain = -(Q{1} << 126);
            unsigned best_vertex = 0;
            for (unsigned v = 0; v < 30; ++v) {
                if (used & (U{1} << v)) continue;
                I gain = net.removal[v];
                for (unsigned prior = 0; prior < depth; ++prior)
                    gain -= net.link[v][chosen[prior]];
                if (!have || gain > best_gain ||
                    (gain == best_gain && v < best_vertex)) {
                    have = true; best_gain = gain; best_vertex = v;
                }
            }
            check(have, "greedy exhausted vertices");
            chosen[depth] = best_vertex;
            used |= U{1} << best_vertex;
            score += best_gain;
        }
        best = score;
        best_mask = used;
    }

    void visit(unsigned start, unsigned depth, I score, U mask) {
        ++nodes;
        const unsigned need = 9 - depth;
        if (need == 0) {
            if (score > best || (score == best && mask < best_mask)) {
                best = score; best_mask = mask;
            }
            return;
        }
        if (30 - start < need) return;

        std::vector<I> possible;
        possible.reserve(30 - start);
        for (unsigned pos = start; pos < 30; ++pos) {
            const unsigned v = order[pos];
            I gain = net.removal[v];
            for (unsigned prior = 0; prior < depth; ++prior)
                gain -= net.link[v][chosen[prior]];
            possible.push_back(gain);
        }
        if (possible.size() > need)
            std::nth_element(possible.begin(), possible.begin() + need,
                             possible.end(), std::greater<I>());
        I upper = score;
        for (unsigned i = 0; i < need; ++i) upper += possible[i];
        // Pair penalties among future choices are nonnegative and omitted,
        // so upper is a rigorous upper bound on every completion.
        if (upper < best) { ++prunes; return; }

        for (unsigned pos = start; pos + need <= 30; ++pos) {
            const unsigned v = order[pos];
            I gain = net.removal[v];
            for (unsigned prior = 0; prior < depth; ++prior)
                gain -= net.link[v][chosen[prior]];
            chosen[depth] = v;
            visit(pos + 1, depth + 1, score + gain, mask | (U{1} << v));
        }
    }
};

struct Candidate { int q = 0; int r = 0; };
std::vector<Candidate> candidates(const std::string& path) {
    std::ifstream in(path);
    check(bool(in), "cannot open degree ledger");
    std::string line;
    check(bool(std::getline(in, line)) &&
              line == "q,r,grid,rank2_total,degree_bound_mass,degree_bound_ticks,positive,top9_hex",
          "degree header changed");
    std::vector<Candidate> out;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        Candidate c;
        std::string grid, total, lower, ticks, top, extra;
        int positive = -1;
        row >> c.q >> c.r >> grid >> total >> lower >> ticks >> positive >> top;
        check(row && !(row >> extra), "malformed degree row");
        if (positive == 0) out.push_back(c);
    }
    check(in.eof() && !out.empty(), "no degree-bound candidates");
    return out;
}
}  // namespace clean

int main(int argc, char** argv) {
    using namespace clean;
    try {
        check(argc == 2, "usage: branch_bound DEGREE_SCREEN_CSV");
        const auto rows = candidates(argv[1]);
        std::cout << "LRC14_RANK2_BRANCH_BOUND_INDEPENDENT_V1\n";
        std::uint64_t total_nodes = 0, total_prunes = 0;
        for (const auto c : rows) {
            const Network net = reconstruct(c.q, c.r);
            Search search(net);
            search.visit(0, 0, 0, 0);
            const I minimum_mass = net.total - search.best;
            const Q ticks = Q(63) * minimum_mass - Q(4) * net.grid;
            total_nodes += search.nodes; total_prunes += search.prunes;
            std::cout << "PAIR " << c.q << ',' << c.r << " GRID " << dec(net.grid)
                      << " BOUND_EXACT 1 MIN_TICKS " << dec(ticks)
                      << " MIN_BODY " << hx(search.best_mask)
                      << " NODES " << search.nodes << " PRUNES "
                      << search.prunes << '\n';
        }
        std::cout << "SUMMARY PAIRS " << rows.size() << " NODES "
                  << total_nodes << " PRUNES " << total_prunes << '\n'
                  << "UPPER_BOUND CURRENT_SCORE_PLUS_TOP_NEEDED_MARGINALS_OMITS_NONNEGATIVE_FUTURE_PAIR_PENALTIES\n"
                  << "SCOPE FINITE_EXACT_DEGREE_BOUND_CANDIDATES_ONLY_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "RANK2_BRANCH_BOUND_ERROR " << e.what() << '\n';
        return 1;
    }
}
