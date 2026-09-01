// Independent wide finite-exact audit of the rank-at-most-two wall skeleton.
//
// This implementation deliberately uses an event-state sweep.  It does not
// sample every speed at every open cell and does not import any theorem packet
// source.  Pool events are built once on the base grid and merged with two
// monotone pair-event streams after exact integral rescaling.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {
using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

constexpr i64 kBaseGrid = INT64_C(18241159416480);
constexpr unsigned kBodyRank = 9;
constexpr std::array<int, 30> kPool = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};

[[noreturn]] void die(const std::string& message) {
    throw std::runtime_error(message);
}
void demand(bool condition, const std::string& message) {
    if (!condition) die(message);
}

std::string dec(i128 value) {
    if (value == 0) return "0";
    const bool negative = value < 0;
    u128 magnitude = static_cast<u128>(value);
    if (negative) magnitude = u128{0} - magnitude;
    std::string text;
    while (magnitude != 0) {
        text.push_back(static_cast<char>('0' + magnitude % 10));
        magnitude /= 10;
    }
    if (negative) text.push_back('-');
    std::reverse(text.begin(), text.end());
    return text;
}

std::string hex8(u32 mask) {
    std::ostringstream output;
    output << std::hex << std::setw(8) << std::setfill('0') << mask;
    return output.str();
}

i128 gcd128(i128 a, i128 b) {
    while (b != 0) {
        const i128 remainder = a % b;
        a = b;
        b = remainder;
    }
    return a;
}

i128 exact_lcm(i128 a, i128 b) {
    demand(a > 0 && b > 0, "nonpositive LCM input");
    const i128 divisor = gcd128(a, b);
    const i128 candidate = static_cast<i128>(a / divisor) * b;
    demand(candidate > 0, "LCM overflow");
    return candidate;
}

struct Pair {
    int q = 0;
    int r = 0;
};

std::vector<Pair> read_pairs(const std::string& path, u64& raw_fnv) {
    std::ifstream input(path, std::ios::binary);
    demand(static_cast<bool>(input), "cannot open pair universe: " + path);
    std::vector<Pair> pairs;
    std::set<std::pair<int, int>> seen;
    raw_fnv = UINT64_C(14695981039346656037);
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        demand(!line.empty(), "blank universe row");
        for (unsigned char byte : line) {
            raw_fnv ^= byte;
            raw_fnv *= UINT64_C(1099511628211);
        }
        raw_fnv ^= static_cast<unsigned char>('\n');
        raw_fnv *= UINT64_C(1099511628211);
        const std::size_t comma = line.find(',');
        demand(comma != std::string::npos &&
                   line.find(',', comma + 1) == std::string::npos,
               "malformed universe row: " + line);
        const int q = std::stoi(line.substr(0, comma));
        const int r = std::stoi(line.substr(comma + 1));
        demand(q > 0 && q < r, "pair ordering/range failure: " + line);
        demand(seen.emplace(q, r).second, "duplicate universe row: " + line);
        pairs.push_back({q, r});
    }
    demand(input.eof(), "pair-universe read failure");
    demand(pairs.size() == 181194, "THM-4231 remainder cardinality changed");
    return pairs;
}

struct BaseEvent {
    i64 position = 0;
    u32 enter = 0;  // bits becoming safe immediately to the right
    u32 leave = 0;  // bits becoming unsafe immediately to the right
};

std::vector<BaseEvent> make_base_events() {
    struct Raw {
        i64 position;
        unsigned bit;
        bool entering;
    };
    std::vector<Raw> raw;
    for (unsigned bit = 0; bit < kPool.size(); ++bit) {
        const int speed = kPool[bit];
        const i64 divisor = 14LL * speed;
        demand(kBaseGrid % divisor == 0, "base grid no longer resolves pool");
        const i64 unit = kBaseGrid / divisor;
        for (int tooth = 0; tooth < speed; ++tooth) {
            raw.push_back({(14LL * tooth + 1) * unit, bit, true});
            raw.push_back({(14LL * tooth + 13) * unit, bit, false});
        }
    }
    std::sort(raw.begin(), raw.end(), [](const Raw& a, const Raw& b) {
        return std::tie(a.position, a.bit, a.entering) <
               std::tie(b.position, b.bit, b.entering);
    });
    std::vector<BaseEvent> events;
    for (std::size_t cursor = 0; cursor < raw.size();) {
        const i64 position = raw[cursor].position;
        BaseEvent event;
        event.position = position;
        while (cursor < raw.size() && raw[cursor].position == position) {
            const u32 bit = u32{1} << raw[cursor].bit;
            if (raw[cursor].entering)
                event.enter |= bit;
            else
                event.leave |= bit;
            ++cursor;
        }
        demand((event.enter & event.leave) == 0,
               "same pool bit enters and leaves at one wall");
        events.push_back(event);
    }
    demand(!events.empty() && events.front().position > 0 &&
               events.back().position < kBaseGrid,
           "bad base event range");
    return events;
}

class PairEvents {
  public:
    PairEvents(int speed, i128 grid) : speed_(speed), grid_(grid) {
        demand(grid % (14LL * speed) == 0, "pair wall not integral");
        unit_ = grid / (14LL * speed);
    }

    i128 next() const {
        if (tooth_ >= speed_) return grid_;
        const i128 numerator = 14LL * tooth_ + (entering_ ? 1 : 13);
        return numerator * unit_;
    }

    void cross(i128 position, bool& safe) {
        demand(next() == position && position < grid_,
               "pair-event merge desynchronized");
        safe = entering_;
        if (entering_) {
            entering_ = false;
        } else {
            entering_ = true;
            ++tooth_;
        }
    }

  private:
    int speed_ = 0;
    int tooth_ = 0;
    i128 grid_ = 0;
    i128 unit_ = 0;
    bool entering_ = true;
};

struct Graph {
    i128 grid = 0;
    i128 rank0 = 0;
    std::array<i128, 30> rank1{};
    std::array<std::array<i128, 30>, 30> rank2{};
    std::array<i128, 30> degree{};
    i128 total = 0;
    u64 cells0 = 0;
    u64 cells1 = 0;
    u64 cells2 = 0;
};

Graph event_sweep_graph(int q, int r,
                        const std::vector<BaseEvent>& base_events) {
    Graph graph;
    graph.grid = exact_lcm(kBaseGrid, 14LL * q);
    graph.grid = exact_lcm(graph.grid, 14LL * r);
    demand(graph.grid % kBaseGrid == 0, "grid lost base divisibility");
    const i128 scale = graph.grid / kBaseGrid;

    std::size_t pool_cursor = 0;
    PairEvents q_events(q, graph.grid), r_events(r, graph.grid);
    bool q_safe = false, r_safe = false;
    u32 failures = (u32{1} << 30) - 1;
    i128 left = 0;

    while (left < graph.grid) {
        const i128 pool_next = pool_cursor < base_events.size()
                                  ? static_cast<i128>(base_events[pool_cursor].position) * scale
                                  : graph.grid;
        const i128 right = std::min({graph.grid, pool_next,
                                    q_events.next(), r_events.next()});
        demand(right > left, "non-increasing merged wall sequence");

        if (q_safe && r_safe) {
            const unsigned rank = std::popcount(failures);
            const i128 width = right - left;
            if (rank == 0) {
                graph.rank0 += width;
                ++graph.cells0;
            } else if (rank == 1) {
                const unsigned vertex = std::countr_zero(failures);
                graph.rank1[vertex] += width;
                ++graph.cells1;
            } else if (rank == 2) {
                const unsigned first = std::countr_zero(failures);
                const unsigned second =
                    std::countr_zero(failures & (failures - 1));
                graph.rank2[first][second] += width;
                graph.rank2[second][first] += width;
                ++graph.cells2;
            }
        }
        left = right;
        if (left == graph.grid) break;

        if (pool_next == left) {
            const BaseEvent& event = base_events[pool_cursor++];
            failures &= ~event.enter;
            failures |= event.leave;
        }
        if (q_events.next() == left) q_events.cross(left, q_safe);
        if (r_events.next() == left) r_events.cross(left, r_safe);
    }

    demand(pool_cursor == base_events.size(), "not all pool events consumed");
    demand(!q_safe && !r_safe, "pair state did not close unsafe");
    demand(failures == (u32{1} << 30) - 1,
           "pool state did not close unsafe");

    graph.total = graph.rank0;
    for (unsigned v = 0; v < 30; ++v) {
        graph.total += graph.rank1[v];
        graph.degree[v] = graph.rank1[v];
    }
    for (unsigned u = 0; u < 30; ++u) {
        for (unsigned v = u + 1; v < 30; ++v) {
            const i128 weight = graph.rank2[u][v];
            graph.total += weight;
            graph.degree[u] += weight;
            graph.degree[v] += weight;
        }
    }
    demand(graph.total >= 0 && graph.total <= graph.grid,
           "rank-two mass outside grid");
    return graph;
}

struct Coarse {
    i128 degree_sum = 0;
    i128 mass = 0;
    i128 ticks = 0;
    u32 top9 = 0;
};

Coarse coarse_bound(const Graph& graph) {
    std::array<unsigned, 30> order{};
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](unsigned a, unsigned b) {
        if (graph.degree[a] != graph.degree[b])
            return graph.degree[a] > graph.degree[b];
        return a < b;
    });
    Coarse coarse;
    for (unsigned i = 0; i < kBodyRank; ++i) {
        coarse.degree_sum += graph.degree[order[i]];
        coarse.top9 |= u32{1} << order[i];
    }
    coarse.mass = static_cast<i128>(graph.total) - coarse.degree_sum;
    coarse.ticks = 63 * coarse.mass - 4 * graph.grid;
    return coarse;
}

i128 reward_of(const Graph& graph, u32 body) {
    demand(std::popcount(body) == kBodyRank, "reward body rank changed");
    i128 reward = 0;
    for (unsigned v = 0; v < 30; ++v)
        if ((body >> v) & 1U) reward += graph.degree[v];
    for (unsigned u = 0; u < 30; ++u)
        if ((body >> u) & 1U)
            for (unsigned v = u + 1; v < 30; ++v)
                if ((body >> v) & 1U) reward -= graph.rank2[u][v];
    return reward;
}

i128 survivor_direct(const Graph& graph, u32 body) {
    demand(std::popcount(body) == kBodyRank, "direct body rank changed");
    i128 mass = graph.rank0;
    for (unsigned v = 0; v < 30; ++v)
        if (((body >> v) & 1U) == 0) mass += graph.rank1[v];
    for (unsigned u = 0; u < 30; ++u)
        if (((body >> u) & 1U) == 0)
            for (unsigned v = u + 1; v < 30; ++v)
                if (((body >> v) & 1U) == 0) mass += graph.rank2[u][v];
    return mass;
}

// Exact cardinality-constrained quadratic optimizer.
//
// For a partial selected set S, the marginal reward of an unprocessed vertex v
// is d(v)-sum_{u in S}w(u,v).  If T is any completion, its true added reward is
// the sum of these current marginals minus the nonnegative edge weight internal
// to T.  Therefore the sum of the largest |T| current marginals is an
// admissible upper bound.  Pruning below the incumbent cannot discard an
// optimum.  Equality is retained so the least-mask tie rule is exact too.
class ExactOptimizer {
  public:
    explicit ExactOptimizer(const Graph& graph) : graph_(graph) {
        std::iota(order_.begin(), order_.end(), 0);
        std::sort(order_.begin(), order_.end(), [&](unsigned a, unsigned b) {
            if (graph_.degree[a] != graph_.degree[b])
                return graph_.degree[a] > graph_.degree[b];
            return a < b;
        });
        seed_greedy();
    }

    void solve() { descend(0, kBodyRank, 0, 0); }
    i128 best_reward() const { return best_reward_; }
    u32 best_body() const { return best_body_; }
    u64 nodes() const { return nodes_; }
    u64 prunes() const { return prunes_; }

  private:
    const Graph& graph_;
    std::array<unsigned, 30> order_{};
    i128 best_reward_ = -((static_cast<i128>(1)) << 120);
    u32 best_body_ = 0;
    u64 nodes_ = 0;
    u64 prunes_ = 0;

    void consider(u32 body, i128 reward) {
        demand(std::popcount(body) == kBodyRank, "optimizer leaf rank changed");
        demand(reward == reward_of(graph_, body), "incremental reward drift");
        if (reward > best_reward_ ||
            (reward == best_reward_ && body < best_body_)) {
            best_reward_ = reward;
            best_body_ = body;
        }
    }

    i128 marginal(unsigned vertex, u32 selected) const {
        i128 gain = graph_.degree[vertex];
        u32 copy = selected;
        while (copy != 0) {
            const unsigned prior = std::countr_zero(copy);
            gain -= graph_.rank2[vertex][prior];
            copy &= copy - 1;
        }
        return gain;
    }

    void seed_greedy() {
        u32 selected = 0;
        i128 reward = 0;
        for (unsigned step = 0; step < kBodyRank; ++step) {
            unsigned chosen = 30;
            i128 chosen_gain = -((static_cast<i128>(1)) << 120);
            for (unsigned vertex = 0; vertex < 30; ++vertex) {
                if ((selected >> vertex) & 1U) continue;
                const i128 gain = marginal(vertex, selected);
                if (chosen == 30 || gain > chosen_gain ||
                    (gain == chosen_gain && vertex < chosen)) {
                    chosen = vertex;
                    chosen_gain = gain;
                }
            }
            demand(chosen < 30, "greedy seed exhausted vertices");
            selected |= u32{1} << chosen;
            reward += chosen_gain;
        }
        consider(selected, reward);
    }

    void take_all(unsigned position, unsigned need, u32 selected,
                  i128 reward) {
        for (unsigned cursor = position; cursor < 30; ++cursor) {
            const unsigned vertex = order_[cursor];
            reward += marginal(vertex, selected);
            selected |= u32{1} << vertex;
        }
        demand(std::popcount(selected) == kBodyRank && need == 30 - position,
               "forced completion rank failure");
        consider(selected, reward);
    }

    void descend(unsigned position, unsigned need, u32 selected, i128 reward) {
        ++nodes_;
        if (need == 0) {
            consider(selected, reward);
            return;
        }
        const unsigned available = 30 - position;
        if (available < need) return;
        if (available == need) {
            take_all(position, need, selected, reward);
            return;
        }

        std::vector<i128> gains;
        gains.reserve(available);
        for (unsigned cursor = position; cursor < 30; ++cursor)
            gains.push_back(marginal(order_[cursor], selected));
        std::nth_element(gains.begin(), gains.begin() + need, gains.end(),
                         std::greater<i128>());
        i128 upper = reward;
        for (unsigned i = 0; i < need; ++i) upper += gains[i];
        if (upper < best_reward_) {
            ++prunes_;
            return;
        }

        const unsigned vertex = order_[position];
        const i128 gain = marginal(vertex, selected);
        descend(position + 1, need - 1, selected | (u32{1} << vertex),
                reward + gain);
        descend(position + 1, need, selected, reward);
    }
};

void fnv_add(u64& hash, const std::string& line) {
    for (unsigned char byte : line) {
        hash ^= byte;
        hash *= UINT64_C(1099511628211);
    }
    hash ^= static_cast<unsigned char>('\n');
    hash *= UINT64_C(1099511628211);
}

std::string fnv_hex(u64 hash) {
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << hash;
    return out.str();
}
}  // namespace

int main(int argc, char** argv) {
    try {
        demand(argc == 4,
               "usage: rank2_event_sweep_cleanroom UNIVERSE SCREEN_CSV EXACT_CSV");
        u64 input_fnv = 0;
        const std::vector<Pair> pairs = read_pairs(argv[1], input_fnv);
        const std::vector<BaseEvent> base_events = make_base_events();

        std::ofstream screen(argv[2]);
        demand(static_cast<bool>(screen), "cannot create screen CSV");
        screen << "q,r,grid,cells0,cells1,cells2,mass0,mass1,mass2,total_mass,top9_degree_sum,coarse_mass,coarse_ticks,coarse_positive,top9_mask_hex\n";
        u64 screen_fnv = UINT64_C(14695981039346656037);
        std::vector<Pair> hostile;
        unsigned coarse_positive = 0;
        int highest_hostile = -1;
        bool have_smallest_positive_coarse = false;
        i128 smallest_positive_coarse_ticks = 0;
        Pair smallest_positive_coarse_pair{};
        u32 smallest_positive_coarse_body = 0;
        unsigned wide_grid_count = 0;
        Pair wide_grid_pair{};
        i128 wide_grid_value = 0;

        for (const Pair pair : pairs) {
            Graph graph = event_sweep_graph(pair.q, pair.r, base_events);
            const Coarse coarse = coarse_bound(graph);
            i128 mass1 = 0, mass2 = 0;
            for (i128 weight : graph.rank1) mass1 += weight;
            for (unsigned u = 0; u < 30; ++u)
                for (unsigned v = u + 1; v < 30; ++v)
                    mass2 += graph.rank2[u][v];
            demand(static_cast<i128>(graph.rank0) + mass1 + mass2 == graph.total,
                   "rank mass partition failure");

            const bool positive = coarse.ticks > 0;
            coarse_positive += positive;
            if (positive && (!have_smallest_positive_coarse ||
                             coarse.ticks < smallest_positive_coarse_ticks)) {
                have_smallest_positive_coarse = true;
                smallest_positive_coarse_ticks = coarse.ticks;
                smallest_positive_coarse_pair = pair;
                smallest_positive_coarse_body = coarse.top9;
            }
            if (!positive) {
                highest_hostile = std::max(highest_hostile, pair.r);
                hostile.push_back(pair);
            }
            if (graph.grid > std::numeric_limits<i64>::max()) {
                ++wide_grid_count;
                wide_grid_pair = pair;
                wide_grid_value = graph.grid;
            }
            std::ostringstream row;
            row << pair.q << ',' << pair.r << ',' << dec(graph.grid) << ','
                << graph.cells0 << ',' << graph.cells1 << ',' << graph.cells2
                << ',' << dec(graph.rank0) << ',' << dec(mass1) << ',' << dec(mass2)
                << ',' << dec(graph.total) << ',' << dec(coarse.degree_sum) << ','
                << dec(coarse.mass) << ',' << dec(coarse.ticks) << ','
                << (positive ? 1 : 0) << ',' << hex8(coarse.top9);
            screen << row.str() << '\n';
            fnv_add(screen_fnv, row.str());
        }
        screen.close();
        demand(static_cast<bool>(screen), "screen CSV write failure");
        demand(wide_grid_count == 1 && wide_grid_pair.q == 713 &&
                   wide_grid_pair.r == 719 &&
                   wide_grid_value == static_cast<i128>(INT64_C(935127565138022256)) * 10,
               "signed-int64 overflow control changed");

        std::ofstream exact(argv[3]);
        demand(static_cast<bool>(exact), "cannot create exact CSV");
        exact << "q,r,grid,coarse_ticks,minimum_mass,minimum_ticks,minimum_body_hex,direct_replay_mass,identity_ok,nodes,prunes\n";
        u64 exact_fnv = UINT64_C(14695981039346656037);
        unsigned exact_positive = 0;
        i128 hostile_exact_minimum = 0;
        Pair hostile_exact_pair{};
        u32 hostile_exact_body = 0;
        u64 total_nodes = 0, total_prunes = 0;

        for (const Pair pair : hostile) {
            const Graph graph = event_sweep_graph(pair.q, pair.r, base_events);
            const Coarse coarse = coarse_bound(graph);
            ExactOptimizer optimizer(graph);
            optimizer.solve();
            const i128 minimum_mass =
                static_cast<i128>(graph.total) - optimizer.best_reward();
            const i128 minimum_ticks =
                63 * minimum_mass - 4 * graph.grid;
            const i128 direct = survivor_direct(graph, optimizer.best_body());
            const bool identity = direct == minimum_mass;
            demand(identity, "direct/degree identity failed");
            demand(minimum_mass >= 0 && minimum_mass <= graph.grid,
                   "minimum mass outside grid");
            exact_positive += minimum_ticks > 0;
            if (total_nodes == 0 || minimum_ticks < hostile_exact_minimum) {
                hostile_exact_minimum = minimum_ticks;
                hostile_exact_pair = pair;
                hostile_exact_body = optimizer.best_body();
            }
            total_nodes += optimizer.nodes();
            total_prunes += optimizer.prunes();

            std::ostringstream row;
            row << pair.q << ',' << pair.r << ',' << dec(graph.grid) << ','
                << dec(coarse.ticks) << ',' << dec(minimum_mass) << ','
                << dec(minimum_ticks) << ',' << hex8(optimizer.best_body())
                << ',' << dec(direct) << ",1," << optimizer.nodes() << ','
                << optimizer.prunes();
            exact << row.str() << '\n';
            fnv_add(exact_fnv, row.str());
        }
        exact.close();
        demand(static_cast<bool>(exact), "exact CSV write failure");

        std::cout << "LRC14_RANK2_WIDE_CLEANROOM_EVENT_SWEEP_AUDIT_V1\n";
        std::cout << "INPUT_PAIRS " << pairs.size()
                  << " INPUT_NORMALIZED_FNV " << fnv_hex(input_fnv)
                  << " BASE_EVENT_COORDINATES " << base_events.size() << "\n";
        std::cout << "COARSE_POSITIVE " << coarse_positive
                  << " COARSE_NONPOSITIVE " << hostile.size()
                  << " HIGHEST_NONPOSITIVE_ENDPOINT " << highest_hostile
                  << " SCREEN_DATA_FNV " << fnv_hex(screen_fnv) << "\n";
        demand(have_smallest_positive_coarse,
               "no strictly positive coarse-bound row");
        std::cout << "SMALLEST_POSITIVE_COARSE_BOUND_TICKS "
                  << dec(smallest_positive_coarse_ticks) << " PAIR "
                  << smallest_positive_coarse_pair.q << ','
                  << smallest_positive_coarse_pair.r << " TOP9_BODY "
                  << hex8(smallest_positive_coarse_body) << "\n";
        std::cout << "SIGNED_INT64_OVERFLOW_GRIDS " << wide_grid_count
                  << " UNIQUE_PAIR " << wide_grid_pair.q << ',' << wide_grid_pair.r
                  << " GRID " << dec(wide_grid_value) << "\n";
        std::cout << "EXACT_PAIRS " << hostile.size()
                  << " EXACT_POSITIVE " << exact_positive
                  << " EXACT_NONPOSITIVE " << hostile.size() - exact_positive
                  << " EXACT_DATA_FNV " << fnv_hex(exact_fnv) << "\n";
        std::cout << "HOSTILE107_EXACT_MINIMUM_TICKS "
                  << dec(hostile_exact_minimum) << " PAIR "
                  << hostile_exact_pair.q << ',' << hostile_exact_pair.r
                  << " BODY " << hex8(hostile_exact_body) << "\n";
        std::cout << "BRANCH_NODES " << total_nodes
                  << " BRANCH_PRUNES " << total_prunes << "\n";
        std::cout << "BOUND L2(B)>=W-SUM_NINE_LARGEST_WEIGHTED_DEGREES\n";
        std::cout << "EXACTNESS CURRENT_MARGINAL_TOP_K_UPPER_BOUND_WITH_NONNEGATIVE_FUTURE_EDGE_PENALTY\n";
        std::cout << "SCOPE FINITE_EXACT_THM4231_REMAINDER181194_RANK2_WALL_AUDIT_NO_PHYSICAL_ENTRY_NO_LRC14\n";
        std::cout << "VERDICT "
                  << ((pairs.size() == 181194 &&
                       exact_positive == hostile.size())
                          ? "PASS"
                          : "FAIL")
                  << "\n";
        return (pairs.size() == 181194 && exact_positive == hostile.size())
                   ? 0
                   : 2;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
