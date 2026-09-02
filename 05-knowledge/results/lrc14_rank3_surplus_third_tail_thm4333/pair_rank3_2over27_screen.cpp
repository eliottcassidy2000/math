// Canonical rank-three strengthening census over THM-4231's pair remainder.
// Reuses the independently audited wide event machinery, but builds a cubic
// failure hypergraph and targets 2/27 rather than THM-4326's 4/63.

#define main rank2_wide_hidden_main
#include "04-computation/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/independent_cleanroom/rank2_event_sweep_wide_cleanroom.cpp"
#undef main

#include <chrono>

namespace {
constexpr unsigned kTripleCount = 4060;

struct Indices {
    std::array<std::array<int, 30>, 30> pair{};
    std::array<std::array<std::array<int, 30>, 30>, 30> triple{};
    std::array<std::pair<unsigned, unsigned>, 435> pair_vertices{};
    std::array<std::array<unsigned, 3>, kTripleCount> triple_vertices{};

    Indices() {
        for (auto& row : pair) row.fill(-1);
        for (auto& slab : triple)
            for (auto& row : slab) row.fill(-1);
        unsigned pi = 0;
        for (unsigned a = 0; a < 30; ++a)
            for (unsigned b = a + 1; b < 30; ++b) {
                pair[a][b] = pair[b][a] = static_cast<int>(pi);
                pair_vertices[pi++] = {a, b};
            }
        demand(pi == pair_vertices.size(), "pair index cardinality");
        unsigned ti = 0;
        for (unsigned a = 0; a < 30; ++a)
            for (unsigned b = a + 1; b < 30; ++b)
                for (unsigned c = b + 1; c < 30; ++c) {
                    triple[a][b][c] = triple[a][c][b] =
                        triple[b][a][c] = triple[b][c][a] =
                        triple[c][a][b] = triple[c][b][a] =
                            static_cast<int>(ti);
                    triple_vertices[ti++] = {a, b, c};
                }
        demand(ti == kTripleCount, "triple index cardinality");
    }
};

const Indices kIndices;

struct Graph3 {
    i128 grid = 0;
    i128 rank0 = 0;
    std::array<i128, 30> rank1{};
    std::array<i128, 435> rank2{};
    std::array<i128, kTripleCount> rank3{};
    std::array<i128, 30> degree{};
    std::array<std::array<i128, 30>, 30> codegree{};
    i128 total = 0;
    u64 cells0 = 0, cells1 = 0, cells2 = 0, cells3 = 0;
};

Graph3 event_sweep_graph3(int q, int r,
                          const std::vector<BaseEvent>& base_events) {
    Graph3 graph;
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
                                   ? i128(base_events[pool_cursor].position) * scale
                                   : graph.grid;
        const i128 right = std::min(
            {graph.grid, pool_next, q_events.next(), r_events.next()});
        demand(right > left, "non-increasing merged wall sequence");

        if (q_safe && r_safe) {
            const unsigned rank = std::popcount(failures);
            const i128 width = right - left;
            if (rank == 0) {
                graph.rank0 += width;
                ++graph.cells0;
            } else if (rank == 1) {
                graph.rank1[std::countr_zero(failures)] += width;
                ++graph.cells1;
            } else if (rank == 2) {
                const unsigned a = std::countr_zero(failures);
                const unsigned b = std::countr_zero(failures & (failures - 1));
                graph.rank2[kIndices.pair[a][b]] += width;
                ++graph.cells2;
            } else if (rank == 3) {
                u32 copy = failures;
                const unsigned a = std::countr_zero(copy);
                copy &= copy - 1;
                const unsigned b = std::countr_zero(copy);
                copy &= copy - 1;
                const unsigned c = std::countr_zero(copy);
                graph.rank3[kIndices.triple[a][b][c]] += width;
                ++graph.cells3;
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
        graph.degree[v] += graph.rank1[v];
    }
    for (unsigned index = 0; index < graph.rank2.size(); ++index) {
        const i128 weight = graph.rank2[index];
        const auto [a, b] = kIndices.pair_vertices[index];
        graph.total += weight;
        graph.degree[a] += weight;
        graph.degree[b] += weight;
        graph.codegree[a][b] += weight;
        graph.codegree[b][a] += weight;
    }
    for (unsigned index = 0; index < graph.rank3.size(); ++index) {
        const i128 weight = graph.rank3[index];
        const auto [a, b, c] = kIndices.triple_vertices[index];
        graph.total += weight;
        graph.degree[a] += weight;
        graph.degree[b] += weight;
        graph.degree[c] += weight;
        graph.codegree[a][b] += weight;
        graph.codegree[b][a] += weight;
        graph.codegree[a][c] += weight;
        graph.codegree[c][a] += weight;
        graph.codegree[b][c] += weight;
        graph.codegree[c][b] += weight;
    }
    demand(graph.total >= 0 && graph.total <= graph.grid,
           "rank-three mass outside grid");
    return graph;
}
i128 triple_weight(const Graph3& graph, unsigned a, unsigned b, unsigned c) {
    return graph.rank3[kIndices.triple[a][b][c]];
}

i128 reward3(const Graph3& graph, u32 body) {
    demand(std::popcount(body) == kBodyRank, "reward body rank changed");
    std::array<unsigned, kBodyRank> vertices{};
    unsigned count = 0;
    u32 copy = body;
    i128 reward = 0;
    while (copy != 0) {
        const unsigned v = std::countr_zero(copy);
        vertices[count++] = v;
        reward += graph.degree[v];
        copy &= copy - 1;
    }
    for (unsigned i = 0; i < count; ++i)
        for (unsigned j = i + 1; j < count; ++j)
            reward -= graph.codegree[vertices[i]][vertices[j]];
    for (unsigned i = 0; i < count; ++i)
        for (unsigned j = i + 1; j < count; ++j)
            for (unsigned k = j + 1; k < count; ++k)
                reward += triple_weight(graph, vertices[i], vertices[j],
                                         vertices[k]);
    return reward;
}

struct Exact3 {
    explicit Exact3(const Graph3& graph) : graph_(graph) {
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
    const Graph3& graph_;
    std::array<unsigned, 30> order_{};
    i128 best_reward_ = -(i128(1) << 120);
    u32 best_body_ = 0;
    u64 nodes_ = 0, prunes_ = 0;

    i128 marginal(unsigned vertex, u32 selected) const {
        i128 gain = graph_.degree[vertex];
        std::array<unsigned, kBodyRank> prior{};
        unsigned count = 0;
        u32 copy = selected;
        while (copy != 0) {
            const unsigned u = std::countr_zero(copy);
            prior[count++] = u;
            gain -= graph_.codegree[vertex][u];
            copy &= copy - 1;
        }
        for (unsigned i = 0; i < count; ++i)
            for (unsigned j = i + 1; j < count; ++j)
                gain += triple_weight(graph_, vertex, prior[i], prior[j]);
        return gain;
    }

    void consider(u32 body, i128 reward) {
        demand(std::popcount(body) == kBodyRank, "optimizer leaf rank");
        demand(reward == reward3(graph_, body), "incremental reward drift");
        if (reward > best_reward_ ||
            (reward == best_reward_ && body < best_body_)) {
            best_reward_ = reward;
            best_body_ = body;
        }
    }

    void seed_greedy() {
        u32 selected = 0;
        i128 reward = 0;
        for (unsigned step = 0; step < kBodyRank; ++step) {
            unsigned best = 30;
            i128 best_gain = -(i128(1) << 120);
            for (unsigned v = 0; v < 30; ++v) {
                if ((selected >> v) & 1U) continue;
                const i128 gain = marginal(v, selected);
                if (best == 30 || gain > best_gain ||
                    (gain == best_gain && v < best)) {
                    best = v;
                    best_gain = gain;
                }
            }
            demand(best < 30, "greedy exhausted vertices");
            selected |= u32{1} << best;
            reward += best_gain;
        }
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
            for (unsigned cursor = position; cursor < 30; ++cursor) {
                const unsigned v = order_[cursor];
                reward += marginal(v, selected);
                selected |= u32{1} << v;
            }
            consider(selected, reward);
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
        // Coverage is monotone submodular: future true marginals can only
        // decrease.  Hence this top-current-marginal sum is admissible.
        if (upper < best_reward_) {
            ++prunes_;
            return;
        }
        const unsigned v = order_[position];
        const i128 gain = marginal(v, selected);
        descend(position + 1, need - 1, selected | (u32{1} << v),
                reward + gain);
        descend(position + 1, need, selected, reward);
    }
};

struct Coarse3 {
    i128 mass = 0, ticks = 0;
    u32 top9 = 0;
};

Coarse3 coarse3(const Graph3& graph) {
    std::array<unsigned, 30> order{};
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](unsigned a, unsigned b) {
        if (graph.degree[a] != graph.degree[b])
            return graph.degree[a] > graph.degree[b];
        return a < b;
    });
    Coarse3 result;
    result.mass = graph.total;
    for (unsigned i = 0; i < kBodyRank; ++i) {
        result.mass -= graph.degree[order[i]];
        result.top9 |= u32{1} << order[i];
    }
    result.ticks = 27 * result.mass - 2 * graph.grid;
    return result;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        demand(argc == 4,
               "usage: pair_rank3_2over27_screen PAIRS BAD_CSV EXACT_CSV");
        u64 input_fnv = 0;
        const auto pairs = read_pairs(argv[1], input_fnv);
        const auto base_events = make_base_events();
        std::ofstream bad_out(argv[2]);
        std::ofstream exact_out(argv[3]);
        demand(bool(bad_out) && bool(exact_out), "cannot create outputs");
        bad_out << "q,r,grid,rank3_total,coarse_mass,coarse_ticks,positive,top9\n";
        exact_out << "q,r,grid,rank3_total,min_mass,min_ticks,least_body,nodes,prunes\n";

        std::vector<Pair> fallback;
        u64 coarse_positive = 0;
        const auto begin = std::chrono::steady_clock::now();
        for (std::size_t index = 0; index < pairs.size(); ++index) {
            const Pair pair = pairs[index];
            const Graph3 graph = event_sweep_graph3(pair.q, pair.r, base_events);
            const Coarse3 coarse = coarse3(graph);
            const bool positive = coarse.ticks > 0;
            if (positive) {
                ++coarse_positive;
            } else {
                fallback.push_back(pair);
            }
            bad_out << pair.q << ',' << pair.r << ',' << dec(graph.grid)
                    << ',' << dec(graph.total) << ',' << dec(coarse.mass)
                    << ',' << dec(coarse.ticks) << ',' << positive << ','
                    << hex8(coarse.top9) << '\n';
            if ((index + 1) % 20000 == 0)
                std::cerr << "SCREEN_PROGRESS " << (index + 1) << '/'
                          << pairs.size() << " FALLBACK " << fallback.size()
                          << '\n';
        }
        demand(bad_out.good(), "bad CSV write failure");

        u64 exact_positive = 0, total_nodes = 0, total_prunes = 0;
        bool have_min = false;
        i128 min_ticks = 0, min_grid = 1;
        Pair min_pair{};
        u32 min_body = 0;
        for (std::size_t index = 0; index < fallback.size(); ++index) {
            const Pair pair = fallback[index];
            const Graph3 graph = event_sweep_graph3(pair.q, pair.r, base_events);
            Exact3 optimizer(graph);
            optimizer.solve();
            const i128 min_mass = graph.total - optimizer.best_reward();
            const i128 ticks = 27 * min_mass - 2 * graph.grid;
            exact_positive += ticks > 0;
            total_nodes += optimizer.nodes();
            total_prunes += optimizer.prunes();
            exact_out << pair.q << ',' << pair.r << ',' << dec(graph.grid) << ','
                      << dec(graph.total) << ',' << dec(min_mass) << ','
                      << dec(ticks) << ',' << hex8(optimizer.best_body()) << ','
                      << optimizer.nodes() << ',' << optimizer.prunes() << '\n';
            if (!have_min || ticks * min_grid < min_ticks * graph.grid ||
                (ticks * min_grid == min_ticks * graph.grid &&
                 std::tuple{pair.q, pair.r, optimizer.best_body()} <
                     std::tuple{min_pair.q, min_pair.r, min_body})) {
                have_min = true;
                min_ticks = ticks;
                min_grid = graph.grid;
                min_pair = pair;
                min_body = optimizer.best_body();
            }
            if ((index + 1) % 1000 == 0)
                std::cerr << "EXACT_PROGRESS " << (index + 1) << '/'
                          << fallback.size() << " POSITIVE " << exact_positive
                          << '\n';
        }
        demand(exact_out.good(), "exact CSV write failure");
        const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - begin);
        std::cout << "LRC14_PAIR_RANK3_2OVER27_SCREEN_V1\n"
                  << "PAIRS " << pairs.size() << " COARSE_POSITIVE "
                  << coarse_positive << " FALLBACK " << fallback.size()
                  << " EXACT_POSITIVE " << exact_positive << '\n'
                  << "EXACT_NODES " << total_nodes << " EXACT_PRUNES "
                  << total_prunes << '\n';
        if (have_min)
            std::cout << "FALLBACK_NORMALIZED_MIN PAIR " << min_pair.q << ','
                      << min_pair.r << " TICKS " << dec(min_ticks) << " GRID "
                      << dec(min_grid) << " BODY " << hex8(min_body) << '\n';
        std::cout << "ELAPSED_MS " << elapsed.count() << '\n'
                  << "IDENTITY L3(B)=W-SUM_D_I+SUM_D_IJ-SUM_W_IJK\n"
                  << "TARGET 27*L3(B)-2*D>0\n"
                  << "SCOPE FINITE_EXACT_THM4231_REMAINDER_RANK3_NINE_BODY_NO_COFINITE_STRENGTHENING_NO_ENTRY_NO_LRC14\n"
                  << "VERDICT "
                  << ((coarse_positive + exact_positive == pairs.size()) ? "PASS"
                                                                         : "FAIL")
                  << '\n';
        return coarse_positive + exact_positive == pairs.size() ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "PAIR_RANK3_ERROR " << error.what() << '\n';
        return 1;
    }
}
