// Rank-four cubic-majorant census over THM-4231's residual-pair universe.
//
// On a retained wall cell let t be the number of selected body labels that
// fail.  For failure rank at most four,
//
//   14 * 1_{t>0} <= g(t) = 14t - 9 C(t,2) + 3 C(t,3),
//
// because g(0..4)=(0,14,19,18,14).  The increments (14,5,-1,-4)
// decrease, so the induced weighted set function is submodular.  We maximize
// this cubic coverage majorant exactly on eight-label bodies and thereby get
// a rigorous lower bound for the rank-four retained mass without storing a
// quartic tensor.

#define main rank2_wide_hidden_main_rank4allpair
#include "lrc14_rank2_wall_graph_complete_pair_closure_thm4326/independent_cleanroom/rank2_event_sweep_wide_cleanroom.cpp"
#undef main

#include <chrono>

namespace {
constexpr unsigned kRank4Body = 8;
constexpr unsigned kTripleCube = 30 * 30 * 30;

unsigned triple_index(unsigned a, unsigned b, unsigned c) {
    return (a * 30 + b) * 30 + c;
}

struct RankCell4 {
    u32 failures = 0;
    i128 width = 0;
};

struct Graph4Majorant {
    i128 grid = 0;
    i128 retained = 0;
    std::array<i128, 5> rank_mass{};
    std::array<u64, 5> rank_cells{};
    std::array<i128, 30> degree{};
    std::array<std::array<i128, 30>, 30> codegree{};
    std::vector<i128> triple_degree = std::vector<i128>(kTripleCube, 0);
    std::vector<RankCell4> cells;
};

i128 triple_incidence(const Graph4Majorant& graph, unsigned a, unsigned b,
                     unsigned c) {
    return graph.triple_degree[triple_index(a, b, c)];
}

Graph4Majorant event_sweep_graph4_majorant(
    int q, int r, const std::vector<BaseEvent>& base_events) {
    Graph4Majorant graph;
    graph.grid = exact_lcm(kBaseGrid, 14LL * q);
    graph.grid = exact_lcm(graph.grid, 14LL * r);
    demand(graph.grid % kBaseGrid == 0, "rank4 grid lost base divisibility");
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
        demand(right > left, "rank4 non-increasing merged wall sequence");

        if (q_safe && r_safe) {
            const unsigned rank = std::popcount(failures);
            if (rank <= 4) {
                const i128 width = right - left;
                graph.retained += width;
                graph.rank_mass[rank] += width;
                ++graph.rank_cells[rank];
                graph.cells.push_back({failures, width});

                std::array<unsigned, 4> vertices{};
                unsigned count = 0;
                u32 copy = failures;
                while (copy != 0) {
                    vertices[count++] = std::countr_zero(copy);
                    copy &= copy - 1;
                }
                demand(count == rank, "rank4 failure extraction changed");
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
                                graph.triple_degree[triple_index(
                                    p[0], p[1], p[2])] += width;
                        }
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
    demand(pool_cursor == base_events.size(),
           "rank4 not all pool events consumed");
    demand(!q_safe && !r_safe, "rank4 pair state did not close unsafe");
    demand(failures == (u32{1} << 30) - 1,
           "rank4 pool state did not close unsafe");
    demand(graph.retained == std::accumulate(graph.rank_mass.begin(),
                                             graph.rank_mass.end(), i128{0}),
           "rank4 retained partition failed");
    return graph;
}

i128 majorant_reward_direct(const Graph4Majorant& graph, u32 body) {
    demand(std::popcount(body) == static_cast<int>(kRank4Body),
           "majorant direct body rank changed");
    constexpr std::array<int, 5> g = {0, 14, 19, 18, 14};
    i128 reward = 0;
    for (const RankCell4& cell : graph.cells) {
        const unsigned hits = std::popcount(cell.failures & body);
        demand(hits <= 4, "majorant cell gained impossible hit count");
        reward += i128(g[hits]) * cell.width;
    }
    return reward;
}

i128 retained_mass_direct(const Graph4Majorant& graph, u32 body) {
    demand(std::popcount(body) == static_cast<int>(kRank4Body),
           "retained direct body rank changed");
    i128 mass = 0;
    for (const RankCell4& cell : graph.cells)
        if ((cell.failures & body) == 0) mass += cell.width;
    return mass;
}

i128 majorant_reward_tensor(const Graph4Majorant& graph, u32 body) {
    demand(std::popcount(body) == static_cast<int>(kRank4Body),
           "majorant tensor body rank changed");
    std::array<unsigned, kRank4Body> vertices{};
    unsigned count = 0;
    u32 copy = body;
    i128 reward = 0;
    while (copy != 0) {
        const unsigned v = std::countr_zero(copy);
        vertices[count++] = v;
        reward += 14 * graph.degree[v];
        copy &= copy - 1;
    }
    for (unsigned i = 0; i < count; ++i)
        for (unsigned j = i + 1; j < count; ++j)
            reward -= 9 * graph.codegree[vertices[i]][vertices[j]];
    for (unsigned i = 0; i < count; ++i)
        for (unsigned j = i + 1; j < count; ++j)
            for (unsigned k = j + 1; k < count; ++k)
                reward += 3 * triple_incidence(
                                  graph, vertices[i], vertices[j], vertices[k]);
    return reward;
}

struct CubicMajorantOptimizer {
    explicit CubicMajorantOptimizer(const Graph4Majorant& graph)
        : graph_(graph) {
        std::iota(order_.begin(), order_.end(), 0);
        std::sort(order_.begin(), order_.end(), [&](unsigned a, unsigned b) {
            if (graph_.degree[a] != graph_.degree[b])
                return graph_.degree[a] > graph_.degree[b];
            return a < b;
        });
        seed_greedy();
    }

    void solve() { descend(0, kRank4Body, 0, 0); }
    i128 best_reward() const { return best_reward_; }
    u32 best_body() const { return best_body_; }
    u64 nodes() const { return nodes_; }
    u64 prunes() const { return prunes_; }

  private:
    const Graph4Majorant& graph_;
    std::array<unsigned, 30> order_{};
    i128 best_reward_ = -(i128{1} << 120);
    u32 best_body_ = 0;
    u64 nodes_ = 0;
    u64 prunes_ = 0;

    i128 marginal(unsigned vertex, u32 selected) const {
        std::array<unsigned, kRank4Body> prior{};
        unsigned count = 0;
        u32 copy = selected;
        i128 gain = 14 * graph_.degree[vertex];
        while (copy != 0) {
            const unsigned u = std::countr_zero(copy);
            prior[count++] = u;
            gain -= 9 * graph_.codegree[vertex][u];
            copy &= copy - 1;
        }
        for (unsigned i = 0; i < count; ++i)
            for (unsigned j = i + 1; j < count; ++j)
                gain += 3 * triple_incidence(
                                graph_, vertex, prior[i], prior[j]);
        return gain;
    }

    void consider(u32 body, i128 reward) {
        demand(std::popcount(body) == static_cast<int>(kRank4Body),
               "majorant optimizer leaf rank changed");
        demand(reward == majorant_reward_tensor(graph_, body),
               "majorant incremental/tensor reward drift");
        demand(reward == majorant_reward_direct(graph_, body),
               "majorant tensor/direct reward drift");
        if (reward > best_reward_ ||
            (reward == best_reward_ && body < best_body_)) {
            best_reward_ = reward;
            best_body_ = body;
        }
    }

    void seed_greedy() {
        u32 selected = 0;
        i128 reward = 0;
        for (unsigned step = 0; step < kRank4Body; ++step) {
            unsigned chosen = 30;
            i128 chosen_gain = -(i128{1} << 120);
            for (unsigned vertex = 0; vertex < 30; ++vertex) {
                if ((selected >> vertex) & 1U) continue;
                const i128 gain = marginal(vertex, selected);
                if (chosen == 30 || gain > chosen_gain ||
                    (gain == chosen_gain && vertex < chosen)) {
                    chosen = vertex;
                    chosen_gain = gain;
                }
            }
            demand(chosen < 30, "majorant greedy exhausted vertices");
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
        demand(std::popcount(selected) == static_cast<int>(kRank4Body) &&
                   need == 30 - position,
               "majorant forced completion rank failure");
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

        std::array<i128, 30> gains{};
        unsigned count = 0;
        for (unsigned cursor = position; cursor < 30; ++cursor)
            gains[count++] = marginal(order_[cursor], selected);
        std::nth_element(gains.begin(), gains.begin() + need,
                         gains.begin() + count, std::greater<i128>());
        i128 upper = reward;
        for (unsigned i = 0; i < need; ++i) upper += gains[i];
        // Each cell contribution has decreasing marginals 14,5,-1,-4.
        // Therefore the objective is submodular, and the sum of the largest
        // current marginals is a sound upper bound for every completion.
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

struct CoarseMajorant {
    i128 lower14 = 0;
    i128 ticks = 0;
    u32 top8 = 0;
};

CoarseMajorant coarse_majorant(const Graph4Majorant& graph) {
    std::array<unsigned, 30> order{};
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](unsigned a, unsigned b) {
        if (graph.degree[a] != graph.degree[b])
            return graph.degree[a] > graph.degree[b];
        return a < b;
    });
    CoarseMajorant result;
    result.lower14 = 14 * graph.retained;
    for (unsigned i = 0; i < kRank4Body; ++i) {
        result.lower14 -= 14 * graph.degree[order[i]];
        result.top8 |= u32{1} << order[i];
    }
    result.ticks = 81 * result.lower14 - 98 * graph.grid;
    return result;
}

std::string body_labels_rank4(u32 body) {
    std::ostringstream out;
    bool first = true;
    for (unsigned bit = 0; bit < 30; ++bit) {
        if (((body >> bit) & 1U) == 0) continue;
        if (!first) out << ';';
        first = false;
        out << kPool[bit];
    }
    return out.str();
}
}  // namespace

int main(int argc, char** argv) {
    try {
        demand(argc >= 4 && argc <= 7,
               "usage: rank4_majorant PAIRS LEDGER_CSV SUMMARY_OUT "
               "[--exact-all] [START COUNT]");
        u64 input_fnv = 0;
        const std::vector<Pair> pairs = read_pairs(argv[1], input_fnv);
        const std::vector<BaseEvent> base_events = make_base_events();
        bool exact_all = false;
        std::size_t start = 0;
        std::size_t count = pairs.size();
        int cursor = 4;
        if (cursor < argc && std::string(argv[cursor]) == "--exact-all") {
            exact_all = true;
            ++cursor;
        }
        demand(argc - cursor == 0 || argc - cursor == 2,
               "rank4 optional arguments must be --exact-all and/or START COUNT");
        if (argc - cursor == 2) {
            start = static_cast<std::size_t>(std::stoull(argv[cursor]));
            count = static_cast<std::size_t>(std::stoull(argv[cursor + 1]));
            demand(start <= pairs.size() && count <= pairs.size() - start,
                   "requested pair slice outside universe");
        }

        std::ofstream ledger(argv[2]);
        std::ofstream summary(argv[3]);
        demand(bool(ledger) && bool(summary), "cannot create rank4 outputs");
        ledger << "index,q,r,grid,rank0,rank1,rank2,rank3,rank4,retained," \
                  "coarse_lower14,coarse_ticks,certificate_type," \
                  "certificate_lower14,certificate_ticks,certificate_body," \
                  "certificate_labels,direct_l4,nodes,prunes\n";

        u64 degree_positive = 0;
        u64 degree_certificates = 0;
        u64 cubic_positive = 0;
        u64 cubic_nonpositive = 0;
        u64 total_nodes = 0;
        u64 total_prunes = 0;
        const auto begin = std::chrono::steady_clock::now();
        for (std::size_t offset = 0; offset < count; ++offset) {
            const std::size_t index = start + offset;
            const Pair pair = pairs[index];
            const Graph4Majorant graph = event_sweep_graph4_majorant(
                pair.q, pair.r, base_events);
            const CoarseMajorant coarse = coarse_majorant(graph);
            degree_positive += coarse.ticks > 0;

            std::string certificate_type = "DEGREE";
            i128 certificate_lower14 = coarse.lower14;
            i128 certificate_ticks = coarse.ticks;
            u32 certificate_body = coarse.top8;
            u64 nodes = 0, prunes = 0;
            if (coarse.ticks > 0 && !exact_all) {
                ++degree_certificates;
            } else {
                certificate_type = "CUBIC_EXACT";
                CubicMajorantOptimizer optimizer(graph);
                optimizer.solve();
                certificate_lower14 =
                    14 * graph.retained - optimizer.best_reward();
                certificate_ticks =
                    81 * certificate_lower14 - 98 * graph.grid;
                certificate_body = optimizer.best_body();
                nodes = optimizer.nodes();
                prunes = optimizer.prunes();
                total_nodes += nodes;
                total_prunes += prunes;
                if (certificate_ticks > 0)
                    ++cubic_positive;
                else
                    ++cubic_nonpositive;
            }
            const i128 direct_l4 =
                retained_mass_direct(graph, certificate_body);
            demand(14 * direct_l4 >= certificate_lower14,
                   "cubic majorant lower bound exceeds direct L4 control");
            ledger << index << ',' << pair.q << ',' << pair.r << ','
                   << dec(graph.grid);
            for (i128 mass : graph.rank_mass) ledger << ',' << dec(mass);
            ledger << ',' << dec(graph.retained) << ','
                   << dec(coarse.lower14) << ',' << dec(coarse.ticks) << ','
                   << certificate_type << ',' << dec(certificate_lower14)
                   << ',' << dec(certificate_ticks) << ','
                   << hex8(certificate_body) << ','
                   << body_labels_rank4(certificate_body) << ','
                   << dec(direct_l4) << ',' << nodes << ',' << prunes << '\n';

            if ((offset + 1) % 10000 == 0)
                std::cerr << "RANK4_PROGRESS " << (offset + 1) << '/'
                          << count << " DEGREE_PASS " << degree_positive
                          << " DEGREE_CERT " << degree_certificates
                          << " CUBIC_PASS " << cubic_positive
                          << " CUBIC_FAIL " << cubic_nonpositive << '\n';
        }
        demand(ledger.good(), "rank4 ledger write failed");
        const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - begin);
        summary << "LRC14_RANK4_CUBIC_MAJORANT_ALLPAIR_PROBE_V1\n"
                << "INPUT_PAIRS " << pairs.size()
                << " INPUT_NORMALIZED_FNV " << fnv_hex(input_fnv)
                << " BASE_EVENT_COORDINATES " << base_events.size() << '\n'
                << "SLICE_START " << start << " SLICE_COUNT " << count
                << " EXACT_ALL " << (exact_all ? 1 : 0) << '\n'
                << "MAJORANT_VALUES g(0..4)=0,14,19,18,14 "
                   "MARGINALS=14,5,-1,-4\n"
                << "TARGET 81*(14W4-MAX_G)-98D>0\n"
                << "DEGREE_POSITIVE " << degree_positive
                << " DEGREE_CERTIFICATES " << degree_certificates
                << " CUBIC_EXACT_POSITIVE " << cubic_positive
                << " CUBIC_EXACT_NONPOSITIVE " << cubic_nonpositive << '\n'
                << "CUBIC_NODES " << total_nodes << " CUBIC_PRUNES "
                << total_prunes << '\n'
                << "ELAPSED_MS " << elapsed.count() << '\n'
                << "SCOPE FINITE_EXACT_REQUESTED_SLICE_THM4231_REMAINDER_"
                   "RANK4_CUBIC_LOWER_BOUND_NOT_FULL_SAFE_MASS_NO_ENTRY_NO_LRC14\n"
                << "VERDICT " << (cubic_nonpositive == 0 ? "PASS" : "INCONCLUSIVE")
                << '\n';
        demand(summary.good(), "rank4 summary write failed");
        return cubic_nonpositive == 0 ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "RANK4_MAJORANT_ERROR " << error.what() << '\n';
        return 1;
    }
}
