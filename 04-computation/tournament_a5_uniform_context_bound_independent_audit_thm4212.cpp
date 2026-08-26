#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Clean-room literal audit for THM-4212.
//
// Deliberately absent:
//   * imports from the THM-4208 transfer or response-jet implementation;
//   * ordinal capacity formulas;
//   * the analytic inequalities used in the proof.
//
// Every composite tournament is built as a labelled adjacency matrix.
// Capacities are then rebuilt from the literal odd-directed-path identity
//
//   c_ij(T) = 2 sum_R H(T-V(R)),
//
// and G_+, R_+, and F_5 are evaluated directly.

using u64 = std::uint64_t;
using i128 = __int128_t;

static std::string decimal(i128 x) {
    if (x == 0) return "0";
    const bool negative = x < 0;
    if (negative) x = -x;
    std::string answer;
    while (x != 0) {
        answer.push_back(static_cast<char>('0' + x % 10));
        x /= 10;
    }
    if (negative) answer.push_back('-');
    std::reverse(answer.begin(), answer.end());
    return answer;
}

static std::ostream& operator<<(std::ostream& stream, i128 x) {
    return stream << decimal(x);
}

struct Tournament {
    int order = 0;
    std::vector<u64> out;

    bool arc(int source, int target) const {
        assert(source != target);
        return ((out.at(source) >> target) & 1ULL) != 0;
    }
};

static Tournament labelled_tournament(int order, u64 label) {
    Tournament answer{order, std::vector<u64>(order, 0)};
    int bit = 0;
    for (int i = 0; i < order; ++i) {
        for (int j = i + 1; j < order; ++j, ++bit) {
            if ((label >> bit) & 1ULL) answer.out[i] |= 1ULL << j;
            else answer.out[j] |= 1ULL << i;
        }
    }
    return answer;
}

static Tournament transitive_tournament(int order) {
    Tournament answer{order, std::vector<u64>(order, 0)};
    for (int i = 0; i < order; ++i)
        for (int j = i + 1; j < order; ++j) answer.out[i] |= 1ULL << j;
    return answer;
}

static Tournament cycle3() {
    Tournament answer{3, std::vector<u64>(3, 0)};
    answer.out[0] |= 1ULL << 1;
    answer.out[1] |= 1ULL << 2;
    answer.out[2] |= 1ULL << 0;
    return answer;
}

static Tournament ordinal(const Tournament& left, const Tournament& right) {
    Tournament answer{left.order + right.order,
                      std::vector<u64>(left.order + right.order, 0)};
    for (int i = 0; i < left.order; ++i) {
        answer.out[i] |= left.out[i];
        for (int j = 0; j < right.order; ++j)
            answer.out[i] |= 1ULL << (left.order + j);
    }
    for (int j = 0; j < right.order; ++j)
        answer.out[left.order + j] |= right.out[j] << left.order;
    return answer;
}

static std::string key_of(const Tournament& tournament) {
    std::string answer;
    answer.push_back(static_cast<char>(tournament.order));
    for (u64 row : tournament.out)
        for (int byte = 0; byte < 8; ++byte)
            answer.push_back(static_cast<char>((row >> (8 * byte)) & 255));
    return answer;
}

struct Analysis {
    u64 hamilton = 0;
    u64 mass = 0;
    i128 gate = 0;
};

static Analysis analyze_uncached(const Tournament& tournament) {
    const int n = tournament.order;
    if (n <= 0 || n >= 63) throw std::runtime_error("unsupported tournament order");
    const u64 states = 1ULL << n;
    const u64 full = states - 1;

    std::vector<unsigned char> population(states, 0);
    for (u64 mask = 1; mask < states; ++mask)
        population[mask] = population[mask >> 1] + static_cast<unsigned char>(mask & 1ULL);

    // Hamilton counts for every induced subtournament.
    std::vector<u64> end_dp(states * n, 0), hamilton(states, 0);
    hamilton[0] = 1;
    for (int vertex = 0; vertex < n; ++vertex)
        end_dp[(1ULL << vertex) * n + vertex] = 1;
    for (u64 mask = 1; mask < states; ++mask) {
        for (int last = 0; last < n; ++last) {
            if (((mask >> last) & 1ULL) == 0) continue;
            if (mask != (1ULL << last)) {
                const u64 previous = mask ^ (1ULL << last);
                u64 count = 0;
                for (int before = 0; before < n; ++before)
                    if (((previous >> before) & 1ULL) && tournament.arc(before, last))
                        count += end_dp[previous * n + before];
                end_dp[mask * n + last] = count;
            }
            hamilton[mask] += end_dp[mask * n + last];
        }
    }

    // For each fixed start, count every directed simple path by its vertex
    // set and endpoint.  Even vertex count means odd arc length.
    std::vector<u64> capacity(n * n, 0), path_dp(states * n, 0);
    for (int start = 0; start < n; ++start) {
        std::fill(path_dp.begin(), path_dp.end(), 0);
        path_dp[(1ULL << start) * n + start] = 1;
        for (u64 mask = 1; mask < states; ++mask) {
            if (((mask >> start) & 1ULL) == 0) continue;
            for (int last = 0; last < n; ++last) {
                if (((mask >> last) & 1ULL) == 0) continue;
                if (!(mask == (1ULL << start) && last == start)) {
                    const u64 previous = mask ^ (1ULL << last);
                    if (((previous >> start) & 1ULL) == 0) continue;
                    u64 count = 0;
                    for (int before = 0; before < n; ++before)
                        if (((previous >> before) & 1ULL) && tournament.arc(before, last))
                            count += path_dp[previous * n + before];
                    path_dp[mask * n + last] = count;
                }
                if (last == start || population[mask] % 2 != 0) continue;
                const u64 weighted = path_dp[mask * n + last] * hamilton[full ^ mask];
                const int low = std::min(start, last), high = std::max(start, last);
                capacity[low * n + high] += 2 * weighted;
            }
        }
    }

    struct Edge {
        int source;
        int target;
        u64 weight;
    };
    std::vector<Edge> edges;
    u64 mass = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const u64 weight = capacity[i * n + j];
            mass += weight;
            if (tournament.arc(i, j)) edges.push_back(Edge{i, j, weight});
            else edges.push_back(Edge{j, i, weight});
        }
    }

    // Direct quadratic-kernel definition of G_+.  Each unordered pair of
    // distinct capacity edges contributes according to its incidence type.
    i128 gate = 0;
    for (std::size_t first = 0; first < edges.size(); ++first) {
        for (std::size_t second = first + 1; second < edges.size(); ++second) {
            const Edge& e = edges[first];
            const Edge& f = edges[second];
            const i128 product = i128(e.weight) * f.weight;
            const bool disjoint = e.source != f.source && e.source != f.target &&
                                  e.target != f.source && e.target != f.target;
            if (disjoint) {
                gate += product;
            } else if (e.source == f.source) {
                gate += 4 * product;
            } else if (e.target == f.target) {
                gate -= 4 * product;
            }
        }
    }
    return Analysis{hamilton[full], mass, gate};
}

static std::map<std::string, Analysis> cache;

static const Analysis& analyze(const Tournament& tournament) {
    const std::string key = key_of(tournament);
    auto found = cache.find(key);
    if (found == cache.end())
        found = cache.emplace(key, analyze_uncached(tournament)).first;
    return found->second;
}

static i128 remainder(const Tournament& left, const Tournament& right) {
    const Analysis& a = analyze(left);
    const Analysis& b = analyze(right);
    const Analysis& child = analyze(ordinal(left, right));
    return child.gate - i128(b.hamilton) * b.hamilton * a.gate
           - i128(a.hamilton) * a.hamilton * b.gate;
}

static i128 response(const Tournament& prefix, const Tournament& b,
                     const Tournament& c) {
    return remainder(ordinal(prefix, b), c) - remainder(b, c);
}

static u64 fnv1a_update(u64 hash, const std::string& line) {
    for (unsigned char byte : line) {
        hash ^= byte;
        hash *= 1099511628211ULL;
    }
    hash ^= '\n';
    hash *= 1099511628211ULL;
    return hash;
}

int main() {
    try {
        std::cout << "THM-4212 clean-room literal audit\n";
        std::cout << "engine=labelled adjacency + Hamilton-subset DP + odd-directed-path capacity DP + direct G_+ kernel\n";
        std::cout << "transfer_or_response_jet_imports=NONE\n";
        std::cout << "label_convention=little-endian pairs; bit 1 means i->j\n";

        const Tournament singleton = labelled_tournament(1, 0);
        const Tournament prefix = ordinal(cycle3(), transitive_tournament(5));
        std::vector<std::pair<std::pair<int, u64>, Tournament>> factors;
        for (int order = 1; order <= 3; ++order) {
            const u64 labels = 1ULL << (order * (order - 1) / 2);
            for (u64 label = 0; label < labels; ++label)
                factors.push_back({{order, label}, labelled_tournament(order, label)});
        }

        u64 rows = 0, equality_rows = 0;
        u64 record_hash = 1469598103934665603ULL;
        i128 minimum_slack = -1, minimum_right_slack = -1, minimum_unary_slack = -1;
        for (const auto& named_b : factors) {
            const int order_b = named_b.first.first;
            const u64 label_b = named_b.first.second;
            const Tournament& b = named_b.second;
            if (order_b > 3) continue;
            const Analysis& analysis_b = analyze(b);
            const i128 unary = response(prefix, b, singleton);
            const i128 unary_slack = unary
                - 10764 * i128(analysis_b.hamilton) * analysis_b.hamilton;
            if (unary_slack < 0) throw std::runtime_error("unary inequality failed");
            if (minimum_unary_slack < 0 || unary_slack < minimum_unary_slack)
                minimum_unary_slack = unary_slack;

            for (const auto& named_c : factors) {
                const int order_c = named_c.first.first;
                const u64 label_c = named_c.first.second;
                if (order_b + order_c > 4) continue;
                const Tournament& c = named_c.second;
                const Analysis& analysis_c = analyze(c);
                const i128 actual = response(prefix, b, c);
                const i128 target = 10764 * i128(analysis_b.hamilton) * analysis_b.hamilton
                                     * analysis_c.hamilton * analysis_c.hamilton;
                const i128 slack = actual - target;
                const i128 right_slack = actual
                    - i128(analysis_c.hamilton) * analysis_c.hamilton * unary;
                if (slack < 0) throw std::runtime_error("uniform lower bound failed");
                if (right_slack < 0) throw std::runtime_error("right-context monotonicity failed");
                if (minimum_slack < 0 || slack < minimum_slack) minimum_slack = slack;
                if (minimum_right_slack < 0 || right_slack < minimum_right_slack)
                    minimum_right_slack = right_slack;
                const bool equality = slack == 0;
                if (equality && !(order_b == 1 && order_c == 1))
                    throw std::runtime_error("unexpected equality row");
                equality_rows += equality;
                ++rows;

                std::ostringstream record;
                record << "B=" << order_b << ':' << label_b
                       << ";C=" << order_c << ':' << label_c
                       << ";HB=" << analysis_b.hamilton
                       << ";HC=" << analysis_c.hamilton
                       << ";F=" << actual
                       << ";slack=" << slack
                       << ";right_slack=" << right_slack;
                const std::string line = record.str();
                record_hash = fnv1a_update(record_hash, line);
                std::cout << "literal_record=" << line << '\n';
            }
        }

        if (rows != 25 || equality_rows != 1 || minimum_slack != 0 ||
            minimum_right_slack != 0 || minimum_unary_slack != 0)
            throw std::runtime_error("literal universe summary mismatch");
        if (response(prefix, singleton, singleton) != 10764 ||
            response(prefix, transitive_tournament(2), singleton) != 68364 ||
            response(prefix, singleton, transitive_tournament(2)) != 79128)
            throw std::runtime_error("named signed controls failed");

        std::cout << "literal_universe=all labelled B,C of orders 1..3 with |B|+|C|<=4\n";
        std::cout << "literal_rows=" << rows << '\n';
        std::cout << "literal_equalities=" << equality_rows << '\n';
        std::cout << "minimum_uniform_slack=" << minimum_slack << '\n';
        std::cout << "minimum_right_context_slack=" << minimum_right_slack << '\n';
        std::cout << "minimum_unary_slack=" << minimum_unary_slack << '\n';
        std::cout << "named_controls=F(1,1):10764,F(P2,1):68364,F(1,P2):79128 PASS\n";
        std::cout << "literal_record_fnv64=" << std::hex << std::setw(16)
                  << std::setfill('0') << record_hash << std::dec << std::setfill(' ') << '\n';
        std::cout << "analysis_cache_entries=" << cache.size() << '\n';
        std::cout << "OVERALL=PASS\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
