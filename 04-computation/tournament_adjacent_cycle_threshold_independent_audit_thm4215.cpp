#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Clean-room literal audit for THM-4215.
//
// Deliberately absent:
//   * response jets and endpoint-energy inequalities;
//   * ordinal capacity-transfer formulas;
//   * imports from any tournament computation in the repository.
//
// Composite tournaments are first built only as labelled adjacency matrices.
// Hamilton counts, odd-directed-path capacities, G_+, R_+, and the contextual
// response are then reconstructed directly from subset dynamic programming.

using u64 = std::uint64_t;
using i128 = __int128_t;

static std::string decimal(i128 value) {
    if (value == 0) return "0";
    const bool negative = value < 0;
    if (negative) value = -value;
    std::string answer;
    while (value != 0) {
        answer.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (negative) answer.push_back('-');
    std::reverse(answer.begin(), answer.end());
    return answer;
}

static std::ostream& operator<<(std::ostream& stream, i128 value) {
    return stream << decimal(value);
}

struct Tournament {
    int order = 0;
    std::vector<u64> out;

    bool arc(int source, int target) const {
        assert(source != target);
        return ((out.at(source) >> target) & 1ULL) != 0;
    }
};

static Tournament transitive(int order) {
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

static Tournament adjacent_cycle_prefix(int tail) {
    Tournament answer = ordinal(cycle3(), cycle3());
    if (tail > 0) answer = ordinal(answer, transitive(tail));
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
    if (n <= 0 || n > 20) throw std::runtime_error("unsupported literal order");
    const u64 states = 1ULL << n;
    const u64 full = states - 1;

    std::vector<unsigned char> population(states, 0);
    for (u64 mask = 1; mask < states; ++mask)
        population[mask] = static_cast<unsigned char>(
            population[mask >> 1] + static_cast<unsigned char>(mask & 1ULL));

    // Hamilton paths of every induced subtournament, sorted by final vertex.
    std::vector<u64> end_dp(states * static_cast<u64>(n), 0);
    std::vector<u64> hamilton(states, 0);
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

    // For each fixed start, count directed simple paths by vertex set and
    // endpoint.  Even vertex count is odd arc length.  Weight the exposed
    // path by the Hamilton count of its complement.
    std::vector<u64> capacity(static_cast<u64>(n) * n, 0);
    std::vector<u64> path_dp(states * static_cast<u64>(n), 0);
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
                        if (((previous >> before) & 1ULL)
                            && tournament.arc(before, last))
                            count += path_dp[previous * n + before];
                    path_dp[mask * n + last] = count;
                }
                if (last == start || population[mask] % 2 != 0) continue;
                const u64 weighted = path_dp[mask * n + last] * hamilton[full ^ mask];
                const int low = std::min(start, last);
                const int high = std::max(start, last);
                capacity[static_cast<u64>(low) * n + high] += 2 * weighted;
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
            const u64 weight = capacity[static_cast<u64>(i) * n + j];
            mass += weight;
            if (tournament.arc(i, j)) edges.push_back(Edge{i, j, weight});
            else edges.push_back(Edge{j, i, weight});
        }
    }

    // Literal incidence-kernel definition of G_+.
    i128 gate = 0;
    for (std::size_t first = 0; first < edges.size(); ++first) {
        for (std::size_t second = first + 1; second < edges.size(); ++second) {
            const Edge& e = edges[first];
            const Edge& f = edges[second];
            const i128 product = i128(e.weight) * f.weight;
            const bool disjoint = e.source != f.source && e.source != f.target
                                  && e.target != f.source && e.target != f.target;
            if (disjoint) gate += product;
            else if (e.source == f.source) gate += 4 * product;
            else if (e.target == f.target) gate -= 4 * product;
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

static i128 response(const Tournament& prefix, const Tournament& middle,
                     const Tournament& right) {
    return remainder(ordinal(prefix, middle), right) - remainder(middle, right);
}

static i128 singleton_formula(int tail) {
    const i128 two = i128(1) << tail;
    const i128 four = two * two;
    return 108 * (2 * four - (12 * tail + 102) * two + 1);
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
        std::cout << "THM-4215 clean-room literal audit\n";
        std::cout << "engine=labelled adjacency + Hamilton-subset DP + odd-directed-path capacity DP + direct G_+ kernel\n";
        std::cout << "transfer_or_response_jet_imports=NONE\n";
        const Tournament one = transitive(1);
        const Tournament p2 = transitive(2);
        const Tournament cycle = cycle3();
        u64 record_hash = 1469598103934665603ULL;

        for (int tail = 0; tail <= 7; ++tail) {
            const Tournament prefix = adjacent_cycle_prefix(tail);
            const i128 actual = response(prefix, one, one);
            const i128 expected = singleton_formula(tail);
            if (actual != expected) throw std::runtime_error("singleton formula mismatch");
            if (tail <= 6 && actual >= 0)
                throw std::runtime_error("missing negative singleton obstruction");
            if (tail == 7 && actual != 967788)
                throw std::runtime_error("wrong positive crossing value");
            std::ostringstream record;
            record << "tail=" << tail << ";F11=" << actual;
            const std::string line = record.str();
            std::cout << "boundary_record=" << line << '\n';
            record_hash = fnv1a_update(record_hash, line);
        }

        const Tournament prefix = adjacent_cycle_prefix(7);
        const Analysis& prefix_analysis = analyze(prefix);
        if (prefix.order != 13 || prefix_analysis.hamilton != 9
            || prefix_analysis.mass != 18414 || prefix_analysis.gate != -734004)
            throw std::runtime_error("Z7 literal prefix invariant mismatch");

        struct Control {
            const char* name;
            Tournament middle;
            Tournament right;
            i128 expected;
        };
        const std::vector<Control> controls = {
            Control{"P1,P1", one, one, 967788},
            Control{"P2,P1", p2, one, 8681580},
            Control{"P1,P2", one, p2, 9649368},
            Control{"C3,P1", cycle, one, 81113472},
            Control{"P1,C3", one, cycle, 58787304432LL},
        };
        u64 equality_rows = 0;
        for (const Control& control : controls) {
            const Analysis& b = analyze(control.middle);
            const Analysis& c = analyze(control.right);
            const i128 actual = response(prefix, control.middle, control.right);
            if (actual != control.expected)
                throw std::runtime_error("named literal response mismatch");
            const i128 floor = 967788 * i128(b.hamilton) * b.hamilton
                               * c.hamilton * c.hamilton;
            const i128 slack = actual - floor;
            if (slack < 0) throw std::runtime_error("named floor failure");
            const bool equality = slack == 0;
            if (equality && std::string(control.name) != "P1,P1")
                throw std::runtime_error("unexpected named equality");
            equality_rows += equality;
            std::ostringstream record;
            record << "contexts=" << control.name << ";HB=" << b.hamilton
                   << ";HC=" << c.hamilton << ";F=" << actual
                   << ";slack=" << slack;
            const std::string line = record.str();
            std::cout << "context_record=" << line << '\n';
            record_hash = fnv1a_update(record_hash, line);
        }
        if (equality_rows != 1) throw std::runtime_error("named equality summary");

        std::cout << "literal_boundary_rows=8\n";
        std::cout << "literal_context_controls=5\n";
        std::cout << "literal_nontransitive_controls=C3,P1 and P1,C3\n";
        std::cout << "Z7_invariants=order:13,H:9,W:18414,G_plus:-734004\n";
        std::cout << "literal_record_fnv64=" << std::hex << std::setw(16)
                  << std::setfill('0') << record_hash << std::dec << std::setfill(' ')
                  << '\n';
        std::cout << "analysis_cache_entries=" << cache.size() << '\n';
        std::cout << "OVERALL=PASS\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
