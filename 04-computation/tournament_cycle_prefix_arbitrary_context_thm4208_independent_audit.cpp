#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

// Clean-room audit for the proposed THM-4208 package.
//
// Deliberately absent:
//   * imports from any tournament transfer/capacity module;
//   * ordinal-sum capacity formulas;
//   * fixed-response jets or the claimed C-finite formula.
//
// Capacities are rebuilt from the literal odd-directed-path definition
//   c_ij = 2 sum_R H(T - V(R)),
// and G_+, R_+, and every composite tournament are then evaluated directly.

using u64 = std::uint64_t;
using i128 = __int128_t;

static std::string dec(i128 x) {
    if (x == 0) return "0";
    const bool neg = x < 0;
    if (neg) x = -x;
    std::string s;
    while (x != 0) {
        s.push_back(static_cast<char>('0' + x % 10));
        x /= 10;
    }
    if (neg) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}

static std::ostream& operator<<(std::ostream& os, i128 x) {
    return os << dec(x);
}

struct Tournament {
    int n = 0;
    std::vector<u64> out;

    bool arc(int i, int j) const {
        assert(i != j);
        return ((out.at(i) >> j) & 1ULL) != 0;
    }
};

static Tournament labelled_tournament(int n, u64 label) {
    Tournament t{n, std::vector<u64>(n, 0)};
    int bit = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j, ++bit) {
            if ((label >> bit) & 1ULL) t.out[i] |= 1ULL << j;
            else t.out[j] |= 1ULL << i;
        }
    }
    return t;
}

static Tournament transitive_tournament(int n) {
    Tournament t{n, std::vector<u64>(n, 0)};
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            t.out[i] |= 1ULL << j;
    return t;
}

static Tournament cycle3() {
    Tournament t{3, std::vector<u64>(3, 0)};
    t.out[0] |= 1ULL << 1;
    t.out[1] |= 1ULL << 2;
    t.out[2] |= 1ULL << 0;
    return t;
}

static Tournament ordinal(const Tournament& a, const Tournament& b) {
    Tournament t{a.n + b.n, std::vector<u64>(a.n + b.n, 0)};
    for (int i = 0; i < a.n; ++i) {
        t.out[i] |= a.out[i];
        for (int j = 0; j < b.n; ++j) t.out[i] |= 1ULL << (a.n + j);
    }
    for (int j = 0; j < b.n; ++j) t.out[a.n + j] |= b.out[j] << a.n;
    return t;
}

static Tournament A(int n) {
    // A_n = C3 ▷ P_n, with P_0 used only as empty-prefix notation.
    return n == 0 ? cycle3() : ordinal(cycle3(), transitive_tournament(n));
}

static std::string key_of(const Tournament& t) {
    std::string s;
    s.reserve(2 + 8 * t.n);
    s.push_back(static_cast<char>(t.n));
    for (u64 x : t.out)
        for (int k = 0; k < 8; ++k) s.push_back(static_cast<char>((x >> (8 * k)) & 255));
    return s;
}

static bool is_path(const Tournament& t, const std::vector<int>& p) {
    for (std::size_t i = 1; i < p.size(); ++i)
        if (!t.arc(p[i - 1], p[i])) return false;
    return true;
}

static bool is_transitive(const Tournament& t) {
    // A tournament is transitive iff it contains no directed triangle.
    for (int i = 0; i < t.n; ++i)
        for (int j = i + 1; j < t.n; ++j)
            for (int k = j + 1; k < t.n; ++k) {
                const bool cyc1 = t.arc(i, j) && t.arc(j, k) && t.arc(k, i);
                const bool cyc2 = t.arc(j, i) && t.arc(i, k) && t.arc(k, j);
                if (cyc1 || cyc2) return false;
            }
    return true;
}

struct Analysis {
    u64 h = 0;
    u64 w = 0;
    i128 g = 0;
    std::vector<u64> vplus;
    std::vector<u64> cap;
};

static Analysis analyze_uncached(const Tournament& t) {
    if (t.n <= 0 || t.n >= 63) throw std::runtime_error("unsupported order");
    const int n = t.n;
    const u64 states = 1ULL << n;
    const u64 full = states - 1;

    std::vector<unsigned char> pop(states, 0);
    for (u64 mask = 1; mask < states; ++mask) pop[mask] = pop[mask >> 1] + (mask & 1ULL);

    // H[mask] is rebuilt by the standard last-vertex subset recurrence.
    std::vector<u64> hp(states * n, 0), h(states, 0);
    h[0] = 1;
    for (int v = 0; v < n; ++v) hp[(1ULL << v) * n + v] = 1;
    for (u64 mask = 1; mask < states; ++mask) {
        for (int last = 0; last < n; ++last) {
            if (((mask >> last) & 1ULL) == 0) continue;
            if (mask != (1ULL << last)) {
                const u64 prev = mask ^ (1ULL << last);
                u64 total = 0;
                for (int u = 0; u < n; ++u)
                    if (((prev >> u) & 1ULL) && t.arc(u, last)) total += hp[prev * n + u];
                hp[mask * n + last] = total;
            }
            h[mask] += hp[mask * n + last];
        }
    }

    // Rebuild rooted path-cover counts and symmetric capacities.  A path with
    // an even number of vertices has odd arc length.
    std::vector<u64> vplus(n, 0), cap(n * n, 0);
    std::vector<u64> dp(states * n, 0);
    for (int start = 0; start < n; ++start) {
        std::fill(dp.begin(), dp.end(), 0);
        dp[(1ULL << start) * n + start] = 1;
        for (u64 mask = 1; mask < states; ++mask) {
            if (((mask >> start) & 1ULL) == 0) continue;
            for (int last = 0; last < n; ++last) {
                if (((mask >> last) & 1ULL) == 0) continue;
                if (!(mask == (1ULL << start) && last == start)) {
                    const u64 prev = mask ^ (1ULL << last);
                    if (((prev >> start) & 1ULL) == 0) continue;
                    u64 total = 0;
                    for (int u = 0; u < n; ++u)
                        if (((prev >> u) & 1ULL) && t.arc(u, last)) total += dp[prev * n + u];
                    dp[mask * n + last] = total;
                }
                const u64 count = dp[mask * n + last];
                if (count == 0) continue;
                const u64 weighted = count * h[full ^ mask];
                vplus[last] += weighted;
                if (last != start && (pop[mask] % 2 == 0)) {
                    const int lo = std::min(start, last), hi = std::max(start, last);
                    cap[lo * n + hi] += 2 * weighted;
                }
            }
        }
    }

    u64 w = 0;
    i128 sum_edge_sq = 0, sum_degree_sq = 0, current = 0;
    std::vector<u64> degree(n, 0);
    std::vector<long long> signed_degree(n, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const u64 c = cap[i * n + j];
            w += c;
            sum_edge_sq += i128(c) * c;
            degree[i] += c;
            degree[j] += c;
            if (t.arc(i, j)) {
                signed_degree[i] += static_cast<long long>(c);
                signed_degree[j] -= static_cast<long long>(c);
            } else {
                signed_degree[i] -= static_cast<long long>(c);
                signed_degree[j] += static_cast<long long>(c);
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        sum_degree_sq += i128(degree[i]) * degree[i];
        current += i128(degree[i]) * signed_degree[i];
    }
    // Evaluate D literally from unordered pairs of vertex-disjoint edges.
    // Also check the quadratic compression, but do not use it as the audit's
    // definition of D.
    i128 disjoint_literal = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            for (int k = 0; k < n; ++k)
                for (int ell = k + 1; ell < n; ++ell) {
                    if (i == k && j == ell) continue;
                    if (!(i < k || (i == k && j < ell))) continue;
                    if (i != k && i != ell && j != k && j != ell)
                        disjoint_literal += i128(cap[i * n + j]) * cap[k * n + ell];
                }
    const i128 disjoint_compressed = (i128(w) * w + sum_edge_sq - sum_degree_sq) / 2;
    if (disjoint_literal != disjoint_compressed)
        throw std::runtime_error("literal disjoint-edge gate mismatch");
    const i128 g = disjoint_literal + 2 * current;

    const u64 vtotal = std::accumulate(vplus.begin(), vplus.end(), u64(0));
    if (vtotal != w + h[full]) throw std::runtime_error("path-cover mass check failed");
    return Analysis{h[full], w, g, std::move(vplus), std::move(cap)};
}

static std::map<std::string, Analysis> analysis_cache;

static const Analysis& analyze(const Tournament& t) {
    const std::string key = key_of(t);
    auto it = analysis_cache.find(key);
    if (it == analysis_cache.end()) it = analysis_cache.emplace(key, analyze_uncached(t)).first;
    return it->second;
}

static i128 Rplus(const Tournament& x, const Tournament& y) {
    const Analysis& ax = analyze(x);
    const Analysis& ay = analyze(y);
    const Analysis& axy = analyze(ordinal(x, y));
    return axy.g - i128(ay.h) * ay.h * ax.g - i128(ax.h) * ax.h * ay.g;
}

static i128 F(int n, const Tournament& b, const Tournament& c) {
    const Tournament an_b = ordinal(A(n), b);
    return Rplus(an_b, c) - Rplus(b, c);
}

static std::vector<int> insert_canonically(const Tournament& t, int x,
                                            const std::vector<int>& q) {
    std::size_t at = 0;
    while (at < q.size() && !t.arc(x, q[at])) ++at;
    std::vector<int> out = q;
    out.insert(out.begin() + static_cast<std::ptrdiff_t>(at), x);
    if (!is_path(t, out)) throw std::runtime_error("Hamilton insertion failed");
    return out;
}

static std::string seq_key(const std::vector<int>& p) {
    std::string s;
    for (int v : p) s.push_back(static_cast<char>('a' + v));
    return s;
}

static u64 fnv1a_update(u64 h, const std::string& s) {
    for (unsigned char c : s) {
        h ^= c;
        h *= 1099511628211ULL;
    }
    h ^= '\n';
    h *= 1099511628211ULL;
    return h;
}

static i128 pow2(int n) { return i128(1) << n; }
static i128 pow4(int n) { return i128(1) << (2 * n); }

int main() {
    try {
        std::cout << "THM-4208 clean-room independent replay\n";
        std::cout << "engine=literal labelled adjacency + Hamilton-subset DP + odd-path capacity DP\n";
        std::cout << "primary_transfer_imports=NONE\n";
        std::cout << "label_convention=little-endian pairs (0,1),(0,2),...; bit 1 means i->j\n";

        // -----------------------------------------------------------------
        // 1. Explicit endpoint injection and its moment consequence.
        // -----------------------------------------------------------------
        u64 tournament_count = 0, vertex_inequalities = 0;
        u64 injection_objects = 0, injection_images = 0;
        u64 literal_capacity_coordinates = 0, literal_capacity_objects = 0;
        u64 empty_first_component_objects = 0;
        u64 transitive_count = 0, equality_count = 0;
        i128 aggregate_E = 0;
        u64 moment_digest = 1469598103934665603ULL;

        for (int n = 1; n <= 5; ++n) {
            const int edges = n * (n - 1) / 2;
            const u64 labels = 1ULL << edges;
            for (u64 label = 0; label < labels; ++label) {
                const Tournament t = labelled_tournament(n, label);
                const Analysis& z = analyze(t);
                ++tournament_count;
                empty_first_component_objects += z.h;

                std::vector<u64> domain(n, 0);
                std::vector<u64> literal_cap(n * n, 0);
                std::vector<std::unordered_set<std::string>> image(n);
                std::vector<int> perm(n);
                std::iota(perm.begin(), perm.end(), 0);
                u64 literal_h = 0;
                do {
                    literal_h += is_path(t, perm);
                    for (int split = 1; split <= n; ++split) {
                        const std::vector<int> p(perm.begin(), perm.begin() + split);
                        const std::vector<int> q(perm.begin() + split, perm.end());
                        if (!is_path(t, p) || !is_path(t, q)) continue;
                        const int i = p.back();
                        ++domain[i];
                        if (p.size() % 2 == 0) {
                            const int lo = std::min(p.front(), p.back());
                            const int hi = std::max(p.front(), p.back());
                            literal_cap[lo * n + hi] += 2;
                            ++literal_capacity_objects;
                        }
                        const std::vector<int> q_with_i = insert_canonically(t, i, q);
                        std::string target;
                        if (p.size() == 1) {
                            if (q_with_i.size() != static_cast<std::size_t>(n))
                                throw std::runtime_error("bad Hamilton target");
                            target = "H|" + seq_key(q_with_i);
                        } else {
                            const int j = p[p.size() - 2];
                            if (!t.arc(j, i)) throw std::runtime_error("bad terminal arc");
                            const std::vector<int> p_short(p.begin(), p.end() - 1);
                            if (!is_path(t, p_short)) throw std::runtime_error("bad shortened path");
                            target = "V" + std::to_string(j) + "|" + seq_key(p_short) + "|" + seq_key(q_with_i);
                        }
                        if (!image[i].insert(target).second)
                            throw std::runtime_error("endpoint insertion is not injective");
                    }
                } while (std::next_permutation(perm.begin(), perm.end()));

                if (literal_h != z.h) throw std::runtime_error("literal Hamilton count mismatch");
                for (int i = 0; i < n; ++i)
                    for (int j = i + 1; j < n; ++j) {
                        if (literal_cap[i * n + j] != z.cap[i * n + j])
                            throw std::runtime_error("literal capacity mismatch");
                        ++literal_capacity_coordinates;
                    }

                i128 m = 0, weighted_slack = 0;
                bool all_slacks_zero = true;
                for (int i = 0; i < n; ++i) {
                    if (domain[i] != z.vplus[i] || image[i].size() != domain[i])
                        throw std::runtime_error("explicit injection/domain mismatch");
                    u64 codomain = z.h;
                    for (int j = 0; j < n; ++j) if (j != i && t.arc(j, i)) codomain += z.vplus[j];
                    if (domain[i] > codomain) throw std::runtime_error("endpoint inequality failed");
                    const u64 slack = codomain - domain[i];
                    all_slacks_zero = all_slacks_zero && (slack == 0);
                    weighted_slack += i128(z.vplus[i]) * slack;
                    m += i128(z.vplus[i]) * z.vplus[i];
                    ++vertex_inequalities;
                    injection_objects += domain[i];
                    injection_images += image[i].size();
                }
                const i128 aa = i128(z.w) + z.h;
                const i128 E = aa * (aa + 2 * z.h) - 3 * m;
                if (E < 0 || E != 2 * weighted_slack)
                    throw std::runtime_error("endpoint energy identity failed");
                const bool trans = is_transitive(t);
                const bool equal = (E == 0);
                if (equal != trans || equal != all_slacks_zero)
                    throw std::runtime_error("endpoint equality classification failed");
                transitive_count += trans;
                equality_count += equal;
                aggregate_E += E;

                std::ostringstream record;
                record << n << ':' << label << ':' << z.h << ':' << z.w << ':' << dec(E) << ':';
                for (u64 v : z.vplus) record << v << ',';
                moment_digest = fnv1a_update(moment_digest, record.str());
            }
        }

        if (tournament_count != 1099 || vertex_inequalities != 5405 ||
            literal_capacity_coordinates != 10650 || transitive_count != 153 ||
            equality_count != 153 || injection_objects != injection_images ||
            injection_objects + empty_first_component_objects != 78418)
            throw std::runtime_error("endpoint universe cardinality check failed");
        std::cout << "endpoint_universe=all labelled tournaments orders 1..5\n";
        std::cout << "endpoint_tournaments=" << tournament_count << "\n";
        std::cout << "endpoint_vertex_inequalities=" << vertex_inequalities << "\n";
        std::cout << "endpoint_injection_objects=" << injection_objects << "\n";
        std::cout << "endpoint_distinct_images=" << injection_images << "\n";
        std::cout << "endpoint_empty_first_component_objects=" << empty_first_component_objects << "\n";
        std::cout << "endpoint_two_path_covers_including_empty="
                  << injection_objects + empty_first_component_objects << "\n";
        std::cout << "endpoint_literal_capacity_coordinates=" << literal_capacity_coordinates << "\n";
        std::cout << "endpoint_literal_odd_path_cover_objects=" << literal_capacity_objects << "\n";
        std::cout << "endpoint_transitive_labels=" << transitive_count << "\n";
        std::cout << "endpoint_energy_equalities=" << equality_count << "\n";
        std::cout << "endpoint_aggregate_E=" << aggregate_E << "\n";
        std::cout << "endpoint_record_fnv64=" << std::hex << std::setw(16) << std::setfill('0')
                  << moment_digest << std::dec << std::setfill(' ') << "\n";
        std::cout << "endpoint_injection=PASS\n";
        std::cout << "literal_permutation_vs_subset_capacity=PASS\n";
        std::cout << "endpoint_energy_iff_transitive=PASS\n";

        // -----------------------------------------------------------------
        // 2. Direct arbitrary-context recurrence, especially the n=0 window.
        // -----------------------------------------------------------------
        struct NamedTournament { int n; u64 label; Tournament t; };
        std::vector<NamedTournament> through3;
        for (int n = 1; n <= 3; ++n) {
            const u64 labels = 1ULL << (n * (n - 1) / 2);
            for (u64 label = 0; label < labels; ++label) {
                NamedTournament x{n, label, labelled_tournament(n, label)};
                through3.push_back(x);
            }
        }

        u64 recurrence_contexts = 0, recurrence_windows = 0;
        u64 recurrence_digest = 1469598103934665603ULL;
        auto check_context = [&](const NamedTournament& nb, const NamedTournament& nc, int max_n) {
            std::vector<i128> fs;
            for (int n = 0; n <= max_n; ++n) fs.push_back(F(n, nb.t, nc.t));
            for (int n = 0; n + 4 <= max_n; ++n) {
                const i128 rhs = 9 * fs[n + 3] - 28 * fs[n + 2] + 36 * fs[n + 1] - 16 * fs[n];
                if (fs[n + 4] != rhs) throw std::runtime_error("C-finite recurrence failed");
                ++recurrence_windows;
            }
            ++recurrence_contexts;
            std::ostringstream record;
            record << "b=" << nb.n << ':' << nb.label << ";c=" << nc.n << ':' << nc.label << ";F=";
            for (const i128& x : fs) record << x << ',';
            const std::string line = record.str();
            recurrence_digest = fnv1a_update(recurrence_digest, line);
            std::cout << "recurrence_record=" << line << "\n";
        };

        // Every one of the 11^2 ordered labelled pairs through factor order
        // three gets the load-bearing n=0 window.  The 3^2 pairs through order
        // two additionally get the consecutive n=1 and n=2 windows.
        for (const auto& b : through3)
            for (const auto& c : through3)
                check_context(b, c, (b.n <= 2 && c.n <= 2) ? 6 : 4);

        if (recurrence_contexts != 121 || recurrence_windows != 139)
            throw std::runtime_error("recurrence universe cardinality check failed");
        std::cout << "recurrence_universe=all 121 ordered labelled factor pairs through order 3 at n=0; 9 pairs through order 2 also at n=1,2\n";
        std::cout << "recurrence_contexts=" << recurrence_contexts << "\n";
        std::cout << "recurrence_windows=" << recurrence_windows << "\n";
        std::cout << "recurrence_n0_windows=121\n";
        std::cout << "recurrence_record_fnv64=" << std::hex << std::setw(16) << std::setfill('0')
                  << recurrence_digest << std::dec << std::setfill(' ') << "\n";
        std::cout << "arbitrary_context_recurrence_including_n0=PASS\n";
        const Tournament singleton = labelled_tournament(1, 0);
        if (Rplus(cycle3(), singleton) != -72 || F(4, singleton, singleton) != -180 ||
            F(5, singleton, singleton) != 10764)
            throw std::runtime_error("signed recurrence controls failed");
        std::cout << "signed_controls=R(C3,1):-72,F4(1,1):-180,F5(1,1):10764 PASS\n";

        // -----------------------------------------------------------------
        // 3. Direct lollipop formula: R_+(C3, P_a ▷ C3).
        // -----------------------------------------------------------------
        u64 lollipop_digest = 1469598103934665603ULL;
        for (int a = 0; a <= 8; ++a) {
            const Tournament right = (a == 0) ? cycle3() : ordinal(transitive_tournament(a), cycle3());
            const i128 actual = Rplus(cycle3(), right);
            const i128 expected = (a == 0)
                ? i128(3204)
                : i128(162) * (36 * pow4(a) - (3 * a + 20) * pow2(a) + 4);
            if (actual != expected) throw std::runtime_error("lollipop formula failed");
            std::ostringstream record;
            record << a << ':' << actual;
            lollipop_digest = fnv1a_update(lollipop_digest, record.str());
            std::cout << "lollipop_a=" << a << " actual=" << actual << " expected=" << expected << "\n";
        }
        std::cout << "lollipop_universe=a=0..8 (direct total orders 6..14)\n";
        std::cout << "lollipop_a0_boundary=direct 3204; formal a>=1 expression at a=0 is 3240\n";
        std::cout << "lollipop_record_fnv64=" << std::hex << std::setw(16) << std::setfill('0')
                  << lollipop_digest << std::dec << std::setfill(' ') << "\n";
        std::cout << "lollipop_formula_with_a0_exception=PASS\n";

        std::cout << "analysis_cache_entries=" << analysis_cache.size() << "\n";
        std::cout << "OVERALL=PASS\n";
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
