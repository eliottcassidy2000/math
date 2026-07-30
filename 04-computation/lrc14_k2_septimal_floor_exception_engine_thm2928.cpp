#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using i64 = std::int64_t;

static i64 gcd64(i64 a, i64 b) {
    a = std::llabs(a);
    b = std::llabs(b);
    while (b) {
        i64 r = a % b;
        a = b;
        b = r;
    }
    return a;
}

static i64 lcm64(i64 a, i64 b) {
    return a / gcd64(a, b) * b;
}

struct Rat {
    i64 n = 0;
    i64 d = 1;

    Rat() = default;
    Rat(i64 numerator) : n(numerator), d(1) {}
    Rat(i64 numerator, i64 denominator) : n(numerator), d(denominator) {
        if (d == 0) throw std::runtime_error("zero rational denominator");
        if (d < 0) {
            n = -n;
            d = -d;
        }
        i64 g = gcd64(n, d);
        n /= g;
        d /= g;
    }
};

static Rat operator+(const Rat& a, const Rat& b) {
    return Rat(a.n * b.d + b.n * a.d, a.d * b.d);
}
static Rat operator-(const Rat& a, const Rat& b) {
    return Rat(a.n * b.d - b.n * a.d, a.d * b.d);
}
static Rat operator*(const Rat& a, const Rat& b) {
    return Rat(a.n * b.n, a.d * b.d);
}
static Rat operator/(const Rat& a, const Rat& b) {
    if (b.n == 0) throw std::runtime_error("division by zero rational");
    return Rat(a.n * b.d, a.d * b.n);
}
static bool operator<(const Rat& a, const Rat& b) {
    return a.n * b.d < b.n * a.d;
}
static bool operator>(const Rat& a, const Rat& b) { return b < a; }

static std::vector<int> divisors_of(int number) {
    std::vector<int> low, high;
    for (int d = 1; static_cast<i64>(d) * d <= number; ++d) {
        if (number % d) continue;
        low.push_back(d);
        if (d * d != number) high.push_back(number / d);
    }
    std::reverse(high.begin(), high.end());
    low.insert(low.end(), high.begin(), high.end());
    return low;
}

static int mobius(int number) {
    int result = 1;
    int remaining = number;
    for (int prime = 2; static_cast<i64>(prime) * prime <= remaining; ++prime) {
        if (remaining % prime) continue;
        remaining /= prime;
        if (remaining % prime == 0) return 0;
        result = -result;
        while (remaining % prime == 0) remaining /= prime;
    }
    if (remaining > 1) result = -result;
    return result;
}

static i64 choose(i64 n, int k) {
    if (k < 0 || n < k) return 0;
    k = std::min<i64>(k, n - k);
    i64 result = 1;
    for (int j = 1; j <= k; ++j) {
        result = result * (n - k + j) / j;
    }
    return result;
}

static i64 multichoose(int types, int copies) {
    if (copies == 0) return 1;
    if (types == 0) return 0;
    return choose(static_cast<i64>(types) + copies - 1, copies);
}

struct Constraint {
    int mask;
    Rat rhs;
};

static bool solve_square(
    const std::vector<Constraint>& constraints,
    const std::vector<int>& chosen,
    int n,
    std::vector<Rat>& answer
) {
    std::array<std::array<Rat, 5>, 4> matrix{};
    for (int row = 0; row < n; ++row) {
        const auto& constraint = constraints[chosen[row]];
        for (int column = 0; column < n; ++column) {
            matrix[row][column] = Rat((constraint.mask >> column) & 1);
        }
        matrix[row][n] = constraint.rhs;
    }
    for (int column = 0; column < n; ++column) {
        int pivot = -1;
        for (int row = column; row < n; ++row) {
            if (matrix[row][column].n != 0) {
                pivot = row;
                break;
            }
        }
        if (pivot < 0) return false;
        std::swap(matrix[column], matrix[pivot]);
        Rat scale = matrix[column][column];
        for (int j = column; j <= n; ++j) matrix[column][j] = matrix[column][j] / scale;
        for (int row = 0; row < n; ++row) {
            if (row == column || matrix[row][column].n == 0) continue;
            scale = matrix[row][column];
            for (int j = column; j <= n; ++j) {
                matrix[row][j] = matrix[row][j] - scale * matrix[column][j];
            }
        }
    }
    answer.resize(n);
    for (int i = 0; i < n; ++i) answer[i] = matrix[i][n];
    return true;
}

struct EventVertices {
    std::vector<int> minimal;
    std::vector<std::vector<Rat>> vertices;
};

static std::unordered_map<std::uint32_t, EventVertices> vertex_cache;
static std::unordered_map<std::uint64_t, bool> pass_cache;

static EventVertices build_event_vertices(int n, int truth) {
    EventVertices result;
    const int states = 1 << n;
    if (truth & 1) throw std::runtime_error("proper threshold event contains empty state");
    for (int state = 1; state < states; ++state) {
        if (!((truth >> state) & 1)) continue;
        bool minimal = true;
        for (int i = 0; i < n; ++i) {
            if ((state >> i) & 1) {
                if ((truth >> (state ^ (1 << i))) & 1) minimal = false;
            }
        }
        if (minimal) result.minimal.push_back(state);
    }
    if (result.minimal.empty()) throw std::runtime_error("nonempty event has no minimal state");

    std::vector<Constraint> constraints;
    for (int mask : result.minimal) constraints.push_back({mask, Rat(1)});
    for (int i = 0; i < n; ++i) constraints.push_back({1 << i, Rat(0)});

    std::vector<int> chosen;
    auto visit = [&](auto&& self, int start) -> void {
        if (static_cast<int>(chosen.size()) == n) {
            std::vector<Rat> point;
            if (!solve_square(constraints, chosen, n, point)) return;
            for (const Rat& value : point) if (value < Rat(0)) return;
            for (int mask : result.minimal) {
                Rat total(0);
                for (int i = 0; i < n; ++i) if ((mask >> i) & 1) total = total + point[i];
                if (total < Rat(1)) return;
            }
            result.vertices.push_back(point);
            return;
        }
        int needed = n - static_cast<int>(chosen.size());
        for (int index = start; index + needed <= static_cast<int>(constraints.size()); ++index) {
            chosen.push_back(index);
            self(self, index + 1);
            chosen.pop_back();
        }
    };
    visit(visit, 0);
    if (result.vertices.empty()) throw std::runtime_error("fractional cover has no vertex");
    return result;
}

static bool event_passes(int n, const std::vector<int>& residues, int truth) {
    int residue_code = 0;
    for (int r : residues) residue_code = residue_code * 7 + r;
    std::uint64_t pass_key =
        (static_cast<std::uint64_t>(n) << 48)
        | (static_cast<std::uint64_t>(residue_code) << 16)
        | static_cast<std::uint16_t>(truth);
    auto known = pass_cache.find(pass_key);
    if (known != pass_cache.end()) return known->second;

    std::uint32_t vertex_key = (static_cast<std::uint32_t>(n) << 16) | truth;
    auto vertex_it = vertex_cache.find(vertex_key);
    if (vertex_it == vertex_cache.end()) {
        vertex_it = vertex_cache.emplace(vertex_key, build_event_vertices(n, truth)).first;
    }
    bool have = false;
    Rat optimum;
    for (const auto& point : vertex_it->second.vertices) {
        Rat cost(0);
        for (int i = 0; i < n; ++i) cost = cost + Rat(residues[i]) * point[i];
        if (!have || cost < optimum) {
            have = true;
            optimum = cost;
        }
    }
    if (!have) throw std::runtime_error("event optimum missing");
    // Actual cost is optimum/7.  Compact/open strictness requires >66/91.
    bool passes = optimum > Rat(462, 91);
    pass_cache.emplace(pass_key, passes);
    return passes;
}

static int status_allowance(
    const std::vector<int>& weights,
    const std::vector<int>& residues
) {
    int n = static_cast<int>(weights.size());
    if (n == 0) return 0;
    std::vector<int> sums(1 << n, 0);
    for (int state = 1; state < (1 << n); ++state) {
        int bit = __builtin_ctz(state);
        sums[state] = sums[state ^ (1 << bit)] + weights[bit];
    }
    std::vector<int> thresholds = sums;
    std::sort(thresholds.begin(), thresholds.end());
    thresholds.erase(std::unique(thresholds.begin(), thresholds.end()), thresholds.end());
    if (thresholds.front() != 0) throw std::runtime_error("zero subset sum missing");
    int low = 0;
    int high = static_cast<int>(thresholds.size());
    while (high - low > 1) {
        int middle = (low + high) / 2;
        int need = thresholds[middle];
        int truth = 0;
        for (int state = 0; state < (1 << n); ++state) {
            if (sums[state] >= need) truth |= 1 << state;
        }
        if (event_passes(n, residues, truth)) low = middle;
        else high = middle;
    }
    return thresholds[low];
}

struct CoefficientKey {
    int c;
    int z;
    int l;
    bool operator<(const CoefficientKey& other) const {
        return std::tie(c, z, l) < std::tie(other.c, other.z, other.l);
    }
};

struct DistributionBuilder {
    int D;
    int q;
    std::vector<int> divisors_D;
    std::vector<int> positive_types;
    std::map<CoefficientKey, i64> coefficient_cache;
    std::array<std::map<int, i64>, 6> distribution;

    explicit DistributionBuilder(int modulus)
        : D(modulus), q(modulus / 7), divisors_D(divisors_of(modulus)) {
        if (D % 7) throw std::runtime_error("nonseptimal resolving modulus");
        for (int d : divisors_of(q)) if (d > 1 && d % 7) positive_types.push_back(d);
    }

    int uniform_types(int E) const {
        int count = 0;
        for (int d : divisors_of(E)) if (d > 1 && q % d) ++count;
        return count;
    }

    int zero_residue_types(int E) const {
        int count = 0;
        int common = gcd64(E, q);
        for (int d : divisors_of(common)) if (d > 1 && d % 7 == 0) ++count;
        return count;
    }

    int spike_types(int E) const {
        return static_cast<int>(divisors_of(gcd64(E, q)).size()) - 1;
    }

    i64 completion_coefficient(int c, int z, int selected_lcm) {
        CoefficientKey key{c, z, selected_lcm};
        auto found = coefficient_cache.find(key);
        if (found != coefficient_cache.end()) return found->second;
        i64 result = 0;
        for (int E : divisors_D) {
            if (E % selected_lcm) continue;
            int sign = mobius(D / E);
            if (!sign) continue;
            result += sign
                * multichoose(uniform_types(E), c)
                * multichoose(zero_residue_types(E), z);
        }
        if (result < 0) throw std::runtime_error("negative exact-lcm completion coefficient");
        coefficient_cache.emplace(key, result);
        return result;
    }

    void accept_positive_multiset(
        int positive_count,
        int selected_lcm,
        int deterministic_base,
        const std::vector<int>& weights,
        const std::vector<int>& residues
    ) {
        int optional = status_allowance(weights, residues);
        for (int c = 1; c <= 5 - positive_count; ++c) {
            int z = 5 - c - positive_count;
            i64 coefficient = completion_coefficient(c, z, selected_lcm);
            if (!coefficient) continue;
            if (z && q % 7) throw std::runtime_error("zero-residue type off septimal q");
            int allowance = deterministic_base + optional + z * (q / 7);
            distribution[c][allowance] += coefficient;
        }
    }

    void enumerate_size(int wanted) {
        std::vector<int> weights, residues;
        auto visit = [&](auto&& self, int start, int left, int selected_lcm, int base) -> void {
            if (left == 0) {
                accept_positive_multiset(wanted, selected_lcm, base, weights, residues);
                return;
            }
            for (int index = start; index < static_cast<int>(positive_types.size()); ++index) {
                int d = positive_types[index];
                int r = d % 7;
                int w = q / d;
                if (q % d) throw std::runtime_error("positive spike type does not divide q");
                weights.push_back(w);
                residues.push_back(r);
                self(self, index, left - 1, lcm64(selected_lcm, d), base + (d / 7) * w);
                residues.pop_back();
                weights.pop_back();
            }
        };
        visit(visit, 0, wanted, 1, 0);
    }

    i64 expected_c_count(int c) const {
        i64 result = 0;
        for (int E : divisors_D) {
            int sign = mobius(D / E);
            if (!sign) continue;
            result += sign
                * multichoose(uniform_types(E), c)
                * multichoose(spike_types(E), 5 - c);
        }
        return result;
    }

    void build() {
        for (int positive_count = 0; positive_count <= 4; ++positive_count) {
            enumerate_size(positive_count);
        }
        i64 total = 0;
        for (int c = 1; c <= 5; ++c) {
            i64 count = 0;
            for (const auto& item : distribution[c]) {
                if (item.second < 0) throw std::runtime_error("negative allowance coefficient");
                count += item.second;
            }
            i64 expected = expected_c_count(c);
            if (count != expected) throw std::runtime_error("allowance GF does not recover c coefficient");
            total += count;
        }
        i64 raw = 0;
        for (int E : divisors_D) {
            int sign = mobius(D / E);
            if (!sign) continue;
            int alphabet = static_cast<int>(divisors_of(E).size()) - 1;
            raw += sign * multichoose(alphabet, 5);
        }
        if (total != raw) throw std::runtime_error("allowance GF does not recover raw shape count");
    }
};

int main() {
    try {
        std::ios::sync_with_stdio(false);
        std::cin.tie(nullptr);
        int case_count;
        if (!(std::cin >> case_count)) throw std::runtime_error("missing case count");
        std::cout << "ENGINE all-c exact floor/exception GF\n";
        for (int case_index = 0; case_index < case_count; ++case_index) {
            int D, query_count;
            std::cin >> D >> query_count;
            std::vector<std::pair<int, int>> queries(query_count);
            for (auto& query : queries) std::cin >> query.first >> query.second;

            DistributionBuilder builder(D);
            builder.build();
            for (int c = 1; c <= 5; ++c) {
                const auto& data = builder.distribution[c];
                i64 raw = 0;
                for (const auto& item : data) raw += item.second;
                int minimum = data.empty() ? -1 : data.begin()->first;
                int maximum = data.empty() ? -1 : data.rbegin()->first;
                std::cout << "C " << D << ' ' << c << ' ' << raw << ' '
                          << data.size() << ' ' << minimum << ' ' << maximum << '\n';
            }

            std::array<std::vector<int>, 6> keys;
            std::array<std::vector<i64>, 6> suffix;
            for (int c = 1; c <= 5; ++c) {
                for (const auto& item : builder.distribution[c]) keys[c].push_back(item.first);
                suffix[c].assign(keys[c].size() + 1, 0);
                for (int index = static_cast<int>(keys[c].size()) - 1; index >= 0; --index) {
                    suffix[c][index] = suffix[c][index + 1]
                        + builder.distribution[c].at(keys[c][index]);
                }
            }
            for (const auto& query : queries) {
                int c = query.first;
                int need = query.second;
                auto it = std::lower_bound(keys[c].begin(), keys[c].end(), need);
                int index = static_cast<int>(it - keys[c].begin());
                i64 count = suffix[c][index];
                int first = it == keys[c].end() ? -1 : *it;
                std::cout << "Q " << D << ' ' << c << ' ' << need << ' '
                          << count << ' ' << first << '\n';
            }
        }
        std::cout << "event_vertex_cache=" << vertex_cache.size() << '\n';
        std::cout << "event_pass_cache=" << pass_cache.size() << '\n';
        std::cout << "all_exact_engine_controls=PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
