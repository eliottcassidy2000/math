// Independent literal-permutation audit for THM-4187.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using i64 = long long;
using u64 = std::uint64_t;
using Out = std::vector<std::uint32_t>;
using Vec2 = std::array<i64, 2>;

static void need(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

static std::string label(const Out& out) {
    std::string answer;
    const int n = static_cast<int>(out.size());
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            answer.push_back((out[static_cast<std::size_t>(i)] &
                              (std::uint32_t{1} << j)) ? '1' : '0');
        }
    }
    return answer;
}

static Out parse(const std::string& bits, int n) {
    need(static_cast<int>(bits.size()) == n * (n - 1) / 2, "bad label length");
    Out out(static_cast<std::size_t>(n), 0);
    int cursor = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (bits[static_cast<std::size_t>(cursor)] == '1') {
                out[static_cast<std::size_t>(i)] |= std::uint32_t{1} << j;
            } else {
                out[static_cast<std::size_t>(j)] |= std::uint32_t{1} << i;
            }
            ++cursor;
        }
    }
    return out;
}

static Out transitive(int n) {
    Out out(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            out[static_cast<std::size_t>(i)] |= std::uint32_t{1} << j;
        }
    }
    return out;
}

static std::vector<Out> all_labelled(int n) {
    const int pairs = n * (n - 1) / 2;
    const u64 total = u64{1} << pairs;
    std::vector<Out> answer;
    answer.reserve(static_cast<std::size_t>(total));
    for (u64 mask = 0; mask < total; ++mask) {
        std::string bits(static_cast<std::size_t>(pairs), '0');
        for (int bit = 0; bit < pairs; ++bit) {
            if (mask & (u64{1} << bit)) bits[static_cast<std::size_t>(bit)] = '1';
        }
        answer.push_back(parse(bits, n));
    }
    return answer;
}

static bool path_valid(const Out& out, const std::vector<int>& path) {
    for (std::size_t i = 0; i + 1 < path.size(); ++i) {
        if (!(out[static_cast<std::size_t>(path[i])] &
              (std::uint32_t{1} << path[i + 1]))) return false;
    }
    return true;
}

static bool has_sink(const Out& out) {
    return std::any_of(out.begin(), out.end(),
                       [](std::uint32_t row) { return row == 0; });
}

static Vec2 convolution(const Vec2& left, const Vec2& right) {
    return {
        left[0] * right[0] + left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    };
}

struct Data {
    Out out;
    i64 hamilton{};
    std::vector<std::vector<i64>> capacities;
    std::vector<Vec2> starts;
    std::vector<Vec2> ends;
    i64 mass{};
    i64 gate{};

    Vec2 optional_state() const {
        Vec2 answer{hamilton, 0};
        for (const Vec2& row : starts) {
            answer[0] += row[0];
            answer[1] += row[1];
        }
        return answer;
    }
};

static i64 gate(const Out& out, const std::vector<std::vector<i64>>& capacities) {
    const int n = static_cast<int>(out.size());
    std::vector<i64> degrees(static_cast<std::size_t>(n), 0);
    std::vector<i64> currents(static_cast<std::size_t>(n), 0);
    i64 mass = 0;
    i64 squares = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const i64 value = capacities[static_cast<std::size_t>(i)]
                                        [static_cast<std::size_t>(j)];
            mass += value;
            squares += value * value;
            degrees[static_cast<std::size_t>(i)] += value;
            degrees[static_cast<std::size_t>(j)] += value;
            if (out[static_cast<std::size_t>(i)] & (std::uint32_t{1} << j)) {
                currents[static_cast<std::size_t>(i)] += value;
                currents[static_cast<std::size_t>(j)] -= value;
            } else {
                currents[static_cast<std::size_t>(i)] -= value;
                currents[static_cast<std::size_t>(j)] += value;
            }
        }
    }
    i64 degree_squares = 0;
    i64 signed_current = 0;
    for (int i = 0; i < n; ++i) {
        degree_squares += degrees[static_cast<std::size_t>(i)] *
                          degrees[static_cast<std::size_t>(i)];
        signed_current += degrees[static_cast<std::size_t>(i)] *
                          currents[static_cast<std::size_t>(i)];
    }
    const i64 numerator = mass * mass + squares - degree_squares;
    need(numerator % 2 == 0, "disjoint-pair parity");
    return numerator / 2 + 2 * signed_current;
}

static Data literal_data(const Out& out) {
    const int n = static_cast<int>(out.size());
    const int size = 1 << n;
    const int full = size - 1;
    std::vector<i64> hamilton(static_cast<std::size_t>(size), 0);
    std::vector<std::vector<i64>> starts_by_subset(
        static_cast<std::size_t>(size),
        std::vector<i64>(static_cast<std::size_t>(n), 0));
    std::vector<std::vector<i64>> ends_by_subset(
        static_cast<std::size_t>(size),
        std::vector<i64>(static_cast<std::size_t>(n), 0));
    hamilton[0] = 1;
    for (int mask = 1; mask < size; ++mask) {
        std::vector<int> vertices;
        for (int vertex = 0; vertex < n; ++vertex) {
            if (mask & (1 << vertex)) vertices.push_back(vertex);
        }
        do {
            if (path_valid(out, vertices)) {
                ++hamilton[static_cast<std::size_t>(mask)];
                ++starts_by_subset[static_cast<std::size_t>(mask)]
                                  [static_cast<std::size_t>(vertices.front())];
                ++ends_by_subset[static_cast<std::size_t>(mask)]
                                [static_cast<std::size_t>(vertices.back())];
            }
        } while (std::next_permutation(vertices.begin(), vertices.end()));
    }

    std::vector<Vec2> starts(static_cast<std::size_t>(n), Vec2{0, 0});
    std::vector<Vec2> ends(static_cast<std::size_t>(n), Vec2{0, 0});
    for (int mask = 1; mask < size; ++mask) {
        const int parity = __builtin_popcount(static_cast<unsigned>(mask)) & 1;
        const i64 complement = hamilton[static_cast<std::size_t>(full ^ mask)];
        for (int vertex = 0; vertex < n; ++vertex) {
            starts[static_cast<std::size_t>(vertex)][static_cast<std::size_t>(parity)] +=
                starts_by_subset[static_cast<std::size_t>(mask)]
                                [static_cast<std::size_t>(vertex)] * complement;
            ends[static_cast<std::size_t>(vertex)][static_cast<std::size_t>(parity)] +=
                ends_by_subset[static_cast<std::size_t>(mask)]
                              [static_cast<std::size_t>(vertex)] * complement;
        }
    }

    std::vector<std::vector<i64>> exposed(
        static_cast<std::size_t>(n),
        std::vector<i64>(static_cast<std::size_t>(n), 0));
    std::vector<int> word(static_cast<std::size_t>(n));
    std::iota(word.begin(), word.end(), 0);
    do {
        for (int cut = 1; cut < n; ++cut) {
            const std::vector<int> left(word.begin(), word.begin() + cut);
            const std::vector<int> right(word.begin() + cut, word.end());
            if (path_valid(out, left) && path_valid(out, right)) {
                ++exposed[static_cast<std::size_t>(left.back())]
                         [static_cast<std::size_t>(right.front())];
            }
        }
    } while (std::next_permutation(word.begin(), word.end()));

    std::vector<std::vector<i64>> capacities(
        static_cast<std::size_t>(n),
        std::vector<i64>(static_cast<std::size_t>(n), 0));
    i64 mass = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const i64 value = exposed[static_cast<std::size_t>(i)]
                                    [static_cast<std::size_t>(j)] +
                              exposed[static_cast<std::size_t>(j)]
                                    [static_cast<std::size_t>(i)];
            capacities[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = value;
            capacities[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] = value;
            mass += value;
        }
    }
    return Data{out, hamilton[static_cast<std::size_t>(full)], capacities,
                starts, ends, mass, gate(out, capacities)};
}

static Out ordinal_out(const Out& left, const Out& right) {
    const int nl = static_cast<int>(left.size());
    const int nr = static_cast<int>(right.size());
    Out answer(static_cast<std::size_t>(nl + nr), 0);
    const std::uint32_t right_mask =
        ((std::uint32_t{1} << nr) - 1) << nl;
    for (int i = 0; i < nl; ++i) {
        answer[static_cast<std::size_t>(i)] =
            left[static_cast<std::size_t>(i)] | right_mask;
    }
    for (int j = 0; j < nr; ++j) {
        answer[static_cast<std::size_t>(nl + j)] =
            right[static_cast<std::size_t>(j)] << nl;
    }
    return answer;
}

static Data compose(const Data& left, const Data& right) {
    const int nl = static_cast<int>(left.out.size());
    const int nr = static_cast<int>(right.out.size());
    const int n = nl + nr;
    std::vector<std::vector<i64>> capacities(
        static_cast<std::size_t>(n),
        std::vector<i64>(static_cast<std::size_t>(n), 0));
    for (int i = 0; i < nl; ++i) {
        for (int j = i + 1; j < nl; ++j) {
            const i64 value = right.hamilton *
                left.capacities[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
            capacities[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = value;
            capacities[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] = value;
        }
    }
    for (int i = 0; i < nr; ++i) {
        for (int j = i + 1; j < nr; ++j) {
            const i64 value = left.hamilton *
                right.capacities[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
            capacities[static_cast<std::size_t>(nl + i)]
                      [static_cast<std::size_t>(nl + j)] = value;
            capacities[static_cast<std::size_t>(nl + j)]
                      [static_cast<std::size_t>(nl + i)] = value;
        }
    }
    for (int i = 0; i < nl; ++i) {
        for (int j = 0; j < nr; ++j) {
            const i64 value = 2 * (
                left.starts[static_cast<std::size_t>(i)][0] *
                    right.ends[static_cast<std::size_t>(j)][0] +
                left.starts[static_cast<std::size_t>(i)][1] *
                    right.ends[static_cast<std::size_t>(j)][1]);
            capacities[static_cast<std::size_t>(i)]
                      [static_cast<std::size_t>(nl + j)] = value;
            capacities[static_cast<std::size_t>(nl + j)]
                      [static_cast<std::size_t>(i)] = value;
        }
    }
    const Vec2 left_optional = left.optional_state();
    const Vec2 right_optional = right.optional_state();
    std::vector<Vec2> starts;
    std::vector<Vec2> ends;
    starts.reserve(static_cast<std::size_t>(n));
    ends.reserve(static_cast<std::size_t>(n));
    for (const Vec2& row : left.starts) starts.push_back(convolution(row, right_optional));
    for (const Vec2& row : right.starts) {
        starts.push_back({left.hamilton * row[0], left.hamilton * row[1]});
    }
    for (const Vec2& row : left.ends) {
        ends.push_back({right.hamilton * row[0], right.hamilton * row[1]});
    }
    for (const Vec2& row : right.ends) ends.push_back(convolution(left_optional, row));
    i64 mass = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            mass += capacities[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        }
    }
    const Out out = ordinal_out(left.out, right.out);
    return Data{out, left.hamilton * right.hamilton, capacities, starts, ends,
                mass, gate(out, capacities)};
}

static bool data_equal(const Data& left, const Data& right) {
    return left.out == right.out && left.hamilton == right.hamilton &&
           left.capacities == right.capacities && left.starts == right.starts &&
           left.ends == right.ends && left.mass == right.mass && left.gate == right.gate;
}

static i64 remainder(const Data& left, const Data& right) {
    const Data child = compose(left, right);
    return child.gate - right.hamilton * right.hamilton * left.gate -
           left.hamilton * left.hamilton * right.gate;
}

static i64 theta(const Data& left, const Data& middle, const Data& right) {
    return remainder(compose(left, middle), right) -
           left.hamilton * left.hamilton * remainder(middle, right);
}

static Out converse_out(const Out& out) {
    const int n = static_cast<int>(out.size());
    Out answer(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (out[static_cast<std::size_t>(i)] & (std::uint32_t{1} << j)) {
                answer[static_cast<std::size_t>(j)] |= std::uint32_t{1} << i;
            } else {
                answer[static_cast<std::size_t>(i)] |= std::uint32_t{1} << j;
            }
        }
    }
    return answer;
}

static i64 minus_gate(const Data& data) {
    return gate(converse_out(data.out), data.capacities);
}

static i64 minus_remainder(const Data& left, const Data& right) {
    const Data child = compose(left, right);
    return minus_gate(child) - right.hamilton * right.hamilton * minus_gate(left) -
           left.hamilton * left.hamilton * minus_gate(right);
}

static i64 incoming_mass(const Data& data, int vertex) {
    i64 answer = 0;
    for (int other = 0; other < static_cast<int>(data.out.size()); ++other) {
        if (data.out[static_cast<std::size_t>(other)] &
            (std::uint32_t{1} << vertex)) {
            answer += data.capacities[static_cast<std::size_t>(vertex)]
                                     [static_cast<std::size_t>(other)];
        }
    }
    return answer;
}

static i64 delta(const Data& one, const Data& data) {
    return compose(one, data).gate - data.gate;
}

static i64 padding_formula(const Data& one, const Data& data) {
    const Data child = compose(one, data);
    const int n = static_cast<int>(data.out.size());
    std::vector<i64> root(static_cast<std::size_t>(n), 0);
    std::vector<i64> degree(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < n; ++i) {
        root[static_cast<std::size_t>(i)] =
            child.capacities[0][static_cast<std::size_t>(1 + i)];
        degree[static_cast<std::size_t>(i)] =
            std::accumulate(data.capacities[static_cast<std::size_t>(i)].begin(),
                            data.capacities[static_cast<std::size_t>(i)].end(), i64{0});
    }
    i64 answer = 0;
    for (int vertex = 0; vertex < n; ++vertex) {
        i64 incoming_root = 0;
        for (int other = 0; other < n; ++other) {
            if (data.out[static_cast<std::size_t>(other)] &
                (std::uint32_t{1} << vertex)) {
                incoming_root += root[static_cast<std::size_t>(other)];
            }
        }
        const i64 fan_defect = incoming_root - incoming_mass(data, vertex);
        need(fan_defect >= 0, "negative literal fan defect");
        answer += root[static_cast<std::size_t>(vertex)] *
                  (data.mass - degree[static_cast<std::size_t>(vertex)]);
        answer += 4 * root[static_cast<std::size_t>(vertex)] * fan_defect;
    }
    need(answer == delta(one, data), "literal padding formula");
    return answer;
}

static std::tuple<i64, i64, i64> supermodularity_split(
    const Data& one, const Data& left, const Data& right) {
    const Data product = compose(left, right);
    const Data product_child = compose(one, product);
    const Data left_child = compose(one, left);
    const int nl = static_cast<int>(left.out.size());
    const int nr = static_cast<int>(right.out.size());
    std::vector<i64> old_root(static_cast<std::size_t>(nl), 0);
    std::vector<i64> all_root(static_cast<std::size_t>(nl + nr), 0);
    for (int i = 0; i < nl; ++i) {
        old_root[static_cast<std::size_t>(i)] =
            left_child.capacities[0][static_cast<std::size_t>(1 + i)];
    }
    for (int i = 0; i < nl + nr; ++i) {
        all_root[static_cast<std::size_t>(i)] =
            product_child.capacities[0][static_cast<std::size_t>(1 + i)];
    }
    for (int i = 0; i < nl; ++i) {
        need(all_root[static_cast<std::size_t>(i)] ==
             right.hamilton * old_root[static_cast<std::size_t>(i)],
             "independent root scaling");
    }
    std::vector<i64> cross_rows(static_cast<std::size_t>(nl), 0);
    for (int i = 0; i < nl; ++i) {
        for (int j = 0; j < nr; ++j) {
            cross_rows[static_cast<std::size_t>(i)] +=
                product.capacities[static_cast<std::size_t>(i)]
                                  [static_cast<std::size_t>(nl + j)];
        }
    }
    const i64 cross_mass =
        std::accumulate(cross_rows.begin(), cross_rows.end(), i64{0});
    i64 left_extra = 0;
    for (int i = 0; i < nl; ++i) {
        left_extra += right.hamilton * old_root[static_cast<std::size_t>(i)] *
            (left.hamilton * right.mass + cross_mass -
             cross_rows[static_cast<std::size_t>(i)]);
    }
    i64 right_terms = 0;
    for (int j = 0; j < nr; ++j) {
        const int vertex = nl + j;
        const i64 root_capacity = all_root[static_cast<std::size_t>(vertex)];
        const i64 degree = std::accumulate(
            product.capacities[static_cast<std::size_t>(vertex)].begin(),
            product.capacities[static_cast<std::size_t>(vertex)].end(), i64{0});
        i64 incoming_root = 0;
        for (int other = 0; other < nl + nr; ++other) {
            if (product.out[static_cast<std::size_t>(other)] &
                (std::uint32_t{1} << vertex)) {
                incoming_root += all_root[static_cast<std::size_t>(other)];
            }
        }
        const i64 fan_defect = incoming_root - incoming_mass(product, vertex);
        need(fan_defect >= 0, "negative independent right fan defect");
        right_terms += root_capacity * (product.mass - degree) +
                       4 * root_capacity * fan_defect;
    }
    const i64 difference = delta(one, product) -
                           right.hamilton * right.hamilton * delta(one, left);
    need(difference == left_extra + right_terms,
         "independent supermodularity identity");
    return {difference, left_extra, right_terms};
}

static i64 negative_remainder_formula(int n) {
    if (n == 1) return -72;
    return 72 - (27 * n + 126) * (i64{1} << (n - 1));
}

static i64 negative_theta_formula(int b, int c) {
    if (b == 1) return 144 - (27 * c + 153) * (i64{1} << c);
    return -(i64{1} << (b - 1)) *
        (((27 * (b + c) + 126) * (i64{1} << c)) - (27 * b + 126));
}

int main() {
    try {
        std::map<int, std::vector<Data>> banks;
        i64 labelled_tournaments = 0;
        for (int n = 1; n <= 4; ++n) {
            for (const Out& out : all_labelled(n)) {
                banks[n].push_back(literal_data(out));
                ++labelled_tournaments;
            }
        }
        need(labelled_tournaments == 75, "labelled factor count through four");
        const Data one = literal_data(transitive(1));
        const Data cycle = literal_data(parse("101", 3));
        need(cycle.hamilton == 3 && cycle.gate == 0, "literal C3 control");

        std::vector<const Data*> factors;
        for (int n = 1; n <= 4; ++n) {
            for (const Data& data : banks[n]) factors.push_back(&data);
        }
        need(factors.size() == 75, "factor pointer total");

        // Directly rebuild small ordinal children by literal permutations.
        std::map<std::string, Data> direct_cache;
        auto direct_cached = [&direct_cache](const Out& out) -> const Data& {
            const std::string key = std::to_string(out.size()) + ":" + label(out);
            auto found = direct_cache.find(key);
            if (found == direct_cache.end()) {
                found = direct_cache.emplace(key, literal_data(out)).first;
            }
            return found->second;
        };
        i64 direct_transfer_checks = 0;
        for (const Data* left : factors) {
            for (const Data* right : factors) {
                if (left->out.size() + right->out.size() > 6) continue;
                const Data child = compose(*left, *right);
                need(data_equal(child, direct_cached(child.out)),
                     "literal ordinal transfer mismatch");
                ++direct_transfer_checks;
            }
        }

        i64 padding_checks = 0;
        for (const Data* data : factors) {
            padding_formula(one, *data);
            ++padding_checks;
        }

        i64 pair_checks = 0;
        i64 pair_zero = 0;
        i64 minimum_positive = std::numeric_limits<i64>::max();
        for (const Data* left : factors) {
            for (const Data* right : factors) {
                const auto [difference, left_extra, right_terms] =
                    supermodularity_split(one, *left, *right);
                need(difference == theta(one, *left, *right),
                     "independent Theta/delta identity");
                need(left_extra >= 0 && right_terms >= 0 && difference >= 0,
                     "independent source supermodularity sign");
                const bool expected_zero = left->out.size() == 1 && right->out.size() == 1;
                need((difference == 0) == expected_zero,
                     "independent source equality boundary");
                pair_zero += difference == 0;
                if (difference > 0) minimum_positive = std::min(minimum_positive, difference);
                ++pair_checks;
            }
        }
        need(pair_checks == 75 * 75 && pair_zero == 1 && minimum_positive == 16,
             "independent pair totals/minimum");

        std::vector<Data> transitive_factors;
        Data current = one;
        for (int n = 1; n <= 8; ++n) {
            if (n > 1) current = compose(one, current);
            need(current.out == transitive(n), "independent transitive recursion");
            transitive_factors.push_back(current);
        }
        i64 remainder_checks = 0;
        i64 remainder_zero = 0;
        i64 interaction_checks = 0;
        i64 interaction_zero = 0;
        i64 os_checks = 0;
        for (const Data& left : transitive_factors) {
            for (const Data* right : factors) {
                const i64 value = remainder(left, *right);
                need(value >= 0, "independent transitive-left remainder");
                remainder_zero += value == 0;
                if (!has_sink(right->out)) {
                    need(value > 0, "independent transitive-left OS+");
                    ++os_checks;
                }
                ++remainder_checks;
            }
            for (const Data* middle : factors) {
                for (const Data* right : factors) {
                    const i64 value = theta(left, *middle, *right);
                    need(value >= 0, "independent transitive-left interaction");
                    interaction_zero += value == 0;
                    ++interaction_checks;
                }
            }
        }
        // The two labelled orientations of P_2 give two copies of the
        // unlabelled (P_1,P_2) equality row.
        need(remainder_zero == 4, "independent remainder zero total");
        need(interaction_zero == 1, "independent interaction zero total");

        i64 dual_remainder_checks = 0;
        i64 dual_remainder_zero = 0;
        for (const Data& right : transitive_factors) {
            for (const Data* left : factors) {
                const Data left_converse = literal_data(converse_out(left->out));
                const i64 value = minus_remainder(*left, right);
                need(value == remainder(right, left_converse),
                     "independent converse remainder identity");
                need(value >= 0, "independent transitive-right dual remainder");
                dual_remainder_zero += value == 0;
                ++dual_remainder_checks;
            }
        }
        need(dual_remainder_zero == 4, "independent dual remainder zero total");

        i64 c3_identity_checks = 0;
        i64 c3_os_checks = 0;
        i64 c3_os_equal = 0;
        for (const Data* data : factors) {
            const Data source_child = compose(one, *data);
            i64 aggregate = data->mass + 2 * data->hamilton;
            i64 sum_a = 0;
            i64 sum_t = 0;
            i64 dot = 0;
            i64 t_squares = 0;
            i64 t_linear = 0;
            for (int i = 0; i < static_cast<int>(data->out.size()); ++i) {
                const i64 a = source_child.capacities[0][static_cast<std::size_t>(1 + i)];
                const i64 t = 2 * data->starts[static_cast<std::size_t>(i)][0] +
                              4 * data->starts[static_cast<std::size_t>(i)][1];
                const i64 degree = std::accumulate(
                    data->capacities[static_cast<std::size_t>(i)].begin(),
                    data->capacities[static_cast<std::size_t>(i)].end(), i64{0});
                i64 outgoing = 0;
                for (int j = 0; j < static_cast<int>(data->out.size()); ++j) {
                    if (data->out[static_cast<std::size_t>(i)] &
                        (std::uint32_t{1} << j)) {
                        outgoing += data->capacities[static_cast<std::size_t>(i)]
                                                    [static_cast<std::size_t>(j)];
                    }
                }
                sum_a += a;
                sum_t += t;
                dot += a * t;
                t_squares += t * t;
                t_linear += t * (data->mass - degree + 4 * outgoing);
            }
            const i64 c3_remainder_formula =
                15 * t_squares + 9 * t_linear - 27 * data->mass * data->mass -
                108 * data->hamilton * data->mass -
                120 * data->hamilton * data->hamilton;
            need(c3_remainder_formula == remainder(*data, cycle),
                 "independent C3 block remainder formula");
            need(sum_a == aggregate && sum_t == 3 * aggregate - 2 * data->hamilton,
                 "independent C3 aggregate identity");
            const i64 exact = 9 *
                (8 * aggregate * (3 * aggregate - data->hamilton) - dot);
            need(exact == theta(one, *data, cycle),
                 "independent C3 interaction identity");
            need(exact >= 648 * data->hamilton * data->hamilton,
                 "independent C3 bound");
            ++c3_identity_checks;
            for (const Data& left : transitive_factors) {
                const i64 bound = 648 * static_cast<i64>(left.out.size()) *
                                  data->hamilton * data->hamilton;
                const i64 value = remainder(left, compose(*data, cycle));
                need(value >= bound, "independent terminal-C3 OS+ bound");
                c3_os_equal += value == bound;
                ++c3_os_checks;
            }
        }
        need(c3_os_equal == 1, "independent terminal-C3 equality total");

        std::vector<Data> long_transitive;
        current = one;
        i64 negative_remainder_checks = 0;
        for (int n = 1; n <= 20; ++n) {
            if (n > 1) current = compose(current, one);
            const i64 value = remainder(cycle, current);
            need(value == negative_remainder_formula(n) && value < 0,
                 "independent negative remainder formula");
            long_transitive.push_back(current);
            ++negative_remainder_checks;
        }
        i64 negative_interaction_checks = 0;
        for (int b = 1; b <= 8; ++b) {
            for (int c = 1; c <= 8; ++c) {
                const i64 value = theta(
                    cycle,
                    long_transitive[static_cast<std::size_t>(b - 1)],
                    long_transitive[static_cast<std::size_t>(c - 1)]);
                need(value == negative_theta_formula(b, c) && value < 0,
                     "independent negative interaction formula");
                ++negative_interaction_checks;
            }
        }
        need(theta(cycle, one, one) == -216, "independent minimal hostile");

        std::cout << "TOURNAMENT_SOURCE_PADDING_SUPERMODULARITY_INDEPENDENT\n";
        std::cout << "literal_labelled_factors " << labelled_tournaments << '\n';
        std::cout << "literal_direct_transfer_checks " << direct_transfer_checks << '\n';
        std::cout << "literal_padding_formula_checks " << padding_checks << '\n';
        std::cout << "labelled_pair_supermodularity_checks " << pair_checks << '\n';
        std::cout << "labelled_pair_zero " << pair_zero << '\n';
        std::cout << "labelled_pair_minimum_positive " << minimum_positive << '\n';
        std::cout << "transitive_left_remainder_checks " << remainder_checks << '\n';
        std::cout << "transitive_left_remainder_zero " << remainder_zero << '\n';
        std::cout << "transitive_left_interaction_checks " << interaction_checks << '\n';
        std::cout << "transitive_left_interaction_zero " << interaction_zero << '\n';
        std::cout << "transitive_right_dual_remainder_checks "
                  << dual_remainder_checks << '\n';
        std::cout << "transitive_right_dual_remainder_zero "
                  << dual_remainder_zero << '\n';
        std::cout << "transitive_left_no_sink_os_checks " << os_checks << '\n';
        std::cout << "terminal_c3_identity_checks " << c3_identity_checks << '\n';
        std::cout << "terminal_c3_os_checks " << c3_os_checks << '\n';
        std::cout << "terminal_c3_os_equal " << c3_os_equal << '\n';
        std::cout << "cycle_first_negative_remainder_checks "
                  << negative_remainder_checks << '\n';
        std::cout << "cycle_first_negative_interaction_checks "
                  << negative_interaction_checks << '\n';
        std::cout << "minimal_cycle_first_theta_hostile "
                  << theta(cycle, one, one) << '\n';
        std::cout << "PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}
