// Independent literal-permutation audit for THM-4184.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
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
            answer.push_back((out[i] & (std::uint32_t{1} << j)) ? '1' : '0');
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

static bool has_sink(const Out& out) {
    return std::any_of(out.begin(), out.end(), [](std::uint32_t row) { return row == 0; });
}

static bool path_valid(const Out& out, const std::vector<int>& path) {
    for (std::size_t i = 0; i + 1 < path.size(); ++i) {
        if (!(out[static_cast<std::size_t>(path[i])] &
              (std::uint32_t{1} << path[i + 1]))) return false;
    }
    return true;
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
        Vec2 state{hamilton, 0};
        for (const auto& row : starts) {
            state[0] += row[0];
            state[1] += row[1];
        }
        return state;
    }
};

static i64 gate(const Out& out, const std::vector<std::vector<i64>>& capacities) {
    const int n = static_cast<int>(out.size());
    std::vector<i64> degree(static_cast<std::size_t>(n), 0);
    std::vector<i64> current(static_cast<std::size_t>(n), 0);
    i64 mass = 0;
    i64 squares = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const i64 value = capacities[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
            mass += value;
            squares += value * value;
            degree[static_cast<std::size_t>(i)] += value;
            degree[static_cast<std::size_t>(j)] += value;
            if (out[static_cast<std::size_t>(i)] & (std::uint32_t{1} << j)) {
                current[static_cast<std::size_t>(i)] += value;
                current[static_cast<std::size_t>(j)] -= value;
            } else {
                current[static_cast<std::size_t>(i)] -= value;
                current[static_cast<std::size_t>(j)] += value;
            }
        }
    }
    i64 degree_squares = 0;
    i64 signed_current = 0;
    for (int i = 0; i < n; ++i) {
        degree_squares += degree[static_cast<std::size_t>(i)] * degree[static_cast<std::size_t>(i)];
        signed_current += degree[static_cast<std::size_t>(i)] * current[static_cast<std::size_t>(i)];
    }
    const i64 numerator = mass * mass + squares - degree_squares;
    need(numerator % 2 == 0, "disjoint parity");
    return numerator / 2 + 2 * signed_current;
}

static Data literal_data(const Out& out) {
    const int n = static_cast<int>(out.size());
    const int size = 1 << n;
    const int full = size - 1;
    std::vector<i64> hamilton(static_cast<std::size_t>(size), 0);
    std::vector<std::vector<i64>> starts_by_subset(
        static_cast<std::size_t>(size), std::vector<i64>(static_cast<std::size_t>(n), 0));
    std::vector<std::vector<i64>> ends_by_subset(
        static_cast<std::size_t>(size), std::vector<i64>(static_cast<std::size_t>(n), 0));
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

    std::vector<Vec2> rooted_start(static_cast<std::size_t>(n), Vec2{0, 0});
    std::vector<Vec2> rooted_end(static_cast<std::size_t>(n), Vec2{0, 0});
    for (int mask = 1; mask < size; ++mask) {
        const int parity = __builtin_popcount(static_cast<unsigned>(mask)) & 1;
        const i64 complement = hamilton[static_cast<std::size_t>(full ^ mask)];
        for (int vertex = 0; vertex < n; ++vertex) {
            rooted_start[static_cast<std::size_t>(vertex)][static_cast<std::size_t>(parity)] +=
                starts_by_subset[static_cast<std::size_t>(mask)][static_cast<std::size_t>(vertex)] * complement;
            rooted_end[static_cast<std::size_t>(vertex)][static_cast<std::size_t>(parity)] +=
                ends_by_subset[static_cast<std::size_t>(mask)][static_cast<std::size_t>(vertex)] * complement;
        }
    }

    std::vector<std::vector<i64>> exposed(
        static_cast<std::size_t>(n), std::vector<i64>(static_cast<std::size_t>(n), 0));
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
        static_cast<std::size_t>(n), std::vector<i64>(static_cast<std::size_t>(n), 0));
    i64 mass = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const i64 value = exposed[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] +
                              exposed[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)];
            capacities[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = value;
            capacities[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] = value;
            mass += value;
        }
    }
    return Data{out, hamilton[static_cast<std::size_t>(full)], capacities,
                rooted_start, rooted_end, mass, gate(out, capacities)};
}

static Out ordinal_out(const Out& left, const Out& right) {
    const int nl = static_cast<int>(left.size());
    const int nr = static_cast<int>(right.size());
    Out answer(static_cast<std::size_t>(nl + nr), 0);
    const std::uint32_t right_mask = ((std::uint32_t{1} << nr) - 1) << nl;
    for (int i = 0; i < nl; ++i) answer[static_cast<std::size_t>(i)] = left[static_cast<std::size_t>(i)] | right_mask;
    for (int j = 0; j < nr; ++j) answer[static_cast<std::size_t>(nl + j)] = right[static_cast<std::size_t>(j)] << nl;
    return answer;
}

static Data compose(const Data& left, const Data& right) {
    const int nl = static_cast<int>(left.out.size());
    const int nr = static_cast<int>(right.out.size());
    const int n = nl + nr;
    std::vector<std::vector<i64>> capacities(
        static_cast<std::size_t>(n), std::vector<i64>(static_cast<std::size_t>(n), 0));
    for (int i = 0; i < nl; ++i) {
        for (int j = i + 1; j < nl; ++j) {
            const i64 value = right.hamilton * left.capacities[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
            capacities[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = value;
            capacities[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] = value;
        }
    }
    for (int i = 0; i < nr; ++i) {
        for (int j = i + 1; j < nr; ++j) {
            const i64 value = left.hamilton * right.capacities[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
            capacities[static_cast<std::size_t>(nl + i)][static_cast<std::size_t>(nl + j)] = value;
            capacities[static_cast<std::size_t>(nl + j)][static_cast<std::size_t>(nl + i)] = value;
        }
    }
    for (int i = 0; i < nl; ++i) {
        for (int j = 0; j < nr; ++j) {
            const i64 value = 2 * (
                left.starts[static_cast<std::size_t>(i)][0] * right.ends[static_cast<std::size_t>(j)][0] +
                left.starts[static_cast<std::size_t>(i)][1] * right.ends[static_cast<std::size_t>(j)][1]);
            capacities[static_cast<std::size_t>(i)][static_cast<std::size_t>(nl + j)] = value;
            capacities[static_cast<std::size_t>(nl + j)][static_cast<std::size_t>(i)] = value;
        }
    }
    const Vec2 left_optional = left.optional_state();
    const Vec2 right_optional = right.optional_state();
    std::vector<Vec2> starts;
    std::vector<Vec2> ends;
    starts.reserve(static_cast<std::size_t>(n));
    ends.reserve(static_cast<std::size_t>(n));
    for (const Vec2& row : left.starts) starts.push_back(convolution(row, right_optional));
    for (const Vec2& row : right.starts) starts.push_back({left.hamilton * row[0], left.hamilton * row[1]});
    for (const Vec2& row : left.ends) ends.push_back({right.hamilton * row[0], right.hamilton * row[1]});
    for (const Vec2& row : right.ends) ends.push_back(convolution(left_optional, row));
    i64 mass = 0;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) mass += capacities[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
    const Out out = ordinal_out(left.out, right.out);
    return Data{out, left.hamilton * right.hamilton, capacities, starts, ends,
                mass, gate(out, capacities)};
}

static i64 remainder(const Data& left, const Data& right) {
    const Data child = compose(left, right);
    return child.gate - right.hamilton * right.hamilton * left.gate -
           left.hamilton * left.hamilton * right.gate;
}

using PathPair = std::pair<std::vector<int>, std::vector<int>>;

static PathPair terminal_move(const Out& out, const PathPair& pair) {
    auto left = pair.first;
    auto right = pair.second;
    if (left.empty()) {
        left.push_back(right.back());
        right.pop_back();
    } else if (right.empty()) {
        right.push_back(left.back());
        left.pop_back();
    } else if (out[static_cast<std::size_t>(left.back())] &
               (std::uint32_t{1} << right.back())) {
        left.push_back(right.back());
        right.pop_back();
    } else {
        right.push_back(left.back());
        left.pop_back();
    }
    return {left, right};
}

static i64 involution_audit(const Out& out) {
    const int n = static_cast<int>(out.size());
    std::vector<int> word(static_cast<std::size_t>(n));
    std::iota(word.begin(), word.end(), 0);
    i64 objects = 0;
    i64 signed_sum = 0;
    do {
        for (int cut = 0; cut <= n; ++cut) {
            PathPair pair{
                std::vector<int>(word.begin(), word.begin() + cut),
                std::vector<int>(word.begin() + cut, word.end()),
            };
            if (!path_valid(out, pair.first) || !path_valid(out, pair.second)) continue;
            const PathPair image = terminal_move(out, pair);
            need(path_valid(out, image.first) && path_valid(out, image.second),
                 "literal terminal move path preservation");
            need(terminal_move(out, image) == pair, "literal terminal move involution");
            need(((pair.first.size() ^ image.first.size()) & 1U) != 0U,
                 "literal terminal move parity");
            signed_sum += (pair.first.size() & 1U) ? -1 : 1;
            ++objects;
        }
    } while (std::next_permutation(word.begin(), word.end()));
    need(signed_sum == 0, "literal two-path-cover balance");
    return objects;
}

static bool data_equal(const Data& left, const Data& right) {
    return left.out == right.out && left.hamilton == right.hamilton &&
           left.capacities == right.capacities && left.starts == right.starts &&
           left.ends == right.ends && left.mass == right.mass && left.gate == right.gate;
}

static i64 transitive_gate_formula(int n) {
    if (n == 1) return 0;
    const i64 p2 = i64{1} << n;
    const i64 p4 = i64{1} << (2 * n);
    return (4 * p4 - 9 * (n + 2) * p2 + 24 * n + 32) / 18;
}

static i64 lollipop_gate_formula(int n) {
    if (n == 0) return 0;
    const i64 p2 = i64{1} << n;
    const i64 p4 = i64{1} << (2 * n);
    return 2 * (37 * p4 - 9 * (n + 4) * p2 + 6 * n - 4);
}

static Data capacity_only_data(const Out& out, std::vector<std::vector<i64>> capacities) {
    const int n = static_cast<int>(out.size());
    i64 mass = 0;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) mass += capacities[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
    return Data{out, 0, capacities, {}, {}, mass, gate(out, capacities)};
}

static Data transitive_capacity_data(int n) {
    std::vector<std::vector<i64>> capacities(
        static_cast<std::size_t>(n), std::vector<i64>(static_cast<std::size_t>(n), 0));
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const int distance = j - i;
            const i64 value = distance == 1 ? 2 : (i64{1} << (distance - 1));
            capacities[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = value;
            capacities[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] = value;
        }
    }
    return capacity_only_data(transitive(n), std::move(capacities));
}

static Data lollipop_capacity_data(int prefix) {
    if (prefix == 0) return literal_data(parse("101", 3));
    const Data p = transitive_capacity_data(prefix);
    const Out cycle = parse("101", 3);
    const int n = prefix + 3;
    std::vector<std::vector<i64>> capacities(
        static_cast<std::size_t>(n), std::vector<i64>(static_cast<std::size_t>(n), 0));
    for (int i = 0; i < prefix; ++i) {
        for (int j = 0; j < prefix; ++j) {
            capacities[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                3 * p.capacities[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            capacities[static_cast<std::size_t>(prefix + i)][static_cast<std::size_t>(prefix + j)] = 2;
            capacities[static_cast<std::size_t>(prefix + j)][static_cast<std::size_t>(prefix + i)] = 2;
        }
    }
    for (int i = 0; i < prefix; ++i) {
        const i64 value = i == prefix - 1 ? 4 : 3 * (i64{1} << (prefix - i - 1));
        for (int j = 0; j < 3; ++j) {
            capacities[static_cast<std::size_t>(i)][static_cast<std::size_t>(prefix + j)] = value;
            capacities[static_cast<std::size_t>(prefix + j)][static_cast<std::size_t>(i)] = value;
        }
    }
    return capacity_only_data(ordinal_out(transitive(prefix), cycle), std::move(capacities));
}

static void fnv_update(u64& digest, const std::string& text) {
    constexpr u64 prime = 1099511628211ULL;
    for (unsigned char byte : text) {
        digest ^= static_cast<u64>(byte);
        digest *= prime;
    }
}

int main() {
    try {
        std::map<int, std::vector<Data>> banks;
        i64 labelled_parity_checks = 0;
        i64 involution_objects = 0;
        for (int n = 1; n <= 5; ++n) {
            for (const Out& out : all_labelled(n)) {
                Data data = literal_data(out);
                Vec2 start_total{0, 0};
                Vec2 end_total{0, 0};
                for (const Vec2& row : data.starts) {
                    start_total[0] += row[0];
                    start_total[1] += row[1];
                }
                for (const Vec2& row : data.ends) {
                    end_total[0] += row[0];
                    end_total[1] += row[1];
                }
                need(start_total == end_total, "literal start/end totals");
                need(start_total[1] - start_total[0] == data.hamilton,
                     "literal all-order parity identity");
                need(data.optional_state()[0] == data.optional_state()[1],
                     "literal optional state balance");
                involution_objects += involution_audit(out);
                ++labelled_parity_checks;
                banks[n].push_back(std::move(data));
            }
        }
        need(labelled_parity_checks == 1099, "labelled tournament count through five");

        std::vector<const Data*> factors4;
        std::vector<const Data*> no_sink4;
        for (int n = 1; n <= 4; ++n) {
            for (const Data& data : banks[n]) {
                factors4.push_back(&data);
                if (n >= 3 && !has_sink(data.out)) no_sink4.push_back(&data);
            }
        }

        i64 pair_direct_checks = 0;
        std::map<std::string, Data> direct_cache;
        auto direct_cached = [&direct_cache](const Out& out) -> const Data& {
            const std::string key = std::to_string(out.size()) + ":" + label(out);
            auto found = direct_cache.find(key);
            if (found == direct_cache.end()) {
                found = direct_cache.emplace(key, literal_data(out)).first;
            }
            return found->second;
        };
        for (const Data* left : factors4) {
            for (const Data* right : factors4) {
                if (left->out.size() + right->out.size() > 6) continue;
                const Data child = compose(*left, *right);
                need(data_equal(child, direct_cached(child.out)), "literal pair transfer");
                ++pair_direct_checks;
            }
        }

        i64 triple_direct_checks = 0;
        i64 cocycle_checks = 0;
        for (const Data* left : factors4) {
            for (const Data* middle : factors4) {
                for (const Data* right : factors4) {
                    if (left->out.size() + middle->out.size() + right->out.size() > 6) continue;
                    const Data lm = compose(*left, *middle);
                    const Data mr = compose(*middle, *right);
                    const Data first = compose(lm, *right);
                    const Data second = compose(*left, mr);
                    need(first.capacities == second.capacities, "literal sidecar associativity");
                    need(data_equal(first, direct_cached(first.out)), "literal triple transfer");
                    const i64 lhs = remainder(lm, *right) +
                                    right->hamilton * right->hamilton * remainder(*left, *middle);
                    const i64 rhs = remainder(*left, mr) +
                                    left->hamilton * left->hamilton * remainder(*middle, *right);
                    need(lhs == rhs, "literal weighted cocycle");
                    ++triple_direct_checks;
                    ++cocycle_checks;
                }
            }
        }

        std::map<std::pair<const Data*, const Data*>, i64> right_remainders;
        for (const Data* middle : factors4) {
            for (const Data* right : no_sink4) {
                right_remainders[{middle, right}] = remainder(*middle, *right);
            }
        }
        i64 theta_checks = 0;
        i64 theta_negative = 0;
        i64 theta_zero = 0;
        i64 minimum_theta = 0;
        i64 minimum_denominator = 1;
        std::tuple<int, int, int, std::string, std::string, std::string, i64, i64>
            minimum_witness;
        bool have_minimum = false;
        u64 theta_digest = 14695981039346656037ULL;
        for (const Data* left : factors4) {
            for (const Data* middle : factors4) {
                const Data lm = compose(*left, *middle);
                for (const Data* right : no_sink4) {
                    const i64 theta = remainder(lm, *right) -
                        left->hamilton * left->hamilton * right_remainders[{middle, right}];
                    const i64 denominator = left->hamilton * left->hamilton *
                        middle->hamilton * middle->hamilton *
                        right->hamilton * right->hamilton;
                    if (!have_minimum ||
                        theta * minimum_denominator < minimum_theta * denominator) {
                        have_minimum = true;
                        minimum_theta = theta;
                        minimum_denominator = denominator;
                        minimum_witness = {
                            static_cast<int>(left->out.size()),
                            static_cast<int>(middle->out.size()),
                            static_cast<int>(right->out.size()),
                            label(left->out), label(middle->out), label(right->out),
                            theta, denominator,
                        };
                    }
                    theta_negative += theta < 0;
                    theta_zero += theta == 0;
                    std::ostringstream row;
                    row << left->out.size() << '|' << middle->out.size() << '|'
                        << right->out.size() << '|' << label(left->out) << '|'
                        << label(middle->out) << '|' << label(right->out) << '|'
                        << theta << '|' << denominator << '\n';
                    fnv_update(theta_digest, row.str());
                    ++theta_checks;
                }
            }
        }
        need(theta_negative == 0 && theta_zero == 0, "labelled no-sink-third Theta");

        std::vector<const Data*> factors3;
        for (int n = 1; n <= 3; ++n) for (const Data& data : banks[n]) factors3.push_back(&data);
        i64 hostile_checks = 0;
        i64 hostile_negative = 0;
        i64 hostile_zero = 0;
        std::tuple<int, int, int, int, std::string, std::string, std::string, i64>
            first_hostile;
        bool have_hostile = false;
        for (const Data* left : factors3) {
            for (const Data* middle : factors3) {
                for (const Data* right : factors3) {
                    const int total = static_cast<int>(left->out.size() + middle->out.size() + right->out.size());
                    if (total > 5) continue;
                    const i64 theta = remainder(compose(*left, *middle), *right) -
                        left->hamilton * left->hamilton * remainder(*middle, *right);
                    const auto record = std::make_tuple(
                        total, static_cast<int>(left->out.size()),
                        static_cast<int>(middle->out.size()), static_cast<int>(right->out.size()),
                        label(left->out), label(middle->out), label(right->out), theta);
                    if (theta < 0) {
                        ++hostile_negative;
                        if (!have_hostile || record < first_hostile) {
                            have_hostile = true;
                            first_hostile = record;
                        }
                    } else if (theta == 0) {
                        ++hostile_zero;
                    }
                    ++hostile_checks;
                }
            }
        }
        need(have_hostile, "overstrong hostile exists");
        need(std::get<0>(first_hostile) == 5 && std::get<1>(first_hostile) == 3 &&
             std::get<2>(first_hostile) == 1 && std::get<3>(first_hostile) == 1 &&
             std::get<7>(first_hostile) == -216, "minimal hostile type/value");

        i64 formula_checks = 0;
        for (int n = 1; n <= 20; ++n) {
            need(transitive_capacity_data(n).gate == transitive_gate_formula(n),
                 "independent transitive formula");
            need(lollipop_capacity_data(n).gate == lollipop_gate_formula(n),
                 "independent lollipop formula");
            formula_checks += 2;
        }
        i64 previous_increment = -1;
        for (int n = 0; n < 20; ++n) {
            const i64 increment = lollipop_gate_formula(n + 1) - lollipop_gate_formula(n);
            need(increment > 0 && increment > previous_increment, "independent strict convexity");
            previous_increment = increment;
        }
        i64 lollipop_checks = 0;
        i64 lollipop_minimum = std::numeric_limits<i64>::max();
        int lollipop_min_a = -1;
        int lollipop_min_b = -1;
        for (int a = 1; a <= 10; ++a) {
            for (int b = 0; b <= 10; ++b) {
                const i64 value = lollipop_gate_formula(a + b) -
                    9 * transitive_gate_formula(a) - lollipop_gate_formula(b);
                need(value > 0, "independent lollipop OS+ remainder");
                if (value < lollipop_minimum) {
                    lollipop_minimum = value;
                    lollipop_min_a = a;
                    lollipop_min_b = b;
                }
                ++lollipop_checks;
            }
        }
        need(lollipop_minimum == 120 && lollipop_min_a == 1 && lollipop_min_b == 0,
             "independent lollipop minimum");

        std::cout << "TOURNAMENT_ORDINAL_COCYCLE_PARITY_INDEPENDENT\n";
        std::cout << "labelled_parity_tournaments " << labelled_parity_checks << '\n';
        std::cout << "terminal_involution_cover_objects " << involution_objects << '\n';
        std::cout << "pair_direct_transfer_checks " << pair_direct_checks << '\n';
        std::cout << "triple_direct_transfer_checks " << triple_direct_checks << '\n';
        std::cout << "weighted_cocycle_checks " << cocycle_checks << '\n';
        std::cout << "labelled_no_sink_third_presentations " << theta_checks << '\n';
        std::cout << "theta_negative " << theta_negative << " zero " << theta_zero << '\n';
        std::cout << "theta_normalized_minimum " << minimum_theta << '/' << minimum_denominator
                  << " witness " << std::get<0>(minimum_witness) << ','
                  << std::get<1>(minimum_witness) << ',' << std::get<2>(minimum_witness)
                  << "," << std::get<3>(minimum_witness) << ','
                  << std::get<4>(minimum_witness) << ',' << std::get<5>(minimum_witness)
                  << ',' << std::get<6>(minimum_witness) << ','
                  << std::get<7>(minimum_witness) << '\n';
        std::cout << "theta_stream_fnv1a64 " << std::hex << theta_digest << std::dec << '\n';
        std::cout << "minimal_hostile_scope_checks " << hostile_checks << '\n';
        std::cout << "minimal_hostile_negative " << hostile_negative << " zero " << hostile_zero << '\n';
        std::cout << "first_hostile " << std::get<0>(first_hostile) << ','
                  << std::get<1>(first_hostile) << ',' << std::get<2>(first_hostile)
                  << ',' << std::get<3>(first_hostile) << ',' << std::get<4>(first_hostile)
                  << ',' << std::get<5>(first_hostile) << ',' << std::get<6>(first_hostile)
                  << ',' << std::get<7>(first_hostile) << '\n';
        std::cout << "closed_form_gate_checks " << formula_checks << '\n';
        std::cout << "lollipop_box_checks " << lollipop_checks << '\n';
        std::cout << "lollipop_remainder_minimum " << lollipop_minimum << ','
                  << lollipop_min_a << ',' << lollipop_min_b << '\n';
        std::cout << "PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}
