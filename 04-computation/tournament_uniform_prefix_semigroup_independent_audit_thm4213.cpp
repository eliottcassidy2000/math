#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using i128 = __int128_t;

struct Tournament {
    std::vector<std::uint64_t> out;
};

struct Data {
    Tournament tournament;
    std::int64_t hamilton;
    std::int64_t mass;
    i128 gate;
};

static void need(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

static std::string decimal(i128 value) {
    if (value == 0) {
        return "0";
    }
    const bool negative = value < 0;
    if (negative) {
        value = -value;
    }
    std::string result;
    while (value > 0) {
        result.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (negative) {
        result.push_back('-');
    }
    std::reverse(result.begin(), result.end());
    return result;
}

static Tournament transitive(int order) {
    Tournament answer{std::vector<std::uint64_t>(static_cast<std::size_t>(order), 0)};
    for (int left = 0; left < order; ++left) {
        for (int right = left + 1; right < order; ++right) {
            answer.out[static_cast<std::size_t>(left)] |= std::uint64_t{1} << right;
        }
    }
    return answer;
}

static Tournament cycle_three() {
    Tournament answer{std::vector<std::uint64_t>(3, 0)};
    answer.out[0] |= std::uint64_t{1} << 1;
    answer.out[1] |= std::uint64_t{1} << 2;
    answer.out[2] |= std::uint64_t{1} << 0;
    return answer;
}

static Tournament ordinal(const Tournament& left, const Tournament& right) {
    const int nl = static_cast<int>(left.out.size());
    const int nr = static_cast<int>(right.out.size());
    need(nl + nr < 63, "ordinal tournament bit width");
    Tournament answer{std::vector<std::uint64_t>(static_cast<std::size_t>(nl + nr), 0)};
    const std::uint64_t right_mask = ((std::uint64_t{1} << nr) - 1) << nl;
    for (int vertex = 0; vertex < nl; ++vertex) {
        answer.out[static_cast<std::size_t>(vertex)] =
            left.out[static_cast<std::size_t>(vertex)] | right_mask;
    }
    for (int vertex = 0; vertex < nr; ++vertex) {
        answer.out[static_cast<std::size_t>(nl + vertex)] =
            right.out[static_cast<std::size_t>(vertex)] << nl;
    }
    return answer;
}

static std::string key(const Tournament& tournament) {
    std::string answer = std::to_string(tournament.out.size()) + ":";
    for (const std::uint64_t row : tournament.out) {
        answer += std::to_string(row) + ",";
    }
    return answer;
}

static Data analyze(const Tournament& tournament) {
    const int n = static_cast<int>(tournament.out.size());
    need(n >= 1 && n <= 20, "literal audit order range");
    const std::size_t size = std::size_t{1} << n;
    const std::size_t full = size - 1;
    std::vector<std::int64_t> starts(size * static_cast<std::size_t>(n), 0);
    std::vector<std::int64_t> ends(size * static_cast<std::size_t>(n), 0);
    const auto at = [n](std::size_t mask, int vertex) {
        return mask * static_cast<std::size_t>(n) + static_cast<std::size_t>(vertex);
    };
    for (int vertex = 0; vertex < n; ++vertex) {
        starts[at(std::size_t{1} << vertex, vertex)] = 1;
        ends[at(std::size_t{1} << vertex, vertex)] = 1;
    }
    for (std::size_t mask = 1; mask < size; ++mask) {
        if ((mask & (mask - 1)) == 0) {
            continue;
        }
        for (int vertex = 0; vertex < n; ++vertex) {
            const std::size_t bit = std::size_t{1} << vertex;
            if ((mask & bit) == 0) {
                continue;
            }
            const std::size_t rest = mask ^ bit;
            std::int64_t start_count = 0;
            std::int64_t end_count = 0;
            for (int other = 0; other < n; ++other) {
                const std::size_t other_bit = std::size_t{1} << other;
                if ((rest & other_bit) == 0) {
                    continue;
                }
                if ((tournament.out[static_cast<std::size_t>(vertex)] & other_bit) != 0) {
                    start_count += starts[at(rest, other)];
                }
                if ((tournament.out[static_cast<std::size_t>(other)] & bit) != 0) {
                    end_count += ends[at(rest, other)];
                }
            }
            starts[at(mask, vertex)] = start_count;
            ends[at(mask, vertex)] = end_count;
        }
    }
    std::int64_t hamilton = 0;
    for (int vertex = 0; vertex < n; ++vertex) {
        hamilton += ends[at(full, vertex)];
    }

    std::vector<std::int64_t> exposed(static_cast<std::size_t>(n * n), 0);
    for (std::size_t mask = 1; mask < full; ++mask) {
        const std::size_t complement = full ^ mask;
        for (int left = 0; left < n; ++left) {
            if ((mask & (std::size_t{1} << left)) == 0) {
                continue;
            }
            const std::int64_t left_count = ends[at(mask, left)];
            if (left_count == 0) {
                continue;
            }
            for (int right = 0; right < n; ++right) {
                if ((complement & (std::size_t{1} << right)) != 0) {
                    exposed[static_cast<std::size_t>(left * n + right)] +=
                        left_count * starts[at(complement, right)];
                }
            }
        }
    }

    std::vector<std::int64_t> capacity(static_cast<std::size_t>(n * n), 0);
    std::vector<std::int64_t> degree(static_cast<std::size_t>(n), 0);
    std::vector<std::int64_t> current(static_cast<std::size_t>(n), 0);
    std::int64_t mass = 0;
    i128 squares = 0;
    for (int left = 0; left < n; ++left) {
        for (int right = left + 1; right < n; ++right) {
            const std::int64_t value =
                exposed[static_cast<std::size_t>(left * n + right)]
                + exposed[static_cast<std::size_t>(right * n + left)];
            capacity[static_cast<std::size_t>(left * n + right)] = value;
            capacity[static_cast<std::size_t>(right * n + left)] = value;
            degree[static_cast<std::size_t>(left)] += value;
            degree[static_cast<std::size_t>(right)] += value;
            mass += value;
            squares += static_cast<i128>(value) * value;
            if ((tournament.out[static_cast<std::size_t>(left)]
                 & (std::uint64_t{1} << right)) != 0) {
                current[static_cast<std::size_t>(left)] += value;
                current[static_cast<std::size_t>(right)] -= value;
            } else {
                current[static_cast<std::size_t>(left)] -= value;
                current[static_cast<std::size_t>(right)] += value;
            }
        }
    }
    i128 degree_squares = 0;
    i128 directed = 0;
    for (int vertex = 0; vertex < n; ++vertex) {
        degree_squares += static_cast<i128>(degree[static_cast<std::size_t>(vertex)])
            * degree[static_cast<std::size_t>(vertex)];
        directed += static_cast<i128>(degree[static_cast<std::size_t>(vertex)])
            * current[static_cast<std::size_t>(vertex)];
    }
    const i128 disjoint =
        (static_cast<i128>(mass) * mass + squares - degree_squares) / 2;
    return Data{tournament, hamilton, mass, disjoint + 2 * directed};
}

static std::map<std::string, Data> cache;

static const Data& data(const Tournament& tournament) {
    const std::string id = key(tournament);
    const auto found = cache.find(id);
    if (found != cache.end()) {
        return found->second;
    }
    const auto inserted = cache.emplace(id, analyze(tournament));
    return inserted.first->second;
}

static i128 remainder(const Tournament& left, const Tournament& right) {
    const Data& left_data = data(left);
    const Data& right_data = data(right);
    const Data& child_data = data(ordinal(left, right));
    return child_data.gate
        - static_cast<i128>(right_data.hamilton) * right_data.hamilton * left_data.gate
        - static_cast<i128>(left_data.hamilton) * left_data.hamilton * right_data.gate;
}

static i128 defect(
    const Tournament& prefix,
    const Tournament& middle,
    const Tournament& right
) {
    return remainder(ordinal(prefix, middle), right) - remainder(middle, right);
}

int main() {
    try {
        const Tournament one = transitive(1);
        const Tournament p2 = transitive(2);
        const Tournament cycle = cycle_three();
        const Tournament q5 = ordinal(cycle, transitive(5));
        const std::vector<Tournament> contexts = {one, p2, cycle};

        int floor_checks = 0;
        int equalities = 0;
        for (const Tournament& middle : contexts) {
            for (const Tournament& right : contexts) {
                const i128 value = defect(q5, middle, right);
                const Data& middle_data = data(middle);
                const Data& right_data = data(right);
                const i128 floor = static_cast<i128>(10764)
                    * middle_data.hamilton * middle_data.hamilton
                    * right_data.hamilton * right_data.hamilton;
                need(value >= floor, "literal Q5 floor control");
                if (value == floor) {
                    ++equalities;
                    need(
                        middle.out.size() == 1 && right.out.size() == 1,
                        "literal Q5 equality boundary"
                    );
                }
                ++floor_checks;
            }
        }
        need(equalities == 1, "unique literal Q5 equality");

        int telescope_checks = 0;
        const Tournament composite = ordinal(cycle, p2);
        for (const Tournament& middle : contexts) {
            for (const Tournament& right : contexts) {
                const i128 lhs = defect(composite, middle, right);
                const i128 rhs = defect(cycle, ordinal(p2, middle), right)
                    + defect(p2, middle, right);
                need(lhs == rhs, "literal contextual-defect telescope");
                ++telescope_checks;
            }
        }

        const Tournament hostile = ordinal(cycle, q5);
        const i128 hostile_value = defect(hostile, one, one);
        need(hostile.out.size() == 11, "literal hostile order");
        need(data(hostile).hamilton == 9, "literal hostile Hamilton count");
        need(hostile_value == static_cast<i128>(-338580), "literal hostile value");

        const Tournament q6 = ordinal(q5, one);
        const i128 osplus_control = remainder(q6, cycle);
        need(osplus_control > 0, "literal OS+ padding control");

        std::cout << "theorem=THM-4213-independent\n";
        std::cout << "engine=literal-labelled-subset-DP-no-ordinal-capacity-transfer\n";
        std::cout << "literal_q5_floor_checks=" << floor_checks << "\n";
        std::cout << "literal_q5_equalities=" << equalities << "\n";
        std::cout << "literal_telescope_checks=" << telescope_checks << "\n";
        std::cout << "good_suffix_hostile_order=" << hostile.out.size() << "\n";
        std::cout << "good_suffix_hostile_hamilton=" << data(hostile).hamilton << "\n";
        std::cout << "good_suffix_hostile_singleton_defect="
                  << decimal(hostile_value) << "\n";
        std::cout << "q6_against_c3_remainder=" << decimal(osplus_control) << "\n";
        std::cout << "status=PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << "\n";
        return 1;
    }
}
