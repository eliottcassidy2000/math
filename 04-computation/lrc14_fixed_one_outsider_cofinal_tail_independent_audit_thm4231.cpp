// Independent exact referee for THM-4231's fixed q=1 outsider ray.
//
// Geometry is reconstructed by an endpoint-event sweep.  The exhaustive
// body census uses a right-half subset transform keyed only by the newly
// reconstructed coefficient ledger; it imports no body rows.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <omp.h>

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

namespace {

constexpr std::array<int, 30> POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
};
constexpr int BODY_SIZE = 9;
constexpr u32 HALF_MASK = (u32{1} << 15) - 1;

[[noreturn]] void fail(const std::string& message) {
    std::cerr << "FAIL " << message << '\n';
    std::exit(1);
}

void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

std::string decimal(u128 value) {
    if (value == 0) return "0";
    std::string answer;
    while (value != 0) {
        answer.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    std::reverse(answer.begin(), answer.end());
    return answer;
}

u64 ceiling(u128 numerator, u128 denominator) {
    require(denominator != 0, "zero denominator");
    const u128 result = (numerator + denominator - 1) / denominator;
    require(result <= std::numeric_limits<u64>::max(), "ceiling overflow");
    return static_cast<u64>(result);
}

u64 splitmix64(u64 x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

struct Digest {
    u64 xor_value = 0;
    u64 sum_value = 0;
    void add(u32 body, u64 threshold, i64 mass, i64 components, i64 surplus) {
        u64 h = splitmix64(body);
        h = splitmix64(h ^ threshold);
        h = splitmix64(h ^ static_cast<u64>(mass));
        h = splitmix64(h ^ static_cast<u64>(components));
        h = splitmix64(h ^ static_cast<u64>(surplus));
        xor_value ^= h;
        sum_value += h;
    }
};

struct Event {
    i64 position;
    int vertex;
    bool enter_safe;
};

struct Cell {
    i64 width;
    u32 failure;
    bool q1_safe;
};

struct Geometry {
    i64 denominator;
    std::vector<Cell> cells;
    std::size_t walls;
};

Geometry endpoint_sweep() {
    i64 denominator = 1;
    for (int p : POOL) denominator = std::lcm(denominator, i64{14} * p);
    denominator = std::lcm(denominator, i64{14});

    // vertex 0..29 denotes a pool label; vertex 30 is the fixed label 1.
    std::vector<Event> events;
    std::vector<i64> walls = {0, denominator};
    auto append_speed = [&](int speed, int vertex) {
        const i64 unit = denominator / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            const i64 enter = (14LL * tooth + 1) * unit;
            const i64 leave = (14LL * tooth + 13) * unit;
            events.push_back({enter, vertex, true});
            events.push_back({leave, vertex, false});
            walls.push_back(enter);
            walls.push_back(leave);
        }
    };
    for (int vertex = 0; vertex < 30; ++vertex) {
        append_speed(POOL[vertex], vertex);
    }
    append_speed(1, 30);

    std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
        return std::tie(a.position, a.vertex, a.enter_safe) <
               std::tie(b.position, b.vertex, b.enter_safe);
    });
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());

    u32 failure = (u32{1} << 30) - 1;
    bool q1_safe = false;
    std::size_t event_index = 0;
    std::vector<Cell> cells;
    cells.reserve(walls.size() - 1);
    for (std::size_t wall_index = 0; wall_index + 1 < walls.size(); ++wall_index) {
        const i64 left = walls[wall_index];
        while (event_index < events.size() && events[event_index].position == left) {
            const Event& event = events[event_index++];
            if (event.vertex == 30) {
                q1_safe = event.enter_safe;
            } else if (event.enter_safe) {
                failure &= ~(u32{1} << event.vertex);
            } else {
                failure |= u32{1} << event.vertex;
            }
        }
        require(walls[wall_index + 1] > left, "nonpositive cell");
        cells.push_back({walls[wall_index + 1] - left, failure, q1_safe});
    }
    require(event_index == events.size(), "unconsumed events");
    return {denominator, std::move(cells), walls.size()};
}

struct Coefficient {
    u32 failure;
    i64 mass;
    i64 components;
};

std::vector<Coefficient> coefficient_ledger(const std::vector<Cell>& cells) {
    std::map<u32, std::pair<i64, i64>> ledger;
    for (std::size_t i = 0; i < cells.size(); ++i) {
        const Cell& current = cells[i];
        const Cell& previous = cells[(i + cells.size() - 1) % cells.size()];
        if (!current.q1_safe) continue;
        if (std::popcount(current.failure) <= 30 - BODY_SIZE) {
            ledger[current.failure].first += current.width;
            ledger[current.failure].second += 1;
        }
        if (previous.q1_safe) {
            const u32 joined = previous.failure | current.failure;
            if (std::popcount(joined) <= 30 - BODY_SIZE) {
                ledger[joined].second -= 1;
            }
        }
    }

    std::vector<Coefficient> answer;
    for (const auto& [failure, pair] : ledger) {
        if (pair.first != 0 || pair.second != 0) {
            answer.push_back({failure, pair.first, pair.second});
        }
    }
    return answer;
}

std::vector<u32> fixed_weight_masks(int bits, int weight) {
    std::vector<u32> answer;
    if (weight < 0 || weight > bits) return answer;
    if (weight == 0) {
        answer.push_back(0);
        return answer;
    }
    u32 mask = (u32{1} << weight) - 1;
    const u32 limit = u32{1} << bits;
    while (mask < limit) {
        answer.push_back(mask);
        const u32 low = mask & (~mask + 1);
        const u32 high = mask + low;
        mask = high | (((high ^ mask) >> 2) / low);
    }
    return answer;
}

std::string labels(u32 body) {
    std::ostringstream out;
    out << '{';
    bool first = true;
    for (int i = 0; i < 30; ++i) {
        if ((body & (u32{1} << i)) == 0) continue;
        if (!first) out << ',';
        out << POOL[i];
        first = false;
    }
    out << '}';
    return out.str();
}

struct ExactSafeSet {
    i64 denominator;
    i64 mass;
    i64 components;
    std::size_t walls;
};

ExactSafeSet literal_geometry(const std::vector<int>& speeds) {
    i64 denominator = 1;
    for (int p : speeds) denominator = std::lcm(denominator, i64{14} * p);
    std::vector<i64> walls = {0, denominator};
    for (int p : speeds) {
        const i64 unit = denominator / (14LL * p);
        for (int tooth = 0; tooth < p; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    }
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    std::vector<unsigned char> safe(walls.size() - 1, 0);
    i64 mass = 0;
    for (std::size_t i = 0; i + 1 < walls.size(); ++i) {
        const i64 twice_midpoint = walls[i] + walls[i + 1];
        const i64 period = 2 * denominator;
        bool good = true;
        for (int p : speeds) {
            const i64 residue = static_cast<i64>((i128{p} * twice_midpoint) % period);
            if (7 * residue < denominator || 7 * residue > 13 * denominator) {
                good = false;
                break;
            }
        }
        safe[i] = static_cast<unsigned char>(good);
        if (good) mass += walls[i + 1] - walls[i];
    }
    i64 components = 0;
    for (std::size_t i = 0; i < safe.size(); ++i) {
        if (safe[i] && !safe[(i + safe.size() - 1) % safe.size()]) ++components;
    }
    return {denominator, mass, components, walls.size()};
}

struct ThreadSummary {
    u64 bodies = 0;
    u64 positive_surpluses = 0;
    u64 minimum_threshold = std::numeric_limits<u64>::max();
    u64 maximum_threshold = 0;
    u64 maximum_count = 0;
    u64 pass_542 = 0;
    u64 pass_543 = 0;
    i64 minimum_mass = std::numeric_limits<i64>::max();
    i64 maximum_mass = 0;
    i64 minimum_components = std::numeric_limits<i64>::max();
    i64 maximum_components = 0;
    u32 maximum_witness = 0;
    u32 minimum_witness = 0;
    std::map<u64, u64> histogram;
    Digest digest;
};

}  // namespace

int main() {
    const Geometry geometry = endpoint_sweep();
    const std::vector<Coefficient> coefficients = coefficient_ledger(geometry.cells);

    std::array<std::vector<u32>, BODY_SIZE + 1> masks_by_weight;
    for (int weight = 0; weight <= BODY_SIZE; ++weight) {
        masks_by_weight[weight] = fixed_weight_masks(15, weight);
    }
    std::vector<u32> left_masks;
    for (int weight = 0; weight <= BODY_SIZE; ++weight) {
        left_masks.insert(left_masks.end(), masks_by_weight[weight].begin(), masks_by_weight[weight].end());
    }
    // A nonnumeric ordering is deliberate: this is a census partition, not a
    // borrowed body-row order.
    std::sort(left_masks.begin(), left_masks.end(), [](u32 a, u32 b) {
        return std::pair{std::popcount(a), a} < std::pair{std::popcount(b), b};
    });

    const int thread_count = omp_get_max_threads();
    std::vector<ThreadSummary> summaries(thread_count);

#pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        ThreadSummary& summary = summaries[tid];
        std::vector<i64> mass(HALF_MASK + 1);
        std::vector<i64> components(HALF_MASK + 1);

#pragma omp for schedule(guided, 8)
        for (i64 left_index = 0; left_index < static_cast<i64>(left_masks.size()); ++left_index) {
            const u32 left = left_masks[left_index];
            const int right_weight = BODY_SIZE - std::popcount(left);
            std::fill(mass.begin(), mass.end(), 0);
            std::fill(components.begin(), components.end(), 0);
            for (const Coefficient& coefficient : coefficients) {
                if ((coefficient.failure & left) != 0) continue;
                const u32 right_failure = coefficient.failure >> 15;
                mass[right_failure] += coefficient.mass;
                components[right_failure] += coefficient.components;
            }
            // Down-set accumulation: entry Q becomes the sum over failure
            // masks contained in Q, hence Q=complement(body-right).
            for (int bit = 0; bit < 15; ++bit) {
                const u32 flag = u32{1} << bit;
                for (u32 upper = 0; upper <= HALF_MASK; ++upper) {
                    if (upper & flag) {
                        mass[upper] += mass[upper ^ flag];
                        components[upper] += components[upper ^ flag];
                    }
                }
            }

            for (u32 right : masks_by_weight[right_weight]) {
                const u32 body = left | (right << 15);
                const u32 query = (~right) & HALF_MASK;
                const i64 m = mass[query];
                const i64 c = components[query];
                require(m > 0 && c > 0, "invalid fixed-ray geometry");
                const i64 surplus = 54 * m - 4 * geometry.denominator;
                ++summary.bodies;
                if (surplus > 0) ++summary.positive_surpluses;
                require(surplus > 0, "nonpositive fixed-ray limit surplus");
                const u64 threshold = ceiling(
                    u128{54} * static_cast<u64>(c) * static_cast<u64>(geometry.denominator),
                    u128{7} * static_cast<u64>(surplus)
                );
                require(
                    u128{7} * static_cast<u64>(surplus) * threshold >=
                        u128{54} * static_cast<u64>(c) * static_cast<u64>(geometry.denominator),
                    "ceiling upper check"
                );
                require(
                    threshold == 0 ||
                    u128{7} * static_cast<u64>(surplus) * (threshold - 1) <
                        u128{54} * static_cast<u64>(c) * static_cast<u64>(geometry.denominator),
                    "ceiling lower check"
                );
                ++summary.histogram[threshold];
                if (threshold <= 542) ++summary.pass_542;
                if (threshold <= 543) ++summary.pass_543;
                if (threshold > summary.maximum_threshold) {
                    summary.maximum_threshold = threshold;
                    summary.maximum_count = 1;
                    summary.maximum_witness = body;
                } else if (threshold == summary.maximum_threshold) {
                    ++summary.maximum_count;
                    summary.maximum_witness = std::min(summary.maximum_witness, body);
                }
                if (threshold < summary.minimum_threshold) {
                    summary.minimum_threshold = threshold;
                    summary.minimum_witness = body;
                } else if (threshold == summary.minimum_threshold) {
                    summary.minimum_witness = std::min(summary.minimum_witness, body);
                }
                summary.minimum_mass = std::min(summary.minimum_mass, m);
                summary.maximum_mass = std::max(summary.maximum_mass, m);
                summary.minimum_components = std::min(summary.minimum_components, c);
                summary.maximum_components = std::max(summary.maximum_components, c);
                summary.digest.add(body, threshold, m, c, surplus);
            }
        }
    }

    ThreadSummary total;
    total.minimum_threshold = std::numeric_limits<u64>::max();
    total.minimum_mass = std::numeric_limits<i64>::max();
    total.minimum_components = std::numeric_limits<i64>::max();
    for (const ThreadSummary& part : summaries) {
        total.bodies += part.bodies;
        total.positive_surpluses += part.positive_surpluses;
        total.pass_542 += part.pass_542;
        total.pass_543 += part.pass_543;
        total.minimum_mass = std::min(total.minimum_mass, part.minimum_mass);
        total.maximum_mass = std::max(total.maximum_mass, part.maximum_mass);
        total.minimum_components = std::min(total.minimum_components, part.minimum_components);
        total.maximum_components = std::max(total.maximum_components, part.maximum_components);
        if (part.maximum_threshold > total.maximum_threshold) {
            total.maximum_threshold = part.maximum_threshold;
            total.maximum_count = part.maximum_count;
            total.maximum_witness = part.maximum_witness;
        } else if (part.maximum_threshold == total.maximum_threshold) {
            total.maximum_count += part.maximum_count;
            total.maximum_witness = std::min(total.maximum_witness, part.maximum_witness);
        }
        if (part.minimum_threshold < total.minimum_threshold) {
            total.minimum_threshold = part.minimum_threshold;
            total.minimum_witness = part.minimum_witness;
        } else if (part.minimum_threshold == total.minimum_threshold) {
            total.minimum_witness = std::min(total.minimum_witness, part.minimum_witness);
        }
        for (const auto& [threshold, count] : part.histogram) total.histogram[threshold] += count;
        total.digest.xor_value ^= part.digest.xor_value;
        total.digest.sum_value += part.digest.sum_value;
    }

    require(total.bodies == 14'307'150, "body universe");
    require(total.positive_surpluses == total.bodies, "positive surplus universe");
    const auto maximum_row = total.histogram.rbegin();
    const auto second_row = std::next(maximum_row);
    require(second_row != total.histogram.rend(), "missing second threshold");

    const u32 extremizer = total.maximum_witness;
    std::vector<int> extremizer_labels;
    for (int i = 0; i < 30; ++i) {
        if (extremizer & (u32{1} << i)) extremizer_labels.push_back(POOL[i]);
    }
    auto exact_body_geometry = [&](u32 body) {
        i64 m = 0;
        i64 c = 0;
        for (const Coefficient& coefficient : coefficients) {
            if ((coefficient.failure & body) == 0) {
                m += coefficient.mass;
                c += coefficient.components;
            }
        }
        return std::pair{i64{m}, i64{c}};
    };
    const auto [extremal_mass, extremal_components] = exact_body_geometry(extremizer);
    const i64 extremal_surplus = 54 * extremal_mass - 4 * geometry.denominator;

    std::vector<int> speeds_542 = extremizer_labels;
    speeds_542.push_back(1);
    speeds_542.push_back(542);
    std::sort(speeds_542.begin(), speeds_542.end());
    std::vector<int> speeds_543 = extremizer_labels;
    speeds_543.push_back(1);
    speeds_543.push_back(543);
    std::sort(speeds_543.begin(), speeds_543.end());
    const ExactSafeSet literal_542 = literal_geometry(speeds_542);
    const ExactSafeSet literal_543 = literal_geometry(speeds_543);

    std::cout << "LRC14_FIXED_Q1_RAY_INDEPENDENT_REFEREE\n";
    std::cout << "GEOMETRY DENOMINATOR " << geometry.denominator
              << " WALLS " << geometry.walls
              << " CELLS " << geometry.cells.size()
              << " COEFFICIENT_KEYS " << coefficients.size() << '\n';
    std::cout << "CENSUS BODIES " << total.bodies
              << " POSITIVE_LIMIT_SURPLUS " << total.positive_surpluses
              << " THRESHOLD_MIN " << total.minimum_threshold
              << " MIN_WITNESS " << labels(total.minimum_witness)
              << " THRESHOLD_MAX " << total.maximum_threshold
              << " MAX_COUNT " << total.maximum_count
              << " MAX_WITNESS " << labels(total.maximum_witness) << '\n';
    std::cout << "SECOND_THRESHOLD " << second_row->first
              << " COUNT " << second_row->second
              << " PASS_AT_542 " << total.pass_542
              << " PASS_AT_543 " << total.pass_543 << '\n';
    std::cout << "EXTREMAL_FIXED_GEOMETRY MASS_TICKS " << extremal_mass
              << " COMPONENTS " << extremal_components
              << " SURPLUS54_TICKS " << extremal_surplus
              << " NUMERATOR54CD "
              << decimal(u128{54} * static_cast<u64>(extremal_components) *
                         static_cast<u64>(geometry.denominator))
              << " DENOMINATOR7S "
              << decimal(u128{7} * static_cast<u64>(extremal_surplus)) << '\n';
    std::cout << "RANGES MASS_TICKS " << total.minimum_mass << ".." << total.maximum_mass
              << " COMPONENTS " << total.minimum_components << ".." << total.maximum_components << '\n';
    std::cout << "COMMUTATIVE_DIGEST XOR " << std::hex << std::setw(16)
              << std::setfill('0') << total.digest.xor_value
              << " SUM " << std::setw(16) << total.digest.sum_value
              << std::dec << std::setfill(' ') << '\n';
    auto print_literal = [](int outsider, const ExactSafeSet& set) {
        const i128 margin = i128{63} * set.mass - i128{4} * set.denominator;
        std::cout << "LITERAL R " << outsider
                  << " MASS " << set.mass << '/' << set.denominator
                  << " MARGIN63_NUM " << decimal(static_cast<u128>(margin))
                  << " COMPONENTS " << set.components
                  << " WALLS " << set.walls << '\n';
    };
    print_literal(542, literal_542);
    print_literal(543, literal_543);
    std::cout << "VERDICT CERTIFICATE_MAX " << total.maximum_threshold
              << " CERTIFICATE_SHARP YES LITERAL_SHARP "
              << ((i128{63} * literal_542.mass - i128{4} * literal_542.denominator > 0)
                      ? "NO_R542_STRICTLY_SAFE" : "NOT_DISPROVED")
              << " LRC14_OPEN\n";
    return 0;
}
