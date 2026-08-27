// Clean-room grouped-event finite audit for the fixed-50 five-petal chart
// S={10,15,20,30,264}: all 10 triple, 5 quadruple, and 1 quintuple faces.

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i32 = std::int32_t;
using i128 = __int128_t;

namespace {

constexpr int NC = 18;
constexpr std::size_t STATES = std::size_t{1} << NC;
constexpr u32 CORE_MASK = (u32{1} << NC) - 1;
constexpr std::array<int, NC> C = {
    8, 16, 40, 42, 80, 84, 85, 88, 95,
    120, 126, 143, 145, 168, 193, 240, 252, 286,
};
constexpr std::array<int, 12> O = {
    10, 15, 20, 30, 60, 63, 132, 170, 176, 190, 264, 290,
};
constexpr std::array<int, 5> S = {10, 15, 20, 30, 264};

struct Lane {
    std::array<int, 5> labels{};
    int petal_count = 0;
    int core_count = 0;
    int extra = 0;
};

std::vector<Lane> make_lanes() {
    std::vector<Lane> answer;
    answer.reserve(16);
    for (int i = 0; i < 5; ++i) for (int j = i + 1; j < 5; ++j)
        for (int k = j + 1; k < 5; ++k) {
            Lane lane;
            lane.labels = {S[i], S[j], S[k], 0, 0};
            lane.petal_count = 3;
            lane.core_count = 6;
            answer.push_back(lane);
        }
    for (int omitted = 0; omitted < 5; ++omitted) {
        Lane lane;
        int slot = 0;
        for (int j = 0; j < 5; ++j) if (j != omitted) lane.labels[slot++] = S[j];
        lane.petal_count = 4;
        lane.core_count = 5;
        answer.push_back(lane);
    }
    Lane lane;
    lane.labels = {10, 15, 20, 30, 264};
    lane.petal_count = 5;
    lane.core_count = 4;
    answer.push_back(lane);
    if (answer.size() != 16) throw std::runtime_error("lane universe");
    return answer;
}

const std::vector<Lane> LANES = make_lanes();

std::string decimal(i128 value) {
    if (value == 0) return "0";
    const bool negative = value < 0;
    if (negative) value = -value;
    std::string answer;
    while (value) {
        answer.push_back(char('0' + value % 10));
        value /= 10;
    }
    if (negative) answer.push_back('-');
    std::reverse(answer.begin(), answer.end());
    return answer;
}

u64 checked_lcm(u64 first, u64 second) {
    const u64 divisor = std::gcd(first, second);
    const i128 value = i128(first / divisor) * second;
    if (value > std::numeric_limits<u64>::max()) throw std::runtime_error("LCM overflow");
    return u64(value);
}

bool excluded_outsider(int value) {
    if (value == 50) return true;
    return std::find(C.begin(), C.end(), value) != C.end() ||
           std::find(O.begin(), O.end(), value) != O.end();
}

std::vector<u32> masks_of_size(int size) {
    std::vector<u32> answer;
    for (u32 mask = 0; mask <= CORE_MASK; ++mask) {
        if (__builtin_popcount(mask) == size) answer.push_back(mask);
    }
    return answer;
}

std::string core_body(u32 mask) {
    std::string answer;
    for (int bit = 0; bit < NC; ++bit) {
        if (!(mask & (u32{1} << bit))) continue;
        if (!answer.empty()) answer += ',';
        answer += std::to_string(C[bit]);
    }
    return answer;
}

std::string lane_labels(const Lane& lane) {
    std::string answer;
    for (int j = 0; j < lane.petal_count; ++j) {
        if (!answer.empty()) answer += ',';
        answer += std::to_string(lane.labels[j]);
    }
    return answer;
}

struct Event { u64 position; u32 toggle; };

std::vector<Event> grouped_events(const std::vector<int>& labels, u64 denominator) {
    std::vector<Event> raw;
    u64 expected = 0;
    for (int speed : labels) expected += 2 * u64(speed);
    raw.reserve(std::size_t(expected));
    for (std::size_t bit = 0; bit < labels.size(); ++bit) {
        const int speed = labels[bit];
        const u64 local_denominator = 14 * u64(speed);
        if (denominator % local_denominator) throw std::runtime_error("wall grid");
        const u64 scale = denominator / local_denominator;
        for (int tooth = 0; tooth < speed; ++tooth) {
            int left = 14 * tooth - 1;
            if (left < 0) left += int(local_denominator);
            const int right = 14 * tooth + 1;
            raw.push_back({u64(left) * scale, u32{1} << bit});
            raw.push_back({u64(right) * scale, u32{1} << bit});
        }
    }
    std::sort(raw.begin(), raw.end(), [](const Event& x, const Event& y) {
        return x.position < y.position;
    });
    std::vector<Event> answer;
    for (const Event& event : raw) {
        if (!answer.empty() && answer.back().position == event.position) {
            answer.back().toggle ^= event.toggle;
        } else {
            answer.push_back(event);
        }
    }
    answer.erase(std::remove_if(answer.begin(), answer.end(),
                                [](const Event& event) { return event.toggle == 0; }),
                 answer.end());
    return answer;
}

template <class Integer>
void subset_zeta(std::vector<Integer>& table, int lanes) {
    for (int bit = 0; bit < NC; ++bit) {
        const std::size_t step = std::size_t{1} << bit;
        for (std::size_t block = 0; block < STATES; block += 2 * step) {
            for (std::size_t offset = 0; offset < step; ++offset) {
                Integer* lower = table.data() + (block + offset) * lanes;
                Integer* upper = table.data() + (block + step + offset) * lanes;
                for (int lane = 0; lane < lanes; ++lane) upper[lane] += lower[lane];
            }
        }
    }
}

struct Geometry { u64 denominator; std::vector<Event> events; u32 initial; };

Geometry geometry(int outsider) {
    std::vector<int> labels;
    labels.insert(labels.end(), C.begin(), C.end());
    labels.insert(labels.end(), O.begin(), O.end());
    labels.push_back(50);
    if (outsider > 0) labels.push_back(outsider);
    u64 denominator = 1;
    for (int speed : labels) denominator = checked_lcm(denominator, 14 * u64(speed));
    const u32 initial = labels.size() == 32
        ? std::numeric_limits<u32>::max()
        : (u32{1} << labels.size()) - 1;
    return {denominator, grouped_events(labels, denominator), initial};
}

std::vector<u32> fixed_masks(int outsider) {
    const u32 bit50 = u32{1} << 30;
    const u32 bitr = outsider > 0 ? (u32{1} << 31) : 0;
    std::vector<u32> answer(LANES.size(), bit50 | bitr);
    for (std::size_t lane = 0; lane < LANES.size(); ++lane) {
        for (int j = 0; j < LANES[lane].petal_count; ++j) {
            const int label = LANES[lane].labels[j];
            const auto it = std::find(O.begin(), O.end(), label);
            if (it == O.end()) throw std::runtime_error("petal outside O");
            answer[lane] |= u32{1} << (NC + int(it - O.begin()));
        }
    }
    return answer;
}

struct MassTable { u64 denominator; std::vector<u64> values; };

MassTable mass_table(int outsider) {
    const Geometry geo = geometry(outsider);
    const int lanes = int(LANES.size());
    const std::vector<u32> fixed = fixed_masks(outsider);
    std::vector<u64> values(STATES * std::size_t(lanes), 0);
    u32 state = geo.initial;
    u64 previous = 0;
    auto add = [&](u64 length) {
        if (length == 0) return;
        const std::size_t row = std::size_t(state & CORE_MASK) * lanes;
        for (int lane = 0; lane < lanes; ++lane) {
            if ((state & fixed[std::size_t(lane)]) == 0) values[row + lane] += length;
        }
    };
    for (const Event& event : geo.events) {
        add(event.position - previous);
        state ^= event.toggle;
        previous = event.position;
    }
    add(geo.denominator - previous);
    if (state != geo.initial) throw std::runtime_error("mass circular mismatch");
    subset_zeta(values, lanes);
    return {geo.denominator, std::move(values)};
}

struct ComponentTable { u64 denominator; std::vector<i32> values; };

ComponentTable component_table() {
    const Geometry geo = geometry(0);
    const int lanes = int(LANES.size());
    const std::vector<u32> fixed = fixed_masks(0);
    std::vector<i32> values(STATES * std::size_t(lanes), 0);
    u32 state = geo.initial;
    for (const Event& event : geo.events) {
        const u32 before = state;
        state ^= event.toggle;
        const std::size_t start = std::size_t(state & CORE_MASK) * lanes;
        const std::size_t bridge = std::size_t((before | state) & CORE_MASK) * lanes;
        for (int lane = 0; lane < lanes; ++lane) {
            const bool before_safe = (before & fixed[std::size_t(lane)]) == 0;
            const bool after_safe = (state & fixed[std::size_t(lane)]) == 0;
            if (after_safe) ++values[start + lane];
            if (before_safe && after_safe) --values[bridge + lane];
        }
    }
    if (state != geo.initial) throw std::runtime_error("component circular mismatch");
    subset_zeta(values, lanes);
    return {geo.denominator, std::move(values)};
}

struct DirectResult { u64 numerator; u64 denominator; u64 components; };

DirectResult direct_body_sweep(std::vector<int> speeds) {
    std::sort(speeds.begin(), speeds.end());
    speeds.erase(std::unique(speeds.begin(), speeds.end()), speeds.end());
    u64 denominator = 1;
    for (int speed : speeds) denominator = checked_lcm(denominator, 14 * u64(speed));
    const std::vector<Event> events = grouped_events(speeds, denominator);
    const u32 initial = (u32{1} << speeds.size()) - 1;
    u32 state = initial;
    u64 previous = 0, numerator = 0, components = 0;
    for (const Event& event : events) {
        if (state == 0) numerator += event.position - previous;
        const u32 before = state;
        state ^= event.toggle;
        if (before != 0 && state == 0) ++components;
        previous = event.position;
    }
    if (state == 0) numerator += denominator - previous;
    if (state != initial) throw std::runtime_error("direct circular mismatch");
    return {numerator, denominator, components};
}

std::vector<int> body_speeds(u32 mask, const Lane& lane, int outsider) {
    std::vector<int> answer;
    for (int bit = 0; bit < NC; ++bit) if (mask & (u32{1} << bit)) answer.push_back(C[bit]);
    for (int j = 0; j < lane.petal_count; ++j) answer.push_back(lane.labels[j]);
    answer.push_back(50);
    if (outsider > 0) answer.push_back(outsider);
    return answer;
}

void require_same_fraction(u64 a, u64 b, u64 c, u64 d, const char* message) {
    if (i128(a) * d != i128(c) * b) throw std::runtime_error(message);
}

i128 ceil_positive(i128 numerator, i128 denominator) {
    return (numerator + denominator - 1) / denominator;
}

struct BaseResult {
    i128 min_delta = -1;
    u32 min_mask = 0;
    u64 min_mass = 0;
    i128 max_cutoff = -1;
    u32 cutoff_mask = 0;
    u64 cutoff_mass = 0;
    u64 cutoff_components = 0;
    u64 nonstrict = 0;
};

struct LiteralRow {
    bool active = false;
    i128 min_margin = -1;
    u64 denominator = 1;
    u64 min_mass = 0;
    u32 min_mask = 0;
    u64 failures = 0;
    u64 equalities = 0;
};

}  // namespace

int main() {
    const std::array<std::vector<u32>, 3> bodies = {
        masks_of_size(4), masks_of_size(5), masks_of_size(6)};
    if (bodies[0].size() != 3060 || bodies[1].size() != 8568 ||
        bodies[2].size() != 18564) throw std::runtime_error("body universes");
    const MassTable mass = mass_table(0);
    const ComponentTable components = component_table();
    if (mass.denominator != UINT64_C(91205797082400) ||
        components.denominator != mass.denominator) throw std::runtime_error("denominator");

    std::vector<BaseResult> results(LANES.size());
    const int lanes = int(LANES.size());
    for (int lane = 0; lane < lanes; ++lane) {
        const std::vector<u32>& universe = bodies[std::size_t(LANES[lane].core_count - 4)];
        for (u32 mask : universe) {
            const std::size_t index = std::size_t(CORE_MASK ^ mask) * lanes + lane;
            const u64 numerator = mass.values[index];
            const i128 delta = 54 * i128(numerator) - 4 * i128(mass.denominator);
            const i32 count = components.values[index];
            if (count <= 0) throw std::runtime_error("components");
            if (delta <= 0) ++results[lane].nonstrict;
            if (results[lane].min_delta < 0 || delta < results[lane].min_delta) {
                results[lane].min_delta = delta;
                results[lane].min_mask = mask;
                results[lane].min_mass = numerator;
            }
            if (delta > 0) {
                const i128 cutoff = ceil_positive(54 * i128(count) * mass.denominator,
                                                   7 * delta);
                if (cutoff > results[lane].max_cutoff) {
                    results[lane].max_cutoff = cutoff;
                    results[lane].cutoff_mask = mask;
                    results[lane].cutoff_mass = numerator;
                    results[lane].cutoff_components = u64(count);
                }
            }
        }
    }

    int global_cutoff = 0;
    for (const BaseResult& result : results) {
        global_cutoff = std::max(global_cutoff, int(result.max_cutoff));
    }
    std::vector<std::vector<LiteralRow>> rows(
        std::size_t(global_cutoff), std::vector<LiteralRow>(LANES.size()));
#pragma omp parallel for schedule(dynamic, 1)
    for (int outsider = 1; outsider < global_cutoff; ++outsider) {
        if (excluded_outsider(outsider)) continue;
        const MassTable table = mass_table(outsider);
        for (int lane = 0; lane < lanes; ++lane) {
            if (outsider >= results[lane].max_cutoff) continue;
            LiteralRow row;
            row.active = true;
            row.denominator = table.denominator;
            const std::vector<u32>& universe =
                bodies[std::size_t(LANES[lane].core_count - 4)];
            for (u32 mask : universe) {
                const u64 numerator =
                    table.values[std::size_t(CORE_MASK ^ mask) * lanes + lane];
                const i128 margin = 63 * i128(numerator) - 4 * i128(table.denominator);
                if (margin < 0) ++row.failures;
                if (margin == 0) ++row.equalities;
                if (row.min_margin < 0 || margin < row.min_margin) {
                    row.min_margin = margin;
                    row.min_mass = numerator;
                    row.min_mask = mask;
                }
            }
            rows[std::size_t(outsider)][std::size_t(lane)] = row;
        }
    }

    std::vector<u64> row_count(LANES.size()), checks(LANES.size());
    std::vector<u64> failures(LANES.size()), equalities(LANES.size());
    std::vector<i128> literal_min(LANES.size(), -1);
    std::vector<u64> literal_den(LANES.size(), 1), literal_mass(LANES.size());
    std::vector<int> literal_r(LANES.size());
    std::vector<u32> literal_mask(LANES.size());
    for (int outsider = 1; outsider < global_cutoff; ++outsider) {
        for (int lane = 0; lane < lanes; ++lane) {
            const LiteralRow& row = rows[std::size_t(outsider)][std::size_t(lane)];
            if (!row.active) continue;
            ++row_count[lane];
            checks[lane] += bodies[std::size_t(LANES[lane].core_count - 4)].size();
            failures[lane] += row.failures;
            equalities[lane] += row.equalities;
            if (literal_min[lane] < 0 ||
                row.min_margin * literal_den[lane] <
                    literal_min[lane] * row.denominator) {
                literal_min[lane] = row.min_margin;
                literal_den[lane] = row.denominator;
                literal_mass[lane] = row.min_mass;
                literal_r[lane] = outsider;
                literal_mask[lane] = row.min_mask;
            }
        }
    }

    for (int lane = 0; lane < lanes; ++lane) {
        const DirectResult direct_min = direct_body_sweep(
            body_speeds(results[lane].min_mask, LANES[lane], 0));
        require_same_fraction(results[lane].min_mass, mass.denominator,
                              direct_min.numerator, direct_min.denominator, "min replay");
        const DirectResult direct_cutoff = direct_body_sweep(
            body_speeds(results[lane].cutoff_mask, LANES[lane], 0));
        require_same_fraction(results[lane].cutoff_mass, mass.denominator,
                              direct_cutoff.numerator, direct_cutoff.denominator,
                              "cutoff replay");
        if (direct_cutoff.components != results[lane].cutoff_components) {
            throw std::runtime_error("component replay");
        }
        const DirectResult direct_literal = direct_body_sweep(
            body_speeds(literal_mask[lane], LANES[lane], literal_r[lane]));
        require_same_fraction(literal_mass[lane], literal_den[lane],
                              direct_literal.numerator, direct_literal.denominator,
                              "literal replay");

        std::cout << "LANE petals=" << LANES[lane].petal_count
                  << " labels=" << lane_labels(LANES[lane])
                  << " base_profiles="
                  << bodies[std::size_t(LANES[lane].core_count - 4)].size()
                  << " nonstrict=" << results[lane].nonstrict
                  << " min_delta=" << decimal(results[lane].min_delta)
                  << " min_body=" << core_body(results[lane].min_mask)
                  << " max_cutoff=" << decimal(results[lane].max_cutoff)
                  << " cutoff_body=" << core_body(results[lane].cutoff_mask)
                  << " cutoff_components=" << results[lane].cutoff_components
                  << " literal_rows=" << row_count[lane]
                  << " literal_checks=" << checks[lane]
                  << " failures=" << failures[lane]
                  << " equalities=" << equalities[lane]
                  << " literal_min_delta=" << decimal(literal_min[lane])
                  << " literal_min_denominator=" << literal_den[lane]
                  << " literal_min_r=" << literal_r[lane]
                  << " literal_min_body=" << core_body(literal_mask[lane])
                  << " direct_base_replays=2 direct_literal_replays=1\n";
    }
    const u64 total_rows = std::accumulate(row_count.begin(), row_count.end(), u64{0});
    const u64 total_checks = std::accumulate(checks.begin(), checks.end(), u64{0});
    const u64 total_failures = std::accumulate(failures.begin(), failures.end(), u64{0});
    const u64 total_equalities = std::accumulate(equalities.begin(), equalities.end(), u64{0});
    std::cout << "SUMMARY denominator=" << mass.denominator
              << " petals=10,15,20,30,264 lanes=16 base_profiles=231540"
              << " literal_rows=" << total_rows
              << " literal_checks=" << total_checks
              << " failures=" << total_failures
              << " equalities=" << total_equalities
              << " direct_base_replays=32 direct_literal_replays=16"
              << " status=FINITE-EXACT-INDEPENDENT LRC14=OPEN\n";
    return total_failures == 0 ? 0 : 3;
}
