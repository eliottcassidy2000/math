// Clean-room grouped-event audit for the fixed-50 quartet
// U={63,132,176,264}.  It simultaneously audits the three new triple
// faces containing 63 and the new four-petal lane.  No floating point.

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
constexpr int NO = 12;
constexpr int NL = 4;
constexpr std::size_t STATES = std::size_t{1} << NC;
constexpr u32 CORE_MASK = (u32{1} << NC) - 1;

constexpr std::array<int, NC> C = {
    8, 16, 40, 42, 80, 84, 85, 88, 95,
    120, 126, 143, 145, 168, 193, 240, 252, 286,
};
constexpr std::array<int, NO> O = {
    10, 15, 20, 30, 60, 63, 132, 170, 176, 190, 264, 290,
};
constexpr std::array<int, 30> P = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
};

struct Lane {
    std::array<int, 4> labels;
    int petal_count;
    int core_count;
};

constexpr std::array<Lane, NL> LANES = {{
    {{{63, 132, 176, 0}}, 3, 6},
    {{{63, 132, 264, 0}}, 3, 6},
    {{{63, 176, 264, 0}}, 3, 6},
    {{{63, 132, 176, 264}}, 4, 5},
}};

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
    if (value > std::numeric_limits<u64>::max()) {
        throw std::runtime_error("LCM overflow");
    }
    return u64(value);
}

bool excluded_outsider(int value) {
    if (value == 50) return true;
    return std::find(P.begin(), P.end(), value) != P.end();
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

std::string lane_labels(int lane) {
    std::string answer;
    for (int j = 0; j < LANES[lane].petal_count; ++j) {
        if (!answer.empty()) answer += ',';
        answer += std::to_string(LANES[lane].labels[j]);
    }
    return answer;
}

struct Event {
    u64 position;
    u32 toggle;
};

std::vector<Event> grouped_events(const std::vector<int>& labels,
                                  u64 denominator) {
    std::vector<Event> raw;
    u64 expected = 0;
    for (int speed : labels) expected += 2 * u64(speed);
    raw.reserve(std::size_t(expected));
    for (std::size_t bit = 0; bit < labels.size(); ++bit) {
        const int speed = labels[bit];
        const u64 local_denominator = 14 * u64(speed);
        if (denominator % local_denominator) {
            throw std::runtime_error("nonintegral wall grid");
        }
        const u64 scale = denominator / local_denominator;
        for (int tooth = 0; tooth < speed; ++tooth) {
            int left = 14 * tooth - 1;
            if (left < 0) left += int(local_denominator);
            const int right = 14 * tooth + 1;
            raw.push_back({u64(left) * scale, u32{1} << bit});
            raw.push_back({u64(right) * scale, u32{1} << bit});
        }
    }
    std::sort(raw.begin(), raw.end(), [](const Event& first,
                                         const Event& second) {
        return first.position < second.position;
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
                                [](const Event& event) {
                                    return event.toggle == 0;
                                }),
                 answer.end());
    return answer;
}

template <class Integer>
void subset_zeta(std::vector<Integer>& table) {
    for (int bit = 0; bit < NC; ++bit) {
        const std::size_t step = std::size_t{1} << bit;
        const std::size_t width = 2 * step;
        for (std::size_t block = 0; block < STATES; block += width) {
            for (std::size_t offset = 0; offset < step; ++offset) {
                Integer* lower = table.data() + (block + offset) * NL;
                Integer* upper = table.data() + (block + step + offset) * NL;
                for (int lane = 0; lane < NL; ++lane) upper[lane] += lower[lane];
            }
        }
    }
}

struct EventGeometry {
    u64 denominator;
    std::vector<Event> events;
    u32 initial;
};

EventGeometry geometry(int outsider) {
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

std::array<u32, NL> fixed_masks(int outsider) {
    const u32 bit50 = u32{1} << 30;
    const u32 bitr = outsider > 0 ? (u32{1} << 31) : 0;
    std::array<u32, NL> answer{};
    for (int lane = 0; lane < NL; ++lane) {
        answer[lane] = bit50 | bitr;
        for (int j = 0; j < LANES[lane].petal_count; ++j) {
            const int label = LANES[lane].labels[j];
            const auto it = std::find(O.begin(), O.end(), label);
            if (it == O.end()) throw std::runtime_error("petal outside O");
            answer[lane] |= u32{1} << (NC + int(it - O.begin()));
        }
    }
    return answer;
}

struct MassTable {
    u64 denominator;
    std::vector<u64> values;
};

MassTable mass_table(int outsider) {
    const EventGeometry geo = geometry(outsider);
    const std::array<u32, NL> fixed = fixed_masks(outsider);
    std::vector<u64> values(STATES * NL, 0);
    u32 state = geo.initial;
    u64 previous = 0;
    auto add_interval = [&](u64 length) {
        if (length == 0) return;
        const std::size_t row = std::size_t(state & CORE_MASK) * NL;
        for (int lane = 0; lane < NL; ++lane) {
            if ((state & fixed[lane]) == 0) values[row + lane] += length;
        }
    };
    for (const Event& event : geo.events) {
        add_interval(event.position - previous);
        state ^= event.toggle;
        previous = event.position;
    }
    add_interval(geo.denominator - previous);
    if (state != geo.initial) throw std::runtime_error("mass circular mismatch");
    subset_zeta(values);
    return {geo.denominator, std::move(values)};
}

struct ComponentTable {
    u64 denominator;
    std::vector<i32> values;
};

ComponentTable component_table() {
    const EventGeometry geo = geometry(0);
    const std::array<u32, NL> fixed = fixed_masks(0);
    std::vector<i32> values(STATES * NL, 0);
    u32 state = geo.initial;
    for (const Event& event : geo.events) {
        const u32 before = state;
        state ^= event.toggle;
        const std::size_t start_row = std::size_t(state & CORE_MASK) * NL;
        const std::size_t bridge_row = std::size_t((before | state) & CORE_MASK) * NL;
        for (int lane = 0; lane < NL; ++lane) {
            const bool before_safe = (before & fixed[lane]) == 0;
            const bool after_safe = (state & fixed[lane]) == 0;
            if (after_safe) ++values[start_row + lane];
            if (before_safe && after_safe) --values[bridge_row + lane];
        }
    }
    if (state != geo.initial) throw std::runtime_error("component circular mismatch");
    subset_zeta(values);
    return {geo.denominator, std::move(values)};
}

struct DirectResult {
    u64 numerator;
    u64 denominator;
    u64 components;
};

DirectResult direct_body_sweep(std::vector<int> speeds) {
    std::sort(speeds.begin(), speeds.end());
    speeds.erase(std::unique(speeds.begin(), speeds.end()), speeds.end());
    u64 denominator = 1;
    for (int speed : speeds) denominator = checked_lcm(denominator, 14 * u64(speed));
    const std::vector<Event> events = grouped_events(speeds, denominator);
    const u32 initial = (u32{1} << speeds.size()) - 1;
    u32 state = initial;
    u64 previous = 0;
    u64 numerator = 0;
    u64 components = 0;
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

std::vector<int> body_speeds(u32 mask, int lane, int outsider) {
    std::vector<int> speeds;
    for (int bit = 0; bit < NC; ++bit) {
        if (mask & (u32{1} << bit)) speeds.push_back(C[bit]);
    }
    for (int j = 0; j < LANES[lane].petal_count; ++j) {
        speeds.push_back(LANES[lane].labels[j]);
    }
    speeds.push_back(50);
    if (outsider > 0) speeds.push_back(outsider);
    return speeds;
}

void require_same_fraction(u64 first_num, u64 first_den,
                           u64 second_num, u64 second_den,
                           const char* message) {
    if (i128(first_num) * second_den != i128(second_num) * first_den) {
        throw std::runtime_error(message);
    }
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
    const std::vector<u32> choose5 = masks_of_size(5);
    const std::vector<u32> choose6 = masks_of_size(6);
    if (choose5.size() != 8568 || choose6.size() != 18564) {
        throw std::runtime_error("body universes");
    }
    const DirectResult one_speed = direct_body_sweep({1});
    if (one_speed.numerator * 7 != one_speed.denominator * 6 ||
        one_speed.components != 1) {
        throw std::runtime_error("one-speed control");
    }

    const MassTable base_mass = mass_table(0);
    const ComponentTable base_components = component_table();
    if (base_mass.denominator != UINT64_C(91205797082400) ||
        base_components.denominator != base_mass.denominator) {
        throw std::runtime_error("base denominator");
    }

    std::array<BaseResult, NL> base{};
    for (int lane = 0; lane < NL; ++lane) {
        const std::vector<u32>& bodies = LANES[lane].core_count == 6 ? choose6 : choose5;
        for (u32 mask : bodies) {
            const std::size_t index = std::size_t(CORE_MASK ^ mask) * NL + lane;
            const u64 numerator = base_mass.values[index];
            const i128 delta = 54 * i128(numerator) - 4 * i128(base_mass.denominator);
            const i32 components = base_components.values[index];
            if (delta <= 0) ++base[lane].nonstrict;
            if (components <= 0) throw std::runtime_error("nonpositive components");
            if (base[lane].min_delta < 0 || delta < base[lane].min_delta) {
                base[lane].min_delta = delta;
                base[lane].min_mask = mask;
                base[lane].min_mass = numerator;
            }
            if (delta > 0) {
                const i128 cutoff = ceil_positive(
                    54 * i128(components) * base_mass.denominator,
                    7 * delta);
                if (cutoff > base[lane].max_cutoff) {
                    base[lane].max_cutoff = cutoff;
                    base[lane].cutoff_mask = mask;
                    base[lane].cutoff_mass = numerator;
                    base[lane].cutoff_components = u64(components);
                }
            }
        }
        if (base[lane].nonstrict != 0) throw std::runtime_error("nonstrict base lane");
    }

    int global_cutoff = 0;
    for (const BaseResult& result : base) {
        global_cutoff = std::max(global_cutoff, int(result.max_cutoff));
    }
    std::vector<std::array<LiteralRow, NL>> rows{std::size_t(global_cutoff)};
#pragma omp parallel for schedule(dynamic, 1)
    for (int outsider = 1; outsider < global_cutoff; ++outsider) {
        if (excluded_outsider(outsider)) continue;
        const MassTable table = mass_table(outsider);
        for (int lane = 0; lane < NL; ++lane) {
            if (outsider >= base[lane].max_cutoff) continue;
            LiteralRow result;
            result.active = true;
            result.denominator = table.denominator;
            const std::vector<u32>& bodies = LANES[lane].core_count == 6 ? choose6 : choose5;
            for (u32 mask : bodies) {
                const u64 numerator =
                    table.values[std::size_t(CORE_MASK ^ mask) * NL + lane];
                const i128 margin = 63 * i128(numerator) - 4 * i128(table.denominator);
                if (margin < 0) ++result.failures;
                if (margin == 0) ++result.equalities;
                if (result.min_margin < 0 || margin < result.min_margin) {
                    result.min_margin = margin;
                    result.min_mass = numerator;
                    result.min_mask = mask;
                }
            }
            rows[std::size_t(outsider)][lane] = result;
        }
    }

    std::array<u64, NL> row_count{}, checks{}, failures{}, equalities{};
    std::array<i128, NL> literal_min{};
    std::array<u64, NL> literal_den{}, literal_mass{};
    std::array<int, NL> literal_r{};
    std::array<u32, NL> literal_mask{};
    literal_min.fill(-1);
    literal_den.fill(1);
    for (int outsider = 1; outsider < global_cutoff; ++outsider) {
        for (int lane = 0; lane < NL; ++lane) {
            const LiteralRow& row = rows[std::size_t(outsider)][lane];
            if (!row.active) continue;
            ++row_count[lane];
            checks[lane] += LANES[lane].core_count == 6 ? choose6.size() : choose5.size();
            failures[lane] += row.failures;
            equalities[lane] += row.equalities;
            if (literal_min[lane] < 0 ||
                row.min_margin * literal_den[lane] < literal_min[lane] * row.denominator) {
                literal_min[lane] = row.min_margin;
                literal_den[lane] = row.denominator;
                literal_mass[lane] = row.min_mass;
                literal_r[lane] = outsider;
                literal_mask[lane] = row.min_mask;
            }
        }
    }

    for (int lane = 0; lane < NL; ++lane) {
        const DirectResult direct_min = direct_body_sweep(body_speeds(base[lane].min_mask, lane, 0));
        require_same_fraction(base[lane].min_mass, base_mass.denominator,
                              direct_min.numerator, direct_min.denominator,
                              "base minimum direct replay");
        const DirectResult direct_cutoff =
            direct_body_sweep(body_speeds(base[lane].cutoff_mask, lane, 0));
        require_same_fraction(base[lane].cutoff_mass, base_mass.denominator,
                              direct_cutoff.numerator, direct_cutoff.denominator,
                              "base cutoff direct replay");
        if (direct_cutoff.components != base[lane].cutoff_components) {
            throw std::runtime_error("cutoff component direct replay");
        }
        const DirectResult direct_literal =
            direct_body_sweep(body_speeds(literal_mask[lane], lane, literal_r[lane]));
        require_same_fraction(literal_mass[lane], literal_den[lane],
                              direct_literal.numerator, direct_literal.denominator,
                              "literal minimum direct replay");

        std::cout << (LANES[lane].petal_count == 3 ? "TRIPLE" : "QUADRUPLE")
                  << " labels=" << lane_labels(lane)
                  << " base_profiles="
                  << (LANES[lane].core_count == 6 ? choose6.size() : choose5.size())
                  << " nonstrict=" << base[lane].nonstrict
                  << " min_delta=" << decimal(base[lane].min_delta)
                  << " min_body=" << core_body(base[lane].min_mask)
                  << " max_cutoff=" << decimal(base[lane].max_cutoff)
                  << " cutoff_body=" << core_body(base[lane].cutoff_mask)
                  << " cutoff_components=" << base[lane].cutoff_components
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
    std::cout << "SUMMARY denominator=" << base_mass.denominator
              << " lanes=4 failures="
              << std::accumulate(failures.begin(), failures.end(), u64{0})
              << " equalities="
              << std::accumulate(equalities.begin(), equalities.end(), u64{0})
              << " status=FINITE-EXACT-INDEPENDENT LRC14=OPEN\n";
    return std::accumulate(failures.begin(), failures.end(), u64{0}) == 0 ? 0 : 3;
}
