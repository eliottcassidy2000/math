#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
__extension__ typedef unsigned __int128 u128;

constexpr std::array<int, 18> C = {
    8, 16, 40, 42, 80, 84, 85, 88, 95,
    120, 126, 143, 145, 168, 193, 240, 252, 286,
};
constexpr std::array<int, 3> TRIPLE = {132, 176, 264};
constexpr std::array<int, 30> P = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
};
constexpr int CENTER = 50;
constexpr std::uint32_t FULL = (std::uint32_t{1} << C.size()) - 1;
constexpr std::size_t STATES = std::size_t{1} << C.size();
constexpr u64 EXPECTED_D = 91'205'797'082'400ULL;

struct Event {
    u64 position;
    std::uint32_t c_toggle;
    std::uint32_t mandatory_toggle;
};

struct Scatter {
    u64 grid;
    std::vector<u64> mass;
    std::vector<u64> starts;
    std::vector<u64> bridges;
    std::size_t raw_events;
    std::size_t grouped_events;
};

struct DirectProfile {
    u64 grid;
    u64 mass;
    u64 components;
};

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

u64 checked_lcm(u64 first, u64 second) {
    const u64 g = std::gcd(first, second);
    const u128 value = u128(first / g) * second;
    assert(value <= std::numeric_limits<u64>::max());
    return static_cast<u64>(value);
}

u64 fixed_grid() {
    u64 lcm = 1;
    for (int speed : C) lcm = checked_lcm(lcm, static_cast<u64>(speed));
    for (int speed : TRIPLE) lcm = checked_lcm(lcm, static_cast<u64>(speed));
    lcm = checked_lcm(lcm, CENTER);
    const u64 answer = 14 * lcm;
    assert(answer == EXPECTED_D);
    return answer;
}

bool in_pool(int value) {
    return std::find(P.begin(), P.end(), value) != P.end();
}

void append_speed_events(std::vector<Event>& events, int speed, u64 grid,
                         int c_index, int mandatory_index) {
    const u64 period_denominator = 14ULL * static_cast<u64>(speed);
    assert(grid % period_denominator == 0);
    const u64 scale = grid / period_denominator;
    const std::uint32_t c_bit = c_index >= 0
        ? (std::uint32_t{1} << c_index) : 0;
    const std::uint32_t mandatory_bit = mandatory_index >= 0
        ? (std::uint32_t{1} << mandatory_index) : 0;
    for (int k = 0; k < speed; ++k) {
        const u64 center = 14ULL * static_cast<u64>(k);
        const u64 left_numerator = (center + period_denominator - 1)
            % period_denominator;
        const u64 right_numerator = (center + 1) % period_denominator;
        events.push_back({left_numerator * scale, c_bit, mandatory_bit});
        events.push_back({right_numerator * scale, c_bit, mandatory_bit});
    }
}

Scatter build_scatter(const std::vector<int>& mandatory_speeds, u64 grid,
                      bool need_components) {
    std::vector<Event> raw;
    std::size_t reserve = 0;
    for (int speed : C) reserve += static_cast<std::size_t>(2 * speed);
    for (int speed : mandatory_speeds) reserve += static_cast<std::size_t>(2 * speed);
    raw.reserve(reserve);
    for (std::size_t index = 0; index < C.size(); ++index) {
        append_speed_events(raw, C[index], grid, static_cast<int>(index), -1);
    }
    for (std::size_t index = 0; index < mandatory_speeds.size(); ++index) {
        append_speed_events(raw, mandatory_speeds[index], grid, -1,
                            static_cast<int>(index));
    }
    std::sort(raw.begin(), raw.end(), [](const Event& first, const Event& second) {
        return first.position < second.position;
    });

    Scatter answer{
        grid,
        std::vector<u64>(STATES, 0),
        need_components ? std::vector<u64>(STATES, 0) : std::vector<u64>(),
        need_components ? std::vector<u64>(STATES, 0) : std::vector<u64>(),
        raw.size(),
        0,
    };
    std::uint32_t c_mask = FULL;
    std::uint32_t mandatory_mask = mandatory_speeds.empty()
        ? 0 : (std::uint32_t{1} << mandatory_speeds.size()) - 1;
    const std::uint32_t initial_c_mask = c_mask;
    const std::uint32_t initial_mandatory_mask = mandatory_mask;
    u64 previous = 0;
    std::size_t cursor = 0;
    while (cursor < raw.size()) {
        const u64 position = raw[cursor].position;
        assert(position > previous || answer.grouped_events == 0);
        if (mandatory_mask == 0) answer.mass[c_mask] += position - previous;

        const std::uint32_t left_c_mask = c_mask;
        const std::uint32_t left_mandatory_mask = mandatory_mask;
        std::uint32_t c_toggle = 0;
        std::uint32_t mandatory_toggle = 0;
        while (cursor < raw.size() && raw[cursor].position == position) {
            c_toggle ^= raw[cursor].c_toggle;
            mandatory_toggle ^= raw[cursor].mandatory_toggle;
            ++cursor;
        }
        c_mask ^= c_toggle;
        mandatory_mask ^= mandatory_toggle;
        if (need_components && mandatory_mask == 0) {
            ++answer.starts[c_mask];
            if (left_mandatory_mask == 0) {
                ++answer.bridges[left_c_mask | c_mask];
            }
        }
        previous = position;
        ++answer.grouped_events;
    }
    if (mandatory_mask == 0) answer.mass[c_mask] += grid - previous;
    assert(c_mask == initial_c_mask);
    assert(mandatory_mask == initial_mandatory_mask);
    assert(answer.raw_events == reserve);
    return answer;
}

DirectProfile direct_profile(const std::vector<int>& speeds, u64 grid) {
    std::vector<std::pair<u64, std::size_t>> events;
    std::size_t reserve = 0;
    for (int speed : speeds) reserve += static_cast<std::size_t>(2 * speed);
    events.reserve(reserve);
    for (std::size_t index = 0; index < speeds.size(); ++index) {
        const u64 denominator = 14ULL * static_cast<u64>(speeds[index]);
        assert(grid % denominator == 0);
        const u64 scale = grid / denominator;
        for (int k = 0; k < speeds[index]; ++k) {
            const u64 center = 14ULL * static_cast<u64>(k);
            const u64 left = (center + denominator - 1) % denominator;
            const u64 right = (center + 1) % denominator;
            events.emplace_back(left * scale, index);
            events.emplace_back(right * scale, index);
        }
    }
    std::sort(events.begin(), events.end());
    std::vector<unsigned char> failing(speeds.size(), 1);
    u64 active = speeds.size();
    u64 mass = 0;
    u64 components = 0;
    u64 previous = 0;
    std::size_t cursor = 0;
    while (cursor < events.size()) {
        const u64 position = events[cursor].first;
        if (active == 0) mass += position - previous;
        const bool left_safe = active == 0;
        while (cursor < events.size() && events[cursor].first == position) {
            const std::size_t index = events[cursor].second;
            if (failing[index] != 0) {
                failing[index] = 0;
                --active;
            } else {
                failing[index] = 1;
                ++active;
            }
            ++cursor;
        }
        const bool right_safe = active == 0;
        if (right_safe && !left_safe) ++components;
        previous = position;
    }
    if (active == 0) mass += grid - previous;
    assert(active == speeds.size());
    assert(events.size() == reserve);
    return {grid, mass, components};
}

void subset_zeta(std::vector<u64>& values) {
    assert(values.size() == STATES);
    for (std::size_t step = 1; step < STATES; step <<= 1) {
        for (std::size_t block = 0; block < STATES; block += 2 * step) {
            for (std::size_t offset = 0; offset < step; ++offset) {
                values[block + step + offset] += values[block + offset];
            }
        }
    }
}

void choose_masks_rec(int next, int remaining, std::uint32_t mask,
                      std::vector<std::uint32_t>& answer) {
    if (remaining == 0) {
        answer.push_back(mask);
        return;
    }
    for (int index = next; index <= static_cast<int>(C.size()) - remaining; ++index) {
        choose_masks_rec(index + 1, remaining - 1,
                         mask | (std::uint32_t{1} << index), answer);
    }
}

std::vector<std::uint32_t> choose_six_masks() {
    std::vector<std::uint32_t> answer;
    answer.reserve(18'564);
    choose_masks_rec(0, 6, 0, answer);
    assert(answer.size() == 18'564);
    return answer;
}

std::string body_string(std::uint32_t mask) {
    std::string answer = "{";
    bool first = true;
    for (std::size_t index = 0; index < C.size(); ++index) {
        if ((mask & (std::uint32_t{1} << index)) == 0) continue;
        if (!first) answer += ',';
        first = false;
        answer += std::to_string(C[index]);
    }
    answer += '}';
    return answer;
}

std::vector<int> body_speeds(std::uint32_t mask, bool include_triple_center,
                             int outsider = 0) {
    std::vector<int> answer;
    for (std::size_t index = 0; index < C.size(); ++index) {
        if (mask & (std::uint32_t{1} << index)) answer.push_back(C[index]);
    }
    if (include_triple_center) {
        answer.insert(answer.end(), TRIPLE.begin(), TRIPLE.end());
        answer.push_back(CENTER);
    }
    if (outsider != 0) answer.push_back(outsider);
    return answer;
}

u64 ceil_ratio(u128 numerator, u128 denominator) {
    assert(denominator != 0);
    const u128 answer = (numerator + denominator - 1) / denominator;
    assert(answer <= std::numeric_limits<u64>::max());
    return static_cast<u64>(answer);
}

void fnv_word(u64& hash, u64 value) {
    constexpr u64 prime = 1'099'511'628'211ULL;
    for (int byte = 0; byte < 8; ++byte) {
        hash ^= (value >> (8 * byte)) & 0xffU;
        hash *= prime;
    }
}

}  // namespace

int main() {
    const u64 D = fixed_grid();
    const std::vector<std::uint32_t> bodies = choose_six_masks();

    std::vector<int> base_mandatory(TRIPLE.begin(), TRIPLE.end());
    base_mandatory.push_back(CENTER);
    Scatter base = build_scatter(base_mandatory, D, true);
    subset_zeta(base.mass);
    subset_zeta(base.starts);
    subset_zeta(base.bridges);

    u64 min_delta = std::numeric_limits<u64>::max();
    std::uint32_t min_delta_body = 0;
    u64 max_cutoff = 0;
    std::uint32_t max_cutoff_body = 0;
    u64 min_mass = std::numeric_limits<u64>::max();
    u64 max_components = 0;
    u64 profile_hash = 1'469'598'103'934'665'603ULL;

    for (std::uint32_t body : bodies) {
        const std::uint32_t complement = FULL ^ body;
        const u64 mass = base.mass[complement];
        const u64 components = base.starts[complement] - base.bridges[complement];
        const u128 delta128 = u128(54) * mass - u128(4) * D;
        assert(delta128 <= std::numeric_limits<u64>::max());
        const u64 delta = static_cast<u64>(delta128);
        assert(delta > 0);
        const u64 cutoff = ceil_ratio(u128(54) * components * D,
                                      u128(7) * delta);
        if (delta < min_delta) {
            min_delta = delta;
            min_delta_body = body;
        }
        if (cutoff > max_cutoff) {
            max_cutoff = cutoff;
            max_cutoff_body = body;
        }
        min_mass = std::min(min_mass, mass);
        max_components = std::max(max_components, components);
        fnv_word(profile_hash, body);
        fnv_word(profile_hash, mass);
        fnv_word(profile_hash, components);
        fnv_word(profile_hash, delta);
        fnv_word(profile_hash, cutoff);
    }

    assert(min_delta == 687'816'435'418'308ULL);
    assert(max_cutoff == 470);

    const DirectProfile direct_one = direct_profile({1}, 14);
    assert(direct_one.mass == 12 && direct_one.components == 1);
    const DirectProfile direct_min = direct_profile(
        body_speeds(min_delta_body, true), D);
    assert(direct_min.mass == base.mass[FULL ^ min_delta_body]);
    assert(direct_min.components
           == base.starts[FULL ^ min_delta_body] - base.bridges[FULL ^ min_delta_body]);
    const DirectProfile direct_cutoff = direct_profile(
        body_speeds(max_cutoff_body, true), D);
    assert(direct_cutoff.mass == base.mass[FULL ^ max_cutoff_body]);
    assert(direct_cutoff.components
           == base.starts[FULL ^ max_cutoff_body] - base.bridges[FULL ^ max_cutoff_body]);

    std::vector<int> literal_rows;
    for (int r = 1; r < static_cast<int>(max_cutoff); ++r) {
        if (r != CENTER && !in_pool(r)) literal_rows.push_back(r);
    }
    assert(literal_rows.size() == 438);

    u64 checks = 0;
    u64 failures = 0;
    u64 equalities = 0;
    bool have_closest = false;
    u128 closest_margin_numerator = 0;
    u64 closest_grid = 1;
    int closest_r = 0;
    std::uint32_t closest_body = 0;
    u64 literal_hash = 1'469'598'103'934'665'603ULL;

    for (int r : literal_rows) {
        const u64 grid = checked_lcm(D, 14ULL * static_cast<u64>(r));
        std::vector<int> mandatory(TRIPLE.begin(), TRIPLE.end());
        mandatory.push_back(CENTER);
        mandatory.push_back(r);
        Scatter current = build_scatter(mandatory, grid, false);
        subset_zeta(current.mass);
        for (std::uint32_t body : bodies) {
            const u64 safe = current.mass[FULL ^ body];
            const u128 lhs = u128(63) * safe;
            const u128 rhs = u128(4) * grid;
            ++checks;
            if (lhs < rhs) {
                ++failures;
            } else if (lhs == rhs) {
                ++equalities;
            } else {
                const u128 margin = lhs - rhs;
                if (!have_closest
                    || margin * closest_grid < closest_margin_numerator * grid) {
                    have_closest = true;
                    closest_margin_numerator = margin;
                    closest_grid = grid;
                    closest_r = r;
                    closest_body = body;
                }
            }
            fnv_word(literal_hash, static_cast<u64>(r));
            fnv_word(literal_hash, body);
            fnv_word(literal_hash, safe);
            fnv_word(literal_hash, grid);
        }
    }

    assert(checks == 8'131'032ULL);
    assert(failures == 0);
    assert(equalities == 0);
    const std::array<int, 6> expected_closest = {8, 95, 120, 145, 168, 286};
    std::uint32_t expected_mask = 0;
    for (int speed : expected_closest) {
        const auto iterator = std::find(C.begin(), C.end(), speed);
        assert(iterator != C.end());
        expected_mask |= std::uint32_t{1} << std::distance(C.begin(), iterator);
    }
    assert(closest_r == 96);
    assert(closest_body == expected_mask);

    const DirectProfile direct_closest = direct_profile(
        body_speeds(closest_body, true, closest_r), closest_grid);
    assert(u128(63) * direct_closest.mass - u128(4) * closest_grid
           == closest_margin_numerator);

    const u128 closest_denominator = u128(63) * closest_grid;
    u64 reduced_numerator = static_cast<u64>(closest_margin_numerator);
    u64 reduced_denominator = static_cast<u64>(closest_denominator);
    const u64 common = std::gcd(reduced_numerator, reduced_denominator);
    reduced_numerator /= common;
    reduced_denominator /= common;

    std::cout << "STATUS=FINITE-EXACT; SELECTED_TRIPLE={132,176,264}; CENTER=50\n";
    std::cout << "METHOD=18BIT_EVENT_COMPLEMENT_SCATTER; PRIMARY_SOURCE_NOT_IMPORTED\n";
    std::cout << "D=" << D << "; PROFILES=" << bodies.size()
              << "; BASE_RAW_EVENTS=" << base.raw_events
              << "; BASE_GROUPED_EVENTS=" << base.grouped_events << '\n';
    std::cout << "MIN_DELTA=" << min_delta
              << "; BODY=" << body_string(min_delta_body) << '\n';
    std::cout << "MAX_CUTOFF=" << max_cutoff
              << "; BODY=" << body_string(max_cutoff_body)
              << "; MAX_COMPONENTS=" << max_components
              << "; MIN_MASS_NUMERATOR=" << min_mass << '\n';
    std::cout << "LITERAL_R_ROWS=" << literal_rows.size()
              << "; CHECKS=" << checks
              << "; FAILURES=" << failures
              << "; EQUALITIES=" << equalities << '\n';
    std::cout << "CLOSEST_R=" << closest_r
              << "; BODY=" << body_string(closest_body)
              << "; MU_MINUS_4_OVER_63=" << reduced_numerator
              << '/' << reduced_denominator
              << "; RAW_MARGIN_OVER_GRID="
              << decimal(closest_margin_numerator) << '/' << closest_grid << '\n';
    std::cout << "PROFILE_FNV64=" << profile_hash
              << "; LITERAL_FNV64=" << literal_hash << '\n';
    std::cout << "DIRECT_BODY_LOCAL_REPLAYS=4; ONE_SPEED+MIN_DELTA+MAX_CUTOFF+CLOSEST\n";
    std::cout << "CONSEQUENCE=CHI_50>=21_RELATIVE_TO_THM4229_AND_ALL_PAIR_LAYER\n";
    std::cout << "PASS\n";
}
