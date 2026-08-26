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
using i64 = std::int64_t;
using i32 = std::int32_t;
using i128 = __int128_t;

namespace {

constexpr int NC = 18;
constexpr int NO = 12;
constexpr int NTRIPLES = 220;
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

struct TripleIndex {
    int first;
    int second;
    int third;
};

std::array<TripleIndex, NTRIPLES> make_triples() {
    std::array<TripleIndex, NTRIPLES> answer{};
    int lane = 0;
    for (int first = 0; first < NO; ++first) {
        for (int second = first + 1; second < NO; ++second) {
            for (int third = second + 1; third < NO; ++third) {
                answer[lane++] = {first, second, third};
            }
        }
    }
    if (lane != NTRIPLES) throw std::runtime_error("triple universe");
    return answer;
}

const std::array<TripleIndex, NTRIPLES> TRIPLES = make_triples();

std::string decimal(i128 value) {
    if (value == 0) return "0";
    bool negative = value < 0;
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
void subset_zeta(std::vector<Integer>& table, int lanes) {
    for (int bit = 0; bit < NC; ++bit) {
        const std::size_t step = std::size_t{1} << bit;
        const std::size_t width = 2 * step;
        for (std::size_t block = 0; block < STATES; block += width) {
            for (std::size_t offset = 0; offset < step; ++offset) {
                Integer* lower = table.data() + (block + offset) * lanes;
                Integer* upper = table.data() + (block + step + offset) * lanes;
                for (int lane = 0; lane < lanes; ++lane) {
                    upper[lane] += lower[lane];
                }
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

std::vector<u32> fixed_masks(const std::vector<int>& active, int outsider) {
    const u32 bit50 = u32{1} << 30;
    const u32 bitr = outsider > 0 ? (u32{1} << 31) : 0;
    std::vector<u32> answer;
    answer.reserve(active.size());
    for (int lane : active) {
        const TripleIndex triple = TRIPLES[lane];
        answer.push_back(bit50 | bitr |
                         (u32{1} << (NC + triple.first)) |
                         (u32{1} << (NC + triple.second)) |
                         (u32{1} << (NC + triple.third)));
    }
    return answer;
}

struct MassTable {
    u64 denominator;
    int lanes;
    std::vector<u64> values;
};

MassTable mass_table(int outsider, const std::vector<int>& active) {
    const EventGeometry geo = geometry(outsider);
    const int lanes = int(active.size());
    const std::vector<u32> fixed = fixed_masks(active, outsider);
    std::vector<u64> values(STATES * std::size_t(lanes), 0);
    u32 state = geo.initial;
    u64 previous = 0;
    auto add_interval = [&](u64 length) {
        if (length == 0) return;
        const std::size_t row = std::size_t(state & CORE_MASK) * lanes;
        for (int lane = 0; lane < lanes; ++lane) {
            if ((state & fixed[lane]) == 0) {
                values[row + std::size_t(lane)] += length;
            }
        }
    };
    for (const Event& event : geo.events) {
        add_interval(event.position - previous);
        state ^= event.toggle;
        previous = event.position;
    }
    add_interval(geo.denominator - previous);
    if (state != geo.initial) throw std::runtime_error("mass circular mismatch");
    subset_zeta(values, lanes);
    return {geo.denominator, lanes, std::move(values)};
}

struct ComponentTable {
    u64 denominator;
    int lanes;
    std::vector<i32> values;
};

ComponentTable component_table(const std::vector<int>& active) {
    const EventGeometry geo = geometry(0);
    const int lanes = int(active.size());
    const std::vector<u32> fixed = fixed_masks(active, 0);
    std::vector<i32> values(STATES * std::size_t(lanes), 0);
    u32 state = geo.initial;
    for (const Event& event : geo.events) {
        const u32 before = state;
        state ^= event.toggle;
        const std::size_t start_row = std::size_t(state & CORE_MASK) * lanes;
        const std::size_t bridge_row =
            std::size_t((before | state) & CORE_MASK) * lanes;
        for (int lane = 0; lane < lanes; ++lane) {
            const bool before_safe = (before & fixed[lane]) == 0;
            const bool after_safe = (state & fixed[lane]) == 0;
            if (after_safe) ++values[start_row + std::size_t(lane)];
            if (before_safe && after_safe) {
                --values[bridge_row + std::size_t(lane)];
            }
        }
    }
    if (state != geo.initial) throw std::runtime_error("component circular mismatch");
    subset_zeta(values, lanes);
    return {geo.denominator, lanes, std::move(values)};
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

std::vector<int> body_speeds(u32 mask, int triple_lane, int outsider) {
    std::vector<int> speeds;
    for (int bit = 0; bit < NC; ++bit) {
        if (mask & (u32{1} << bit)) speeds.push_back(C[bit]);
    }
    const TripleIndex triple = TRIPLES[triple_lane];
    speeds.push_back(O[triple.first]);
    speeds.push_back(O[triple.second]);
    speeds.push_back(O[triple.third]);
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
    const std::vector<u32> choose6 = masks_of_size(6);
    if (choose6.size() != 18564) throw std::runtime_error("choose(18,6)");
    const DirectResult one_speed = direct_body_sweep({1});
    if (one_speed.numerator * 7 != one_speed.denominator * 6 ||
        one_speed.components != 1) {
        throw std::runtime_error("one-speed control");
    }

    std::vector<int> all_triples(NTRIPLES);
    std::iota(all_triples.begin(), all_triples.end(), 0);
    std::array<BaseResult, NTRIPLES> base_results;
    std::vector<i64> deltas(choose6.size() * NTRIPLES);

    {
        const MassTable mass = mass_table(0, all_triples);
        if (mass.denominator != UINT64_C(91205797082400) ||
            mass.lanes != NTRIPLES) {
            throw std::runtime_error("triple base mass table");
        }
        for (std::size_t body_index = 0; body_index < choose6.size(); ++body_index) {
            const u32 mask = choose6[body_index];
            const std::size_t row = std::size_t(CORE_MASK ^ mask) * NTRIPLES;
            for (int lane = 0; lane < NTRIPLES; ++lane) {
                BaseResult& result = base_results[lane];
                const u64 numerator = mass.values[row + std::size_t(lane)];
                const i128 wide_delta =
                    54 * i128(numerator) - 4 * i128(mass.denominator);
                if (wide_delta < std::numeric_limits<i64>::min() ||
                    wide_delta > std::numeric_limits<i64>::max()) {
                    throw std::runtime_error("delta storage range");
                }
                const i64 delta = i64(wide_delta);
                deltas[body_index * NTRIPLES + std::size_t(lane)] = delta;
                if (delta <= 0) ++result.nonstrict;
                if (result.min_delta < 0 || delta < result.min_delta) {
                    result.min_delta = delta;
                    result.min_mask = mask;
                    result.min_mass = numerator;
                }
            }
        }
    }

    {
        const ComponentTable components = component_table(all_triples);
        if (components.denominator != UINT64_C(91205797082400) ||
            components.lanes != NTRIPLES) {
            throw std::runtime_error("triple component table");
        }
        for (std::size_t body_index = 0; body_index < choose6.size(); ++body_index) {
            const u32 mask = choose6[body_index];
            const std::size_t row = std::size_t(CORE_MASK ^ mask) * NTRIPLES;
            for (int lane = 0; lane < NTRIPLES; ++lane) {
                BaseResult& result = base_results[lane];
                const i64 delta = deltas[body_index * NTRIPLES + std::size_t(lane)];
                const i32 component_count = components.values[row + std::size_t(lane)];
                if (delta <= 0 || component_count <= 0) {
                    throw std::runtime_error("triple strict/component gate");
                }
                const i128 cutoff = ceil_positive(
                    54 * i128(component_count) * components.denominator,
                    7 * i128(delta));
                if (cutoff > result.max_cutoff) {
                    result.max_cutoff = cutoff;
                    result.cutoff_mask = mask;
                    result.cutoff_components = u64(component_count);
                    const i128 mass_numerator =
                        (i128(delta) + 4 * i128(components.denominator)) / 54;
                    result.cutoff_mass = u64(mass_numerator);
                }
            }
        }
    }
    deltas.clear();
    deltas.shrink_to_fit();

    int best_lane = -1;
    for (int lane = 0; lane < NTRIPLES; ++lane) {
        if (base_results[lane].nonstrict != 0 || base_results[lane].min_delta <= 0) {
            throw std::runtime_error("all-triple limiting strictness");
        }
        if (best_lane < 0 ||
            base_results[lane].min_delta > base_results[best_lane].min_delta) {
            best_lane = lane;
        }
    }
    const TripleIndex best = TRIPLES[best_lane];
    if (O[best.first] != 132 || O[best.second] != 176 || O[best.third] != 264 ||
        base_results[best_lane].min_delta != i128(UINT64_C(687816435418308)) ||
        base_results[best_lane].max_cutoff != 470) {
        throw std::runtime_error("selected strongest-minimum triple");
    }

    const int selected_cutoff = int(base_results[best_lane].max_cutoff);
    std::vector<LiteralRow> rows{std::size_t(selected_cutoff)};
#pragma omp parallel for schedule(dynamic, 1)
    for (int outsider = 1; outsider < selected_cutoff; ++outsider) {
        if (excluded_outsider(outsider)) continue;
        const MassTable table = mass_table(outsider, std::vector<int>{best_lane});
        LiteralRow result;
        result.active = true;
        result.denominator = table.denominator;
        for (u32 mask : choose6) {
            const u64 numerator = table.values[std::size_t(CORE_MASK ^ mask)];
            const i128 margin = 63 * i128(numerator) - 4 * i128(table.denominator);
            if (margin < 0) ++result.failures;
            if (margin == 0) ++result.equalities;
            if (result.min_margin < 0 || margin < result.min_margin) {
                result.min_margin = margin;
                result.min_mass = numerator;
                result.min_mask = mask;
            }
        }
        rows[std::size_t(outsider)] = result;
    }

    u64 literal_row_count = 0;
    u64 literal_checks = 0;
    u64 literal_failures = 0;
    u64 literal_equalities = 0;
    i128 literal_min_margin = -1;
    u64 literal_min_denominator = 1;
    u64 literal_min_mass = 0;
    int literal_min_r = 0;
    u32 literal_min_mask = 0;
    for (int outsider = 1; outsider < selected_cutoff; ++outsider) {
        const LiteralRow& row = rows[std::size_t(outsider)];
        if (!row.active) continue;
        ++literal_row_count;
        literal_checks += choose6.size();
        literal_failures += row.failures;
        literal_equalities += row.equalities;
        if (literal_min_margin < 0 ||
            row.min_margin * literal_min_denominator <
                literal_min_margin * row.denominator) {
            literal_min_margin = row.min_margin;
            literal_min_denominator = row.denominator;
            literal_min_mass = row.min_mass;
            literal_min_r = outsider;
            literal_min_mask = row.min_mask;
        }
    }
    if (literal_row_count != 438 || literal_checks != UINT64_C(8131032) ||
        literal_failures != 0 || literal_equalities != 0) {
        throw std::runtime_error("selected literal totals");
    }

    // Direct body-local controls for both limiting extremizers of all 220
    // triples, plus the selected triple's closest literal body.
    for (int lane = 0; lane < NTRIPLES; ++lane) {
        const BaseResult& result = base_results[lane];
        const DirectResult direct_min =
            direct_body_sweep(body_speeds(result.min_mask, lane, 0));
        require_same_fraction(result.min_mass, UINT64_C(91205797082400),
                              direct_min.numerator, direct_min.denominator,
                              "triple min direct replay");
        const DirectResult direct_cutoff =
            direct_body_sweep(body_speeds(result.cutoff_mask, lane, 0));
        require_same_fraction(result.cutoff_mass, UINT64_C(91205797082400),
                              direct_cutoff.numerator, direct_cutoff.denominator,
                              "triple cutoff mass direct replay");
        if (direct_cutoff.components != result.cutoff_components) {
            throw std::runtime_error("triple cutoff component replay");
        }
    }
    const DirectResult direct_literal =
        direct_body_sweep(body_speeds(literal_min_mask, best_lane, literal_min_r));
    require_same_fraction(literal_min_mass, literal_min_denominator,
                          direct_literal.numerator, direct_literal.denominator,
                          "selected triple literal replay");

    for (int lane = 0; lane < NTRIPLES; ++lane) {
        const TripleIndex triple = TRIPLES[lane];
        const BaseResult& result = base_results[lane];
        std::cout << "TRIPLE p=" << O[triple.first]
                  << " q=" << O[triple.second]
                  << " u=" << O[triple.third]
                  << " nonstrict=" << result.nonstrict
                  << " min_delta=" << decimal(result.min_delta)
                  << " min_body=" << core_body(result.min_mask)
                  << " max_cutoff=" << decimal(result.max_cutoff)
                  << " cutoff_body=" << core_body(result.cutoff_mask)
                  << '\n';
    }
    std::cout << "TRIPLE_SUMMARY denominator=91205797082400"
              << " triples=" << NTRIPLES
              << " profiles=" << u64(NTRIPLES) * choose6.size()
              << " strict=" << NTRIPLES
              << " nonstrict=0"
              << " best=" << O[best.first] << ',' << O[best.second] << ','
              << O[best.third]
              << " best_min_delta=" << decimal(base_results[best_lane].min_delta)
              << " best_cutoff=" << decimal(base_results[best_lane].max_cutoff)
              << '\n';
    std::cout << "SELECTED_TRIPLE labels=" << O[best.first] << ','
              << O[best.second] << ',' << O[best.third]
              << " limiting_profiles=" << choose6.size()
              << " limiting_min_delta=" << decimal(base_results[best_lane].min_delta)
              << " cutoff=" << selected_cutoff
              << " literal_rows=" << literal_row_count
              << " literal_checks=" << literal_checks
              << " failures=" << literal_failures
              << " equalities=" << literal_equalities
              << " min_delta=" << decimal(literal_min_margin)
              << " min_denominator=" << literal_min_denominator
              << " min_r=" << literal_min_r
              << " min_body=" << core_body(literal_min_mask)
              << " direct_base_replays=" << 2 * NTRIPLES
              << " direct_literal_replays=1"
              << " status=FINITE-EXACT-INDEPENDENT"
              << " LRC14=OPEN\n";
    return 0;
}
