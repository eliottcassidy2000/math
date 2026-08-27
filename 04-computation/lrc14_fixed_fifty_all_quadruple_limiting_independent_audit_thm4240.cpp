// Clean-room exact limiting audit of every fixed-50 four-petal chart.
// Universe: T in C(O,4), L in C(C,5), V=G_{L union T union {50}}.
// Grouped simultaneous wall events + subset-zeta transforms; no floating point.

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
constexpr int NQUADS = 495;
constexpr int BATCH = 45;
constexpr std::size_t STATES = std::size_t{1} << NC;
constexpr u32 CORE_MASK = (u32{1} << NC) - 1;
constexpr std::array<int, NC> C = {
    8, 16, 40, 42, 80, 84, 85, 88, 95,
    120, 126, 143, 145, 168, 193, 240, 252, 286,
};
constexpr std::array<int, NO> O = {
    10, 15, 20, 30, 60, 63, 132, 170, 176, 190, 264, 290,
};

struct Quad { int a; int b; int c; int d; };

std::array<Quad, NQUADS> make_quads() {
    std::array<Quad, NQUADS> answer{};
    int lane = 0;
    for (int a = 0; a < NO; ++a) for (int b = a + 1; b < NO; ++b)
        for (int c = b + 1; c < NO; ++c) for (int d = c + 1; d < NO; ++d) {
            answer[std::size_t(lane++)] = {a, b, c, d};
        }
    if (lane != NQUADS) throw std::runtime_error("quadruple universe");
    return answer;
}

const std::array<Quad, NQUADS> QUADS = make_quads();

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

std::string quad_labels(int lane) {
    const Quad q = QUADS[std::size_t(lane)];
    return std::to_string(O[q.a]) + ',' + std::to_string(O[q.b]) + ',' +
           std::to_string(O[q.c]) + ',' + std::to_string(O[q.d]);
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
#pragma omp parallel for schedule(static)
        for (std::int64_t mask = 0; mask < std::int64_t(STATES); ++mask) {
            if ((u32(mask) & (u32{1} << bit)) == 0) continue;
            Integer* upper = table.data() + std::size_t(mask) * lanes;
            const Integer* lower =
                table.data() + std::size_t(u32(mask) ^ (u32{1} << bit)) * lanes;
            for (int lane = 0; lane < lanes; ++lane) upper[lane] += lower[lane];
        }
    }
}

struct Geometry { u64 denominator; std::vector<Event> events; u32 initial; };

Geometry geometry() {
    std::vector<int> labels;
    labels.insert(labels.end(), C.begin(), C.end());
    labels.insert(labels.end(), O.begin(), O.end());
    labels.push_back(50);
    u64 denominator = 1;
    for (int speed : labels) denominator = checked_lcm(denominator, 14 * u64(speed));
    return {denominator, grouped_events(labels, denominator), (u32{1} << 31) - 1};
}

std::vector<u32> fixed_masks(const std::vector<int>& active) {
    const u32 bit50 = u32{1} << 30;
    std::vector<u32> answer;
    answer.reserve(active.size());
    for (int lane : active) {
        const Quad q = QUADS[std::size_t(lane)];
        answer.push_back(bit50 | (u32{1} << (NC + q.a)) |
                         (u32{1} << (NC + q.b)) | (u32{1} << (NC + q.c)) |
                         (u32{1} << (NC + q.d)));
    }
    return answer;
}

struct MassTable { u64 denominator; int lanes; std::vector<u64> values; };

MassTable mass_table(const std::vector<int>& active) {
    const Geometry geo = geometry();
    const int lanes = int(active.size());
    const std::vector<u32> fixed = fixed_masks(active);
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
    return {geo.denominator, lanes, std::move(values)};
}

struct ComponentTable { u64 denominator; int lanes; std::vector<i32> values; };

ComponentTable component_table(const std::vector<int>& active) {
    const Geometry geo = geometry();
    const int lanes = int(active.size());
    const std::vector<u32> fixed = fixed_masks(active);
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
    return {geo.denominator, lanes, std::move(values)};
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

std::vector<int> body_speeds(u32 mask, int lane) {
    std::vector<int> answer;
    for (int bit = 0; bit < NC; ++bit) if (mask & (u32{1} << bit)) answer.push_back(C[bit]);
    const Quad q = QUADS[std::size_t(lane)];
    answer.push_back(O[q.a]); answer.push_back(O[q.b]);
    answer.push_back(O[q.c]); answer.push_back(O[q.d]);
    answer.push_back(50);
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

}  // namespace

int main() {
    const std::vector<u32> choose5 = masks_of_size(5);
    if (choose5.size() != 8568) throw std::runtime_error("choose(18,5)");
    std::array<BaseResult, NQUADS> results{};
    u64 denominator = 0;

    for (int first = 0; first < NQUADS; first += BATCH) {
        const int last = std::min(first + BATCH, NQUADS);
        std::vector<int> active(std::size_t(last - first));
        std::iota(active.begin(), active.end(), first);
        const MassTable mass = mass_table(active);
        const ComponentTable components = component_table(active);
        if (mass.denominator != UINT64_C(91205797082400) ||
            components.denominator != mass.denominator ||
            mass.lanes != last - first || components.lanes != last - first) {
            throw std::runtime_error("batch geometry");
        }
        denominator = mass.denominator;
        for (u32 mask : choose5) {
            const std::size_t row = std::size_t(CORE_MASK ^ mask) * active.size();
            for (int local = 0; local < last - first; ++local) {
                BaseResult& result = results[std::size_t(first + local)];
                const u64 numerator = mass.values[row + std::size_t(local)];
                const i128 delta = 54 * i128(numerator) - 4 * i128(denominator);
                const i32 count = components.values[row + std::size_t(local)];
                if (count <= 0) throw std::runtime_error("nonpositive components");
                if (delta <= 0) ++result.nonstrict;
                if (result.min_delta < 0 || delta < result.min_delta) {
                    result.min_delta = delta;
                    result.min_mask = mask;
                    result.min_mass = numerator;
                }
                if (delta > 0) {
                    const i128 cutoff = ceil_positive(54 * i128(count) * denominator,
                                                       7 * delta);
                    if (cutoff > result.max_cutoff) {
                        result.max_cutoff = cutoff;
                        result.cutoff_mask = mask;
                        result.cutoff_mass = numerator;
                        result.cutoff_components = u64(count);
                    }
                }
            }
        }
    }

    int strict = 0, nonstrict = 0, global_cutoff = 0;
    int global_cutoff_lane = -1, best_lane = -1;
    for (int lane = 0; lane < NQUADS; ++lane) {
        BaseResult& result = results[std::size_t(lane)];
        if (result.nonstrict == 0) ++strict; else ++nonstrict;
        if (result.max_cutoff > global_cutoff) {
            global_cutoff = int(result.max_cutoff);
            global_cutoff_lane = lane;
        }
        if (best_lane < 0 || result.min_delta > results[std::size_t(best_lane)].min_delta) {
            best_lane = lane;
        }
        const DirectResult direct_min = direct_body_sweep(body_speeds(result.min_mask, lane));
        require_same_fraction(result.min_mass, denominator,
                              direct_min.numerator, direct_min.denominator, "min replay");
        const DirectResult direct_cutoff = direct_body_sweep(body_speeds(result.cutoff_mask, lane));
        require_same_fraction(result.cutoff_mass, denominator,
                              direct_cutoff.numerator, direct_cutoff.denominator, "cutoff replay");
        if (direct_cutoff.components != result.cutoff_components) {
            throw std::runtime_error("component replay");
        }
        std::cout << "QUADRUPLE labels=" << quad_labels(lane)
                  << " nonstrict=" << result.nonstrict
                  << " min_delta=" << decimal(result.min_delta)
                  << " min_body=" << core_body(result.min_mask)
                  << " max_cutoff=" << decimal(result.max_cutoff)
                  << " cutoff_body=" << core_body(result.cutoff_mask)
                  << " cutoff_components=" << result.cutoff_components
                  << " direct_replays=2\n";
    }
    std::cout << "SUMMARY denominator=" << denominator
              << " quadruples=" << NQUADS
              << " profiles=" << u64(NQUADS) * choose5.size()
              << " strict=" << strict << " nonstrict=" << nonstrict
              << " global_max_cutoff=" << global_cutoff
              << " global_cutoff_quadruple=" << quad_labels(global_cutoff_lane)
              << " best=" << quad_labels(best_lane)
              << " best_min_delta=" << decimal(results[std::size_t(best_lane)].min_delta)
              << " best_cutoff=" << decimal(results[std::size_t(best_lane)].max_cutoff)
              << " direct_replays=" << 2 * NQUADS
              << " status=FINITE-EXACT-INDEPENDENT-LIMITING LRC14=OPEN\n";
    return nonstrict == 0 ? 0 : 3;
}
