#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

namespace {

constexpr int NC = 18;
constexpr int NO = 12;
constexpr int NPAIRS = 66;
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

struct PairIndex {
    int first;
    int second;
};

std::array<PairIndex, NPAIRS> make_pairs() {
    std::array<PairIndex, NPAIRS> answer{};
    int lane = 0;
    for (int first = 0; first < NO; ++first) {
        for (int second = first + 1; second < NO; ++second) {
            answer[lane++] = {first, second};
        }
    }
    if (lane != NPAIRS) throw std::runtime_error("pair universe");
    return answer;
}

const std::array<PairIndex, NPAIRS> PAIRS = make_pairs();

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
    u64 divisor = std::gcd(first, second);
    i128 value = i128(first / divisor) * second;
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
            throw std::runtime_error("nonintegral event grid");
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

struct Tables {
    u64 denominator = 0;
    int lanes = 0;
    std::vector<u64> mass;
    std::vector<u32> starts;
    std::vector<u32> bridges;
};

Tables sweep_all_pairs(int outsider, const std::vector<int>& active_pairs,
                       bool components) {
    std::vector<int> labels;
    labels.insert(labels.end(), C.begin(), C.end());
    labels.insert(labels.end(), O.begin(), O.end());
    labels.push_back(50);
    if (outsider > 0) labels.push_back(outsider);
    if ((outsider == 0 && labels.size() != 31) ||
        (outsider > 0 && labels.size() != 32)) {
        throw std::runtime_error("label-bit universe");
    }

    u64 denominator = 1;
    for (int speed : labels) denominator = checked_lcm(denominator, 14 * u64(speed));
    const std::vector<Event> events = grouped_events(labels, denominator);
    const int lanes = int(active_pairs.size());
    if (lanes == 0) throw std::runtime_error("empty active pair set");

    std::vector<u32> fixed_masks;
    fixed_masks.reserve(std::size_t(lanes));
    const u32 bit50 = u32{1} << 30;
    const u32 bitr = outsider > 0 ? (u32{1} << 31) : 0;
    for (int pair_lane : active_pairs) {
        const PairIndex pair = PAIRS[pair_lane];
        fixed_masks.push_back(bit50 | bitr |
                              (u32{1} << (NC + pair.first)) |
                              (u32{1} << (NC + pair.second)));
    }

    Tables answer;
    answer.denominator = denominator;
    answer.lanes = lanes;
    answer.mass.assign(STATES * std::size_t(lanes), 0);
    if (components) {
        answer.starts.assign(STATES * std::size_t(lanes), 0);
        answer.bridges.assign(STATES * std::size_t(lanes), 0);
    }

    const u32 initial = labels.size() == 32
        ? std::numeric_limits<u32>::max()
        : (u32{1} << labels.size()) - 1;
    u32 state = initial;
    u64 previous = 0;

    auto add_interval = [&](u64 length) {
        if (length == 0) return;
        const std::size_t base = std::size_t(state & CORE_MASK) * lanes;
        for (int lane = 0; lane < lanes; ++lane) {
            if ((state & fixed_masks[lane]) == 0) {
                answer.mass[base + std::size_t(lane)] += length;
            }
        }
    };

    for (const Event& event : events) {
        if (event.position < previous) throw std::runtime_error("event order");
        add_interval(event.position - previous);
        const u32 before = state;
        state ^= event.toggle;
        if (components) {
            const std::size_t start_base = std::size_t(state & CORE_MASK) * lanes;
            const std::size_t bridge_base =
                std::size_t((before | state) & CORE_MASK) * lanes;
            for (int lane = 0; lane < lanes; ++lane) {
                const bool before_safe = (before & fixed_masks[lane]) == 0;
                const bool after_safe = (state & fixed_masks[lane]) == 0;
                if (after_safe) ++answer.starts[start_base + std::size_t(lane)];
                if (before_safe && after_safe) {
                    ++answer.bridges[bridge_base + std::size_t(lane)];
                }
            }
        }
        previous = event.position;
    }
    add_interval(denominator - previous);
    if (state != initial) throw std::runtime_error("circular toggle mismatch");

    subset_zeta(answer.mass, lanes);
    if (components) {
        subset_zeta(answer.starts, lanes);
        subset_zeta(answer.bridges, lanes);
    }
    return answer;
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

struct LiteralOne {
    bool active = false;
    u64 failures = 0;
    u64 equalities = 0;
    i128 min_margin = -1;
    u64 denominator = 1;
    u64 min_mass = 0;
    u32 min_mask = 0;
};

struct LiteralRow {
    std::array<LiteralOne, NPAIRS> pairs;
};

struct LiteralAggregate {
    u64 rows = 0;
    u64 checks = 0;
    u64 failures = 0;
    u64 equalities = 0;
    i128 min_margin = -1;
    u64 denominator = 1;
    u64 min_mass = 0;
    int min_r = 0;
    u32 min_mask = 0;
};

std::vector<int> body_speeds(u32 core_mask, int pair_lane, int outsider) {
    std::vector<int> speeds;
    for (int bit = 0; bit < NC; ++bit) {
        if (core_mask & (u32{1} << bit)) speeds.push_back(C[bit]);
    }
    const PairIndex pair = PAIRS[pair_lane];
    speeds.push_back(O[pair.first]);
    speeds.push_back(O[pair.second]);
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

}  // namespace

int main() {
    const std::vector<u32> choose7 = masks_of_size(7);
    if (choose7.size() != 31824) throw std::runtime_error("choose(18,7)");

    const DirectResult one_speed = direct_body_sweep({1});
    if (one_speed.numerator * 7 != one_speed.denominator * 6 ||
        one_speed.components != 1) {
        throw std::runtime_error("one-speed positive control");
    }

    std::vector<int> all_pairs(NPAIRS);
    std::iota(all_pairs.begin(), all_pairs.end(), 0);
    const Tables base = sweep_all_pairs(0, all_pairs, true);
    if (base.denominator != UINT64_C(91205797082400) || base.lanes != NPAIRS) {
        throw std::runtime_error("base chart denominator/lanes");
    }

    std::array<BaseResult, NPAIRS> base_results;
    u64 profiles = 0;
    for (u32 mask : choose7) {
        const std::size_t row = std::size_t(CORE_MASK ^ mask) * NPAIRS;
        for (int lane = 0; lane < NPAIRS; ++lane) {
            BaseResult& result = base_results[lane];
            const u64 mass = base.mass[row + std::size_t(lane)];
            const u64 components =
                u64(base.starts[row + std::size_t(lane)]) -
                u64(base.bridges[row + std::size_t(lane)]);
            if (components == 0) throw std::runtime_error("zero base components");
            const i128 delta = 54 * i128(mass) - 4 * i128(base.denominator);
            ++profiles;
            if (delta <= 0) ++result.nonstrict;
            if (result.min_delta < 0 || delta < result.min_delta) {
                result.min_delta = delta;
                result.min_mask = mask;
                result.min_mass = mass;
            }
            if (delta > 0) {
                const i128 cutoff = ceil_positive(
                    54 * i128(components) * base.denominator, 7 * delta);
                if (cutoff > result.max_cutoff) {
                    result.max_cutoff = cutoff;
                    result.cutoff_mask = mask;
                    result.cutoff_mass = mass;
                    result.cutoff_components = components;
                }
            }
        }
    }
    if (profiles != UINT64_C(2100384)) throw std::runtime_error("profile count");

    int global_cutoff = 0;
    for (const BaseResult& result : base_results) {
        if (result.nonstrict != 0 || result.min_delta <= 0) {
            throw std::runtime_error("nonstrict base profile");
        }
        global_cutoff = std::max(global_cutoff, int(result.max_cutoff));
    }
    if (global_cutoff != 563) throw std::runtime_error("global cutoff");

    std::vector<LiteralRow> literal_rows{std::size_t(global_cutoff)};

#pragma omp parallel for schedule(dynamic, 1)
    for (int outsider = 1; outsider < global_cutoff; ++outsider) {
        if (excluded_outsider(outsider)) continue;
        std::vector<int> active;
        for (int lane = 0; lane < NPAIRS; ++lane) {
            if (outsider < base_results[lane].max_cutoff) active.push_back(lane);
        }
        if (active.empty()) continue;
        const Tables table = sweep_all_pairs(outsider, active, false);
        LiteralRow row_result;
        for (int pair_lane : active) row_result.pairs[pair_lane].active = true;

        for (u32 mask : choose7) {
            const std::size_t row = std::size_t(CORE_MASK ^ mask) * table.lanes;
            for (int local = 0; local < table.lanes; ++local) {
                const int pair_lane = active[local];
                LiteralOne& result = row_result.pairs[pair_lane];
                const u64 mass = table.mass[row + std::size_t(local)];
                const i128 margin = 63 * i128(mass) - 4 * i128(table.denominator);
                if (margin < 0) ++result.failures;
                if (margin == 0) ++result.equalities;
                if (result.min_margin < 0 || margin < result.min_margin) {
                    result.min_margin = margin;
                    result.denominator = table.denominator;
                    result.min_mass = mass;
                    result.min_mask = mask;
                }
            }
        }
        literal_rows[std::size_t(outsider)] = std::move(row_result);
    }

    std::array<LiteralAggregate, NPAIRS> literal;
    for (int outsider = 1; outsider < global_cutoff; ++outsider) {
        const LiteralRow& row = literal_rows[std::size_t(outsider)];
        for (int lane = 0; lane < NPAIRS; ++lane) {
            const LiteralOne& source = row.pairs[lane];
            if (!source.active) continue;
            LiteralAggregate& target = literal[lane];
            ++target.rows;
            target.checks += choose7.size();
            target.failures += source.failures;
            target.equalities += source.equalities;
            if (target.min_margin < 0 ||
                source.min_margin * target.denominator <
                    target.min_margin * source.denominator) {
                target.min_margin = source.min_margin;
                target.denominator = source.denominator;
                target.min_mass = source.min_mass;
                target.min_r = outsider;
                target.min_mask = source.min_mask;
            }
        }
    }

    // A body-local event path replays all three extremizers for every pair.
    for (int lane = 0; lane < NPAIRS; ++lane) {
        const BaseResult& b = base_results[lane];
        DirectResult direct_min = direct_body_sweep(body_speeds(b.min_mask, lane, 0));
        require_same_fraction(b.min_mass, base.denominator,
                              direct_min.numerator, direct_min.denominator,
                              "base minimum direct replay");

        DirectResult direct_cutoff =
            direct_body_sweep(body_speeds(b.cutoff_mask, lane, 0));
        require_same_fraction(b.cutoff_mass, base.denominator,
                              direct_cutoff.numerator, direct_cutoff.denominator,
                              "cutoff mass direct replay");
        if (direct_cutoff.components != b.cutoff_components) {
            throw std::runtime_error("cutoff component direct replay");
        }

        const LiteralAggregate& l = literal[lane];
        DirectResult direct_literal =
            direct_body_sweep(body_speeds(l.min_mask, lane, l.min_r));
        require_same_fraction(l.min_mass, l.denominator,
                              direct_literal.numerator, direct_literal.denominator,
                              "literal minimum direct replay");
    }

    u64 total_rows = 0;
    u64 total_checks = 0;
    u64 total_failures = 0;
    u64 total_equalities = 0;
    for (int lane = 0; lane < NPAIRS; ++lane) {
        const PairIndex pair = PAIRS[lane];
        const BaseResult& b = base_results[lane];
        const LiteralAggregate& l = literal[lane];
        if (l.checks != l.rows * choose7.size()) {
            throw std::runtime_error("literal row/check product");
        }
        total_rows += l.rows;
        total_checks += l.checks;
        total_failures += l.failures;
        total_equalities += l.equalities;
        std::cout << "PAIR p=" << O[pair.first]
                  << " q=" << O[pair.second]
                  << " nonstrict=" << b.nonstrict
                  << " min_delta=" << decimal(b.min_delta)
                  << " min_body=" << core_body(b.min_mask)
                  << " max_cutoff=" << decimal(b.max_cutoff)
                  << " cutoff_body=" << core_body(b.cutoff_mask)
                  << " literal_rows=" << l.rows
                  << " literal_checks=" << l.checks
                  << " literal_failures=" << l.failures
                  << " literal_equalities=" << l.equalities
                  << " literal_min_delta=" << decimal(l.min_margin)
                  << " literal_min_denominator=" << l.denominator
                  << " literal_min_r=" << l.min_r
                  << " literal_min_body=" << core_body(l.min_mask)
                  << '\n';
    }

    if (total_rows != UINT64_C(29042) ||
        total_checks != UINT64_C(924232608) ||
        total_failures != 0 || total_equalities != 0) {
        throw std::runtime_error("literal aggregate totals");
    }
    std::cout << "SUMMARY denominator=" << base.denominator
              << " pairs=" << NPAIRS
              << " profiles=" << profiles
              << " strict_pairs=" << NPAIRS
              << " nonstrict_pairs=0"
              << " global_max_cutoff=" << global_cutoff
              << " literal_rows=" << total_rows
              << " literal_checks=" << total_checks
              << " literal_failures=" << total_failures
              << " literal_equalities=" << total_equalities
              << " direct_extremal_replays=" << 3 * NPAIRS
              << " twenty_label_charts=66"
              << " all_omitted_union_claim=0"
              << " status=FINITE-EXACT-INDEPENDENT"
              << " LRC14=OPEN\n";
    return 0;
}
