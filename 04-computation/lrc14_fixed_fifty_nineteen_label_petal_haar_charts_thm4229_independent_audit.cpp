#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using u64 = std::uint64_t;
using u32 = std::uint32_t;
using i128 = __int128_t;

namespace {

constexpr int NC = 18;
constexpr int NP = 12;
constexpr int LANES = NP + 1;  // twelve petals plus inherited C-chart control
constexpr u32 CORE_MASK = (u32{1} << NC) - 1;
constexpr std::size_t STATES = std::size_t{1} << NC;

const std::array<int, NC> C = {
    8, 16, 40, 42, 80, 84, 85, 88, 95,
    120, 126, 143, 145, 168, 193, 240, 252, 286,
};
const std::array<int, NP> O = {
    10, 15, 20, 30, 60, 63, 132, 170, 176, 190, 264, 290,
};
const std::array<int, 30> P = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
};

std::string to_string_i128(i128 value) {
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

bool in_pool_or_50(int value) {
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

std::string body_string(u32 mask, int petal) {
    std::string answer;
    for (int index = 0; index < NC; ++index) {
        if (!(mask & (u32{1} << index))) continue;
        if (!answer.empty()) answer += ',';
        answer += std::to_string(C[index]);
    }
    if (petal >= 0) {
        if (!answer.empty()) answer += ',';
        answer += std::to_string(O[petal]);
    }
    return answer;
}

struct Event {
    u64 position;
    u32 toggle;
};

std::vector<Event> grouped_events(const std::vector<int>& labels, u64 denominator) {
    std::vector<Event> raw;
    u64 event_bound = 0;
    for (int speed : labels) event_bound += 2 * u64(speed);
    raw.reserve(std::size_t(event_bound));
    for (std::size_t bit = 0; bit < labels.size(); ++bit) {
        int speed = labels[bit];
        u64 local_denominator = 14 * u64(speed);
        if (denominator % local_denominator != 0) {
            throw std::runtime_error("nonintegral event grid");
        }
        u64 scale = denominator / local_denominator;
        for (int k = 0; k < speed; ++k) {
            int left = 14 * k - 1;
            if (left < 0) left += int(local_denominator);
            int right = 14 * k + 1;
            raw.push_back({u64(left) * scale, u32{1} << bit});
            raw.push_back({u64(right) * scale, u32{1} << bit});
        }
    }
    std::sort(raw.begin(), raw.end(), [](const Event& first, const Event& second) {
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
    return answer;
}

void zeta_complement(std::vector<u64>& table) {
    for (int bit = 0; bit < NC; ++bit) {
        std::size_t step = std::size_t{1} << bit;
        std::size_t width = 2 * step;
        for (std::size_t block = 0; block < STATES; block += width) {
            for (std::size_t offset = 0; offset < step; ++offset) {
                std::size_t source = (block + offset) * LANES;
                std::size_t target = (block + step + offset) * LANES;
                for (int lane = 0; lane < LANES; ++lane) {
                    table[target + lane] += table[source + lane];
                }
            }
        }
    }
}

struct SweepTables {
    u64 denominator = 0;
    std::vector<u64> mass;
    std::vector<u64> starts;
    std::vector<u64> bridges;
};

SweepTables chart_sweep(int outsider, bool components) {
    std::vector<int> labels;
    labels.insert(labels.end(), C.begin(), C.end());
    labels.insert(labels.end(), O.begin(), O.end());
    labels.push_back(50);
    if (outsider > 0) labels.push_back(outsider);
    if (labels.size() > 32) throw std::runtime_error("mask width exceeded");

    u64 denominator = 1;
    for (int speed : labels) denominator = checked_lcm(denominator, 14 * u64(speed));
    std::vector<Event> events = grouped_events(labels, denominator);

    SweepTables answer;
    answer.denominator = denominator;
    answer.mass.assign(STATES * LANES, 0);
    if (components) {
        answer.starts.assign(STATES * LANES, 0);
        answer.bridges.assign(STATES * LANES, 0);
    }

    const u32 initial = labels.size() == 32
        ? std::numeric_limits<u32>::max()
        : (u32{1} << labels.size()) - 1;
    u32 state = initial;
    u64 previous = 0;
    const u32 bit50 = u32{1} << 30;
    const u32 bitr = outsider > 0 ? (u32{1} << 31) : 0;

    auto add_interval = [&](u64 length) {
        if (length == 0) return;
        u32 core_failure = state & CORE_MASK;
        std::size_t base = std::size_t(core_failure) * LANES;
        for (int petal = 0; petal < NP; ++petal) {
            u32 fixed = bit50 | bitr | (u32{1} << (NC + petal));
            if ((state & fixed) == 0) answer.mass[base + petal] += length;
        }
        if ((state & (bit50 | bitr)) == 0) {
            answer.mass[base + NP] += length;
        }
    };

    for (const Event& event : events) {
        if (event.position < previous) throw std::runtime_error("event order");
        add_interval(event.position - previous);
        u32 before = state;
        state ^= event.toggle;
        if (components) {
            u32 before_core = before & CORE_MASK;
            u32 after_core = state & CORE_MASK;
            std::size_t after_base = std::size_t(after_core) * LANES;
            std::size_t union_base = std::size_t(before_core | after_core) * LANES;
            for (int petal = 0; petal < NP; ++petal) {
                u32 fixed = bit50 | (u32{1} << (NC + petal));
                bool before_safe = (before & fixed) == 0;
                bool after_safe = (state & fixed) == 0;
                if (after_safe) ++answer.starts[after_base + petal];
                if (before_safe && after_safe) ++answer.bridges[union_base + petal];
            }
            bool before_safe = (before & bit50) == 0;
            bool after_safe = (state & bit50) == 0;
            if (after_safe) ++answer.starts[after_base + NP];
            if (before_safe && after_safe) ++answer.bridges[union_base + NP];
        }
        previous = event.position;
    }
    add_interval(denominator - previous);
    if (state != initial) throw std::runtime_error("circular toggle mismatch");

    zeta_complement(answer.mass);
    if (components) {
        zeta_complement(answer.starts);
        zeta_complement(answer.bridges);
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
    std::vector<Event> events = grouped_events(speeds, denominator);
    u32 initial = (u32{1} << speeds.size()) - 1;
    u32 state = initial;
    u64 previous = 0;
    u64 numerator = 0;
    u64 components = 0;
    for (const Event& event : events) {
        if (state == 0) numerator += event.position - previous;
        u32 before = state;
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
    i128 max_cutoff = -1;
    u32 cutoff_mask = 0;
    u64 cutoff_components = 0;
    u64 cutoff_mass = 0;
};

struct LiteralPetal {
    u64 checks = 0;
    u64 outsider_rows = 0;
    u64 failures = 0;
    u64 equalities = 0;
    i128 min_margin = -1;
    u64 min_denominator = 1;
    int min_r = 0;
    u32 min_mask = 0;
    u64 min_mass = 0;
};

struct LiteralRow {
    std::array<LiteralPetal, NP> petals;
    u64 old_checks = 0;
    u64 old_failures = 0;
    u64 old_equalities = 0;
};

std::vector<int> speeds_from_mask(u32 mask, int petal, int outsider) {
    std::vector<int> speeds;
    for (int index = 0; index < NC; ++index) {
        if (mask & (u32{1} << index)) speeds.push_back(C[index]);
    }
    if (petal >= 0) speeds.push_back(O[petal]);
    speeds.push_back(50);
    if (outsider > 0) speeds.push_back(outsider);
    return speeds;
}

void require_equal_fraction(u64 first_num, u64 first_den,
                            u64 second_num, u64 second_den,
                            const std::string& message) {
    if (i128(first_num) * second_den != i128(second_num) * first_den) {
        throw std::runtime_error(message);
    }
}

}  // namespace

int main() {
    const std::vector<u32> choose8 = masks_of_size(8);
    const std::vector<u32> choose9 = masks_of_size(9);
    if (choose8.size() != 43758 || choose9.size() != 48620) {
        throw std::runtime_error("binomial universes");
    }

    DirectResult one_speed = direct_body_sweep({1});
    if (one_speed.numerator * 7 != one_speed.denominator * 6
        || one_speed.components != 1) {
        throw std::runtime_error("one-speed 6/7 control");
    }

    SweepTables base = chart_sweep(0, true);
    if (base.denominator != 91205797082400ULL) {
        throw std::runtime_error("base denominator");
    }

    std::array<BaseResult, NP> base_results;
    u64 analytic_profiles = 0;
    for (int petal = 0; petal < NP; ++petal) {
        BaseResult result;
        for (u32 mask : choose8) {
            u32 complement = CORE_MASK ^ mask;
            std::size_t index = std::size_t(complement) * LANES + petal;
            u64 mass = base.mass[index];
            u64 components = base.starts[index] - base.bridges[index];
            i128 delta = 54 * i128(mass) - 4 * i128(base.denominator);
            ++analytic_profiles;
            if (delta <= 0) throw std::runtime_error("nonstrict limiting margin");
            i128 cutoff = ceil_positive(
                54 * i128(components) * base.denominator,
                7 * delta
            );
            if (result.min_delta < 0 || delta < result.min_delta) {
                result.min_delta = delta;
                result.min_mask = mask;
            }
            if (cutoff > result.max_cutoff) {
                result.max_cutoff = cutoff;
                result.cutoff_mask = mask;
                result.cutoff_components = components;
                result.cutoff_mass = mass;
            }
        }
        base_results[petal] = result;
    }
    if (analytic_profiles != 525096) throw std::runtime_error("analytic count");

    // Independent inherited C-chart control.
    i128 inherited_min_delta = -1;
    i128 inherited_max_cutoff = -1;
    u32 inherited_min_mask = 0;
    u32 inherited_cutoff_mask = 0;
    for (u32 mask : choose9) {
        u32 complement = CORE_MASK ^ mask;
        std::size_t index = std::size_t(complement) * LANES + NP;
        u64 mass = base.mass[index];
        u64 components = base.starts[index] - base.bridges[index];
        i128 delta = 54 * i128(mass) - 4 * i128(base.denominator);
        if (delta <= 0) throw std::runtime_error("inherited limiting margin");
        i128 cutoff = ceil_positive(
            54 * i128(components) * base.denominator,
            7 * delta
        );
        if (inherited_min_delta < 0 || delta < inherited_min_delta) {
            inherited_min_delta = delta;
            inherited_min_mask = mask;
        }
        if (cutoff > inherited_max_cutoff) {
            inherited_max_cutoff = cutoff;
            inherited_cutoff_mask = mask;
        }
    }

    int maximum_cutoff = 0;
    for (const BaseResult& result : base_results) {
        maximum_cutoff = std::max(maximum_cutoff, int(result.max_cutoff));
    }
    std::vector<LiteralRow> literal_rows(maximum_cutoff);

    // Each r is independent.  OpenMP only changes scheduling; every result
    // is written to its unique row and aggregated deterministically below.
#pragma omp parallel
    {
#pragma omp for schedule(dynamic, 1)
        for (int r = 1; r < maximum_cutoff; ++r) {
            if (in_pool_or_50(r)) continue;
            SweepTables table = chart_sweep(r, false);
            LiteralRow row;
            for (int petal = 0; petal < NP; ++petal) {
                if (r >= base_results[petal].max_cutoff) continue;
                LiteralPetal& result = row.petals[petal];
                result.outsider_rows = 1;
                for (u32 mask : choose8) {
                    u32 complement = CORE_MASK ^ mask;
                    u64 mass = table.mass[std::size_t(complement) * LANES + petal];
                    i128 margin = 63 * i128(mass) - 4 * i128(table.denominator);
                    ++result.checks;
                    if (margin < 0) ++result.failures;
                    if (margin == 0) ++result.equalities;
                    if (result.min_margin < 0 || margin * result.min_denominator
                                                < result.min_margin * table.denominator) {
                        result.min_margin = margin;
                        result.min_denominator = table.denominator;
                        result.min_r = r;
                        result.min_mask = mask;
                        result.min_mass = mass;
                    }
                }
            }
            // Reconstruct each inherited body once for this literal r.  The
            // same values serve as controls for every petal whose cutoff uses r.
            for (u32 mask : choose9) {
                u32 complement = CORE_MASK ^ mask;
                u64 mass = table.mass[std::size_t(complement) * LANES + NP];
                i128 margin = 63 * i128(mass) - 4 * i128(table.denominator);
                ++row.old_checks;
                if (margin < 0) ++row.old_failures;
                if (margin == 0) ++row.old_equalities;
            }
            literal_rows[r] = row;
        }
    }

    std::array<LiteralPetal, NP> literal;
    u64 unique_old_checks = 0;
    u64 unique_old_failures = 0;
    u64 unique_old_equalities = 0;
    for (int r = 1; r < maximum_cutoff; ++r) {
        const LiteralRow& row = literal_rows[r];
        unique_old_checks += row.old_checks;
        unique_old_failures += row.old_failures;
        unique_old_equalities += row.old_equalities;
        for (int petal = 0; petal < NP; ++petal) {
            const LiteralPetal& source = row.petals[petal];
            LiteralPetal& target = literal[petal];
            target.checks += source.checks;
            target.outsider_rows += source.outsider_rows;
            target.failures += source.failures;
            target.equalities += source.equalities;
            if (source.checks && (target.min_margin < 0
                || source.min_margin * target.min_denominator
                   < target.min_margin * source.min_denominator)) {
                target.min_margin = source.min_margin;
                target.min_denominator = source.min_denominator;
                target.min_r = source.min_r;
                target.min_mask = source.min_mask;
                target.min_mass = source.min_mass;
            }
        }
    }

    // Body-local direct event sweeps independently replay every extremizer.
    for (int petal = 0; petal < NP; ++petal) {
        const BaseResult& b = base_results[petal];
        u32 min_complement = CORE_MASK ^ b.min_mask;
        u64 min_mass = base.mass[std::size_t(min_complement) * LANES + petal];
        DirectResult direct_min = direct_body_sweep(
            speeds_from_mask(b.min_mask, petal, 0));
        require_equal_fraction(min_mass, base.denominator,
                               direct_min.numerator, direct_min.denominator,
                               "direct limiting-min replay");

        u32 cutoff_complement = CORE_MASK ^ b.cutoff_mask;
        std::size_t cutoff_index = std::size_t(cutoff_complement) * LANES + petal;
        DirectResult direct_cutoff = direct_body_sweep(
            speeds_from_mask(b.cutoff_mask, petal, 0));
        require_equal_fraction(base.mass[cutoff_index], base.denominator,
                               direct_cutoff.numerator, direct_cutoff.denominator,
                               "direct cutoff-mass replay");
        if (direct_cutoff.components != base.starts[cutoff_index]
                                       - base.bridges[cutoff_index]) {
            throw std::runtime_error("direct component replay");
        }

        const LiteralPetal& l = literal[petal];
        DirectResult direct_literal = direct_body_sweep(
            speeds_from_mask(l.min_mask, petal, l.min_r));
        require_equal_fraction(l.min_mass, l.min_denominator,
                               direct_literal.numerator, direct_literal.denominator,
                               "direct literal replay");
    }

    DirectResult inherited_direct = direct_body_sweep(
        speeds_from_mask(inherited_min_mask, -1, 0));
    u64 inherited_mass = base.mass[
        std::size_t(CORE_MASK ^ inherited_min_mask) * LANES + NP];
    require_equal_fraction(inherited_mass, base.denominator,
                           inherited_direct.numerator, inherited_direct.denominator,
                           "direct inherited replay");

    u64 total_rows = 0;
    u64 total_new_checks = 0;
    u64 total_failures = 0;
    u64 total_equalities = 0;
    for (int petal = 0; petal < NP; ++petal) {
        const BaseResult& b = base_results[petal];
        const LiteralPetal& l = literal[petal];
        total_rows += l.outsider_rows;
        total_new_checks += l.checks;
        total_failures += l.failures;
        total_equalities += l.equalities;
        std::cout << "PETAL p=" << O[petal]
                  << " min_limit_delta=" << to_string_i128(b.min_delta)
                  << " min_limit_body=" << body_string(b.min_mask, petal)
                  << " max_cutoff=" << to_string_i128(b.max_cutoff)
                  << " max_cutoff_body=" << body_string(b.cutoff_mask, petal)
                  << " literal_rows=" << l.outsider_rows
                  << " literal_checks=" << l.checks
                  << " literal_failures=" << l.failures
                  << " literal_equalities=" << l.equalities
                  << " literal_min_margin=" << to_string_i128(l.min_margin)
                  << '/' << l.min_denominator
                  << " literal_min_r=" << l.min_r
                  << " literal_min_body=" << body_string(l.min_mask, petal)
                  << '\n';
    }
    std::cout << "INHERITED min_limit_delta=" << to_string_i128(inherited_min_delta)
              << " min_limit_body=" << body_string(inherited_min_mask, -1)
              << " max_cutoff=" << to_string_i128(inherited_max_cutoff)
              << " max_cutoff_body=" << body_string(inherited_cutoff_mask, -1)
              << " unique_literal_checks=" << unique_old_checks
              << " failures=" << unique_old_failures
              << " equalities=" << unique_old_equalities << '\n';
    std::cout << "SUMMARY analytic_profiles=" << analytic_profiles
              << " literal_rows=" << total_rows
              << " literal_new_checks=" << total_new_checks
              << " failures=" << total_failures
              << " equalities=" << total_equalities
              << " direct_extremal_replays=" << (3 * NP + 1)
              << " chart_union_claim=0"
              << " status=FINITE-EXACT-INDEPENDENT"
              << " LRC14=OPEN\n";
    return 0;
}
