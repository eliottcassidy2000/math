// Exact hostile probe for extending THM-4211's 18-label fixed-50 chart
// by one of the twelve omitted labels of the ambient 30-label pool.
//
// Scope: for each p in P\C, scan outsider labels r in increasing order and
// test every newly created rank-nine body {p} union L, L in binom(C,8).
// Record the first literal direct-mass failure, together with the labelled
// singleton petal recovered by deleting p.  This is scratch evidence only.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using i64 = std::int64_t;
using u32 = std::uint32_t;
using i128 = __int128_t;

namespace {

constexpr std::array<int, 18> CHART = {
    8,16,40,42,80,84,85,88,95,120,126,143,145,168,193,240,252,286};
constexpr std::array<int, 30> POOL = {
    8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
    120,126,132,143,145,168,170,176,190,193,240,252,264,286,290};
constexpr std::array<int, 12> OMITTED = {
    10,15,20,30,60,63,132,170,176,190,264,290};
constexpr int FIXED_Q = 50;
constexpr int SCAN_MAX = 600;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

i64 exact_lcm(i64 left, i64 right) {
    const i64 common = std::gcd(left, right);
    const i128 value = static_cast<i128>(left / common) * right;
    require(value <= std::numeric_limits<i64>::max(), "lcm overflow");
    return static_cast<i64>(value);
}

bool safe_mid(int speed, i64 left, i64 right, i64 denominator) {
    i128 residue = static_cast<i128>(speed) * (left + right) %
                   (static_cast<i128>(2) * denominator);
    if (residue < 0) residue += static_cast<i128>(2) * denominator;
    return static_cast<i128>(7) * residue >= denominator &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * denominator;
}

std::string decimal(i128 value) {
    if (value == 0) return "0";
    bool negative = value < 0;
    if (negative) value = -value;
    std::string answer;
    while (value > 0) {
        answer.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (negative) answer.push_back('-');
    std::reverse(answer.begin(), answer.end());
    return answer;
}

bool in_pool(int value) {
    return std::find(POOL.begin(), POOL.end(), value) != POOL.end();
}

std::string labels(u32 mask, const std::array<int, 19>& vertices) {
    std::string answer;
    for (int bit = 0; bit < 19; ++bit) {
        if ((mask & (u32{1} << bit)) == 0) continue;
        if (!answer.empty()) answer += ',';
        answer += std::to_string(vertices[bit]);
    }
    return answer;
}

struct ExactRow {
    i64 denominator = 0;
    std::vector<i64> safe_mass_by_allowed_failure;
};

struct BaseProfileRow {
    i64 denominator = 0;
    std::vector<i64> mass_by_allowed_failure;
    std::vector<i64> starts_by_allowed_failure;
    std::vector<i64> continuations_by_allowed_failure;
};

BaseProfileRow build_base_profile_row(const std::array<int, 19>& vertices) {
    std::vector<int> speeds(vertices.begin(), vertices.end());
    speeds.push_back(FIXED_Q);
    std::sort(speeds.begin(), speeds.end());
    speeds.erase(std::unique(speeds.begin(), speeds.end()), speeds.end());

    i64 denominator = 1;
    for (int speed : speeds) denominator = exact_lcm(denominator, 14LL * speed);
    std::vector<i64> walls{0, denominator};
    for (int speed : speeds) {
        const i64 unit = denominator / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    }
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());

    struct CellState { i64 length; u32 failed; bool q_safe; };
    std::vector<CellState> cells;
    cells.reserve(walls.size() - 1);
    for (std::size_t index = 0; index + 1 < walls.size(); ++index) {
        u32 failed = 0;
        for (int bit = 0; bit < 19; ++bit) {
            if (!safe_mid(vertices[bit], walls[index], walls[index + 1], denominator)) {
                failed |= u32{1} << bit;
            }
        }
        cells.push_back({walls[index + 1] - walls[index], failed,
                         safe_mid(FIXED_Q, walls[index], walls[index + 1], denominator)});
    }

    constexpr u32 count = u32{1} << 19;
    std::vector<i64> mass(count, 0), starts(count, 0), continuations(count, 0);
    for (std::size_t index = 0; index < cells.size(); ++index) {
        const CellState& current = cells[index];
        const CellState& previous = cells[(index + cells.size() - 1) % cells.size()];
        if (!current.q_safe) continue;
        mass[current.failed] += current.length;
        starts[current.failed] += 1;
        if (previous.q_safe) continuations[current.failed | previous.failed] += 1;
    }
    for (int bit = 0; bit < 19; ++bit) {
        for (u32 mask = 0; mask < count; ++mask) {
            if ((mask & (u32{1} << bit)) == 0) continue;
            const u32 lower = mask ^ (u32{1} << bit);
            mass[mask] += mass[lower];
            starts[mask] += starts[lower];
            continuations[mask] += continuations[lower];
        }
    }
    return {denominator, std::move(mass), std::move(starts),
            std::move(continuations)};
}

ExactRow build_zeta_row(const std::array<int, 19>& vertices, int outsider) {
    std::vector<int> speeds(vertices.begin(), vertices.end());
    speeds.push_back(FIXED_Q);
    speeds.push_back(outsider);
    std::sort(speeds.begin(), speeds.end());
    speeds.erase(std::unique(speeds.begin(), speeds.end()), speeds.end());

    i64 denominator = 1;
    for (int speed : speeds) denominator = exact_lcm(denominator, 14LL * speed);

    std::vector<i64> walls{0, denominator};
    for (int speed : speeds) {
        const i64 unit = denominator / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    }
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());

    constexpr u32 count = u32{1} << 19;
    std::vector<i64> zeta(count, 0);
    for (std::size_t index = 0; index + 1 < walls.size(); ++index) {
        const i64 left = walls[index];
        const i64 right = walls[index + 1];
        if (!safe_mid(FIXED_Q, left, right, denominator) ||
            !safe_mid(outsider, left, right, denominator)) {
            continue;
        }
        u32 failed = 0;
        for (int bit = 0; bit < 19; ++bit) {
            if (!safe_mid(vertices[bit], left, right, denominator)) {
                failed |= u32{1} << bit;
            }
        }
        zeta[failed] += right - left;
    }

    for (int bit = 0; bit < 19; ++bit) {
        for (u32 mask = 0; mask < count; ++mask) {
            if ((mask & (u32{1} << bit)) != 0) {
                zeta[mask] += zeta[mask ^ (u32{1} << bit)];
            }
        }
    }
    return {denominator, std::move(zeta)};
}

struct Minimum {
    bool initialized = false;
    u32 body = 0;
    i64 mass = 0;
    i128 delta = 0;
};

Minimum minimum_new_body(const ExactRow& row) {
    constexpr u32 full = (u32{1} << 19) - 1;
    constexpr u32 pbit = u32{1} << 18;
    Minimum result;
    const u32 chart_full = (u32{1} << 18) - 1;
    for (u32 subset = chart_full;; subset = (subset - 1) & chart_full) {
        if (std::popcount(subset) == 8) {
            const u32 body = subset | pbit;
            const i64 mass = row.safe_mass_by_allowed_failure[full ^ body];
            const i128 delta = static_cast<i128>(63) * mass -
                               static_cast<i128>(4) * row.denominator;
            if (!result.initialized || delta < result.delta) {
                result = {true, body, mass, delta};
            }
        }
        if (subset == 0) break;
    }
    return result;
}

Minimum minimum_old_chart_body(const ExactRow& row) {
    constexpr u32 full = (u32{1} << 19) - 1;
    Minimum result;
    const u32 chart_full = (u32{1} << 18) - 1;
    for (u32 subset = chart_full;; subset = (subset - 1) & chart_full) {
        if (std::popcount(subset) == 9) {
            const i64 mass = row.safe_mass_by_allowed_failure[full ^ subset];
            const i128 delta = static_cast<i128>(63) * mass -
                               static_cast<i128>(4) * row.denominator;
            if (!result.initialized || delta < result.delta) {
                result = {true, subset, mass, delta};
            }
        }
        if (subset == 0) break;
    }
    return result;
}

}  // namespace

int main() {
    std::cout << "SCOPE fixed_q=" << FIXED_Q
              << " omitted_pool_labels=" << OMITTED.size()
              << " new_bodies_per_extension=43758 scan_max=" << SCAN_MAX << '\n';
    int refuted = 0;
    int completed = 0;
    std::uint64_t total_literal_outsiders = 0;
    std::uint64_t total_new_body_checks = 0;
    std::uint64_t total_old_control_checks = 0;
    for (int added : OMITTED) {
        std::array<int, 19> vertices{};
        std::copy(CHART.begin(), CHART.end(), vertices.begin());
        vertices[18] = added;
        const BaseProfileRow base = build_base_profile_row(vertices);
        constexpr u32 full = (u32{1} << 19) - 1;
        constexpr u32 pbit = u32{1} << 18;
        const u32 chart_full = (u32{1} << 18) - 1;
        int nonstrict = 0;
        int maximum_cutoff = 0;
        i128 minimum_limit_delta = 0;
        u32 minimum_limit_body = 0;
        u32 maximum_cutoff_body = 0;
        for (u32 subset = chart_full;; subset = (subset - 1) & chart_full) {
            if (std::popcount(subset) == 8) {
                const u32 body = subset | pbit;
                const u32 allowed = full ^ body;
                const i64 mass = base.mass_by_allowed_failure[allowed];
                const i64 components = base.starts_by_allowed_failure[allowed] -
                                       base.continuations_by_allowed_failure[allowed];
                require(components > 0, "base component formula failed");
                const i128 limit_delta = static_cast<i128>(54) * mass -
                                         static_cast<i128>(4) * base.denominator;
                if (minimum_limit_body == 0 || limit_delta < minimum_limit_delta) {
                    minimum_limit_delta = limit_delta;
                    minimum_limit_body = body;
                }
                if (limit_delta <= 0) {
                    ++nonstrict;
                } else {
                    const i128 numerator = static_cast<i128>(54) * components *
                                           base.denominator;
                    const i128 divisor = static_cast<i128>(7) * limit_delta;
                    const int cutoff = static_cast<int>((numerator + divisor - 1) /
                                                        divisor);
                    if (cutoff > maximum_cutoff) {
                        maximum_cutoff = cutoff;
                        maximum_cutoff_body = body;
                    }
                }
            }
            if (subset == 0) break;
        }
        i128 old_minimum_limit_delta = 0;
        u32 old_minimum_limit_body = 0;
        int old_maximum_cutoff = 0;
        for (u32 subset = chart_full;; subset = (subset - 1) & chart_full) {
            if (std::popcount(subset) == 9) {
                const u32 allowed = full ^ subset;
                const i64 mass = base.mass_by_allowed_failure[allowed];
                const i64 components = base.starts_by_allowed_failure[allowed] -
                                       base.continuations_by_allowed_failure[allowed];
                require(components > 0, "old-chart component formula failed");
                const i128 limit_delta = static_cast<i128>(54) * mass -
                                         static_cast<i128>(4) * base.denominator;
                if (old_minimum_limit_body == 0 ||
                    limit_delta < old_minimum_limit_delta) {
                    old_minimum_limit_delta = limit_delta;
                    old_minimum_limit_body = subset;
                }
                require(limit_delta > 0, "old chart lost strict limiting reserve");
                const i128 numerator = static_cast<i128>(54) * components *
                                       base.denominator;
                const i128 divisor = static_cast<i128>(7) * limit_delta;
                const int cutoff = static_cast<int>((numerator + divisor - 1) /
                                                    divisor);
                old_maximum_cutoff = std::max(old_maximum_cutoff, cutoff);
            }
            if (subset == 0) break;
        }
        require(base.denominator == 91205797082400LL,
                "fixed-50 base denominator changed");
        require(old_minimum_limit_delta == static_cast<i128>(517215373867152LL) &&
                    labels(old_minimum_limit_body, vertices) ==
                        "8,80,85,88,95,143,145,168,193" &&
                    old_maximum_cutoff == 448,
                "THM-4211 mass/component/cutoff control changed");
        std::cout << "BASE added=" << added
                  << " denominator=" << base.denominator
                  << " nonstrict=" << nonstrict
                  << " min_limit_delta=" << decimal(minimum_limit_delta)
                  << " min_limit_body=" << labels(minimum_limit_body, vertices)
                  << " max_cutoff=" << maximum_cutoff
                  << " max_cutoff_body=" << labels(maximum_cutoff_body, vertices)
                  << '\n' << std::flush;
        bool found = false;
        const int scan_ceiling = std::min(SCAN_MAX, std::max(0, maximum_cutoff - 1));
        int scanned_outsiders = 0;
        bool literal_min_initialized = false;
        i128 literal_min_delta = 0;
        i64 literal_min_denominator = 1;
        int literal_min_outsider = 0;
        u32 literal_min_body = 0;
        for (int outsider = 1; outsider <= scan_ceiling; ++outsider) {
            if (outsider == FIXED_Q || in_pool(outsider)) continue;
            ++scanned_outsiders;
            const ExactRow row = build_zeta_row(vertices, outsider);
            const Minimum old_control = minimum_old_chart_body(row);
            require(old_control.initialized && old_control.delta > 0,
                    "THM-4211 old-chart positive control failed");
            const Minimum fresh = minimum_new_body(row);
            require(fresh.initialized, "new-body universe empty");
            if (!literal_min_initialized ||
                fresh.delta * literal_min_denominator <
                    literal_min_delta * row.denominator) {
                literal_min_initialized = true;
                literal_min_delta = fresh.delta;
                literal_min_denominator = row.denominator;
                literal_min_outsider = outsider;
                literal_min_body = fresh.body;
            }
            if (fresh.delta >= 0) continue;

            const u32 reduced = fresh.body ^ pbit;
            const i64 reduced_mass = row.safe_mass_by_allowed_failure[full ^ reduced];
            const i128 reduced_delta = static_cast<i128>(63) * reduced_mass -
                                       static_cast<i128>(4) * row.denominator;
            const i64 singleton_petal = reduced_mass - fresh.mass;
            require(singleton_petal >= 0 && reduced_delta > 0,
                    "labelled deletion petal control failed");

            // Independent consequence extraction: re-evaluate the selected body
            // directly from the pre-zeta exact cells by rebuilding the same literal
            // wall row with a singleton query encoded as an allowed-failure mask.
            const i64 replay_mass =
                row.safe_mass_by_allowed_failure[full ^ fresh.body];
            require(replay_mass == fresh.mass, "selected-body replay mismatch");

            std::cout << "EXTENSION added=" << added
                      << " first_failure_r=" << outsider
                      << " body=" << labels(fresh.body, vertices)
                      << " delta=" << decimal(fresh.delta)
                      << " denominator=" << row.denominator
                      << " delete_added_delta=" << decimal(reduced_delta)
                      << " singleton_petal_numerator=" << singleton_petal
                      << " old_chart_min_delta=" << decimal(old_control.delta)
                      << " old_chart_min_body=" << labels(old_control.body, vertices)
                      << '\n';
            found = true;
            ++refuted;
            break;
        }
        if (!found) {
            const bool cofinal_complete =
                nonstrict == 0 && scan_ceiling >= maximum_cutoff - 1;
            if (cofinal_complete) ++completed;
            total_literal_outsiders += static_cast<std::uint64_t>(scanned_outsiders);
            total_new_body_checks += static_cast<std::uint64_t>(scanned_outsiders) *
                                     UINT64_C(43758);
            total_old_control_checks += static_cast<std::uint64_t>(scanned_outsiders) *
                                        UINT64_C(48620);
            std::cout << "EXTENSION added=" << added
                      << " no_failure_through=" << scan_ceiling
                      << " analytic_cutoff=" << maximum_cutoff
                      << " literal_outsiders=" << scanned_outsiders
                      << " new_body_checks="
                      << static_cast<std::uint64_t>(scanned_outsiders) *
                             UINT64_C(43758)
                      << " literal_min_delta=" << decimal(literal_min_delta)
                      << " literal_min_denominator=" << literal_min_denominator
                      << " literal_min_r=" << literal_min_outsider
                      << " literal_min_body=" << labels(literal_min_body, vertices)
                      << " cofinal_complete="
                      << (cofinal_complete ? 1 : 0)
                      << '\n' << std::flush;
        }
    }
    std::cout << "SUMMARY refuted_extensions=" << refuted
              << " completed_extensions=" << completed
              << " unresolved_extensions=" << (OMITTED.size() - refuted - completed)
              << " analytic_new_profiles=" << OMITTED.size() * UINT64_C(43758)
              << " literal_outsiders=" << total_literal_outsiders
              << " literal_new_body_checks=" << total_new_body_checks
              << " literal_old_control_checks=" << total_old_control_checks
              << '\n';
    return 0;
}
