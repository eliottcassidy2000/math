// Independent fixed-prefix/refined-atom audit for THM-4211.
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

#ifdef _OPENMP
#include <omp.h>
#endif

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

namespace {

constexpr std::array<int, 18> CHART = {
    8,16,40,42,80,84,85,88,95,120,126,143,145,168,193,240,252,286};
constexpr std::array<int, 30> POOL = {
    8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
    120,126,132,143,145,168,170,176,190,193,240,252,264,286,290};
constexpr int Q = 50;
constexpr u32 FULL = (u32{1} << 18) - 1;
constexpr std::array<int, 9> HOSTILE_BODY_VALUES =
    {16,85,88,95,143,145,168,193,240};

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

i64 exact_lcm(i64 a, i64 b) {
    const i64 g = std::gcd(a, b);
    const i128 product = static_cast<i128>(a / g) * b;
    require(product <= std::numeric_limits<i64>::max(), "lcm overflow");
    return static_cast<i64>(product);
}

std::string decimal(i128 value) {
    if (value == 0) return "0";
    bool negative = value < 0;
    if (negative) value = -value;
    std::string text;
    while (value != 0) {
        text.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (negative) text.push_back('-');
    std::reverse(text.begin(), text.end());
    return text;
}

std::string labels(u32 mask) {
    std::string text;
    for (int i = 0; i < 18; ++i) {
        if ((mask & (u32{1} << i)) == 0) continue;
        if (!text.empty()) text += ',';
        text += std::to_string(CHART[i]);
    }
    return text;
}

u32 mask_from_values(const std::array<int, 9>& values) {
    u32 mask = 0;
    for (int value : values) {
        const auto found = std::find(CHART.begin(), CHART.end(), value);
        require(found != CHART.end(), "requested body is outside chart");
        mask |= u32{1} << std::distance(CHART.begin(), found);
    }
    require(std::popcount(mask) == 9, "requested body is not a nine-set");
    return mask;
}

i128 exact_gcd(i128 a, i128 b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b != 0) {
        const i128 remainder = a % b;
        a = b;
        b = remainder;
    }
    return a;
}

bool safe_mid(int speed, i64 left, i64 right, i64 denominator) {
    const i128 period = static_cast<i128>(2) * denominator;
    i128 residue = static_cast<i128>(speed) * (left + right) % period;
    if (residue < 0) residue += period;
    return static_cast<i128>(7) * residue >= denominator &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * denominator;
}

i128 safe_prefix(i64 tick, int speed, i64 lattice) {
    const i128 product = static_cast<i128>(speed) * tick;
    const i128 whole = product / lattice;
    const i64 remainder = static_cast<i64>(product % lattice);
    const i128 scaled = static_cast<i128>(14) * remainder;
    i128 partial = 0;
    if (scaled >= static_cast<i128>(13) * lattice) {
        partial = static_cast<i128>(12) * lattice;
    } else if (scaled > lattice) {
        partial = scaled - lattice;
    }
    return static_cast<i128>(12) * whole * lattice + partial;
}

struct Cell {
    i64 left;
    i64 right;
    u32 failed;
    bool q_safe;
    int support;
};

struct Geometry {
    i64 denominator = 1;
    std::vector<Cell> cells;
    std::vector<u32> support_masks;
};

Geometry fixed_geometry() {
    Geometry geometry;
    std::vector<int> speeds(CHART.begin(), CHART.end());
    speeds.push_back(Q);
    for (int speed : speeds) {
        geometry.denominator = exact_lcm(geometry.denominator, 14LL * speed);
    }
    std::vector<i64> walls{0, geometry.denominator};
    for (int speed : speeds) {
        const i64 unit = geometry.denominator / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    }
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    std::vector<int> support_index(std::size_t{1} << 18, -1);
    geometry.cells.reserve(walls.size() - 1);
    for (std::size_t i = 0; i + 1 < walls.size(); ++i) {
        u32 failed = 0;
        for (int vertex = 0; vertex < 18; ++vertex) {
            if (!safe_mid(CHART[vertex], walls[i], walls[i + 1],
                          geometry.denominator)) {
                failed |= u32{1} << vertex;
            }
        }
        const bool q_safe = safe_mid(Q, walls[i], walls[i + 1],
                                     geometry.denominator);
        int index = -1;
        if (q_safe) {
            index = support_index[failed];
            if (index < 0) {
                index = static_cast<int>(geometry.support_masks.size());
                support_index[failed] = index;
                geometry.support_masks.push_back(failed);
            }
        }
        geometry.cells.push_back({walls[i], walls[i + 1], failed, q_safe, index});
    }
    return geometry;
}

struct Profile {
    u32 body;
    i64 base_mass = 0;
    int components = 0;
    i128 delta27 = 0;
    int cutoff = 0;
};

struct Row {
    int r = 0;
    int failures = 0;
    int equalities = 0;
    i128 minimum_delta = 0;
    u32 minimum_body = 0;
    i128 hostile_delta = 0;
};

Row evaluate_r(int r, const Geometry& geometry,
               const std::vector<Profile>& profiles,
               const std::vector<std::vector<int>>& targets) {
    std::vector<i128> atoms(geometry.support_masks.size(), 0);
    i128 total_qr_safe = 0;
    for (const Cell& cell : geometry.cells) {
        if (!cell.q_safe) continue;
        const i128 contribution =
            safe_prefix(cell.right, r, geometry.denominator) -
            safe_prefix(cell.left, r, geometry.denominator);
        require(contribution >= 0, "negative prefix contribution");
        atoms[cell.support] += contribution;
        total_qr_safe += contribution;
    }
    require(total_qr_safe > 0 &&
                total_qr_safe < static_cast<i128>(12) * r * geometry.denominator,
            "q/r safe total is degenerate");
    std::vector<i128> masses(profiles.size(), 0);
    for (std::size_t support = 0; support < atoms.size(); ++support) {
        if (atoms[support] == 0) continue;
        for (int body : targets[support]) masses[body] += atoms[support];
    }
    Row row;
    row.r = r;
    const i128 threshold_scale = static_cast<i128>(r) * geometry.denominator;
    const u32 hostile_body = mask_from_values(HOSTILE_BODY_VALUES);
    for (std::size_t body = 0; body < profiles.size(); ++body) {
        const i128 delta = static_cast<i128>(9) * masses[body] -
                           static_cast<i128>(8) * threshold_scale;
        if (body == 0 || delta < row.minimum_delta) {
            row.minimum_delta = delta;
            row.minimum_body = profiles[body].body;
        }
        if (delta < 0) ++row.failures;
        if (delta == 0) ++row.equalities;
        if (profiles[body].body == hostile_body) row.hostile_delta = delta;
    }
    return row;
}

bool excluded_r(int r) {
    if (r == Q) return true;
    return std::find(POOL.begin(), POOL.end(), r) != POOL.end();
}

}  // namespace

int main() {
    const Geometry geometry = fixed_geometry();
    std::vector<Profile> profiles;
    for (u32 mask = 0; mask <= FULL; ++mask) {
        if (std::popcount(mask) != 9) continue;
        profiles.push_back({mask});
    }
    require(profiles.size() == 48620, "body universe changed");

    int minimum_cutoff = 0;
    int maximum_cutoff = 0;
    u32 minimum_cutoff_body = 0;
    u32 maximum_cutoff_body = 0;
    int minimum_cutoff_count = 0;
    int maximum_cutoff_count = 0;
    int base_equalities = 0;
    int base_failures = 0;
    i128 minimum_delta27 = 0;
    u32 minimum_delta_body = 0;
    int minimum_delta_components = 0;
    i64 minimum_delta_mass = 0;
    for (std::size_t body_index_value = 0;
         body_index_value < profiles.size(); ++body_index_value) {
        Profile& profile = profiles[body_index_value];
        std::vector<bool> good;
        good.reserve(geometry.cells.size());
        for (const Cell& cell : geometry.cells) {
            const bool value = cell.q_safe && (cell.failed & profile.body) == 0;
            good.push_back(value);
            if (value) profile.base_mass += cell.right - cell.left;
        }
        bool previous = good.back();
        for (bool current : good) {
            if (current && !previous) ++profile.components;
            previous = current;
        }
        profile.delta27 = static_cast<i128>(27) * profile.base_mass -
                          static_cast<i128>(2) * geometry.denominator;
        if (body_index_value == 0 || profile.delta27 < minimum_delta27) {
            minimum_delta27 = profile.delta27;
            minimum_delta_body = profile.body;
            minimum_delta_components = profile.components;
            minimum_delta_mass = profile.base_mass;
        }
        if (profile.delta27 < 0) ++base_failures;
        if (profile.delta27 == 0) ++base_equalities;
        if (profile.delta27 > 0) {
            const i128 numerator = static_cast<i128>(27) * geometry.denominator *
                                   profile.components;
            const i128 denominator = static_cast<i128>(7) * profile.delta27;
            profile.cutoff = static_cast<int>((numerator + denominator - 1) /
                                              denominator);
            if (minimum_cutoff == 0 || profile.cutoff < minimum_cutoff) {
                minimum_cutoff = profile.cutoff;
                minimum_cutoff_body = profile.body;
                minimum_cutoff_count = 1;
            } else if (profile.cutoff == minimum_cutoff) {
                ++minimum_cutoff_count;
            }
            if (profile.cutoff > maximum_cutoff) {
                maximum_cutoff = profile.cutoff;
                maximum_cutoff_body = profile.body;
                maximum_cutoff_count = 1;
            } else if (profile.cutoff == maximum_cutoff) {
                ++maximum_cutoff_count;
            }
        }
    }

    std::vector<std::vector<int>> targets(geometry.support_masks.size());
    u64 target_memberships = 0;
    for (std::size_t support = 0; support < geometry.support_masks.size(); ++support) {
        const u32 failed = geometry.support_masks[support];
        for (std::size_t body = 0; body < profiles.size(); ++body) {
            if ((failed & profiles[body].body) == 0) {
                targets[support].push_back(static_cast<int>(body));
            }
        }
        target_memberships += targets[support].size();
    }

    std::vector<int> finite_r;
    for (int r = 1; r < maximum_cutoff; ++r) {
        if (!excluded_r(r)) finite_r.push_back(r);
    }
    std::vector<Row> rows(finite_r.size());

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (std::size_t i = 0; i < finite_r.size(); ++i) {
        rows[i] = evaluate_r(finite_r[i], geometry, profiles, targets);
    }

    int total_failures = 0;
    int total_equalities = 0;
    i128 minimum_normalized_delta = 0;
    int minimum_normalized_r = 0;
    u32 minimum_normalized_body = 0;
    for (const Row& row : rows) {
        total_failures += row.failures;
        total_equalities += row.equalities;
        if (minimum_normalized_r == 0 ||
            row.minimum_delta * minimum_normalized_r <
                minimum_normalized_delta * row.r) {
            minimum_normalized_delta = row.minimum_delta;
            minimum_normalized_r = row.r;
            minimum_normalized_body = row.minimum_body;
        }
    }
    const i128 minimum_fraction_denominator =
        static_cast<i128>(2) * minimum_normalized_r * geometry.denominator;
    const i128 minimum_fraction_gcd =
        exact_gcd(minimum_normalized_delta, minimum_fraction_denominator);

    const std::array<int, 7> spot_r = {6,24,25,49,105,447,448};
    std::array<Row, spot_r.size()> spots;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (std::size_t i = 0; i < spot_r.size(); ++i) {
        spots[i] = evaluate_r(spot_r[i], geometry, profiles, targets);
    }

    std::cout << "METHOD fixed-prefix atom-to-body complement scatter; no zeta\n";
    std::cout << "GEOMETRY DEN " << geometry.denominator
              << " CELLS " << geometry.cells.size()
              << " SUPPORT " << geometry.support_masks.size()
              << " TARGET_MEMBERSHIPS " << target_memberships << '\n';
    std::cout << "BASE BODIES " << profiles.size()
              << " FAILURES " << base_failures
              << " EQUALITIES " << base_equalities
              << " MIN_DELTA27 " << decimal(minimum_delta27)
              << " MIN_LIMIT_DELTA54 " << decimal(2 * minimum_delta27)
              << " MIN_BODY " << labels(minimum_delta_body)
              << " MIN_BASE_MASS " << minimum_delta_mass
              << " MIN_COMPONENTS " << minimum_delta_components
              << " MIN_CUTOFF " << minimum_cutoff
              << " MIN_CUTOFF_BODY " << labels(minimum_cutoff_body)
              << " MIN_CUTOFF_COUNT " << minimum_cutoff_count
              << " MAX_CUTOFF " << maximum_cutoff
              << " MAX_CUTOFF_BODY " << labels(maximum_cutoff_body)
              << " MAX_CUTOFF_COUNT " << maximum_cutoff_count << '\n';
    std::cout << "FINITE R_COUNT " << finite_r.size()
              << " FIRST " << finite_r.front()
              << " LAST " << finite_r.back()
              << " COMPARISONS " << static_cast<u64>(finite_r.size()) * profiles.size()
              << " FAILURES " << total_failures
              << " EQUALITIES " << total_equalities
              << " MIN_NORMALIZED_R " << minimum_normalized_r
              << " MIN_NORMALIZED_DELTA " << decimal(minimum_normalized_delta)
              << " MIN_NORMALIZED_DEN " << decimal(minimum_fraction_denominator)
              << " MIN_NORMALIZED_REDUCED "
              << decimal(minimum_normalized_delta / minimum_fraction_gcd) << '/'
              << decimal(minimum_fraction_denominator / minimum_fraction_gcd)
              << " MIN_NORMALIZED_BODY " << labels(minimum_normalized_body) << '\n';
    for (const Row& row : spots) {
        std::cout << "SPOT R " << row.r
                  << " MIN_DELTA " << decimal(row.minimum_delta)
                  << " MIN_BODY " << labels(row.minimum_body)
                  << " FAILURES " << row.failures
                  << " EQUALITIES " << row.equalities << '\n';
    }
    const Row& hostile_row = spots[0];
    const i128 hostile_denominator =
        static_cast<i128>(2) * hostile_row.r * geometry.denominator;
    const i128 hostile_gcd = exact_gcd(hostile_row.hostile_delta,
                                       hostile_denominator);
    std::cout << "HOSTILE R 6 BODY " << labels(mask_from_values(HOSTILE_BODY_VALUES))
              << " DELTA " << decimal(hostile_row.hostile_delta)
              << " DEN " << decimal(hostile_denominator)
              << " REDUCED " << decimal(hostile_row.hostile_delta / hostile_gcd)
              << '/' << decimal(hostile_denominator / hostile_gcd) << '\n';
}
