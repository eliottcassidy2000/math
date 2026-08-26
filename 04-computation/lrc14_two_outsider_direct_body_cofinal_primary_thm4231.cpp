// Exact direct/composite two-outsider certificate for THM-4231.
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
constexpr int DEPTH = 6;
constexpr int BODY_SIZE = 9;
constexpr u64 INF = std::numeric_limits<u64>::max();

struct Cell {
    i64 width;
    u32 failed;
};

struct RepairRow {
    u32 deletion;
    u64 activation;
    i64 mass;
    i64 components;
    i64 surplus;
};

struct BodyRow {
    u64 mass;
    u64 repair_activation;
    u32 body;
    u32 direct_activation;
    u32 components;
    u32 padding;
};
static_assert(sizeof(BodyRow) == 32);

[[noreturn]] void fail(const std::string& label) {
    std::cerr << "FAIL " << label << '\n';
    std::exit(1);
}

void require(bool condition, const std::string& label) {
    if (!condition) fail(label);
}

u64 ceil_div(u128 numerator, u128 denominator) {
    require(denominator > 0, "zero ceiling denominator");
    const u128 answer = (numerator + denominator - 1) / denominator;
    require(answer <= std::numeric_limits<u64>::max(), "ceiling overflow");
    return static_cast<u64>(answer);
}

u32 mask_of(std::initializer_list<int> values) {
    u32 answer = 0;
    for (int value : values) {
        const auto iterator = std::find(POOL.begin(), POOL.end(), value);
        require(iterator != POOL.end(), "label outside pool");
        answer |= u32{1} << std::distance(POOL.begin(), iterator);
    }
    return answer;
}

std::string labels(u32 mask) {
    std::string answer = "{";
    bool first = true;
    for (int vertex = 0; vertex < 30; ++vertex) {
        if ((mask & (u32{1} << vertex)) == 0) continue;
        if (!first) answer += ",";
        answer += std::to_string(POOL[vertex]);
        first = false;
    }
    answer += "}";
    return answer;
}

struct Fnv1a64 {
    u64 value = 0xcbf29ce484222325ULL;
    void add(u64 word) {
        for (int shift = 0; shift < 64; shift += 8) {
            value ^= (word >> shift) & 0xffULL;
            value *= 0x100000001b3ULL;
        }
    }
};

std::pair<i64, std::vector<Cell>> build_cells() {
    i64 common = 1;
    for (int speed : POOL) common = std::lcm(common, i64{14} * speed);
    require(common == 18'241'159'416'480LL, "common denominator changed");

    std::vector<i64> walls = {0, common};
    for (int speed : POOL) {
        const i64 unit = common / (14 * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14 * tooth + 1) * unit);
            walls.push_back((14 * tooth + 13) * unit);
        }
    }
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    require(walls.size() == 7'134, "wall count changed");

    std::vector<Cell> cells;
    cells.reserve(walls.size() - 1);
    const i64 period = 2 * common;
    for (std::size_t index = 0; index + 1 < walls.size(); ++index) {
        const i64 left = walls[index];
        const i64 right = walls[index + 1];
        const i64 midpoint_twice = left + right;
        u32 failed = 0;
        for (int vertex = 0; vertex < 30; ++vertex) {
            const i64 residue = static_cast<i64>(
                (i128{POOL[vertex]} * midpoint_twice) % period
            );
            if (!(7 * residue >= common && 7 * residue <= 13 * common)) {
                failed |= u32{1} << vertex;
            }
        }
        cells.push_back({right - left, failed});
    }
    require(cells.size() == 7'133, "cell count changed");
    return {common, cells};
}

using Coefficients = std::vector<std::pair<u32, i64>>;

std::pair<Coefficients, Coefficients> build_coefficients(
    const std::vector<Cell>& cells,
    int maximum_failure_size
) {
    std::map<u32, i64> mass;
    std::map<u32, i64> components;
    for (std::size_t index = 0; index < cells.size(); ++index) {
        const u32 current = cells[index].failed;
        const u32 previous = cells[(index + cells.size() - 1) % cells.size()].failed;
        if (std::popcount(current) <= maximum_failure_size) {
            mass[current] += cells[index].width;
            components[current] += 1;
        }
        const u32 joined = previous | current;
        if (std::popcount(joined) <= maximum_failure_size) {
            components[joined] -= 1;
        }
    }

    Coefficients mass_rows;
    Coefficients component_rows;
    for (const auto& [mask, value] : mass) {
        if (value != 0) mass_rows.push_back({mask, value});
    }
    for (const auto& [mask, value] : components) {
        if (value != 0) component_rows.push_back({mask, value});
    }
    return {mass_rows, component_rows};
}

i64 subset_sum(u32 mask, const std::map<u32, i64>& table) {
    i64 answer = 0;
    u32 subset = mask;
    while (true) {
        const auto iterator = table.find(subset);
        if (iterator != table.end()) answer += iterator->second;
        if (subset == 0) return answer;
        subset = (subset - 1) & mask;
    }
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
        mask = (((high ^ mask) >> 2) / low) | high;
    }
    return answer;
}

std::vector<RepairRow> build_repairs(
    i64 common,
    const Coefficients& mass_coefficients,
    const Coefficients& component_coefficients
) {
    std::map<u32, i64> mass_table(mass_coefficients.begin(), mass_coefficients.end());
    std::map<u32, i64> component_table(
        component_coefficients.begin(), component_coefficients.end()
    );
    std::vector<RepairRow> rows;
    rows.reserve(140'082);
    i64 nonpositive = 0;
    i64 equalities = 0;
    for (u32 deletion : fixed_weight_masks(30, DEPTH)) {
        const i64 mass = subset_sum(deletion, mass_table);
        const i64 components = subset_sum(deletion, component_table);
        require(mass > 0 && components > 0, "invalid repair geometry");
        const i64 surplus = 45 * mass - 4 * common;
        if (surplus == 0) ++equalities;
        if (surplus <= 0) {
            ++nonpositive;
            continue;
        }
        const u64 activation = ceil_div(
            u128{108} * components * common,
            u128{7} * surplus
        );
        require(
            u128{7} * surplus * activation >= u128{108} * components * common,
            "repair ceiling upper audit"
        );
        require(
            u128{7} * surplus * (activation - 1) <
                u128{108} * components * common,
            "repair ceiling lower audit"
        );
        rows.push_back({deletion, activation, mass, components, surplus});
    }
    require(rows.size() == 140'082, "strict repair count changed");
    require(nonpositive == 453'693 && equalities == 0, "repair sign census changed");

    Fnv1a64 fingerprint;
    for (const RepairRow& row : rows) {
        fingerprint.add(row.deletion);
        fingerprint.add(row.activation);
        fingerprint.add(row.mass);
        fingerprint.add(row.components);
        fingerprint.add(row.surplus);
    }
    require(fingerprint.value == 0xa8b79ad77ad91a62ULL, "repair fingerprint changed");
    std::sort(rows.begin(), rows.end(), [](const RepairRow& left, const RepairRow& right) {
        return std::tie(left.activation, left.deletion) <
               std::tie(right.activation, right.deletion);
    });
    require(rows.front().activation == 3'077, "minimum repair activation changed");
    return rows;
}

u64 first_disjoint_activation(u32 body, const std::vector<RepairRow>& repairs) {
    for (const RepairRow& row : repairs) {
        if ((body & row.deletion) == 0) return row.activation;
    }
    return INF;
}

struct ThresholdCensus {
    u64 cutoff;
    u64 covers = 0;
    u64 direct_good_covers = 0;
    u64 uncertified = 0;
    u64 maximum_cover_activation = 0;
    u32 first_uncertified = 0;
    u32 first_cover = 0;
    Fnv1a64 cover_fingerprint;
};

ThresholdCensus census_at(const std::vector<BodyRow>& rows, u64 cutoff) {
    ThresholdCensus result;
    result.cutoff = cutoff;
    for (const BodyRow& row : rows) {
        if (row.repair_activation <= cutoff) continue;
        ++result.covers;
        result.cover_fingerprint.add(row.body);
        if (result.first_cover == 0) result.first_cover = row.body;
        result.maximum_cover_activation = std::max<u64>(
            result.maximum_cover_activation, row.direct_activation
        );
        if (row.direct_activation <= cutoff) {
            ++result.direct_good_covers;
        } else {
            ++result.uncertified;
            if (result.first_uncertified == 0) result.first_uncertified = row.body;
        }
    }
    return result;
}

void print_census(const ThresholdCensus& census) {
    std::cout << "CENSUS Q " << census.cutoff
              << " COVERS " << census.covers
              << " DIRECT_GOOD " << census.direct_good_covers
              << " UNCERTIFIED " << census.uncertified
              << " MAX_COVER_K " << census.maximum_cover_activation
              << " FIRST_COVER " << labels(census.first_cover)
              << " FIRST_UNCERTIFIED " << labels(census.first_uncertified)
              << " COVER_FNV " << std::hex << std::setw(16) << std::setfill('0')
              << census.cover_fingerprint.value << std::dec << std::setfill(' ') << '\n';
}

}  // namespace

int main() {
    const auto [common, cells] = build_cells();
    const auto [repair_mass_coefficients, repair_component_coefficients] =
        build_coefficients(cells, DEPTH);
    require(repair_mass_coefficients.size() == 2'213, "repair mass coefficient count");
    require(repair_component_coefficients.size() == 1'083, "repair component coefficient count");
    const std::vector<RepairRow> repairs = build_repairs(
        common, repair_mass_coefficients, repair_component_coefficients
    );

    // A nine-body can be disjoint only from failure masks of size at most 21.
    const auto [body_mass_coefficients, body_component_coefficients] =
        build_coefficients(cells, 30 - BODY_SIZE);
    require(body_mass_coefficients.size() == 2'939, "body mass coefficient count");
    require(body_component_coefficients.size() == 1'457, "body component coefficient count");

    std::array<std::vector<u32>, BODY_SIZE + 1> right_by_size;
    std::vector<u32> left_bodies;
    for (int size = 0; size <= BODY_SIZE; ++size) {
        right_by_size[size] = fixed_weight_masks(15, size);
        std::vector<u32> layer = fixed_weight_masks(15, size);
        left_bodies.insert(left_bodies.end(), layer.begin(), layer.end());
    }
    std::sort(left_bodies.begin(), left_bodies.end());
    require(left_bodies.size() == 27'824, "left-body universe changed");

    std::vector<std::size_t> offsets(left_bodies.size() + 1, 0);
    for (std::size_t index = 0; index < left_bodies.size(); ++index) {
        const int right_size = BODY_SIZE - std::popcount(left_bodies[index]);
        offsets[index + 1] = offsets[index] + right_by_size[right_size].size();
    }
    require(offsets.back() == 14'307'150, "body universe changed");
    std::vector<BodyRow> body_rows(offsets.back());

#pragma omp parallel
    {
        std::vector<i64> mass_zeta(u32{1} << 15);
        std::vector<i64> component_zeta(u32{1} << 15);
#pragma omp for schedule(dynamic, 1)
        for (i64 left_index = 0; left_index < static_cast<i64>(left_bodies.size()); ++left_index) {
            const u32 left_body = left_bodies[left_index];
            const int right_size = BODY_SIZE - std::popcount(left_body);
            std::fill(mass_zeta.begin(), mass_zeta.end(), 0);
            std::fill(component_zeta.begin(), component_zeta.end(), 0);
            for (const auto& [failure, value] : body_mass_coefficients) {
                if ((failure & left_body) == 0) mass_zeta[failure >> 15] += value;
            }
            for (const auto& [failure, value] : body_component_coefficients) {
                if ((failure & left_body) == 0) component_zeta[failure >> 15] += value;
            }
            for (int bit = 0; bit < 15; ++bit) {
                const u32 flag = u32{1} << bit;
                for (u32 mask = 0; mask < (u32{1} << 15); ++mask) {
                    if ((mask & flag) == 0) continue;
                    mass_zeta[mask] += mass_zeta[mask ^ flag];
                    component_zeta[mask] += component_zeta[mask ^ flag];
                }
            }

            std::size_t output = offsets[left_index];
            for (u32 right_body : right_by_size[right_size]) {
                const u32 body = left_body | (right_body << 15);
                const u32 query = (~right_body) & 0x7fffU;
                const i64 mass = mass_zeta[query];
                const i64 components = component_zeta[query];
                require(mass > 0 && components > 0, "invalid body geometry");
                const i64 surplus = 45 * mass - 4 * common;
                require(surplus > 0, "nonpositive direct body surplus");
                const u64 direct_activation = ceil_div(
                    u128{108} * components * common,
                    u128{7} * surplus
                );
                require(direct_activation <= std::numeric_limits<u32>::max(), "body activation overflow");
                require(
                    u128{7} * surplus * direct_activation >=
                        u128{108} * components * common,
                    "body ceiling upper audit"
                );
                require(
                    u128{7} * surplus * (direct_activation - 1) <
                        u128{108} * components * common,
                    "body ceiling lower audit"
                );
                const u64 repair_activation = first_disjoint_activation(body, repairs);
                body_rows[output++] = {
                    static_cast<u64>(mass), repair_activation, body,
                    static_cast<u32>(direct_activation), static_cast<u32>(components), 0
                };
            }
            require(output == offsets[left_index + 1], "body block size changed");
        }
    }

    std::sort(body_rows.begin(), body_rows.end(), [](const BodyRow& left, const BodyRow& right) {
        return left.body < right.body;
    });
    require(
        std::adjacent_find(
            body_rows.begin(), body_rows.end(),
            [](const BodyRow& left, const BodyRow& right) { return left.body == right.body; }
        ) == body_rows.end(),
        "duplicate body"
    );

    const u32 witness_w = mask_of({85, 88, 143, 168, 193, 240, 252, 264, 290});
    const auto w_iterator = std::lower_bound(
        body_rows.begin(), body_rows.end(), witness_w,
        [](const BodyRow& row, u32 mask) { return row.body < mask; }
    );
    require(w_iterator != body_rows.end() && w_iterator->body == witness_w, "W missing");
    require(
        w_iterator->mass == 4'802'564'195'362ULL &&
        w_iterator->components == 506 &&
        w_iterator->direct_activation == 995 &&
        w_iterator->repair_activation == 17'548,
        "W geometry or thresholds changed"
    );

    u64 direct_maximum = 0;
    u64 direct_minimum = INF;
    u64 composite_maximum = 0;
    u64 repair_route_maximum = 0;
    u32 direct_witness = 0;
    u32 direct_minimum_witness = 0;
    u32 composite_witness = 0;
    u32 repair_route_witness = 0;
    u64 minimum_mass = INF;
    u64 maximum_mass = 0;
    u32 minimum_components = std::numeric_limits<u32>::max();
    u32 maximum_components = 0;
    std::map<u32, u64> direct_histogram;
    std::map<u64, u64> composite_histogram;
    Fnv1a64 body_fingerprint;
    for (const BodyRow& row : body_rows) {
        const u64 composite = std::min<u64>(row.direct_activation, row.repair_activation);
        if (row.direct_activation > direct_maximum) {
            direct_maximum = row.direct_activation;
            direct_witness = row.body;
        }
        if (row.direct_activation < direct_minimum) {
            direct_minimum = row.direct_activation;
            direct_minimum_witness = row.body;
        }
        if (composite > composite_maximum) {
            composite_maximum = composite;
            composite_witness = row.body;
        }
        if (row.repair_activation > repair_route_maximum) {
            repair_route_maximum = row.repair_activation;
            repair_route_witness = row.body;
        }
        minimum_mass = std::min(minimum_mass, row.mass);
        maximum_mass = std::max(maximum_mass, row.mass);
        minimum_components = std::min(minimum_components, row.components);
        maximum_components = std::max(maximum_components, row.components);
        ++direct_histogram[row.direct_activation];
        ++composite_histogram[composite];
        const i64 surplus = 45 * static_cast<i64>(row.mass) - 4 * common;
        body_fingerprint.add(row.body);
        body_fingerprint.add(row.direct_activation);
        body_fingerprint.add(row.mass);
        body_fingerprint.add(row.components);
        body_fingerprint.add(surplus);
        body_fingerprint.add(row.repair_activation);
        body_fingerprint.add(composite);
    }
    require(composite_maximum <= direct_maximum, "composite exceeds direct maximum");
    require(repair_route_maximum != INF, "body without strict repair");
    const auto direct_second = std::next(direct_histogram.rbegin());
    require(direct_second != direct_histogram.rend(), "direct histogram too short");

    const auto direct_iterator = std::lower_bound(
        body_rows.begin(), body_rows.end(), direct_witness,
        [](const BodyRow& row, u32 mask) { return row.body < mask; }
    );
    const auto composite_iterator = std::lower_bound(
        body_rows.begin(), body_rows.end(), composite_witness,
        [](const BodyRow& row, u32 mask) { return row.body < mask; }
    );
    require(direct_iterator != body_rows.end(), "direct witness missing");
    require(composite_iterator != body_rows.end(), "composite witness missing");

    std::cout << "LRC14_TWO_OUTSIDER_COMPOSITE_DESCENT_PRIMARY_THM4231\n";
    std::cout << "UNIVERSES CELLS " << cells.size()
              << " D6 " << repairs.size()
              << " BODIES " << body_rows.size() << '\n';
    std::cout << "COEFFICIENTS REPAIR_MASS " << repair_mass_coefficients.size()
              << " REPAIR_COMPONENT " << repair_component_coefficients.size()
              << " BODY_MASS " << body_mass_coefficients.size()
              << " BODY_COMPONENT " << body_component_coefficients.size() << '\n';
    std::cout << "DIRECT_MAX Q " << direct_maximum
              << " COUNT " << direct_histogram[static_cast<u32>(direct_maximum)]
              << " WITNESS " << labels(direct_witness)
              << " MASS " << direct_iterator->mass
              << " COMPONENTS " << direct_iterator->components
              << " SURPLUS45 "
              << 45 * static_cast<i64>(direct_iterator->mass) - 4 * common
              << " REPAIR_K " << direct_iterator->repair_activation << '\n';
    std::cout << "DIRECT_SECOND Q " << direct_second->first
              << " COUNT " << direct_second->second << '\n';
    std::cout << "DIRECT_MIN Q " << direct_minimum
              << " COUNT " << direct_histogram[static_cast<u32>(direct_minimum)]
              << " WITNESS " << labels(direct_minimum_witness) << '\n';
    std::cout << "COMPOSITE_MAX Q " << composite_maximum
              << " COUNT " << composite_histogram[composite_maximum]
              << " WITNESS " << labels(composite_witness)
              << " DIRECT_K " << composite_iterator->direct_activation
              << " REPAIR_K " << composite_iterator->repair_activation << '\n';
    std::cout << "REPAIR_ROUTE_MAX Q " << repair_route_maximum
              << " COUNT " << std::count_if(
                     body_rows.begin(), body_rows.end(),
                     [repair_route_maximum](const BodyRow& row) {
                         return row.repair_activation == repair_route_maximum;
                     }
                 )
              << " WITNESS " << labels(repair_route_witness)
              << " MIN_ACTIVE_Q " << repairs.front().activation << '\n';
    std::cout << "BODY_GEOMETRY MASS_MIN " << minimum_mass
              << " MASS_MAX " << maximum_mass
              << " COMPONENTS_MIN " << minimum_components
              << " COMPONENTS_MAX " << maximum_components << '\n';
    std::cout << "W_CONTROL MASS " << w_iterator->mass
              << " COMPONENTS " << w_iterator->components
              << " SURPLUS45 "
              << 45 * static_cast<i64>(w_iterator->mass) - 4 * common
              << " DIRECT_K " << w_iterator->direct_activation
              << " REPAIR_K " << w_iterator->repair_activation << '\n';

    std::vector<u64> cutoffs = {
        composite_maximum - 1,
        composite_maximum,
        direct_maximum - 1,
        direct_maximum,
        3'076,
        3'077,
        17'547,
        17'548,
    };
    std::sort(cutoffs.begin(), cutoffs.end());
    cutoffs.erase(std::unique(cutoffs.begin(), cutoffs.end()), cutoffs.end());
    for (u64 cutoff : cutoffs) print_census(census_at(body_rows, cutoff));

    std::cout << "FINGERPRINT BODY_LEDGER " << std::hex << std::setw(16)
              << std::setfill('0') << body_fingerprint.value << std::dec
              << std::setfill(' ') << '\n';
    std::cout << "VERDICT EXACT_COMPOSITE_THRESHOLD " << composite_maximum
              << " DIRECT_THRESHOLD " << direct_maximum
              << " ACTIVATION_ONLY_TRANSITION 17548 LRC14_OPEN\n";
    return 0;
}
