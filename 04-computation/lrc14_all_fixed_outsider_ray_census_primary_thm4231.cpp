#define FIXED_Q_TAIL_NO_MAIN
#include "lrc14_fixed_one_outsider_cofinal_tail_primary_thm4231.cpp"

namespace {

struct RayResult {
    int q = 0;
    i64 common = 0;
    u64 cells = 0;
    u64 mass_coefficients = 0;
    u64 component_coefficients = 0;
    u64 bodies = 0;
    u64 nonpositive = 0;
    i64 minimum_surplus = std::numeric_limits<i64>::max();
    u32 minimum_surplus_body = 0;
    u64 worst_activation = 0;
    u64 worst_count = 0;
    u32 worst_body = 0;
    u64 worst_mass = 0;
    u32 worst_components = 0;
    i64 worst_surplus = 0;
    u64 fingerprint = 0;
};

RayResult scan_ray(
    int q,
    const std::vector<u32>& left_bodies,
    const std::array<std::vector<u32>, BODY_SIZE + 1>& right_by_size
) {
    const auto [common, cells] = build_q_cells(q);
    const auto [mass_coefficients, component_coefficients] =
        build_q_coefficients(cells);

    RayResult result;
    result.q = q;
    result.common = common;
    result.cells = cells.size();
    result.mass_coefficients = mass_coefficients.size();
    result.component_coefficients = component_coefficients.size();
    Fnv1a64 fingerprint;

    std::vector<i64> mass_zeta(u32{1} << 15);
    std::vector<i64> component_zeta(u32{1} << 15);
    Coefficients filtered_mass;
    Coefficients filtered_components;
    filtered_mass.reserve(mass_coefficients.size());
    filtered_components.reserve(component_coefficients.size());

    const auto accept = [&](u32 body, i64 mass, i64 components) {
        require(mass > 0 && components > 0, "invalid ray geometry");
        const i64 surplus = 54 * mass - 4 * common;
        u64 activation = INF;
        if (surplus > 0) {
            activation = ceil_div(
                u128{54} * components * common,
                u128{7} * surplus
            );
            require(
                u128{7} * surplus * activation >=
                    u128{54} * components * common,
                "ray ceiling upper audit"
            );
            require(
                u128{7} * surplus * (activation - 1) <
                    u128{54} * components * common,
                "ray ceiling lower audit"
            );
        } else {
            ++result.nonpositive;
        }
        ++result.bodies;
        if (surplus < result.minimum_surplus) {
            result.minimum_surplus = surplus;
            result.minimum_surplus_body = body;
        }
        if (activation > result.worst_activation) {
            result.worst_activation = activation;
            result.worst_count = 1;
            result.worst_body = body;
            result.worst_mass = mass;
            result.worst_components = components;
            result.worst_surplus = surplus;
        } else if (activation == result.worst_activation) {
            ++result.worst_count;
        }
        fingerprint.add(q);
        fingerprint.add(body);
        fingerprint.add(activation);
        fingerprint.add(mass);
        fingerprint.add(components);
        fingerprint.add(static_cast<u64>(surplus));
    };

    for (u32 left_body : left_bodies) {
        const int left_size = std::popcount(left_body);
        const int right_size = BODY_SIZE - left_size;
        if (left_size <= 5) {
            std::fill(mass_zeta.begin(), mass_zeta.end(), 0);
            std::fill(component_zeta.begin(), component_zeta.end(), 0);
            for (const auto& [failure, value] : mass_coefficients) {
                if ((failure & left_body) == 0) mass_zeta[failure >> 15] += value;
            }
            for (const auto& [failure, value] : component_coefficients) {
                if ((failure & left_body) == 0) component_zeta[failure >> 15] += value;
            }
            for (int bit = 0; bit < 15; ++bit) {
                const u32 flag = u32{1} << bit;
                for (u32 base = 0; base < (u32{1} << 15); base += 2 * flag) {
                    for (u32 offset = 0; offset < flag; ++offset) {
                        const u32 source = base + offset;
                        const u32 target = source + flag;
                        mass_zeta[target] += mass_zeta[source];
                        component_zeta[target] += component_zeta[source];
                    }
                }
            }
            for (u32 right_body : right_by_size[right_size]) {
                const u32 query = (~right_body) & 0x7fffU;
                accept(
                    left_body | (right_body << 15),
                    mass_zeta[query], component_zeta[query]
                );
            }
        } else {
            filtered_mass.clear();
            filtered_components.clear();
            for (const auto& [failure, value] : mass_coefficients) {
                if ((failure & left_body) == 0) {
                    filtered_mass.push_back({failure >> 15, value});
                }
            }
            for (const auto& [failure, value] : component_coefficients) {
                if ((failure & left_body) == 0) {
                    filtered_components.push_back({failure >> 15, value});
                }
            }
            for (u32 right_body : right_by_size[right_size]) {
                i64 mass = 0;
                i64 components = 0;
                for (const auto& [failure, value] : filtered_mass) {
                    if ((failure & right_body) == 0) mass += value;
                }
                for (const auto& [failure, value] : filtered_components) {
                    if ((failure & right_body) == 0) components += value;
                }
                accept(left_body | (right_body << 15), mass, components);
            }
        }
    }
    require(result.bodies == 14'307'150, "ray body count changed");
    result.fingerprint = fingerprint.value;
    return result;
}

bool rational_less(i64 a_num, i64 a_den, i64 b_num, i64 b_den) {
    return i128{a_num} * b_den < i128{b_num} * a_den;
}

void print_ray(const RayResult& row) {
    std::cout << "RAY Q " << row.q
              << " D " << row.common
              << " CELLS " << row.cells
              << " COEFF " << row.mass_coefficients << "," << row.component_coefficients
              << " NONPOSITIVE " << row.nonpositive
              << " MIN_SURPLUS " << row.minimum_surplus
              << " MIN_BODY " << labels(row.minimum_surplus_body)
              << " WORST_K ";
    if (row.worst_activation == INF) std::cout << "INF";
    else std::cout << row.worst_activation;
    std::cout << " WORST_COUNT " << row.worst_count
              << " WORST_BODY " << labels(row.worst_body)
              << " WORST_MASS " << row.worst_mass
              << " WORST_C " << row.worst_components
              << " WORST_SURPLUS " << row.worst_surplus
              << " FNV " << std::hex << std::setw(16) << std::setfill('0')
              << row.fingerprint << std::dec << std::setfill(' ') << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    const int lower = argc >= 2 ? std::stoi(argv[1]) : 1;
    const int upper = argc >= 3 ? std::stoi(argv[2]) : 1'306;
    require(1 <= lower && lower <= upper, "invalid q range");

    std::array<std::vector<u32>, BODY_SIZE + 1> right_by_size;
    std::vector<u32> left_bodies;
    for (int size = 0; size <= BODY_SIZE; ++size) {
        right_by_size[size] = fixed_weight_masks(15, size);
        std::vector<u32> layer = fixed_weight_masks(15, size);
        left_bodies.insert(left_bodies.end(), layer.begin(), layer.end());
    }
    std::sort(left_bodies.begin(), left_bodies.end());

    std::vector<int> q_values;
    for (int q = lower; q <= upper; ++q) {
        if (std::find(POOL.begin(), POOL.end(), q) == POOL.end()) q_values.push_back(q);
    }
    std::vector<RayResult> rows(q_values.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (i64 index = 0; index < static_cast<i64>(q_values.size()); ++index) {
        rows[index] = scan_ray(q_values[index], left_bodies, right_by_size);
    }

    u64 positive_rays = 0;
    u64 nonpositive_rays = 0;
    u64 nonpositive_body_cases = 0;
    int first_positive = 0;
    int first_hostile = 0;
    u32 first_hostile_body = 0;
    u64 global_worst_activation = 0;
    int global_worst_q = 0;
    u32 global_worst_body = 0;
    std::size_t minimum_margin_index = 0;
    Fnv1a64 range_fingerprint;
    for (std::size_t index = 0; index < rows.size(); ++index) {
        const RayResult& row = rows[index];
        if (row.nonpositive == 0) {
            ++positive_rays;
            if (first_positive == 0) first_positive = row.q;
            if (row.worst_activation > global_worst_activation) {
                global_worst_activation = row.worst_activation;
                global_worst_q = row.q;
                global_worst_body = row.worst_body;
            }
        } else {
            ++nonpositive_rays;
            if (first_hostile == 0) {
                first_hostile = row.q;
                first_hostile_body = row.minimum_surplus_body;
            }
        }
        nonpositive_body_cases += row.nonpositive;
        if (
            index == 0 || rational_less(
                row.minimum_surplus, row.common,
                rows[minimum_margin_index].minimum_surplus,
                rows[minimum_margin_index].common
            )
        ) {
            minimum_margin_index = index;
        }
        range_fingerprint.add(row.q);
        range_fingerprint.add(row.common);
        range_fingerprint.add(row.nonpositive);
        range_fingerprint.add(static_cast<u64>(row.minimum_surplus));
        range_fingerprint.add(row.minimum_surplus_body);
        range_fingerprint.add(row.worst_activation);
        range_fingerprint.add(row.worst_body);
        range_fingerprint.add(row.fingerprint);
    }

    std::cout << "LRC14_FIXED_FIRST_OUTSIDER_RANGE_SCAN\n";
    std::cout << "RANGE " << lower << " " << upper
              << " RAYS " << rows.size()
              << " BODY_CASES " << u64{rows.size()} * 14'307'150ULL << '\n';
    for (const RayResult& row : rows) {
        if (row.q == lower || row.q == upper || row.q == 1 || row.q == 2 || row.nonpositive != 0) {
            print_ray(row);
        }
    }
    std::vector<std::size_t> by_tail(rows.size());
    std::iota(by_tail.begin(), by_tail.end(), 0);
    std::sort(by_tail.begin(), by_tail.end(), [&](std::size_t left, std::size_t right) {
        return std::tie(rows[right].worst_activation, rows[left].q) <
               std::tie(rows[left].worst_activation, rows[right].q);
    });
    std::cout << "TOP_TAIL_RAYS";
    for (std::size_t rank = 0; rank < std::min<std::size_t>(10, by_tail.size()); ++rank) {
        const RayResult& row = rows[by_tail[rank]];
        std::cout << " " << row.q << ":";
        if (row.worst_activation == INF) std::cout << "INF";
        else std::cout << row.worst_activation;
    }
    std::cout << '\n';

    std::vector<std::size_t> by_margin(rows.size());
    std::iota(by_margin.begin(), by_margin.end(), 0);
    std::sort(by_margin.begin(), by_margin.end(), [&](std::size_t left, std::size_t right) {
        if (rational_less(
                rows[left].minimum_surplus, rows[left].common,
                rows[right].minimum_surplus, rows[right].common
            )) return true;
        if (rational_less(
                rows[right].minimum_surplus, rows[right].common,
                rows[left].minimum_surplus, rows[left].common
            )) return false;
        return rows[left].q < rows[right].q;
    });
    std::cout << "SMALLEST_MARGIN_RAYS";
    for (std::size_t rank = 0; rank < std::min<std::size_t>(10, by_margin.size()); ++rank) {
        const RayResult& row = rows[by_margin[rank]];
        std::cout << " " << row.q << ":" << row.minimum_surplus << "/" << row.common;
    }
    std::cout << '\n';
    const RayResult& minimum_margin = rows[minimum_margin_index];
    std::cout << "SUMMARY POSITIVE_RAYS " << positive_rays
              << " NONPOSITIVE_RAYS " << nonpositive_rays
              << " NONPOSITIVE_BODY_CASES " << nonpositive_body_cases
              << " FIRST_POSITIVE " << first_positive
              << " FIRST_HOSTILE " << first_hostile
              << " FIRST_HOSTILE_BODY " << labels(first_hostile_body) << '\n';
    std::cout << "GLOBAL_MIN_MARGIN Q " << minimum_margin.q
              << " TICKS " << minimum_margin.minimum_surplus
              << " D " << minimum_margin.common
              << " BODY " << labels(minimum_margin.minimum_surplus_body) << '\n';
    std::cout << "GLOBAL_WORST_TAIL K " << global_worst_activation
              << " Q " << global_worst_q
              << " BODY " << labels(global_worst_body) << '\n';
    std::cout << "RANGE_FNV " << std::hex << std::setw(16) << std::setfill('0')
              << range_fingerprint.value << std::dec << std::setfill(' ') << '\n';
    std::cout << "VERDICT "
              << (nonpositive_rays == 0 ? "ALL_RAYS_POSITIVE" : "HOSTILE_RAYS_PRESENT")
              << '\n';
    return 0;
}
