#define main composite_primary_main
#include "lrc14_two_outsider_direct_body_cofinal_primary_thm4231.cpp"
#undef main

// Follow-up scout: fix a first newcomer q and exhaust all nine-bodies B.
// This deliberately reuses only the primary's basic exact types and helpers;
// the q-dependent arrangement and coefficient ledger are rebuilt below.

namespace {

struct QCell {
    i64 width;
    u32 failed_pool;
    bool q_safe;
};

struct QBodyRow {
    u64 mass;
    u32 body;
    u32 tail_activation;
    u32 components;
    i64 surplus;
};

std::pair<i64, std::vector<QCell>> build_q_cells(int q) {
    i64 common = 1;
    for (int speed : POOL) common = std::lcm(common, i64{14} * speed);
    common = std::lcm(common, i64{14} * q);

    std::vector<i64> walls = {0, common};
    const auto add_walls = [&](int speed) {
        const i64 unit = common / (14 * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14 * tooth + 1) * unit);
            walls.push_back((14 * tooth + 13) * unit);
        }
    };
    for (int speed : POOL) add_walls(speed);
    add_walls(q);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());

    std::vector<QCell> cells;
    cells.reserve(walls.size() - 1);
    const i64 period = 2 * common;
    for (std::size_t index = 0; index + 1 < walls.size(); ++index) {
        const i64 midpoint_twice = walls[index] + walls[index + 1];
        const auto is_safe = [&](int speed) {
            const i64 residue = static_cast<i64>(
                (i128{speed} * midpoint_twice) % period
            );
            return 7 * residue >= common && 7 * residue <= 13 * common;
        };
        u32 failed = 0;
        for (int vertex = 0; vertex < 30; ++vertex) {
            if (!is_safe(POOL[vertex])) failed |= u32{1} << vertex;
        }
        cells.push_back({
            walls[index + 1] - walls[index], failed, is_safe(q)
        });
    }
    return {common, cells};
}

std::pair<Coefficients, Coefficients> build_q_coefficients(
    const std::vector<QCell>& cells
) {
    std::map<u32, i64> mass;
    std::map<u32, i64> components;
    for (std::size_t index = 0; index < cells.size(); ++index) {
        const QCell& current = cells[index];
        const QCell& previous = cells[(index + cells.size() - 1) % cells.size()];
        if (current.q_safe && std::popcount(current.failed_pool) <= 21) {
            mass[current.failed_pool] += current.width;
            components[current.failed_pool] += 1;
        }
        const u32 joined = previous.failed_pool | current.failed_pool;
        if (
            previous.q_safe && current.q_safe &&
            std::popcount(joined) <= 21
        ) {
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

}  // namespace

#ifndef FIXED_Q_TAIL_NO_MAIN
int main(int argc, char** argv) {
    require(argc == 2, "usage: fixed_q_tail Q");
    const int q = std::stoi(argv[1]);
    require(q > 0, "q must be positive");
    require(std::find(POOL.begin(), POOL.end(), q) == POOL.end(), "q lies in P");

    const auto [common, cells] = build_q_cells(q);
    const auto [mass_coefficients, component_coefficients] =
        build_q_coefficients(cells);

    std::array<std::vector<u32>, BODY_SIZE + 1> right_by_size;
    std::vector<u32> left_bodies;
    for (int size = 0; size <= BODY_SIZE; ++size) {
        right_by_size[size] = fixed_weight_masks(15, size);
        std::vector<u32> layer = fixed_weight_masks(15, size);
        left_bodies.insert(left_bodies.end(), layer.begin(), layer.end());
    }
    std::sort(left_bodies.begin(), left_bodies.end());
    std::vector<std::size_t> offsets(left_bodies.size() + 1, 0);
    for (std::size_t index = 0; index < left_bodies.size(); ++index) {
        const int right_size = BODY_SIZE - std::popcount(left_bodies[index]);
        offsets[index + 1] = offsets[index] + right_by_size[right_size].size();
    }
    require(offsets.back() == 14'307'150, "q-body universe changed");
    std::vector<QBodyRow> rows(offsets.back());

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
            for (const auto& [failure, value] : mass_coefficients) {
                if ((failure & left_body) == 0) mass_zeta[failure >> 15] += value;
            }
            for (const auto& [failure, value] : component_coefficients) {
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
                require(mass > 0 && components > 0, "invalid fixed-q geometry");
                const i64 surplus = 54 * mass - 4 * common;
                u64 activation = std::numeric_limits<u32>::max();
                if (surplus > 0) {
                    activation = ceil_div(
                        u128{54} * components * common,
                        u128{7} * surplus
                    );
                    require(activation < std::numeric_limits<u32>::max(), "fixed-q activation overflow");
                }
                rows[output++] = {
                    static_cast<u64>(mass), body, static_cast<u32>(activation),
                    static_cast<u32>(components), surplus
                };
            }
            require(output == offsets[left_index + 1], "fixed-q block changed");
        }
    }

    std::sort(rows.begin(), rows.end(), [](const QBodyRow& left, const QBodyRow& right) {
        return left.body < right.body;
    });
    u64 nonpositive = 0;
    u64 worst_activation = 0;
    u64 worst_count = 0;
    u32 worst_body = 0;
    i64 minimum_surplus = std::numeric_limits<i64>::max();
    u32 minimum_surplus_body = 0;
    Fnv1a64 fingerprint;
    for (const QBodyRow& row : rows) {
        if (row.surplus <= 0) ++nonpositive;
        if (row.tail_activation > worst_activation) {
            worst_activation = row.tail_activation;
            worst_count = 1;
            worst_body = row.body;
        } else if (row.tail_activation == worst_activation) {
            ++worst_count;
        }
        if (row.surplus < minimum_surplus) {
            minimum_surplus = row.surplus;
            minimum_surplus_body = row.body;
        }
        fingerprint.add(q);
        fingerprint.add(row.body);
        fingerprint.add(row.tail_activation);
        fingerprint.add(row.mass);
        fingerprint.add(row.components);
        fingerprint.add(static_cast<u64>(row.surplus));
    }
    const auto worst = std::lower_bound(
        rows.begin(), rows.end(), worst_body,
        [](const QBodyRow& row, u32 mask) { return row.body < mask; }
    );
    require(worst != rows.end(), "fixed-q worst witness missing");

    std::cout << "LRC14_FIXED_FIRST_OUTSIDER_TAIL_SCOUT\n";
    std::cout << "Q " << q << " COMMON " << common
              << " CELLS " << cells.size()
              << " MASS_COEFF " << mass_coefficients.size()
              << " COMPONENT_COEFF " << component_coefficients.size()
              << " BODIES " << rows.size() << '\n';
    std::cout << "LIMIT_SURPLUS NONPOSITIVE " << nonpositive
              << " MIN_TICKS " << minimum_surplus
              << " WITNESS " << labels(minimum_surplus_body) << '\n';
    std::cout << "TAIL_MAX K " << worst_activation
              << " COUNT " << worst_count
              << " WITNESS " << labels(worst_body)
              << " MASS_TICKS " << worst->mass
              << " COMPONENTS " << worst->components
              << " SURPLUS54 " << worst->surplus << '\n';
    std::cout << "FINGERPRINT " << std::hex << std::setw(16)
              << std::setfill('0') << fingerprint.value << std::dec
              << std::setfill(' ') << '\n';
    std::cout << "VERDICT "
              << (nonpositive == 0 ? "UNIFORM_COFINAL_TAIL" : "LIMIT_OBSTRUCTION")
              << '\n';
    return 0;
}
#endif
