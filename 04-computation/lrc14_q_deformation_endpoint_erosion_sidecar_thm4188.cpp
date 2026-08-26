#pragma GCC diagnostic push
// Renaming the included program's main removes C++'s implicit `return 0` rule;
// the renamed function is never called by this sidecar.
#pragma GCC diagnostic ignored "-Wreturn-type"
#define main lrc14_q_deformation_full_scout_main
#include "lrc14_all_newcomer_zero_original_anchor_hierarchy_thm4188.cpp"
#undef main
#pragma GCC diagnostic pop

// Exact structural sidecar for THM-4188.
//
// For a circular union U of components I_j and the one-runner safe comb G_q,
// interval discrepancy gives
//
//   mu(U cap G_q) >= (6/7) sum_j max(|I_j|-1/(7q),0).
//
// Hence an erosion length at least 2/27 at Q certifies the 4/63 threshold
// for every q>=Q.  The sliding-window failure-mask path below audits all
// q=50 repair edges; a direct component sum independently checks the named
// hostile and minimizing edges.  Q=2479 is sharp only for this erosion
// certificate.  THM-4188's exact inclusion cutoff remains a separate result.
//
// Reproduce (from the repository root):
//   c++ -std=c++20 -O2 -DNDEBUG -Wall -Wextra -pedantic \
//     04-computation/lrc14_q_deformation_endpoint_erosion_sidecar_thm4188.cpp \
//     -o /tmp/lrc14_endpoint_erosion && /tmp/lrc14_endpoint_erosion
// A matching UBSan replay replaces `-O2 -DNDEBUG` by
// `-O1 -g -fsanitize=undefined -fno-sanitize-recover=all`.

namespace {

struct ErosionAtoms {
    i64 scale;
    std::unordered_map<u32, i64> mass;
    i64 ignored_mass;
    std::size_t events;
};

std::size_t containing_cell(const std::vector<i64>& scaled_walls,
                            i64 point_numerator) {
    const auto found = std::upper_bound(scaled_walls.begin(),
                                        scaled_walls.end(),
                                        point_numerator);
    require(found != scaled_walls.begin(), "scaled point precedes first wall");
    const std::size_t index =
        static_cast<std::size_t>(found - scaled_walls.begin() - 1);
    require(index + 1 < scaled_walls.size(), "scaled point reached period end");
    return index;
}

ErosionAtoms build_erosion_atoms(const std::vector<Cell>& cells,
                                 int q_bound,
                                 int max_arity) {
    require(q_bound >= 1, "nonpositive erosion q bound");
    const i64 scale = 7LL * q_bound;
    const i128 period_wide = static_cast<i128>(scale) * COMMON;
    require(period_wide <= std::numeric_limits<i64>::max(),
            "scaled erosion period overflow");
    const i64 period = static_cast<i64>(period_wide);

    std::vector<i64> walls;
    walls.reserve(cells.size() + 1);
    for (const Cell& cell : cells) walls.push_back(cell.left);
    walls.push_back(COMMON);

    std::vector<i64> events;
    events.reserve(2 * walls.size());
    for (i64 wall : walls) {
        const i64 left_event = static_cast<i64>(static_cast<i128>(scale) * wall %
                                                period);
        i64 right_event = left_event - COMMON;
        if (right_event < 0) right_event += period;
        events.push_back(left_event);
        events.push_back(right_event);
    }
    std::sort(events.begin(), events.end());
    events.erase(std::unique(events.begin(), events.end()), events.end());
    require(!events.empty() && events.front() == 0,
            "erosion event partition lost zero");

    const i64 doubled_scale = 2 * scale;
    std::vector<i64> doubled_scaled_walls;
    doubled_scaled_walls.reserve(walls.size());
    for (i64 wall : walls) {
        const i128 value = static_cast<i128>(doubled_scale) * wall;
        require(value <= std::numeric_limits<i64>::max(),
                "doubled scaled wall overflow");
        doubled_scaled_walls.push_back(static_cast<i64>(value));
    }
    const i64 doubled_period = 2 * period;

    ErosionAtoms atoms{scale, {}, 0, events.size()};
    atoms.mass.reserve(4096);
    i128 total = 0;
    for (std::size_t event_index = 0; event_index < events.size(); ++event_index) {
        const i64 left_event = events[event_index];
        const i64 right_event = event_index + 1 < events.size()
                                    ? events[event_index + 1]
                                    : period;
        require(left_event < right_event, "empty erosion event interval");
        const i64 midpoint_numerator = left_event + right_event;
        const i64 right_midpoint_numerator =
            (midpoint_numerator + 2 * COMMON) % doubled_period;
        const std::size_t first =
            containing_cell(doubled_scaled_walls, midpoint_numerator);
        const std::size_t last =
            containing_cell(doubled_scaled_walls, right_midpoint_numerator);

        u32 required = 0;
        std::size_t cell_index = first;
        while (true) {
            required |= cells[cell_index].failed_pool;
            if (cell_index == last) break;
            cell_index = (cell_index + 1) % cells.size();
        }
        const i64 interval_mass = right_event - left_event;
        if (std::popcount(required) <= max_arity) {
            atoms.mass[required] += interval_mass;
        } else {
            atoms.ignored_mass += interval_mass;
        }
        total += interval_mass;
    }
    require(total == period_wide, "erosion atoms do not partition circle");
    return atoms;
}

struct ErosionLayerResult {
    u64 failures = 0;
    i128 minimum_margin = 0;
    u32 minimum_edge = 0;
};

struct ComponentInterval {
    i64 left;
    i64 right;
};

std::vector<ComponentInterval> deletion_components(
    const std::vector<Cell>& cells,
    u32 deletion) {
    std::vector<bool> active(cells.size());
    std::size_t inactive = cells.size();
    for (std::size_t index = 0; index < cells.size(); ++index) {
        active[index] = (cells[index].failed_pool & ~deletion) == 0;
        if (!active[index]) inactive = index;
    }
    require(inactive != cells.size(), "whole circle unexpectedly active");
    std::vector<ComponentInterval> components;
    std::size_t step = 1;
    while (step <= cells.size()) {
        const std::size_t index = (inactive + step) % cells.size();
        if (!active[index]) {
            ++step;
            continue;
        }
        const i64 left = cells[index].left;
        std::size_t end_step = step;
        while (end_step + 1 <= cells.size() &&
               active[(inactive + end_step + 1) % cells.size()]) {
            ++end_step;
        }
        const std::size_t end_index = (inactive + end_step) % cells.size();
        components.push_back({left, cells[end_index].right});
        step = end_step + 1;
    }
    return components;
}

struct EndpointAudit {
    i64 components = 0;
    i64 length = 0;
    i128 surplus = 0;
    i128 endpoint_error = 0;
    i64 residue_support = 0;
    i64 residue_l1 = 0;
    i64 cyclic_boundary_height = 0;
    i128 exact_delta = 0;
};

EndpointAudit audit_endpoints(const std::vector<Cell>& cells,
                              u32 deletion,
                              int q) {
    EndpointAudit audit;
    const std::vector<ComponentInterval> components =
        deletion_components(cells, deletion);
    audit.components = static_cast<i64>(components.size());
    std::unordered_map<i64, i64> oriented_residues;
    oriented_residues.reserve(2 * components.size());
    i128 left_error_sum = 0;
    i128 right_error_sum = 0;
    for (const ComponentInterval& component : components) {
        i64 length = component.right - component.left;
        if (length < 0) length += COMMON;
        audit.length += length;
        const i128 left_error =
            static_cast<i128>(safe_prefix(component.left, q)) -
            static_cast<i128>(12) * q * component.left;
        const i128 right_error =
            static_cast<i128>(safe_prefix(component.right, q)) -
            static_cast<i128>(12) * q * component.right;
        left_error_sum += left_error;
        right_error_sum += right_error;
        audit.endpoint_error += right_error - left_error;
        const i64 left_residue = static_cast<i64>(
            (static_cast<i128>(q) * component.left) % COMMON);
        const i64 right_residue = static_cast<i64>(
            (static_cast<i128>(q) * component.right) % COMMON);
        --oriented_residues[left_residue];
        ++oriented_residues[right_residue];
    }
    require(left_error_sum == -right_error_sum,
            "reflection endpoint-error identity failed");
    audit.surplus = static_cast<i128>(27) * audit.length -
                    static_cast<i128>(2) * COMMON;
    audit.exact_delta = static_cast<i128>(4) * q * audit.surplus +
                        static_cast<i128>(9) * audit.endpoint_error;

    std::vector<std::pair<i64, i64>> ordered(oriented_residues.begin(),
                                              oriented_residues.end());
    std::sort(ordered.begin(), ordered.end());
    i64 cumulative = 0;
    i64 cumulative_minimum = 0;
    i64 cumulative_maximum = 0;
    for (const auto& [residue, coefficient] : ordered) {
        (void)residue;
        if (coefficient != 0) {
            ++audit.residue_support;
            audit.residue_l1 += std::abs(coefficient);
        }
        cumulative += coefficient;
        cumulative_minimum = std::min(cumulative_minimum, cumulative);
        cumulative_maximum = std::max(cumulative_maximum, cumulative);
    }
    require(cumulative == 0, "oriented endpoint mass did not cancel");
    audit.cyclic_boundary_height = cumulative_maximum - cumulative_minimum;

    const AtomMass q_mass = build_atom_mass(cells, q, 7);
    require(audit.exact_delta == delta(deletion, q_mass, q),
            "endpoint formula disagrees with atom integration");
    return audit;
}

ErosionLayerResult audit_erosion_layer(const CompressedLayer& layer,
                                       const std::vector<i64>& dense_mass,
                                       i64 scale) {
    ErosionLayerResult result;
    bool first = true;
    const i128 target = static_cast<i128>(2) * COMMON * scale;
    for (std::size_t edge_index = 0; edge_index < layer.edges.size();
         ++edge_index) {
        i64 eroded_scaled_mass = 0;
        for (u32 position = layer.offsets[edge_index];
             position < layer.offsets[edge_index + 1]; ++position) {
            eroded_scaled_mass += dense_mass[layer.atom_indices[position]];
        }
        const i128 margin = static_cast<i128>(27) * eroded_scaled_mass - target;
        if (first || margin < result.minimum_margin) {
            result.minimum_margin = margin;
            result.minimum_edge = layer.edges[edge_index];
        }
        if (margin < 0) ++result.failures;
        first = false;
    }
    return result;
}

u32 mask_from_values(std::initializer_list<int> values) {
    u32 mask = 0;
    for (int value : values) {
        const auto found = std::find(POOL.begin(), POOL.end(), value);
        require(found != POOL.end(), "requested value is outside pool");
        mask |= u32{1} << std::distance(POOL.begin(), found);
    }
    return mask;
}

i64 direct_eroded_scaled_mass(const std::vector<Cell>& cells,
                              u32 deletion,
                              i64 scale) {
    i128 result = 0;
    for (const ComponentInterval& component :
         deletion_components(cells, deletion)) {
        i64 length = component.right - component.left;
        if (length < 0) length += COMMON;
        const i128 contribution = static_cast<i128>(scale) * length - COMMON;
        if (contribution > 0) result += contribution;
    }
    require(result <= std::numeric_limits<i64>::max(),
            "direct eroded mass overflow");
    return static_cast<i64>(result);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::vector<int> q_bounds;
        for (int index = 1; index < argc; ++index) {
            q_bounds.push_back(std::stoi(argv[index]));
        }
        if (q_bounds.empty()) q_bounds = {2478, 2479};
        const std::vector<Cell> cells = build_pool_cells();
        const AtomMass mass50 = build_atom_mass(cells, 50, 7);
        const std::array<EdgeLayer, 3> layers = {
            build_layer(mass50, 50, 5),
            build_layer(mass50, 50, 6),
            build_layer(mass50, 50, 7)};
        std::vector<ErosionAtoms> atom_tables;
        atom_tables.reserve(q_bounds.size());
        AtomSupport support;
        for (int q_bound : q_bounds) {
            atom_tables.push_back(build_erosion_atoms(cells, q_bound, 7));
            for (const auto& [mask, mass] : atom_tables.back().mass) {
                (void)mass;
                support.masks.push_back(mask);
            }
        }
        std::sort(support.masks.begin(), support.masks.end());
        support.masks.erase(
            std::unique(support.masks.begin(), support.masks.end()),
            support.masks.end());
        support.index.reserve(2 * support.masks.size());
        for (u32 index = 0; index < support.masks.size(); ++index) {
            support.index.emplace(support.masks[index], index);
        }
        const std::array<CompressedLayer, 3> compressed = {
            compress_layer(layers[0], support),
            compress_layer(layers[1], support),
            compress_layer(layers[2], support)};

        const u32 q366_hostile =
            mask_from_values({8, 15, 145, 193, 290});
        const u32 erosion_minimizer =
            mask_from_values({8, 10, 15, 20, 40, 95, 145});
        std::cout << "AUDIT q_deformation_endpoint_erosion_sidecar\n";
        std::cout << "SCOPE VERIFIED_EXACT_EROSION_CERTIFICATE_Q_GE_2479 "
                     "TRUE_INCLUSION_CUTOFF_NOT_CLAIMED LRC14_OPEN\n";
        for (int q : {366, 367}) {
            const EndpointAudit endpoint =
                audit_endpoints(cells, q366_hostile, q);
            std::cout << "ENDPOINT q=" << q
                      << " EDGE " << labels(q366_hostile)
                      << " COMPONENTS " << endpoint.components
                      << " LENGTH " << endpoint.length
                      << " SURPLUS " << decimal(endpoint.surplus)
                      << " ERROR " << decimal(endpoint.endpoint_error)
                      << " RESIDUE_SUPPORT " << endpoint.residue_support
                      << " RESIDUE_L1 " << endpoint.residue_l1
                      << " CYCLIC_HEIGHT " << endpoint.cyclic_boundary_height
                      << " DELTA " << decimal(endpoint.exact_delta) << '\n';
        }
        std::cout << "UNION_ATOM_SUPPORT " << support.masks.size() << '\n';
        for (std::size_t q_index = 0; q_index < q_bounds.size(); ++q_index) {
            const int q_bound = q_bounds[q_index];
            const ErosionAtoms& atoms = atom_tables[q_index];
            const i64 hostile_zeta_mass =
                subset_sum(q366_hostile, atoms.mass);
            const i64 hostile_direct_mass =
                direct_eroded_scaled_mass(cells, q366_hostile, atoms.scale);
            require(hostile_zeta_mass == hostile_direct_mass,
                    "erosion zeta/direct hostile disagreement");
            std::vector<i64> dense_mass(support.masks.size(), 0);
            for (const auto& [mask, mass] : atoms.mass) {
                const auto found = support.index.find(mask);
                require(found != support.index.end(),
                        "erosion introduced an unsupported atom mask");
                dense_mass[found->second] = mass;
            }

            std::cout << "Q_BOUND " << q_bound
                      << " SCALE " << atoms.scale
                      << " EVENTS " << atoms.events
                      << " ATOMS " << atoms.mass.size()
                      << " IGNORED_MASS " << atoms.ignored_mass
                      << " HOSTILE_MASS " << hostile_zeta_mass << '\n';
            ErosionLayerResult d7_result;
            for (std::size_t index = 0; index < layers.size(); ++index) {
                const ErosionLayerResult result = audit_erosion_layer(
                    compressed[index], dense_mass, atoms.scale);
                if (layers[index].arity == 7) d7_result = result;
                std::cout << "D" << layers[index].arity
                          << " EDGES " << layers[index].edges.size()
                          << " FAILURES " << result.failures
                          << " MIN_MARGIN " << decimal(result.minimum_margin)
                          << " MIN_EDGE " << labels(result.minimum_edge) << '\n';
            }
            const i64 minimizer_zeta_mass =
                subset_sum(erosion_minimizer, atoms.mass);
            const i64 minimizer_direct_mass = direct_eroded_scaled_mass(
                cells, erosion_minimizer, atoms.scale);
            require(minimizer_zeta_mass == minimizer_direct_mass,
                    "erosion zeta/direct minimizing-edge disagreement");
            const i128 minimizer_margin =
                static_cast<i128>(27) * minimizer_zeta_mass -
                static_cast<i128>(2) * COMMON * atoms.scale;
            if (q_bound == 2478 || q_bound == 2479) {
                require(d7_result.minimum_edge == erosion_minimizer,
                        "declared erosion minimizer changed");
                require(d7_result.minimum_margin == minimizer_margin,
                        "declared erosion minimum margin changed");
            }
            std::cout << "DIRECT_MIN_EDGE " << labels(erosion_minimizer)
                      << " ZETA_MASS " << minimizer_zeta_mass
                      << " COMPONENT_MASS " << minimizer_direct_mass
                      << " MARGIN " << decimal(minimizer_margin) << '\n';
        }
        std::cout << "VERDICT ENDPOINT_EROSION_SIDECAR_COMPLETE\n";
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
