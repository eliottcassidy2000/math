#!/usr/bin/env python3
"""Exact global census and analytic cutoff for THM-4166.

Fraction constructs the one global wall arrangement of the fixed THM-4156
pool.  Every wall lies on a proved common integer lattice.  NumPy is used only
to vectorize signed int64 arithmetic on that lattice; all overflow bounds and
the exact threshold numerator are checked explicitly.  This is independent of
the primary script's sequential leave-two Fraction intersections.
"""

from fractions import Fraction as F
from itertools import combinations
from math import comb, floor, gcd, lcm

import numpy as np


THRESHOLD = F(4, 63)
DENSITY = F(6, 7)
OSCILLATION = F(6, 49)
ANCHORS = (120, 126, 143)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
    264, 286, 290,
)
OPTIONAL = tuple(value for value in POOL if value not in ANCHORS)
GLOBAL_LIMIT = 49_493
POST_CUTOFF = 49_494
COMMON = 18_241_159_416_480
BLOCK = 128
UINT64_MASK = (1 << 64) - 1


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def safe_at(point, speed):
    residue = (speed * point) % 1
    return F(1, 14) <= residue <= F(13, 14)


def maximum_independent_set(adjacency):
    order = len(adjacency)
    universe = (1 << order) - 1
    complement = tuple(
        universe & ~(1 << vertex) & ~adjacency[vertex]
        for vertex in range(order)
    )
    best_size = 0
    best_mask = 0

    def search(chosen, candidates, excluded):
        nonlocal best_size, best_mask
        if chosen.bit_count() + candidates.bit_count() <= best_size:
            return
        if not candidates and not excluded:
            size = chosen.bit_count()
            if size > best_size:
                best_size = size
                best_mask = chosen
            return
        pivot_pool = candidates | excluded
        if pivot_pool:
            pivot = max(
                (v for v in range(order) if pivot_pool >> v & 1),
                key=lambda v: (candidates & complement[v]).bit_count(),
            )
            extensions = candidates & ~complement[pivot]
        else:
            extensions = candidates
        while extensions:
            bit = extensions & -extensions
            vertex = bit.bit_length() - 1
            search(
                chosen | bit,
                candidates & complement[vertex],
                excluded & complement[vertex],
            )
            candidates &= ~bit
            excluded |= bit
            extensions &= ~bit
            if chosen.bit_count() + candidates.bit_count() <= best_size:
                return

    search(0, universe, 0)
    return best_size, best_mask


def graph_from_pairs(pairs):
    adjacency = [0] * 27
    for first, second in pairs:
        adjacency[first] |= 1 << second
        adjacency[second] |= 1 << first
    return tuple(adjacency)


def fnv_update(state, value):
    return ((state ^ int(value)) * 1_099_511_628_211) & UINT64_MASK


def format_fraction(value):
    return f"{value.numerator}/{value.denominator}"


def main():
    require(len(OPTIONAL) == 27, "optional size changed")
    computed_common = 1
    for speed in POOL:
        computed_common = lcm(computed_common, 14 * speed)
    require(computed_common == COMMON, "common wall lattice changed")

    walls = {F(0), F(1)}
    for speed in POOL:
        for tooth in range(speed):
            walls.add(F(14 * tooth + 1, 14 * speed))
            walls.add(F(14 * tooth + 13, 14 * speed))
    walls = tuple(sorted(walls))
    require(len(walls) == 7_134, "global wall count changed")
    wall_ticks = np.array(
        [int(point * COMMON) for point in walls], dtype=np.int64
    )
    require(
        all(F(int(tick), COMMON) == point for tick, point in zip(wall_ticks, walls)),
        "wall did not embed exactly in common lattice",
    )

    pair_list = tuple(combinations(range(27), 2))
    pair_index = {pair: index for index, pair in enumerate(pair_list)}
    require(len(pair_list) == 351, "pair count changed")
    group_ids = []
    failure_masks = []
    class_histogram = [0, 0, 0, 0]
    for left, right in zip(walls, walls[1:]):
        point = (left + right) / 2
        if any(not safe_at(point, anchor) for anchor in ANCHORS):
            group_ids.append(-1)
            failure_masks.append(None)
            class_histogram[3] += 1
            continue
        failures = tuple(
            index for index, speed in enumerate(OPTIONAL)
            if not safe_at(point, speed)
        )
        if len(failures) == 0:
            group_ids.append(0)
            failure_masks.append(0)
            class_histogram[0] += 1
        elif len(failures) == 1:
            group_ids.append(1 + failures[0])
            failure_masks.append(1 << failures[0])
            class_histogram[1] += 1
        elif len(failures) == 2:
            first, second = failures
            group_ids.append(28 + pair_index[(first, second)])
            failure_masks.append((1 << first) | (1 << second))
            class_histogram[2] += 1
        else:
            group_ids.append(-1)
            failure_masks.append(None)
            class_histogram[3] += 1
    require(class_histogram == [150, 328, 518, 6_137], "cell classes changed")

    cell_ticks = np.diff(wall_ticks)
    require(np.all(cell_ticks > 0), "nonpositive wall cell")

    # Exact base masses, component counts, and the analytic supergraph cutoff.
    base_mass_ticks = {}
    component_counts = {}
    stable_pairs = []
    transient_rows = []
    for first, second in pair_list:
        allowed = (1 << first) | (1 << second)
        selected = np.array(
            [mask is not None and mask & ~allowed == 0 for mask in failure_masks],
            dtype=np.bool_,
        )
        mass_ticks = int(cell_ticks[selected].sum(dtype=np.int64))
        starts = selected & np.concatenate((np.array([True]), ~selected[:-1]))
        components = int(starts.sum())
        base_mass_ticks[(first, second)] = mass_ticks
        component_counts[(first, second)] = components
        if 27 * mass_ticks >= 2 * COMMON:
            require(27 * mass_ticks != 2 * COMMON, "limiting equality appeared")
            stable_pairs.append((first, second))
        else:
            deficit = THRESHOLD - DENSITY * F(mass_ticks, COMMON)
            bound = OSCILLATION * components / deficit
            transient_rows.append((bound, first, second))
    require(len(stable_pairs) == 39, "stable edge count changed")
    require(len(transient_rows) == 312, "transient edge count changed")
    stable_adjacency = graph_from_pairs(stable_pairs)
    stable_alpha, stable_witness = maximum_independent_set(stable_adjacency)
    require((stable_alpha, 27 - stable_alpha) == (21, 6), "stable graph changed")

    def supergraph_at(q):
        pairs = list(stable_pairs)
        pairs.extend(
            (first, second) for bound, first, second in transient_rows
            if floor(bound) >= q
        )
        adjacency = graph_from_pairs(pairs)
        alpha, witness = maximum_independent_set(adjacency)
        return tuple(pairs), alpha, 27 - alpha, witness

    super_493 = supergraph_at(GLOBAL_LIMIT)
    super_494 = supergraph_at(POST_CUTOFF)
    require(
        (len(super_493[0]), super_493[1], super_493[2]) == (53, 19, 8),
        "q=49493 analytic supergraph changed",
    )
    require(
        (len(super_494[0]), super_494[1], super_494[2]) == (52, 20, 7),
        "q=49494 analytic supergraph changed",
    )
    transition = next(
        row for row in transient_rows
        if (OPTIONAL[row[1]], OPTIONAL[row[2]]) == (170, 240)
    )
    require(
        transition[0] == F(2_997_437_002_560, 60_562_157),
        "cutoff transition bound changed",
    )
    require(floor(transition[0]) == GLOBAL_LIMIT, "cutoff floor changed")

    # Arrange every relevant cell by its base/single/pair accumulation group.
    selected_cells = np.flatnonzero(np.array(group_ids, dtype=np.int64) >= 0)
    selected_groups = np.array(group_ids, dtype=np.int64)[selected_cells]
    order = np.argsort(selected_groups, kind="stable")
    ordered_cells = selected_cells[order]
    ordered_groups = selected_groups[order]
    starts = np.concatenate((
        np.array([0], dtype=np.int64),
        np.flatnonzero(ordered_groups[1:] != ordered_groups[:-1]) + 1,
    ))
    present_groups = ordered_groups[starts]
    require(
        len(set(int(value) for value in present_groups)) == len(present_groups),
        "duplicate accumulation group",
    )
    require(0 in present_groups, "base accumulation group missing")

    adjacency = np.zeros((POST_CUTOFF + 1, 27), dtype=np.uint32)
    edge_counts = np.zeros(POST_CUTOFF + 1, dtype=np.int16)
    equality_counts = np.zeros(POST_CUTOFF + 1, dtype=np.int16)
    first_indices = np.array([pair[0] for pair in pair_list], dtype=np.int64)
    second_indices = np.array([pair[1] for pair in pair_list], dtype=np.int64)
    max_int64 = np.iinfo(np.int64).max

    for q_start in range(1, POST_CUTOFF + 1, BLOCK):
        q_stop = min(POST_CUTOFF + 1, q_start + BLOCK)
        q_values = np.arange(q_start, q_stop, dtype=np.int64)
        products = q_values[:, None] * wall_ticks[None, :]
        require(np.all(products >= 0), "q-wall product overflow")
        whole = products // COMMON
        remainder = products - whole * COMMON
        scaled = 14 * remainder
        partial = 14 * remainder - COMMON
        partial[scaled <= COMMON] = 0
        partial[scaled >= 13 * COMMON] = 12 * COMMON
        whole_difference = whole[:, 1:] - whole[:, :-1]
        require(
            int(np.max(np.abs(whole_difference))) * 12 * COMMON < max_int64,
            "cell whole-part overflow",
        )
        contributions = (
            12 * whole_difference * COMMON + partial[:, 1:] - partial[:, :-1]
        )
        require(np.all(contributions >= 0), "negative safe-cell contribution")
        reduced = np.add.reduceat(
            contributions[:, ordered_cells], starts, axis=1
        )
        grouped = np.zeros((len(q_values), 379), dtype=np.int64)
        grouped[:, present_groups] = reduced
        require(grouped.shape == (len(q_values), 379), "grouped shape changed")
        numerators = (
            grouped[:, 0, None]
            + grouped[:, 1 + first_indices]
            + grouped[:, 1 + second_indices]
            + grouped[:, 28 + np.arange(351)]
        )
        require(int(numerators.max()) <= max_int64 // 9, "9*numerator overflow")
        rhs = 8 * q_values[:, None] * COMMON
        require(np.all(rhs >= 0), "threshold rhs overflow")
        differences = 9 * numerators - rhs
        edge_matrix = differences >= 0
        equality_matrix = differences == 0
        edge_counts[q_start:q_stop] = edge_matrix.sum(axis=1, dtype=np.int16)
        equality_counts[q_start:q_stop] = equality_matrix.sum(axis=1, dtype=np.int16)
        for pair_number, (first, second) in enumerate(pair_list):
            local = edge_matrix[:, pair_number]
            adjacency[q_start:q_stop, first] |= (
                local.astype(np.uint32) << np.uint32(second)
            )
            adjacency[q_start:q_stop, second] |= (
                local.astype(np.uint32) << np.uint32(first)
            )

    labels = tuple(q for q in range(1, GLOBAL_LIMIT + 1) if q not in POOL)
    require(len(labels) == 49_463, "global q universe changed")
    require(sum(int(equality_counts[q]) for q in labels) == 0, "finite equality")

    cache = {}
    rows = []
    tau_histogram = {}
    semantic_hash = 14_695_981_039_346_656_037
    for q in labels:
        graph = tuple(int(mask) for mask in adjacency[q])
        if graph not in cache:
            cache[graph] = maximum_independent_set(graph)
        alpha, witness = cache[graph]
        tau = 27 - alpha
        row = (q, int(edge_counts[q]), alpha, tau, witness, graph)
        rows.append(row)
        tau_histogram[tau] = tau_histogram.get(tau, 0) + 1
        semantic_hash = fnv_update(semantic_hash, q)
        semantic_hash = fnv_update(semantic_hash, row[1])
        semantic_hash = fnv_update(semantic_hash, alpha)
        for mask in graph:
            semantic_hash = fnv_update(semantic_hash, mask)

    admitted = tuple(row for row in rows if row[3] > 7)
    expected_histogram = {
        0: 45, 1: 127, 2: 124, 3: 596, 4: 793, 5: 6241, 6: 38003,
        7: 2502, 8: 377, 9: 435, 10: 81, 11: 45, 12: 32, 13: 19,
        14: 11, 15: 13, 16: 4, 17: 5, 18: 3, 19: 5, 20: 2,
    }
    require(tau_histogram == expected_histogram, "global tau histogram changed")
    require(len(admitted) == 1_032, "global admitted count changed")
    require(admitted[-1][0] == 8_265, "last admitted q changed")
    require(semantic_hash == 0x995AA971AF1069E4, "semantic hash changed")
    maximum_tau_rows = tuple((row[0], row[3]) for row in rows if row[3] == 20)
    require(maximum_tau_rows == ((380, 20), (386, 20)), "maximum tau rows changed")
    one_deletion_ten = (5, 66, 182, 298, 336, 340, 380, 386, 528, 572)
    admitted_set = {row[0] for row in admitted}
    require(set(one_deletion_ten) <= admitted_set, "one-deletion ten not subsumed")

    def row_at(q):
        graph = tuple(int(mask) for mask in adjacency[q])
        alpha, witness = maximum_independent_set(graph)
        return (q, int(edge_counts[q]), alpha, 27 - alpha, witness)

    boundary_positive = row_at(8_265)
    boundary_hostile = row_at(8_266)
    cutoff_actual = row_at(GLOBAL_LIMIT)
    post_cutoff_actual = row_at(POST_CUTOFF)
    require(boundary_positive[3] == 8, "last positive is not admitted")
    require(boundary_hostile[3] <= 7, "immediate hostile unexpectedly admitted")
    require(cutoff_actual[3] <= 7 and post_cutoff_actual[3] <= 7,
            "cutoff controls changed")

    extension_bodies = len(admitted) * comb(27, 7)
    old_bodies = comb(27, 8)
    require(extension_bodies == 916_446_960, "extension count changed")
    require(extension_bodies + old_bodies == 918_667_035, "total count changed")

    print("THM4166_EXACT_FRACTION_LATTICE_GLOBAL_CENSUS")
    print("walls", len(walls), "cells", len(walls) - 1, "common", COMMON)
    print("cell_class_histogram", tuple(class_histogram))
    print("q_universe", f"1..{GLOBAL_LIMIT} outside P", len(labels))
    print("threshold", format_fraction(THRESHOLD), "comparison", ">=")
    print("threshold_equalities", sum(int(equality_counts[q]) for q in labels))
    print("stable_edges", len(stable_pairs), "stable_alpha", stable_alpha,
          "stable_tau", 27 - stable_alpha)
    print("stable_edge_labels", tuple(
        (OPTIONAL[first], OPTIONAL[second]) for first, second in stable_pairs
    ))
    print("stable_independent_witness", tuple(
        OPTIONAL[index] for index in range(27) if stable_witness >> index & 1
    ))
    print("transition_bound_170_240", format_fraction(transition[0]),
          "floor", floor(transition[0]))
    print("analytic_supergraph_q49493", len(super_493[0]), super_493[1], super_493[2])
    print("analytic_supergraph_q49494", len(super_494[0]), super_494[1], super_494[2])
    print("analytic_supergraph_q49493_edge_labels", tuple(
        (OPTIONAL[first], OPTIONAL[second]) for first, second in super_493[0]
    ))
    print("analytic_supergraph_q49494_edge_labels", tuple(
        (OPTIONAL[first], OPTIONAL[second]) for first, second in super_494[0]
    ))
    print("analytic_supergraph_q49494_cover", tuple(
        OPTIONAL[index] for index in range(27) if not (super_494[3] >> index & 1)
    ))
    print("tau_histogram", tuple(sorted(tau_histogram.items())))
    print("admitted_count", len(admitted), "last_admitted", admitted[-1][0])
    print("admitted_q", tuple(row[0] for row in admitted))
    print("admitted_rows_q_edges_alpha_tau", tuple(row[:4] for row in admitted))
    print("maximum_tau_rows", maximum_tau_rows)
    print("boundary_q8265", boundary_positive[:4])
    print("hostile_q8266", boundary_hostile[:4])
    print("actual_q49493", cutoff_actual[:4])
    print("actual_q49494", post_cutoff_actual[:4])
    print("one_deletion_ten_subsumed", set(one_deletion_ten) <= admitted_set)
    print("bodies_per_q", comb(27, 7))
    print("extension_bodies", extension_bodies)
    print("with_old_thm4156", extension_bodies + old_bodies)
    print("semantic_fnv1a64", f"{semantic_hash:016x}")
    print("unique_graphs", len(cache))
    print("GLOBAL_CENSUS_PASS")


if __name__ == "__main__":
    main()
