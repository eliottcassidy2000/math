#!/usr/bin/env python3
"""Exact finite three-deletion census for THM-4170.

Fraction arithmetic constructs the fixed-pool walls.  Every wall is embedded
in an audited common integer lattice, and NumPy only vectorizes exact signed
int64 prefix differences.  A separate Python transversal recursion classifies
each resulting 3-uniform hypergraph.  Direct Fraction primitives independently
rebuild the q=924/925 boundary rows.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import lcm

import numpy as np


ANCHORS = (120, 126, 143)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
    264, 286, 290,
)
OPTIONAL = tuple(value for value in POOL if value not in ANCHORS)
POOL_SET = frozenset(POOL)
COMMON = 18_241_159_416_480
LIMIT = 9_699
BLOCK = 128
THRESHOLD = F(4, 63)
FNV_OFFSET = 14_695_981_039_346_656_037
FNV_PRIME = 1_099_511_628_211
UINT64_MASK = (1 << 64) - 1
EXPECTED_HOSTILE = (
    3, 6, 22, 24, 25, 46, 48, 50, 55, 64, 70, 72, 75, 83, 93,
    96, 100, 103, 105, 110, 122, 127, 128, 140, 147, 153, 158, 166,
    172, 173, 183, 186, 192, 206, 210, 220, 256, 260, 270, 282, 294,
    306, 313, 320, 325, 332, 346, 366, 372, 384, 416, 440, 462, 512,
    516, 520, 550, 567, 744, 768, 924,
)
EXPECTED_TAU_HIST = Counter({
    0: 4, 1: 6, 2: 5, 3: 8, 4: 5, 5: 9, 6: 10, 7: 14,
})
CONTROLS = frozenset((1, 7, 200, 924, 925, 4959, 4960, 8266, 9699))


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def safe_at(point, speed):
    residue = (speed * point) % 1
    return F(1, 14) <= residue <= F(13, 14)


def labels_text(mask):
    return tuple(label for index, label in enumerate(OPTIONAL) if mask >> index & 1)


def greedy_matching(edges):
    used = 0
    witness = []
    for edge in edges:
        if edge & used:
            continue
        used |= edge
        witness.append(edge)
    return tuple(witness)


def find_cover(edges, budget):
    failed = set()

    def search(chosen, remaining):
        key = (chosen, remaining)
        if key in failed:
            return None
        uncovered = 0
        matching_used = 0
        matching = 0
        for edge in edges:
            if edge & chosen:
                continue
            if not uncovered:
                uncovered = edge
            if not edge & matching_used:
                matching_used |= edge
                matching += 1
                if matching > remaining:
                    failed.add(key)
                    return None
        if not uncovered:
            return chosen
        if not remaining:
            failed.add(key)
            return None
        branch = uncovered
        while branch:
            bit = branch & -branch
            result = search(chosen | bit, remaining - 1)
            if result is not None:
                return result
            branch ^= bit
        failed.add(key)
        return None

    return search(0, budget)


def minimum_cover_through_seven(edges):
    for budget in range(8):
        witness = find_cover(edges, budget)
        if witness is not None:
            return budget, witness
    return None, None


def safe_prefix(point, speed):
    phase = speed * point
    whole = phase.numerator // phase.denominator
    residue = phase - whole
    partial = max(F(0), min(residue, F(13, 14)) - F(1, 14))
    return (F(6, 7) * whole + partial) / speed


def direct_fraction_row(speed, walls, failure_masks, triples):
    buckets = defaultdict(F)
    previous = safe_prefix(walls[0], speed)
    for index, failure in enumerate(failure_masks):
        current = safe_prefix(walls[index + 1], speed)
        if failure is not None:
            buckets[failure] += current - previous
        previous = current
    measures = []
    for edge in triples:
        measures.append(sum((mass for failure, mass in buckets.items()
                             if failure & ~edge == 0), F(0)))
    return tuple(measures)


def main():
    require(len(OPTIONAL) == 27, "optional size")
    common = 1
    for speed in POOL:
        common = lcm(common, 14 * speed)
    require(common == COMMON, "common lattice")

    walls = {F(0), F(1)}
    for speed in POOL:
        for tooth in range(speed):
            walls.add(F(14 * tooth + 1, 14 * speed))
            walls.add(F(14 * tooth + 13, 14 * speed))
    walls = tuple(sorted(walls))
    require(len(walls) == 7_134, "wall count")
    wall_ticks = np.array([int(point * COMMON) for point in walls], dtype=np.int64)
    require(all(F(int(tick), COMMON) == point
                for tick, point in zip(wall_ticks, walls)), "wall embedding")

    pair_list = tuple(combinations(range(27), 2))
    pair_index = {pair: index for index, pair in enumerate(pair_list)}
    triple_list = tuple(combinations(range(27), 3))
    triple_masks = tuple(sum(1 << vertex for vertex in triple)
                         for triple in triple_list)
    triple_index = {triple: index for index, triple in enumerate(triple_list)}
    require((len(pair_list), len(triple_list)) == (351, 2925), "bank counts")

    group_ids = []
    failure_masks = []
    cell_hist = Counter()
    for left, right in zip(walls, walls[1:]):
        point = (left + right) / 2
        if any(not safe_at(point, anchor) for anchor in ANCHORS):
            group_ids.append(-1)
            failure_masks.append(None)
            cell_hist[4] += 1
            continue
        failures = tuple(index for index, speed in enumerate(OPTIONAL)
                         if not safe_at(point, speed))
        if len(failures) > 3:
            group_ids.append(-1)
            failure_masks.append(None)
            cell_hist[4] += 1
            continue
        mask = sum(1 << vertex for vertex in failures)
        failure_masks.append(mask)
        cell_hist[len(failures)] += 1
        if not failures:
            group_ids.append(0)
        elif len(failures) == 1:
            group_ids.append(1 + failures[0])
        elif len(failures) == 2:
            group_ids.append(28 + pair_index[failures])
        else:
            group_ids.append(379 + triple_index[failures])
    require(cell_hist == Counter({0: 150, 1: 328, 2: 518, 3: 678, 4: 5459}),
            "cell histogram")

    selected_cells = np.flatnonzero(np.array(group_ids, dtype=np.int64) >= 0)
    selected_groups = np.array(group_ids, dtype=np.int64)[selected_cells]
    order = np.argsort(selected_groups, kind="stable")
    ordered_cells = selected_cells[order]
    ordered_groups = selected_groups[order]
    starts = np.concatenate((np.array([0], dtype=np.int64),
                             np.flatnonzero(ordered_groups[1:] != ordered_groups[:-1]) + 1))
    present_groups = ordered_groups[starts]
    first = np.array([row[0] for row in triple_list], dtype=np.int64)
    second = np.array([row[1] for row in triple_list], dtype=np.int64)
    third = np.array([row[2] for row in triple_list], dtype=np.int64)
    pair_first_second = np.array([pair_index[(a, b)] for a, b, _ in triple_list], dtype=np.int64)
    pair_first_third = np.array([pair_index[(a, c)] for a, _, c in triple_list], dtype=np.int64)
    pair_second_third = np.array([pair_index[(b, c)] for _, b, c in triple_list], dtype=np.int64)

    hostile_rows = []
    tau_hist = Counter()
    qualifier_count = 0
    greedy_qualifiers = 0
    branch_qualifiers = 0
    equality_total = 0
    fnv = FNV_OFFSET
    controls = {}
    max_int64 = np.iinfo(np.int64).max

    for q_start in range(1, LIMIT + 1, BLOCK):
        q_stop = min(LIMIT + 1, q_start + BLOCK)
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
        contributions = 12 * whole_difference * COMMON + partial[:, 1:] - partial[:, :-1]
        require(np.all(contributions >= 0), "negative cell contribution")
        reduced = np.add.reduceat(contributions[:, ordered_cells], starts, axis=1)
        grouped = np.zeros((len(q_values), 3304), dtype=np.int64)
        grouped[:, present_groups] = reduced
        numerators = (
            grouped[:, 0, None]
            + grouped[:, 1 + first] + grouped[:, 1 + second] + grouped[:, 1 + third]
            + grouped[:, 28 + pair_first_second]
            + grouped[:, 28 + pair_first_third]
            + grouped[:, 28 + pair_second_third]
            + grouped[:, 379 + np.arange(2925)]
        )
        require(int(numerators.max()) <= max_int64 // 9, "numerator overflow")
        differences = 9 * numerators - 8 * q_values[:, None] * COMMON
        edge_matrix = differences >= 0
        equality_matrix = differences == 0

        for local, q_value in enumerate(q_values):
            q = int(q_value)
            if q in POOL_SET:
                continue
            equality_count = int(equality_matrix[local].sum())
            equality_total += equality_count
            active_indices = np.flatnonzero(edge_matrix[local])
            edges = tuple(triple_masks[int(index)] for index in active_indices)
            matching = greedy_matching(edges)
            cover = None if len(matching) >= 8 else find_cover(edges, 7)
            if cover is None:
                qualifier_count += 1
                if len(matching) >= 8:
                    greedy_qualifiers += 1
                else:
                    branch_qualifiers += 1
                fnv = ((fnv ^ q) * FNV_PRIME) & UINT64_MASK
            else:
                tau, minimum = minimum_cover_through_seven(edges)
                require(tau is not None, "hostile minimum cover")
                hostile_rows.append((q, tau, minimum))
                tau_hist[tau] += 1
                cover = minimum
            if q in CONTROLS:
                controls[q] = (edges, matching, cover, equality_count)

    hostile_q = tuple(row[0] for row in hostile_rows)
    require(hostile_q == EXPECTED_HOSTILE, "hostile q ledger")
    require(tau_hist == EXPECTED_TAU_HIST, "hostile tau histogram")
    require((qualifier_count, greedy_qualifiers, branch_qualifiers) == (9608, 8837, 771),
            "qualifier census")
    require(equality_total == 0, "finite equality audit")
    require(fnv == 0x02784121A66537AC, "qualifier standard FNV")
    require(labels_text(controls[924][2]) == (16, 85, 88, 145, 168, 252),
            "q924 cover")
    require(len(controls[925][1]) == 9 and controls[925][2] is None,
            "q925 matching control")

    direct_rows = {}
    for q in (924, 925):
        measures = direct_fraction_row(q, walls, failure_masks, triple_masks)
        active = tuple(triple_masks[index] for index, measure in enumerate(measures)
                       if measure >= THRESHOLD)
        equalities = tuple(triple_masks[index] for index, measure in enumerate(measures)
                           if measure == THRESHOLD)
        require(active == controls[q][0], "direct Fraction active ledger")
        require(not equalities, "direct Fraction equality")
        if q == 924:
            require(all(edge & controls[q][2] for edge in active), "direct q924 cover")
        witness = controls[q][1][:8]
        witness_rows = tuple((labels_text(edge), measures[triple_masks.index(edge)],
                              measures[triple_masks.index(edge)] - THRESHOLD)
                             for edge in witness)
        direct_rows[q] = witness_rows

    print("LRC14_TRIPLE_DELETION_THM4170_FINITE_CENSUS_20260826")
    print(f"universe=1<=q<9700,q_notin_P;count={LIMIT-len(POOL)}")
    print(f"walls={len(walls)};cells={len(failure_masks)};cell_hist={tuple(sorted(cell_hist.items()))}")
    print(f"qualifiers={qualifier_count};hostile={len(hostile_rows)};last_hostile={hostile_q[-1]};"
          f"greedy_matching8_qualifiers={greedy_qualifiers};branch_only_qualifiers={branch_qualifiers};"
          f"equality_total={equality_total}")
    print(f"hostile_q={hostile_q}")
    print(f"hostile_tau_hist={tuple(sorted(tau_hist.items()))}")
    print("hostile_rows=" + ";".join(
        f"{q}:tau{tau}:{labels_text(cover)}" for q, tau, cover in hostile_rows
    ))
    for q in sorted(CONTROLS):
        edges, matching, cover, equalities = controls[q]
        print(f"control_q={q};edges={len(edges)};greedy={len(matching)};"
              f"cover7={cover is not None};cover={labels_text(cover or 0)};equalities={equalities}")
    for q in (924, 925):
        print(f"direct_fraction_q={q};matching_rows={direct_rows[q]}")
    print(f"qualifier_q_fnv1a64={fnv:016x}")
    print("certificate=all_q>=925_by_finite_925_to_9699_plus_eventual_q>=9700")
    print("PASS")


if __name__ == "__main__":
    main()
