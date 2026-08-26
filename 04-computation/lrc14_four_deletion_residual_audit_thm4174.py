#!/usr/bin/env python3
"""Exact four-deletion scout on the 61 THM-4170 triple-certificate failures."""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import lcm

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp


ANCHORS = (120, 126, 143)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
    264, 286, 290,
)
OPTIONAL = tuple(v for v in POOL if v not in ANCHORS)
BAD3 = (
    3, 6, 22, 24, 25, 46, 48, 50, 55, 64, 70, 72, 75, 83, 93,
    96, 100, 103, 105, 110, 122, 127, 128, 140, 147, 153, 158, 166,
    172, 173, 183, 186, 192, 206, 210, 220, 256, 260, 270, 282, 294,
    306, 313, 320, 325, 332, 346, 366, 372, 384, 416, 440, 462, 512,
    516, 520, 550, 567, 744, 768, 924,
)
COMMON = 18_241_159_416_480


def safe_at(point, speed):
    residue = (speed * point) % 1
    return F(1, 14) <= residue <= F(13, 14)


def labels(mask):
    return tuple(v for i, v in enumerate(OPTIONAL) if mask >> i & 1)


def find_cover_exact(edges, budget):
    """Exact bounded transversal recursion, independent of the MILP path."""
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

    answer = search(0, budget)
    return answer, len(failed)


def main():
    common = 1
    for speed in POOL:
        common = lcm(common, 14 * speed)
    assert common == COMMON
    walls = {F(0), F(1)}
    for speed in POOL:
        for tooth in range(speed):
            walls.add(F(14 * tooth + 1, 14 * speed))
            walls.add(F(14 * tooth + 13, 14 * speed))
    walls = tuple(sorted(walls))
    assert len(walls) == 7134
    wall_ticks = np.array([int(w * COMMON) for w in walls], dtype=np.int64)

    singles = tuple((i,) for i in range(27))
    pairs = tuple(combinations(range(27), 2))
    triples = tuple(combinations(range(27), 3))
    quads = tuple(combinations(range(27), 4))
    pair_idx = {x: i for i, x in enumerate(pairs)}
    triple_idx = {x: i for i, x in enumerate(triples)}
    quad_idx = {x: i for i, x in enumerate(quads)}
    offsets = (0, 1, 28, 379, 3304)
    assert len(quads) == 17550

    group_ids = []
    hist = Counter()
    for left, right in zip(walls, walls[1:]):
        mid = (left + right) / 2
        if any(not safe_at(mid, a) for a in ANCHORS):
            group_ids.append(-1)
            hist[">4/anchor"] += 1
            continue
        failure = tuple(i for i, v in enumerate(OPTIONAL) if not safe_at(mid, v))
        if len(failure) > 4:
            group_ids.append(-1)
            hist[">4/anchor"] += 1
            continue
        hist[len(failure)] += 1
        if len(failure) == 0:
            group_ids.append(0)
        elif len(failure) == 1:
            group_ids.append(1 + failure[0])
        elif len(failure) == 2:
            group_ids.append(28 + pair_idx[failure])
        elif len(failure) == 3:
            group_ids.append(379 + triple_idx[failure])
        else:
            group_ids.append(3304 + quad_idx[failure])

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

    subset_groups = []
    quad_masks = []
    incidence = np.zeros((len(quads), 27), dtype=float)
    for row, quad in enumerate(quads):
        subset_groups.append((
            0,
            *(1 + i for i in quad),
            *(28 + pair_idx[p] for p in combinations(quad, 2)),
            *(379 + triple_idx[t] for t in combinations(quad, 3)),
            3304 + row,
        ))
        mask = sum(1 << i for i in quad)
        quad_masks.append(mask)
        incidence[row, list(quad)] = 1.0
    subset_groups = np.array(subset_groups, dtype=np.int64)
    assert subset_groups.shape == (17550, 16)

    rows = []
    for q in BAD3:
        products = q * wall_ticks
        whole = products // COMMON
        rem = products - whole * COMMON
        scaled = 14 * rem
        partial = 14 * rem - COMMON
        partial[scaled <= COMMON] = 0
        partial[scaled >= 13 * COMMON] = 12 * COMMON
        contributions = 12 * (whole[1:] - whole[:-1]) * COMMON + partial[1:] - partial[:-1]
        assert np.all(contributions >= 0)
        reduced = np.add.reduceat(contributions[ordered_cells], starts)
        grouped = np.zeros(20854, dtype=np.int64)
        grouped[present_groups] = reduced
        numerators = grouped[subset_groups].sum(axis=1)
        differences = 9 * numerators - 8 * q * COMMON
        active = np.flatnonzero(differences >= 0)
        equalities = int(np.count_nonzero(differences == 0))
        assert len(active)

        # Minimum transversal as a 27-variable binary MILP.  The returned
        # integer witness is checked directly against every exact active edge.
        constraint = LinearConstraint(incidence[active], np.ones(len(active)), np.inf)
        result = milp(
            c=np.ones(27),
            integrality=np.ones(27),
            bounds=Bounds(np.zeros(27), np.ones(27)),
            constraints=constraint,
            options={"presolve": True, "time_limit": 30},
        )
        assert result.success, (q, result.message)
        cover = sum(1 << i for i, x in enumerate(result.x) if x > 0.5)
        assert all(quad_masks[i] & cover for i in active)
        tau = cover.bit_count()
        assert abs(result.fun - tau) < 1e-7
        exact_cover7, states = find_cover_exact(
            tuple(quad_masks[i] for i in active), 7
        )
        assert (exact_cover7 is None) == (tau > 7), (q, tau, exact_cover7)
        if exact_cover7 is not None:
            assert exact_cover7.bit_count() <= 7
            assert all(quad_masks[i] & exact_cover7 for i in active)
        rows.append((q, len(active), tau, labels(cover), equalities, states))
        print("ROW", rows[-1], flush=True)

    print("TAU_HIST", tuple(sorted(Counter(row[2] for row in rows).items())))
    print("FAIL_THROUGH_7", tuple(row for row in rows if row[2] <= 7))
    print("RESCUED", tuple(row[0] for row in rows if row[2] > 7))
    print("EQUALITY_TOTAL", sum(row[4] for row in rows))
    print("EXACT_STATE_RANGE", min(row[5] for row in rows), max(row[5] for row in rows))


if __name__ == "__main__":
    main()
