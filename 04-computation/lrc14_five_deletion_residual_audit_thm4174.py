#!/usr/bin/env python3
"""Exact five-deletion scout on the nine four-deletion residual labels."""

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
RESIDUAL4 = (25, 50, 96, 100, 105, 128, 192, 210, 256)
COMMON = 18_241_159_416_480


def safe_at(point, speed):
    residue = (speed * point) % 1
    return F(1, 14) <= residue <= F(13, 14)


def labels(mask):
    return tuple(v for i, v in enumerate(OPTIONAL) if mask >> i & 1)


def find_cover_exact(edges, budget):
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
        if remaining == 0:
            failed.add(key)
            return None
        branch = uncovered
        while branch:
            bit = branch & -branch
            answer = search(chosen | bit, remaining - 1)
            if answer is not None:
                return answer
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

    levels = [tuple(combinations(range(27), d)) for d in range(6)]
    offsets = []
    total = 0
    indices = []
    for level in levels:
        offsets.append(total)
        indices.append({subset: i for i, subset in enumerate(level)})
        total += len(level)
    assert total == 101584
    quints = levels[5]
    assert len(quints) == 80730

    group_ids = []
    hist = Counter()
    for left, right in zip(walls, walls[1:]):
        mid = (left + right) / 2
        if any(not safe_at(mid, a) for a in ANCHORS):
            group_ids.append(-1)
            hist[">5/anchor"] += 1
            continue
        failure = tuple(i for i, v in enumerate(OPTIONAL) if not safe_at(mid, v))
        if len(failure) > 5:
            group_ids.append(-1)
            hist[">5/anchor"] += 1
            continue
        hist[len(failure)] += 1
        group_ids.append(offsets[len(failure)] + indices[len(failure)][failure])

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

    subset_groups = np.empty((len(quints), 32), dtype=np.int64)
    quint_masks = np.empty(len(quints), dtype=np.uint32)
    incidence = np.zeros((len(quints), 27), dtype=float)
    for row, quint in enumerate(quints):
        col = 0
        for d in range(6):
            for subset in combinations(quint, d):
                subset_groups[row, col] = offsets[d] + indices[d][subset]
                col += 1
        assert col == 32
        quint_masks[row] = sum(1 << i for i in quint)
        incidence[row, list(quint)] = 1.0

    rows = []
    for q in RESIDUAL4:
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
        grouped = np.zeros(total, dtype=np.int64)
        grouped[present_groups] = reduced
        numerators = grouped[subset_groups].sum(axis=1)
        differences = 9 * numerators - 8 * q * COMMON
        active = np.flatnonzero(differences >= 0)
        equalities = int(np.count_nonzero(differences == 0))
        assert len(active)

        result = milp(
            c=np.ones(27),
            integrality=np.ones(27),
            bounds=Bounds(np.zeros(27), np.ones(27)),
            constraints=LinearConstraint(
                incidence[active], np.ones(len(active)), np.inf
            ),
            options={"presolve": True, "time_limit": 60},
        )
        assert result.success, (q, result.message)
        cover = sum(1 << i for i, x in enumerate(result.x) if x > 0.5)
        active_masks = tuple(int(quint_masks[i]) for i in active)
        assert all(edge & cover for edge in active_masks)
        tau = cover.bit_count()
        assert abs(result.fun - tau) < 1e-7
        cover7, states = find_cover_exact(active_masks, 7)
        assert (cover7 is None) == (tau > 7), (q, tau, cover7)
        if cover7 is not None:
            assert all(edge & cover7 for edge in active_masks)
        row = (q, len(active), tau, labels(cover), equalities, states)
        rows.append(row)
        print("ROW", row, flush=True)

    print("CELL_HIST", tuple(sorted(hist.items(), key=lambda x: str(x[0]))))
    print("TAU_HIST", tuple(sorted(Counter(row[2] for row in rows).items())))
    print("FAIL_THROUGH_7", tuple(row for row in rows if row[2] <= 7))
    print("RESCUED", tuple(row[0] for row in rows if row[2] > 7))
    print("EQUALITY_TOTAL", sum(row[4] for row in rows))
    print("EXACT_STATE_RANGE", min(row[5] for row in rows), max(row[5] for row in rows))


if __name__ == "__main__":
    main()
