#!/usr/bin/env python3
"""Exact controls for THM-2798.

The truth-bearing checks use explicit ``require`` calls so ``python -O``
executes the same verification.
"""

from collections import deque
from fractions import Fraction
from itertools import permutations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def defect(m, centre):
    """Number of missing distinct edge differences on K_{1,m}."""
    differences = {abs(j - centre) for j in range(m + 1) if j != centre}
    return m - len(differences)


def predicted_barrier(m, radius):
    return max(0, (m - radius + 1) // 2)


def position_minimax(m, radius):
    """Minimax defect from 0 to m in the radius-r position graph."""
    best = [m + 1] * (m + 1)
    best[0] = defect(m, 0)
    pending = {0}
    while pending:
        c = min(pending, key=lambda x: best[x])
        pending.remove(c)
        for nxt in range(max(0, c - radius), min(m, c + radius) + 1):
            if nxt == c:
                continue
            candidate = max(best[c], defect(m, nxt))
            if candidate < best[nxt]:
                best[nxt] = candidate
                pending.add(nxt)
    return best[m]


def canonical_path(m, radius):
    """A matching path: left bank, one bridge, right bank."""
    barrier = predicted_barrier(m, radius)
    path = [0]
    while path[-1] < barrier:
        path.append(min(barrier, path[-1] + radius))
    right = m - barrier
    if path[-1] != right:
        path.append(right)
    while path[-1] < m:
        path.append(min(m, path[-1] + radius))
    return path


def rotate_centre_right(order, step):
    """Left-rotate the block beginning at the centre by ``step``."""
    centre = order.index(0)
    require(1 <= step < len(order) - centre, "illegal centre rotation")
    return (
        order[:centre]
        + order[centre + 1 : centre + step + 1]
        + (0,)
        + order[centre + step + 1 :]
    )


def realize_position_path(m, positions):
    order = tuple(range(m + 1))
    require(order[0] == 0, "bad initial order")
    for left, right in zip(positions, positions[1:]):
        require(order.index(0) == left, "position/path mismatch")
        order = rotate_centre_right(order, right - left)
    require(order == tuple(range(1, m + 1)) + (0,), "leaf order was not preserved")
    return order


def full_move_neighbours(order, radius):
    """All nontrivial permutations inside one consecutive window of width <=r+1."""
    n = len(order)
    seen = set()
    for width in range(2, min(radius + 1, n) + 1):
        for left in range(n - width + 1):
            block = order[left : left + width]
            for perm in permutations(block):
                if perm == block:
                    continue
                nxt = order[:left] + perm + order[left + width :]
                if nxt not in seen:
                    seen.add(nxt)
                    yield nxt


def full_permutation_minimax(m, radius):
    start = tuple(range(m + 1))
    target = tuple(range(1, m + 1)) + (0,)
    best = {start: 0}
    buckets = [deque() for _ in range(m + 1)]
    buckets[0].append(start)
    for level in range(m + 1):
        while buckets[level]:
            order = buckets[level].popleft()
            if best.get(order) != level:
                continue
            if order == target:
                return level
            for nxt in full_move_neighbours(order, radius):
                nxt_level = max(level, defect(m, nxt.index(0)))
                if nxt_level < best.get(nxt, m + 1):
                    best[nxt] = nxt_level
                    buckets[nxt_level].append(nxt)
    raise RuntimeError("target unreachable")


def main():
    print("THM-2798 STAR GRACEFUL-DEFECT MINIMAX AUDIT")

    rows = []
    checked_pairs = 0
    for m in range(1, 81):
        require(defect(m, 0) == 0 and defect(m, m) == 0, "extreme not graceful")
        for centre in range(m + 1):
            require(defect(m, centre) == min(centre, m - centre), "defect formula")
        for radius in range(1, m + 2):
            expected = predicted_barrier(m, radius)
            actual = position_minimax(m, radius)
            require(actual == expected, "position minimax formula")
            path = canonical_path(m, radius)
            require(path[0] == 0 and path[-1] == m, "bad canonical endpoints")
            require(
                all(1 <= b - a <= radius for a, b in zip(path, path[1:])),
                "canonical move too wide",
            )
            require(max(defect(m, c) for c in path) == expected, "path misses barrier")
            realize_position_path(m, path)
            checked_pairs += 1
        if m <= 12:
            rows.append(
                (
                    m,
                    predicted_barrier(m, 1),
                    predicted_barrier(m, 2),
                    predicted_barrier(m, 3),
                )
            )

    full_cases = 0
    for m in range(1, 7):
        for radius in range(1, min(3, m) + 1):
            require(
                full_permutation_minimax(m, radius) == predicted_barrier(m, radius),
                "full permutation graph disagrees with projection",
            )
            full_cases += 1

    print("formula: defect_m(c)=min(c,m-c)")
    print("barrier: B_r(m)=max(0,ceil((m-r)/2))")
    print("m | C2 radius1 | C2+C3 radius2 | radius3")
    for row in rows:
        print("%2d | %10d | %14d | %7d" % row)
    print(f"position pairs checked: {checked_pairs}")
    print(f"full permutation minimax cases through m=6: {full_cases}")
    print("single-seed bounded-defect barrier is unbounded: PASS")


if __name__ == "__main__":
    main()
