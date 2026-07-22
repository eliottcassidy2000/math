#!/usr/bin/env python3
"""Exact referee for THM-2081's relative Hunter rank-seven gate.

All three-frequency edge weights are computed by complete rational-boundary
atomization. Runtime checks remain active under ``python -O``.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from math import gcd


MAX_CORE_SPEED = 24
FULL_DIVISOR_MASK = (1 << 13) - 1
DIVISOR_MASK = {
    q: sum(1 << (d - 2) for d in range(2, 15) if q % d == 0)
    for q in range(1, MAX_CORE_SPEED + 1)
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_distance(x: F) -> F:
    residue = x % 1
    return min(residue, 1 - residue)


def divisor_complete(Q: tuple[int, ...]) -> bool:
    mask = 0
    for q in Q:
        mask |= DIVISOR_MASK[q]
    return mask == FULL_DIVISOR_MASK


def hereditarily_primitive(Q: tuple[int, ...]) -> bool:
    return all(gcd(*(Q[:i] + Q[i + 1 :])) == 1 for i in range(len(Q)))


def mixed_overlap(h: int, q: int) -> F:
    g = gcd(h, q)
    a, b = h // g, q // g
    x, y = F(a % 14, 14), F(b % 7, 7)
    correction = min(x, y) + max(F(0), x + y - 1) - 2 * x * y
    return F(2, 49) + F(2, a * b) * correction


def restricted_pair(h: int, p: int, q: int, cache: dict[tuple[int, int, int], F]) -> F:
    key = (h, min(p, q), max(p, q))
    if key in cache:
        return cache[key]

    boundaries = {F(0), F(1)}
    for speed, radius in ((h, F(1, 7)), (p, F(1, 14)), (q, F(1, 14))):
        for k in range(speed + 1):
            center = F(k, speed)
            for sign in (-1, 1):
                point = center + sign * radius / speed
                if 0 <= point <= 1:
                    boundaries.add(point)

    points = sorted(boundaries)
    measure = F(0)
    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        if (
            circle_distance(h * midpoint) >= F(1, 7)
            and circle_distance(p * midpoint) < F(1, 14)
            and circle_distance(q * midpoint) < F(1, 14)
        ):
            measure += right - left
    cache[key] = measure
    return measure


def maximum_spanning_tree(
    Q: tuple[int, ...], h: int, cache: dict[tuple[int, int, int], F]
) -> tuple[F, tuple[tuple[F, int, int], ...]]:
    used = {0}
    total = F(0)
    edges: list[tuple[F, int, int]] = []
    while len(used) < len(Q):
        best = None
        for i in used:
            for j in range(len(Q)):
                if j in used:
                    continue
                weight = restricted_pair(h, Q[i], Q[j], cache)
                candidate = (weight, i, j)
                if best is None or candidate > best:
                    best = candidate
        require(best is not None, "Prim search lost a vertex")
        total += best[0]
        edges.append(best)
        used.add(best[2])
    return total, tuple(edges)


def main() -> None:
    cache: dict[tuple[int, int, int], F] = {}
    core_count = 0
    pair_count = 0
    scalar_survivors = 0
    relative_survivors = 0
    worst = None

    for Q in combinations(range(1, MAX_CORE_SPEED + 1), 7):
        if not divisor_complete(Q) or not hereditarily_primitive(Q):
            continue
        core_count += 1
        # THM-2077 at terminal rank seven: 6h < 16 max(Q).
        for h in range(1, (8 * max(Q) - 1) // 3 + 1, 2):
            if 6 * h >= 16 * max(Q):
                continue
            pair_count += 1
            mixed_sum = sum((mixed_overlap(h, q) for q in Q), F(0))
            deficit = F(2, 7) - mixed_sum
            if deficit < 0:
                continue
            scalar_survivors += 1
            tree, edges = maximum_spanning_tree(Q, h, cache)
            margin = tree - deficit
            if margin <= 0:
                relative_survivors += 1
            row = (margin, Q, h, mixed_sum, deficit, tree, edges)
            if worst is None or row < worst:
                worst = row

    require(core_count == 131, f"core count changed: {core_count}")
    require(pair_count == 4120, f"pair count changed: {pair_count}")
    require(relative_survivors == 0, "relative Hunter left a bounded survivor")
    require(worst is not None, "no scalar survivor was audited")
    require(worst[0] == F(561797, 8288280), f"worst margin changed: {worst[0]}")
    require(worst[1] == (1, 9, 10, 11, 13, 14, 24), "worst core changed")
    require(worst[2] == 23, "worst guard changed")

    print("THM-2081 RELATIVE HUNTER RANK-SEVEN REFEREE")
    print(f"hereditary divisor-complete cores through 24: {core_count}")
    print(f"allowed odd core/guard pairs: {pair_count}")
    print(f"mixed-star scalar survivors: {scalar_survivors}")
    print(f"relative-tree survivors: {relative_survivors}")
    print(f"three-frequency weights cached: {len(cache)}")
    print(f"minimum relative-tree margin: {worst[0]}")
    print(f"hostile core/guard: {worst[1]}, {worst[2]}")
    print(f"hostile mixed sum: {worst[3]}")
    print(f"hostile deficit: {worst[4]}")
    print(f"hostile tree weight: {worst[5]}")
    print(f"hostile tree edges: {worst[6]}")
    print("PASS")


if __name__ == "__main__":
    main()
