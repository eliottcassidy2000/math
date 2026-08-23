#!/usr/bin/env python3
"""Exact pair-obstruction widths for the surviving THM-3878 11+2 seams.

This is a scratch research companion.  It builds the 5,855 THM-3818 pair
atlas, computes the largest connected component of the scale-one pair-danger
set and of the odd scale-two two-lift obstruction quotient, and audits the
result by exact wall-cell topology.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")
H = Fraction(1, 14)
GATES = 0


def require(condition: bool, label: object) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def distance(x: Fraction) -> Fraction:
    x %= 1
    return min(x, 1 - x)


def factor_trial(n: int) -> tuple[tuple[int, int], ...]:
    result = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            exponent = 0
            while n % p == 0:
                exponent += 1
                n //= p
            result.append((p, exponent))
        p += 1
    if n > 1:
        result.append((n, 1))
    return tuple(result)


def admissible_sum(n: int) -> bool:
    factors = factor_trial(n)
    return bool(factors) and all(p % 3 == 2 and exponent <= 2 for p, exponent in factors)


def atlas() -> tuple[tuple[int, int], ...]:
    return tuple(
        (p, q)
        for p in range(1, 356)
        for q in range(p + 1, 357 - p)
        if gcd(p, q) == 1 and admissible_sum(p + q)
    )


def pair_bad(p: int, q: int, z: Fraction) -> bool:
    return min(distance(p * z), distance(q * z)) < H


def two_lift_bad(p: int, q: int, z: Fraction) -> bool:
    return pair_bad(p, q, z) and pair_bad(p, q, z + Fraction(1, 2))


def danger_walls(p: int, q: int) -> tuple[Fraction, ...]:
    walls = {Fraction(0)}
    for speed in (p, q):
        for k in range(speed):
            for sign in (-1, 1):
                walls.add(((Fraction(k) + sign * H) / speed) % 1)
    return tuple(sorted(walls))


def doubled_two_lift_walls(p: int, q: int) -> tuple[Fraction, ...]:
    # The branch shift disappears after doubling, so these are exactly the
    # quotient-circle walls w=2z mod 1.
    walls = {Fraction(0)}
    for speed in (p, q):
        for k in range(speed):
            for sign in (-1, 1):
                walls.add((2 * (Fraction(k) + sign * H) / speed) % 1)
    return tuple(sorted(walls))


def max_open_component(
    walls: tuple[Fraction, ...], predicate
) -> tuple[Fraction, int, Fraction]:
    """Return (max component length, component count, total open length).

    Cells are joined across a wall only when the wall itself satisfies the
    open predicate.  This handles coincident boundaries covered by another
    danger interval without spuriously splitting the component.
    """
    n = len(walls)
    require(n > 0 and walls[0] == 0, ("wall normalization", walls[:2]))
    lengths: list[Fraction] = []
    active: list[bool] = []
    for i, left in enumerate(walls):
        right = walls[(i + 1) % n] + (1 if i + 1 == n else 0)
        length = right - left
        require(length > 0, ("positive cell", i, left, right))
        midpoint = (left + length / 2) % 1
        lengths.append(length)
        active.append(bool(predicate(midpoint)))

    parent = list(range(n))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri

    for i, wall in enumerate(walls):
        left_cell = (i - 1) % n
        right_cell = i
        wall_active = bool(predicate(wall))
        if wall_active:
            require(active[left_cell] and active[right_cell],
                    ("open wall neighborhood", wall, left_cell, right_cell))
            union(left_cell, right_cell)

    masses: dict[int, Fraction] = {}
    for i in range(n):
        if active[i]:
            root = find(i)
            masses[root] = masses.get(root, Fraction(0)) + lengths[i]
    if not masses:
        return Fraction(0), 0, Fraction(0)
    total = sum(masses.values(), Fraction(0))
    return max(masses.values()), len(masses), total


def widths(p: int, q: int) -> tuple[Fraction, int, Fraction, Fraction, int, Fraction]:
    beta1, components1, measure1 = max_open_component(
        danger_walls(p, q), lambda z: pair_bad(p, q, z)
    )
    beta2, components2, measure2 = max_open_component(
        doubled_two_lift_walls(p, q),
        lambda w: two_lift_bad(p, q, w / 2),
    )
    return beta1, components1, measure1, beta2, components2, measure2


def fmt(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def top(entries: list[tuple[Fraction, int, int]], count: int = 12):
    return sorted(entries, key=lambda item: (-item[0], item[1], item[2]))[:count]


def main() -> None:
    pairs = atlas()
    require(len(pairs) == 5855, "atlas size")

    table = []
    for p, q in pairs:
        beta1, c1, m1, beta2, c2, m2 = widths(p, q)
        require(Fraction(1, 7) <= m1 <= Fraction(2, 7), ("danger measure", p, q, m1))
        require(Fraction(0) <= m2 <= m1, ("two-lift measure", p, q, m1, m2))
        if p % 2 and q % 2:
            require((beta2 == 0) == (p + q <= 7), ("odd universal criterion", p, q, beta2))
        table.append((p, q, beta1, c1, m1, beta2, c2, m2))

    odd = [row for row in table if row[0] % 2 and row[1] % 2]
    odd_residual = [row for row in odd if row[5] > 0]
    require(len(odd) == 1651, "odd scale-two pre-residual")
    require(len(odd_residual) == 1650, "odd scale-two residual")

    b1_entries = [(row[2], row[0], row[1]) for row in table]
    b2_entries = [(row[5], row[0], row[1]) for row in odd_residual]
    max_b1 = max(x[0] for x in b1_entries)
    max_b2 = max(x[0] for x in b2_entries)
    max_b1_pairs = tuple((p, q) for beta, p, q in b1_entries if beta == max_b1)
    max_b2_pairs = tuple((p, q) for beta, p, q in b2_entries if beta == max_b2)

    # K_s=42 beta_s is the exact coefficient in the necessary scale-ratio
    # inequality t < K_s U furnished by an 11-speed 1/12 witness interval.
    k1 = [(42 * beta, p, q) for beta, p, q in b1_entries]
    k2 = [(42 * beta, p, q) for beta, p, q in b2_entries]
    thresholds = (Fraction(1), Fraction(2), Fraction(3), Fraction(4), Fraction(5), Fraction(6))
    k1_counts = tuple((fmt(a), sum(k <= a for k, _, _ in k1)) for a in thresholds)
    k2_counts = tuple((fmt(a), sum(k <= a for k, _, _ in k2)) for a in thresholds)

    # Positive and hostile topology controls.
    controls = {}
    for pair in ((1, 3), (3, 7), (1, 5), (1, 7), (1, 355), (2, 3), (175, 181)):
        if pair in pairs or pair == (1, 5):
            controls[pair] = widths(*pair)
    require(controls[(1, 3)][3] == 0, "positive two-lift control")
    require(controls[(3, 7)][3] > 0, "minimal two-lift hostile")
    require(controls[(1, 5)][3] == 0, "off-atlas positive control")

    # THM-3743 saturation: every THM-3818 ratio already has its primitive
    # pair Graver relation of norm p+q<=356, independent of both scales.
    require(all(gcd(p, q) == 1 and p + q <= 356 for p, q in pairs),
            "flatness pair relation saturation")
    pair_norm_hist = Counter(p + q for p, q in pairs)

    semantic = {
        "atlas": len(pairs),
        "odd": len(odd),
        "odd_residual": len(odd_residual),
        "max_beta1": fmt(max_b1),
        "max_beta1_pairs": max_b1_pairs,
        "max_beta2": fmt(max_b2),
        "max_beta2_pairs": max_b2_pairs,
        "k1_counts": k1_counts,
        "k2_counts": k2_counts,
        "controls": {
            f"{p},{q}": tuple(fmt(value) if isinstance(value, Fraction) else value for value in data)
            for (p, q), data in controls.items()
        },
        "pair_norm_hist": tuple(sorted(pair_norm_hist.items())),
    }
    digest = sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()

    print("LRC14_REMAINING_CYCLIC_INTERVAL_PROBE_20260823")
    print("scope=THM3818_exact_rank11_two_component_11+2_branch;LRC14=OPEN")
    print(f"atlas={len(pairs)};odd_scale2={len(odd)};odd_residual={len(odd_residual)}")
    print(f"scale1_max_beta={fmt(max_b1)};K1max={fmt(42*max_b1)};pairs={max_b1_pairs}")
    print(f"scale2_max_beta={fmt(max_b2)};K2max={fmt(42*max_b2)};pairs={max_b2_pairs}")
    print(f"scale1_counts_K_le_1_to_6={k1_counts}")
    print(f"scale2_counts_K_le_1_to_6={k2_counts}")
    print("top_scale1=" + repr(tuple((fmt(42*b), p, q) for b, p, q in top(b1_entries))))
    print("top_scale2=" + repr(tuple((fmt(42*b), p, q) for b, p, q in top(b2_entries))))
    for pair, data in sorted(controls.items()):
        b1, c1, m1, b2, c2, m2 = data
        print(
            f"control_{pair[0]}_{pair[1]}="
            f"beta1:{fmt(b1)},components1:{c1},measure1:{fmt(m1)},"
            f"beta2:{fmt(b2)},components2:{c2},measure2:{fmt(m2)}"
        )
    print("flatness_join=all_5855_pairs_already_carry_internal_primitive_Graver_(q,-p)_with_norm_p+q<=356;no_scale_filter")
    print(f"semantic_sha256={digest}")
    print(f"gates={GATES}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
