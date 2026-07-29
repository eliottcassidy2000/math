#!/usr/bin/env python3
"""Exact audit for THM-2928's literal six-body drift reductions.

All arithmetic is rational.  The phase-window masses are computed twice:
first by explicit interval construction and then by an independent
breakpoint/midpoint integration.  The fixed two-drift constants and the
body-uniform endpoint ruler are also checked exactly.
"""

from fractions import Fraction
from math import gcd, lcm


Q = Fraction
TARGET = Q(61, 273)
FIVE_SAFE = Q(478, 1365)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def ceil_fraction(x: Fraction) -> int:
    return -((-x.numerator) // x.denominator)


def merge(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    merged: list[tuple[Fraction, Fraction]] = []
    for left, right in sorted(intervals):
        if left >= right:
            continue
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    return merged


def phase_intervals(c: int, d: int, residue: int) -> list[tuple[Fraction, Fraction]]:
    """Return {u in [0,1]: ||c(u+residue)/d||<1/14}, up to endpoints."""

    intervals: list[tuple[Fraction, Fraction]] = []
    # The affine argument lies in [0,c].  The two spare integers on either
    # side make the endpoint handling manifestly harmless.
    for integer in range(-2, c + 3):
        left = Q(d, c) * (integer - Q(1, 14)) - residue
        right = Q(d, c) * (integer + Q(1, 14)) - residue
        left = max(left, Q(0))
        right = min(right, Q(1))
        if left < right:
            intervals.append((left, right))
    return merge(intervals)


def intersect(
    first: list[tuple[Fraction, Fraction]],
    second: list[tuple[Fraction, Fraction]],
) -> list[tuple[Fraction, Fraction]]:
    answer: list[tuple[Fraction, Fraction]] = []
    i = 0
    j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            answer.append((left, right))
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return answer


def interval_mass(intervals: list[tuple[Fraction, Fraction]]) -> Fraction:
    return sum((right - left for left, right in intervals), Q(0))


def distance_to_integer(x: Fraction) -> Fraction:
    floor_x = x.numerator // x.denominator
    fractional = x - floor_x
    return min(fractional, 1 - fractional)


def breakpoint_mass(c: int, d: int, residues: tuple[int, ...]) -> Fraction:
    """Independent cell-decomposition integration of all phase predicates."""

    breakpoints = {Q(0), Q(1)}
    for residue in residues:
        for integer in range(-2, c + 3):
            for sign in (-1, 1):
                point = Q(d, c) * (integer + sign * Q(1, 14)) - residue
                if 0 < point < 1:
                    breakpoints.add(point)
    ordered = sorted(breakpoints)
    total = Q(0)
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        if all(
            distance_to_integer(Q(c, d) * (midpoint + residue)) < Q(1, 14)
            for residue in residues
        ):
            total += right - left
    return total


def simultaneous_mass(c: int, d: int, residues: tuple[int, ...]) -> Fraction:
    intervals = [(Q(0), Q(1))]
    for residue in residues:
        intervals = intersect(intervals, phase_intervals(c, d, residue))
    direct = interval_mass(intervals)
    independent = breakpoint_mass(c, d, residues)
    require(direct == independent, f"phase mass mismatch for {(c, d, residues)}")
    return direct


def main() -> None:
    safe_floors = {
        1: Q(6, 7),
        2: Q(66, 91),
        3: Q(55, 91),
        4: Q(558, 1183),
        5: FIVE_SAFE,
        6: TARGET,
    }
    phase_tariffs = tuple(
        Q(49, 6) * (safe_floors[k] - Q(7 - k, 7)) for k in range(1, 7)
    )
    require(
        phase_tariffs == (Q(0), Q(7, 78), Q(7, 26), Q(119, 338), Q(308, 585), Q(77, 117)),
        "phase-load tariff table",
    )
    require(Q(2, 1) / phase_tariffs[4] == Q(585, 154), "two-drift first slope")
    require(Q(585, 154) < 5, "first drift lies below the largest aligned speed")
    require(Q(585, 154) + Q(8, 7) == Q(5327, 1078) < 5, "four-tooth cutoff")
    require(Q(91, 7) == 13, "second-clock BV coefficient")
    require(14 * lcm(*range(1, 15)) == 5_045_040, "universal six-body ruler")

    # From g>L(h-1/7), h>=61/273, and d=L/g.
    require(TARGET - Q(1, 7) == Q(22, 273), "safe-density difference")
    require(Q(273, 22) < 13, "coindex cutoff")
    coindices = tuple(range(1, 13))

    # Since c*u/nu>0, an occupied d-clock has at most ceil(d/7)
    # residue classes.  Its density is therefore at most ceil(d/7)/d.
    density_candidates = tuple(
        d for d in coindices if Q(ceil_fraction(Q(d, 7)), d) >= TARGET
    )
    require(density_candidates == (1, 2, 3, 4, 8), "density coindices")

    # Reflection sends r to -1-r mod d.  The d=2 and d=4 candidates require
    # one occupied class but have no reflection-fixed class.
    fixed_residues = {
        d: tuple(r for r in range(d) if (2 * r + 1) % d == 0)
        for d in (1, 2, 3, 4)
    }
    require(fixed_residues == {1: (0,), 2: (), 3: (1,), 4: ()}, "reflection")

    # THM-1094 on an interval of length 1/d gives
    # u_A <= 1/7+6d/(49c).  With u_A>=61/273 this is
    # c<=819d/539.
    require(Q(6 * 273, 49 * 22) == Q(819, 539), "phase cutoff constant")
    c_candidates = {
        d: tuple(c for c in range(1, (819 * d) // 539 + 1) if gcd(c, d) == 1)
        for d in (1, 3, 8)
    }
    require(
        c_candidates == {1: (1,), 3: (1, 2, 4), 8: (1, 3, 5, 7, 9, 11)},
        "phase numerator candidates",
    )

    expected_d1 = {1: Q(1, 7)}
    expected_d3 = {1: Q(0), 2: Q(3, 14), 4: Q(3, 28)}
    expected_d8 = {
        1: (Q(1, 7), Q(0), Q(0), Q(0)),
        3: (Q(0), Q(0), Q(1, 21), Q(0)),
        5: (Q(0), Q(1, 35), Q(0), Q(0)),
        7: (Q(0), Q(0), Q(0), Q(1, 49)),
        9: (Q(2, 63), Q(0), Q(0), Q(1, 63)),
        11: (Q(0), Q(1, 77), Q(2, 77), Q(0)),
    }

    actual_d1 = {c: simultaneous_mass(c, 1, (0,)) for c in c_candidates[1]}
    actual_d3 = {c: simultaneous_mass(c, 3, (1,)) for c in c_candidates[3]}
    actual_d8 = {
        c: tuple(simultaneous_mass(c, 8, (r, 7 - r)) for r in range(4))
        for c in c_candidates[8]
    }
    require(actual_d1 == expected_d1, "d=1 phase table")
    require(actual_d3 == expected_d3, "d=3 phase table")
    require(actual_d8 == expected_d8, "d=8 phase table")

    every_mass = list(actual_d1.values()) + list(actual_d3.values())
    every_mass += [mass for row in actual_d8.values() for mass in row]
    require(max(every_mass) == Q(3, 14), "global phase maximum")
    require(max(every_mass) < TARGET, "strict phase contradiction")

    # With only five aligned combs, the inner survivor has the stronger mass
    # 478/1365.  If one of two drifts covered that survivor by itself on
    # every selected body cell, the same outer density/reflection reduction
    # would leave only the following clocks.
    require(FIVE_SAFE - Q(1, 7) == Q(283, 1365), "five-safe excess")
    require(Q(6 * 1365, 49 * 283) == Q(8190, 13867), "five-safe phase cutoff")
    five_c_candidates = {
        d: tuple(c for c in range(1, (8190 * d) // 13867 + 1) if gcd(c, d) == 1)
        for d in (1, 3, 8)
    }
    require(
        five_c_candidates == {1: (), 3: (1,), 8: (1, 3)},
        "five-aligned numerator candidates",
    )
    five_d3 = {1: simultaneous_mass(1, 3, (1,))}
    five_d8 = {
        c: tuple(simultaneous_mass(c, 8, (r, 7 - r)) for r in range(4))
        for c in five_c_candidates[8]
    }
    require(five_d3 == {1: Q(0)}, "five-aligned d=3 table")
    require(
        five_d8
        == {
            1: (Q(1, 7), Q(0), Q(0), Q(0)),
            3: (Q(0), Q(0), Q(1, 21), Q(0)),
        },
        "five-aligned d=8 table",
    )
    five_phase_max = max([*five_d3.values(), *(x for row in five_d8.values() for x in row)])
    require(five_phase_max == Q(1, 7), "five-aligned maximum")
    require(five_phase_max < FIVE_SAFE, "five-aligned strict phase contradiction")

    print("THM-2928 literal one-drift clock audit")
    print(f"safe target: {TARGET}")
    print(f"phase-load tariffs k=1..6: {phase_tariffs}")
    print("two-drift finite box: c1/d1<585/154, c2/d2<13(max A)")
    print("direct six-body universal ruler: L divides 5,045,040")
    print("coindex cutoff: d<=12")
    print(f"density candidates: {density_candidates}")
    print("reflection survivors: (1, 3, 8)")
    print(f"phase numerators: {c_candidates}")
    print(f"d=1 masses: {actual_d1}")
    print(f"d=3 masses: {actual_d3}")
    print(f"d=8 reflected-pair masses: {actual_d8}")
    print(f"largest candidate phase mass: {max(every_mass)} < {TARGET}")
    print(f"five-aligned safe target: {FIVE_SAFE}")
    print(f"five-aligned one-drift numerators: {five_c_candidates}")
    print(f"five-aligned d=3 masses: {five_d3}")
    print(f"five-aligned d=8 reflected-pair masses: {five_d8}")
    print(f"largest five-aligned phase mass: {five_phase_max} < {FIVE_SAFE}")
    print("PASS: no one drift covers the five-aligned residual on the full body carrier")
    print("PASS: the literal six-body/six-aligned-plus-one-drift branch is empty")


if __name__ == "__main__":
    main()
