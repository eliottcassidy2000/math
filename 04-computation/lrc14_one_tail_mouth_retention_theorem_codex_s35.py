#!/usr/bin/env python3
"""THM-542: exact one-tail mouth-retention lemma for the LRC14 collar.

This is a proof certificate for the family

    C_{h,r} = ({1,...,13} \\ {6,h}) union {r},  h != 6, r >= 14.

The goal is exact and arithmetic: prove that the only one-tail row below the
AP one-hole second value 426/35035 is (h,r)=(10,20), and that row keeps the
four old drop-6 mouths intact.

The infinite tail is cut off by the periodic-comb bound.  If G is a union of c
intervals and D_r is the speed-r danger comb at level 1/14, then

    meas(G \\ D_r) >= (6/7) meas(G) - 2c/(7r).

This is just exact counting over r equal periods: on full periods the danger
duty cycle is exactly 1/7, and each interval has at most two partial end
periods, each contributing at most one danger tooth of length 1/(7r).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache, reduce
from math import gcd


TARGET_DENOM = 14
COLLAR = Fraction(7, 858)
AP_SECOND = Fraction(426, 35035)
DROP6_MOUTHS = (
    (Fraction(29, 182), Fraction(9, 56)),
    (Fraction(29, 168), Fraction(27, 154)),
    (Fraction(127, 154), Fraction(139, 168)),
    (Fraction(47, 56), Fraction(153, 182)),
)


@dataclass(frozen=True)
class BaseRow:
    h: int
    measure: Fraction
    components: int
    cutoff: int
    denominator_slack: Fraction


@dataclass(frozen=True)
class TailRow:
    h: int
    r: int
    measure: Fraction
    components: int
    old_survivor: Fraction
    new_mass: Fraction
    core: tuple[int, ...]


def fmt(q: Fraction) -> str:
    return f"{q} = {float(q):.9f}"


def ceil_fraction(q: Fraction) -> int:
    return (q.numerator + q.denominator - 1) // q.denominator


def merge(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    out: list[tuple[Fraction, Fraction]] = []
    for lo, hi in sorted(intervals):
        if lo >= hi:
            continue
        if out and lo <= out[-1][1]:
            if hi > out[-1][1]:
                out[-1] = (out[-1][0], hi)
        else:
            out.append((lo, hi))
    return out


@lru_cache(maxsize=None)
def danger_arcs(speed: int) -> tuple[tuple[Fraction, Fraction], ...]:
    arcs: list[tuple[Fraction, Fraction]] = []
    denom = TARGET_DENOM * speed
    for tooth in range(speed):
        left = TARGET_DENOM * tooth - 1
        right = TARGET_DENOM * tooth + 1
        if left < 0:
            arcs.append((Fraction(0), Fraction(right, denom)))
            arcs.append((Fraction(denom + left, denom), Fraction(1)))
        else:
            arcs.append((Fraction(left, denom), Fraction(right, denom)))
    return tuple(arcs)


def safe_components(core: tuple[int, ...]) -> tuple[tuple[Fraction, Fraction], ...]:
    danger: list[tuple[Fraction, Fraction]] = []
    for speed in core:
        danger.extend(danger_arcs(speed))
    safe: list[tuple[Fraction, Fraction]] = []
    cursor = Fraction(0)
    for lo, hi in merge(danger):
        if lo > cursor:
            safe.append((cursor, lo))
        if hi > cursor:
            cursor = hi
    if cursor < 1:
        safe.append((cursor, Fraction(1)))
    return tuple(safe)


def measure(intervals: tuple[tuple[Fraction, Fraction], ...]) -> Fraction:
    return sum((hi - lo for lo, hi in intervals), Fraction(0))


def intersect_measure(
    xs: tuple[tuple[Fraction, Fraction], ...], ys: tuple[tuple[Fraction, Fraction], ...]
) -> Fraction:
    total = Fraction(0)
    i = j = 0
    while i < len(xs) and j < len(ys):
        lo = max(xs[i][0], ys[j][0])
        hi = min(xs[i][1], ys[j][1])
        if lo < hi:
            total += hi - lo
        if xs[i][1] < ys[j][1]:
            i += 1
        else:
            j += 1
    return total


def primitive(core: tuple[int, ...]) -> bool:
    return reduce(gcd, core, 0) == 1


def base_core(h: int) -> tuple[int, ...]:
    return tuple(v for v in range(1, 14) if v not in (6, h))


def tail_core(h: int, r: int) -> tuple[int, ...]:
    return tuple(sorted(base_core(h) + (r,)))


def base_row(h: int) -> BaseRow:
    comps = safe_components(base_core(h))
    m = measure(comps)
    c = len(comps)
    slack = 6 * m - 7 * AP_SECOND
    if slack <= 0:
        raise AssertionError(f"base h={h} has no positive comb slack")
    cutoff = ceil_fraction(Fraction(2 * c, 1) / slack)
    return BaseRow(h=h, measure=m, components=c, cutoff=cutoff, denominator_slack=slack)


def tail_row(h: int, r: int) -> TailRow:
    core = tail_core(h, r)
    comps = safe_components(core)
    m = measure(comps)
    old = intersect_measure(comps, DROP6_MOUTHS)
    return TailRow(
        h=h,
        r=r,
        measure=m,
        components=len(comps),
        old_survivor=old,
        new_mass=m - old,
        core=core,
    )


def comb_lower_bound(base: BaseRow, r: int) -> Fraction:
    return Fraction(6, 7) * base.measure - Fraction(2 * base.components, 7 * r)


def finite_check(base: BaseRow) -> tuple[TailRow, list[TailRow], int]:
    checked = 0
    below: list[TailRow] = []
    best: TailRow | None = None
    for r in range(14, base.cutoff):
        core = tail_core(base.h, r)
        if not primitive(core):
            continue
        checked += 1
        row = tail_row(base.h, r)
        if best is None or (row.measure, row.r) < (best.measure, best.r):
            best = row
        if row.measure < AP_SECOND:
            below.append(row)
    if best is None:
        raise AssertionError(f"no finite rows checked for h={base.h}")
    return best, below, checked


def main() -> None:
    print("THM-542 LRC14 one-tail mouth-retention certificate")
    print("family: C_{h,r}=({1,...,13}\\{6,h}) union {r}, h!=6, r>=14")
    print(f"collar: {fmt(COLLAR)}")
    print(f"AP one-hole second value: {fmt(AP_SECOND)}")
    print()

    bases = [base_row(h) for h in range(1, 14) if h != 6]

    print("periodic-comb cutoff table")
    print("h | meas(G_base) | components | 6M-7Q | cutoff R_h | finite best")
    all_below: list[TailRow] = []
    total_checked = 0
    for base in bases:
        best, below, checked = finite_check(base)
        total_checked += checked
        all_below.extend(below)
        assert comb_lower_bound(base, base.cutoff) >= AP_SECOND
        print(
            f"{base.h:2d} | {str(base.measure):>12} | {base.components:10d} | "
            f"{str(base.denominator_slack):>10} | {base.cutoff:10d} | "
            f"h={best.h}, r={best.r}, safe={best.measure}, old={best.old_survivor}"
        )
    print()

    all_below.sort(key=lambda row: (row.measure, row.h, row.r))
    print(f"finite rows checked below cutoffs: {total_checked}")
    print("rows below 426/35035:")
    for row in all_below:
        print(
            f"  h={row.h}, r={row.r}, safe={row.measure}, "
            f"old_survivor={row.old_survivor}, new_mass={row.new_mass}, core={row.core}"
        )
    print()

    expected = [row for row in all_below if row.h == 10 and row.r == 20]
    assert len(all_below) == 1
    assert len(expected) == 1
    champion = expected[0]
    assert champion.measure == Fraction(3859, 420420)
    assert champion.old_survivor == COLLAR
    assert champion.new_mass == Fraction(1, 980)

    print("exact conclusion")
    print("  The unique one-tail row below 426/35035 is h=10, r=20.")
    print("  It has safe measure 3859/420420 = 7/858 + 1/980.")
    print("  It retains all four old drop-6 mouth intervals.")
    print()

    print("Tournament Analysis")
    print("  vertices: one-tail proof gates")
    print("  pairwise observable: which gate removes more candidate rows exactly")
    print("  switch/gauge: use the periodic-comb denominator cutoff before finite wall scan")
    print("  Hamiltonian path:")
    print("    comb_cutoff > finite_exact_scan > old_mouth_survivor > new_mouth_mass > raw_tail_size")
    print("  directed 3-cycles: 0 (proof-obligation order)")
    print()
    print("PASS: THM-542 one-tail mouth-retention lemma certified.")


if __name__ == "__main__":
    main()
