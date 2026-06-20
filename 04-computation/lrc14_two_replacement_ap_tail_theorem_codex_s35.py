#!/usr/bin/env python3
"""THM-544: exact two-replacement AP-tail theorem for LRC14.

This certificate proves the two-replacement AP-tail layer

    C_{a,b,c,r,s} = ({1,...,13} \\ {a,b,c}) union {r,s},
    1 <= a < b < c <= 13, 14 <= r < s.

The proof is exact rational arithmetic.  It first removes all sufficiently
large pairs (r,s) by a two-comb union bound, then fixes the smaller tail r and
removes all sufficiently large s by a one-comb bound on the exact row after r.
Only the remaining finite pairs are scanned by exact interval subtraction.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache, reduce
from itertools import combinations
from math import gcd


TARGET_DENOM = 14
AP_SECOND = Fraction(426, 35035)
DROP6_MOUTHS = (
    (Fraction(29, 182), Fraction(9, 56)),
    (Fraction(29, 168), Fraction(27, 154)),
    (Fraction(127, 154), Fraction(139, 168)),
    (Fraction(47, 56), Fraction(153, 182)),
)


@dataclass(frozen=True)
class BaseTriple:
    holes: tuple[int, int, int]
    measure: Fraction
    components: int
    two_tail_slack: Fraction
    cutoff: int


@dataclass(frozen=True)
class FixedTail:
    holes: tuple[int, int, int]
    r: int
    measure: Fraction
    components: int
    one_tail_slack: Fraction
    cutoff: int


@dataclass(frozen=True)
class PairRow:
    holes: tuple[int, int, int]
    r: int
    s: int
    measure: Fraction
    old_survivor: Fraction
    core: tuple[int, ...]


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


def subtract_arcs(
    safe: tuple[tuple[Fraction, Fraction], ...],
    arcs: tuple[tuple[Fraction, Fraction], ...],
) -> tuple[tuple[Fraction, Fraction], ...]:
    out: list[tuple[Fraction, Fraction]] = []
    j = 0
    for lo, hi in safe:
        cursor = lo
        while j < len(arcs) and arcs[j][1] <= lo:
            j += 1
        k = j
        while k < len(arcs) and arcs[k][0] < hi:
            arc_lo, arc_hi = arcs[k]
            if arc_lo > cursor:
                out.append((cursor, min(arc_lo, hi)))
            if arc_hi > cursor:
                cursor = max(cursor, min(arc_hi, hi))
            if cursor >= hi:
                break
            k += 1
        if cursor < hi:
            out.append((cursor, hi))
    return tuple(out)


@lru_cache(maxsize=None)
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
    xs: tuple[tuple[Fraction, Fraction], ...],
    ys: tuple[tuple[Fraction, Fraction], ...],
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


def base_core(holes: tuple[int, int, int]) -> tuple[int, ...]:
    return tuple(v for v in range(1, 14) if v not in holes)


def base_triple(holes: tuple[int, int, int]) -> tuple[BaseTriple, tuple[tuple[Fraction, Fraction], ...]]:
    comps = safe_components(base_core(holes))
    m = measure(comps)
    c = len(comps)
    slack = 5 * m - 7 * AP_SECOND
    if slack <= 0:
        raise AssertionError(f"base holes={holes} has no positive two-tail slack")
    cutoff = ceil_fraction(Fraction(4 * c, 1) / slack)
    return BaseTriple(holes, m, c, slack, cutoff), comps


def fixed_tail(
    base: BaseTriple,
    base_comps: tuple[tuple[Fraction, Fraction], ...],
    r: int,
) -> tuple[FixedTail, tuple[tuple[Fraction, Fraction], ...]]:
    comps = subtract_arcs(base_comps, danger_arcs(r))
    m = measure(comps)
    c = len(comps)
    slack = 6 * m - 7 * AP_SECOND
    if slack <= 0:
        raise AssertionError(f"fixed row holes={base.holes}, r={r} has no slack")
    cutoff = ceil_fraction(Fraction(2 * c, 1) / slack)
    return FixedTail(base.holes, r, m, c, slack, cutoff), comps


def pair_row(
    fixed: FixedTail,
    fixed_comps: tuple[tuple[Fraction, Fraction], ...],
    s: int,
) -> PairRow:
    comps = subtract_arcs(fixed_comps, danger_arcs(s))
    m = measure(comps)
    old = intersect_measure(comps, DROP6_MOUTHS)
    core = tuple(sorted(base_core(fixed.holes) + (fixed.r, s)))
    return PairRow(fixed.holes, fixed.r, s, m, old, core)


def two_comb_lower_bound(base: BaseTriple, r: int, s: int) -> Fraction:
    return (
        Fraction(5, 7) * base.measure
        - Fraction(2 * base.components, 7 * r)
        - Fraction(2 * base.components, 7 * s)
    )


def one_comb_lower_bound(fixed: FixedTail, s: int) -> Fraction:
    return Fraction(6, 7) * fixed.measure - Fraction(2 * fixed.components, 7 * s)


def main() -> None:
    print("THM-544 LRC14 two-replacement AP-tail certificate")
    print(
        "family: C=({1,...,13}\\{a,b,c}) union {r,s}, "
        "1<=a<b<c<=13, 14<=r<s"
    )
    print(f"AP one-hole second value Q: {AP_SECOND}")
    print()

    bases: list[BaseTriple] = []
    fixed_rows = 0
    finite_pairs = 0
    below: list[PairRow] = []
    best: PairRow | None = None
    fixed_cutoff_rows: list[FixedTail] = []

    for holes in combinations(range(1, 14), 3):
        base, base_comps = base_triple(holes)
        bases.append(base)
        assert two_comb_lower_bound(base, base.cutoff, base.cutoff) >= AP_SECOND
        for r in range(14, base.cutoff):
            first_core = tuple(sorted(base_core(holes) + (r,)))
            if not primitive(first_core):
                continue
            fixed, fixed_comps = fixed_tail(base, base_comps, r)
            fixed_rows += 1
            fixed_cutoff_rows.append(fixed)
            if fixed.cutoff > r + 1:
                assert one_comb_lower_bound(fixed, fixed.cutoff) >= AP_SECOND
            for s in range(r + 1, fixed.cutoff):
                core = tuple(sorted(first_core + (s,)))
                if not primitive(core):
                    continue
                finite_pairs += 1
                row = pair_row(fixed, fixed_comps, s)
                if best is None or (row.measure, row.holes, row.r, row.s) < (
                    best.measure,
                    best.holes,
                    best.r,
                    best.s,
                ):
                    best = row
                if row.measure < AP_SECOND:
                    below.append(row)

    if best is None:
        raise AssertionError("finite scan unexpectedly empty")
    assert not below
    assert best.holes == (4, 6, 10)
    assert best.r == 14
    assert best.s == 15
    assert best.measure == Fraction(14249, 252252)

    bases_by_cutoff = sorted(bases, key=lambda row: (-row.cutoff, row.holes))
    fixed_active = [row for row in fixed_cutoff_rows if row.cutoff > row.r + 1]
    fixed_by_cutoff = sorted(fixed_active, key=lambda row: (-row.cutoff, row.holes, row.r))

    print("largest two-comb base cutoffs")
    print("holes | M_base | components | 5M-7Q | cutoff R")
    for base in bases_by_cutoff[:20]:
        print(
            f"{base.holes!s:>11} | {str(base.measure):>12} | "
            f"{base.components:10d} | {str(base.two_tail_slack):>14} | {base.cutoff:8d}"
        )
    print()

    print("largest active fixed-tail cutoffs")
    print("holes, r | M_after_r | components | 6M-7Q | cutoff S")
    for fixed in fixed_by_cutoff[:20]:
        print(
            f"{fixed.holes!s:>11}, {fixed.r:3d} | {str(fixed.measure):>12} | "
            f"{fixed.components:10d} | {str(fixed.one_tail_slack):>14} | {fixed.cutoff:8d}"
        )
    print()

    print(f"three-hole bases checked: {len(bases)}")
    print(f"fixed smaller-tail rows checked for cutoff: {fixed_rows}")
    print(f"finite two-tail pairs checked exactly: {finite_pairs}")
    print(f"rows below 426/35035: {len(below)}")
    print()

    print("exact minimum finite row")
    print(
        f"  holes={best.holes}, tails=({best.r},{best.s}), safe={best.measure}, "
        f"gap_above_Q={best.measure - AP_SECOND}, old_survivor={best.old_survivor}, "
        f"core={best.core}"
    )
    print()

    worst_base = bases_by_cutoff[0]
    worst_fixed = fixed_by_cutoff[0]
    print("exact extremal cutoff data")
    print(
        f"  largest two-comb cutoff: holes={worst_base.holes}, "
        f"M={worst_base.measure}, c={worst_base.components}, "
        f"5M-7Q={worst_base.two_tail_slack}, R={worst_base.cutoff}"
    )
    print(
        f"  largest active fixed-tail cutoff: holes={worst_fixed.holes}, r={worst_fixed.r}, "
        f"M={worst_fixed.measure}, c={worst_fixed.components}, "
        f"6M-7Q={worst_fixed.one_tail_slack}, S={worst_fixed.cutoff}"
    )
    print()

    print("exact conclusion")
    print("  Every two-replacement AP-tail row is at least 426/35035.")
    print("  The exact finite minimum is 14249/252252, exceeding Q by 1141/25740.")
    print()

    print("Tournament Analysis")
    print("  vertices: two-tail proof gates")
    print("  pairwise observable: exact eliminations before threshold Q")
    print("  switch/gauge: two-comb cutoff, then fixed-tail one-comb cutoff, then scan")
    print("  Hamiltonian path:")
    print("    base_two_tail_slack > fixed_tail_slack > finite_scan > best_gap > mouth_survivor")
    print("  directed 3-cycles: 0 (transitive proof-obligation order)")
    print("  challenged assumption: the proof vertices are gates, not runners or arcs.")
    print()
    print("PASS: THM-544 two-replacement AP-tail theorem certified.")


if __name__ == "__main__":
    main()
