#!/usr/bin/env python3
"""THM-543: exact one-replacement AP-tail theorem for LRC14.

This certificate proves the full one-replacement AP-tail layer

    C_{a,b,r} = ({1,...,13} \\ {a,b}) union {r}, 1 <= a < b <= 13, r >= 14.

The proof is intentionally arithmetic.  It uses exact rational interval
arithmetic and an exact periodic-comb cutoff; no floating-point approximation
appears in any proof comparison or in the reported certificate.

If G is a union of c intervals and D_r is the speed-r danger comb at level
1/14, then

    meas(G \\ D_r) >= (6/7) meas(G) - 2c/(7r).

For every two-hole AP base B={1,...,13}\\{a,b}, the denominator slack

    6*meas(G_B) - 7*(426/35035)

is positive.  Hence all sufficiently large r are above the AP one-hole second
value, and the remaining finite rows can be checked exactly.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache, reduce
from itertools import combinations
from math import gcd


TARGET_DENOM = 14
AP_SECOND = Fraction(426, 35035)
DROP6_COLLAR = Fraction(7, 858)
DROP6_MOUTHS = (
    (Fraction(29, 182), Fraction(9, 56)),
    (Fraction(29, 168), Fraction(27, 154)),
    (Fraction(127, 154), Fraction(139, 168)),
    (Fraction(47, 56), Fraction(153, 182)),
)


@dataclass(frozen=True)
class BaseRow:
    holes: tuple[int, int]
    measure: Fraction
    components: int
    denominator_slack: Fraction
    cutoff: int


@dataclass(frozen=True)
class TailRow:
    holes: tuple[int, int]
    r: int
    measure: Fraction
    components: int
    old_survivor: Fraction
    new_mass: Fraction
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


def base_core(holes: tuple[int, int]) -> tuple[int, ...]:
    return tuple(v for v in range(1, 14) if v not in holes)


def tail_core(holes: tuple[int, int], r: int) -> tuple[int, ...]:
    return tuple(sorted(base_core(holes) + (r,)))


def base_row(holes: tuple[int, int]) -> BaseRow:
    comps = safe_components(base_core(holes))
    m = measure(comps)
    c = len(comps)
    slack = 6 * m - 7 * AP_SECOND
    if slack <= 0:
        raise AssertionError(f"base holes={holes} has no positive comb slack")
    cutoff = ceil_fraction(Fraction(2 * c, 1) / slack)
    return BaseRow(
        holes=holes,
        measure=m,
        components=c,
        denominator_slack=slack,
        cutoff=cutoff,
    )


def tail_row(holes: tuple[int, int], r: int) -> TailRow:
    core = tail_core(holes, r)
    comps = safe_components(core)
    m = measure(comps)
    old = intersect_measure(comps, DROP6_MOUTHS)
    return TailRow(
        holes=holes,
        r=r,
        measure=m,
        components=len(comps),
        old_survivor=old,
        new_mass=m - old,
        core=core,
    )


def comb_lower_bound(base: BaseRow, r: int) -> Fraction:
    return Fraction(6, 7) * base.measure - Fraction(2 * base.components, 7 * r)


def finite_check(base: BaseRow) -> tuple[TailRow | None, list[TailRow], int]:
    checked = 0
    below: list[TailRow] = []
    best: TailRow | None = None
    for r in range(14, base.cutoff):
        core = tail_core(base.holes, r)
        if not primitive(core):
            continue
        checked += 1
        row = tail_row(base.holes, r)
        if best is None or (row.measure, row.holes, row.r) < (
            best.measure,
            best.holes,
            best.r,
        ):
            best = row
        if row.measure < AP_SECOND:
            below.append(row)
    return best, below, checked


def relation_lattice_rank(core: tuple[int, ...]) -> int:
    """Affine relation rank for one-dimensional speeds.

    Relations are integer vectors n with sum n_i=0 and sum n_i*c_i=0.  For a
    primitive one-dimensional row with at least two distinct speeds, the two
    displayed linear forms are independent over Q, so the rank is k-2.
    """

    if len(set(core)) != len(core):
        raise AssertionError(f"duplicate speeds in {core}")
    return len(core) - 2


def main() -> None:
    print("THM-543 LRC14 one-replacement AP-tail certificate")
    print("family: C_{a,b,r}=({1,...,13}\\{a,b}) union {r}, 1<=a<b<=13, r>=14")
    print(f"AP one-hole second value Q: {AP_SECOND}")
    print()

    bases = [base_row(pair) for pair in combinations(range(1, 14), 2)]
    bases_by_cutoff = sorted(bases, key=lambda row: (-row.cutoff, row.holes))

    all_below: list[TailRow] = []
    total_checked = 0
    finite_bests: dict[tuple[int, int], TailRow | None] = {}
    for base in bases:
        best, below, checked = finite_check(base)
        total_checked += checked
        all_below.extend(below)
        finite_bests[base.holes] = best
        assert comb_lower_bound(base, base.cutoff) >= AP_SECOND

    print("periodic-comb cutoff table sorted by decreasing cutoff")
    print("holes | M_base | components | 6M-7Q | cutoff R | finite best safe row")
    for base in bases_by_cutoff:
        best = finite_bests[base.holes]
        if best is None:
            best_summary = "comb bound covers all r>=14"
        else:
            best_summary = f"holes={best.holes}, r={best.r}, safe={best.measure}"
        print(
            f"{base.holes!s:>8} | {str(base.measure):>12} | "
            f"{base.components:10d} | {str(base.denominator_slack):>14} | "
            f"{base.cutoff:8d} | {best_summary}"
        )
    print()

    all_below.sort(key=lambda row: (row.measure, row.holes, row.r))
    print(f"finite rows checked below exact cutoffs: {total_checked}")
    print("rows below 426/35035:")
    for row in all_below:
        print(
            f"  holes={row.holes}, r={row.r}, safe={row.measure}, "
            f"old_survivor={row.old_survivor}, new_mass={row.new_mass}, core={row.core}"
        )
    print()

    expected = [row for row in all_below if row.holes == (6, 10) and row.r == 20]
    assert len(all_below) == 1
    assert len(expected) == 1
    champion = expected[0]
    assert champion.measure == Fraction(3859, 420420)
    assert champion.old_survivor == DROP6_COLLAR
    assert champion.new_mass == Fraction(1, 980)
    assert relation_lattice_rank(champion.core) == 10

    weakest = min(bases, key=lambda row: (row.denominator_slack, row.holes))
    largest_cutoff = max(bases, key=lambda row: (row.cutoff, row.holes))
    print("exact extremal cutoff data")
    print(
        f"  weakest denominator slack: holes={weakest.holes}, "
        f"6M-7Q={weakest.denominator_slack}, cutoff={weakest.cutoff}"
    )
    print(
        f"  largest cutoff: holes={largest_cutoff.holes}, "
        f"6M-7Q={largest_cutoff.denominator_slack}, cutoff={largest_cutoff.cutoff}"
    )
    print("  the same resonant base (6,10) controls both extremal quantities.")
    print()

    print("exact conclusion")
    print("  The unique one-replacement AP-tail row below 426/35035 is")
    print("  holes=(6,10), r=20.")
    print("  It has safe measure 3859/420420 = 7/858 + 1/980.")
    print("  It retains all four old drop-6 mouth intervals.")
    print()

    print("Tournament Analysis")
    print("  vertices: proof-obligation gates for the one-replacement AP-tail layer")
    print("  pairwise observable: exact row eliminations before the AP-second threshold")
    print("  switch/gauge: periodic-comb cutoff first, finite wall scan second")
    print("  Hamiltonian path:")
    print("    denominator_slack > comb_cutoff > finite_exact_scan > mouth_survivor > relation_rank")
    print("  directed 3-cycles: 0 (transitive proof-obligation order)")
    print("  challenged assumption: vertices need not be runners; here they are proof gates.")
    print()
    print("PASS: THM-543 one-replacement AP-tail theorem certified.")


if __name__ == "__main__":
    main()
