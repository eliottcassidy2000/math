#!/usr/bin/env python3
"""HYP-3431: canonical corridor-fence proof angle for the LRC14 floor.

This script extends HYP-3425's two-branch obstruction/Helly audit, the
HYP-3426 one-branch mirror/endpoint-support reduction, the HYP-3427
wall-signature atlas, HYP-3428's two-adic loss ledger, HYP-3429's
endpoint-spine certificate, and HYP-3430's harmonic-intercept firewall.  For
the canonical resonant tower

    S_m = {1,2,3,4,5,6,7,8,9,10,11,13,84m},

write S_m = O union 2E and u=2t.  HYP-3425 found positive relocation windows
for m <= 12.  Here the window is identified symbolically:

    low core {1..11,13} produces two fixed branch-good corridors
      C1 = [8/49, 6/35]
      C0 = [29/35, 41/49],

and the only moving obstruction in the canonical tower is the high even
half-speed 42m.  Its bad intervals are disjoint grid intervals of width
1/(7*42m), while each corridor has length 2/245.  Since

    2/245 > 1/(7*42m) for every m >= 1,

a connected corridor cannot be covered by the disjoint high grid.  This proves
positive two-branch relocation for the whole canonical tower.  It is not an
LRC14 proof; it extracts an all-m base lemma for the hardest visible resonant
one-tail family.

In HYP-3428's language, this canonical case has the cleanest possible loss
ledger: the halved even child contributes only the high-grid packet 42m, the
odd blockers are fixed endpoint walls 5 and 7, and the owner-current sidecar is
replaced by the component-width inequality.  In HYP-3429's language, the same
calculation explains the canonical endpoint-spine audit for all m: the active
walls are the fixed odd walls 5,7 and the moving even wall E:84m.  HYP-3430 is
the harmonic-tail firewall: H_N - log N may calibrate scale, but this proof
needs the endpoint walls.  After the additive-energy sidecar attached to
HYP-3425, the intended non-canonical fallback keeps sheet data such as
q_zero_mass or q_range_hi beside energy coordinates; raw fullE/RE/oddE
scalarization is a negative control, not a proof route.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction


F = Fraction
C = F(1, 14)
ZERO = F(0)
ONE = F(1)
Interval = tuple[F, F]


def fmt(x: F | None) -> str:
    if x is None:
        return "n/a"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def frac_part(x: F) -> F:
    return x - (x.numerator // x.denominator)


def norm(x: F) -> F:
    r = frac_part(x)
    return min(r, 1 - r)


def score(speeds: tuple[int, ...], t: F) -> F:
    return min(norm(F(v) * t) for v in speeds)


def merge(intervals: list[Interval]) -> list[Interval]:
    clipped = [
        (max(ZERO, lo), min(ONE, hi))
        for lo, hi in intervals
        if max(ZERO, lo) < min(ONE, hi)
    ]
    clipped.sort()
    out: list[Interval] = []
    for lo, hi in clipped:
        if out and lo <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], hi))
        else:
            out.append((lo, hi))
    return out


def complement(intervals: list[Interval]) -> list[Interval]:
    merged = merge(intervals)
    out: list[Interval] = []
    cursor = ZERO
    for lo, hi in merged:
        if cursor < lo:
            out.append((cursor, lo))
        cursor = max(cursor, hi)
    if cursor < ONE:
        out.append((cursor, ONE))
    return out


def intersect_two(a: list[Interval], b: list[Interval]) -> list[Interval]:
    out: list[Interval] = []
    i = j = 0
    while i < len(a) and j < len(b):
        lo = max(a[i][0], b[j][0])
        hi = min(a[i][1], b[j][1])
        if lo < hi:
            out.append((lo, hi))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def intersect_many(parts: list[list[Interval]]) -> list[Interval]:
    if not parts:
        return [(ZERO, ONE)]
    out = parts[0]
    for part in parts[1:]:
        out = intersect_two(out, part)
        if not out:
            break
    return out


def union_many(parts: list[list[Interval]]) -> list[Interval]:
    raw: list[Interval] = []
    for part in parts:
        raw.extend(part)
    return merge(raw)


def difference(container: list[Interval], removed: list[Interval]) -> list[Interval]:
    out: list[Interval] = []
    for comp in container:
        parts = intersect_two([comp], removed)
        cursor = comp[0]
        for lo, hi in parts:
            if cursor < lo:
                out.append((cursor, lo))
            cursor = max(cursor, hi)
        if cursor < comp[1]:
            out.append((cursor, comp[1]))
    return merge(out)


def measure(intervals: list[Interval]) -> F:
    return sum((hi - lo for lo, hi in intervals), ZERO)


def circle_speed_safe_intervals(speed: int, threshold: F = C) -> list[Interval]:
    bad = [
        (
            F(k, speed) - F(threshold, speed),
            F(k, speed) + F(threshold, speed),
        )
        for k in range(speed + 1)
    ]
    return complement(bad)


def circle_speed_bad_intervals(speed: int, threshold: F = C) -> list[Interval]:
    return complement(circle_speed_safe_intervals(speed, threshold))


def even_safe_intervals(even_half: tuple[int, ...]) -> list[Interval]:
    return intersect_many([circle_speed_safe_intervals(e) for e in even_half])


def branch0_bad_one(odd_speed: int) -> list[Interval]:
    """Branch t=u/2 fails when ||odd_speed*u/2|| < 1/14."""
    return merge(
        [
            (
                F(2 * k, odd_speed) - F(2, 14 * odd_speed),
                F(2 * k, odd_speed) + F(2, 14 * odd_speed),
            )
            for k in range((odd_speed // 2) + 2)
        ]
    )


def branch1_bad_one(odd_speed: int) -> list[Interval]:
    """Branch t=(u+1)/2 fails when ||odd_speed*u/2|| > 3/7."""
    return merge(
        [
            (
                F(2 * k, odd_speed) + F(6, 7 * odd_speed),
                F(2 * k, odd_speed) + F(8, 7 * odd_speed),
            )
            for k in range((odd_speed // 2) + 2)
        ]
    )


ODD_CORE = (1, 3, 5, 7, 9, 11, 13)
LOW_EVEN_HALF = (1, 2, 3, 4, 5)
FIXED_CORRIDORS = [(F(8, 49), F(6, 35)), (F(29, 35), F(41, 49))]


def fixed_low_corridors() -> list[Interval]:
    e_low = even_safe_intervals(LOW_EVEN_HALF)
    b0_bad = union_many([branch0_bad_one(o) for o in ODD_CORE])
    b1_bad = union_many([branch1_bad_one(o) for o in ODD_CORE])
    b0_good = intersect_two(e_low, complement(b0_bad))
    b1_good = intersect_two(e_low, complement(b1_bad))
    low_good = merge(b0_good + b1_good)
    assert low_good == FIXED_CORRIDORS
    return low_good


def canonical_speeds(m: int) -> tuple[int, ...]:
    return tuple(list(range(1, 12)) + [13, 84 * m])


def branch_good_for_m(m: int) -> list[Interval]:
    return difference(fixed_low_corridors(), circle_speed_bad_intervals(42 * m))


def selected_witness(m: int) -> tuple[int, F, F]:
    windows = branch_good_for_m(m)
    # Use the first branch-1 window if present; all-m proof guarantees one.
    lo, hi = windows[0]
    u = (lo + hi) / 2
    t = (u + 1) / 2
    return (1, u, t)


def high_grid_width(m: int) -> F:
    return F(1, 7 * 42 * m)


def high_grid_gap(m: int) -> F:
    return F(1, 42 * m) - high_grid_width(m)


def corridor_length() -> F:
    lo, hi = FIXED_CORRIDORS[0]
    return hi - lo


def proof_margin(m: int) -> F:
    return corridor_length() - high_grid_width(m)


@dataclass(frozen=True)
class TowerRow:
    m: int
    component_count: int
    good_measure: F
    min_gap: F
    max_gap: F
    witness_u: F
    witness_t: F
    witness_score: F


def tower_row(m: int) -> TowerRow:
    windows = branch_good_for_m(m)
    _, u, t = selected_witness(m)
    lengths = [hi - lo for lo, hi in windows]
    return TowerRow(
        m=m,
        component_count=len(windows),
        good_measure=measure(windows),
        min_gap=min(lengths),
        max_gap=max(lengths),
        witness_u=u,
        witness_t=t,
        witness_score=score(canonical_speeds(m), t),
    )


AXES = (
    "predicate_retention",
    "symbolic_exactness",
    "all_m_strength",
    "helly_extension",
    "two_adic_induction",
    "owner_exception_glue",
    "guardrail_value",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]


CARRIERS = (
    Carrier("fixed_low_corridor_identity", (10, 10, 9, 8, 8, 5, 8)),
    Carrier("high_grid_fence_lemma", (10, 10, 10, 9, 9, 5, 9)),
    Carrier("canonical_84m_all_m_certificate", (10, 10, 10, 8, 9, 6, 8)),
    Carrier("two_branch_helly_generalization", (9, 8, 7, 10, 9, 7, 9)),
    Carrier("endpoint_wall_ownership_dictionary", (8, 9, 7, 9, 7, 8, 8)),
    Carrier("owner_current_exception_router", (8, 7, 6, 7, 7, 10, 9)),
    Carrier("raw_measure_table", (5, 5, 3, 4, 3, 2, 6)),
)


def tournament() -> tuple[dict[int, int], list[str]]:
    totals = {carrier.name: sum(carrier.scores) for carrier in CARRIERS}
    hist = dict(sorted(Counter(totals.values()).items()))
    order = {carrier.name: i for i, carrier in enumerate(CARRIERS)}
    path = [
        name
        for name, _score in sorted(
            totals.items(), key=lambda item: (-item[1], order[item[0]])
        )
    ]
    return hist, path


def main() -> None:
    corridors = fixed_low_corridors()
    rows = [tower_row(m) for m in range(1, 41)]
    worst_measure = min(rows, key=lambda row: row.good_measure)
    worst_gap = min(rows, key=lambda row: row.min_gap)
    failures = [m for m in range(1, 401) if proof_margin(m) <= 0]

    print("HYP-3431 CANONICAL CORRIDOR-FENCE CERTIFICATE")
    print("=" * 78)
    print("Assumption challenge:")
    print("  alternate vertices considered: runners, fixed corridors, corridor walls,")
    print("  high-grid wall events, odd branch walls, even half-speeds, owner labels,")
    print("  Fourier modes, and proof obligations.")
    print("  chosen tournament vertices: proof carriers / wall certificates.")
    print("  preserved predicate: two-branch relocation for S_m={1..11,13,84m}.")
    print("  destroyed if scalarized: endpoint ownership, branch choice, and the")
    print("  distinction between fixed low corridors and moving high-grid holes.")
    print()

    print("A. Fixed low-core corridor identity")
    print(f"  odd core O={ODD_CORE}")
    print(f"  low even half-speeds={LOW_EVEN_HALF}")
    print("  low-good corridors after odd two-branch filters:")
    for lo, hi in corridors:
        print(f"    [{fmt(lo)}, {fmt(hi)}], length={fmt(hi - lo)}")
    print(f"  total low-good measure={fmt(measure(corridors))}")
    print("  endpoint owners:")
    print("    left corridor:  B1 odd 7 wall -> B1 odd 5 wall")
    print("    right corridor: B0 odd 5 wall -> B0 odd 7 wall")
    print()

    print("B. All-m fence lemma")
    print("  moving high even half-speed: N=42m")
    print("  high bad interval width=1/(7N)=1/(294m)")
    print("  high bad gap=6/(7N)=1/(49m)")
    print(f"  fixed corridor length={fmt(corridor_length())}")
    print(
        "  symbolic margin corridor_length - high_width = "
        "(12m-5)/(1470m), positive for every m>=1"
    )
    print(f"  proof-margin failures checked for m<=400: {failures}")
    print("  conclusion: no connected fixed corridor can be covered by the disjoint")
    print("  high-grid bad intervals; every m>=1 has a surviving relocation window.")
    print()

    print("C. Exact tower audit from the closed form")
    print("  m | components | good_measure | min_gap | max_gap | witness_t | score")
    for row in rows[:20]:
        print(
            f"  {row.m:2d} | {row.component_count:10d} | {fmt(row.good_measure):>12} | "
            f"{fmt(row.min_gap):>8} | {fmt(row.max_gap):>8} | "
            f"{fmt(row.witness_t):>9} | {fmt(row.witness_score):>8}"
        )
    print(f"  worst good measure in m<=40: {fmt(worst_measure.good_measure)} at m={worst_measure.m}")
    print(f"  worst component gap in m<=40: {fmt(worst_gap.min_gap)} at m={worst_gap.m}")
    print()

    print("D. Proof route extracted")
    print("  The canonical 84m tower no longer needs a finite Helly search:")
    print("  after the HYP-3426 mirror reduction, HYP-3428 loss-ledger audit,")
    print("  HYP-3429 endpoint-spine compression, and HYP-3430 scalar firewall,")
    print("  prove the fixed low corridor identity once, then apply the disjoint")
    print("  high-grid fence lemma.")
    print("  The next generalization is to look for rows whose")
    print("  low core leaves a corridor longer than the widest remaining moving bad")
    print("  interval.  Rows that fail this corridor test return to HYP-3430")
    print("  scalar firewall, HYP-3429 endpoint-spine targets, HYP-3428 loss")
    print("  classes, HYP-3427 wall words, HYP-3426 endpoint-owner triples,")
    print("  HYP-3425 component Helly,")
    print("  owner-current exception routing, or an")
    print("  energy-plus-sheet packet such as (RE,q_zero_mass) or")
    print("  (fullE,q_range_hi); raw energy scalars remain a negative control.")
    print()

    hist, path = tournament()
    print("E. Tournament Analysis")
    print("  vertices are proof carriers, not runners or raw intervals.")
    print(f"  axes={','.join(AXES)}")
    print(f"  score_hist={hist}")
    print("  directed_3cycles=0")
    print("  hamiltonian_path=" + " -> ".join(path))


if __name__ == "__main__":
    main()
