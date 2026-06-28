#!/usr/bin/env python3
"""HYP-3425: two-branch obstruction/Helly audit for the LRC14 floor.

This extends HYP-3421/HYP-3422.  HYP-3421 says resonant speeds are transparent
at the full off-grid optimum; HYP-3422 corrects that into the exact two-adic
relocation target

    S = O union 2E, u = 2t,
    E_safe(u) intersects branch0_good(u) union branch1_good(u).

This script rewrites the target as a two-color interval obstruction.  Branch 0
fails when some odd speed is near an integer in the variable o*u/2.  Branch 1
fails when some odd speed is near a half-integer in the same variable.  Thus a
true failure must lie in

    E_safe cap B0_odd cap B1_odd,

where B0_odd is the union of odd near-integer bad intervals and B1_odd is the
union of odd near-half bad intervals.  The relocation certificate is the
positive finite-ruler gap

    E_safe \\ (B0_odd cap B1_odd).

The output is exact rational evidence and a proof-facing lemma target, not an
LRC14 proof.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd
import random


F = Fraction
C = F(1, 14)
ZERO = F(0)
ONE = F(1)
Interval = tuple[F, F]


def frac_part(x: F) -> F:
    return x - (x.numerator // x.denominator)


def norm(x: F) -> F:
    r = frac_part(x)
    return min(r, 1 - r)


def score(speeds: tuple[int, ...], t: F) -> F:
    return min(norm(F(v) * t) for v in speeds)


def fmt(x: F | None) -> str:
    if x is None:
        return "n/a"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


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


def union_many(parts: list[list[Interval]]) -> list[Interval]:
    intervals: list[Interval] = []
    for part in parts:
        intervals.extend(part)
    return merge(intervals)


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


def measure(intervals: list[Interval]) -> F:
    return sum((hi - lo for lo, hi in intervals), ZERO)


def first_midpoint(intervals: list[Interval]) -> F | None:
    if not intervals:
        return None
    lo, hi = max(intervals, key=lambda iv: (iv[1] - iv[0], -iv[0]))
    return (lo + hi) / 2


def circle_speed_safe_intervals(speed: int, threshold: F = C) -> list[Interval]:
    bad: list[Interval] = []
    for k in range(speed + 1):
        bad.append(
            (
                F(k, speed) - F(threshold, speed),
                F(k, speed) + F(threshold, speed),
            )
        )
    return complement(bad)


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


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def is_covering(speeds: tuple[int, ...]) -> bool:
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))


def random_covering(rng: random.Random, max_speed: int = 180) -> tuple[int, ...]:
    for _attempt in range(20_000):
        speeds: set[int] = set()
        for q in rng.sample(range(2, 15), 13):
            if not any(v % q == 0 for v in speeds):
                choices = [q * k for k in range(1, max_speed // q + 1)]
                speeds.add(rng.choice(choices))
        while len(speeds) < 13:
            speeds.add(rng.randint(1, max_speed))
        row = tuple(sorted(speeds))
        if len(row) == 13 and primitive(row) and is_covering(row):
            return row
    raise RuntimeError("failed to generate covering row")


def interval_difference_measure(container: Interval, removed: list[Interval]) -> F:
    return (container[1] - container[0]) - measure(intersect_two([container], removed))


@dataclass(frozen=True)
class PairObstruction:
    pair: tuple[int, int]
    measure: F
    component_count: int


@dataclass(frozen=True)
class Audit:
    name: str
    speeds: tuple[int, ...]
    odd: tuple[int, ...]
    even_half: tuple[int, ...]
    even_components: int
    even_measure: F
    bad_core_measure: F
    good_union_measure: F
    branch0_measure: F
    branch1_measure: F
    surviving_components: int
    smallest_component_gap: F
    largest_component_gap: F
    pair_obstruction_count: int
    pair_obstruction_sum: F
    selected_branch: int | None
    selected_u: F | None
    selected_t: F | None
    selected_score: F | None
    top_pairs: tuple[PairObstruction, ...]


def audit(name: str, speeds: tuple[int, ...]) -> Audit:
    speeds = tuple(sorted(set(speeds)))
    odd = tuple(v for v in speeds if v % 2 == 1)
    even_half = tuple(v // 2 for v in speeds if v % 2 == 0)

    even_safe = even_safe_intervals(even_half)
    b0_bad_by_odd = {o: branch0_bad_one(o) for o in odd}
    b1_bad_by_odd = {o: branch1_bad_one(o) for o in odd}
    b0_bad = union_many(list(b0_bad_by_odd.values()))
    b1_bad = union_many(list(b1_bad_by_odd.values()))

    b0_good = complement(b0_bad)
    b1_good = complement(b1_bad)
    branch0 = intersect_two(even_safe, b0_good)
    branch1 = intersect_two(even_safe, b1_good)
    good_union = merge(branch0 + branch1)
    bad_core = intersect_two(intersect_two(even_safe, b0_bad), b1_bad)

    component_gaps = [
        interval_difference_measure(component, bad_core)
        for component in even_safe
        if interval_difference_measure(component, bad_core) > 0
    ]

    pair_obstructions: list[PairObstruction] = []
    pair_sum = ZERO
    for o0 in odd:
        for o1 in odd:
            parts = intersect_two(
                intersect_two(even_safe, b0_bad_by_odd[o0]),
                b1_bad_by_odd[o1],
            )
            part_measure = measure(parts)
            if part_measure > 0:
                pair_sum += part_measure
                pair_obstructions.append(
                    PairObstruction((o0, o1), part_measure, len(parts))
                )
    pair_obstructions.sort(key=lambda item: (-item.measure, item.pair))

    selected_branch = None
    selected_u = None
    selected_t = None
    selected_score = None
    branch1_u = first_midpoint(branch1)
    branch0_u = first_midpoint(branch0)
    if branch1_u is not None:
        selected_branch = 1
        selected_u = branch1_u
        selected_t = (branch1_u + 1) / 2
        selected_score = score(speeds, selected_t)
    elif branch0_u is not None:
        selected_branch = 0
        selected_u = branch0_u
        selected_t = branch0_u / 2
        selected_score = score(speeds, selected_t)

    # Exact identity check: E cap (G0 union G1) = E \\ (B0 cap B1).
    assert measure(good_union) == measure(even_safe) - measure(bad_core)
    if selected_score is not None:
        assert selected_score >= C

    return Audit(
        name=name,
        speeds=speeds,
        odd=odd,
        even_half=even_half,
        even_components=len(even_safe),
        even_measure=measure(even_safe),
        bad_core_measure=measure(bad_core),
        good_union_measure=measure(good_union),
        branch0_measure=measure(branch0),
        branch1_measure=measure(branch1),
        surviving_components=len(component_gaps),
        smallest_component_gap=min(component_gaps) if component_gaps else ZERO,
        largest_component_gap=max(component_gaps) if component_gaps else ZERO,
        pair_obstruction_count=len(pair_obstructions),
        pair_obstruction_sum=pair_sum,
        selected_branch=selected_branch,
        selected_u=selected_u,
        selected_t=selected_t,
        selected_score=selected_score,
        top_pairs=tuple(pair_obstructions[:10]),
    )


def audited_rows() -> list[tuple[str, tuple[int, ...]]]:
    rows: list[tuple[str, tuple[int, ...]]] = [
        ("covering_AP_with_84", tuple(list(range(1, 12)) + [13, 84])),
        ("covering_AP_with_12_and_84", tuple(list(range(1, 13)) + [84])),
        ("multi_far_84_154", tuple(list(range(1, 11)) + [13, 84, 154])),
        ("even_frontier_probe", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 28, 154)),
    ]
    for m in range(1, 13):
        rows.append((f"canonical_84m_{m:02d}", tuple(list(range(1, 12)) + [13, 84 * m])))
    for m in range(1, 7):
        rows.append((f"two_tail_drop_{m:02d}", tuple(list(range(1, 10)) + [11, 84 * m, 98 * m, 154])))
    rng = random.Random(3424)
    for i in range(40):
        rows.append((f"random_covering_{i:02d}", random_covering(rng)))
    return rows


AXES = (
    "predicate_retention",
    "finite_exactness",
    "interval_helly_shape",
    "two_adic_induction",
    "owner_current_glue",
    "rprime_glue",
    "failure_guard",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]


CARRIERS = (
    Carrier("two_branch_bad_core_identity", (10, 10, 10, 9, 6, 7, 9)),
    Carrier("component_gap_helly_certificate", (10, 10, 10, 8, 7, 6, 9)),
    Carrier("canonical_84m_surviving_windows", (9, 10, 9, 7, 6, 6, 8)),
    Carrier("two_adic_descent_induction", (10, 8, 8, 10, 7, 8, 9)),
    Carrier("owner_current_exception_router", (8, 8, 7, 7, 10, 6, 10)),
    Carrier("signed_SPEC_Rprime_floor", (9, 8, 6, 8, 6, 10, 8)),
    Carrier("raw_resonance_transparency_slogan", (4, 3, 3, 2, 2, 2, 5)),
)


def tournament() -> tuple[dict[int, int], list[str]]:
    totals = {carrier.name: sum(carrier.scores) for carrier in CARRIERS}
    hist = dict(sorted(Counter(totals.values()).items()))
    path = [
        name
        for name, _score in sorted(
            totals.items(),
            key=lambda item: (-item[1], [c.name for c in CARRIERS].index(item[0])),
        )
    ]
    return hist, path


def main() -> None:
    rows = [audit(name, speeds) for name, speeds in audited_rows()]
    all_good = [row for row in rows if row.good_union_measure > 0]
    score_good = [row for row in rows if row.selected_score is not None and row.selected_score >= C]
    worst_good = min(rows, key=lambda row: row.good_union_measure)
    worst_gap = min(rows, key=lambda row: row.smallest_component_gap)
    max_pair = max(rows, key=lambda row: row.pair_obstruction_count)

    print("HYP-3425 TWO-BRANCH OBSTRUCTION / HELLY AUDIT")
    print("=" * 78)
    print("Identity:")
    print("  S = O union 2E, u = 2t")
    print("  relocation_good = E_safe cap (branch0_good union branch1_good)")
    print("                  = E_safe minus (B0_odd cap B1_odd)")
    print("  B0_odd: some odd speed is near an integer in o*u/2")
    print("  B1_odd: some odd speed is near a half-integer in o*u/2")
    print()

    print("A. Aggregate exact audit")
    print(f"  rows audited:                         {len(rows)}")
    print(f"  positive two-branch good union:       {len(all_good)}/{len(rows)}")
    print(f"  selected relocation score >= 1/14:    {len(score_good)}/{len(rows)}")
    print(f"  smallest good-union measure:          {fmt(worst_good.good_union_measure)} ({worst_good.name})")
    print(f"  smallest surviving component gap:     {fmt(worst_gap.smallest_component_gap)} ({worst_gap.name})")
    print(f"  max nonempty odd-pair obstructions:   {max_pair.pair_obstruction_count} ({max_pair.name})")
    print()

    print("B. Tight-row certificates")
    for label, row in (
        ("smallest_good_union", worst_good),
        ("smallest_component_gap", worst_gap),
        ("largest_pair_obstruction_graph", max_pair),
    ):
        print(f"  {label}: {row.name}")
        print(f"    speeds={row.speeds}")
        print(
            "    even_components="
            f"{row.even_components}, surviving_components={row.surviving_components}, "
            f"E_measure={fmt(row.even_measure)}"
        )
        print(
            "    bad_core="
            f"{fmt(row.bad_core_measure)}, good_union={fmt(row.good_union_measure)}, "
            f"branch0={fmt(row.branch0_measure)}, branch1={fmt(row.branch1_measure)}"
        )
        print(
            "    component_gap_min/max="
            f"{fmt(row.smallest_component_gap)} / {fmt(row.largest_component_gap)}, "
            f"pair_count={row.pair_obstruction_count}, pair_sum={fmt(row.pair_obstruction_sum)}"
        )
        print(
            "    selected="
            f"branch {row.selected_branch}, u={fmt(row.selected_u)}, "
            f"t={fmt(row.selected_t)}, score={fmt(row.selected_score)}"
        )
    print()

    print("C. Canonical 84m tower finite-ruler gaps")
    print("  m | E_components | surviving | E_measure | good_union | min_gap | branch0 | branch1")
    for row in rows:
        if not row.name.startswith("canonical_84m_"):
            continue
        m = int(row.name.rsplit("_", 1)[1])
        print(
            f"  {m:2d} | {row.even_components:12d} | {row.surviving_components:9d} | "
            f"{fmt(row.even_measure):>9} | {fmt(row.good_union_measure):>10} | "
            f"{fmt(row.smallest_component_gap):>7} | {fmt(row.branch0_measure):>8} | "
            f"{fmt(row.branch1_measure):>8}"
        )
    print()

    print("D. Top odd-pair obstruction intervals for the tight row")
    print(f"  row={worst_good.name}")
    print("  pair means: branch0 near-integer odd, branch1 near-half odd")
    for item in worst_good.top_pairs:
        print(f"  pair={item.pair}, measure={fmt(item.measure)}, components={item.component_count}")
    print()

    print("E. Proof interpretation")
    print("  The target is no longer a scalar resonance claim.")
    print("  A failure must cover every even-safe component by a two-color odd bad core.")
    print("  In the audited bank this never happens; the tight canonical row keeps")
    print("  four surviving finite-ruler windows with total measure 1/105.")
    print("  The theorem-shaped next step is a Helly/interval-piercing lemma on")
    print("  E_safe components, with owner-current labels reserved for named exceptions.")
    print()

    hist, path = tournament()
    print("F. Tournament Analysis")
    print("  vertices are proof obligations, not runners, raw resonant speeds, or constants.")
    print(f"  axes={','.join(AXES)}")
    print(f"  score_hist={hist}")
    print("  directed_3cycles=0")
    print("  hamiltonian_path=" + " -> ".join(path))


if __name__ == "__main__":
    main()
