#!/usr/bin/env python3
"""S258 two-frontier observer gluing scout for LRC14.

This script works two proof obligations back and forth:

A. HYP-3096 polynomial-method witness route:
   direct lonely measure, component count, largest direct arc, finite
   bad-denominator budget.

B. HYP-3097 Pascal/equivalence/scissors route:
   pair-normalized cap skeleton, cap defect, Farey pair-mass signal, and
   branch/K33 handoff debt.

The point is not to prove LRC14.  The point is to expose which observer chart
has to pay when another chart fails.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import comb, gcd


THRESHOLD = Fraction(1, 14)
AP = tuple(range(1, 14))


@dataclass(frozen=True)
class Row:
    name: str
    speeds: tuple[int, ...]
    route: str
    chart_hint: str


def lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


def lcm_upto(n: int) -> int:
    out = 1
    for v in range(1, n + 1):
        out = lcm(out, v)
    return out


def mod1(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def norm(x: Fraction) -> Fraction:
    r = mod1(x)
    return r if r <= Fraction(1, 2) else 1 - r


def fmt(x: Fraction) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def gcd_all(values: tuple[int, ...]) -> int:
    out = 0
    for v in values:
        out = gcd(out, v)
    return out


def is_witness(speeds: tuple[int, ...], t: Fraction, strict: bool = False) -> bool:
    md = min(norm(v * t) for v in speeds)
    return md > THRESHOLD if strict else md >= THRESHOLD


def has_grid_witness(speeds: tuple[int, ...], d: int) -> bool:
    return any(is_witness(speeds, Fraction(a, d)) for a in range(d))


def danger_intervals(speed: int) -> list[tuple[Fraction, Fraction]]:
    radius = THRESHOLD / speed
    intervals: list[tuple[Fraction, Fraction]] = []
    for k in range(speed):
        center = Fraction(k, speed)
        lo = center - radius
        hi = center + radius
        if lo < 0:
            intervals.append((lo + 1, Fraction(1)))
            intervals.append((Fraction(0), hi))
        elif hi > 1:
            intervals.append((lo, Fraction(1)))
            intervals.append((Fraction(0), hi - 1))
        else:
            intervals.append((lo, hi))
    return intervals


def merge_intervals(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    intervals = sorted(intervals)
    if not intervals:
        return []
    merged = [intervals[0]]
    for lo, hi in intervals[1:]:
        old_lo, old_hi = merged[-1]
        if lo <= old_hi:
            if hi > old_hi:
                merged[-1] = (old_lo, hi)
        else:
            merged.append((lo, hi))
    return merged


def safe_components(speeds: tuple[int, ...]) -> list[tuple[Fraction, Fraction]]:
    intervals: list[tuple[Fraction, Fraction]] = []
    for v in speeds:
        intervals.extend(danger_intervals(v))
    comps: list[tuple[Fraction, Fraction]] = []
    cursor = Fraction(0)
    for lo, hi in merge_intervals(intervals):
        if cursor < lo:
            comps.append((cursor, lo))
        if hi > cursor:
            cursor = hi
    if cursor < 1:
        comps.append((cursor, Fraction(1)))
    return comps


def strict_grid_bound(length: Fraction) -> int | None:
    if length <= 0:
        return None
    return length.denominator // length.numerator + 1


def candidate_taus(speeds: tuple[int, ...]) -> set[Fraction]:
    speeds = tuple(sorted(set(speeds)))
    cands: set[Fraction] = {Fraction(1, 2)}
    for v in speeds:
        k = 0
        while True:
            t = Fraction(2 * k + 1, 2 * v)
            if t > Fraction(1, 2):
                break
            cands.add(t)
            k += 1
    for a, b in combinations(speeds, 2):
        for denom in (a + b, abs(b - a)):
            if denom <= 0:
                continue
            k = 1
            while True:
                t = Fraction(k, denom)
                if t > Fraction(1, 2):
                    break
                cands.add(t)
                k += 1
    return cands


def m_exact(speeds: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    best = Fraction(0)
    arg = Fraction(0)
    for tau in candidate_taus(speeds):
        val = min(norm(v * tau) for v in speeds)
        if val > best:
            best = val
            arg = tau
    return best, arg


def binder_word(speeds: tuple[int, ...], tau: Fraction, m_value: Fraction) -> str:
    binders = [v for v in speeds if norm(v * tau) == m_value]
    return ",".join(f"{v}:r{v % 14}:g{gcd(v, 14)}" for v in binders)


def rows() -> list[Row]:
    n6 = 84 * lcm_upto(6)
    return [
        Row("AP tight 1..13", AP, "boundary", "equality_atom"),
        Row("cover tail 12->84", tuple(list(range(1, 12)) + [13, 84]), "O2", "nested_refinement"),
        Row("cover tail 12->168", tuple(list(range(1, 12)) + [13, 168]), "O2", "nested_refinement"),
        Row("divisor-loaded B6", tuple(list(range(1, 12)) + [13, n6]), "A", "raw_time_failure"),
        Row("near/K33 12->36", tuple(list(range(1, 12)) + [13, 36]), "O3", "cross_handoff"),
        Row("P10+K33", tuple(sorted((set(AP) - {10, 12}) | {20, 36})), "O3", "cross_handoff"),
    ]


def print_witness_route() -> None:
    print("[A] WITNESS-ROUTE LEDGER")
    print(
        "row|route|chart|gcd|mult7|mult14|M|tau|mu_direct|components|"
        "largest_arc|D_from_largest|grid_hits_d<=14|binders"
    )
    for row in rows():
        comps = safe_components(row.speeds)
        lengths = [hi - lo for lo, hi in comps]
        measure = sum(lengths, Fraction(0))
        largest = max(lengths, default=Fraction(0))
        m_value, tau = m_exact(row.speeds)
        hits = [d for d in range(2, 15) if has_grid_witness(row.speeds, d)]
        print(
            "|".join(
                [
                    row.name,
                    row.route,
                    row.chart_hint,
                    str(gcd_all(row.speeds)),
                    str(sum(v % 7 == 0 for v in row.speeds)),
                    str(sum(v % 14 == 0 for v in row.speeds)),
                    fmt(m_value),
                    fmt(tau),
                    fmt(measure),
                    str(len(comps)),
                    fmt(largest),
                    str(strict_grid_bound(largest)),
                    ",".join(map(str, hits)) or "none",
                    binder_word(row.speeds, tau, m_value),
                ]
            )
        )


def cap_rows() -> list[tuple[int, Fraction, Fraction, str]]:
    c142 = Fraction(comb(14, 2), 1)
    actual = {
        1: Fraction(comb(13, 2), 91),
        2: Fraction(comb(12, 2), 91),
        3: Fraction(comb(11, 2), 91),
        4: Fraction(1979, 4004),
        5: Fraction(2243, 5880),
    }
    minimizers = {
        1: "{1}",
        2: "{1,13}",
        3: "{1,12,13}",
        4: "{1,11,12,13}",
        5: "{1,5,7,8,9}",
    }
    return [
        (j, Fraction(comb(14 - j, 2), c142), actual[j], minimizers[j])
        for j in range(1, 6)
    ]


def print_pascal_scissors() -> None:
    print("[B] PASCAL-SCISSORS LEDGER")
    print("j|smooth_pair_cap|actual_cap|minimizer|defect|defect_pair_units|status")
    for j, smooth, actual, minimizer in cap_rows():
        defect = smooth - actual
        pair_units = defect * 4004
        status = "pure_pair_mass" if defect == 0 else "scissors_debt"
        print(
            "|".join(
                [
                    str(j),
                    fmt(smooth),
                    fmt(actual),
                    minimizer,
                    fmt(defect),
                    fmt(pair_units),
                    status,
                ]
            )
        )
    print("  4004-readout: j=4 is exactly one affine pair-mass unit of debt.")
    print("  j=5-readout: the defect is no longer a one-pair correction; it asks for S3/S4 or Perron data.")


def print_back_and_forth() -> None:
    print("[C] BACK-AND-FORTH")
    print("  A progress: finite rows have exact direct largest arcs, hence finite bad-denominator budgets.")
    print("  A failure -> B: divisor-loaded B6 already pushes D_from_largest to 5881, so raw time cannot be the invariant; the missing coordinate is the apex/ruler normalization plus the cap/scissors payload.")
    print("  B progress: the cap chart is pure pair mass through j<=3, and j=4 is a one-unit 1/4004 defect.")
    print("  B failure -> A: pair mass does not distinguish positive covering rows from positive K33 cross-handoff rows; the witness ledger must retain branch chart, binders, and endpoint owners.")
    print("  Synthesis: the hidden object is not the measure, the count, or the branch label. It is the overlap map saying which chart reconstructs the coordinate another chart legally forgot.")


def tournament() -> None:
    vertices = [
        ("observer_gluing_packet", (7, 7, 7, 7, 7)),
        ("normalized_arc_chart", (7, 6, 7, 6, 7)),
        ("pascal_scissors_chart", (6, 7, 5, 7, 6)),
        ("branch_k33_chart", (6, 5, 6, 7, 7)),
        ("level7_crt_chart", (6, 6, 7, 5, 6)),
        ("safe_mass_scalar", (4, 4, 4, 3, 3)),
        ("raw_denominator_floor", (3, 2, 3, 2, 1)),
        ("raw_pair_count", (2, 3, 1, 2, 2)),
    ]
    wins = Counter({name: 0 for name, _ in vertices})
    flips = []
    for i, (a, av) in enumerate(vertices):
        for b, bv in vertices[i + 1 :]:
            if av >= bv:
                wins[a] += 1
            elif bv >= av:
                wins[b] += 1
                flips.append((b, a))
            else:
                # Incomparable charts force gluing; break by input path.
                wins[a] += 1
    ordered = sorted((name for name, _ in vertices), key=lambda n: wins[n], reverse=True)
    print("[D] TOURNAMENT_ANALYSIS")
    print("  vertices_are: observer charts and proof obligations, not runners")
    print("  pairwise_observable: predicate retention, denominator-net survival, CRT/lift debt, scissors payload, branch-handoff debt")
    print("  score_histogram:", dict(sorted(Counter(wins.values()).items())))
    print("  hamiltonian_path:", " > ".join(ordered))
    print("  directed_cycles_known: 0 under this lexicographic gauge; incomparability is resolved by chart gluing")
    print("  edge_flips_against_input_order:", len(flips))
    if flips:
        print("  first_flips:", "; ".join(f"{a}>{b}" for a, b in flips[:5]))


def assumption_challenge() -> None:
    print("[E] ASSUMPTION_CHALLENGE")
    print("  considered_vertices: runners, residues, denominators, direct arcs, cap sectors, Pascal entries, endpoint owners, branch handoffs, observer charts, proof obligations")
    print("  selected_vertices: observer charts / proof obligations")
    print("  preserved_predicate: existence of a legal level-1/14 witness or named residual debt")
    print("  destroyed_by_A_alone: apex scale and pair/scissors defect")
    print("  destroyed_by_B_alone: direct denominator-net topology and K33 cross-handoff")


def main() -> None:
    print("S258 LRC14 TWO-FRONTIER OBSERVER GLUING")
    print("=" * 78)
    print_witness_route()
    print()
    print_pascal_scissors()
    print()
    print_back_and_forth()
    print()
    tournament()
    print()
    assumption_challenge()


if __name__ == "__main__":
    main()
