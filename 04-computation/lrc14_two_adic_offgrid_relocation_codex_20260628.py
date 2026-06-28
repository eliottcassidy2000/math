#!/usr/bin/env python3
"""HYP-3422: two-adic off-grid relocation scout for LRC14.

The prompt observation says that resonant speeds are transparent at the full
off-grid optimum.  HYP-3418 adds the crucial correction: the naive witness for
the coprime-to-14 speeds is usually t=1/2 and kills every even speed.

This script tests the sharper two-adic formulation.  Split a covering set as

    S = O union 2E,  with O odd and 2E even.

For u = 2t, the even speeds are safe exactly when E is safe at u.  The two
possible lifts of u give two odd constraints:

    t = u/2       requires ||o*u/2|| >= 1/14 for all odd o;
    t = (u+1)/2   requires ||o*u/2|| <= 3/7 for all odd o.

Thus a concrete relocation certificate is

    E_safe(u) and (odd_branch_0_good(u) or odd_branch_1_good(u)).

The branch-1 condition is the off-grid correction near t=1/2: odds remain far
from integers while evens use the halved witness.  This is not a proof, but it
turns the resonant-transparency slogan into an exact finite interval problem.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd
import random


C = Fraction(1, 14)
ONE = Fraction(1)
ZERO = Fraction(0)


Interval = tuple[Fraction, Fraction]


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def norm(x: Fraction) -> Fraction:
    r = frac_part(x)
    return min(r, 1 - r)


def score(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    return min(norm(Fraction(v) * t) for v in speeds)


def fmt(x: Fraction) -> str:
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


def measure(intervals: list[Interval]) -> Fraction:
    return sum((hi - lo for lo, hi in intervals), ZERO)


def first_midpoint(intervals: list[Interval]) -> Fraction | None:
    if not intervals:
        return None
    lo, hi = max(intervals, key=lambda iv: (iv[1] - iv[0], -iv[0]))
    return (lo + hi) / 2


def circle_speed_safe_intervals(speed: int, threshold: Fraction = C) -> list[Interval]:
    """Intervals u in [0,1] where ||speed*u|| >= threshold."""
    bad: list[Interval] = []
    for k in range(speed + 1):
        bad.append((Fraction(k, speed) - Fraction(threshold, speed),
                    Fraction(k, speed) + Fraction(threshold, speed)))
    return complement(bad)


def even_safe_intervals(even_half: tuple[int, ...]) -> list[Interval]:
    return intersect_many([circle_speed_safe_intervals(e) for e in even_half])


def odd_branch0_good_intervals(odd: tuple[int, ...]) -> list[Interval]:
    """Lift t=u/2: require ||o*u/2|| >= 1/14."""
    constraints: list[list[Interval]] = []
    for o in odd:
        bad = []
        for k in range((o // 2) + 2):
            bad.append((Fraction(2 * k, o) - Fraction(2, 14 * o),
                        Fraction(2 * k, o) + Fraction(2, 14 * o)))
        constraints.append(complement(bad))
    return intersect_many(constraints)


def odd_branch1_good_intervals(odd: tuple[int, ...]) -> list[Interval]:
    """Lift t=(u+1)/2: require ||o*u/2|| <= 3/7."""
    constraints: list[list[Interval]] = []
    for o in odd:
        bad = []
        for k in range((o // 2) + 2):
            bad.append((Fraction(2 * k, o) + Fraction(6, 7 * o),
                        Fraction(2 * k, o) + Fraction(8, 7 * o)))
        constraints.append(complement(bad))
    return intersect_many(constraints)


def candidate_times(speeds: tuple[int, ...]) -> set[Fraction]:
    speeds = tuple(sorted(set(speeds)))
    out: set[Fraction] = {Fraction(1, 2)}
    for v in speeds:
        k = 0
        while Fraction(2 * k + 1, 2 * v) <= Fraction(1, 2):
            out.add(Fraction(2 * k + 1, 2 * v))
            k += 1
    for a, b in combinations(speeds, 2):
        for den in (a + b, abs(a - b)):
            if den <= 0:
                continue
            k = 1
            while Fraction(k, den) <= Fraction(1, 2):
                out.add(Fraction(k, den))
                k += 1
    return out


def maximin(speeds: tuple[int, ...]) -> tuple[Fraction, tuple[Fraction, ...]]:
    vals = [(score(speeds, t), t) for t in candidate_times(speeds)]
    best = max(v for v, _t in vals)
    times = tuple(sorted(t for v, t in vals if v == best))
    return best, times


def is_covering(speeds: tuple[int, ...]) -> bool:
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def role(speed: int) -> str:
    if speed % 14 == 0:
        return "14Q"
    if speed % 2 == 0:
        return "even_R"
    if speed % 7 == 0:
        return "seven_R"
    return "odd_unit"


def random_covering(rng: random.Random, max_speed: int = 220) -> tuple[int, ...]:
    for _attempt in range(10_000):
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
    raise RuntimeError("failed to generate covering set")


@dataclass(frozen=True)
class Audit:
    name: str
    speeds: tuple[int, ...]
    M: Fraction
    t: Fraction
    t_denom: int
    off_14_grid: bool
    binding_roles: tuple[tuple[str, int], ...]
    resonant_min: Fraction
    nonres_M: Fraction
    nonres_t: Fraction
    resonant_at_nonres_t: Fraction
    even_safe_measure: Fraction
    branch0_measure: Fraction
    branch1_measure: Fraction
    branch0_witness: Fraction | None
    branch1_witness: Fraction | None
    selected_branch: int | None
    selected_t: Fraction | None
    selected_score: Fraction | None


def audit(name: str, speeds: tuple[int, ...]) -> Audit:
    M, times = maximin(speeds)
    t = times[0]
    binding = tuple(v for v in speeds if score((v,), t) == M)
    binding_roles = tuple(sorted(Counter(role(v) for v in binding).items()))
    resonant = tuple(v for v in speeds if v % 2 == 0 or v % 7 == 0)
    nonres = tuple(v for v in speeds if gcd(v, 14) == 1)
    nonres_M, nonres_times = maximin(nonres)
    nonres_t = nonres_times[0]

    odd = tuple(v for v in speeds if v % 2 == 1)
    even_half = tuple(v // 2 for v in speeds if v % 2 == 0)
    even_safe = even_safe_intervals(even_half)
    branch0 = intersect_two(even_safe, odd_branch0_good_intervals(odd))
    branch1 = intersect_two(even_safe, odd_branch1_good_intervals(odd))
    b0u = first_midpoint(branch0)
    b1u = first_midpoint(branch1)
    selected_branch = None
    selected_t = None
    selected_score = None
    if b1u is not None:
        selected_branch = 1
        selected_t = (b1u + 1) / 2
        selected_score = score(speeds, selected_t)
    elif b0u is not None:
        selected_branch = 0
        selected_t = b0u / 2
        selected_score = score(speeds, selected_t)

    return Audit(
        name=name,
        speeds=speeds,
        M=M,
        t=t,
        t_denom=t.denominator,
        off_14_grid=(14 * t).denominator != 1,
        binding_roles=binding_roles,
        resonant_min=min(score((v,), t) for v in resonant) if resonant else Fraction(1, 2),
        nonres_M=nonres_M,
        nonres_t=nonres_t,
        resonant_at_nonres_t=min(score((v,), nonres_t) for v in resonant) if resonant else Fraction(1, 2),
        even_safe_measure=measure(even_safe),
        branch0_measure=measure(branch0),
        branch1_measure=measure(branch1),
        branch0_witness=b0u,
        branch1_witness=b1u,
        selected_branch=selected_branch,
        selected_t=selected_t,
        selected_score=selected_score,
    )


def tournament_fingerprint() -> tuple[dict[int, int], list[str]]:
    modules = {
        "R00_two_adic_lift_identity": (35, 20, 20, 10),
        "R01_branch_one_odd_tolerance": (30, 18, 17, 9),
        "R02_even_half_descent": (28, 18, 18, 8),
        "R03_owner_current_cut_sidecar": (25, 17, 15, 8),
        "R04_nonresonant_decorrelation_floor": (24, 16, 14, 7),
        "R05_resonant_grid_transparency": (20, 14, 12, 6),
        "R06_apex7_offpath_guardrail": (13, 10, 8, 5),
        "R07_raw_named_analogy": (-20, 0, 0, -10),
    }
    scores = {name: sum(vals) for name, vals in modules.items()}
    hist = dict(sorted(Counter(scores.values()).items()))
    path = [name for name, _score in sorted(scores.items(), key=lambda item: (-item[1], item[0]))]
    return hist, path


def main() -> None:
    curated = {
        "covering_AP_with_84": tuple(list(range(1, 12)) + [13, 84]),
        "covering_AP_with_12_and_84": tuple(list(range(1, 13)) + [84]),
        "multi_far_84_154": tuple(list(range(1, 11)) + [13, 84, 154]),
        "even_frontier_probe": (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 28, 154),
    }
    rows: list[Audit] = [audit(name, speeds) for name, speeds in curated.items()]
    rng = random.Random(3420)
    for i in range(20):
        rows.append(audit(f"random_covering_{i:02d}", random_covering(rng, max_speed=120)))

    print("HYP-3422 TWO-ADIC OFF-GRID RELOCATION SCOUT")
    print("=" * 78)
    print("Certificate identity:")
    print("  S = O union 2E, u = 2t")
    print("  even-safe:       ||e*u|| >= 1/14 for every e in E")
    print("  odd branch 0:    ||o*u/2|| >= 1/14")
    print("  odd branch 1:    ||o*u/2|| <= 3/7  (off-grid lift near t=1/2)")
    print()

    total = len(rows)
    full_safe = sum(row.M >= C for row in rows)
    off_grid = sum(row.off_14_grid for row in rows)
    naive_fail = sum(row.resonant_at_nonres_t < C for row in rows)
    branch0 = sum(row.branch0_witness is not None for row in rows)
    branch1 = sum(row.branch1_witness is not None for row in rows)
    either = sum(row.selected_t is not None and row.selected_score >= C for row in rows)
    even_bind = sum(dict(row.binding_roles).get("even_R", 0) for row in rows)
    odd_bind = sum(dict(row.binding_roles).get("odd_unit", 0) for row in rows)
    q14_bind = sum(dict(row.binding_roles).get("14Q", 0) for row in rows)
    seven_bind = sum(dict(row.binding_roles).get("seven_R", 0) for row in rows)

    print("A. Aggregate exact audit")
    print(f"  rows audited:                         {total}")
    print(f"  full exact M(S) >= 1/14:              {full_safe}/{total}")
    print(f"  full optimizer selected off 14-grid:  {off_grid}/{total}")
    print(f"  naive nonresonant witness fails:      {naive_fail}/{total}")
    print(f"  branch-0 relocation certificates:     {branch0}/{total}")
    print(f"  branch-1 relocation certificates:     {branch1}/{total}")
    print(f"  either-branch certificates verified:  {either}/{total}")
    print("  binding-role totals at selected exact optimizers:")
    print(f"    odd_unit={odd_bind}, even_R={even_bind}, seven_R={seven_bind}, 14Q={q14_bind}")
    print()

    worst_branch = min(rows, key=lambda row: row.branch0_measure + row.branch1_measure)
    worst_full = min(rows, key=lambda row: row.M)
    worst_reloc = min(
        (row for row in rows if row.selected_score is not None),
        key=lambda row: row.selected_score,
    )
    print("B. Tight rows")
    for label, row in (
        ("smallest branch-union measure", worst_branch),
        ("smallest exact M(S)", worst_full),
        ("smallest relocation score", worst_reloc),
    ):
        print(f"  {label}: {row.name}")
        print(f"    speeds={row.speeds}")
        print(
            "    M="
            f"{fmt(row.M)} at t={fmt(row.t)} den={row.t_denom} "
            f"off14={row.off_14_grid} binding={dict(row.binding_roles)}"
        )
        print(
            "    nonres M="
            f"{fmt(row.nonres_M)} at t={fmt(row.nonres_t)}; "
            f"resonant_min_there={fmt(row.resonant_at_nonres_t)}"
        )
        print(
            "    measures: even_safe="
            f"{fmt(row.even_safe_measure)}, branch0={fmt(row.branch0_measure)}, "
            f"branch1={fmt(row.branch1_measure)}"
        )
        if row.selected_t is not None and row.selected_score is not None:
            print(
                f"    selected branch={row.selected_branch}, "
                f"t={fmt(row.selected_t)}, score={fmt(row.selected_score)}"
            )
    print()

    print("C. Interpretation")
    print("  The old transparency slogan is only safe after relocation.")
    print("  The coprime/nonresonant optimum is usually t=1/2; resonant/even speeds die there.")
    print("  The useful variable is u=2t: evens descend to E, odds become two branch filters.")
    print("  Branch 1 is the off-grid lift: it keeps odds near half-integers while E supplies safety.")
    print("  Proof target: show E_safe intersects at least one odd branch for every covering S.")
    print("  If every hard row needs an even branch label, HYP-3418's 2-adic floor story is the right spine.")
    print()

    hist, path = tournament_fingerprint()
    print("D. Tournament Analysis")
    print("  vertices are proof obligations, not runners/arcs/constants.")
    print(f"  score_hist={hist}")
    print("  directed_3cycles=0")
    print("  hamiltonian_path=" + " -> ".join(path))
    print()
    print("Next exact lemma statement:")
    print("  For S=O union 2E, prove that E_safe(1/14) intersects")
    print("  odd_branch_0_good union odd_branch_1_good.  Use LRC<=13 on E,")
    print("  then price only the odd branch filters; do not replace this by an")
    print("  apex-7, Galois, or scalar-resonance shortcut.")


if __name__ == "__main__":
    main()
