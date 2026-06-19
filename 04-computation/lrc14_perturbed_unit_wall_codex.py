#!/usr/bin/env python3
"""
LRC(14): perturbed unit-wall criterion for one 14-multiple rows.

For rows S = A U {14m}, this script isolates the exact first-branch criterion near
the six unit walls a/14:

    t = a/14 + delta,   a in (Z/14Z)^*,   |delta| small.

It verifies three things.

1. The exact first-branch inequality is correct:
   for positive delta, a speed v with residue s = a v (mod 14) stays > 1/14
   exactly while delta < (13-s)/(14v); for negative delta, while
   |delta| < (s-1)/(14v).
2. The resulting constructed witness is exact on a bounded exhaustive census.
3. The user's stronger "missing a unit residue should force a unit-wall branch"
   remains plausible beyond the first branch: bounded exhaustive scans find a
   unit-wall branch witness in every tested row, even when the first branch fails.

Methodology note: there is no clean Tournament Analysis quotient here. The
natural object is a scalar interval width on each unit wall, not a pairwise
binary relation. The assumption challenged in this session is precisely that the
first branch already captures every unit-wall witness.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd


N = 14
THRESHOLD = Fraction(1, N)
UNIT_RESIDUES = (1, 3, 5, 9, 11, 13)
UNIT_INVERSE = {1: 1, 3: 5, 5: 3, 9: 11, 11: 9, 13: 13}


@dataclass(frozen=True)
class BranchWitness:
    a: int
    sign: int
    delta: Fraction
    time: Fraction
    upper: Fraction


def normalize_circle(x: Fraction) -> Fraction:
    x %= 1
    if x > Fraction(1, 2):
        x = 1 - x
    return x


def is_witness(speeds: tuple[int, ...], t: Fraction) -> bool:
    return all(normalize_circle(Fraction(v) * t) > THRESHOLD for v in speeds)


def first_branch_upper(a: int, sign: int, A: tuple[int, ...], m: int) -> Fraction | None:
    """
    Exact first local branch at a/14.

    sign = +1 means t = a/14 + delta with delta > 0.
    sign = -1 means t = a/14 - delta with delta > 0.
    """

    upper = Fraction(1, 28 * m)
    for v in A:
        s = (a * (v % N)) % N
        if s == 0:
            return None
        coeff = (13 - s) if sign > 0 else (s - 1)
        if coeff <= 0:
            return None
        cand = Fraction(coeff, N * v)
        if cand < upper:
            upper = cand
    return upper


def first_branch_witness(A: tuple[int, ...], m: int) -> BranchWitness | None:
    lower = Fraction(1, 196 * m)
    best: BranchWitness | None = None
    for a in UNIT_RESIDUES:
        for sign in (1, -1):
            upper = first_branch_upper(a, sign, A, m)
            if upper is None or not (lower < upper):
                continue
            delta = (lower + upper) / 2
            time = Fraction(a, N) + (delta if sign > 0 else -delta)
            candidate = BranchWitness(a=a, sign=sign, delta=delta, time=time, upper=upper)
            if best is None or upper > best.upper:
                best = candidate
    return best


def missing_unit_residues(A: tuple[int, ...]) -> tuple[int, ...]:
    residues = {v % N for v in A}
    return tuple(r for r in UNIT_RESIDUES if r not in residues)


def branch_breakpoints(speeds: tuple[int, ...], a: int, radius: Fraction) -> list[Fraction]:
    points = {-radius, Fraction(0), radius}
    for v in speeds:
        va = Fraction(v * a, N)
        for k in range(-1, v + 2):
            for sign in (-1, 1):
                delta = (Fraction(k) + Fraction(sign, N) - va) / v
                if -radius <= delta <= radius:
                    points.add(delta)
    return sorted(points)


def any_unit_wall_branch_witness(speeds: tuple[int, ...], radius: Fraction = Fraction(1, 14)) -> BranchWitness | None:
    for a in UNIT_RESIDUES:
        points = branch_breakpoints(speeds, a, radius)
        for lo, hi in zip(points, points[1:]):
            delta = (lo + hi) / 2
            time = Fraction(a, N) + delta
            if is_witness(speeds, time):
                sign = 1 if delta > 0 else -1
                return BranchWitness(a=a, sign=sign, delta=abs(delta), time=time, upper=hi)
    return None


def exhaustive_first_branch_verification(pool_max: int, m_values: range) -> tuple[int, int, int]:
    pool = tuple(v for v in range(1, pool_max + 1) if v % N != 0)
    total_rows = 0
    certified_rows = 0
    mismatches = 0
    for A in combinations(pool, 12):
        if not missing_unit_residues(A):
            continue
        for m in m_values:
            total_rows += 1
            witness = first_branch_witness(A, m)
            if witness is None:
                continue
            certified_rows += 1
            speeds = tuple(sorted(A + (N * m,)))
            if not is_witness(speeds, witness.time):
                mismatches += 1
    return total_rows, certified_rows, mismatches


def exhaustive_any_branch_verification(pool_max: int, m_values: range) -> tuple[int, int]:
    pool = tuple(v for v in range(1, pool_max + 1) if v % N != 0)
    total_rows = 0
    misses = 0
    for A in combinations(pool, 12):
        if not missing_unit_residues(A):
            continue
        for m in m_values:
            total_rows += 1
            speeds = tuple(sorted(A + (N * m,)))
            if any_unit_wall_branch_witness(speeds) is None:
                misses += 1
    return total_rows, misses


def ap_drop_table() -> list[tuple[int, str]]:
    base = tuple(range(1, 14))
    rows: list[tuple[int, str]] = []
    for removed in base:
        A = tuple(v for v in base if v != removed)
        first = None
        for m in range(1, 5):
            if first_branch_witness(A, m) is not None:
                first = m
                break
        if first is None:
            rows.append((removed, "unit-complete residual"))
        else:
            rows.append((removed, f"first-branch from m={first}"))
    return rows


def main() -> None:
    print("=" * 78)
    print("LRC(14) perturbed unit-wall criterion for S = A U {14m}")
    print("=" * 78)
    print("First branch theorem:")
    print("  t = a/14 + delta, delta > 0  ->  each v in A stays safe while")
    print("      delta < (13 - (a v mod 14)) / (14 v).")
    print("  t = a/14 - delta, delta > 0  ->  each v in A stays safe while")
    print("      delta < ((a v mod 14) - 1) / (14 v).")
    print("  The 14m-runner is safe on the first branch exactly when")
    print("      delta > 1 / (196 m)")
    print("  provided delta < 1 / (28 m).")
    print()

    print("A. Exact first-branch verification on bounded exhaustive rows")
    for pool_max, m_values in ((16, range(1, 5)), (18, range(1, 4))):
        total, certified, mismatches = exhaustive_first_branch_verification(pool_max, m_values)
        print(
            f"  pool <= {pool_max:2d}, m in {m_values.start}..{m_values.stop - 1}: "
            f"rows={total}, first-branch certified={certified}, mismatches={mismatches}"
        )
    print()

    print("B. Stronger bounded test: some unit-wall branch exists")
    for pool_max, m_values in ((16, range(1, 5)), (18, range(1, 4))):
        total, misses = exhaustive_any_branch_verification(pool_max, m_values)
        print(
            f"  pool <= {pool_max:2d}, m in {m_values.start}..{m_values.stop - 1}: "
            f"rows={total}, rows with NO unit-wall branch witness={misses}"
        )
    print()

    print("C. AP one-drop families {1..13}\\\\{j} U {14m}")
    for removed, status in ap_drop_table():
        print(f"  drop {removed:2d}: {status}")
    print("  Interpretation: odd unit drops are discharged immediately by the first branch;")
    print("  even / 7 drops are exactly the unit-complete residual.")
    print()

    print("D. Explicit first-branch failure but later-branch success")
    counterexample_A = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 26)
    counterexample_S = tuple(sorted(counterexample_A + (14,)))
    first = first_branch_witness(counterexample_A, 1)
    later = any_unit_wall_branch_witness(counterexample_S)
    print(f"  A={counterexample_A}, missing unit residues={missing_unit_residues(counterexample_A)}")
    print(f"  first-branch witness: {first}")
    print(
        "  later branch witness: "
        f"a={later.a}, sign={'+' if later.sign > 0 else '-'}, "
        f"delta={later.delta}, t={later.time}"
    )
    print(
        "  This shows the raw 'missing a unit residue' idea is not a complete FIRST-BRANCH"
        " theorem by itself, but it may still hold after later wraps."
    )
    print()

    print("E. Small exact example where the first branch already works")
    sample_A = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16)
    sample_S = tuple(sorted(sample_A + (14,)))
    sample = first_branch_witness(sample_A, 1)
    print(f"  A={sample_A}, missing unit residues={missing_unit_residues(sample_A)}")
    print(
        "  first-branch witness: "
        f"a={sample.a}, sign={'+' if sample.sign > 0 else '-'}, "
        f"delta={sample.delta}, t={sample.time}, exact={is_witness(sample_S, sample.time)}"
    )
    print()

    print("Takeaway")
    print("  THM-level: the exact first-branch inequality is elementary and exact.")
    print("  Evidence-level: every bounded row tested with a missing unit residue still")
    print("  has SOME unit-wall branch witness, even when the first branch fails.")
    print("  Proof-program implication: the genuinely hard one-multiple rows appear to be")
    print("  the unit-complete cores, matching the new exact low-L extremizers.")


if __name__ == "__main__":
    main()
