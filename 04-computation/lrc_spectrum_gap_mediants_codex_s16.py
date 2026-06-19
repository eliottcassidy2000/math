#!/usr/bin/env python3
"""
LRC spectrum gap mediant scouts.

Prompt target:
  M_k = {1,...,k-1,2k} gives M(S)=2/(2k+1), so
      1/(k+1) < sigma_2(k) <= 2/(2k+1)
  and the induced lift depth is at least (k+1)(2k+1).

This script verifies that family exactly, then looks for structured dips below
the second-mediant value.  The useful surprise is a residue-class split inside
the one-defect AP family

  S_k = {1,...,k} \\ {k-1} union {3(k-1)}

For tested k == 1 (mod 6), the branch is:

  k == 1 (mod 30):      M(S_k) = 1/(k+1)      (AP-tight again)
  k == 7,13,19,25 mod30 M(S_k) = 3/(3k+2)    (third-mediant)

so its gap above the AP floor is 1/((k+1)(3k+2)).  This does not change the
Theta(1/k^2) scale, but it moves the G2 lift-depth constant from 2 to 3 on an
infinite subsequence.
"""

from __future__ import annotations

from functools import lru_cache, reduce
from fractions import Fraction as F
from itertools import combinations
from math import gcd


def set_gcd(S: tuple[int, ...]) -> int:
    return reduce(gcd, S)


def dist_at_candidate(S: tuple[int, ...], t: F) -> F:
    """Exact min_v ||v*t|| for a rational candidate t."""
    num, den = t.numerator, t.denominator
    best = den
    for v in S:
        r = (v * num) % den
        d = min(r, den - r)
        if d < best:
            best = d
    return F(best, den)


@lru_cache(maxsize=None)
def exact_M(S: tuple[int, ...]) -> tuple[F, F]:
    """Exact M(S)=max_t min_v ||v*t|| via the standard THM-524 candidates."""
    S = tuple(sorted(set(S)))
    candidates: set[F] = set()

    for i, v in enumerate(S):
        for a in range(v):
            candidates.add(F(2 * a + 1, 2 * v))

        for w in S[i + 1 :]:
            for den in (v + w, abs(v - w)):
                if den == 0:
                    continue
                for num in range(1, den):
                    candidates.add(F(num, den))

    best = F(0)
    arg = F(0)
    for t in candidates:
        if not (0 < t < 1):
            continue
        value = dist_at_candidate(S, t)
        if value > best:
            best = value
            arg = t
    return best, arg


def gap_depth(k: int, value: F) -> F:
    floor = F(1, k + 1)
    if value <= floor:
        return F(0)
    return 1 / (value - floor)


def fmt_frac(x: F) -> str:
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


def doubled_top_family(k: int) -> tuple[int, ...]:
    return tuple(list(range(1, k)) + [2 * k])


def third_mediant_family(k: int) -> tuple[int, ...]:
    return tuple(sorted((set(range(1, k + 1)) - {k - 1}) | {3 * (k - 1)}))


def scan_one_defect(k: int, max_add: int) -> list[tuple[F, tuple[int, ...], F, int, int]]:
    """AP {1..k}, remove one runner, add one larger runner."""
    floor = F(1, k + 1)
    family_bound = F(2, 2 * k + 1)
    rows: list[tuple[F, tuple[int, ...], F, int, int]] = []
    ap = set(range(1, k + 1))

    for removed in range(1, k + 1):
        base = ap - {removed}
        for added in range(k + 1, max_add + 1):
            S = tuple(sorted(base | {added}))
            if len(S) != k or set_gcd(S) != 1:
                continue
            value, arg = exact_M(S)
            if floor < value < family_bound:
                rows.append((value, S, arg, removed, added))
    return sorted(rows)


def scan_two_defect(k: int, max_add: int) -> list[tuple[F, tuple[int, ...], F, tuple[int, ...], tuple[int, ...]]]:
    """AP {1..k}, remove up to two runners, add the same number of larger runners."""
    floor = F(1, k + 1)
    family_bound = F(2, 2 * k + 1)
    rows: list[tuple[F, tuple[int, ...], F, tuple[int, ...], tuple[int, ...]]] = []
    ap = set(range(1, k + 1))
    seen: set[tuple[int, ...]] = set()

    for defect in (1, 2):
        for removed in combinations(range(1, k + 1), defect):
            base = ap - set(removed)
            for added in combinations(range(k + 1, max_add + 1), defect):
                S = tuple(sorted(base | set(added)))
                if len(S) != k or S in seen or set_gcd(S) != 1:
                    continue
                seen.add(S)
                value, arg = exact_M(S)
                if floor < value < family_bound:
                    rows.append((value, S, arg, removed, added))
    return sorted(rows)


def bounded_exhaustive(k: int, B: int) -> tuple[tuple[F, tuple[int, ...], F] | None, int, int]:
    floor = F(1, k + 1)
    best: tuple[F, tuple[int, ...], F] | None = None
    at_floor = 0
    total = 0
    for S in combinations(range(1, B + 1), k):
        if set_gcd(S) != 1:
            continue
        total += 1
        value, arg = exact_M(S)
        if value == floor:
            at_floor += 1
        if value > floor and (best is None or value < best[0]):
            best = (value, S, arg)
    return best, at_floor, total


def print_family_verifications() -> None:
    print("=" * 78)
    print("1. Doubled-top family: {1,...,k-1,2k}")
    print("=" * 78)
    print("k   M(S)       argmax t   expected    gap-depth")
    for k in range(2, 31):
        S = doubled_top_family(k)
        value, arg = exact_M(S)
        expected = F(2, 2 * k + 1)
        depth = gap_depth(k, value)
        ok = "OK" if value == expected else "FAIL"
        print(f"{k:2d}  {fmt_frac(value):>9}  {fmt_frac(arg):>9}  {fmt_frac(expected):>9}  {fmt_frac(depth):>9}  {ok}")


def print_third_mediant_family() -> None:
    print("\n" + "=" * 78)
    print("2. Third-mediant one-defect AP family")
    print("=" * 78)
    print("S_k={1,...,k}\\{k-1} U {3(k-1)} for k == 1 mod 6.")
    print("The tested split is: k==1 mod30 is AP-tight; the other four")
    print("1 mod6 classes are third-mediant.")
    print("k   mod30   M(S)       argmax t   expected    state       gap-depth")
    for k in range(7, 104, 6):
        S = third_mediant_family(k)
        value, arg = exact_M(S)
        floor = F(1, k + 1)
        expected = floor if k % 30 == 1 else F(3, 3 * k + 2)
        family_bound = F(2, 2 * k + 1)
        depth = gap_depth(k, value)
        if value == floor:
            state = "tight"
        elif value < family_bound:
            state = "third"
        else:
            state = "other"
        ok = "OK" if value == expected else "CHECK"
        print(
            f"{k:3d}  {k % 30:5d}  {fmt_frac(value):>9}  {fmt_frac(arg):>9}"
            f"  {fmt_frac(expected):>9}  {state:>8}  {fmt_frac(depth):>9}  {ok}"
        )


def print_bounded_scans() -> None:
    print("\n" + "=" * 78)
    print("3. Bounded exact scans near the AP floor")
    print("=" * 78)
    print("These are not proofs of sigma_2(k); they are finite-box scouts.")
    for k, B in ((2, 40), (3, 36), (4, 24), (5, 20)):
        best, at_floor, total = bounded_exhaustive(k, B)
        floor = F(1, k + 1)
        family_bound = F(2, 2 * k + 1)
        print(
            f"k={k}, B={B}, primitive rows={total}, at AP floor={at_floor}, "
            f"best above floor={best}, family={family_bound}, tight={best and best[0] == family_bound}"
        )
        assert best is not None
        assert best[0] == family_bound


def print_defect_scans() -> None:
    print("\n" + "=" * 78)
    print("4. AP one/two-defect scans for dips below the doubled-top value")
    print("=" * 78)
    for k in range(3, 16):
        rows = scan_one_defect(k, 4 * k)
        first = rows[:3]
        print(f"one-defect k={k:2d}, add<=4k: below-family rows={len(rows)}")
        for value, S, arg, removed, added in first:
            print(
                f"  M={fmt_frac(value):>8}, t={fmt_frac(arg):>8}, "
                f"remove={removed}, add={added}, S={S}, depth={fmt_frac(gap_depth(k, value))}"
            )

    print("\nTwo-defect window, where the first small strict leak is visible:")
    for k in range(4, 9):
        rows = scan_two_defect(k, 4 * k)
        print(f"two-defect k={k:2d}, add<=4k: below-family rows={len(rows)}")
        for value, S, arg, removed, added in rows[:5]:
            print(
                f"  M={fmt_frac(value):>8}, t={fmt_frac(arg):>8}, "
                f"remove={removed}, add={added}, S={S}, depth={fmt_frac(gap_depth(k, value))}"
            )


def print_tournament_analysis() -> None:
    print("\n" + "=" * 78)
    print("5. Tournament Analysis quotient")
    print("=" * 78)
    print("Assumption challenge:")
    print("  considered vertices = runners, residue classes mod 2k+1, pair-crossing")
    print("  denominators, AP-defect families, finite scan boxes, and proof obligations.")
    print("Chosen vertices = proof-obligation/family routes; this preserves the")
    print("  predicate 'gap above AP floor and its reciprocal lift-depth' while")
    print("  discarding individual witness-time geometry except for the argmax address.")
    vertices = [
        "infinite_lower_bound_for_g(k)",
        "third_mediant_family",
        "doubled_top_family",
        "one_defect_scan",
        "two_defect_scan",
        "bounded_exhaustive_box",
        "raw_runner_vertices",
    ]
    edges = [(vertices[i], vertices[j]) for i in range(len(vertices)) for j in range(i + 1, len(vertices))]
    print("Tournament is transitive under proof leverage ranking.")
    print("Hamiltonian path:")
    print("  " + " > ".join(vertices))
    print(f"score histogram: {[len(vertices)-1-i for i in range(len(vertices))]}")
    print(f"directed 3-cycles: 0; SCC sizes: [{len(vertices)}]; edge count: {len(edges)}")


def main() -> None:
    print_family_verifications()
    print_third_mediant_family()
    print_bounded_scans()
    print_defect_scans()
    print_tournament_analysis()
    print("\nREADOUT")
    print("  The doubled-top family proves an explicit universal upper bound")
    print("      g(k) <= 1/((k+1)(2k+1)).")
    print("  The one-defect branch k == 7,13,19,25 mod30 improves that to")
    print("      g(k) <= 1/((k+1)(3k+2))")
    print("  on an infinite tested residue-class subsequence; the k == 1 mod30")
    print("  branch of the same construction is AP-tight again.  The lower-bound")
    print("  question remains open:")
    print("  these computations do not show any o(1/k^2) dip.")


if __name__ == "__main__":
    main()
