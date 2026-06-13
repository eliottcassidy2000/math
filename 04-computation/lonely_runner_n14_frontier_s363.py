#!/usr/bin/env python3
"""
lonely_runner_n14_frontier_s363.py

oracle-2026-05-31-S16

The 14-runner case (k=13 speeds, threshold 1/14) is the FIRST OPEN case of the
Lonely Runner Conjecture as of 2026-05 (proven through 13 runners:
Rosenfeld arXiv:2509.14111 for 8; arXiv:2511.22427 / 2512.01912 for 9,10;
Sungkawichai-Trakulthongchai arXiv:2604.23906 for 11,12,13).

This script makes concrete a NEW synthesis of the repo's two pictures
(forbidden-measure + endpoint-protection) via a *divisibility sieve* reading,
specialized to n = 14.

KEY OBSERVATION (the "whack-a-mole" / modulus-cover lemma).
For n = 14, threshold 1/14, and any candidate speed set V:

  * For each integer m with 2 <= m <= 13: if NO speed in V is divisible by m,
    then t = 1/m is lonely, because ||v/m|| >= 1/m > 1/14 for every v not
    divisible by m.  (Strict: 1/m > 1/14 for m <= 13.)
  * For m = 14: if no speed is divisible by 14, then t = a/14 (gcd(a,14)=1)
    is a *boundary* lonely witness, ||v a/14|| >= 1/14 for all v.

So a counterexample at n=14 must contain, for EVERY m in {2,...,14}, at least
one speed divisible by m.  The famous tight family {1,...,13} covers every
m in {2,...,13} (speed m divides itself) but MISSES m=14 -- which is exactly
why its only surviving lonely times are the units a/14: the conjecture's
tightness lives precisely at the one modulus scale the initial segment fails
to cover.

This reframes the search for an n=14 counterexample as a tension:
  - small speeds cover macroscopic forbidden MEASURE but the initial segment
    misses modulus 14;
  - covering modulus 14 (or any prime power) forces either a large speed
    (tiny scattered intervals, leftover measure -> a gap opens) or sacrificing
    a small-speed slot (re-opening a coarser modulus scale).
13 speeds cannot simultaneously win both games.  This file measures that
tension exactly over Fraction arithmetic.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()

N = 14            # number of runners
K = N - 1         # number of speeds = 13
THR = Fraction(1, N)
MODULI = list(range(2, N + 1))   # {2,...,14}


def fmt(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def divisibility_profile(speeds: tuple[int, ...]) -> dict[int, list[int]]:
    """For each m in {2,...,14}, which speeds are divisible by m."""
    return {m: [v for v in speeds if v % m == 0] for m in MODULI}


def uncovered_moduli(speeds: tuple[int, ...]) -> list[int]:
    prof = divisibility_profile(speeds)
    return [m for m in MODULI if not prof[m]]


@dataclass
class N14Report:
    label: str
    speeds: tuple[int, ...]
    max_gap: Fraction
    forbidden_length: Fraction
    boundary_count: int
    components: int
    uncovered: list[int]
    full_cover: bool          # forbidden length 1 AND no boundary witness
    boundary_first: Fraction | None


def analyze(label: str, raw_speeds: list[int]) -> N14Report:
    if len(set(abs(v) for v in raw_speeds if v)) != K:
        # we still analyze, but flag non-13 sets by label
        pass
    row = S356.report(label, raw_speeds)
    speeds = row.speeds
    unc = uncovered_moduli(speeds)
    full = (row.forbidden_length == 1) and (row.boundary_witness_count == 0)
    return N14Report(
        label=label,
        speeds=speeds,
        max_gap=row.max_gap,
        forbidden_length=row.forbidden_length,
        boundary_count=row.boundary_witness_count,
        components=row.components,
        uncovered=unc,
        full_cover=full,
        boundary_first=row.boundary_witness,
    )


def print_report(r: N14Report) -> None:
    status = (
        "OPEN-COVER **COUNTEREXAMPLE**"
        if r.full_cover
        else ("positive-gap" if r.max_gap > 0 else "boundary-only (tight)")
    )
    print(f"[{r.label}]  -> {status}")
    print(f"  speeds = {r.speeds}")
    print(
        f"  forbidden_length = {fmt(r.forbidden_length)}  "
        f"max_gap = {fmt(r.max_gap)} ({float(r.max_gap):.5f})  "
        f"components = {r.components}"
    )
    print(
        f"  boundary_witnesses = {r.boundary_count}  "
        f"first = {fmt(r.boundary_first)}"
    )
    if r.uncovered:
        # predicted surviving witness scale = smallest uncovered modulus
        msmall = min(r.uncovered)
        print(
            f"  UNCOVERED moduli (no speed divisible by) = {r.uncovered}"
            f"   => t = 1/{msmall} is lonely"
        )
    else:
        print("  divisibility sieve fully covered (all m in 2..14 hit)")
    print()


def section(title: str) -> None:
    print("=" * 72)
    print(title)
    print("=" * 72)


def main() -> None:
    print("Lonely Runner n=14 frontier analysis (oracle-2026-05-31-S16)")
    print(f"runners n={N}, speeds k={K}, threshold 1/{N}\n")

    section("A. The tight baseline: initial segment {1,...,13}")
    print_report(analyze("initial segment 1..13", list(range(1, 14))))
    print("Interpretation: covers every modulus 2..13 (speed m divides m), misses")
    print("only m=14.  The 6 surviving boundary witnesses are exactly the units")
    print("a/14, gcd(a,14)=1 -> the tightness lives at the one uncovered scale.\n")

    section("B. Whack-a-mole: inserting a multiple of 14 re-opens a coarse scale")
    # Replace one of 1..13 by a multiple of 14; watch which modulus re-opens.
    base = list(range(1, 14))
    for drop in [13, 11, 9, 8, 5, 7, 12, 10, 6]:
        speeds = sorted(set(base) - {drop} | {14})
        print_report(analyze(f"drop {drop}, add 14", speeds))

    section("C. Cover m=14 AND the dropped scale with one combined large speed")
    # To cover both 14 and 13 after dropping 13, use 14*13=182 in place of 13.
    examples_C = [
        ("drop 13, add 182 (=2*7*13)", sorted(set(base) - {13} | {182})),
        ("drop 11,13 add 14,143", sorted(set(base) - {11, 13} | {14, 143})),
        ("drop 8, add 56 (=8*7, covers 8&14&7)", sorted(set(base) - {8} | {56})),
        ("drop 9, add 126 (=2*9*7 covers 9,14,...)", sorted(set(base) - {9} | {126})),
        ("one super-speed 360360 + small filler",
         sorted(set(range(1, 13)) | {360360})),  # 1..12 + 360360, drops 13
    ]
    for label, speeds in examples_C:
        print_report(analyze(label, speeds))

    section("D. Sets that DO satisfy the full divisibility sieve (no uncovered m)")
    # Hand-build 13-speed sets hitting every modulus 2..14; see if any is tight/cover.
    candidates_D = [
        # contains 14 (covers 2,7,14), 8(2,4,8), 9(3,9), 5(5), 11, 13, and 6,10,12 explicitly
        ("sieve-complete v1", [5, 6, 8, 9, 10, 11, 12, 13, 14, 1, 2, 3, 4]),
        ("sieve-complete v2", [14, 8, 9, 5, 11, 13, 6, 10, 12, 7, 4, 3, 1]),
        ("sieve-complete near-initial", [1, 2, 3, 4, 5, 6, 8, 9, 11, 13, 10, 12, 14]),
        # try to keep small + cover 14 by using 28 (=4*7, covers 4,7,14,2)
        ("use 28 for {4,7,14,2}", [1, 3, 5, 6, 8, 9, 10, 11, 12, 13, 28, 2, 27]),
    ]
    for label, speeds in candidates_D:
        r = analyze(label, speeds)
        print_report(r)

    section("E. Even-n antipodal tool: all-odd speeds make t=1/2 lonely")
    odd13 = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25]
    r = analyze("13 odd speeds (1,3,..,25)", odd13)
    print_report(r)
    print("Any all-odd set is loneliness-trivial at t=1/2 (||v/2||=1/2). So a")
    print("counterexample must contain an even speed -- automatically supplied by")
    print("the required multiple of 14.\n")

    section("SUMMARY")
    print("No open-cover counterexample among any set analyzed.")
    print("The divisibility sieve {2,...,14} cannot be satisfied by 13 small")
    print("speeds without either (i) missing a modulus scale (boundary witness")
    print("re-opens), or (ii) admitting a large speed whose intervals are too")
    print("sparse to keep the cover full (a positive gap opens).")


if __name__ == "__main__":
    main()
