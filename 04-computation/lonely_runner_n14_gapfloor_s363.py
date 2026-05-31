#!/usr/bin/env python3
"""
lonely_runner_n14_gapfloor_s363.py

oracle-2026-05-31-S16

Stress-test the n=14 "gap-floor" hypothesis:

  HYP (n=14 gap floor).  Every PRIMITIVE 13-speed set V (gcd 1) that contains a
  multiple of 14 has max_gap > 0 -- i.e. it is not even tight, let alone a
  full-cover counterexample.

Why this matters.  A counterexample to LRC at n=14 must contain a multiple of
14 (else t=1/14 is a boundary lonely witness; cf. modulus-sieve lemma).  If the
gap floor held, no such set could be a cover, so LRC@14 would follow.  We do NOT
expect to prove it here; we search aggressively for a violation (a tight, gap=0
set containing a multiple of 14) and record the smallest gaps found.

Strategy:
  1. Structured perturbations of the initial segment that keep many small speeds
     (these are the only plausibly-tight neighborhoods) while forcing a
     multiple of 14 in.
  2. Random primitive 13-subsets of [1..B] that contain a multiple of 14,
     biased toward small speeds + full divisibility-sieve coverage.
Report the minimum max_gap encountered and any tight (gap=0) hits.
"""

from __future__ import annotations

import random
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()

N = 14
THR = Fraction(1, N)


def primitive(vs) -> bool:
    g = 0
    for v in vs:
        g = gcd(g, v)
    return g == 1


def has_mult_of_14(vs) -> bool:
    return any(v % 14 == 0 for v in vs)


def max_gap(vs) -> Fraction:
    return S356.report("x", list(vs)).max_gap


def gap_and_status(vs):
    r = S356.report("x", list(vs))
    return r.max_gap, r.forbidden_length, r.boundary_witness_count


def main() -> None:
    print("n=14 gap-floor stress test (oracle-2026-05-31-S16)")
    print("Searching primitive 13-sets containing a multiple of 14 for max_gap=0.\n")

    best = None          # (gap, speeds)
    tight_hits = []      # gap == 0
    cover_hits = []      # gap == 0 and 0 boundary witnesses (a TRUE counterexample)
    tested = 0

    def consider(vs):
        nonlocal best, tested
        vs = tuple(sorted(set(vs)))
        if len(vs) != 13 or not primitive(vs) or not has_mult_of_14(vs):
            return
        tested += 1
        g, flen, bw = gap_and_status(vs)
        if best is None or g < best[0]:
            best = (g, vs)
        if g == 0:
            tight_hits.append((vs, flen, bw))
            if bw == 0 and flen == 1:
                cover_hits.append(vs)

    # ---- 1. Structured perturbations of the initial segment -------------
    base = list(range(1, 14))
    mults14 = [14, 28, 42, 56, 70, 84, 98, 112, 126, 140, 154, 168, 182]
    # single drop, add one multiple of 14
    for drop in base:
        for m14 in mults14:
            consider(set(base) - {drop} | {m14})
    # double drop, add a multiple of 14 plus a nearby small replacement
    for d1, d2 in combinations(base, 2):
        for m14 in [14, 28, 42]:
            for repl in range(14, 41):
                if repl % 14 == 0:
                    continue
                consider(set(base) - {d1, d2} | {m14, repl})

    print(f"after structured perturbations: tested={tested}, "
          f"min_gap={float(best[0]):.6g} at {best[1]}", flush=True)

    # ---- 2. Randomized search, small-speed biased, sieve-complete -------
    rng = random.Random(363)
    B = 42
    pool = list(range(1, B + 1))
    weights = [1.0 / v for v in pool]  # bias toward small speeds
    NITER = 40000
    for it in range(NITER):
        # force at least one multiple of 14
        m14 = rng.choice([14, 28, 42])
        rest = set()
        # weighted sample without replacement
        while len(rest) < 12:
            v = rng.choices(pool, weights=weights, k=1)[0]
            if v != m14:
                rest.add(v)
        consider(rest | {m14})
        if (it + 1) % 10000 == 0:
            print(f"   ...{it+1}/{NITER} random sets, tested={tested}, "
                  f"min_gap={float(best[0]):.6g}", flush=True)

    print(f"after randomized search: tested={tested}, "
          f"min_gap={float(best[0]):.8g} at {best[1]}")

    print()
    print(f"TOTAL primitive 13-sets w/ multiple of 14 tested: {tested}")
    print(f"minimum max_gap found: {best[0]}  (= {float(best[0]):.8g})  at {best[1]}")
    print(f"tight (gap=0) hits: {len(tight_hits)}")
    for vs, flen, bw in tight_hits[:20]:
        print(f"   TIGHT {vs}  forbidden={flen} boundary_witnesses={bw}")
    print(f"full-cover COUNTEREXAMPLES (gap=0, no boundary witness): {len(cover_hits)}")
    for vs in cover_hits[:20]:
        print(f"   ***COUNTEREXAMPLE*** {vs}")
    if not tight_hits:
        print("\nNo tight set with a multiple of 14 was found: the gap floor held")
        print("across the entire search (evidence FOR LRC@14, not a proof).")


if __name__ == "__main__":
    main()
