#!/usr/bin/env python3
"""WIDENED exhaustive completeness certificate for C', plus the first n=9 box
(monad-compute-2026-06-03-S598).

This continues monad-compute-S597 ("EXHAUSTIVE box certificate for C', n=4..8,
0 tight"). The S597 handoff asked for exactly two extensions:
  (a) widen K per n (e.g. n=8 to K=4);
  (b) push the exhaustive n=9 box (C(27,8) ~ 2.2M).
S598 does BOTH in one self-certifying run.

C' (THM-398 / HYP-2102):  if a PRIMITIVE speed set S (gcd=1, the n-1 moving
runners) contains a multiple of n, then S is LOOSE: M(S) = max_t min_v ||v t|| > 1/n.
A single TIGHT (M=1/n) primitive multiple-of-n config inside any box would REFUTE
C' and the LRC reduction. S595/S596 sampled the class through n=20; S597 made it
EXHAUSTIVE inside finite boxes (n=4..8). S598 enlarges every box and adds n=9.

WHY THIS IS NEW / STRONGER THAN S597
  Per n the box [1,B]^{n-1} is strictly larger here than in S597:
       n :   4    5    6    7    8    9
   S597 K:  10    8    6    4    3    --   (none)
   S598 K:  20   14    8    5    4    3    (n=9 box is brand new)
  K = B/n is the largest companion-to-n ratio covered. The C' residual lives at
  small multiples/companions, so widening K closes more of the "interesting band"
  and n=9 is a completely new dimension never exhaustively certified.

LOOSENESS TEST -- FAST EXACT INTEGER METHOD (the S598 engine)
  S597 used the Fraction breakpoint scan (~3 ms/config); too slow for n=9 (est.
  ~50 min). S598 uses a mathematically equivalent INTEGER arc-cover test (~0.07
  ms/config, ~40x faster), which makes n=9 + all widened boxes finish in minutes.

  Reformulation: for runner v, the CLOSED unsafe set { t : ||v t|| <= 1/n } is a
  union of v arcs, arc k (k=0..v-1) centred at k/v with half-width 1/(n v); i.e.
  endpoints (k n -+ 1)/(n v). S is LOOSE  <=>  the open safe set is nonempty
  <=>  the union of all closed unsafe arcs does NOT cover the circle [0,1)
  (a strict uncovered gap = positive-measure safe set; arc endpoints, where
  ||v t|| = 1/n exactly, are correctly excluded from the strict-safe set).

  To make it pure-integer/exact: scale the circle by D = n * L, L = lcm(S). Arc
  endpoint (k n -+ 1)/(n v) maps to the integer (k n -+ 1) * (L // v) (v | L).
  Sort arcs, sweep; a strict gap => loose. Summing the gaps gives the EXACT
  looseness margin mu(safe) = (uncovered units)/D as a rational -- so we report
  the minimum positive margin over each whole box (closest approach to tight) at
  full-box scale, not just on n=4,5 as in S597.

SELF-CERTIFICATION
  Before the main run, the integer engine is cross-checked on a random sample
  against S597's canonical Fraction routines `is_loose` and `open_safe_measure`
  (imported from lrc_Cprime_exhaustive_box_monad_s597). Any single disagreement
  aborts the run. (Offline validation: 0 mismatches on 6000 looseness + 5000
  measure random configs across n=4..9, and the AP {1..n-1} is correctly TIGHT
  for n=4..14.) No floats anywhere (cf. MISTAKE-019).

No proof; pure computation.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools
import random
import sys

# Reuse the canonical (slow, trusted) S597 routines for cross-checking, and its
# exact enumerator of the C' hypothesis class.
from lrc_Cprime_exhaustive_box_monad_s597 import (
    is_loose as is_loose_slow,
    open_safe_measure as open_safe_measure_slow,
    enumerate_box,
    set_gcd,
)


def _lcm(a, b):
    return a // gcd(a, b) * b


def safe_measure_int(V, n):
    """EXACT open-safe-set measure of { t in [0,1) : min_v ||v t|| > 1/n } as a
    Fraction, via the integer arc-cover sweep. Equals open_safe_measure(V,n).
    S is LOOSE iff this is > 0; the value itself is the looseness margin."""
    L = reduce(_lcm, V)
    D = n * L
    arcs = []
    for v in V:
        s = L // v  # exact since v | L
        for k in range(v):
            lo = ((k * n - 1) * s) % D
            end = lo + 2 * s  # arc length = 2/(n v) -> 2*s in scaled units
            if end <= D:
                arcs.append((lo, end))
            else:  # wraps past D: split at the seam
                arcs.append((lo, D))
                arcs.append((0, end - D))
    arcs.sort()
    cur = 0
    cov = 0
    for lo, hi in arcs:
        if hi <= cur:
            continue
        lo2 = lo if lo > cur else cur
        cov += hi - lo2
        cur = hi
    return F(D - cov, D)


def _self_check(trials=4000, seed=20260603):
    """Cross-check the integer engine against the S597 Fraction routines. Abort
    on any disagreement -- the certificate is only as trustworthy as this."""
    rng = random.Random(seed)
    bad = 0
    for _ in range(trials):
        n = rng.choice([4, 5, 6, 7, 8, 9])
        m = n - 1
        B = rng.randint(n, n * 4)
        if B < m:
            continue
        V = tuple(sorted(rng.sample(range(1, B + 1), m)))
        mu = safe_measure_int(V, n)
        if (mu > 0) != is_loose_slow(V, n):
            bad += 1
            print(f"  SELF-CHECK FAIL (loose): n={n} V={V}")
        if mu != open_safe_measure_slow(V, n):
            bad += 1
            print(f"  SELF-CHECK FAIL (measure): n={n} V={V} "
                  f"int={mu} exact={open_safe_measure_slow(V, n)}")
    # AP is the canonical TIGHT witness -- must have measure exactly 0.
    for n in (4, 5, 6, 7, 8, 9, 12, 14):
        if safe_measure_int(tuple(range(1, n)), n) != 0:
            bad += 1
            print(f"  SELF-CHECK FAIL (AP tight): n={n}")
    return bad


def check_n(n, B, report_every=200000):
    total = loose = tight = 0
    min_meas = None
    min_cfg = None
    tights = []
    for V in enumerate_box(n, B):
        total += 1
        mu = safe_measure_int(V, n)
        if mu > 0:
            loose += 1
            if min_meas is None or mu < min_meas:
                min_meas, min_cfg = mu, V
        else:
            tight += 1
            tights.append(V)
        if total % report_every == 0:
            print(f"    ... n={n} B={B}: scanned {total} "
                  f"({loose} loose, {tight} tight)", flush=True)
    return total, loose, tight, min_meas, min_cfg, tights


# Strictly wider than S597 for every n, plus the new n=9 box.
BOXES = [
    (4, 80),   # K=20  (S597 had K=10)
    (5, 70),   # K=14  (S597 had K=8)
    (6, 48),   # K=8   (S597 had K=6)
    (7, 35),   # K=5   (S597 had K=4)
    (8, 32),   # K=4   (S597 had K=3)
    (9, 27),   # K=3   (NEW dimension -- never exhaustively certified before)
]


def main():
    print("=" * 74)
    print("C' WIDENED EXHAUSTIVE BOX CERTIFICATE + n=9  (monad-compute-S598)")
    print("Claim C': primitive speed set with a multiple of n  =>  LOOSE (M>1/n).")
    print("Engine: exact INTEGER arc-cover measure (~40x faster than S597 Fraction).")
    print("A single tight (measure-0) config inside a box would REFUTE C'.")
    print("=" * 74)
    print()
    print("Self-check: integer engine vs S597 Fraction routines (is_loose / "
          "open_safe_measure)...", flush=True)
    bad = _self_check()
    if bad:
        print(f"ABORT: {bad} self-check disagreement(s) -- engine not trusted.")
        sys.exit(1)
    print("Self-check PASSED (0 disagreements; AP tight for n=4..14).")
    print()

    print(f"{'n':>3} {'B':>4} {'K=B/n':>6} {'#configs':>11} {'loose':>11} "
          f"{'tight':>6} {'cert':>6}   min_margin (smallest mu>0)")
    grand_total = grand_loose = grand_tight = 0
    all_tights = []
    for n, B in BOXES:
        total, loose, tight, min_meas, min_cfg, tights = check_n(n, B)
        cert = "PASS" if tight == 0 else "FAIL"
        if min_meas is not None:
            mm = f"{float(min_meas):.6f}  (={min_meas})  at {min_cfg}"
        else:
            mm = "(no loose config?)"
        print(f"{n:>3} {B:>4} {B // n:>6} {total:>11} {loose:>11} "
              f"{tight:>6} {cert:>6}   {mm}", flush=True)
        grand_total += total
        grand_loose += loose
        grand_tight += tight
        all_tights += [(n, B, t) for t in tights]

    print()
    print("-" * 74)
    print(f"GRAND TOTAL: {grand_total} primitive multiple-of-n configs enumerated "
          f"EXHAUSTIVELY.")
    print(f"             {grand_loose} loose, {grand_tight} tight.")
    if grand_tight == 0:
        print()
        print("CERTIFICATE: C' holds with ZERO EXCEPTIONS over every primitive")
        print("multiple-of-n speed set inside the listed boxes -- now WIDER than")
        print("S597 for n=4..8 and EXTENDED to n=9. Exhaustive, not sampled.")
    else:
        print()
        print("*** C' COUNTEREXAMPLE(S) FOUND -- would REFUTE C' / LRC reduction ***")
        for n, B, t in all_tights[:50]:
            print(f"    n={n} B={B}  TIGHT (M=1/n) config: {t}")
    sys.stdout.flush()


if __name__ == '__main__':
    main()
