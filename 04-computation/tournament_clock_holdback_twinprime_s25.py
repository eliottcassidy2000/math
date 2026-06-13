#!/usr/bin/env python3
"""
tournament_clock_holdback_twinprime_s25.py

oracle-2026-06-01-S25

HOLDBACK, twin primes, and the LRC — through the tournament clock (S24).

Setup recap.  Speeds s_0<...<s_{n-1}.  Phase-comparator edge (i,j) flips at the
walls t = m / (2 d_ij), d_ij = s_j - s_i.  So edge (i,j) HOLDS one orientation
for exactly

    holdback(i,j) = 1 / (2 d_ij)

of the lap -- its persistence.  The stickiest (most "held back") edges are the
SMALL-DIFFERENCE pairs.  This script connects:

  * HOLDBACK  = 1/(2*difference)          -> persistence of a clock edge
  * min difference = MAX holdback         -> the "stickiest" relation
  * twin-prime pairs (difference 2)       -> the rate-2 held-back edges; the
       closest two primes can sit, the prime analogue of consecutive integers
  * SYNCHRONY: equal differences share identical walls (lockstep flips), so the
       clock complexity ~ #DISTINCT differences (Sidon = max, AP = min)
  * LRC / ADMISSIBILITY DUALITY:  by the sieve lemma (THM-369) a speed set has a
       lonely time t=1/m whenever it MISSES residue 0 mod m (is "admissible mod
       m").  Prime/twin-prime sets miss almost all small residues -> LRC-easy;
       the consecutive set {1..n-1} covers every residue -> anti-admissible ->
       LRC-extremal (and, S24, the MINIMAL clock with MAX holdback).
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parent
S356 = SourceFileLoader("s356", str(ROOT / "lonely_runner_residue_probe_s356.py")).load_module()


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def differences(speeds):
    return Counter(b - a for a, b in combinations(sorted(speeds), 2))


def holdback_spectrum(speeds):
    diffs = differences(speeds)
    # holdback = 1/(2d); group by value
    return {Fraction(1, 2 * d): cnt for d, cnt in diffs.items()}


def admissibility_profile(speeds, mmax):
    """For each m in 2..mmax: does the speed set cover residue 0 mod m?
    'covered' = some speed divisible by m (anti-admissible at m, LRC-hard at m).
    'missed'  = no speed divisible by m  => t=1/m is lonely (admissible)."""
    covered, missed = [], []
    for m in range(2, mmax + 1):
        if any(s % m == 0 for s in speeds if s != 0):
            covered.append(m)
        else:
            missed.append(m)
    return covered, missed


def lrc_gap(speeds):
    nz = tuple(s for s in speeds if s != 0)
    if len(nz) < 2:
        return None
    r = S356.report("x", list(nz))
    return r.max_gap, r.threshold, r.forbidden_length


def fmt(x):
    if x is None:
        return "-"
    if isinstance(x, Fraction):
        return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)
    return str(x)


def report(label, speeds):
    speeds = tuple(sorted(set(speeds)))
    diffs = differences(speeds)
    distinct = len(diffs)
    total_pairs = sum(diffs.values())
    min_d = min(diffs)
    twin_edges = diffs.get(2, 0)                 # difference-2 pairs (twin-prime-like)
    consec_edges = diffs.get(1, 0)               # difference-1 pairs (consecutive)
    n = len(speeds)
    covered, missed = admissibility_profile(speeds, n)  # moduli up to n
    g = lrc_gap(speeds)
    print(f"[{label}] speeds={speeds}")
    print(f"   pairs={total_pairs} distinct_diffs={distinct}  "
          f"(Sidon-ratio={distinct}/{total_pairs})  min_diff={min_d} "
          f"=> MAX holdback={fmt(Fraction(1,2*min_d))}")
    print(f"   consecutive(d=1) edges={consec_edges}   twin(d=2) edges={twin_edges}")
    print(f"   LRC mod-cover: covered(anti-adm)={covered}  missed(adm->lonely t=1/m)={missed}")
    if g:
        gap, thr, forb = g
        status = ("LONELY easy (positive gap)" if gap > 0
                  else "tight/boundary")
        print(f"   LRC gap={fmt(gap)} thr={fmt(thr)} forbidden_len={fmt(forb)}  -> {status}")
    print()


def main():
    print("Holdback / twin primes / LRC via the tournament clock (oracle-S25)\n")

    print("=" * 70)
    print("A. The extremal axis: consecutive integers vs primes vs twin primes")
    print("=" * 70)
    report("consecutive 1..7 (LRC-extremal)", range(1, 8))
    report("primes <=17 (7 of them)", [2, 3, 5, 7, 11, 13, 17])
    report("odd primes 3..19 (7)", [3, 5, 7, 11, 13, 17, 19])
    # twin-prime members only: 3,5,11,13,17,19,29,31...
    report("twin-prime members (7)", [3, 5, 11, 13, 17, 19, 29])
    report("a pure twin cluster 3,5 / 11,13 / 17,19 / 29 (=7)", [3, 5, 11, 13, 17, 19, 29])

    print("=" * 70)
    print("B. Synchrony: equal differences flip in lockstep (clock complexity)")
    print("=" * 70)
    print(" #distinct differences controls how many independent flip-rhythms the")
    print(" clock has. Arithmetic progressions (all diffs multiples of a gap) and")
    print(" twin clusters (many d=2) synchronise; Sidon sets maximise complexity.\n")
    report("AP gap 1: 1..7", range(1, 8))
    report("AP gap 3: 3,6,9,12,15,18,21", [3, 6, 9, 12, 15, 18, 21])
    report("Sidon-ish 1,2,5,11,22,33,50", [1, 2, 5, 11, 22, 33, 50])

    print("=" * 70)
    print("C. Holdback spectrum (persistence = 1/(2d)) for a few sets")
    print("=" * 70)
    for label, s in {
        "consecutive 1..6": list(range(1, 7)),
        "twin members 3,5,11,13,17,19": [3, 5, 11, 13, 17, 19],
        "geometric 1,2,4,8,16,32": [1, 2, 4, 8, 16, 32],
    }.items():
        spec = holdback_spectrum(s)
        top = sorted(spec.items(), key=lambda kv: -kv[0])[:5]
        print(f" [{label}] holdbacks (value:count, stickiest first): "
              + ", ".join(f"{fmt(v)}:{c}" for v, c in top))
    print()

    print("SUMMARY")
    print(" holdback = 1/(2*difference): small speed-gap = sticky edge.")
    print(" Consecutive integers: min_diff=1 => MAX holdback, covers all moduli")
    print("   => anti-admissible => LRC-extremal (and minimal clock, S24).")
    print(" Twin primes (d=2) are the stickiest pairs the primes allow -- the")
    print("   prime analogue of consecutive integers -- but prime sets MISS the")
    print("   small residues (admissible) and so are LRC-easy. The tension:")
    print("   holdback wants LRC-hardness; admissibility (primes) grants ease.")


if __name__ == "__main__":
    main()
