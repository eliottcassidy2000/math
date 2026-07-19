#!/usr/bin/env python3
"""
boxeph-2026-07-19-S122 — the COVERING-restricted rigidity margin, through the determinant lens.

CONTEXT. S121 (witness-blocking cascade): M(C)=1/13 forces C to contain a multiple of every
q in {2,...,12} (covering) -- this blocks EVERY witness t=p/q with q<=12 at once (a multiple of
q kills t=p/q for all p). So covering is the necessary "block the small sieves" condition, and
the residual tightness is decided by q>=13 witnesses, which by the Pinch Lemma sit at pairwise
SUMS q=v_i+v_j >= 13. opus THM-1210: M = D/(v_i+v_j) with D = |v_i a_j - v_j a_i| the DETERMINANT;
D=1 is exactly the classical sieve at modulus s=v_i+v_j.

So the sharpest probe of Tao n=12 uniqueness is NOT the raw spectrum (dominated by non-covering
sets) but the COVERING-restricted spectrum: among 12-sets that already block all small sieves,
how much lonelier than 1/13 is the tightest one that is not {1,...,12}?

This script:
 (A) enumerates primitive COVERING 12-subsets of {1,...,N} (contain a multiple of each q=2..12),
     computes M (Pinch maximizer over pairwise sums, all numerators), and reports the tightest ones
     != {1,...,12}, with the maximizing pair, its sum s, and determinant D = M*s;
 (B) the RESIDUE-LIFT family: {1,...,12} with a subset of residues r lifted to r+13 (still a
     complete residue system mod 13, still covering-ish) -- to see the mechanism by which lifting
     an element past 13 opens a D>=2 witness at a larger pair-sum and pushes M above 1/13.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations


def fd(x: F) -> F:
    r = x - (x.numerator // x.denominator)
    return min(r, 1 - r)


def fmin(sp, t: F) -> F:
    return min(fd(F(v) * t) for v in sp)


def M_pair_det(sp):
    """M via Pinch maximizer; return (M, s=v_i+v_j, D=det, pair). Search all m/q, q a pair-sum."""
    sums = {}
    n = len(sp)
    for i in range(n):
        for j in range(i + 1, n):
            sums.setdefault(sp[i] + sp[j], (sp[i], sp[j]))
    best, bq, bm = F(0), None, None
    for q in sums:
        for m in range(1, q):
            v = fmin(sp, F(m, q))
            if v > best:
                best, bq, bm = v, q, m
    D = best.numerator * bq // best.denominator if best.numerator else 0
    # D = best * bq  as integer (best = D/bq)
    D = (best * bq)
    return best, bq, D, sums[bq]


def primitive(sp):
    g = 0
    for v in sp:
        g = gcd(g, v)
    return g == 1


def covering(sp):
    s = set(sp)
    return all(any(v % q == 0 for v in sp) for q in range(2, 13))


THIRT = F(1, 13)
BASE = tuple(range(1, 13))


def run_cover(N):
    print("=" * 78, flush=True)
    print(f"(A) tightest primitive COVERING 12-subsets of {{1,...,{N}}} (!= {{1..12}})", flush=True)
    print("=" * 78, flush=True)
    rows = []
    ncov = 0
    for combo in combinations(range(1, N + 1), 12):
        if not primitive(combo) or not covering(combo):
            continue
        ncov += 1
        M, s, D, pair = M_pair_det(list(combo))
        rows.append((M, combo, s, D, pair))
    rows.sort(key=lambda r: (r[0], r[1]))
    print(f"  primitive covering 12-subsets: {ncov}", flush=True)
    print(f"  {'M':>8} {'~dec':>7} {'pair':>10} {'s':>4} {'D':>3}  set", flush=True)
    shown = 0
    for M, combo, s, D, pair in rows:
        tag = "  <-- {1..12}" if combo == BASE else ""
        exs = "{" + ",".join(map(str, combo)) + "}"
        print(f"  {str(M):>8} {float(M):>7.4f} {str(pair):>10} {s:>4} {str(D):>3}  {exs}{tag}", flush=True)
        shown += 1
        if shown >= 14:
            break
    # the covering-rigidity gap: tightest covering set that is not {1..12}
    non_base = [r for r in rows if r[1] != BASE]
    if non_base and rows[0][1] == BASE:
        m0, m1 = rows[0][0], non_base[0][0]
        print(f"\n  tightest covering = {{1..12}} at {m0}; tightest covering competitor = {m1} "
              f"({non_base[0][1]}); covering-rigidity gap = {m1 - m0} ~ {float(m1-m0):.5f}", flush=True)
    return rows


rows18 = run_cover(18)

print("", flush=True)
print("=" * 78, flush=True)
print("(B) RESIDUE-LIFT family: {1..12} with residues in L lifted r -> r+13", flush=True)
print("=" * 78, flush=True)
print("    (each stays a complete residue system mod 13; measures how lifting opens a witness)", flush=True)
print(f"  {'lifted set L':>18} {'M':>8} {'~dec':>7} {'pair':>10} {'s':>4} {'D':>3}  set(short)", flush=True)
import itertools
# lift small subsets of residues
base = list(range(1, 13))
tests = [[], [12], [11], [1], [6], [12, 11], [1, 12], [11, 12, 10], [1, 2, 3]]
for L in tests:
    C = sorted((r + 13) if r in L else r for r in base)
    M, s, D, pair = M_pair_det(C)
    exs = "{" + ",".join(map(str, C[:2])) + ",..," + str(C[-1]) + "}"
    print(f"  {str(L):>18} {str(M):>8} {float(M):>7.4f} {str(pair):>10} {s:>4} {str(D):>3}  {exs}",
          flush=True)

print("\n  reading: L=[] is {1..12} (M=1/13, D=1, s=13); any lift raises M and (typically) D,", flush=True)
print("  and moves the maximizer to a larger pair-sum -- the D>=2 mechanism opus THM-1210 predicts.", flush=True)
print("\nDONE.", flush=True)
