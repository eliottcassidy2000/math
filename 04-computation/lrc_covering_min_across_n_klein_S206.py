#!/usr/bin/env python3
"""
klein-2026-07-09-S206: the covering/non-covering dichotomy across runner counts n.

THM-523 (n=14): a 13-set OMITTING a multiple of some q in {2..14} has the explicit lonely witness
tau = 1/q, so M(S) >= 1/q >= 1/14.  Hence a counterexample must be a COVERING set (contains a multiple
of EVERY q in {2..n}).  The tight minimizers (M = 1/14 exactly) are all in the trivial/non-covering
case; the covering core is empirically LOOSE (min M = 7/89 ~ 10% above 1/14; the proved lower bound is
n/Phi6(n) = 14/183, gap 35/16287).

TEST the same dichotomy for general n:
  (A) is every TIGHT set (M(S) = 1/n exactly) NON-covering?          [the equality locus dodges covering]
  (B) is the covering-min strictly > 1/n, and how big is the cushion? [does a margin argument have room?]
  (C) does the cushion track n^2/Phi6(n) -> 1 (the Phi6 resonance)?   [why hardness grows with n]

EXACT M(S) = max_t min_i ||v_i t||: f(t) is piecewise linear with slopes +-v_i, so every local max sits
at t = p/(v_i+v_j) (including q = 2 v_i, the single-piece apex).  So scan q in {v_i+v_j}, p in 1..q-1.
"""
import sys
from math import gcd
from fractions import Fraction
from itertools import combinations

def phi6(n): return n*n - n + 1

def M_exact(v):
    """exact (M, t) with M = max_t min_i ||v_i t||, comparing m/q by cross-multiplication."""
    qs = sorted({v[i] + v[j] for i in range(len(v)) for j in range(i, len(v))})
    bm, bq, bp = 0, 1, 0
    for q in qs:
        for p in range(1, q):
            m = q
            for vk in v:
                r = (vk * p) % q
                d = r if r <= q - r else q - r
                if d < m:
                    m = d
                    if m * bq <= bm * q: break     # prune: cannot beat current best
            if m * bq > bm * q:
                bm, bq, bp = m, q, p
    return Fraction(bm, bq), Fraction(bp, bq)

def is_covering(S, n):
    return all(any(s % q == 0 for s in S) for q in range(2, n + 1))

print("Dichotomy: TIGHT sets vs COVERING sets, exact, per n.\n", flush=True)
print(f"{'n':>3} {'cap':>4} {'#prim':>8} {'#cov':>6} {'#tight':>7} {'#tight&cov':>11} "
      f"{'min M(cov)':>12} {'1/n':>8} {'cushion':>9} {'n^2/Phi6':>9}", flush=True)

for n, cap in [(4, 20), (5, 22), (6, 24), (7, 26)]:
    k = n - 1
    tgt = Fraction(1, n)
    lb = Fraction(n, phi6(n))
    nprim = ncov = ntight = ntightcov = 0
    best = None
    tight_examples = []
    for S in combinations(range(1, cap + 1), k):
        g = 0
        for s in S: g = gcd(g, s)
        if g != 1: continue                      # primitive only
        nprim += 1
        cov = is_covering(S, n)
        if cov: ncov += 1
        M, t = M_exact(list(S))
        if M == tgt:
            ntight += 1
            if len(tight_examples) < 4: tight_examples.append((S, str(t)))
            if cov: ntightcov += 1
        if cov and (best is None or M < best[0]):
            best = (M, t, S)
    mM, mt, mS = best if best else (None, None, None)
    cushion = float(mM / tgt) if mM else float('nan')
    print(f"{n:>3} {cap:>4} {nprim:>8} {ncov:>6} {ntight:>7} {ntightcov:>11} "
          f"{str(mM):>12} {str(tgt):>8} {cushion:>9.4f} {float(Fraction(n*n,phi6(n))):>9.4f}", flush=True)
    print(f"      covering-min at S={mS}, t={mt};  proved LB n/Phi6(n)={lb}={float(lb):.6f}, "
          f"min>=LB? {mM>=lb}", flush=True)
    print(f"      tight examples (all should be NON-covering): {tight_examples}", flush=True)

print("\nREADING:", flush=True)
print(" (A) #tight&cov == 0  =>  the EQUALITY locus M=1/n is entirely NON-covering: it is certified by the", flush=True)
print("     explicit witness tau=1/q (THM-523), never by a margin/averaging argument.", flush=True)
print(" (B) min M(cov) > 1/n  =>  the covering branch has a STRICT cushion: margin arguments have room.", flush=True)
print(" (C) cushion vs n^2/Phi6(n): the Phi6 resonance is the proved floor n/Phi6(n); the cushion shrinks", flush=True)
print("     to 1 as n grows -- that, not any n=14-specific coincidence, is why LRC gets harder.", flush=True)
