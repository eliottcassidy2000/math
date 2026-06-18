#!/usr/bin/env python3
"""
LRC(14) S3 residual -- THE MULTI-BAND-FIT LEMMA, made precise and tested.

Structure of S3:  S = P  ∪  L,   P ⊆ {1..13} (small part), L = cluster (all >13),
                  with max(S) = Vmax ∈ L, internal spread s = max(L)-min(L).

THE LEMMA WE WANT (multi-band-fit / generalized LEMMA 2):
  The small part P has level-1/14 safe set with a WIDE arc I_P (width W_P >= mu_P/N_P,
  pigeonhole, with N_P <= sum P <= 91). Inside I_P, the cluster L plants teeth: runner
  u ∈ L has teeth spaced 1/u ≈ 1/V apart, each of width 1/(7u) ≈ 1/(7V). Across the
  whole cluster, the union of teeth inside I_P has measure <= |L| * (#teeth per runner
  in I_P) * 1/(7V). If I_P is WIDE (width ~ const) it contains ~ V*W_P teeth-slots per
  runner, so a naive union bound is useless. BUT the cluster teeth are NOT independent:
  because L is tight (spread s), the c=|L| runners' teeth nearly COINCIDE -- they form c
  near-parallel comb, and the LARGEST GAP between consecutive merged cluster-teeth inside
  I_P is ~ 1/(c'V) for some effective count, which we must show exceeds 1/(7Vmax).

  So the precise quantity to bound below is:
     G* := the largest gap, inside the small-part safe arc I_P, that is free of ALL
           cluster teeth (i.e. a point safe for P AND for all of L).
     CLAIM: G* > 1/(7*Vmax).   [then its midpoint is a global witness, M(S)>=1/14]

THIS SCRIPT:
  1. Computes G* exactly for many S3 sets and checks G* > 1/(7*Vmax) (the *global
     witness* form -- stronger than C(S) via v=max, since it needs NO runner removed!).
     Wait: a global witness means tau safe for ALL of S incl Vmax. That proves M>=1/14
     OUTRIGHT. Let's test whether such a global witness exists in I_P.
  2. If global-witness-in-I_P fails sometimes, fall back to the leave-one-out (remove
     Vmax) version: largest cluster-tooth-free gap inside I_P for the cluster MINUS Vmax,
     vs 1/(7*Vmax). That is exactly W(S\\{Vmax}) restricted to I_P -- and equals C via max.
  3. Decompose G* lower bound: relate it to s (spread) and c (cluster size) to get a
     closed-form provable bound.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random
from collections import Counter

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def danger_intervals(A, lo, hi, h=F(1, 14)):
    """All danger teeth of runners in A that intersect [lo,hi] (lo,hi in [0,1), lo<hi).
    Returns list of (a,b) clipped to [lo,hi]."""
    out = []
    for u in A:
        # teeth centered j/u, half width h/u. Find j with center near [lo,hi].
        jmin = int(lo * u) - 1; jmax = int(hi * u) + 1
        for j in range(jmin, jmax + 1):
            c = F(j, u); a = c - h/u; b = c + h/u
            # clip to [lo,hi]
            aa = max(a, lo); bb = min(b, hi)
            if aa < bb: out.append((aa, bb))
    return out

def safe_components(A, h=F(1, 14)):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h/u) % 1; b = (c + h/u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
    safe = []; prev = F(0)
    for a, b in merged:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe

def safe_arcs(A): return safe_components(A)

def largest_free_gap_in(interval, danger):
    """interval=(lo,hi). danger=list of (a,b) sub-intervals. Return width of largest
    sub-interval of [lo,hi] free of all danger (closed teeth: free = open complement)."""
    lo, hi = interval
    ivs = sorted((max(a, lo), min(b, hi)) for a, b in danger if min(b, hi) > max(a, lo))
    merged = []
    for a, b in ivs:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
    best = F(0); prev = lo; argbest = None
    for a, b in merged:
        if a - prev > best: best = a - prev; argbest = (prev, a)
        prev = max(prev, b)
    if hi - prev > best: best = hi - prev; argbest = (prev, hi)
    return best, argbest

def Wwidth(A):
    sc = safe_components(A)
    if not sc: return F(0)
    ws = [b - a for a, b in sc]
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1:
        ws.append((sc[0][1]) + (1 - sc[-1][0]))
    return max(ws)

def gcd_all(S): return reduce(gcd, S, 0)
def is_covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def classify(S):
    S = sorted(set(S)); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13 * Vmin: return 'S2'
    return 'S3'

def gen_S3(seed=0, target=2000, vmax_cap=1500):
    rng = random.Random(seed); out = []; tries = 0
    while len(out) < target and tries < target * 400:
        tries += 1
        c = rng.choice([2, 2, 3, 3, 4, 5, 6]); nsmall = 13 - c
        if nsmall < 1: continue
        V = rng.choice([rng.randint(40, 80), rng.randint(80, 200),
                        rng.randint(200, 600), rng.randint(600, vmax_cap)])
        spread = rng.choice([14, 20, 28, 35, 42, 56, 70])
        cluster = set()
        while len(cluster) < c: cluster.add(V + rng.randint(0, spread))
        small = set(); pool = list(range(1, 14)); rng.shuffle(pool)
        for x in pool:
            if len(small) >= nsmall: break
            small.add(x)
        if len(small) < nsmall: continue
        S = sorted(small | cluster)
        if len(S) != 13 or gcd_all(S) != 1: continue
        if not is_covering(S) or classify(S) != 'S3': continue
        out.append(S)
    return out

if __name__ == "__main__":
    S3 = gen_S3(seed=5, target=600)
    print(f"n S3 = {len(S3)}")
    glob_fires = 0; loo_fires = 0
    glob_margins = []; loo_margins = []
    glob_fail = []
    for S in S3:
        Vmax = S[-1]
        small = [u for u in S if u <= 13]
        cluster = [u for u in S if u > 13]
        thr = F(1, 7*Vmax)
        # small-part safe arcs
        sc = safe_arcs(small)
        # widest small arc (no wrap-handling needed for global witness existence; just scan all)
        # (1) GLOBAL witness: largest gap free of ALL cluster teeth, inside any small arc
        bestG = F(0)
        for (lo, hi) in sc:
            dgr = danger_intervals(cluster, lo, hi)
            g, _ = largest_free_gap_in((lo, hi), dgr)
            if g > bestG: bestG = g
        glob_margins.append(float(bestG * 7 * Vmax))
        if bestG > thr: glob_fires += 1
        else: glob_fail.append((S, float(bestG*7*Vmax)))
        # (2) leave-one-out (remove Vmax): cluster' = cluster\{Vmax}
        cl2 = [u for u in cluster if u != Vmax]
        bestL = F(0)
        for (lo, hi) in sc:
            dgr = danger_intervals(cl2, lo, hi)
            g, _ = largest_free_gap_in((lo, hi), dgr)
            if g > bestL: bestL = g
        loo_margins.append(float(bestL * 7 * Vmax))
        if bestL > thr: loo_fires += 1
    n = len(S3)
    def st(a): return f"min={min(a):.3f} mean={sum(a)/len(a):.3f} max={max(a):.3f}"
    print(f"\n(1) GLOBAL witness inside small arc, free of ALL cluster teeth:")
    print(f"    G*·7·Vmax: {st(glob_margins)}  -> M(S)>=1/14 OUTRIGHT in {glob_fires}/{n}")
    print(f"(2) Leave-one-out (drop Vmax), free of cluster\\Vmax teeth:")
    print(f"    ·7·Vmax: {st(loo_margins)}  -> fires in {loo_fires}/{n}")
    if glob_fail:
        print(f"\n  GLOBAL-witness FAILS on {len(glob_fail)} sets (need leave-one-out):")
        for S, m in glob_fail[:10]:
            print(f"     S={S}  G*·7·max={m:.3f}")
