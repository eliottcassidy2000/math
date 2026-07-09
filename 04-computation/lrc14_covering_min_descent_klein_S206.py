#!/usr/bin/env python3
"""
klein-2026-07-09-S206: HOW LOW DOES THE COVERING-MIN GO at n=14?  (the prove/disprove frontier)

THM-523: a 13-set omitting a multiple of some q in {2..14} has the explicit witness tau=1/q, so
M>=1/q>=1/14.  Hence LRC(14) reduces to COVERING sets (a multiple of every q in {2..14}).  The tight
(M=1/14) sets are all NON-covering (verified S206 for n=4..7 and here for {1..13}).  So:

    LRC(14) is TRUE  <=>  inf{ M(S) : S primitive covering 13-set }  >=  1/14.
    LRC(14) is FALSE <=>  some covering 13-set has M(S) < 1/14.

The repo records min M ~ 7/89 = 0.078652 (THM-523, "10% above 1/14"), with UNIFORM LOOSENESS only
CONJECTURAL (HYP-2566).  This script measures the DESCENT of the covering-min with the speed cap, and
hill-climbs at larger speeds, to see whether the infimum plateaus above 1/14 or descends toward it.

EXACT M(S) = max_t min_i ||v_i t||; local maxima of the piecewise-linear min sit at t=p/(v_i+v_j).
"""
import sys, random
from math import gcd
from fractions import Fraction
from itertools import combinations

def M_exact(v):
    qs = sorted({v[i] + v[j] for i in range(len(v)) for j in range(i, len(v))})
    bm, bq = 0, 1
    for q in qs:
        for p in range(1, q):
            m = q
            for vk in v:
                r = (vk * p) % q
                d = r if r <= q - r else q - r
                if d < m:
                    m = d
                    if m * bq <= bm * q: break
            if m * bq > bm * q: bm, bq = m, q
    return Fraction(bm, bq)

def is_cov(S): return all(any(s % q == 0 for s in S) for q in range(2, 15))
def prim(S):
    g = 0
    for s in S: g = gcd(g, s)
    return g == 1

TGT = Fraction(1, 14); REPO = Fraction(7, 89)
print(f"1/14 = {float(TGT):.6f}    repo covering-min 7/89 = {float(REPO):.6f}\n", flush=True)

print("(1) EXHAUSTIVE covering-min vs speed cap  (does it descend?)", flush=True)
print(f"{'cap':>4} {'#cov prim':>10} {'min M':>10} {'float':>10} {'/(1/14)':>9} {'argmin':>0}", flush=True)
for cap in [15, 16, 17, 18, 19, 20]:
    best = None; cnt = 0
    for S in combinations(range(1, cap + 1), 13):
        if not prim(S) or not is_cov(S): continue
        cnt += 1
        M = M_exact(list(S))
        if best is None or M < best[0]: best = (M, S)
    if best:
        M, S = best
        print(f"{cap:>4} {cnt:>10} {str(M):>10} {float(M):>10.6f} {float(M/TGT):>9.4f}  {S}", flush=True)

print("\n(2) HILL-CLIMB at larger speeds (minimise M over primitive covering 13-sets, speeds<=CAP)", flush=True)
random.seed(206)
def rand_cov(cap):
    while True:
        S = sorted(random.sample(range(1, cap + 1), 13))
        if prim(S) and is_cov(S): return S
for CAP in [40, 90, 200]:
    best = None
    for restart in range(6):
        S = rand_cov(CAP)
        M = M_exact(S)
        improved = True
        it = 0
        while improved and it < 400:
            improved = False; it += 1
            for idx in range(13):
                for _ in range(14):
                    T = list(S); T[idx] = random.randint(1, CAP)
                    T = sorted(set(T))
                    if len(T) != 13 or not prim(T) or not is_cov(T): continue
                    MT = M_exact(T)
                    if MT < M:
                        S, M = T, MT; improved = True; break
        if best is None or M < best[0]: best = (M, list(S))
    M, S = best
    flag = ""
    if M < REPO: flag = "  <-- BEATS repo 7/89"
    if M < TGT: flag = "  <-- BELOW 1/14: COUNTEREXAMPLE"
    print(f"  CAP={CAP:>3}: best M = {str(M):>12} = {float(M):.6f}  ratio/(1/14) = {float(M/TGT):.4f}{flag}", flush=True)
    print(f"           S = {S}", flush=True)

print("\nREADING: if the covering-min DESCENDS toward 1/14 as speeds grow, the covering branch is also", flush=True)
print("'barely true' and uniform looseness (HYP-2566) is the crux. If it PLATEAUS above 1/14, a margin", flush=True)
print("argument suffices there and LRC(14) = [tau=1/q witness on non-covering] + [margin on covering].", flush=True)
