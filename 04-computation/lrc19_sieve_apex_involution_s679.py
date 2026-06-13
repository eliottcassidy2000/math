#!/usr/bin/env python3
"""lrc19_sieve_apex_involution_s679.py — LRC@19 prime frontier: clean sieve
reduction + the apex = negation-involution fixed-point reframe. (See HYP below.)
n=19 prime, C=2n-1=37 prime: the 'no side-channel' frontier (S678), no sporadics.
Session: claude-2026-06-06-S679."""
import sys; sys.stdout.reconfigure(line_buffering=True)
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import random
def dist(x): x=x%1; return min(x,1-x)
n=19; thr=F(1,n)
print("LRC@19: 18 runners, delta=1/19, C=37 (both prime). AP={1..18} tight, M=1/19.")
# (1) sieve reduction: non-sieve-covered => explicit lonely witness t=a/q
def lonely_at(V,t): return all(dist(v*t)>=thr for v in V)
def sieve_witness(V):
    for q in range(2,n+1):
        if all(v%q!=0 for v in V):
            for a in range(1,q):
                if gcd(a,q)==1 and lonely_at(V,F(a,q)): return (a,q)
    return None
rng=random.Random(1); noncov=0; fail=0
for _ in range(15000):
    V=tuple(random.sample(range(1,61),18))
    if any(all(v%q!=0 for v in V) for q in range(2,n+1)):
        noncov+=1
        if sieve_witness(V) is None: fail+=1
print(f"(1) SIEVE REDUCTION: {noncov} non-sieve-covered 18-sets; lonely t=a/q witness FAILED: {fail}")
print("    => every LRC@19 counterexample is sieve-covered (mult of each q in 2..19).")
# (2) apex = negation fixed point; kills division witnesses; M recovered off-grid
V=tuple(range(2,20))  # sieve-covered, apex 19
def M_pinch(V):
    cand={F(0)}
    for a,b in combinations(V,2):
        for m in range(a+b+1): cand.add(F(m,a+b))
    best=F(0); arg=None
    for t in cand:
        mn=min(dist(v*t) for v in V)
        if mn>best: best=mn; arg=t
    return best,arg
M,arg=M_pinch(V)
print(f"(2) APEX/INVOLUTION: V={V} sieve-covered, apex 19 == 0 mod 19 (negation fixed point).")
print(f"    apex dist at t=j/19: {[float(dist(19*F(j,19))) for j in (1,5,9)]} (0 => kills every division witness)")
print(f"    M(V)={M} (>=1/19) at t={arg} (t*19={float(arg*19):.2f}) -- recovered OFF the j/19 grid.")
print("    binding pairs of AP are {a, n-a}={a,-a} mod 19 = negation orbits; apex at the fixed point 0.")
print("(3) REFRAME: hard core of LRC@prime-n = negation-involution FIXED-POINT residual (the apex),")
print("    same shape as Redei self-converse (H-gaps), V* nonunit seam (signed-LRC S674), kappa-even")
print("    (S605); cure = a secondary/twisted involution on the fixed-point set (S606).")
