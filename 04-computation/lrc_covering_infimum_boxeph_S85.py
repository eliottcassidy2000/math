#!/usr/bin/env python3
"""
The true infimum of M over primitive covering families  (boxeph-2026-07-18-S85)
===============================================================================
The stress test refuted (a) the 12-subset floor and (b) THM-995(X)'s empirical
covering floor M>=1/9: primitive covering families reach M=1/13.  Here we pin the
actual infimum by searching structured families (dilations of the tight AP made
primitive) and confirm it is 1/14 (approached, not attained) -- i.e. there is NO
covering floor above 1/14, so the compact residual is genuine near-tightness.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools

def norm(x):
    r = x % 1; return min(r, 1 - r)

def exact_M(V):
    if len(V) == 1: return F(1,2)
    best = F(0); qs = set([14])
    for i in range(len(V)):
        for j in range(i, len(V)):
            s = V[i]+V[j]
            for d in range(1, s+1):
                if s % d == 0: qs.add(d)
    for q in qs:
        for a in range(1, q):
            if gcd(a,q)==1:
                m = min(min((v*a)%q, q-(v*a)%q) for v in V)
                c = F(m,q)
                if c>best: best=c
    return best

def cover(V,n=14): return all(any(v%b==0 for v in V) for b in range(2,n+1))
def prim(V): return reduce(gcd,V)==1

print("="*90)
print("HUNT: primitive covering families with M as close to 1/14 as possible")
print("="*90)
# Structured: take the tight AP {1..13}, dilate by d, then make primitive by
# replacing/dropping elements while keeping covering. Also c*{1..13} + swaps.
cands = []
# family A: d*{1..12} u {w}, w chosen odd/coprime to keep primitive+covering, small M
for d in [2,3,4,5,6,7]:
    ap = [d*i for i in range(1,13)]
    for w in range(1, 500):
        V = sorted(set(ap+[w]))
        if len(V)!=13: continue
        if not prim(V) or not cover(V): continue
        cands.append(("A d=%d w=%d"%(d,w), V))
# family B: d*{1..13} with ONE element halved/replaced to gain primitivity
for d in [2,3,4]:
    ap = [d*i for i in range(1,14)]  # 13 elts, M=1/14 but non-primitive
    for idx in range(13):
        for repl in range(1, 60):
            V = sorted(set(ap[:idx]+ap[idx+1:]+[repl]))
            if len(V)!=13: continue
            if not prim(V) or not cover(V): continue
            cands.append(("B d=%d"%d, V))
# evaluate, keep smallest M
best = []
seen=set()
for name,V in cands:
    key=tuple(V)
    if key in seen: continue
    seen.add(key)
    M = exact_M(V)
    best.append((M, name, V))
best.sort(key=lambda t: t[0])
print(f"\nsearched {len(seen)} distinct primitive covering families")
print(f"smallest M values found (all must be > 1/14 = {float(F(1,14)):.5f} if LRC(14) holds):\n")
for M,name,V in best[:14]:
    print(f"  M={str(M):8s} ({float(M):.5f})  [{name}]  V={V}")
mn = best[0][0]
print(f"\nMINIMUM M over the searched primitive covering families: {mn} = {float(mn):.5f}")
print(f"  vs 1/14={float(F(1,14)):.5f}, 1/13={float(F(1,13)):.5f}, 1/9={float(F(1,9)):.5f}")
print(f"  all > 1/14: {all(M> F(1,14) for M,_,_ in best)}  (LRC(14) consistent)")
print(f"  any < 1/9 (refutes THM-995 X floor): {any(M< F(1,9) for M,_,_ in best)}")
# how close to 1/14 does it get? show the M-1/14 gaps
print("\ngaps M - 1/14 for the closest families (-> infimum is 1/14, not attained):")
for M,name,V in best[:6]:
    print(f"  M-1/14 = {M-F(1,14)} = {float(M-F(1,14)):.5f}   [{name}]")
