#!/usr/bin/env python3
"""
lrc14_disprove_counterexample_kps  (kind-pasteur 2026-06-16, DISPROVE side, part B)

Hunt for a GENUINE LRC(14) counterexample:
  a 13-set S with max-min := max_tau min_v ||v tau|| < 1/14.
LRC(14) holds  <=>  every 13-set has max-min >= 1/14.
Known tight configs (AP, sporadics) have max-min = 1/14 EXACTLY (equality).
A counterexample would have max-min STRICTLY < 1/14.

max-min is computed EXACTLY: min_v ||v tau|| is piecewise linear; its max over
[0,1) is attained at a breakpoint (kink j/v) or a crossing of two ||.|| graphs.
We enumerate the exact candidate tau set and evaluate.

We scan:
  - the AP and all single/double perturbations (tight-like neighborhood),
  - the tight sporadics and their neighbors,
  - dense low-value sets {1..k, ...} which maximize coverage.
"""
from fractions import Fraction as Fr
from itertools import combinations
from math import gcd
from functools import reduce

THRESH=Fr(1,14)

def frac_norm(x):
    f=x-(x.numerator//x.denominator)
    return f if f<=1-f else 1-f

def maxmin_exact(S):
    S=sorted(set(S))
    cands=set()
    for v in S:
        for j in range(v+1):
            cands.add(Fr(j,v))
    n=len(S)
    for i in range(n):
        a=S[i]
        for jx in range(i+1,n):
            b=S[jx]
            for p in range(a+1):
                ap=Fr(p,1)
                for q in range(b+1):
                    if a!=b:
                        t=Fr(p-q,a-b)
                        if 0<=t<=1: cands.add(t)
                    s=a+b
                    t2=Fr(p+q,s)
                    if 0<=t2<=1: cands.add(t2)
    best=Fr(-1); bt=None
    for t in cands:
        m=min(frac_norm(v*t) for v in S)
        if m>best: best=m; bt=t
    return best,bt

def lcm(S): return reduce(lambda a,b:a*b//gcd(a,b),S,1)
def is_prim(S): return reduce(gcd,S)==1

BASE=list(range(1,14))
results=[]   # (maxmin, S, tag)

def consider(S,tag):
    Sf=tuple(sorted(set(S)))
    if len(Sf)!=13: return
    mm,t=maxmin_exact(Sf)
    results.append((mm,Sf,tag,t))

print("="*78)
print("[1] AP single perturbations e->w (w<=120) and the AP itself")
print("="*78)
consider(BASE,"AP")
for e in BASE:
    for w in range(14,121):
        if w in BASE and w!=e: continue
        consider([x for x in BASE if x!=e]+[w], f"single {e}->{w}")

print("[2] AP double perturbations e1,e2->w1,w2 (w<=60)")
WR=list(range(14,61))
for e1,e2 in combinations(BASE,2):
    rest=[x for x in BASE if x not in (e1,e2)]
    for w1,w2 in combinations(WR,2):
        if w1 in rest or w2 in rest: continue
        consider(rest+[w1,w2], f"double {e1},{e2}->{w1},{w2}")

print("[3] tight sporadics + their single perturbations")
SPOR=[1,2,3,4,5,6,7,8,9,10,11,13,24]
consider(SPOR,"sporadic")
for e in SPOR:
    for w in range(14,80):
        if w in SPOR and w!=e: continue
        consider([x for x in SPOR if x!=e]+[w], f"spor {e}->{w}")

print("[4] dense low-value sets: {1..12, w} maximizing coverage, w<=120")
for w in range(13,121):
    consider(list(range(1,13))+[w], f"dense 1..12,{w}")
# {1..11, a, b}
for a,b in combinations(range(12,40),2):
    consider(list(range(1,12))+[a,b], f"1..11,{a},{b}")

results.sort()
print("\n"+"="*78)
print("SMALLEST max-min found (counterexample iff < 1/14):")
print("="*78)
shown=0; seen=set()
for mm,S,tag,t in results:
    if S in seen: continue
    seen.add(S)
    flag=" <<< COUNTEREXAMPLE (max-min < 1/14) " if mm<THRESH else ("  (= 1/14, tight)" if mm==THRESH else "")
    print(f"  maxmin={float(mm):.9g}={mm}{flag}  tau={t}")
    print(f"        S={S} (lcm={lcm(S)}, prim={is_prim(S)}) [{tag}]")
    shown+=1
    if shown>=20: break

print("\nGlobal min max-min:", min(r[0] for r in results), "= 1/14?" , min(r[0] for r in results)==THRESH)
n_below=sum(1 for mm,_,_,_ in results if mm<THRESH)
print("configs with max-min < 1/14 (GENUINE COUNTEREXAMPLES):", n_below)
