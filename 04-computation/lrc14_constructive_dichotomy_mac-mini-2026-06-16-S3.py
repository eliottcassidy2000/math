#!/usr/bin/env python3
"""
lrc14_constructive_dichotomy — mac-mini-2026-06-16-S3

A CONSTRUCTIVE proof of half of LRC(14), and the isolation of the hard half.

LEMMA (constructive, elementary). For ANY 13-set S of distinct positive integers
with NO multiple of 14, tau* = 1/14 is a lonely point:
    ||v/14|| = min(v mod 14, 14-(v mod 14))/14 >= 1/14   <=>   v not= 0 (mod 14),
so min_v ||v/14|| >= 1/14, hence M(S) >= 1/14.  => LRC(14) holds for all such S.

By scale-invariance M(S)=M(cS) (tau->c tau measure-preserving), the only shapes
NOT covered are primitive sets containing a multiple of 14 = exactly C'(14).
So LRC(14)  <=>  M(S) >= 1/14 for primitive S containing a multiple of 14.
This recovers the repo reduction LRC(14) <= C'(14) CONSTRUCTIVELY (explicit witness),
and tells the disproof search to live ENTIRELY in C'(14).

This script: (1) verifies the lemma over a large random sample (0 failures expected);
(2) searches C'(14) (S contains 14m) for the min gap M -- is it always >= 1/14
(no counterexample) and is the inf > 1/14 (strict looseness)?; (3) records the
lonely points tau* of C'(14) extremizers to look for a SECOND explicit construction.
"""
from fractions import Fraction as F
from itertools import combinations
import random

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r; return r if r<=F(1,2) else 1-r
def g(S,t): return min(nrm(v*t) for v in S)
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def M(S):
    b=F(0); at=None
    for t in cand(S):
        v=g(S,t)
        if v>b: b=v; at=t
    return b,at

ONE14=F(1,14)
print("="*72)
print("[1] LEMMA CHECK: S with no multiple of 14  =>  tau=1/14 lonely, M>=1/14")
print("="*72)
random.seed(3); fails=0; n=4000
for _ in range(n):
    S=random.sample([x for x in range(1,120) if x%14!=0],13)
    # witness tau=1/14
    wit=g(S,F(1,14))
    if wit<ONE14: fails+=1; print("  WITNESS FAIL", S, wit)
print(f"  random no-mult-14 13-sets tested: {n};  tau=1/14 witness failures: {fails}")
print(f"  (0 failures => M(S)>=1/14 PROVED constructively for all no-mult-14 S)")
# also confirm via exact M on a few that M>=1/14
for _ in range(5):
    S=sorted(random.sample([x for x in range(1,60) if x%14!=0],13))
    m,t=M(S); assert m>=ONE14
print("  exact M>=1/14 confirmed on samples.")

print("\n"+"="*72)
print("[2] THE HARD HALF: C'(14) = primitive S containing a multiple of 14.")
print("    Search for min gap M over C'(14): any M<1/14 (counterexample)? inf>1/14?")
print("="*72)
from math import gcd
from functools import reduce
def isprim(S): return reduce(gcd,S)==1
# interior-drop family {1..13}\{j} u {14m} (all C'(14)); + broader: 12 small + one 14m
best=(F(1),None); ce=[]
# (a) interior-drop
for j in range(1,14):
    A=[x for x in range(1,14) if x!=j]
    for m in range(1,16):
        S=A+[14*m]
        if len(set(S))!=13: continue
        if not isprim(S): continue
        mm,t=M(S)
        if mm<best[0]: best=(mm,(tuple(sorted(S)),t))
        if mm<ONE14: ce.append((mm,tuple(sorted(S)),t))
# (b) broader: 12-subset of [1,20]\{14} plus one multiple of 14
pool=[x for x in range(1,21) if x%14!=0]
for core in combinations(pool,12):
    for w in (14,28,42,56):
        if w in core: continue
        S=list(core)+[w]
        if not isprim(S): continue
        mm,t=M(S)
        if mm<best[0]: best=(mm,(tuple(sorted(S)),t))
        if mm<ONE14: ce.append((mm,tuple(sorted(S)),t))
mm,(S,t)=best
print(f"  min gap M over C'(14) search = {mm} = {float(mm):.6f}  (1/14={float(ONE14):.6f})")
print(f"    at S={list(S)}, lonely tau*={t}")
print(f"  M = 1/14 exactly? {mm==ONE14}   M > 1/14 (strictly loose)? {mm>ONE14}")
print(f"  COUNTEREXAMPLES (M<1/14): {ce if ce else 'NONE'}")

print("\n"+"="*72)
print("[3] Lonely points tau* of C'(14) interior-drop extremizers (seek 2nd construction)")
print("="*72)
for j in (4,6,12):
    A=[x for x in range(1,14) if x!=j]
    print(f"  j={j}:")
    for m in (1,2,4,7):
        S=A+[14*m]
        if len(set(S))!=13: continue
        mm,t=M(S)
        print(f"    +{14*m}: M={mm}={float(mm):.5f}  tau*={t}  (14*tau*={14*t}, denom={t.denominator})")
print("\n  => the no-mult-14 case is DONE (tau=1/14). The whole LRC(14) difficulty is C'(14),")
print("     where the lonely point MOVES off 1/14 (denominators above) — find its law.")
