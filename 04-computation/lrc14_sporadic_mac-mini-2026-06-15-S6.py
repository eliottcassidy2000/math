#!/usr/bin/env python3
"""LRC(14): sporadic tight family {1..11,13,x} and {1..12,x} sweep + residue analysis.
mac-mini-2026-06-15-S6. Float-screened, exact-confirmed. stdlib only."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1,2) else 1 - r
def g(S, t): return min(nrm(v*t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1,2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1,2): C.add(F(k, d)); k += 1
    C.add(F(1,2)); return C
def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at
def gf(S, t):
    m = 2.0
    for v in S:
        r = (v*t) % 1.0; r = r if r <= 0.5 else 1.0 - r
        if r < m: m = r
    return m
def Mf_at_5_14(S):  # cheap: is M possibly tight? screen the known minimizer tau=5/14 plus generic
    # we only need a quick upper screen; use full float Mfloat-lite via 1/(2v) and pair taus
    return gf(S, 5.0/14.0)
def is_prim(S): return reduce(gcd, S) == 1
THR = F(1,14)

# Sweep {1..11,13,x}: exact M only when float at 5/14 == ~1/14 (necessary for tightness via this vertex)
print("--- {1..11,13,x} x=14..600: tight (M==1/14) ---")
t1=[]
for x in range(14, 601):
    S=sorted(set(list(range(1,12))+[13,x]))
    if len(S)!=13: continue
    m,_=M(S)
    if m<THR:
        print("  COUNTEREXAMPLE", S, m);
    if m==THR and is_prim(S): t1.append(x)
print("  tight x:", t1)

print("--- {1..12,x} x=13..600: tight ---")
t2=[]
for x in range(13, 601):
    S=sorted(set(list(range(1,13))+[x]))
    if len(S)!=13: continue
    m,_=M(S)
    if m<THR: print("  COUNTEREXAMPLE", S, m)
    if m==THR and is_prim(S): t2.append(x)
print("  tight x:", t2)

# residue patterns at tau=5/14
print("\n--- residues at tau=5/14 ---")
for S in [list(range(1,14)), list(range(1,12))+[13,24]]:
    print(" S=",sorted(S))
    for v in sorted(S):
        print(f"    v={v:2d}: ||5v/14|| = {nrm(v*F(5,14))}")
