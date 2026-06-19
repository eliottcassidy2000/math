#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL part 3: exhaustive search over primitive sets E (0 in E, |E|=k,
bounded spread) computing L_y EXACTLY. Verify:
  (a) consec is the global maximizer over bounded-spread sets (HYP-2607 finite part),
  (b) consec L_y < cap_k,
  (c) wide primitive sets stay well below cap (Step 3),
  (d) ANY set found with L_y >= cap would BREAK the claim.
"""
import sys, itertools
from fractions import Fraction as F
from functools import reduce
from math import gcd, comb

def N_at(E, x):
    hit=set()
    for e in E:
        v=e*x; v=v-(v.numerator//v.denominator)
        hit.add((v.numerator*7)//v.denominator)
    return sum(1 for j in range(1,7) if j not in hit)

def dist_p(E):
    E=sorted(set(E))
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(F(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1)
    p=[F(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        p[N_at(E,(lo+hi)/2)]+=(hi-lo)
    return p

def g_poly(k):
    g=[]
    for t in range(7):
        if k==8: g.append(F((t-1)*(t-2)*(t-4)*(t-5),40))
        elif k in (9,10): g.append(F(-(t-2)*(t-3)*(t-6),36))
        else: g.append(F((t-3)*(t-4),12))
    return g

def L_y(E,k):
    p=dist_p(E); g=g_poly(k)
    return sum(p[t]*g[t] for t in range(7))

caps={8:F(3815,10000),9:F(49426,100000),10:F(6044,10000)}

def is_primitive(E):
    return reduce(gcd, E)==1

# k=8 exhaustive over spread<=14 (claim says consec UNIQUE over spread<=14)
k=8
g=g_poly(k); cap=caps[k]
print(f"=== k={k} exhaustive primitive sets, 0 in E, |E|={k}, max<=14 ===")
consecE=list(range(k))
Lc=L_y(consecE,k)
print(f"consec L_y = {Lc} = {float(Lc):.6f}  (cap={float(cap):.6f})  consec<cap? {Lc<cap}")
best=Lc; bestE=consecE; ties=[]; over_cap=[]
cnt=0
# choose 7 more elements from 1..14
for combo in itertools.combinations(range(1,15),7):
    E=[0]+list(combo)
    if not is_primitive(E): continue
    cnt+=1
    L=L_y(E,k)
    if L>=cap: over_cap.append((E,L))
    if L>best:
        best=L; bestE=E
    elif L==best and E!=consecE:
        ties.append(E)
print(f"primitive sets checked: {cnt}")
print(f"global best L_y = {best} = {float(best):.6f}  at E={bestE}")
print(f"best == consec? {bestE==consecE}")
print(f"ties with best (besides consec): {ties[:5]}{'...' if len(ties)>5 else ''}  count={len(ties)}")
print(f"any set with L_y >= cap? {len(over_cap)} -> {over_cap[:3]}")
