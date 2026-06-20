#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GENERAL 12-CORE carry test v2 (OPEN-Q-108).  Larger caps: exhaustive to 22, sample to 40.
Tests HYP-2661 generalized: meas(G_C)<426/35035 ==> {1,2,4,8} subset C.
EXACT rationals.  kind-pasteur-2026-06-19.
"""
import sys, itertools, random
from fractions import Fraction
from functools import reduce
from math import gcd, comb

if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')

THETA=Fraction(1,14); THR2=Fraction(426,35035); THR1=Fraction(7,858)

def lonely_measure(C, theta=THETA):
    segs=[]
    for d in C:
        if d==0: continue
        w=theta/d
        for m in range(0,d+1):
            lo=Fraction(m,d)-w; hi=Fraction(m,d)+w
            for s in (-1,0,1):
                a=max(lo+s,Fraction(0)); b=min(hi+s,Fraction(1))
                if a<b: segs.append((a,b))
    segs.sort()
    cur=Fraction(-1); U=[]
    for a,b in segs:
        if a>cur: U.append([a,b]); cur=b
        elif b>cur: U[-1][1]=b; cur=b
    return Fraction(1)-sum(b-a for a,b in U)

def carry_profile(C):
    d={}
    for e in C:
        m=e;a=0
        while m%2==0: m//=2;a+=1
        d[m]=d.get(m,0)+2**a
    return dict(sorted(d.items()))
def tower_full(C): return all(b in C for b in (1,2,4,8))

# QUICK PRE-FILTER: meas(G_C) <= meas(G_{C'}) whenever C' subset C (more constraints = smaller set).
# So meas(G_C) >= product/independent lower bound? Use a cheap necessary screen to skip far rows:
# meas(G_C) <= meas(G_{single smallest few}) won't help. Instead just compute; but prune by
# noting meas only depends on C. We compute fully (exact) but skip non-primitive early.

def run(cap, sample=None, seed=0, label=""):
    rng=random.Random(seed)
    below=[]; bnt=[]; total=0
    if sample is None:
        it=itertools.combinations(range(1,cap+1),12)
        for C in it:
            if reduce(gcd,C)!=1: continue
            total+=1
            L=lonely_measure(C)
            if L<THR2:
                below.append((L,C))
                if not tower_full(C): bnt.append((L,C))
    else:
        seen=set(); tries=0
        while len(seen)<sample and tries<sample*80:
            tries+=1
            C=tuple(sorted(rng.sample(range(1,cap+1),12)))
            if C in seen: continue
            if reduce(gcd,C)!=1: continue
            seen.add(C); total+=1
            L=lonely_measure(C)
            if L<THR2:
                below.append((L,C))
                if not tower_full(C): bnt.append((L,C))
    below.sort(); bnt.sort()
    print(f"=== {label} cap={cap} {'EXHAUSTIVE' if sample is None else f'SAMPLE n={sample}'}: primitive={total} below_THR2={len(below)} below_NO_tower={len(bnt)} ===")
    for L,C in below[:12]:
        print(f"   meas={float(L):.6f}={L} C={C} {'1248' if tower_full(C) else '*BROKEN*'}")
    if bnt:
        print("   *** COUNTEREXAMPLES below THR2 without full tower: ***")
        for L,C in bnt[:30]:
            print(f"      meas={float(L):.6f}={L} C={C} carry={carry_profile(C)}")
    print()
    return below, bnt

if __name__=="__main__":
    print(f"THR2=426/35035={float(THR2):.6f}  THR1=7/858={float(THR1):.6f}\n")
    run(22, label="exhaustive22")
    run(30, sample=40000, seed=1, label="rand30")
    run(40, sample=40000, seed=2, label="rand40")
    run(50, sample=30000, seed=3, label="rand50")
