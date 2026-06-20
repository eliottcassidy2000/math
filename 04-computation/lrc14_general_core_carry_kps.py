#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GENERAL 12-CORE carry-conservation test (OPEN-Q-108, beyond AP-tail).
codex's lonely-measure route.  EXACT rationals (fractions.Fraction).

Question (HYP-2661 generalized): for ARBITRARY primitive positive 12-subsets C (a "12-core",
elements up to some cap), is it true that
        meas(G_C) < 426/35035  ==>  {1,2,4,8} subset C   (full dyadic-1 tower / carry[1]=15) ?

meas(G_C) = meas{ t in [0,1) : ||d t|| > 1/14  for all d in C }   (lonely measure, gap theta=1/14).

If TRUE generally: the dyadic-1 tower is a UNIVERSAL mouth-owner -> sub-threshold census reduces
to {1,2,4,8}-containing cores.  If FALSE: find the counterexample + its alternative mouth-owner.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd

if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')

THETA = Fraction(1,14)
THR1  = Fraction(7,858)      # the mouth (drop-6 global min)
THR2  = Fraction(426,35035)  # AP one-hole second value

def lonely_measure(C, theta=THETA):
    """meas{t in [0,1): ||d t|| > theta for all d in C}, exact rationals."""
    segs=[]
    for d in C:
        if d==0: continue
        w=theta/d
        for m in range(0, d+1):
            lo=Fraction(m,d)-w; hi=Fraction(m,d)+w
            for shift in (-1,0,1):
                a=lo+shift; b=hi+shift
                a2=max(a,Fraction(0)); b2=min(b,Fraction(1))
                if a2<b2: segs.append((a2,b2))
    segs.sort()
    cur=Fraction(-1); union=[]
    for a,b in segs:
        if a>cur: union.append([a,b]); cur=b
        elif b>cur: union[-1][1]=b; cur=b
    covered=sum(b-a for a,b in union)
    return Fraction(1)-covered

def carry_profile(C):
    d={}
    for e in C:
        if e==0: continue
        m=e; a=0
        while m%2==0: m//=2; a+=1
        d[m]=d.get(m,0)+2**a
    return dict(sorted(d.items()))

def tower_full(C):
    """is the full dyadic-1 tower {1,2,4,8} subset C  (carry[1]==15)?"""
    return all(b in C for b in (1,2,4,8))

def scan(cap, sample_only=None, seed=0):
    """All primitive 12-subsets of [1,cap]; if too many, random sample."""
    universe=list(range(1,cap+1))
    total=0
    below=[]                  # rows with meas<THR2
    below_no_tower=[]         # the dangerous ones: below but tower NOT full -> counterexamples
    import random
    rng=random.Random(seed)
    if sample_only is None:
        combos=itertools.combinations(universe,12)
        for C in combos:
            if reduce(gcd,C)!=1: continue
            total+=1
            L=lonely_measure(C)
            if L<THR2:
                below.append((L,C))
                if not tower_full(C):
                    below_no_tower.append((L,C))
    else:
        seen=set()
        tries=0
        while len(seen)<sample_only and tries<sample_only*50:
            tries+=1
            C=tuple(sorted(rng.sample(universe,12)))
            if C in seen: continue
            if reduce(gcd,C)!=1: continue
            seen.add(C); total+=1
            L=lonely_measure(C)
            if L<THR2:
                below.append((L,C))
                if not tower_full(C):
                    below_no_tower.append((L,C))
    below.sort(); below_no_tower.sort()
    return total, below, below_no_tower

if __name__=="__main__":
    print(f"THR1 (mouth)={THR1}={float(THR1):.6f}   THR2 (AP 2nd)={THR2}={float(THR2):.6f}")
    print("Testing HYP-2661 generalized: meas<THR2 ==> {1,2,4,8} subset C, on GENERAL 12-cores.\n")

    # 1) EXHAUSTIVE over [1, cap] for caps where C(cap,12) is feasible.
    for cap in [13,14,15,16,17,18]:
        from math import comb
        nC=comb(cap,12)
        print(f"--- EXHAUSTIVE cap={cap}: C({cap},12)={nC} subsets ---")
        total,below,bnt=scan(cap)
        print(f"    primitive scanned={total}; below THR2: {len(below)}; below & NO full tower: {len(bnt)}")
        for L,C in below[:8]:
            print(f"      meas={float(L):.6f}={L}  C={C}  tower={'1248' if tower_full(C) else 'BROKEN'} carry1={carry_profile(C).get(1)}")
        if bnt:
            print("    *** COUNTEREXAMPLES (below THR2, tower NOT full): ***")
            for L,C in bnt[:20]:
                print(f"      meas={float(L):.6f}={L}  C={C}  carry={carry_profile(C)}")
        print()
