#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The CARRY CONSERVATION law for LRC14 AP-tail mouth-retention (HYP-2654/2659, user hint).
kind-pasteur-2026-06-19-S16.  EXACT rationals (lonely measure at gap 1/14).

Glaisher carry profile of a speed set C: shell(m) = sum of 2^a over speeds 2^a*m in C (m odd).
meas(G_C) = meas{t in [0,1): ||d t|| > 1/14 for all d in C}  (the lonely set; codex's object).

User's claim (to make precise + prove): the one-tail exception {5:+2} (pure carry in an existing
shell) keeps the FULL drop-6 mouth (old_survivor=7/858) so stays >= the AP one-hole 2nd value;
the two-tail min {1:-4,5:+2,23:+2} DAMAGES the mouth (the {1:-4} removes a shell-1 tower bit) but
its measure already PAYS the 2nd threshold 426/35035. The conservation: you cannot damage the
mouth without paying. GOAL: find the exact carry->measure law that forces this.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def lonely_measure(C, theta=Fraction(1,14)):
    """meas{t in [0,1): ||d t|| > theta for all d in C}, exact."""
    # danger arcs for speed d: union_m [m/d - theta/d, m/d + theta/d]; complement, intersect all.
    arcs=[]
    for d in C:
        if d==0: continue
        w=theta/d
        for m in range(0, d+1):
            lo=Fraction(m,d)-w; hi=Fraction(m,d)+w
            arcs.append((lo,hi))
    # union of danger arcs on the circle [0,1): merge (handle wrap by working on [0,1) with mod)
    # represent each arc clipped/split at 0 and 1
    segs=[]
    for lo,hi in arcs:
        lo=lo; hi=hi
        # split into [0,1) pieces
        # shift to canonical: consider integer translates intersecting [0,1]
        k0=lo.numerator//lo.denominator if lo>=0 else -((-lo.numerator)//lo.denominator+ (1 if (-lo.numerator)%lo.denominator else 0))
        # simpler: just add the raw interval and all its +-1 translates that meet [0,1]
        for shift in (-1,0,1):
            a=lo+shift; b=hi+shift
            a2=max(a,Fraction(0)); b2=min(b,Fraction(1))
            if a2<b2: segs.append((a2,b2))
    segs.sort()
    # union length
    covered=Fraction(0); cur=Fraction(-1)
    union=[]
    for a,b in segs:
        if a>cur:
            union.append([a,b]); cur=b
        elif b>cur:
            union[-1][1]=b; cur=b
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

def ap_tail_core(holes, tails):
    return tuple(sorted([d for d in range(1,14) if d not in holes]+list(tails)))

if __name__ == "__main__":
    print("=== KNOWN ROWS: carry profile + exact lonely measure (verify the hint data) ===")
    rows = {
      "drop-6 (THM-541)": ap_tail_core({6}, []),
      "one-tail exc (THM-543)": ap_tail_core({6,10}, [20]),
      "two-tail min (THM-544)": ap_tail_core({4,6,10}, [20,46]),
      "drop-12 (AP 2nd)": ap_tail_core({12}, []),
    }
    base=carry_profile(ap_tail_core({6},[]))
    thr1=Fraction(7,858); thr2=Fraction(426,35035)
    for name,C in rows.items():
        L=lonely_measure(C); cp=carry_profile(C)
        delta={m:cp.get(m,0)-base.get(m,0) for m in set(cp)|set(base) if cp.get(m,0)!=base.get(m,0)}
        print(f"  {name}: meas={L}={float(L):.6f}  carry={cp}")
        print(f"      delta_vs_drop6={ {m:('+' if v>0 else '')+str(v) for m,v in delta.items()} }  >=thr2(426/35035)? {L>=thr2}")
    print(f"\n  drop-6 mouth 7/858={float(thr1):.6f}; AP-2nd 426/35035={float(thr2):.6f}")

    print("\n=== EXPLORE: does the carry profile determine whether meas < 426/35035 ? ===")
    print("    scan two-tail rows (remove 2-3 holes, add 2 tails<=30), group by carry-delta signature")
    below=[]; cnt=0
    for nh in [2,3]:
        for holes in itertools.combinations(range(1,14), nh):
            for tails in itertools.combinations(range(14,31), nh-1+1 if False else 2):
                C=ap_tail_core(set(holes), tails)
                if len(C)!=12: continue
                if reduce(gcd,C)!=1: continue
                cnt+=1
                L=lonely_measure(C)
                if L<thr2:
                    cp=carry_profile(C); delta={m:cp.get(m,0)-base.get(m,0) for m in set(cp)|set(base) if cp.get(m,0)!=base.get(m,0)}
                    below.append((L,holes,tails,delta))
    below.sort()
    print(f"    scanned {cnt} two-tail rows; {len(below)} below 426/35035:")
    for L,h,t,d in below[:12]:
        print(f"      meas={float(L):.6f}={L} holes={h} tails={t} carry-delta={d}")
    if not below: print("      NONE below threshold (consistent with THM-544).")
