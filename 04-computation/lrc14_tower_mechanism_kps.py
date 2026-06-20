#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WHY the dyadic-1 tower {1,2,4,8} owns the mouth.  Mechanism diagnostic.
For a 12-core C the lonely set lives in [0,1).  The element d=1 forbids ||t||>1/14 fails only
near 0,1 -> arc [0,1/14] u [1-1/14,1] removed.  d=2 removes around 0,1/2,1.  d=4 around k/4.
d=8 around k/8.  The tower {1,2,4,8} = all dyadic rationals of denominator | 8, each carrying a
radius theta/d = 1/(14d).  Together they remove a comb at every k/8 with shrinking radii.

We measure: how much of [0,1) does the tower ALONE remove, vs the surviving set, and where the
drop-6 survivors sit (THM-541: 4 components owners R(13,2)->L(12,2) etc).  Then we test the
LOCAL claim: any core achieving meas<426/35035 must have its surviving mass pinned in the gaps
that ONLY the full tower leaves open.  Removing any of {1,2,4,8} widens a gap -> mass jumps up.
EXACT.  kind-pasteur-2026-06-19.
"""
import sys
from fractions import Fraction
if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')
THETA=Fraction(1,14)

def lonely_set(C, theta=THETA):
    """return surviving components as list of [a,b] in [0,1)."""
    segs=[]
    for d in C:
        w=theta/d
        for m in range(0,d+1):
            lo=Fraction(m,d)-w; hi=Fraction(m,d)+w
            for s in (-1,0,1):
                a=max(lo+s,Fraction(0)); b=min(hi+s,Fraction(1))
                if a<b: segs.append((a,b))
    segs.sort(); cur=Fraction(-1); U=[]
    for a,b in segs:
        if a>cur: U.append([a,b]); cur=b
        elif b>cur: U[-1][1]=b; cur=b
    # complement of U in [0,1)
    surv=[]; prev=Fraction(0)
    for a,b in U:
        if a>prev: surv.append([prev,a])
        prev=max(prev,b)
    if prev<1: surv.append([prev,Fraction(1)])
    return surv

def meas(surv): return sum(b-a for a,b in surv)

if __name__=="__main__":
    drop6=(1,2,3,4,5,7,8,9,10,11,12,13)
    print("=== drop-6 surviving components (THM-541 says 4) ===")
    S=lonely_set(drop6)
    print(f"   meas={meas(S)}={float(meas(S)):.6f}  #components={len(S)}")
    for a,b in S:
        print(f"     [{a}, {b}]  width={b-a}={float(b-a):.6f}  center~{float((a+b)/2):.5f}")

    print("\n=== tower {1,2,4,8} ALONE: surviving set ===")
    St=lonely_set((1,2,4,8))
    print(f"   meas={meas(St)}={float(meas(St)):.6f}  #comp={len(St)}")
    for a,b in St[:20]:
        print(f"     [{float(a):.5f},{float(b):.5f}] w={float(b-a):.5f}")

    print("\n=== drop-6 WITH one tower bit removed (widen test) ===")
    for kill in (1,2,4,8):
        C=tuple(x for x in drop6 if x!=kill)
        # not primitive necessarily; just measure
        S2=lonely_set(C)
        print(f"   drop-6 minus {kill}: meas={meas(S2)}={float(meas(S2)):.6f}  #comp={len(S2)}  (vs 7/858={float(Fraction(7,858)):.6f})")

    print("\n=== drop-6 survivors: which dyadic comb teeth bound each survivor? ===")
    # report the endpoints as m/d to see which element is the 'owner'
    for a,b in S:
        def owners(x):
            r=[]
            for d in drop6:
                for m in range(0,d+1):
                    if abs(Fraction(m,d)-x)==THETA/d: r.append(f"{m}/{d}")
            return r
        print(f"   [{a},{b}]: left-owner {owners(a)} right-owner {owners(b)}")
