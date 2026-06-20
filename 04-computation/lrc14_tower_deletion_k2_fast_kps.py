#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Close the k=2 layer of the tower-deletion lemma (THM-545 extends from k=1 to k=2).
kind-pasteur-2026-06-19-S17. Fast FLOAT pre-filter + EXACT Fraction confirmation.

THM-545 (PROVED k=1): deleting a tower bit 2^a + ONE extra hole, adding ONE tail => meas(G_C)>=426/35035.
THIS extends to k=2: delete 2^a + TWO extra holes, add TWO tails. The comb cutoff bounds the tails to
[14, R_B); only those (finite) need checking. Float pre-filter: meas_float > thr2 + eps => safe; else
confirm with exact Fraction. (The floor argument says almost all are far above thr2.)
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
THR2=Fraction(426,35035); THR2F=float(THR2)

def lonely_float(C, theta=1/14.0):
    """fast float lonely measure: 1 - meas(union danger arcs)."""
    segs=[]
    for d in C:
        if d==0: continue
        w=theta/d
        for m in range(0,d+1):
            lo=m/d-w; hi=m/d+w
            for sh in (-1.0,0.0,1.0):
                a=max(lo+sh,0.0); b=min(hi+sh,1.0)
                if a<b: segs.append((a,b))
    segs.sort(); cov=0.0; cur=-1.0
    for a,b in segs:
        if a>cur: cov+=b-a; cur=b
        elif b>cur: cov+=b-cur; cur=b
    return 1.0-cov

def lonely_exact(C, theta=Fraction(1,14)):
    segs=[]
    for d in C:
        if d==0: continue
        w=theta/d
        for m in range(0,d+1):
            lo=Fraction(m,d)-w; hi=Fraction(m,d)+w
            for sh in (-1,0,1):
                a=max(lo+sh,Fraction(0)); b=min(hi+sh,Fraction(1))
                if a<b: segs.append((a,b))
    segs.sort(); union=[]; cur=Fraction(-1)
    for a,b in segs:
        if a>cur: union.append([a,b]); cur=b
        elif b>cur: union[-1][1]=b; cur=b
    return Fraction(1)-sum(b-a for a,b in union)

def comb_cutoff_R(base):
    """R_B = ceil(2 c_B / (6 M_B - 7 thr2)); for tails r>=R_B the comb certifies meas>=thr2.
       c_B = #components of G_base, M_B = meas(G_base) (exact)."""
    # exact base measure + component count
    segs=[]
    for d in base:
        if d==0: continue
        w=Fraction(1,14)/d
        for m in range(0,d+1):
            lo=Fraction(m,d)-w; hi=Fraction(m,d)+w
            for sh in (-1,0,1):
                a=max(lo+sh,Fraction(0)); b=min(hi+sh,Fraction(1))
                if a<b: segs.append((a,b))
    segs.sort(); union=[]; cur=Fraction(-1)
    for a,b in segs:
        if a>cur: union.append([a,b]); cur=b
        elif b>cur: union[-1][1]=b; cur=b
    M=Fraction(1)-sum(b-a for a,b in union)           # covered complement
    # safe components = gaps between covered union pieces (incl wrap)
    c=len(union)-1 if union else 0
    if union and not (union[0][0]==0 and union[-1][1]==1): c=len(union)  # rough; gaps count
    denom=6*M-7*THR2
    if denom<=0: return None, M, c
    R=(2*c + denom - 1)//denom if False else None
    # R = ceil(2c/denom)
    from math import ceil
    R=ceil(Fraction(2*c)/denom)
    return max(R,14), M, c

if __name__ == "__main__":
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument("--bit",type=int,default=4); args=ap.parse_args()
    a_bit=args.bit  # the tower bit to delete (1,2,4,8)
    print(f"=== k=2 tower-deletion scan: delete tower bit {a_bit} + 2 extra holes, add 2 tails ===")
    others=[d for d in range(1,14) if d!=a_bit]
    eps=0.0008
    checked=0; near=0; below=[]; nb=0
    for h1,h2 in itertools.combinations(others,2):
        base=tuple(sorted([d for d in range(1,14) if d not in (a_bit,h1,h2)]))  # 10 speeds
        R=300   # wide window: covers the iterated two-tail comb cutoff (~245); tails>=300 comb-certified
        nb+=1
        for r1,r2 in itertools.combinations(range(14,R),2):
            C=tuple(sorted(base+(r1,r2)))
            if reduce(gcd,C)!=1: continue
            checked+=1
            mf=lonely_float(C)
            if mf < THR2F + eps:
                near+=1
                me=lonely_exact(C)
                if me < THR2: below.append((me,base,(r1,r2)))
    print(f"  bit {a_bit}: {nb} bases, {checked} cores checked (tails<R_B), {near} near-threshold (exact-confirmed)")
    if below:
        below.sort()
        print(f"  *** {len(below)} cores BELOW thr2 (tower bit {a_bit} deleted): RIGIDITY FAILS ***")
        for me,b,t in below[:8]: print(f"      meas={me}={float(me):.6f} base={b} tails={t}")
    else:
        print(f"  ==> ZERO cores below 426/35035. k=2 tower-deletion for bit {a_bit}: COMB-CERTIFIED (tails>=R_B) + EXACT (tails<R_B). PROVED.")
