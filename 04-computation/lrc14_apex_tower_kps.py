#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WHY 4 levels {1,2,4,8} and the apex-prime-7 / n=14 connection.
theta=1/N with N=14=2*7.  Tower = powers of 2 that are < N: {1,2,4,8} since 8<14<16.
Depth = floor(log2(N-1))+1 = floor(log2 14)= 3 -> exponents 0,1,2,3 -> 4 levels.

Test the META-claim across N: for the lonely measure at gap 1/N over k-cores from [1,N-1],
is the global-min core ("the mouth") always the (N-1)\{apex-ish hole}, and is its mouth owned
by the dyadic tower {2^a < N}?  Run N=2p for small primes p (the LRC(2p) family) and inspect
the tower depth + whether dropping a tower bit always blows the measure above the 2nd value.
EXACT.  kind-pasteur-2026-06-19.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd, log2, floor
if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')

def lonely_measure(C, N):
    theta=Fraction(1,N); segs=[]
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
    return Fraction(1)-sum(b-a for a,b in U)

def tower(N): return [1<<a for a in range(0, N.bit_length()) if (1<<a)<N]

if __name__=="__main__":
    print("apex tower per N (powers of 2 < N):")
    for N in [6,10,14,22,26]:
        print(f"  N={N}=2*{N//2}: tower={tower(N)} depth={len(tower(N))} (floor log2(N-1)+1={floor(log2(N-1))+1})")
    print()
    # For the genuine LRC family the core size is N-2 (drop one hole from {1..N-1}).
    # LRC(14): cores of size 12 from {1..13}. Generalize: size N-2 from {1..N-1}, drop 1.
    for N in [6,10,14]:
        k=N-2
        print(f"=== N={N}: one-hole {k}-cores from {{1..{N-1}}} (drop 1), exact lonely measure ===")
        rows=[]
        for hole in range(1,N):
            C=tuple(d for d in range(1,N) if d!=hole)
            if reduce(gcd,C)!=1: continue
            L=lonely_measure(C,N)
            tw=all(b in C for b in tower(N))
            rows.append((L,hole,tw))
        rows.sort()
        thr1=rows[0][0]; thr2=rows[1][0] if len(rows)>1 else None
        for L,hole,tw in rows:
            mark='<-MOUTH' if L==thr1 else ('<-2nd' if L==thr2 else '')
            print(f"    drop {hole}: meas={L}={float(L):.6f} tower_full={tw} {mark}")
        # does dropping a tower bit from the MOUTH core blow it above thr2?
        mouth_hole=rows[0][1]
        Cm=tuple(d for d in range(1,N) if d!=mouth_hole)
        print(f"    MOUTH=drop {mouth_hole}; thr1={float(thr1):.6f} thr2={float(thr2):.6f}")
        for b in tower(N):
            if b in Cm:
                C2=tuple(x for x in Cm if x!=b)
                L2=lonely_measure(C2,N)
                print(f"      remove tower bit {b}: meas={float(L2):.6f} {'> thr2 (blows up)' if L2>=thr2 else 'still < thr2 !!'}")
        print()
