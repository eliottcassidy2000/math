#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fast version of the Candidate-2 envelope test using meet-in-the-middle on the
relation constraint sum_j n_j e_j = 0. Splits the 7 nonzero offsets into two halves,
enumerates partial sums, and matches.  Computes:
  ENV(L)  = sum over n in Lambda(E), |n|_inf<=L, all coords nonzero&not-7-mult,
            of prod c1/|n_j|   (support exactly 7 here since k=8 -> 7 nonzero offsets;
            support>=6 also captured by allowing ONE coord = 0).
We track ENV(L) for growing L to test the claimed (pi^2/6)^3 < 17 convergence.
"""
import sys, math
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True) if hasattr(sys.stdout,'reconfigure') else None
from fractions import Fraction as F

c1 = 0.6973026

def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*abs(e)+1): bps.add(F(m,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set()
        for e in E:
            y=(e*xm)%1; secs.add((y.numerator*7)//y.denominator)
        if len(secs)==7: total+=x1-x0
    return total

def M7(k):
    s=F(0)
    for t in range(0,7):
        s += F((-1)**t * math.comb(6,t)) * F(7-t,7)**(k-1)
    return s

def coord_weight(n):
    """factor for one coordinate value n in [-L,L]: 0-coord contributes weight 1 (and 'missing'
       support); nonzero & 7-mult contributes 0 (kills term); else c1/|n|."""
    if n==0: return 1.0, 0      # weight, support-contribution
    if n%7==0: return 0.0, 0    # kills the term
    return c1/abs(n), 1

def env_partial_mitm(Enz, L, min_support=6):
    """meet in the middle over Enz (the nonzero offsets). Sum prod weights over n with
       sum n_j e_j = 0 and support>=min_support (support = #coords nonzero&not-7-mult)."""
    h = len(Enz)//2
    A = Enz[:h]; B = Enz[h:]
    # left: map partial-sum -> dict over support-count -> accumulated weight
    def build(part):
        table = defaultdict(lambda: defaultdict(float))
        # iterate cartesian product of coord values for this part
        ranges = [range(-L,L+1)]*len(part)
        import itertools
        for ns in itertools.product(*ranges):
            s = sum(n*e for n,e in zip(ns,part))
            w = 1.0; sup = 0; dead=False
            for n in ns:
                cw,sc = coord_weight(n)
                if cw==0.0: dead=True; break
                w*=cw; sup+=sc
            if dead: continue
            table[s][sup]+=w
        return table
    TA = build(A); TB = build(B)
    # match s_A + s_B = 0
    total=0.0; cnt=0
    for sA, dA in TA.items():
        dB = TB.get(-sA)
        if not dB: continue
        for supA,wA in dA.items():
            for supB,wB in dB.items():
                if supA+supB>=min_support:
                    total += wA*wB
    # subtract the all-zero term (n=0): support 0, excluded by min_support>=6 anyway
    return total

if __name__=="__main__":
    print("="*80)
    print("CANDIDATE 2 (fast): support>=6 lattice envelope vs claimed (pi^2/6)^3<17")
    print("="*80)
    tests = [
        ("AP k=8 {0..7}", list(range(8))),
        ("wide k=8 {0,1,2,3,4,5,6,30}", [0,1,2,3,4,5,6,30]),
        ("Sidon k=8 {0,1,3,7,12,20,30,44}", [0,1,3,7,12,20,30,44]),
        ("dissoc k=8 {0,1,3,7,15,31,63,127}", [0,1,3,7,15,31,63,127]),
    ]
    for name,E in tests:
        Enz=[e for e in E if e!=0]
        g=measS7(E); corr=g-M7(len(E))
        print(f"\n  {name}: TRUE signed correction = {float(corr):+.5f}")
        for L in [10,15,20,25,30]:
            ev=env_partial_mitm(Enz,L)
            print(f"     |n|_inf<={L:>3}: support>=6 envelope = {ev:12.4f}")
    print("\n  Claimed convergent bound (pi^2/6)^3 =", round((math.pi**2/6)**3,4))
