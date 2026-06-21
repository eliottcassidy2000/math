#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A — the CLEAN B_2 rung, fast & exact (no relation-box enumeration).

B_2 = 1 - S_1 + S_2,  S_1=E[N], S_2=E[C(N,2)].
Equivalent: B_2 = E[(N-1)(N-2)/2] = E[g_2(N)], g_2 convex.

EXACT decomposition into single- and pair-sector miss measures (occupancy law,
fast):
   S_1 = sum_{s=1}^{6} P(s missed)               (sector 0 never missed)
   S_2 = sum_{1<=s<t<=6} P(s,t both missed)
       + sum_{t=1}^{6} P(0,t both missed)         (but 0 never missed -> =0)
   so all miss events live on inner sectors {1..6}.

THE SHARP REFORMULATION (this is the proof-bearing identity):
   B_2 = 1 - sum_s P_s + sum_{s<t} P_{st}
   where P_s = P(sector s empty), P_{st}=P(sectors s,t both empty).
   By inclusion-exclusion truncation, B_2 = E[ 1 - N + C(N,2) ] and since
   g_2(N)=(N-1)(N-2)/2 vanishes at N=1,2 and is convex,
   B_2 = P(N=0) + sum_{N>=3} C(N-1,2)... -> measS7 + (deg>=3 convex tail).
   B_2 is an UPPER BOUND on measS7 (Bonferroni, even truncation), and the GAP
   B_2 - measS7 = sum_{N>=3} [g_2(N)-[N=0]] pi = sum_{N>=3} g_2(N) pi(N) >= 0.

CLAIM (clean rung): consec maximizes B_2.  Test fast & exact over LARGE boxes,
and report the exact margin to the runner-up (the certificate).
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
from collections import defaultdict

def occupancy_pi(E):
    E=sorted(set(E))
    bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*abs(e)+1): bps.add(F(a,7*abs(e)))
    bps=sorted(b for b in bps if 0<=b<=1)
    pi=[F(0)]*8
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        xm=(lo+hi)/2; hit=set()
        for e in E:
            v=e*xm; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        pi[len(hit)]+=hi-lo
    return pi

def B2_from_pi(pi):
    # B_2 = sum_h pi[h] * g2(7-h), g2(N)=(N-1)(N-2)/2
    return sum(pi[h]*F((7-h-1)*(7-h-2),2) for h in range(8))
def measS7_from_pi(pi): return pi[7]
def primitive(E): return reduce(gcd,[e for e in E if e!=0],0)==1
def consec(k): return list(range(k))

if __name__=="__main__":
    print("="*78)
    print("CLEAN RUNG: consec maximizes B_2 = E[(N-1)(N-2)/2]  (fast exact)")
    print("="*78)
    for k,W in [(8,13),(8,15),(9,13),(10,12)]:
        C=consec(k); piC=occupancy_pi(C); B2c=B2_from_pi(piC)
        bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
        bank=[E for E in bank if primitive(E)]
        beat=0; second=F(-1); secondE=None
        for E in bank:
            if list(E)==C: continue
            b2=B2_from_pi(occupancy_pi(list(E)))
            if b2>B2c+F(1,10**18): beat+=1
            if b2>second: second=b2; secondE=list(E)
        margin=B2c-second
        print(f" k={k} span<= {W}: {len(bank)} shapes; B_2(consec)={float(B2c):.7f}")
        print(f"    beaters={beat}  ({'CONSEC MAX' if beat==0 else 'FAIL'})  runner-up={secondE} B_2={float(second):.7f}")
        print(f"    EXACT margin B_2(consec)-B_2(runnerup) = {margin} = {float(margin):.2e}")
