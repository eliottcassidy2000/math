#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A — HONEST CORRECTION: consec does NOT maximize the UNSIGNED relation
weight-enumerator (A_2 short-count maximizer is [0,1,2,3,4,6,8,12], not consec).
So the MacWilliams link is NOT 'consec = most short relations'.  The correct
statement must be SIGNED: corr = sum_n K(n) with K(n) SIGN-VARYING.

This script pins down the TRUE structure by computing, for consec vs the A_2
champion vs a Sidon (MDS) set:
  - measS7, corr, B_2, B_4
  - the EXACT depth law pi and the even-Krawtchouk bands
  - sector-pair co-miss EXACT (the S_2 = E[C(N,2)] association term)
and shows WHICH object consec extremizes: the CONVEX g_2 read of the depth law,
NOT the relation count.  The MacWilliams duality is then:
  measS7 = even-Krawtchouk (MacWilliams) read of the OCCUPANCY weight law pi_E,
  and the relation code is the DUAL coordinate (Poisson summation / Weyl), but the
  EXTREMALITY is carried by the CONVEX (even-band) read, which is a SIGNED
  combination of relation shells, NOT the unsigned enumerator.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
from collections import defaultdict

def occupancy_pi(E):
    E=sorted(set(E)); bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*abs(e)+1): bps.add(F(a,7*abs(e)))
    bps=sorted(b for b in bps if 0<=b<=1); pi=[F(0)]*8
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        xm=(lo+hi)/2; hit=set()
        for e in E:
            v=e*xm; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        pi[len(hit)]+=hi-lo
    return pi
def gJ(h,J): return sum((-1)**j*comb(7-h,j) for j in range(J+1))
def bands(E):
    pi=occupancy_pi(E)
    return {J:sum(pi[h]*gJ(h,J) for h in range(8)) for J in range(7)},pi
def iid_measS7(k):
    Snk=sum((-1)**(7-i)*comb(7,i)*i**k for i in range(8)); return F(Snk,7**k)
def primitive(E): return reduce(gcd,[e for e in E if e!=0],0)==1

if __name__=="__main__":
    k=8
    shapes={
        "consec [0..7]":list(range(8)),
        "A_2-champ [0,1,2,3,4,6,8,12]":[0,1,2,3,4,6,8,12],
        "single-shift [0,2,3,4,5,6,7,8]":[0,2,3,4,5,6,7,8],
        "Sidon/MDS [0,1,3,7,12,20,30,44]":[0,1,3,7,12,20,30,44],
    }
    print("="*92)
    print(f"{'shape':38s} {'measS7':>9s} {'corr':>9s} {'B_1':>8s} {'B_2':>8s} {'B_3':>8s} {'B_4':>8s} {'B_5':>8s}")
    print("="*92)
    ik=iid_measS7(k)
    rows={}
    for nm,E in shapes.items():
        b,pi=bands(E); m=pi[7]; rows[nm]=(b,pi,m)
        print(f"{nm:38s} {float(m):9.5f} {float(m-ik):+9.5f} {float(b[1]):8.4f} {float(b[2]):8.4f} {float(b[3]):8.4f} {float(b[4]):8.4f} {float(b[5]):8.4f}")
    print(f"\n iid_measS7(k=8) = {float(ik):.6f}")
    print("\n DEPTH LAW pi[h] (h=0..7), and N=7-h occupancy:")
    for nm,E in shapes.items():
        _,pi,_=rows[nm]
        print(f"   {nm:38s} pi={[round(float(x),4) for x in pi]}")
    print("\n OBSERVATION: among these, the measS7 (=B_6) ORDER and the B_2,B_4 ORDER")
    print(" should AGREE (even bands consec-max) but B_3 (odd) may not.  And the relation")
    print(" A_2-champion does NOT win measS7 -> unsigned relation count is the wrong enumerator.")
    # Rank by each band
    print("\n RANKING (descending) by each EVEN band vs measS7:")
    for J in [2,4,6]:
        order=sorted(rows.items(),key=lambda kv:-(kv[1][0][J] if J<6 else kv[1][2]))
        print(f"   B_{J}: "+" > ".join(nm.split()[0] for nm,_ in order))
    order=sorted(rows.items(),key=lambda kv:-(kv[1][0][3]))
    print(f"   B_3(odd): "+" > ".join(nm.split()[0] for nm,_ in order))
