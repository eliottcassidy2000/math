#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The POCKET STRUCTURE of L_y via SUMSET EXCESS / Freiman dimension. kind-pasteur-2026-06-19-S12.

HYP-2607 crux: consec maximizes L_y(E)=E[g(N_E)] (N=#missed of 6 sectors at 1/7-resolution).
NEW FRAME (from S560 summand/multiplicand geometry + Freiman-Vosper):
  excess(E) = |E+E| - (2k-1) >= 0,  excess=0  <=>  E is an AP (dilate of consec)  [Freiman-Vosper].
  All APs have the SAME L_y (= L_y(consec)) by scale-invariance (THM-531).
HYPOTHESIS: L_y(E) is controlled by excess; the "pockets" are indexed by Freiman dimension d:
  d=1 (excess 0)   = AP                 = pocket 1 (the extremal; finite/exact)
  d=2 (small exc)  = 2-dim GAP          = pocket 3
  d>=3             = higher GAP         = pocket 4+
  dissociated (excess ~ C(k,2)) = pocket 2 (-> independent limit L_y^inf)
GOAL: is L_y(E) <= L_y(consec) with the deficit increasing in excess? Look for pocket 4.
"""
import sys, itertools, random
from fractions import Fraction
from functools import reduce
from math import gcd, comb
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def N_at(E,x):
    hit=set()
    for e in E:
        v=e*x; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
    return sum(1 for j in range(1,7) if j not in hit)
def dist_p(E):
    E=sorted(set(E)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); p=[Fraction(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        p[N_at(E,(lo+hi)/2)]+=hi-lo
    return p
def g_poly(k,t):
    if k==8: return Fraction((t-1)*(t-2)*(t-4)*(t-5),40)
    if k in(9,10): return Fraction(-(t-2)*(t-3)*(t-6),36)
    return Fraction((t-3)*(t-4),12)
def L_y(E,k):
    p=dist_p(E); return sum(p[t]*g_poly(k,t) for t in range(7))
def excess(E):
    E=set(E); ss={a+b for a in E for b in E}
    return len(ss)-(2*len(E)-1)
def primitive(E): return reduce(gcd,E)==1

def gap2(base_len, step, k):
    """a 2-dim GAP: {i + step*j} arranged to give k elements. {0..a-1} + step*{0..b-1}."""
    a=base_len; b=(k+a-1)//a
    S=sorted({i+step*j for i in range(a) for j in range(b)})[:k]
    return S

if __name__ == "__main__":
    print("=== L_y vs sumset excess: is L_y controlled by additive structure? ===\n")
    for k in [8,9]:
        C=list(range(k)); Lc=L_y(C,k); exc_c=excess(C)
        capf={8:0.38153,9:0.49426}[k]
        print(f"--- k={k}: consec L_y={float(Lc):.5f} excess={exc_c} (must be 0=AP) cap={capf} ---")
        rows=[]
        # AP dilate (excess 0)
        for d in [2,3,5]:
            E=[d*i for i in range(k)]
            if primitive(E): continue  # gcd=d, skip (not primitive) -- but check L_y invariance
            rows.append((f"dilate x{d}", E, L_y(E,k), excess(E)))
        # explicit AP-dilate primitive? APs are {a*i}; primitive only if a=1. So all primitive APs = consec.
        # 2-dim GAPs (small excess)
        for (a,st) in [(2,7),(2,5),(3,7),(2,9),(4,7),(3,11),(2,11)]:
            E=gap2(a,st,k)
            if len(E)==k and primitive(E):
                rows.append((f"GAP2 base{a} step{st}", tuple(E), L_y(E,k), excess(E)))
        # the known third-pocket examples
        for E in [[0,3,5,16,28,30,33,35],[0,4,12,15,20,21,25,31]] if k==8 else []:
            if len(set(E))==k: rows.append(("3rd-pocket ex", tuple(E), L_y(E,k), excess(E)))
        # random dissociated (large excess)
        rng=random.Random(7+k)
        for _ in range(6):
            E=sorted(rng.sample(range(1,5*k),k-1)); E=[0]+E
            if primitive(E): rows.append(("random", tuple(E), L_y(E,k), excess(E)))
        rows.sort(key=lambda r:-float(r[2]))
        for name,E,L,ex in rows:
            flag=" >consec!" if L>Lc else ""
            print(f"   L_y={float(L):.5f} exc={ex:3d}  {name:22s} {E}{flag}")
        print()
    print("=== TEST: over ALL primitive k=8 sets (spread<=16), is L_y monotone-bounded by excess? ===")
    k=8; C=list(range(k)); Lc=L_y(C,k)
    by_exc={}  # excess -> max L_y
    cnt=0
    for tail in itertools.combinations(range(1,17),k-1):
        E=(0,)+tail
        if not primitive(E): continue
        cnt+=1
        ex=excess(E); L=L_y(E,k)
        if ex not in by_exc or L>by_exc[ex][0]: by_exc[ex]=(L,E)
    print(f"   scanned {cnt} primitive sets. max L_y by excess (excess: maxL_y  witness):")
    for ex in sorted(by_exc):
        L,E=by_exc[ex]
        print(f"     excess={ex:3d}:  maxL_y={float(L):.5f}  {E}")
    print(f"\n   consec L_y={float(Lc):.5f} at excess=0. If maxL_y STRICTLY DECREASES with excess => excess is the right order.")
