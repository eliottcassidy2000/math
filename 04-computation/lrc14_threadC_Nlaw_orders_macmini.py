#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD C (mac-mini) -- a FRESH probabilistic angle on the empty-count law N.

measS7(E) = P(N_E = 0), where N_E(x) = # empty sectors of {floor(7 frac(e x)) mod 7 : e in E}.
The occupancy h = 7 - N.  occupancy_pi(E)[h] = P(occupancy = h) EXACTLY (rationals).

CRUX target (sector route, VERIFIED not PROVED): consec_k maximizes measS7 = P(N=0)
for k = 8,9,10 over all offset sets E.

HYP-2722c names the goal as "convex-order dominance of N".  But P(N=0)=E[1[N=0]]
is governed by ONLY the LEFT tail of N.  Pure convex order is two-sided and is
KNOWN to be wrong as a clean ladder (B_3 dirty, k=12 p0 anomaly).  This script
SCANS a *menu* of probabilistic orderings of the law of N to find WHICH one
consec actually wins on the binding rows k=8,9,10, focusing on LEFT-TAIL / Laplace
orders that control P(N=0) directly:

  (A) usual stochastic order  N_consec <=_st N_E   (i.e. P(N>=b) smaller for consec, all b)
      -- this would say consec is "smallest" N.  EXPECTED FALSE (consec has FAT right tail).
  (B) reversed:  P(N>=b) LARGER for consec, all b  (consec stochastically LARGEST)
      -- EXPECTED FALSE on left (consec has SMALL P(N>=1)=1-p0).
  (C) Laplace / MGF order at theta<0:  L(theta)=E[e^{theta N}] for theta<=0.
      Since 1[N=0] = lim_{theta->-inf} e^{theta N}, P(N=0) is the theta->-inf limit
      of L(theta).  HYP: consec MAXIMIZES L(theta)=E[e^{theta N}] for ALL theta<=0
      (uniformly, not just the limit).  This is "consec is smallest in Laplace order".
      THIS is the natural one-sided strengthening of P(N=0)-max.  TEST it exactly.
  (D) entropy of the occupancy law:  HYP: consec MINIMIZES Shannon entropy H(pi)
      (most concentrated occupancy).  TEST.
  (E) E[z^N] (prob gen fn) for z in [0,1]:  P(N=0)=pgf(0).  HYP consec max pgf(z) all z in[0,1].
  (F) the death-chain / decorrelated comparison: gap = true measS7 - decorrelated measS7.

We report, for each order, on how LARGE a box consec wins, and the FIRST violation.
A clean WIN on (C) or (E) for ALL theta/z would be a genuinely NEW, stronger,
one-sided certificate that implies P(N=0)-max as a limit -- the deliverable.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd, exp, log
from functools import reduce

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

def N_law(E):
    """pmf of N=7-occupancy as list indexed by N=0..7 (exact Fractions)."""
    pi=occupancy_pi(E)
    nl=[F(0)]*8
    for h in range(8):
        nl[7-h]+=pi[h]
    return nl

def measS7(nl): return nl[0]               # P(N=0)
def survival(nl):                          # G_b = P(N>=b), b=0..7 (G_0=1)
    G=[F(0)]*8
    acc=F(0)
    for b in range(7,-1,-1):
        acc+=nl[b]; G[b]=acc
    return G

def entropy_float(nl):
    s=0.0
    for p in nl:
        pf=float(p)
        if pf>0: s-= pf*log(pf)
    return s

def mgf_at(nl, theta):
    # E[e^{theta N}] exact requires exp; use high-precision float since we only
    # need a robust comparison.  But to keep it EXACT-comparable we use rationals
    # for theta = log of a rational?  Instead: compare via pgf E[z^N] at rational z
    # (E (z in [0,1])), which IS exact.  mgf at theta<=0 <=> pgf at z=e^theta in (0,1].
    s=0.0
    for n,p in enumerate(nl):
        s+=float(p)*exp(theta*n)
    return s

def pgf_at(nl, z):  # exact rational pgf E[z^N], z rational in [0,1]
    return sum(p*(z**n) for n,p in enumerate(nl))

def primitive(E): return reduce(gcd,[e for e in E if e!=0],0)==1
def consec(k): return list(range(k))

def run(k, W):
    C=consec(k); nlC=N_law(C); GC=survival(nlC); p0C=measS7(nlC)
    bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
    bank=[E for E in bank if primitive(E)]
    # menu of rational z for pgf order (z=e^theta, theta<=0)
    zs=[F(j,10) for j in range(0,11)]            # 0,0.1,...,1.0
    zs+=[F(1,100),F(1,1000),F(99,100)]
    res={
        'p0_beaters':0,
        'st_smallest_fail':0,   # (A) consec smallest in st: G_consec<=G_E all b -> fails if some E has smaller G at some b? we test consec smallest
        'st_largest_fail':0,    # (B) consec largest
        'pgf_fail':0,           # (E) consec max pgf(z) for all z in zs
        'pgf_fail_zlist':set(),
        'entropy_min_fail':0,   # (D) consec min entropy
        'n_shapes':len(bank),
    }
    HC=entropy_float(nlC)
    pgfC={z:pgf_at(nlC,z) for z in zs}
    for E in bank:
        if list(E)==C: continue
        nl=N_law(list(E)); G=survival(nl)
        if measS7(nl)>p0C: res['p0_beaters']+=1
        # (A) consec st-smallest: need G_consec[b] <= G[b] for all b
        if not all(GC[b]<=G[b] for b in range(8)): res['st_smallest_fail']+=1
        # (B) consec st-largest: G_consec[b] >= G[b] all b
        if not all(GC[b]>=G[b] for b in range(8)): res['st_largest_fail']+=1
        # (E) pgf: consec max pgf at all z in [0,1]
        bad=False
        for z in zs:
            if pgf_at(nl,z)>pgfC[z]:
                res['pgf_fail_zlist'].add(float(z)); bad=True
        if bad: res['pgf_fail']+=1
        # (D) entropy: consec min entropy
        if entropy_float(nl)<HC-1e-12: res['entropy_min_fail']+=1
    return C, nlC, p0C, HC, res

if __name__=="__main__":
    print("="*78)
    print("THREAD C: menu of probabilistic orderings of the empty-count law N")
    print("="*78)
    for k,W in [(8,13),(8,15),(9,13),(10,12)]:
        C,nlC,p0C,HC,res=run(k,W)
        print(f"\n--- k={k}, span<={W}, primitive shapes={res['n_shapes']} ---")
        print(f"  consec N-law (P(N=n), n=0..7): {[str(x) for x in nlC]}")
        print(f"  measS7=P(N=0)={float(p0C):.6f}   entropy(occ)={HC:.6f}")
        print(f"  (p0) beaters of consec on P(N=0): {res['p0_beaters']}   "
              f"({'CONSEC MAX' if res['p0_beaters']==0 else 'FAIL'})")
        print(f"  (A) consec stochastically SMALLEST N: fails on {res['st_smallest_fail']} shapes")
        print(f"  (B) consec stochastically LARGEST  N: fails on {res['st_largest_fail']} shapes")
        print(f"  (E) consec MAX pgf E[z^N] all z in [0,1]: fails on {res['pgf_fail']} shapes")
        if res['pgf_fail']:
            print(f"      pgf fails at z-values: {sorted(res['pgf_fail_zlist'])}")
        print(f"  (D) consec MIN entropy(occ): fails on {res['entropy_min_fail']} shapes")
    print("\nDONE.")
