#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE 5 (Markoff/CF, sharpened via THREE-GAP) -- opus 2026-06-21.

The prompt's deepest claim: AP = simplest continued fraction; the three-gap
(Steinhaus) theorem is the CF/Markoff bridge. For consec, the orbit at x is the
AP {0,x,...,(k-1)x} whose THREE gaps g1>=g2>=g3=g1-g2 (Steinhaus: 3rd gap is the
difference) are determined by the CF of x. measS7(consec) = meas{x : maxgap<=1/7}.

MARKOFF/THREE-GAP HYPOTHESIS (falsifiable):
  The three-gap lengths (L,S,L+S form -- a Stern-Brocot/Fibonacci triple) make the
  AP orbit the MOST UNIFORM k-point set reachable by a single frequency. Define the
  three-gap "Markoff balance"  B(x) = L*S*(L+S) ... no. The right invariant from CF:
  the orbit is a 1/7-net iff maxgap L <= 1/7. The measure of net x is governed by
  the FAREY/STERN-BROCOT intervals where L <= 1/7.

TEST 1: For consec_k, decompose measS7 as a sum over FAREY fractions p/q (q<=...)
  of local net-widths, and check each local width equals a THREE-GAP closed form
  (the "g-leg" of the verified G_P).
TEST 2 (the Markoff extremality): for a GENERAL E, the orbit is NOT an AP; it has
  up to k gaps. The net region {maxgap<=1/7} shrinks when the orbit is LESS uniform.
  Falsifiable: measS7(E) <= measS7(consec) because the AP achieves the three-gap
  (CF) optimum -- i.e. consec MINIMIZES the maximum gap pointwise? TEST whether
  for ALL x, maxgap_consec(x) <= min over residue-matching E of maxgap_E(x).
  (If TRUE pointwise, the wall is a pointwise three-gap optimum -- a clean proof!)
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def frac(x): return x-(x.numerator//x.denominator)
def consec(k): return list(range(k))
def is_full_residue(E): return frozenset(e%7 for e in E)==frozenset(range(7))
def primitive(E): return reduce(gcd,[abs(e) for e in E if e!=0],0)==1

def maxgap(E, x):
    pts=sorted(set(frac(F(e)*x) for e in E))
    if len(pts)==1: return F(1)
    g=F(0)
    for i in range(len(pts)):
        nxt=pts[(i+1)%len(pts)] if i<len(pts)-1 else pts[0]+1
        g=max(g, nxt-pts[i])
    return g

def measS7(E):
    E=sorted(set(int(e) for e in E)); bps={F(0),F(1)}
    for e in E:
        ae=abs(e)
        if ae==0: continue
        for m in range(7*ae+1): bps.add(F(m,7*ae))
    bps=sorted(b for b in bps if 0<=b<=1); tot=F(0)
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        if maxgap(E,(lo+hi)/2)<=F(1,7): tot+=hi-lo
    return tot

if __name__=="__main__":
    print("#"*78); print("# ANGLE 5: THREE-GAP / CF MARKOFF EXTREMALITY (POINTWISE TEST)"); print("#"*78)
    k=8; C=consec(k); mC=measS7(C)
    print(f"\nconsec(8): measS7 = {mC} = {float(mC):.6f}")

    # TEST 2: POINTWISE maxgap comparison. For each full-residue primitive E,
    # is maxgap_consec(x) <= maxgap_E(x) for ALL x? If yes pointwise, then the
    # net{maxgap<=1/7} of consec CONTAINS that of E => measS7 max, a CLEAN proof.
    print("\nTEST 2 (POINTWISE three-gap optimum): does consec minimize maxgap(x) for all x?")
    W=13
    bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
    full=[list(E) for E in bank if primitive(E) and is_full_residue(E)]
    print(f"  testing {len(full)} full-res shapes against consec, pointwise on fine grid...")
    # build common breakpoint grid from all E plus consec
    allden=set([7])
    for E in [C]+full[:0]:  # consec only for grid base; sample midpoints densely
        pass
    # sample at midpoints of consec's own breakpoints AND a dense rational grid
    grid=set()
    for e in C:
        ae=abs(e)
        if ae==0: continue
        for m in range(7*ae+1): grid.add(F(m,7*ae))
    # add finer grid q up to 20
    for q in range(1,21):
        for p in range(q+1): grid.add(F(p,q))
    grid=sorted(grid); mids=[(grid[i]+grid[i+1])/2 for i in range(len(grid)-1)]

    pointwise_dominators=0; violators=[]
    for E in full:
        dom=True
        for x in mids:
            if maxgap(C,x) > maxgap(E,x):
                dom=False; break
        if dom: pointwise_dominators+=1
        else: violators.append(E)
    print(f"  shapes where consec maxgap(x) <= E maxgap(x) for all sampled x: {pointwise_dominators}/{len(full)}")
    print(f"  => POINTWISE optimum holds: {pointwise_dominators==len(full)}")
    if violators:
        print(f"  REFUTED pointwise. Example violator: {violators[0]}")
        # show an x where consec has a LARGER maxgap
        E=violators[0]
        for x in mids:
            if maxgap(C,x)>maxgap(E,x):
                print(f"    at x={x}: maxgap_consec={maxgap(C,x)}={float(maxgap(C,x)):.4f} > "
                      f"maxgap_E={maxgap(E,x)}={float(maxgap(E,x)):.4f}")
                break
        # but does consec still win the MEASURE despite losing pointwise?
        print(f"    measS7(consec)={float(mC):.6f}  measS7(violator)={float(measS7(E)):.6f}  consec still bigger: {measS7(E)<mC}")
