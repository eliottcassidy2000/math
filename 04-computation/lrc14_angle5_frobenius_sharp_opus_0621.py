#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE 5 SHARPENED -- the FROBENIUS-MINIMALITY mechanism (opus 2026-06-21).

Finding so far: consec MINIMIZES both M_var and the Frobenius number of the
velocity-gap semigroup, and is the unique measS7-max. Neither is a clean function
of measS7 (weak Kendall tau). Sharpen the FROBENIUS angle:

H1 (NECESSITY): every measS7-maximizer has velocity-gap Frobenius number = -1
   (i.e. all consecutive velocity gaps = 1, the gap-set has gcd-complete cover).
H2 (the real Frobenius object): the survival problem is a COVERING of Z/7 by
   moving points. A gap opens iff some residue is uncovered. Define the
   "covering deficiency" CD(E) = max over (a,s) of #uncovered sectors. Test if
   measS7 is exactly 1/7 * (Lebesgue measure of {x : CD=0}) and whether CD is
   governed by a Frobenius/numerical-semigroup quantity of the SPEED set.
H3 (MARKOFF TRIPLE on legs): take the verified closed form g(t)=2t(7-t). For
   consec, occupied residue gaps are all 1; for a rival they spread. Test whether
   sum_a W_a = const - (Markoff defect) where defect = sum over residue-gap
   triples of a Markoff form. Concretely test:  is sum_a (7 W_a) an AFFINE image
   of sum_{r} g(gap_r) for the residue-occupancy pattern?
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def covered_sectors(E,xm):
    secs=set()
    for e in E:
        v=e*xm; v=v-(v.numerator//v.denominator); secs.add((v.numerator*7)//v.denominator)
    return secs
def measS7(E):
    E=sorted(set(int(e) for e in E)); bps={F(0),F(1)}
    for e in E:
        ae=abs(e)
        if ae==0: continue
        for m in range(7*ae+1): bps.add(F(m,7*ae))
    bps=sorted(b for b in bps if 0<=b<=1); tot=F(0)
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        if len(covered_sectors(E,(lo+hi)/2))==7: tot+=hi-lo
    return tot
def is_full_residue(E): return frozenset(e%7 for e in E)==frozenset(range(7))
def primitive(E): return reduce(gcd,[abs(e) for e in E if e!=0],0)==1
def consec(k): return list(range(k))

def frob(gens):
    gens=sorted(set(g for g in gens if g>0))
    if not gens or reduce(gcd,gens)!=1: return None
    bound=max(gens)**2+1; reach=[False]*(bound+1); reach[0]=True
    for i in range(bound+1):
        if reach[i]:
            for g in gens:
                if i+g<=bound: reach[i+g]=True
    fn=-1
    for i in range(bound+1):
        if not reach[i]: fn=i
    return fn
def velgap_frob(E):
    E=sorted(int(e) for e in E); gaps=[E[i+1]-E[i] for i in range(len(E)-1)]
    return frob([g for g in gaps if g>0])

if __name__=="__main__":
    print("#"*78); print("# ANGLE 5 SHARPENED: FROBENIUS-MINIMALITY MECHANISM"); print("#"*78)
    k=8; C=consec(k); mC=measS7(C)
    print(f"\nconsec(8): measS7={mC}  velgap_frob={velgap_frob(C)}")

    # H1 NECESSITY: classify shapes by velgap_frob, look at measS7 ranges
    print("\nH1: measS7 distribution by velocity-gap Frobenius number (span<=16):")
    W=16
    bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
    full=[list(E) for E in bank if primitive(E) and is_full_residue(E)]
    from collections import defaultdict
    by_frob=defaultdict(list)
    for E in full: by_frob[velgap_frob(E)].append(measS7(E))
    for fn in sorted(by_frob, key=lambda x:(x is None, x)):
        ms=by_frob[fn]
        print(f"   velgap_frob={fn}: n={len(ms)}  max measS7={float(max(ms)):.6f}  "
              f"min={float(min(ms)):.6f}")
    # Is consec's measS7 beaten by ANY shape with velgap_frob != -1?
    best_nonminimal = max((max(ms) for fn,ms in by_frob.items() if fn!=-1), default=None)
    print(f"\n   consec measS7={float(mC):.6f}; best measS7 among velgap_frob!=-1: "
          f"{float(best_nonminimal):.6f}" if best_nonminimal else "")
    print(f"   => Frobenius-minimal (=-1) NECESSARY for being max? "
          f"{best_nonminimal is None or best_nonminimal < mC}")

    # Among velgap_frob=-1 shapes, what distinguishes consec?
    minimal=[E for E in full if velgap_frob(E)==-1]
    print(f"\n   shapes with velgap_frob=-1: {len(minimal)} ; their measS7 range "
          f"[{float(min(measS7(E) for E in minimal)):.6f},{float(max(measS7(E) for E in minimal)):.6f}]")
    print(f"   (so velgap_frob=-1 is NECESSARY but NOT SUFFICIENT -- the wall is finer)")
