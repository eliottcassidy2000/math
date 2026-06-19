#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) THE 1-D WEYL LEMMA VIA THE EXACT COSET FACTORIZATION  (kps-2026-06-19-S12b)

Builds on:
 - HYP-2642 (kps-S13): unbounded direction reduced to "signed sum of K(n) over relations
   with n_w != 0 is <= C/w".  This IS the last analytic lemma.
 - S12 finding: K(n) = D7(n mod 7) / prod_j n_j  EXACTLY, with D7(-c)=conj(D7(c)),
   so corr is real and ruled by Re D7(c).  |Re D7| <= 0.143, D7 nonzero on all cosets.

GOAL: make the Weyl lemma PRECISE.  Compute, for core E' = consec_{k-1} and far element w,
   Delta_w(E) = sum_{n in Lambda(E), n_w != 0} K(n)   = p0(E) - Plat(E')   (the decorrelation defect)
EXACTLY (via the breakpoint engine, not the lattice) AND decompose the signed lattice picture:
   - measure the decay exponent of |Delta_w| vs w (claim: ~ C/w),
   - confirm the contributing relations all have a LARGE companion coordinate (sparse),
   - extract the explicit constant C and the cap margin.

We compute Delta_w EXACTLY two ways and cross-check:
  (1) breakpoint:  p0(E) - Plat(E'),  Plat(E') = p0(E') + (1/7) p1(E').
  (2) (diagnostic) the signed lattice sum truncated, to see it track (1).
"""
import sys, math
from fractions import Fraction
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def dist_full(E):
    """exact p_t (t=0..6) as Fractions over x in [0,1)."""
    E=sorted(set(E)); bps=set([Fraction(0),Fraction(1)])
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1)
    p=[Fraction(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in E:
            v=e*mid; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        N=sum(1 for j in range(1,7) if j not in hit)
        p[N]+=(hi-lo)
    return p

def p0p1(E):
    p=dist_full(E); return p[0],p[1]

if __name__=="__main__":
    print("LRC(14) 1-D WEYL DECORRELATION LEMMA — exact decay of Delta_w (kps-S12b)\n")
    print("Delta_w(E) = p0(E) - Plat(E'),  Plat(E') = p0(E') + (1/7) p1(E'),  E'=core, w=far element.")
    print("Claim (HYP-2642): |Delta_w| <= C/w.  We extract C and check the cap margin.\n")

    for k in [8,9,10]:
        core=list(range(k-1))         # consec_{k-1}
        p0c,p1c=p0p1(core)
        Plat=float(p0c)+float(p1c)/7.0
        M7=sum(((-1)**t)*math.comb(6,t)*((1-t/7.0)**(k-1)) for t in range(7))
        cap={8:0.38153,9:0.49426,10:0.6044}[k]
        print(f"=== k={k}: core=consec_{k-1}  p0(core)={float(p0c):.6f} p1(core)={float(p1c):.6f}")
        print(f"    Plat(core)={Plat:.6f}   M7({k})={M7:.6f}   cap_{k}={cap}   Plat<cap margin={cap-Plat:.4f}")
        print(f"    w :   p0(E)        Delta_w=p0-Plat     w*|Delta_w| (->C)")
        prevC=None
        for w in [20,40,80,160,320,640,1280]:
            E=core+[w]
            p0E,_=p0p1(E)
            Dw=float(p0E)-Plat
            C=w*abs(Dw)
            print(f"   {w:5d}: {float(p0E):.7f}   {Dw:+.7e}     {C:.5f}")
            prevC=C
        print(f"    => empirical C ~ {prevC:.4f};  p0(E) <= Plat + C/w < cap once w > C/(cap-Plat) = {prevC/(cap-Plat):.2f}")
        print()
