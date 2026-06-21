#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE 5 -- the EXACT meas(G_{p,q}) closed form & its Markoff/CF extremality.
opus 2026-06-21.

DISCOVERY: meas(G_{1,q}) = 5/7 + 1/(7q), monotone DEC in q, minimized at q=13.
The binding pair {1,13} has the simplest CF [13] (single partial quotient) and the
largest q. This is the CF/Markoff extremal pair: simplest CF + extreme denominator.

Now find the GENERAL closed form meas(G_{p,q}) and verify the Markoff structure:
the binding (minimal-meas) pair over ALL {p,q} is {1,13}, and we conjecture it is
the unique pair achieving the Frobenius/CF extreme. Test:
  C1: meas(G_{p,q}) closed form via the THREE legs (p, q, and the overlap pq-resonance).
  C2: the minimal-meas pair is exactly the {1, q_max} family endpoint -- i.e. the
      'most linear' pair (gcd 1, smallest p, largest q): a degenerate AP/CF.
  C3: connect to Frobenius: meas deficit (6/7 - meas) = sum of overlap arc lengths;
      for {1,q} it is (q-1)/(7q) = (1/7)(1 - 1/q), the Frobenius 'density' of the
      semigroup <1,q> hitting the 1/14-window. (Frobenius number of <1,q> is -1 but
      the WEIGHTED overlap density is (q-1)/q -- the 'genus fraction'.)
"""
import sys
from fractions import Fraction as F
from itertools import combinations
from math import gcd
sys.stdout.reconfigure(line_buffering=True)

def measGP(P):
    bad=[]
    for p in P:
        p=int(p)
        for kk in range(0,p+1):
            lo=max((F(kk)-F(1,14))/p,F(0)); hi=min((F(kk)+F(1,14))/p,F(1))
            if lo<hi: bad.append((lo,hi))
    bad.sort(); tot=F(0); cur=cb=None
    for lo,hi in bad:
        if cur is None: cur,cb=lo,hi
        elif lo<=cb: cb=max(cb,hi)
        else: tot+=cb-cur; cur,cb=lo,hi
    if cur is not None: tot+=cb-cur
    return F(1)-tot

if __name__=="__main__":
    print("#"*78); print("# ANGLE 5: EXACT meas(G_{p,q}) & MARKOFF/CF EXTREMALITY"); print("#"*78)

    # C1: verify the {1,q} closed form 5/7 + 1/(7q)
    print("\nC1: meas(G_{1,q}) == 5/7 + 1/(7q) ?")
    allok=True
    for q in range(2,14):
        pred=F(5,7)+F(1,7*q); act=measGP([1,q])
        ok=pred==act; allok&=ok
        print(f"   q={q}: pred={pred} act={act} match={ok}")
    print(f"   ALL match: {allok}")

    # general meas(G_{p,q}). For gcd(p,q)=1, the two BAD-arc systems (around k/p and
    # j/q) overlap. Without overlap, meas(G) = 1 - p*(1/(7p)) - q*(1/(7q)) + overlaps
    # = 1 - 1/7 - 1/7 + overlap = 5/7 + overlap_correction.
    # Each p has p arcs of width 1/(7p) (total 1/7); plus endpoints. Overlap = where
    # ||px||<1/14 AND ||qx||<1/14 simultaneously.
    print("\nC1b: meas(G_{p,q}) = 5/7 + (overlap of the two bad-arc systems). Check deficit:")
    for p,q in [(1,13),(1,12),(2,13),(3,13),(5,8),(2,11),(3,11)]:
        m=measGP([p,q]); overlap=m-F(5,7)
        print(f"   {{{p},{q}}}: meas={m}={float(m):.5f}  overlap=meas-5/7={overlap}={float(overlap):.5f}")

    # C2: the global pair-min. Confirm {1,13}. Then test: among pairs with gcd 1,
    # is meas monotone increasing in p (for fixed q) and dec in q (for fixed p)?
    print("\nC2: monotonicity -- meas(G_{p,q}) vs p (fixed q=13) and vs q (fixed p=1):")
    print("   fixed q=13, vary p:")
    for p in range(1,13):
        if gcd(p,13)==1: print(f"     p={p}: meas={float(measGP([p,13])):.5f}")
    # C3: Frobenius/CF reading
    print("\nC3: FROBENIUS/CF reading of the {1,q} family.")
    print("   meas(G_{1,q}) = 5/7 + 1/(7q). deficit from 6/7 = (q-1)/(7q) = (1/7)(1-1/q).")
    print("   (1-1/q) = density of <1,q>-overlap; q->inf gives the pure 1/7 deficit (lonely-runner floor).")
    print("   The binding q=13 is the LARGEST q in P={1..13}: the extremal CF endpoint.")
    print("   => binding G_P = the single-partial-quotient (simplest CF) pair at max denominator.")
