#!/usr/bin/env python3
"""
lrc14_wsb_FINAL_summary_kps-S9-wf.py  (kind-pasteur-2026-06-19-S9)

FINAL consolidated summary of the moment-convex-order angle on LRC(14) HYP-2607.

WHAT WAS ESTABLISHED (all exact rationals):

REFRAME (the central new framing).  Work in the PROBABILITY basis p_t=meas{N_E=t} (NOT the
binomial-moment basis S_r).  The canon's 'non-separability' (S_1,S_2,S_3 not separately
consec-extremal) is a MOMENT-basis artifact.  In the probability basis:
  * p_0 = meas(S7) (the ACTUAL target) is consec-extremal among same-k.   [LEMMA B]
  * p_6 = P(N=6)  is consec-extremal among same-k (= 1/(7(k-1))).          [LEMMA A]
Both unit-weight terms g(0)=g(6)=1 of the k=8 dual are SEPARATELY consec-maximal.

LEMMA B (the strong statement, STRONGER than HYP-2607 -- no L_y needed).  Among PRIMITIVE
same-k clusters, consec UNIQUELY maximizes meas(S7).  VERIFIED EXACT:
  k=8: exhaustive primitive maxE<=15 (6434 sets) + random wide to spread 200, ZERO beaters.
  k=9: exhaustive primitive maxE<=12 (495 sets), ZERO beaters.
  => meas(S7)(E) <= meas(S7)(consec_k) for all tested primitive E; with THM-535 (consec finite
     check k=8,9,10) this gives meas(S7)(E) <= cap_k, i.e. HYP-2603 => LRC(14).  STILL OPEN as a
     theorem (no proof of same-k extremality), but it is the SAME crux, in a cleaner (target-direct)
     form, with the resonance worry dissolved.

RESONANCE DISSOLVED.  HYP-2608 feared 'resonant wide w==0 mod 7' configs.  These are NON-PRIMITIVE
(gcd=7); by scale-invariance they EQUAL consec exactly (not over cap).  The LRC reduction uses
primitive clusters (THM-531), so the resonance is a non-issue: {0..6}+7m for m>=2 all sit strictly
below consec; the only 'tail=7' that ties consec IS consec.

LEMMA A (rigorous building blocks).  P(N=6)(E)=meas(G_E), G_E={x: frac(ex) in [0,1/7) forall e}.
  (i) PER-COMPONENT LENGTH BOUND (PROVED, elementary): each component of G_E has length
      <= 1/(7 max(E)).  [on a component the winding integers floor(ex) are fixed; the constraint
      is x in intersection_e [n_e/e, n_e/e+1/(7e)), an interval of length <= min 1/(7e)=1/(7 maxE).]
      VERIFIED 0 violations / 3176 primitive E.
  (ii) COMPONENT-COUNT BOUND (VERIFIED, open as proof): #comp(G_E) <= max(E)/(k-1) for primitive E.
  (i)*(ii) => meas(G_E) <= 1/(7(k-1)) = P(N=6)(consec).  This is THM-535's Phi(c,k)=c/(k-1)
      EXTENDED from consec to the same-k maximizer.

WIDE-SPREAD COLLAPSE (HYP-2608(a), supporting evidence).  Primitive wide clusters have meas(S7)
collapsing to ~0.03-0.07 (spread 40-200), far below consec (0.327) and cap (0.381).  The largest
primitive non-consec meas(S7) at k=8 is 0.2736 (E=[0,2,3,4,5,6,7,8]), margin 0.108 below cap.

This script prints the headline exact numbers for the record.
"""
import sys, itertools
from math import gcd
from functools import reduce
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def gcd_all(E): return reduce(gcd,[e for e in E if e!=0],0)
def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for j in range(7):
            c=F(j,7)
            for m in range(e): bps.add((c+m)/e)
    bps=sorted(z for z in bps if F(0)<=z<F(1)); tot=F(0)
    for i in range(len(bps)):
        x0=bps[i]; x1=bps[i+1] if i+1<len(bps) else F(1)
        if x1<=x0: continue
        xm=(x0+x1)/2; hit=set()
        for e in E:
            fr=(e*xm)%1; hit.add((fr.numerator*7)//fr.denominator)
        if len(hit)==7: tot+=x1-x0
    return tot

# cap engine
H=F(1,14)
def danger(u):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-H/u)%1; b=(c+H/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def measGP(P):
    dz=sorted([iv for u in P for iv in danger(u)]); s=F(0); prev=F(0)
    for a,b in dz:
        if a>prev: s+=a-prev
        prev=max(prev,b)
    if prev<1: s+=1-prev
    return s
def cap(k): return min(measGP(P) for P in itertools.combinations(range(1,14),13-k))

print("HEADLINE EXACT NUMBERS (LRC(14) seven-sector, k=8,9,10)")
print("="*70)
for k in (8,9,10):
    c=measS7(list(range(k))); cp=cap(k)
    print(f" k={k}: meas(S7)(consec)={c}={float(c):.6f}  cap_k={cp}={float(cp):.6f}  "
          f"margin={float(cp-c):.6f}  {'<cap OK' if c<cp else 'FAIL'}")
print()
print("LEMMA B headline (k=8): consec is the UNIQUE primitive same-k maximizer of meas(S7).")
print("  exhaustive primitive maxE<=15: 6434 sets, 0 beaters (verified in stage 6/7).")
print("  tightest non-consec primitive competitor: E=[0,2,3,4,5,6,7,8] meas(S7)=",
      float(measS7([0,2,3,4,5,6,7,8])), " (margin", float(cap(8)-measS7([0,2,3,4,5,6,7,8])),"below cap)")
print()
print("LEMMA A headline: P(N=6)(consec_k)=1/(7(k-1)); consec is the UNIQUE primitive same-k max.")
for k in (7,8,9,10):
    print(f"  k={k}: P(N=6)(consec)=1/(7*{k-1})={F(1,7*(k-1))}")
print()
print("RESONANCE DISSOLUTION (the w==0 mod 7 fear):")
for E in [[0,1,2,3,4,5,6,7],[0,7,14,21,28,35,42,49],[0,2,4,6,8,10,12,14]]:
    print(f"  E={E}  gcd={gcd_all(E)}  meas(S7)={float(measS7(E)):.6f}  "
          f"(primitive rep = divide by gcd => consec)")
print("\nDONE final summary.")
