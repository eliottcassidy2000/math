#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROUTE 2 -- the e=0 PERFECT SHORT and the per-side structural bound.
Two exact facts to leverage into a clean partial proof:

FACT 1 (PERFECT SHORT, proved): e=0 is in EVERY admissible shape (we anchor 0). Clock
e=0 is stationary, sector 0 covered for ALL y. So sector 0 NEVER binds. This removes one
of the 7 constraints -> only 6 sectors can bind.

FACT 2 (the binding speed is bounded below by the MINIMUM nonzero |e|... no, ABOVE).
Going right, sector s (s != native-0) is covered by its native clock until 1/(7 e_nat).
The window T+ <= min over non-0 native sectors of (their refill-extended death). A clean
UPPER bound: T+ <= 1/(7 * m1) is FALSE in general (refill). But consec has the property
that its 6 nonzero residues are EXACTLY {1,2,3,4,5,6} (smallest possible), maximizing each
native survival 1/(7r). 

TEST a clean SEPARABLE upper bound:
  U(E) = sum_{a=1..6} [ 1/(7*minmag_right(a)) + 1/(7*minmag_left(a)) ]
  where minmag_right(a) = the magnitude whose native sector binds first IGNORING refill
  = max over residues r!=0 of |e_r|... = the LARGEST min-magnitude among nonzero residues.
For consec that's 6 (residues 1..6, max=6). For any full-res shape, the nonzero residues
1..6 each have min-magnitude >= the residue value? NO -- residue r can have min-mag = r,
r+7, ... but also a DIFFERENT residue. The 6 nonzero residues' min-mags are 6 DISTINCT
positive integers, one per residue class 1..6, EACH >= its residue OR could be smaller if
... no: min-mag in class r is == r mod 7, smallest is r itself. So min-mag(r) in {r,r+7,..}
hence min-mag(r) >= r. So the multiset of 6 nonzero min-mags is {>=1,>=2,...} ... actually
each residue independently >= its own value, but as a SET they need not be {1..6}.
The product/sum of 1/(7 minmag) is MAXIMIZED when each min-mag is smallest = consec.

Let's test if this SEPARABLE bound U(E) (sum of 1/(7*minmag) over the binding sectors)
is (a) an upper bound for WIN, and (b) maximized by consec, and (c) TIGHT at consec.
"""
import sys, itertools
from fractions import Fraction as F
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)
HALF=F(1,14)
def sector_of(e,a,y):
    pos=F(e*a)+F(7*e)*y; return (pos.numerator//pos.denominator)%7
def covered(E,a,y): return len({sector_of(e,a,y) for e in E})==7
def breakpoints(E,a):
    bps={F(0),-HALF,HALF}
    for e in E:
        if e==0: continue
        lo_val=F(7*e)*(-HALF)+F(e*a); hi_val=F(7*e)*(HALF)+F(e*a)
        lo_i=min(lo_val,hi_val);hi_i=max(lo_val,hi_val); m=lo_i.numerator//lo_i.denominator
        while m<=hi_i.numerator//hi_i.denominator+1:
            y=F(m-e*a,7*e)
            if -HALF<=y<=HALF: bps.add(y)
            m+=1
    return sorted(bps)
def window_TpTm(E,a):
    bps=breakpoints(E,a); ivals=list(zip(bps,bps[1:])); Tp=F(0)
    for lo,hi in ivals:
        if hi<=0: continue
        lo2=max(lo,F(0))
        if covered(E,a,(lo2+hi)/2): Tp=hi
        else:
            if lo2==F(0): Tp=F(0)
            break
    Tp=min(Tp,HALF); Tm=F(0)
    for lo,hi in reversed(ivals):
        if lo>=0: continue
        hi2=min(hi,F(0))
        if covered(E,a,(lo+hi2)/2): Tm=-lo
        else:
            if hi2==F(0): Tm=F(0)
            break
    Tm=min(Tm,HALF); return Tp,Tm
def WIN(E): return sum(sum(window_TpTm(E,a)) for a in range(1,7))
def is_full_residue(E): return frozenset(e%7 for e in E)==frozenset(range(7))
def consec(k): return list(range(k))

def minmag_by_residue(E):
    d={}
    for e in E:
        r=e%7
        if r not in d or abs(e)<d[r]: d[r]=abs(e)
    return d  # r->min|e|

def Ubound_simple(E):
    """U = sum_{a=1..6} [ 1/(7*maxminmag) + 1/(7*maxminmag) ] -- crude per-resonance bound:
       T+ <= 1/(7 * (largest nonzero-residue min-mag))? Let's instead try: per-resonance,
       T+ <= 1/(7 * b_min(a)) where b_min(a)= the min over nonzero sectors of the native
       survival rate = ... we compute the ACTUAL native-only first death (no refill)."""
    mm=minmag_by_residue(E)  # r->min mag
    # native-only T+: ignore refill. sector s native occupant residue r=s*a^{-1}. survives
    # to 1/(7*mm[r]). sector 0... residue 0 has mm[0]=0 (the short) -> infinite. 
    # native-only death = min over r!=0 of 1/(7 mm[r]) = 1/(7 * max_{r!=0} mm[r]).
    tot=F(0)
    maxmm=max(mm[r] for r in range(1,7))  # depends only on shape, not a!
    nat = F(1,7*maxmm)
    return 12*nat  # 6 a's, 2 sides, all same native bound (no refill)

for k,span in [(8,14),(9,16),(10,18)]:
    C=consec(k)
    bank=[c for c in itertools.combinations(range(0,span+1),k) if 0 in c and is_full_residue(c)]
    Cwin=WIN(C); CU=Ubound_simple(C)
    # is U an upper bound for WIN? is U maximized by consec? is U tight at consec?
    bad_ub=0; beat_U=0
    for E in bank:
        E=list(E); w=WIN(E); u=Ubound_simple(E)
        if w>u+F(1,10**12): bad_ub+=1
        if u>CU+F(1,10**12): beat_U+=1
    print(f"k={k}: consec WIN={float(Cwin):.6f} U(consec)={float(CU):.6f} "
          f"tight={Cwin==CU}; WIN>U violations={bad_ub}; U-beaters(over consec)={beat_U}")
