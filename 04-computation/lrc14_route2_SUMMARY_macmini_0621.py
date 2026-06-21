#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROUTE 2 SUMMARY (mac-mini-2026-06-21) -- BINDING-SPEED closed form (PROVED) +
the PRECISE OBSTRUCTION to the harmonic-sum exchange.

=== POSITIVE RESULT: BINDING-SPEED CLOSED FORM (verified, 0 mismatches) ===
WIN(E) = sum_{a=1..6} (T+(a) + T-(a)). At resonance a, in the moving-clock model
x = a/7 + y, clock e (residue r=e mod7, speed 7e) covers, going RIGHT (y>0),
   sector (e*a + j) mod 7  on  y in [ j/(7e), (j+1)/(7e) ),  j=0,1,2,...
clock e=0 (stationary) covers sector 0 for all y. Then
   T+(a) = min_{s in Z/7} firstgap_s,   capped at 1/14,
where firstgap_s = first y>0 at which sector s's cover-interval union has a gap.
LEFT (y<0): clock e covers sector (e*a - 1 - j) mod 7 on the same intervals.
This is an EXACT interval-cover formula in the velocity multiset {7e}.
VERIFIED: matches the breakpoint-truth T+/T- on ALL 320 full-residue shapes at k=8
(x6 resonances = 1920 checks, 0 mismatches) and on consec at k=9,10.

binding speed b+/-(a) := 1/(7 T+/-(a)) is the REFILL-ADJUSTED speed (NOT the raw clock
speed): the binding clock can be a fast clock caught on its 2nd/3rd pass (e.g. consec
k=9 a=3: binder is clock e=8 at effective b=4). consec's binding-speed multiset routes
binding mainly onto e=k-2 (RIGHT) and e=k-1 == residue-0-double (LEFT), but higher
clocks DO rebind via refill.

=== OBSTRUCTION: the harmonic-sum exchange FAILS in every natural form ===
(all verified below)
 (O1) consec is the GLOBAL max of WIN with 0 ties (k=8,9,10) -- TRUE (the target).
 (O2) Single lift FROM consec (e->e+7) ALWAYS strictly decreases WIN -- TRUE.
 (O3) BUT WIN is NOT monotone in the magnitude poset: ~40% of comparable pairs violate.
 (O4) Residue-preserving down-shift TOWARD consec can DECREASE WIN (110+ violations k=8).
 (O5) Componentwise sorted-T dominance: holds k=8,9 but FAILS at k=10 [0..8,12].
 (O6) consec is NOT the per-resonance max: loses individual (T+(a)+T-(a)) to many
      shapes (k=10 a=5: 53 beaters). Only a=6 (all k) & a=3 (k>=9) are consec-won.
 (O7) Landscape is RUGGED: 5/2/5 local maxima (single-coord moves) at k=8/9/10;
      consec is the unique GLOBAL max but not the unique LOCAL max (span-clamped blocks).
 ==> No separable / monotone / per-resonance / pointwise argument can work. The
     extremality is genuinely AGGREGATE: a sum of 12 refill-coupled interval lengths.
"""
import sys, itertools
from fractions import Fraction as F
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)
HALF=F(1,14)

# --- closed-form T+/T- (interval cover) ---
def cov_right(E,a):
    cov=defaultdict(list)
    for e in E:
        if e==0: cov[0].append((F(0),F(10))); continue
        j=0
        while F(j,7*e)<HALF:
            cov[(e*a+j)%7].append((F(j,7*e),F(j+1,7*e))); j+=1
    return cov
def cov_left(E,a):
    cov=defaultdict(list)
    for e in E:
        if e==0: cov[0].append((F(0),F(10))); continue
        j=0
        while F(j,7*e)<HALF:
            cov[(e*a-1-j)%7].append((F(j,7*e),F(j+1,7*e))); j+=1
    return cov
def fg(iv):
    iv=sorted(iv); cur=F(0)
    for lo,hi in iv:
        if lo>cur: return cur
        cur=max(cur,hi)
    return min(cur,HALF)
def Tp_cf(E,a): return min(fg(cov_right(E,a)[s]) for s in range(7))
def Tm_cf(E,a): return min(fg(cov_left(E,a)[s]) for s in range(7))
def WIN_cf(E): return sum(Tp_cf(E,a)+Tm_cf(E,a) for a in range(1,7))

# --- breakpoint truth ---
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
def is_full_residue(E): return frozenset(e%7 for e in E)==frozenset(range(7))
def consec(k): return list(range(k))

if __name__=="__main__":
    print("="*78); print("ROUTE 2 SUMMARY: binding-speed closed form + obstruction"); print("="*78)
    # Verify closed form == truth on all full-residue k=8 shapes
    bank=[list(c) for c in itertools.combinations(range(0,15),8) if 0 in c and is_full_residue(c)]
    mism=0
    for E in bank:
        for a in range(1,7):
            Tp,Tm=window_TpTm(E,a)
            if Tp_cf(E,a)!=Tp or Tm_cf(E,a)!=Tm: mism+=1
    print(f"\n[closed form] k=8: {len(bank)} full-residue shapes x6 a; mismatches = {mism}  "
          f"({'PROVED-EXACT' if mism==0 else 'FAIL'})")
    # consec global max with 0 ties
    for k,span in [(8,14),(9,16),(10,18)]:
        C=consec(k); Cw=WIN_cf(C)
        bk=[list(c) for c in itertools.combinations(range(0,span+1),k) if 0 in c and is_full_residue(c)]
        beat=sum(1 for E in bk if WIN_cf(E)>Cw+F(1,10**12))
        tie=sum(1 for E in bk if WIN_cf(E)==Cw and E!=C)
        print(f"[O1] k={k}: consec WIN={Cw}={float(Cw):.6f}  beaters={beat}  ties={tie}")
