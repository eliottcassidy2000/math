#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROUTE 2 -- the MECHANISM: why e=6 binds right and e=7 binds left for consec.

THE CLEAN MODEL (going RIGHT, y in (0, 1/14)):
  At y=0, x=a/7, clock e sits on a sector BOUNDARY (e*a is an integer, so 7 frac(e x)
  = e*a mod 7 = integer => boundary between sector (e*a-1 mod7) and (e*a mod7)).
  Convention: floor of an integer-valued point = that integer, so at y=0+epsilon the
  clock is in sector s_e = (e*a) mod 7, having JUST entered from below.
  Going right it drifts UP; it will LEAVE sector s_e and enter s_e+1 at the y where
  7*e*y = 1, i.e. y = 1/(7e).  So clock e occupies sector (e*a mod 7) on the interval
  y in [0, 1/(7e)).

  Residues = Z/7 => the 7 residue classes hit all 7 sectors at y=0 (a invertible).
  EACH sector s = (e*a mod 7) is occupied at y=0+ by EXACTLY the clocks e' with
  e'*a == s (mod 7), i.e. e' == s*a^{-1} (mod 7). Among those, the LONGEST occupant
  (last to leave) is the SMALLEST magnitude |e'| (slowest, leaves at 1/(7|e'|), largest).
  ==> sector s stays covered until y = 1/(7 * minmag(s)) where minmag(s)=min{|e'|:e'a==s},
  UNLESS a NEW clock drifts INTO sector s before then (refill).

  REFILL: clock e'' enters sector s at y where it crosses from s-1 to s:
  7*e''*y = 1 + (height needed). e'' started in sector (e''*a mod7)=s-1 requires... 
  Actually e'' enters s at y = 1/(7 e'') if (e''*a mod7) = s-1. So refill of s comes
  from the clock e'' with residue (s-1)*a^{-1}, arriving at 1/(7 e'').

  So sector s DIES (becomes empty with no refill) at the first y where its last native
  occupant left and no refiller present.

  THE BINDING SECTOR = the sector whose (native-departure) happens earliest AND whose
  refill (if any) arrives later. T+ = min over sectors of that death time.

  consec right: native occupant of sector s=(e*a mod7) is the clock e=s*a^{-1} mod 7
  in {0..6} PLUS possibly e=7 (residue 0). The FASTEST native departure (smallest T)
  is the sector occupied by the LARGEST-magnitude slowest... no: the sector with the
  LARGEST min-mag dies first. For consec, residues 1..6 have min-mag = the residue
  itself (1..6); residue 0 has min-mag 0 (e=0 NEVER moves -> sector (0) always covered
  by e=0!). So the residue-0 sector is a PERFECT SHORT (e=0 stationary). The other
  6 sectors have min-mag = their residue 1..6. The sector dying first is residue with
  the LARGEST min-mag among 1..6 = 6, dying at y=1/(7*6)=1/42. THAT's b+=6. checks out.

This script verifies the closed form:
  T+(a) = min over sectors s != 0-sector of [ death(s) ], death(s)= time last native
          occupant leaves OR is refilled. And derives b+ = max usable min-mag.
"""
import sys, itertools
from fractions import Fraction as F
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)
HALF=F(1,14)

def sector_of(e,a,y):
    pos=F(e*a)+F(7*e)*y
    return (pos.numerator//pos.denominator)%7

def covered(E,a,y): return len({sector_of(e,a,y) for e in E})==7

def breakpoints(E,a):
    bps={F(0),-HALF,HALF}
    for e in E:
        if e==0: continue
        lo_val=F(7*e)*(-HALF)+F(e*a); hi_val=F(7*e)*(HALF)+F(e*a)
        lo_i=min(lo_val,hi_val); hi_i=max(lo_val,hi_val)
        m=lo_i.numerator//lo_i.denominator
        while m<=hi_i.numerator//hi_i.denominator+1:
            y=F(m-e*a,7*e)
            if -HALF<=y<=HALF: bps.add(y)
            m+=1
    return sorted(bps)

def window_TpTm(E,a):
    bps=breakpoints(E,a); ivals=list(zip(bps,bps[1:]))
    Tp=F(0)
    for lo,hi in ivals:
        if hi<=0: continue
        lo2=max(lo,F(0))
        if covered(E,a,(lo2+hi)/2): Tp=hi
        else:
            if lo2==F(0): Tp=F(0)
            break
    Tp=min(Tp,HALF)
    Tm=F(0)
    for lo,hi in reversed(ivals):
        if lo>=0: continue
        hi2=min(hi,F(0))
        if covered(E,a,(lo+hi2)/2): Tm=-lo
        else:
            if hi2==F(0): Tm=F(0)
            break
    Tm=min(Tm,HALF)
    return Tp,Tm

def consec(k): return list(range(k))
def is_full_residue(E): return frozenset(e%7 for e in E)==frozenset(range(7))

# --- CLOSED FORM via interval-cover per sector (exact, no breakpoint scan) ----
def Tplus_closedform(E,a):
    """Going right. Each clock e>0 covers sector (e*a + j) mod 7 on
       y in [j/(7e), (j+1)/(7e)) for j=0,1,2,...; clock e=0 covers sector 0 forever.
       T+ = sup{ Y<=1/14 : every sector covered on [0,Y] continuously from 0 }.
       Compute as min over sectors of the first GAP-start where sector loses cover."""
    # build, per sector, the union of covering y-intervals within [0, HALF]
    cover=defaultdict(list)  # sector -> list of (lo,hi)
    for e in E:
        if e==0:
            cover[0].append((F(0),HALF)); continue
        # clock e (assume e>0 for consec) in sector (e*a + j)mod7 for y in [j/(7e),(j+1)/(7e))
        j=0
        while F(j,7*e) < HALF:
            lo=F(j,7*e); hi=min(F(j+1,7*e),HALF)
            cover[(e*a+j)%7].append((lo,hi))
            j+=1
    # for each sector, find first y in (0,HALF] where it is NOT covered
    def first_gap(secivs):
        # merge intervals, find first uncovered point >0
        ivs=sorted(secivs)
        cur=F(0)
        for lo,hi in ivs:
            if lo>cur:  # gap (cur,lo)
                return cur
            cur=max(cur,hi)
        return cur  # covered up to cur (then gap if cur<HALF, else full)
    deaths=[first_gap(cover[s]) for s in range(7)]
    return min(deaths)

def Tminus_closedform(E,a):
    """Going left, y<0. clock e covers sector (e*a + j) for y in ((j-1)/(7e) ... wait
       for y<0 going down, sector decreases. Use symmetry: y->-y is e->-e, i.e. the
       same as Tplus with E negated. sector_of(-e,a,y)=( -e*a -7e y)... = sector_of(e,a,-y)
       reversed. Easiest: mirror via reusing Tplus on reflected config."""
    # going left = going right with all velocities reversed. The cover pattern:
    cover=defaultdict(list)
    for e in E:
        if e==0:
            cover[0].append((F(0),HALF)); continue
        # for y<0, let t=-y>0; clock e at y=-t sits at e*a -7e t. covers sector
        # (e*a - 1 - j)mod7 ... at t=0+ it's leaving sector e*a downward -> in (e*a-1)mod7.
        # crosses down each 1/(7e). So sector (e*a-1-j)mod7 on t in [j/(7e),(j+1)/(7e)).
        j=0
        while F(j,7*e)<HALF:
            lo=F(j,7*e); hi=min(F(j+1,7*e),HALF)
            cover[(e*a-1-j)%7].append((lo,hi))
            j+=1
    def first_gap(secivs):
        ivs=sorted(secivs); cur=F(0)
        for lo,hi in ivs:
            if lo>cur: return cur
            cur=max(cur,hi)
        return cur
    deaths=[first_gap(cover[s]) for s in range(7)]
    return min(deaths)

if __name__=="__main__":
    print("="*80)
    print("ROUTE 2: closed-form T+/T- (sector interval-cover) vs breakpoint truth")
    print("="*80)
    allok=True
    for k in (8,9,10):
        C=consec(k)
        print(f"\n### k={k} consec")
        for a in range(1,7):
            Tp,Tm=window_TpTm(C,a)
            Tpc=Tplus_closedform(C,a); Tmc=Tminus_closedform(C,a)
            ok = (min(Tpc,HALF)==Tp) and (min(Tmc,HALF)==Tm)
            allok = allok and ok
            print(f"  a={a}: T+ truth={Tp} closed={min(Tpc,HALF)} {'OK' if min(Tpc,HALF)==Tp else 'XX'}"
                  f"   T- truth={Tm} closed={min(Tmc,HALF)} {'OK' if min(Tmc,HALF)==Tm else 'XX'}")
    print(f"\nALL closed-form match: {allok}")

    # Now verify on RANDOM full-residue shapes (not just consec) to validate the formula
    print("\n--- validate closed form on full-residue rivals (k=8) ---")
    import random
    bank=[(0,)+c for c in itertools.combinations(range(1,15),7)]
    bank=[E for E in bank if is_full_residue(E)]
    bad=0
    for E in bank:
        E=list(E)
        for a in range(1,7):
            Tp,Tm=window_TpTm(E,a)
            Tpc=min(Tplus_closedform(E,a),HALF); Tmc=min(Tminus_closedform(E,a),HALF)
            if Tpc!=Tp or Tmc!=Tm:
                bad+=1
                if bad<=5:
                    print(f"  MISMATCH E={E} a={a}: T+({Tp} vs {Tpc}) T-({Tm} vs {Tmc})")
    print(f"  total mismatches over {len(bank)} full-residue shapes x6 a: {bad}")
