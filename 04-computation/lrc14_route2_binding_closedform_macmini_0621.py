#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROUTE 2 (mac-mini-2026-06-21) -- BINDING-SPEED CLOSED FORM + harmonic-sum exchange.

WIN(E) = sum_{a=1..6} (T+(a) + T-(a)), the CONTIGUOUS-through-center survival
window at each resonance a/7.  At resonance a, set x = a/7 + y, y in (-1/14,1/14).
Each clock e (residue r_e = e mod 7, magnitude |e|) sits at sector position
   pos_e(y) = e*a + 7*e*y   (real); its sector = floor(frac of that / ... ) but in
   the moving-clock model: sector_e(y) = floor( (e*a + 7*e*y) mod 7 ).
At y=0, sector_e = (e*a) mod 7 = (r_e * a) mod 7. Since residues are ALL of Z/7,
the 7 clocks' base sectors {r*a mod 7 : r in Z/7} = Z/7 exactly (a invertible).
Plus the DOUBLED residue (consec: residue 0 by e=0 and e=7) gives an extra clock.

GOAL 1: closed form for T+(a) (and T-(a)) as an explicit function of the velocity
        multiset {7e} and the refill order.

DERIVATION (moving-clock / first-empty-sector model):
  Going RIGHT (y>0 increasing), each clock e drifts up at speed 7e (down if e<0; here
  all e>=0 so all drift up, e=0 stationary). A sector s is covered as long as >=1
  clock sits in it. Sector s loses its last occupant at the smallest y where the
  last clock in s has drifted OUT of s. The window dies at T+ = min over sectors s
  of (time for sector s to become empty).

  Clock e enters sector at offset o_e(0) = frac(e*a) ... within its sector it sits at
  height h_e in [0,1) (the fractional part of (e*a + 7*e*y) ). It LEAVES its current
  sector when h_e reaches 1, i.e. at y where 7*e*y carries it across: the clock
  leaves at y = (1 - h_e(0)) / (7 e)   [for e>0]. Then it ENTERS the next sector.

  A sector s becomes EMPTY at the moment its last occupant leaves AND no new clock
  has entered. With residues = Z/7 and a invertible, at y=0 every sector has exactly
  the clocks e with (e*a mod 7) = s. The DEPARTURE time of clock e from its sector is
     d_e = (1 - frac(e*a)) / (7 e)     (e>0)   [time to exit current sector going up]
  but frac(e*a) here = (e*a mod 7)/7 ... wait e*a is an integer. Need the SUB-sector
  height. At x=a/7 exactly, 7*frac(e*x) = 7*frac(e*a/7) = (e*a mod 7), an INTEGER ->
  the clock sits EXACTLY on a sector BOUNDARY. So the sub-sector height is 0 (just
  entered) going one way and 1 (about to leave) the other way. This boundary
  degeneracy is why we look at the two SIDES separately.

We will COMPUTE T+/T- exactly and FIT the closed form, then verify on all shapes.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

HALF = F(1,14)

def sector_of(e, a, y):
    pos = F(e*a) + F(7*e)*y
    fl = pos.numerator // pos.denominator
    return fl % 7

def covered(E, a, y):
    return len({sector_of(e,a,y) for e in E}) == 7

def breakpoints(E, a):
    bps = {F(0), -HALF, HALF}
    for e in E:
        if e == 0: continue
        for sign in (-1,1):
            yend = sign*HALF
            lo_val = F(7*e)*(-HALF)+F(e*a); hi_val = F(7*e)*(HALF)+F(e*a)
            lo_i = min(lo_val,hi_val); hi_i=max(lo_val,hi_val)
            m = lo_i.numerator//lo_i.denominator
            while m <= hi_i.numerator//hi_i.denominator + 1:
                y = F(m - e*a, 7*e)
                if -HALF <= y <= HALF: bps.add(y)
                m += 1
    return sorted(bps)

def window_TpTm(E, a):
    bps = breakpoints(E,a); ivals=list(zip(bps,bps[1:]))
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

# ---- Closed-form CANDIDATE for the binding speed -----------------------------
# Going RIGHT from y=0: every clock e>0 sits on the upper boundary of its sector
# (since e*a integer => about to step UP into sector (e*a+1) mod 7). The sector it
# is LEAVING is (e*a) mod 7. It leaves at y_e = 1/(7e) (needs to move height 1).
# Actually at y=0+ it's ABOUT to leave the sector below; the sector (e*a mod 7) is
# occupied on (-1/(7e), 0) going down, and on the RIGHT it's already moved up... 
# Let's just directly read off, per side, which sector dies first & at what y, and
# call b = 1/(7*T) the binding speed; then express b in terms of magnitudes.

def first_death(E, a, side):
    """Return (T, binding_sector, binding_clocks) for the given side ('+'/'-')."""
    Tp,Tm = window_TpTm(E,a)
    T = Tp if side=='+' else Tm
    if T==0: return (F(0), None, None)
    if T==HALF: return (HALF, 'CAP', None)
    # the sector that dies: just past T (on that side), which sector is missing?
    eps = F(1, 7*max(abs(e) for e in E if e)*1000)
    ytest = (T+eps) if side=='+' else -(T+eps)
    inside = {sector_of(e,a,F(0)) for e in E} if False else None
    secs_just_inside = {sector_of(e,a, (T-eps) if side=='+' else -(T-eps)) for e in E}
    secs_just_outside = {sector_of(e,a, ytest) for e in E}
    dead = secs_just_inside - secs_just_outside
    # which clock(s) were the last occupant(s) of the dead sector just inside?
    if not dead:
        return (T, 'CAP?', None)
    s = sorted(dead)[0]
    yin = (T-eps) if side=='+' else -(T-eps)
    occupants = [e for e in E if sector_of(e,a,yin)==s]
    return (T, s, occupants)

def is_full_residue(E): return frozenset(e%7 for e in E)==frozenset(range(7))
def consec(k): return list(range(k))
def WIN(E): return sum(sum(window_TpTm(E,a)) for a in range(1,7))

if __name__=="__main__":
    print("="*80)
    print("ROUTE 2: BINDING-SPEED CLOSED FORM (first-death sector & clock per side)")
    print("="*80)
    for k in (8,9,10):
        C=consec(k)
        print(f"\n### k={k}, consec={C}, WIN={WIN(C)}={float(WIN(C)):.6f}")
        for a in range(1,7):
            Tp,Tm=window_TpTm(C,a)
            fp=first_death(C,a,'+'); fm=first_death(C,a,'-')
            bp = (1/(7*Tp)) if (Tp not in (F(0),HALF)) else ('CAP' if Tp==HALF else 'X')
            bm = (1/(7*Tm)) if (Tm not in (F(0),HALF)) else ('CAP' if Tm==HALF else 'X')
            print(f"  a={a}: T+={Tp} (b+={bp}, deadsec={fp[1]}, occ={fp[2]})  "
                  f"T-={Tm} (b-={bm}, deadsec={fm[1]}, occ={fm[2]})")
