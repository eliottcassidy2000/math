#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROUTE 2 -- EXPLICIT binding speed b+(a), b-(a) in terms of the magnitude profile.

We proved: going right, sector s is first covered by NATIVE clocks (those with
e*a==s mod7) starting at y=0; native occupant e leaves at y=1/(7e). Sector s gets
REFILLED by clock e' with e'*a==s-1 (mod7) arriving at y=1/(7e') (it leaves its own
sector s-1 and enters s). More generally clock e enters sector (e*a+j) at j/(7e).

For the FULL-RESIDUE stratum, residue r=e mod7 has a set of magnitudes M_r={|e|:e==r}.
Define for each residue r the SORTED magnitudes. The native occupant of sector s at
y=0 are clocks with residue r=s*a^{-1}. 

CLAIM (binding-speed law, right side):
  T+(a) = min_{s in Z/7}  death_right(s)
  where death_right(s) is computed by the per-sector interval cover:
   - contributions to sector s: from each residue r, the clock(s) e (|e| in M_r,
     e==r) contribute interval [ ((s - r*a) * inv ... )]. Cleaner: clock e covers s
     at the unique j in {0,..} with (e*a+j)==s mod7, i.e. j = (s - e*a) mod 7, on
     [j/(7e),(j+1)/(7e)) and again at j+7, j+14,... 
  The FIRST GAP in sector s determines death_right(s).

We want: b+(a)=1/(7 T+(a)). For consec, T+ = 1/42 (b=6) except a=6 where T+=2/49.
Let's express b+ and b- as explicit functions and CHECK the "two fastest clocks"
mechanism: consec routes ~all binding to e=6 (right) and e=k-1 (the doubled-resid-0
clock, left).  Then study the harmonic sum.
"""
import sys, itertools
from fractions import Fraction as F
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)
HALF=F(1,14)

def Tplus_cf(E,a):
    cover=defaultdict(list)
    for e in E:
        if e==0: cover[0].append((F(0),HALF)); continue
        ae=abs(e)
        # sign of drift: e>0 up. assume e>=0 here.
        j=0
        while F(j,7*e)<HALF:
            lo=F(j,7*e); hi=min(F(j+1,7*e),HALF)
            cover[(e*a+j)%7].append((lo,hi)); j+=1
    def fg(iv):
        iv=sorted(iv); cur=F(0)
        for lo,hi in iv:
            if lo>cur: return cur
            cur=max(cur,hi)
        return cur
    return min(fg(cover[s]) for s in range(7))

def Tminus_cf(E,a):
    cover=defaultdict(list)
    for e in E:
        if e==0: cover[0].append((F(0),HALF)); continue
        j=0
        while F(j,7*e)<HALF:
            lo=F(j,7*e); hi=min(F(j+1,7*e),HALF)
            cover[(e*a-1-j)%7].append((lo,hi)); j+=1
    def fg(iv):
        iv=sorted(iv); cur=F(0)
        for lo,hi in iv:
            if lo>cur: return cur
            cur=max(cur,hi)
        return cur
    return min(fg(cover[s]) for s in range(7))

def consec(k): return list(range(k))

# Which sector binds, and the binding magnitude, for consec
def analyze_consec(k):
    C=consec(k)
    print(f"\n=== consec k={k} : binding sector & magnitude per (a,side) ===")
    for a in range(1,7):
        # right
        coverR=defaultdict(list)
        for e in C:
            if e==0: coverR[0].append((F(0),HALF,'STAT-0')); continue
            j=0
            while F(j,7*e)<HALF:
                coverR[(e*a+j)%7].append((F(j,7*e),min(F(j+1,7*e),HALF),f'e{e}')); j+=1
        deaths={}
        for s in range(7):
            iv=sorted(coverR[s]); cur=F(0); lastclock=None
            for lo,hi,tag in iv:
                if lo>cur: break
                if hi>cur: lastclock=tag
                cur=max(cur,hi)
            deaths[s]=(cur,lastclock)
        Tp=min(d[0] for d in deaths.values())
        binds=[(s,d[1]) for s,d in deaths.items() if d[0]==Tp]
        b = (1/(7*Tp)) if Tp not in (F(0),HALF) else ('CAP' if Tp==HALF else 'X')
        print(f"  a={a} RIGHT: T+={Tp} b+={b} binding sector(s)={binds}")

if __name__=="__main__":
    print("="*80); print("ROUTE 2: explicit binding speed & sector for consec"); print("="*80)
    for k in (8,9,10): analyze_consec(k)
