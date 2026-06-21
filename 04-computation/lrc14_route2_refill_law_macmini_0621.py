#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROUTE 2 -- THE REFILL LAW. We have the EXACT closed form:
  going RIGHT, sector s is covered on union of [j/(7e),(j+1)/(7e)) over clocks e with
  (e*a+j)==s (mod 7). T+ = min_s first-gap(s).

KEY SIMPLIFICATION (full-residue): at y=0 sector s is occupied by clocks e with
  e*a==s (mod7), i.e. residue r=s*a^{-1}. The native occupant LEAVES at y=1/(7*e_native)
  where e_native = the clock(s) of residue r. For consec, residue r in 1..6 has a UNIQUE
  clock e=r; residue 0 has e=0 (stationary, never leaves -> sector (0) PERMANENTLY covered
  -- the PERFECT SHORT) AND e=7 (the doubled clock). Actually for consec k=8 residue 0 has
  {0,7}.

  After native leaves, sector s is refilled by the clock e' with (e'*a) == s-1 (mod7) i.e.
  residue r'=(s-1)*a^{-1}, which ENTERS s at y=1/(7 e') ... but only if 1/(7e') <= 1/(7 e_nat)
  (arrive before native leaves) OR shortly after. The gap in s = interval where s empty.

  CLAIM: for consec, EVERY sector except the binding one is refilled in time, and the
  binding sector is the one whose native clock is FASTEST (largest |e| among the slow legs).

Let's compute, for consec, the FULL per-sector cover timeline going right, and confirm:
T+ = 1/(7 * b+) where b+ = the magnitude of the BINDING native clock = the one sector
that empties first. Then express b+(a) in closed form.
"""
import sys
from fractions import Fraction as F
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)
HALF=F(1,14)

def sector_timeline_right(E,a):
    """For each sector s, the sorted list of (entry,exit,clock) covering intervals in [0,HALF]."""
    cov=defaultdict(list)
    for e in E:
        if e==0:
            cov[0].append((F(0),F(10),'0*')); continue   # stationary, covers sector 0 forever
        j=0
        while F(j,7*e)<HALF:
            cov[(e*a+j)%7].append((F(j,7*e),F(j+1,7*e),f'{e}')); j+=1
    return cov

def first_gap(intervals):
    iv=sorted(intervals); cur=F(0); last=None
    for lo,hi,c in iv:
        if lo>cur: return cur,last
        if hi>cur: last=c
        cur=max(cur,hi)
    return min(cur,HALF),last

def Tplus_and_binding(E,a):
    cov=sector_timeline_right(E,a)
    deaths={s:first_gap(cov[s]) for s in range(7)}
    T=min(d[0] for d in deaths.values()); T=min(T,HALF)
    binders=[(s,deaths[s][1]) for s in range(7) if deaths[s][0]==min(d[0] for d in deaths.values())]
    return T,binders

def consec(k): return list(range(k))

for k in (8,9,10):
    C=consec(k)
    print(f"\n### consec k={k}: per-a binding native clock (RIGHT)")
    for a in range(1,7):
        T,binders=Tplus_and_binding(C,a)
        b = 1/(7*T) if T not in (F(0),HALF) else ('CAP' if T==HALF else 0)
        print(f"  a={a}: T+={T} b+={b} binders(sector,clock)={binders}")
