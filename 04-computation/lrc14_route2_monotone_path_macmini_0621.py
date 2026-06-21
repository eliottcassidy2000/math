#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROUTE 2 -- does a WIN-INCREASING PATH to consec exist from every shape?
Local moves: any single magnitude change e->e' (same residue, stay distinct & full-res,
within span). From each shape, is there a neighbor with STRICTLY higher WIN, until consec?
i.e. is consec the unique 'peak' reachable by greedy ascent? If YES from EVERY shape,
then consec is the unique local max AND globally attracting -> strong structural result
(even though WIN isn't monotone in a fixed poset, it has a unique basin).
We test: (1) is consec the UNIQUE local max (no other shape has all neighbors <= it)?
         (2) does greedy steepest ascent from every shape reach consec?
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

def neighbors(E, span):
    """single magnitude change e->e' (any residue, any value 0..span), keep distinct & full-res."""
    Es=set(E); out=[]
    for e in list(E):
        for v in range(0, span+1):
            if v==e or v in Es: continue
            cand=sorted((Es-{e})|{v})
            if len(cand)==len(E) and is_full_residue(cand):
                out.append(tuple(cand))
    return out

def test(k, span):
    C=tuple(consec(k))
    bank=set(c for c in itertools.combinations(range(0,span+1),k) if 0 in c and is_full_residue(c))
    wins={E:WIN(list(E)) for E in bank}
    Cwin=wins[C]
    # (1) unique local max?
    localmax=[]
    for E in bank:
        nb=[n for n in neighbors(list(E),span) if n in bank]
        if all(wins[n] <= wins[E] for n in nb):
            localmax.append(E)
    print(f"k={k} span<={span}: bank={len(bank)} #local-maxima(single-coord, full-res)= {len(localmax)}")
    for lm in localmax[:10]:
        print(f"   localmax {list(lm)} WIN={float(wins[lm]):.6f} {'<-- CONSEC' if lm==C else ''}")
    return localmax,C

for k,span in [(8,14),(9,16),(10,16)]:
    test(k,span)
