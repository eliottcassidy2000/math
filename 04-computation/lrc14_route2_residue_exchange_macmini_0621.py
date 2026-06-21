#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROUTE 2 -- RESIDUE-PRESERVING EXCHANGE. The right poset is NOT raw magnitude. The
correct exchange keeps the RESIDUE assignment and the DOUBLING pattern, only changing
the magnitude WITHIN a residue class. consec uses the AP refill order: residue r gets
the SMALLEST available magnitudes.

DEFINITION (refill order). A full-residue shape assigns to each residue r a SET M_r of
magnitudes (>=1 each, with 0 allowed for residue 0). Full-residue: every r in Z/7 has
>=1 magnitude. k = total count. consec: M_0={0,..} M_r={r, r+7,...} taking smallest k.
Precisely consec_k = the k smallest non-negative integers; its residue->magset is the
'greedy AP' filling.

EXCHANGE STEP (the real one): take a shape E, pick a residue r and a magnitude e in M_r,
and DECREASE it to the smallest unused magnitude e' < e with e'==r (mod7) (a 'down-shift'
toward consec). Claim: this NEVER decreases WIN (moving toward consec increases WIN).
Equivalently every shape can be reduced to consec by down-shifts, each non-decreasing WIN.

TEST: build the down-shift graph; check WIN is monotone NON-decreasing along down-shifts
(toward consec). Down-shift = pick e in E, replace by e-7 if e-7>=0 (or >=1 for r!=0),
e-7 not in E, stays full-residue (it will, same residue)."""
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

def downshifts(E):
    """All single down-shifts e -> e-7 (same residue), keeping distinct & full-residue & >=0."""
    res=[]
    Es=set(E)
    for e in E:
        e2=e-7
        if e2<0: continue
        if e2 in Es: continue
        cand=sorted((Es-{e})|{e2})
        if is_full_residue(cand) and len(cand)==len(E):
            res.append(cand)
    return res

def test(k,span):
    C=consec(k)
    bank=set(c for c in itertools.combinations(range(0,span+1),k) if 0 in c and is_full_residue(c))
    wins={E:WIN(list(E)) for E in bank}
    viol=0; edges=0; ex=[]
    for E in bank:
        for D in downshifts(list(E)):
            td=tuple(D)
            if td not in bank: continue   # restrict to bank
            edges+=1
            # down-shift toward consec should NOT decrease WIN
            if wins[td] < wins[E] - F(1,10**15):
                viol+=1
                if len(ex)<8: ex.append((E,td,wins[E],wins[td]))
    print(f"k={k} span<={span}: bank={len(bank)} downshift-edges={edges} "
          f"WIN-DECREASE-on-downshift violations={viol}")
    for E,td,w,wd in ex:
        print(f"   VIOL {list(E)} -down-> {list(td)}: WIN {float(w):.6f} -> {float(wd):.6f} (DECREASED toward consec)")
    return viol

for k,span in [(8,15),(9,17),(10,19),(11,20)]:
    test(k,span)
