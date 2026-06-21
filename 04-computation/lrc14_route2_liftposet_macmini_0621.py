#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ROUTE 2 -- LIFT-POSET monotonicity (efficient). The poset on full-residue shapes:
E1 <= E2 iff sorted(E1) <= sorted(E2) componentwise. consec is the unique minimum.
We test WIN monotone-DECREASING along COVER (Hasse) edges only: E -> E' where E' is
obtained from E by a minimal upward step (increase ONE sorted-coordinate to the next
admissible value keeping full-residue & distinct). If WIN strictly decreases along
EVERY cover edge, then WIN is monotone on the whole poset (transitivity), hence consec
(the bottom) is the unique max. This is the exchange lemma in poset form."""
import sys, itertools
from fractions import Fraction as F
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

def cover_edges_test(k, span):
    C=consec(k)
    bank=set(c for c in itertools.combinations(range(0,span+1),k)
             if 0 in c and is_full_residue(c))
    wins={E:WIN(list(E)) for E in bank}
    # cover edge: E' obtained from E by raising one coordinate to a larger value, 
    # with NO admissible intermediate in bank. Equivalently for each E, for each i,
    # raise E[i] to the smallest v>E[i] such that result is in bank (distinct,full-res),
    # and v is a Hasse cover if no w in (E[i],v) gives an in-bank shape <= E'.
    # Simpler & still valid: check ALL single-coordinate raises that land in bank
    # (these generate the order); monotone on these => monotone overall is NOT guaranteed
    # unless they're covers. We check the MINIMAL raise per coordinate (true Hasse cover
    # in the product order restricted to bank is the min raise).
    viol=0; edges=0; ex=[]
    blist=sorted(bank)
    for E in blist:
        Es=list(E)
        for i in range(k):
            # raise coordinate i to next admissible value keeping sorted & distinct & full-res
            for v in range(Es[i]+1, span+1):
                cand=sorted(Es[:i]+[v]+Es[i+1:])
                if len(set(cand))!=k: continue
                t=tuple(cand)
                if t in bank:
                    edges+=1
                    if wins[t] >= wins[E]:   # should strictly DEcrease
                        viol+=1
                        if len(ex)<6: ex.append((E,t,wins[E],wins[t]))
                    break  # minimal raise = cover
    print(f"k={k} span<={span}: bank={len(bank)} cover-edges={edges} "
          f"WIN-increase-or-tie violations={viol}")
    for E,t,w,wt in ex:
        print(f"   VIOL {list(E)} -> {list(t)}: WIN {float(w):.6f} -> {float(wt):.6f} (not decreasing)")
    return viol

for k,span in [(8,14),(9,16),(10,16)]:
    cover_edges_test(k,span)
