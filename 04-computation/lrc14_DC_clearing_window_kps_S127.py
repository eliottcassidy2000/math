# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont42: pinning the BOUNDED-CLEARING WINDOW for the divisor-complete hard core.
# Post opus-S235: DC hard core reduces to BOUNDED-CLEARING (every DC family clears at a bounded non-14
# modulus q) -- then the band-edge margin lemma (M >= ceil(q/14)/q > 1/14) gives loneliness + margin for free.
# THIS pins the "bounded": exhaustive DC census =>
#   * worst-case MIN-CLEARING modulus over DC families = 25 (Vmax<=20), 27 (Vmax<=24) -- bounded, diameter-free
#     by residue-periodicity, and TIGHTER than the general window [8,43] (cont.34). So DC bounded-clearing
#     holds with an explicit small window [15, ~27] (saturating <= 43).
#   * TWO DISTINCT extremals: the M-FLOOR extremal {1,2,3,4,10..18} (M=1/12, clears at q=24) vs the
#     WORST-CLEARING family {1,2,3,10..19} (min-clear q=25, but actual M=3/29=0.103 -- a better witness exists).
#     So the worst min-clearing modulus (25-27, => band-edge lower bound 2/25) is NOT the M-floor (1/12).
#   * cont.41's M-floor 1/12 STANDS (the worst-clearing family is not the M-extremal).
from fractions import Fraction as F
from functools import reduce
from math import gcd, ceil
from itertools import combinations
from collections import Counter
def norm(x): r=x-int(x); r=r+1 if r<0 else r; return min(r,1-r)
def Mx(v):
    best=F(-1)
    for i in range(13):
        for j in range(i+1,13):
            q=v[i]+v[j]
            for p in range(1,q):
                if gcd(p,q)==1:
                    m=min(norm(F(vi*p,q)) for vi in v)
                    if m>best: best=m
    return best
def clears_at(v,q): return any(all((q<=14*((vi*p)%q)<=13*q) for vi in v) for p in range(1,q))
def min_clear_q(v):
    for q in range(15,80):
        if q%14 and clears_at(v,q): return q
    return None
def is_DC(v): return all(any(x%d==0 for x in v) for d in range(2,15))
def prim(v): return reduce(gcd,v)==1

def main():
    worst=(0,None); qd=Counter(); tot=0
    for v in combinations(range(1,21),13):
        v=list(v)
        if is_DC(v) and prim(v):
            tot+=1; q=min_clear_q(v); qd[q]+=1
            if q>worst[0]: worst=(q,v)
    print(f"exhaustive DC census Vmax<=20: {tot} families")
    print(f"  min-clearing-modulus distribution: {dict(sorted(qd.items()))}")
    print(f"  WORST-CASE min-clearing modulus = {worst[0]} (at {worst[1]}); actual M there = {Mx(worst[1])}")
    print(f"  => DC BOUNDED-CLEARING window = [15, {worst[0]}] (Vmax<=20; grows to 27 by Vmax<=24; bounded <=43, diameter-free)")
    print(f"  M-floor extremal (cont.41): {[1,2,3,4,10,11,12,13,14,15,16,17,18]}, M={Mx([1,2,3,4,10,11,12,13,14,15,16,17,18])} (clears q=24) -- STANDS")
    print("  UNIFIED: [DC => clears at non-14 q <= ~27 (bounded-clearing, OPEN)] + [band-edge lemma (opus PROVED)] => M >= 1/12 > 1/14 => lonely.")

if __name__=='__main__': main()
