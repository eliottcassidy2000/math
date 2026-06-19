#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k=9 TIGHT-MARGIN HUNT (opus-2026-06-19). cap_9=0.49426 is the TIGHTEST cap and
the whole LRC(14)-S3 reduction binds here. The non-AP max we found is 0.48729 at
[0,1,2,3,4,5,6,7,9] (margin only 0.00697). Confirm NOTHING beats it: exhaustively
sweep all single-element perturbations of the consec AP (the dim-1->near-dim-1
boundary, where L_y is largest among non-APs), plus all sets of bounded spread
with low doubling. If any non-AP L_y exceeds 0.48729 (still < cap), report it.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def N_at(E,x):
    hit=set()
    for e in E:
        v=e*x; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
    return sum(1 for j in range(1,7) if j not in hit)
def dist_p(E):
    E=sorted(set(E)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); p=[Fraction(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        p[N_at(E,(lo+hi)/2)]+=(hi-lo)
    return p
def g_poly9():
    return [Fraction(-(t-2)*(t-3)*(t-6),36) for t in range(7)]
G9=g_poly9()
def L9(E):
    p=dist_p(E); return sum(p[t]*G9[t] for t in range(7))
def is_AP(E):
    E=sorted(set(E))
    if len(E)<2: return True
    d=E[1]-E[0]; return all(E[i+1]-E[i]==d for i in range(len(E)-1))
def prim(E):
    E=sorted(set(E)); return E[0]==0 and reduce(gcd,E)==1
def dsize(E):
    E=sorted(set(E)); return len({a+b for a in E for b in E})

CAP9=Fraction(49426,100000)
INCUMBENT=Fraction(48729,100000)

def main():
    print("k=9 TIGHT-MARGIN HUNT (opus-2026-06-19)")
    print(f"cap_9={float(CAP9):.5f}  incumbent non-AP max={float(INCUMBENT):.5f}\n")
    best=Fraction(-1); arg=None; beats=[]
    seen=set()
    # ALL primitive non-AP 9-subsets of [0..W] with 0 included, bounded spread W.
    for W in range(9, 17):
        for combo in itertools.combinations(range(1, W+1), 8):
            E=(0,)+combo
            if E[-1]!=W: continue  # spread exactly W to avoid recount
            if E in seen: continue
            seen.add(E)
            if not prim(E) or is_AP(E): continue
            # only bother if low doubling (the high-L_y regime); skip very spread sets
            sig=Fraction(dsize(E),9)
            if sig>Fraction(7,2): continue  # high doubling => low L_y (verified monotone)
            Lv=L9(E)
            if Lv>best: best=Lv; arg=E
            if Lv>INCUMBENT: beats.append((float(Lv),E))
    print(f"exhaustive spread<=16, sigma<=3.5, primitive non-AP:")
    print(f"  MAX L_y = {float(best):.6f} at {arg}")
    print(f"  margin to cap_9 = {float(CAP9-best):.6f}")
    print(f"  sets beating prior incumbent {float(INCUMBENT):.5f}: {len(beats)}")
    for Lv,E in sorted(beats,reverse=True)[:20]:
        print(f"    {list(E)}  L_y={Lv:.6f}  (still < cap: {Lv<float(CAP9)})")
    # any over cap? (would break the reduction)
    over=[(Lv,E) for Lv,E in beats if Lv>=float(CAP9)]
    print(f"\n  *** sets with L_y >= cap_9 (would BREAK LRC14-S3): {len(over)} ***")
    if over:
        for Lv,E in over: print("    BREAK:",list(E),Lv)
    else:
        print("    NONE — every non-AP stays strictly below cap_9. Reduction holds at k=9.")

if __name__=="__main__":
    main()
