#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
LRC(14) OPEN-Q-108 -- boundary-window-finite ANGLE -- the DEFINITIVE EXHAUSTIVE
window check (companion to lrc_q108_boundary-window-finite_kps-Sx-wf.py).

This script does the FULLY EXHAUSTIVE true-wide check that the structured
families in the main script approximate: it enumerates EVERY primitive
true-wide k-set (0 in E, 2nd-largest>14, gcd=1) with bounded span, exactly,
for the tight cases k=8,9,10, and reports max p0 and any violation.

It is the rigorous core of LINK 5 (the window 15<=span<=B').

EXACT arithmetic (fractions.Fraction).  kind-pasteur workflow.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd

try:
    sys.stdout.reconfigure(encoding='utf-8', line_buffering=True)
except Exception:
    pass

CAPS = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91), 11: F(66,91), 12: F(6,7)}

def bp(E):
    s={F(0),F(1)}
    for e in E:
        if e==0: continue
        for j in range(7):
            m=0
            while True:
                xv=(F(j,7)+m)/e
                if xv>=1: break
                if xv>=0: s.add(xv)
                m+=1
    return sorted(b for b in s if 0<=b<1)

def p0(E):
    E=sorted(set(E)); B=bp(E); a=F(0)
    for lo,hi in zip(B,B[1:]+[F(1)]):
        if hi<=lo: continue
        mid=(lo+hi)/2
        miss=set(range(1,7))-set(int((e*mid)%1*7) for e in E)
        if len(miss)==0: a+=hi-lo
    return a

def exhaustive(k, SPANMAX):
    """ALL primitive true-wide k-sets: 0 in E, max<=SPANMAX, 2nd-largest>14, gcd=1."""
    cap=CAPS[k]
    maxp0=F(0); argmax=None; viol=0; cnt=0; nearrows=[]
    for combo in itertools.combinations(range(1,SPANMAX+1), k-1):
        if combo[-2] <= 14:        # 2nd largest of full set
            continue
        g=0
        for e in combo: g=gcd(g,e)
        if g!=1: continue          # primitive
        cnt+=1
        v=p0((0,)+combo)
        if v>maxp0:
            maxp0=v; argmax=(0,)+combo
        if v>cap:
            viol+=1
            print(f"  *** VIOLATION k={k}: p0={v}={float(v):.6f} > cap={float(cap):.6f} E={(0,)+combo}")
        elif cap-v < F(1,25):
            nearrows.append(((0,)+combo, v))
    return cnt, maxp0, argmax, viol, nearrows

def main():
    print("#"*78)
    print("# LRC(14) OPEN-Q-108 boundary-window-finite -- DEFINITIVE EXHAUSTIVE window")
    print("#"*78)
    print()
    print("Every primitive true-wide k-set (0 in E, 2nd-largest>14, gcd=1, span<=SPANMAX).")
    print("This is the rigorous LINK-5 finite check (no sampling).")
    print()
    cases = [(8,24),(9,22),(10,21)]
    allok=True
    for (k,SP) in cases:
        cnt,mx,am,vi,near=exhaustive(k,SP)
        if vi>0: allok=False
        print(f"k={k} span<={SP}: #true-wide(prim)={cnt}  maxp0={float(mx):.6f}  "
              f"cap={float(CAPS[k]):.6f}  margin={float(CAPS[k]-mx):.6f}  viol={vi}")
        print(f"   argmax E={am}")
        for (E,v) in sorted(near,key=lambda z:-z[1])[:4]:
            print(f"   near: p0={float(v):.6f} margin={float(CAPS[k]-v):+.6f} E={E}")
        print()
    print("="*78)
    print("RESULT:", "0 violations across ALL exhaustive true-wide windows" if allok
          else "VIOLATION FOUND (see above)")
    print("="*78)

if __name__=="__main__":
    main()
