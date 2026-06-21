#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THREAD 2 -- the EXACT p0-argmax pattern across k, and WHY k=12 breaks consec-extremality.

Established (lrc14_canonical_p0_k12_kpswf7):
  CANONICAL p0=measS7_full over [0,1]:
   - k=8,9,10,11 : consec is the UNIQUE argmax on full-residue/primitive/span<=14.
   - k=12        : consec BEATEN by E*=(0..10,12) and 3 others, but ALL << cap_12.

This script:
  (1) For each k=8..14, find the TOP-5 p0-shapes on full-residue/primitive/span<=14,
      EXACT. Print the residue multiset & which residue is "doubled" beyond consec.
      Conjecture: the argmax always DOUBLES residues 0 and (the small ones) -- consec
      doubles {0,1,2,3} (= residues of 0,7,1,8,...). At k=12 the optimum prefers
      doubling residue 5 (via e=12) over residue 4. Characterize the crossover.
  (2) Confirm the wall p0 <= cap_k holds for ALL with margin; report the worst
      (max p0, min margin) per k. cap_k EXACT from canon:
        cap_8=2243/5880, cap_9, ..., cap_12=6/7.  (THM-538: each >= (k-6)/7.)
  (3) The DECISIVE strategic readout: is the FINITE layer's needed statement
      "consec maximizes p0" (FALSE at k>=12) or "p0 <= cap_k" (TRUE, margin)?
"""
import itertools, math, sys
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)
P=7

def measS7_full(E):
    E = sorted(set(int(e) for e in E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        d = P*abs(e)
        for m in range(0, d+1): bps.add(F(m, d))
    bps = sorted(b for b in bps if 0 <= b <= 1); tot = F(0)
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi <= lo: continue
        mid = (lo+hi)/2
        if len({int(((e*mid)%1)*P) for e in E})==P: tot += hi-lo
    return tot

def full_residue(E): return set(e % P for e in E)==set(range(P))
def primitive(E): return reduce(gcd,[e for e in E if e!=0],0)==1
def consec(k): return tuple(range(k))

# EXACT caps from canon THM-538 (cap_8=2243/5880; cap_12=6/7). Fill known + (k-6)/7 floor.
CAP = {8:F(2243,5880), 12:F(6,7)}
def cap_floor(k): return F(k-6,7)  # THM-538 lower bound; use as conservative wall when exact cap unknown

def stratum(k, box):
    out=[]
    for combo in itertools.combinations(range(1,box+1), k-1):
        E=(0,)+combo
        if full_residue(E) and primitive(E): out.append(E)
    return out

if __name__=="__main__":
    print("="*92)
    print("EXACT p0-argmax pattern + wall margin across k (full-residue/primitive/span<=14)")
    print("="*92)
    for k in range(8,15):
        S=stratum(k,14)
        if not S:
            print(f"k={k}: empty stratum at span<=14"); continue
        scored=sorted(((measS7_full(E),E) for E in S), key=lambda z:-z[0])
        pC=measS7_full(consec(k))
        top=scored[0]
        is_consec = top[1]==consec(k)
        capk = CAP.get(k, None)
        floor = cap_floor(k)
        # wall margin against floor (conservative) and exact cap if known
        worst_p0 = top[0]
        marg_floor = floor - worst_p0
        print(f"\nk={k}: |strat|={len(S)}  consec p0={float(pC):.6f}  argmax={'CONSEC' if is_consec else 'NOT consec'}")
        for v,E in scored[:5]:
            resmult = sorted(e%P for e in E)
            # doubled residues beyond a single copy:
            from collections import Counter
            cnt=Counter(e%P for e in E)
            doubled = sorted([r for r in cnt if cnt[r]>=2])
            tag = " <-- CONSEC" if E==consec(k) else ""
            print(f"    p0={float(v):.6f}={v}  E={E}  doubled_res={doubled}{tag}")
        capstr = (f"cap_k={capk}={float(capk):.5f}" if capk else f"cap_k=? (floor (k-6)/7={float(floor):.5f})")
        print(f"    {capstr}  max p0={float(worst_p0):.6f}  margin-to-floor={float(marg_floor):+.6f}  WALL(vs floor) {'OK' if worst_p0<=floor else 'VIOLATED'}")
        if capk:
            print(f"    margin-to-EXACTcap = {float(capk-worst_p0):+.6f}  WALL(vs cap) {'OK' if worst_p0<=capk else 'VIOLATED'}")
