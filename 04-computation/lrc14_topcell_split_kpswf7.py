#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THREAD 2 -- THE TOP-CELL SPLIT (the joint coupling, final form).

Robust finding (lrc14_matched_coupling_kpswf7.py): the sorted W-cell partial-sum
deficit D_t = sum_{i<=t}(WC_(i) - WE_(i)) is >= 0 for ALL t>=2 over EVERY rival
(k=8..11). It can only go negative at t=1 (a single rival cell beats consec's
single best cell). At t=6, D_6 = measS7(consec)-measS7(E) >= 0 = THE WALL.

So define the SPLIT:
   measS7(consec) - measS7(E)
     = [WC_(1) - WE_(1)]            (top-cell term, can be NEGATIVE)
     + sum_{t=2..6}[WC_(t)-WE_(t)]  (tail term, must overcompensate)

The clean structural claim to PROVE/refute:
   (P1)  TAIL DOMINANCE:  for all t in {2,3,4,5,6},
            sum_{i=1..t} WC_(i)  >=  sum_{i=1..t} WE_(i).
         (cumulative sorted dominance from t>=2; equiv. consec's "all but the
          single largest cell" already wins at every depth.)  <-- if universal,
          this is a JOINT 5-inequality certificate STRICTLY STRONGER than the wall
          minus the single hardest cell. The wall (t=6) is one of them.

This script:
 (1) Verify P1 exhaustively at k=8..12, LARGE boxes, EXACT. Report any t<6 failure.
 (2) Characterize the t=1 beaters: what is W_(1)(E) - W_(1)(consec) and which E?
     Conjecture: the ONLY t=1 beaters are the "shift-one" shapes
     consec\{k-1} u {k+1}, ... i.e. residue-doubling near the top. Confirm.
 (3) The MIN-CUT / bottleneck reading: is consec the unique maximizer of the
     SMALLEST cell min_a W_a (the bottleneck)? (max-min). If consec maxes the
     bottleneck AND the tail-cumulative, that is a max-min + majorization combo.
"""
import itertools, math, sys
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def W_a(E, a):
    E = sorted(set(E))
    lo = F(a, 7) - F(1, 14); hi = F(a, 7) + F(1, 14)
    bps = {lo, hi}
    for e in E:
        if e == 0: continue
        d = 7 * abs(e)
        j0 = math.floor(lo * d); j1 = math.ceil(hi * d)
        for j in range(j0 - 1, j1 + 2):
            x = F(j, d)
            if lo <= x <= hi: bps.add(x)
    bps = sorted(bps); tot = F(0)
    for l, h in zip(bps, bps[1:]):
        if h <= l: continue
        xm = (l + h) / 2; hit = set()
        for e in E:
            v = e * xm; v = v - (v.numerator // v.denominator)
            hit.add((v.numerator * 7) // v.denominator)
        if len(hit) == 7: tot += h - l
    return tot

def Wraw(E): return [W_a(E,a) for a in range(1,7)]
def measS7(E): return sum(Wraw(E))
def full_residue(E): return set(e % 7 for e in E) == set(range(7))
def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return tuple(range(k))

def stratum(k, box):
    out=[]
    for combo in itertools.combinations(range(1,box+1), k-1):
        E=(0,)+combo
        if full_residue(E) and primitive(E): out.append(E)
    return out

if __name__=="__main__":
    print("="*88)
    print("TOP-CELL SPLIT: cumulative sorted W-cell dominance for t>=2 (joint certificate)")
    print("="*88)
    for k, box in [(8,18),(9,16),(10,15),(11,14),(12,13)]:
        C=consec(k); WCraw=Wraw(C); mC=measS7(C); WCs=sorted(WCraw,reverse=True)
        S=stratum(k,box)
        print(f"\n# k={k} box=[0,{box}] | consec sorted W={[round(float(x),5) for x in WCs]} measS7={float(mC):.6f} | stratum={len(S)}")
        # P1: tail-cumulative dominance for t=2..6
        tail_fail={2:0,3:0,4:0,5:0,6:0}
        tail_fail_ex={2:None,3:None,4:None,5:None,6:None}
        # bottleneck (min cell) max-min
        bott_C=min(WCraw); bott_beaten=0; bott_ex=None
        # t=1 beaters
        t1_beaters=[]
        for E in S:
            if E==C: continue
            WEraw=Wraw(E); WEs=sorted(WEraw,reverse=True)
            pu=pv=F(0)
            for t in range(6):
                pu+=WCs[t]; pv+=WEs[t]
                if t>=1 and pu<pv:
                    tail_fail[t+1]+=1
                    if tail_fail_ex[t+1] is None: tail_fail_ex[t+1]=(E,float(pu),float(pv))
            # t=1 beaters
            if WEs[0]>WCs[0]:
                t1_beaters.append((E, float(WEs[0]-WCs[0])))
            # bottleneck
            if min(WEraw)>bott_C:
                bott_beaten+=1
                if bott_ex is None: bott_ex=(E,float(min(WEraw)),float(bott_C))
        print(f"  TAIL cumulative dominance failures (t=2..6): {tail_fail}")
        for t in (2,3,4,5,6):
            if tail_fail_ex[t]: print(f"     t={t} FAIL example: {tail_fail_ex[t]}")
        print(f"  bottleneck min_a W_a(consec)={float(bott_C):.6f}; #rivals beating bottleneck={bott_beaten}",
              (f" ex={bott_ex}" if bott_ex else ""))
        print(f"  #t=1 beaters (single cell beats consec top cell) = {len(t1_beaters)}")
        for E,d in sorted(t1_beaters,key=lambda z:-z[1])[:6]:
            # characterize: is E = consec with the largest element shifted up by gap?
            diff_from_consec = (set(E)-set(C), set(C)-set(E))
            print(f"     t1-beater {E} delta_topcell=+{d:.5f} added={sorted(diff_from_consec[0])} removed={sorted(diff_from_consec[1])}")
