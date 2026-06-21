#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THREAD 2 -- THE W-CELL WEAK-MAJORIZATION CERTIFICATE (the joint coupling, made precise).

Discovery (lrc14_wholecircle_coupling_kpswf7.py):
  Sort each config's 6 cell-coverages W_a desc. consec's sorted W-vector
  WEAKLY MAJORIZES (dominates) every rival's at k=8 (318/318), and all but ONE
  rival at k=9 (686/687) and k=10 (398/399).

Weak majorization u >=_w v (u,v desc-sorted, length 6):
   sum_{i<=t} u_i >= sum_{i<=t} v_i  for ALL t=1..6.
At t=6 this is sum u >= sum v, i.e. measS7(consec) >= measS7(E) -- THE WALL.
So weak-majorization of the W-cell vector IMPLIES the wall (take t=6).
=> if weak-maj holds for ALL rivals, the wall is PROVED.

This script:
  (1) Re-verify weak-maj exhaustively, EXACT, larger boxes; pin down the exceptions.
  (2) For each exception E*: is it actually a wall violation (measS7(E*)>consec)?
      If measS7(E*) < consec but weak-maj fails at some t<6, then weak-maj is
      STRICTLY STRONGER than the wall and the exception is harmless to the wall
      but breaks THIS certificate -- report the precise t and cell where it fails.
  (3) Test the SHARPER ordering: is consec's sorted W-vector >= rival's
      ENTRYWISE (sorted-pointwise dominance, the strongest)? That would be a
      one-line coupling. Count entrywise failures.
  (4) Identify the structure of the exception(s): are they near-consec
      (drop-one-add-one) shapes? characterize them.
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

def Wvec(E): return sorted([W_a(E,a) for a in range(1,7)], reverse=True)
def measS7(E): return sum(W_a(E,a) for a in range(1,7))
def full_residue(E): return set(e % 7 for e in E) == set(range(7))
def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return tuple(range(k))

def weak_maj(u, v):
    """u weakly majorizes v? both desc-sorted length 6. Return (ok, first_fail_t)."""
    pu=pv=F(0)
    for t in range(6):
        pu+=u[t]; pv+=v[t]
        if pu < pv: return False, t+1
    return True, None

def entrywise(u,v):
    return all(u[i] >= v[i] for i in range(6))

def stratum(k, box):
    out=[]
    for combo in itertools.combinations(range(1,box+1), k-1):
        E=(0,)+combo
        if full_residue(E) and primitive(E): out.append(E)
    return out

if __name__=="__main__":
    print("="*88)
    print("W-CELL WEAK-MAJORIZATION CERTIFICATE for the LAYER-3 wall")
    print("="*88)
    for k, box in [(8,16),(9,15),(10,14),(11,13)]:
        print(f"\n{'#'*84}\n# k={k}, box=[0,{box}]\n{'#'*84}")
        C=consec(k); WC=Wvec(C); mC=measS7(C)
        print(f"consec sorted W-vec = {[round(float(x),5) for x in WC]}  measS7={float(mC):.6f}")
        S=stratum(k,box)
        print(f"stratum size = {len(S)}")
        wm_ok=0; wm_fail=0; ent_ok=0; ent_fail=0
        wall_viol=0
        fail_details=[]
        ent_fail_examples=[]
        for E in S:
            if E==C: continue
            WE=Wvec(E); mE=measS7(E)
            if mE>mC: wall_viol+=1
            ok,t=weak_maj(WC,WE)
            if ok: wm_ok+=1
            else:
                wm_fail+=1
                fail_details.append((E,t,float(mE),[round(float(x),5) for x in WE]))
            if entrywise(WC,WE): ent_ok+=1
            else:
                ent_fail+=1
                if len(ent_fail_examples)<6: ent_fail_examples.append((E,[round(float(x),5) for x in WE]))
        print(f"  WALL violations (measS7(E)>consec): {wall_viol}   <-- must be 0")
        print(f"  WEAK-MAJ holds: {wm_ok}/{len(S)-1}   FAILS: {wm_fail}")
        for E,t,mE,WE in fail_details[:8]:
            print(f"     WEAK-MAJ FAIL at t={t}: {E} measS7={mE:.6f} (consec {float(mC):.6f}) Wvec={WE}")
        print(f"  ENTRYWISE (sorted-pointwise) holds: {ent_ok}/{len(S)-1}   FAILS: {ent_fail}")
        for E,WE in ent_fail_examples:
            print(f"     ENTRYWISE FAIL: {E} Wvec={WE}")
