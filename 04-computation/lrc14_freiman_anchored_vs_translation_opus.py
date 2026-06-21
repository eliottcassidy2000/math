#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE 3 -- THE RECONCILIATION: anchored leg-profile vs translation-invariant energy.
(opus 2026-06-21)

The k=10 inversion E1=[0,2,3,4,5,6,7,8,9,10] vs E2=[0,1,2,3,4,5,6,7,9,10]:
  - IDENTICAL doubled-residue structure {0:[0,7],2:[2,9],3:[3,10]}.
  - DIFFER only in residue-1's leg: E2 has min|e|=1, E1 has min|e|=8.
  - E2 wins measS7 (anchored leg better); E1 wins additive energy (longer run).
=> additive energy is TRANSLATION-INVARIANT (rewards runs anywhere), measS7 is
   ANCHORED at 0 (rewards small min|e| per residue). They agree at the apex,
   diverge off it. The Freiman route proves the APEX but not the order.

FALSIFIABLE HYPOTHESIS (the real lever, replacing raw Freiman):
  HYP-AnchoredSchur: Define the LEG PROFILE L(E) = sorted multiset {min|e| over
  e in E with e≡r (mod 7)} for r=0..6. Among full-residue shapes, the leg profile
  is MAJORIZED-minimized by consec (legs {0,1,2,3,4,5,6}), and measS7 is
  SCHUR-CONCAVE-decreasing in L... BUT this is Layer 1/2 only. The TEST: is measS7
  monotone (Schur) in the leg profile ALONE? If yes, Layer 1 + Schur-concavity =>
  consec max with NO finite check (since consec uniquely minimizes the leg profile
  in majorization order). If NO (measS7 depends on MORE than legs), identify the
  smallest extra statistic.

TEST:
  (1) Group shapes by IDENTICAL leg profile. Within a group, does measS7 vary?
      If measS7 is CONSTANT on leg-profile classes => leg profile is a sufficient
      statistic and Schur-concavity is the whole story (HUGE).
      If it varies => legs are necessary but not sufficient; report the residual.
  (2) Confirm consec uniquely minimizes the leg profile (majorization) over stratum.
  (3) Across leg-profile classes, is measS7 Schur-monotone (smaller-majorized leg
      profile => higher measS7)? Find any violation.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def occupancy_pi(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7*ae+1): bps.add(F(m, 7*ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    pi = [F(0)]*8
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        xm = (lo+hi)/2
        hit = set()
        for e in E:
            v = e*xm; v = v - (v.numerator // v.denominator)
            hit.add((v.numerator*7)//v.denominator)
        pi[len(hit)] += hi - lo
    return pi
def measS7(E): return occupancy_pi(E)[7]

def leg_profile(E):
    byres=defaultdict(list)
    for e in E: byres[e%7].append(abs(e))
    return tuple(sorted(min(byres[r]) for r in range(7) if r in byres))

def second_profile(E):
    # 2nd-smallest magnitude per residue (the LAYER 2 statistic), -1 if singleton
    byres=defaultdict(list)
    for e in E: byres[e%7].append(abs(e))
    return tuple(sorted((sorted(byres[r])[1] if len(byres[r])>1 else -1) for r in range(7) if r in byres))

def residues(E): return frozenset(e%7 for e in E)
def is_full_residue(E): return residues(E)==frozenset(range(7))
def primitive(E): return reduce(gcd,[abs(e) for e in E if e!=0],0)==1
def consec(k): return list(range(k))

def majorizes(a,b):
    """does sorted-desc a majorize b? (a more spread). Both same length & sum."""
    if sum(a)!=sum(b): return None
    A=sorted(a,reverse=True); B=sorted(b,reverse=True)
    pa=pb=0
    for i in range(len(A)):
        pa+=A[i]; pb+=B[i]
        if pa<pb: return False
    return True

if __name__=="__main__":
    print("#"*78); print("# ANCHORED LEG-PROFILE vs measS7 -- is leg profile sufficient?"); print("#"*78)
    for k,W in [(8,13),(9,13),(10,14)]:
        C=consec(k)
        bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
        full=[list(E) for E in bank if primitive(E) and is_full_residue(E)]
        if not full:
            print(f"\nk={k} W={W}: empty"); continue
        print(f"\n{'='*70}\nk={k} span<= {W}: {len(full)} shapes  (consec legs={leg_profile(C)})")
        # (1) group by leg profile; does measS7 vary within a group?
        groups=defaultdict(list)
        for E in full: groups[leg_profile(E)].append(E)
        varying=0; constant=0; max_spread=F(0); worst=None
        for lp,Es in groups.items():
            ms=set(measS7(E) for E in Es)
            if len(ms)>1:
                varying+=1; sp=max(ms)-min(ms)
                if sp>max_spread: max_spread=sp; worst=(lp,Es)
            else: constant+=1
        print(f"  leg-profile classes: {len(groups)}  (constant measS7: {constant}, varying: {varying})")
        print(f"  => leg profile SUFFICIENT statistic? {varying==0}")
        if worst:
            lp,Es=worst
            print(f"  worst within-class measS7 spread = {float(max_spread):.6f}, leg profile {lp}, {len(Es)} shapes:")
            for E in sorted(Es,key=lambda e:-measS7(e))[:4]:
                print(f"      {E}  measS7={float(measS7(E)):.6f}  2nd-profile={second_profile(E)}")
        # (2) consec uniquely minimizes leg profile (majorization)?
        lpC=leg_profile(C)
        not_dominated=[E for E in full if majorizes(leg_profile(E),lpC) is False]
        # consec leg sum is minimal? check sums
        sums=set(sum(leg_profile(E)) for E in full)
        print(f"  consec leg-sum={sum(lpC)}; leg-sums in stratum range {min(sums)}..{max(sums)}; "
              f"consec minimal-sum? {sum(lpC)==min(sums)}")
        # (3) Schur monotone: across classes with consec's leg-sum, is consec's the
        # majorization-minimal AND measS7-max?
        same_sum=[E for E in full if sum(leg_profile(E))==sum(lpC)]
        bm_ss=max(same_sum,key=measS7)
        print(f"  among leg-sum=={sum(lpC)} shapes ({len(same_sum)}): measS7-max = "
              f"{'consec' if tuple(bm_ss)==tuple(C) else bm_ss} (={float(measS7(bm_ss)):.6f})")
        # global apex still consec?
        gm=max(full,key=measS7)
        print(f"  GLOBAL measS7-max = {'consec' if tuple(gm)==tuple(C) else gm}")
