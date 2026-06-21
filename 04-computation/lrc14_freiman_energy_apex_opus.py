#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE 3 -- THE APEX COINCIDENCE (opus 2026-06-21).

Integer additive ENERGY is the strongest additive predictor of measS7 (Spearman
~+0.65..+0.70), and at the EXTREME, consec is simultaneously:
    argmax measS7 == argmax additive_energy == argmin doubling.
This is the falsifiable apex hypothesis. The interior ranking disagrees (energy is
not rank-faithful), so a correlation cannot prove Layer 3. BUT the APEX coincidence,
if it holds at all k, gives a clean BRIDGE:

   HYP-FreimanApex:  over the full-residue primitive stratum, the measS7-maximizer
   EQUALS the additive-energy-maximizer, and additive-energy is maximized (among
   k-sets of integers, with the residue-cover constraint) UNIQUELY by the AP -- a
   THEOREM (additive energy is maximized by APs; e(A) <= ... with equality iff AP,
   classical inverse theory). Hence measS7-max => AP => consec.

DECISIVE TESTS:
  (A) k=9, k=10: does argmax measS7 == argmax additive_energy == consec persist?
  (B) Is additive energy maximized UNIQUELY by consec over the full-residue stratum?
      (i.e. is the apex coincidence really a *coincidence of two maxima*, both at AP?)
  (C) The KEY potential REFUTATION: find a shape with HIGHER additive energy than
      consec but LOWER measS7 (or vice versa) -- if the two maxima ever separate,
      the bridge is broken. Report the closest competitors on BOTH functionals.
  (D) Within the stratum, restrict to APs (gcd-normalized translates/dilates):
      is consec the measS7-max AMONG APs, and is it the energy-max among APs?
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
def add_energy(E):
    c=defaultdict(int)
    for a in E:
        for b in E: c[a+b]+=1
    return sum(v*v for v in c.values())
def doubling(E): return F(len({a+b for a in E for b in E}), len(E))

def residues(E): return frozenset(e%7 for e in E)
def is_full_residue(E): return residues(E)==frozenset(range(7))
def primitive(E): return reduce(gcd,[abs(e) for e in E if e!=0],0)==1
def consec(k): return list(range(k))
def is_AP(E):
    s=sorted(set(E))
    if len(s)<2: return True
    d=s[1]-s[0]
    return all(s[i+1]-s[i]==d for i in range(len(s)-1))

if __name__=="__main__":
    print("#"*78); print("# APEX COINCIDENCE: measS7-max == energy-max == AP ?"); print("#"*78)
    for k,W in [(8,13),(9,12),(9,13),(10,13)]:
        C=consec(k)
        bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
        full=[list(E) for E in bank if primitive(E) and is_full_residue(E)]
        if not full:
            print(f"\nk={k} span<= {W}: EMPTY stratum"); continue
        rows=[(E,measS7(E),add_energy(E),doubling(E),is_AP(E)) for E in full]
        bm=max(rows,key=lambda r:r[1])   # max measS7
        be=max(rows,key=lambda r:r[2])   # max energy
        print(f"\n{'='*70}\nk={k} span<= {W}: {len(full)} shapes  (consec={C})")
        print(f"  measS7-max : {bm[0]}  measS7={float(bm[1]):.6f} energy={bm[2]} doub={float(bm[3]):.3f} AP={bm[4]}")
        print(f"  energy-max : {be[0]}  measS7={float(be[1]):.6f} energy={be[2]} doub={float(be[3]):.3f} AP={be[4]}")
        print(f"  >>> APEX COINCIDENCE (measS7-max == energy-max)? {tuple(bm[0])==tuple(be[0])}")
        print(f"  >>> both == consec? measS7max:{tuple(bm[0])==tuple(C)}  energymax:{tuple(be[0])==tuple(C)}")
        # uniqueness of energy max
        maxe=be[2]; at_maxe=[r for r in rows if r[2]==maxe]
        print(f"  energy-max value={maxe}; #achieving it={len(at_maxe)}; all AP? {all(r[4] for r in at_maxe)}")
        for r in at_maxe[:6]: print(f"      {r[0]}  measS7={float(r[1]):.6f} AP={r[4]}")
        # (C) closest competitors: 2nd by measS7 and 2nd by energy
        s_meas=sorted(rows,key=lambda r:-r[1]); s_en=sorted(rows,key=lambda r:-r[2])
        print(f"  2nd measS7: {s_meas[1][0]} measS7={float(s_meas[1][1]):.6f} energy={s_meas[1][2]} AP={s_meas[1][4]}")
        print(f"  2nd energy: {s_en[1][0]} measS7={float(s_en[1][1]):.6f} energy={s_en[1][2]} AP={s_en[1][4]}")
        # (D) AP-only sub-stratum
        aps=[r for r in rows if r[4]]
        if aps:
            apm=max(aps,key=lambda r:r[1]); ape=max(aps,key=lambda r:r[2])
            print(f"  AMONG APs ({len(aps)}): measS7-max={apm[0]}(={float(apm[1]):.6f}) "
                  f"energy-max={ape[0]}(en={ape[2]})  both consec? "
                  f"{tuple(apm[0])==tuple(C) and tuple(ape[0])==tuple(C)}")
