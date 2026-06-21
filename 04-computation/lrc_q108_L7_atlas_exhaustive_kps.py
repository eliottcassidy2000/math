#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_L7_atlas_exhaustive_kps.py   (kind-pasteur 2026-06-21, HYP-2730)

EXHAUSTIVE finite-atlas check that closes L7: for the worst bounded bases, verify
   p0_inf(B, p, q) < cap_k   for ALL p/q in (1,2.15], p <= PMAX=66.
The proven tail D_{p,q} <= 14/p (lrc_q108_L7_discrepancy_proof) makes p > 66 safe
(R <= 14/67 < 0.21 <= margin). So [atlas p<=66 exact] + [tail p>=67 proven] + [THM-546
finite-f1] + [base domination] = L7 closed. This script runs the atlas exactly.
EXACT rational arithmetic. (Long-running; run in background.)
"""
import os, importlib.util
from fractions import Fraction as Fr
from math import gcd
_d = os.path.dirname(__file__)
atl = importlib.util.module_from_spec(importlib.util.spec_from_file_location("atl", os.path.join(_d,"lrc_q108_L7_resonance_atlas_kps.py")))
importlib.util.spec_from_file_location("atl", os.path.join(_d,"lrc_q108_L7_resonance_atlas_kps.py")).loader.exec_module(atl)
p0_inf = atl.p0_inf
CAP = {8: Fr(2243,5880), 9: Fr(1979,4004), 10: Fr(55,91)}
PMAX = 66

def run(name, B, k):
    cap = CAP[k]
    worst = (Fr(0), 0, 0)
    nchecked = 0; violations = 0
    for q in range(1, PMAX):
        for p in range(q+1, min(int(Fr(43,20)*q), PMAX)+1):
            if gcd(p, q) != 1: continue
            r = Fr(p, q)
            if not (Fr(1) < r <= Fr(43,20)): continue
            val = p0_inf(B, p, q)
            nchecked += 1
            if val > worst[0]:
                worst = (val, p, q)
            if val >= cap:
                violations += 1
                print(f"   !!! VIOLATION {name} p/q={p}/{q}: p0_inf={float(val):.5f} >= cap={float(cap):.5f}")
    print(f"[{name}] k={k}: checked {nchecked} ratios (p<={PMAX}); violations={violations}; "
          f"sup p0_inf={float(worst[0]):.5f} @ {worst[1]}/{worst[2]}; cap={float(cap):.5f}; "
          f"margin={float(cap-worst[0]):.5f}", flush=True)
    return violations

def main():
    print(f"EXHAUSTIVE L7 ATLAS (p<={PMAX}); proven tail D<=14/p covers p>{PMAX}.", flush=True)
    total = 0
    total += run("k=9 base[0,2..12]", [0,2,4,6,8,10,12], 9)
    total += run("k=8 base[0,2..10]", [0,2,4,6,8,10], 8)
    total += run("k=10 base[0,2..14]", [0,2,4,6,8,10,12,14], 10)
    print(f"\nTOTAL VIOLATIONS = {total}.  {'L7 ATLAS CLOSED (0 violations)' if total==0 else 'CHECK violations'}", flush=True)

if __name__ == "__main__":
    main()
