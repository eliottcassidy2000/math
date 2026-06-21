#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_L7_finite_convergence_kps.py   (kind-pasteur 2026-06-21, HYP-2729)

VALIDATION of the L7 limit: does the EXACT finite-f1 two-far cover
   p0(B u {f1, f2}),  f2 = round(gamma*f1)
converge to p0_inf(B, gamma) (the f1->inf limit from lrc_q108_L7_resonance_atlas_kps.py)
at the THM-546 rate O(1/f1) ?  If yes, then
   L7 bound = [p0_inf <= 0.25 << cap, atlas+decay]  +  [finite comb correction <=(6/49)V/f, THM-546 PROVED]
is a COMPLETE argument (same shape as the single-far closure, one dimension up).
EXACT rational arithmetic.
"""
import os, importlib.util
from fractions import Fraction as Fr
from math import gcd

_d = os.path.dirname(__file__)
def _load(name, path):
    s = importlib.util.spec_from_file_location(name, os.path.join(_d, path))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
rcm = _load("rcm", "lrc_q108_relation_code_mds_kps.py")          # measS7 (exact, general int set)
atl = _load("atl", "lrc_q108_L7_resonance_atlas_kps.py")          # p0_inf, P2_decorrelated
measS7 = rcm.measS7

def main():
    B = [0,2,4,6,8,10,12]   # worst balanced base (k=9), 7 base + 2 far
    print("="*80)
    print("L7 finite->limit convergence  (base k=9 [0,2..12], cap_9=0.49426)")
    print("="*80)
    for p, q, name in [(2,1,"gamma=2/1 (top resonance)"), (5,3,"gamma=5/3 (2nd peak)"),
                       (89,61,"gamma~1.459 (non-resonant proxy)")]:
        pinf = float(atl.p0_inf(B, p, q))
        print(f"\n--- {name}: p0_inf = {pinf:.6f} ---")
        print(f"   {'f1':>5} {'f2':>6} {'p0(finite)':>12} {'p0-p0_inf':>12} {'(p0-pinf)*f1':>14}")
        gamma = Fr(p, q)
        for f1 in (20, 40, 80, 160, 320):
            f2 = int(gamma*f1)             # exact for q | f1; choose f1 multiple of q
            f1m = (f1 // q) * q
            if f1m == 0: f1m = q
            f2m = int(gamma*f1m)
            E = sorted(set(B) | {f1m, f2m})
            if len(E) != len(B)+2:         # avoid collision with base
                continue
            val = float(measS7(E))
            d = val - pinf
            print(f"   {f1m:>5} {f2m:>6} {val:>12.6f} {d:>+12.6f} {d*f1m:>+14.4f}")
    print("\nINTERPRETATION: if (p0-p0_inf)*f1 stays BOUNDED, the correction is O(1/f1) =")
    print("the single-far comb rate (6/49)V/f1 (THM-546 PROVED) => for f1 > W* the finite")
    print("value is within margin of p0_inf << cap; f1 <= W* is the finite window. L7 closes.")
    print("\nDONE.")

if __name__ == "__main__":
    main()
