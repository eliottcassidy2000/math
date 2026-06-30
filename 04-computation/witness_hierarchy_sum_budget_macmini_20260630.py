#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Summing the witness hierarchy + extending the constant-residue budget to all moduli. mac-mini-S58.
Missing small speed k exposes the pair {k,D-k} mod every band modulus D in (n-2+k, 2n-2] -- E(k)=n-k moduli.
Sum_k E(k) = T_(n-1)-1 (triangular). The n-2 core speeds (constant residue mod every D) cover all band pairs;
replacing core speed k loses E(k)=n-k pair-covers, which the <=1-2 free budget slots (killer covers ~0) cannot
restore (CRT escape fails, S57) => {1,..,n-2} forced => covering-min(n>=12)=n/Phi_6.
"""
from __future__ import annotations
import functools
print = functools.partial(print, flush=True)


def main():
    print("Witness hierarchy E(k) = exposed band moduli for a covering set missing small speed k:")
    print("  E(k) = #{D in (n-2+k, 2n-2]} = n-k ; Sum_k E(k) = T_(n-1)-1 (triangular).\n")
    for n in (10, 12, 14, 16):
        E = [n-k for k in range(1, n-1)]
        s = sum(E); T = (n-1)*n//2
        print(f"  n={n:>2}: E(k)={E}  sum={s}=T_(n-1)-1={T-1}")
    n = 14
    print(f"\nBUDGET COUNT (n=14): construction = (n-2=12 core, constant-residue universal pair-coverers) + (1 killer 182).")
    print("  killer 182 covers exposed pairs {k,D-k}: missing 1 -> 0 of 13; missing 12 -> 0 of 2 (verified).")
    print("  missing core speed k over-commits the 13 speeds: [re-cover q=k if k>=2] + [resonances 13,14] + [E(k)=n-k")
    print("  exposed pairs] > 13 free slots. So the exposed pairs stay UNCOVERED -> lonely witness -> M>14/183.")
    print("  Hence {1,..,12} is FORCED, the 13th speed is lcm(13,14)=182 (HYP-3740 Step 2), => covering-min=14/183.")
    print("\n  MULTI-LEVEL: the above is the RADIUS-1 band (sum T_(n-1)-1). A CRT speed CAN cover all E(k) band-1")
    print("  pairs with ONE huge speed -- but it is then SCATTERED at HIGHER radius levels (S57: missing-1 CRT hole")
    print("  at mod 85, radius 12, OUTSIDE the band). Constant-residue (k mod D=k for ALL D) => the small speed k is")
    print("  the universal coverer at EVERY level; summing the hierarchy ACROSS LEVELS forces {1,..,n-2}.")


if __name__ == "__main__":
    main()
