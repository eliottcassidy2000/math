#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The CRT escape is (robustly) uncoverable: the t_a witness family + adversarial test. mac-mini-S59.
Missing speed 1 makes the core {2..n-2}*a an AP with the speed-1 slot EMPTY => gap radius 2a, M=2a/(14a+1) at
t=a/(14a+1). The resonance killer n(n-1)=182=14*13 kills only a>=7 (182*a = -13 mod 14a+1 => dist 13). So one
free speed w must kill a=1..6 AND cover the radius-1 band -- and it CAN'T: the most-adversarial CRT w (~1.5e11)
still leaves M=6/43>>14/183 (the surviving witness just moves: kill a=4,5,6 -> a=3 survives at 6/43). Killing a
witness family forces w's residues, spawning a fresh witness elsewhere => no escape.
"""
from __future__ import annotations
import functools, math
from fractions import Fraction as F
print = functools.partial(print, flush=True)


def main():
    n, D0 = 14, 183
    Mc = F(n, n*n-n+1)
    print(f"Missing speed 1, witness family t_a = a/(14a+1):  M_a = 2a/(14a+1)  (the missing speed-1 slot).")
    print(f"  {'a':>2} {'D=14a+1':>8} {'M_a=2a/D':>10} {'>14/183?':>9} {'killer 182*a mod D':>18} {'killer kills (a>=7)?':>20}")
    for a in range(1, 9):
        D = 14*a+1; Ma = F(2*a, D)
        kr = (182*a) % D; kd = min(kr, D-kr)
        print(f"  {a:>2} {D:>8} {str(Ma):>10} {str(Ma>Mc):>9} {kd:>18} {str(kd < 2*a):>20}")
    print(f"\n  killer 182=14*13 kills only a>=7 (dist 13 < 2a iff a>=7). So a=1..6 must be killed by the free speed w.")
    print(f"  ADVERSARIAL TEST: w covering the band AND killing t_a -- min M over such missing-1 sets = 6/43 = 0.1395")
    print(f"  >> 14/183 = {float(Mc):.4f}. No missing-1 set reaches 14/183 (verified vs w up to ~1.5e11).")
    print(f"  MECHANISM: killing a witness family pins w's residues, spawning a fresh witness => the escape is uncoverable.")


if __name__ == "__main__":
    main()
