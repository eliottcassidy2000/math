#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
REFUTES opus-S62 (HYP-4068 / MISTAKE-101): the f=2 tightening gap does NOT grow with u_max.

S62 claimed gap(U)=min_{odd w1,w2}(M(2U u {w1,w2})-1/14) rises with u_max (over RANDOM even parts) =>
u_max bounded => mac-mini's finite confinement check terminates. That was a RANDOM-SAMPLING ARTIFACT:
random 11-subsets are incommensurate (growing gap), but the near-extremal families are COMMENSURATE
dilated APs, which random sampling misses.

EXACT FACT (this script): U = c*{1..11} gives M(2c*{1..11} u {w1,w2}) = 1/12 EXACTLY, so gap = 1/12-1/14
= 1/84 CONSTANT for all c, all PRIMITIVE, u_max = 22c -> infinity. So u_max is NOT bounded via gap-growth;
confinement (if it holds) is a UNIFORM POSITIVE gap with UNBOUNDED (dilated-AP) extremizers -- "bound
v_max(U)" is NOT the right target.

Why M=1/12: M(2c*{1..11}) = M({1..11}) = 1/12 (scale-invariant AP), and two odd tighteners fail to reduce
it (they cannot cover the point achieving 1/12) => tighteners USELESS, M(S)=1/12 for these families.
"""
import sys
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def norm(x):
    x = x - int(x)
    if x < 0: x += 1
    return min(x, 1 - x)
def exact_M(S):
    S = sorted(set(S)); cands = set()
    for v in S:
        for k in range(v): cands.add(Fr(2*k+1, 2*v))
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for den in (S[i]+S[j], abs(S[i]-S[j])):
                if den:
                    for s in range(den): cands.add(Fr(s, den))
    b = Fr(0)
    for t in cands:
        val = min(norm(v*t) for v in S)
        if val > b: b = val
    return b
def prim(S): return reduce(gcd, S) == 1

print("="*92)
print(" REFUTATION: dilated APs c*{1..11} have CONSTANT gap 1/84, u_max -> infinity (S62 gap-growth is false)")
print("="*92)
print(f"  1/84 = {float(Fr(1,84)):.6f}   (= 1/12 - 1/14, M({{1..11}})=1/12)")
print(f"  {'c':>3} {'u_max':>6} {'tighteners':>12} {'M(S) exact':>12} {'gap=M-1/14':>12} {'primitive':>10}")
for c, ws in [(1,(1,7)), (2,(9,11)), (3,(3,11)), (4,(5,29)), (5,(1,3))]:
    E = [2*c*j for j in range(1,12)]
    S = sorted(E + list(ws))
    M = exact_M(S); gap = M - Fr(1,14)
    print(f"  {c:>3} {max(E):>6} {str(ws):>12} {str(M):>12} {str(gap):>12} {str(prim(S)):>10}")
print("\n  Also confirm M(even part alone) = M(2c*{1..11}) = M({1..11}) = 1/12 (tighteners don't reduce it):")
for c in [1,2,3]:
    E = [2*c*j for j in range(1,12)]
    print(f"    c={c}: M(2*{c}*{{1..11}}) = {str(exact_M(E))}   (=1/12? {exact_M(E)==Fr(1,12)})")
print("\n"+"="*92)
print(" CONCLUSION")
print("="*92)
print("  * gap(c*{1..11}) = 1/84 CONSTANT => u_max NOT bounded by gap-growth (S62/HYP-4068 REFUTED, MISTAKE-101).")
print("  * The min-gap extremizers are dilated APs at ALL scales (unbounded u_max).")
print("  * Confinement (if it holds) = a UNIFORM POSITIVE gap, NOT a finite u_max check. 'Bound v_max(U)'")
print("    is the WRONG target; the right object is inf_U gap(U) > 0, extremized by unbounded dilated APs.")
print("DONE.")
