#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
BELOW-V* WINDOW WITNESS GENERATOR for the LRCWitnessCert shapes (kps HYP-3958/3959).

For each certificate shape (P, offs) with tail threshold V*, and each V in the finite window
[max(offs)+1, V*), the configuration S_V = P u {V - o : o in offs} is a fixed <=13-runner set.
We find an exact rational witness t = a/b with  min_{v in S_V} ||v t|| >= 1/14, verified by the
CLOSED Nat criterion  b <= 14 * min((v*a) % b, b - (v*a) % b)  for every v -- the same check the
Lean side re-runs by kernel `decide`.  Output: Lean rows + per-branch tactic lines.

Witness search = exact max-min over the complete breakpoint set (kinks k/(2v), crossings
k/(v+w), k/(w-v)), as in covering_min_court_verification (MISTAKE-086-safe).

opus-2026-07-02-S35.
"""
import sys
from fractions import Fraction as Fr
from math import gcd
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

H = Fr(1, 14)

def dist_frac(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def nat_check(S, a, b):
    """the exact Lean-side criterion: b <= 14*min((v*a)%b, b-(v*a)%b) for all v in S."""
    for v in S:
        m = (v * a) % b
        if 14 * min(m, b - m) < b:
            return False
    return True

def best_witness(S):
    """exact max-min witness over the complete breakpoint set; returns Fraction."""
    dens = set()
    for v in S:
        dens.add(v); dens.add(2 * v)
    L = sorted(S)
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            dens.add(L[i] + L[j])
            if L[j] > L[i]: dens.add(L[j] - L[i])
    dens.discard(0)
    best, argt = Fr(-1), None
    for d in sorted(dens):
        for k in range(0, d + 1):
            t = Fr(k, d)
            m = min(dist_frac(v * t) for v in S)
            if m > best:
                best, argt = m, t
    return best, argt

SHAPES = [
    ("AP", [], [0,1,2,3,4,5,6,7,8,9,10,11,12], 15),
    ("3",  [1,2,3], [0,1,2,3,5,8,11,14,17,20], 47),
    ("4",  [2,3,4,5,6], [0,2,4,6,8,10,12,14], 40),
]

print("shape | V | witness a/b | max-min | Nat-check")
rows = {}
allok = True
for name, P, offs, Vstar in SHAPES:
    lo = max(offs) + 1
    out = []
    for V in range(lo, Vstar):
        S = sorted(set(P) | {V - o for o in offs})
        mm, t = best_witness(S)
        a, b = t.numerator, t.denominator
        ok = mm >= H and nat_check(S, a, b)
        allok = allok and ok
        print(f"  {name:>3} | {V:>3} | {a}/{b} | {float(mm):.6f} | {'OK' if ok else '*** FAIL ***'}")
        out.append((V, a, b))
    rows[name] = (lo, out)
print(f"\nALL WINDOW ROWS PASS: {allok}")

if allok:
    print("\n" + "=" * 60 + "\n LEAN EMISSION\n" + "=" * 60)
    for name, P, offs, Vstar in SHAPES:
        lo, out = rows[name]
        lst = ", ".join(f"({V}, {a}, {b})" for V, a, b in out)
        print(f"\ndef window{name} : List (ℤ × ℕ × ℕ) := [{lst}]")
        print(f"-- window{name}: V ∈ [{lo}, {Vstar})  ({len(out)} rows)")
        for V, a, b in out:
            print(f"  -- V={V}: exact lonely_of_row (a := {a}) (b := {b}) (by decide)")
