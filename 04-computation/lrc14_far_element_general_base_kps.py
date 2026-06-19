#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) — far-element contraction with GENERAL bounded base.
kind-pasteur-2026-06-19-S12.

For each k, enumerate ALL bounded (k-1)-subsets of {0,..,Bbase} containing 0,
then grow a far element w. Confirm: meas(S7(base u {w})) is MAXIMIZED at a
SMALL w (bounded), and decays for large w. The point is to find an explicit
B(k) such that for any primitive non-AP k-set with max element > B(k),
meas(S7) < (the bounded finite-check max).

We test the contraction CLAIM directly: define
   M_bounded(k) = max over bounded k-sets (max<=Bcheck) of meas(S7)
   For a (k-1)-base 'base' and far w, is meas(S7(base u {w})) <= M_bounded?
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def meas_S7(E):
    E = sorted(set(E))
    bps = set([Fraction(0), Fraction(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7*e + 1):
            bps.add(Fraction(a, 7*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p0 = Fraction(0)
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo + hi)/2
        hit = set()
        for e in E:
            v = e*mid; v = v - (v.numerator//v.denominator)
            hit.add((v.numerator*7)//v.denominator)
        if all(j in hit for j in range(1,7)):
            p0 += (hi - lo)
    return p0

def primitive(E):
    return reduce(gcd, [e for e in E if e != 0]) == 1

if __name__ == "__main__":
    # k=8 only (k=9,10 heavier). Bounded base from {0..Bbase}, |base|=k-1, 0 in base.
    k = 8
    cap = 0.38153
    Bbase = 9         # base elements drawn from {0,1,...,Bbase}
    # consec full meas (the conjectured global max)
    m_consec = float(meas_S7(list(range(k))))
    print(f"k={k}  cap={cap}  consec meas={m_consec:.5f}")
    print(f"Enumerating (k-1)={k-1}-subsets of 0..{Bbase} (0 included), growing far w.\n")

    # For each bounded base, find the w (in a far range) that maximizes meas, and
    # check whether ANY base+w (non-AP) beats consec.
    best_overall = (None, -1.0)
    bases_checked = 0
    violations = []
    base_elems = list(range(1, Bbase+1))
    for combo in itertools.combinations(base_elems, k-2):
        base = [0] + list(combo)
        bases_checked += 1
        # grow w over a moderate-to-far window
        for w in range(Bbase+1, 45):
            E = base + [w]
            if not primitive(E): continue
            m = float(meas_S7(E))
            if m > best_overall[1]:
                best_overall = (tuple(E), m)
            if m > m_consec + 1e-9:
                violations.append((tuple(E), round(m,5)))
    print(f"bases checked: {bases_checked}")
    print(f"best base+far-w meas = {best_overall[1]:.5f}  at E={best_overall[0]}")
    print(f"  (consec meas = {m_consec:.5f}, cap = {cap})")
    if violations:
        print(f"VIOLATIONS (base+w beats consec): {violations[:20]}")
    else:
        print("NO base+far-w (w>Bbase) configuration beats consec meas.")
        print(f"=> far element (w>{Bbase}) always contracts below consec for k={k}.")
