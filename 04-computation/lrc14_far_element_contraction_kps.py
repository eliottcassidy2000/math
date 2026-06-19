#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) ENDGAME — the UNBOUNDED->BOUNDED reduction for meas(S7).
kind-pasteur-2026-06-19-S12.

GOAL: prove that for fixed k, the sup of meas(S7(E)) over primitive non-AP
k-sets is attained at a BOUNDED set (max element <= explicit B(k)).

STRATEGY (far-element contraction): take a bounded base B0 of k-1 elements and
grow one element w. Track meas(S7(B0 u {w})) as w -> infinity. By Weyl
equidistribution, frac(w x) decorrelates from the bounded orbit, so meas(S7)
should converge to a limiting value that is the meas of the (k-1)-set under an
INDEPENDENT extra sector-hit. We measure the convergence rate and find B(k).

Definitions match lrc14_empty_sector_distribution_kps.py:
  meas(S7)(E) = measure{x in [0,1): orbit {frac(e x): e in E} hits ALL of sectors 1..6}.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def meas_S7(E):
    """exact meas(S7)(E) = P(N=0) as a Fraction."""
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

def is_AP(E):
    E = sorted(set(E))
    if len(E) < 2: return True
    d = E[1]-E[0]
    return all(E[i+1]-E[i] == d for i in range(len(E)-1))

if __name__ == "__main__":
    print("=== FAR-ELEMENT CONTRACTION: meas(S7) as one element grows ===\n")
    # For each k, take bounded base of (k-1) elements, grow the k-th element w.
    for k in [8, 9, 10]:
        cap = {8:0.38153, 9:0.49426, 10:0.6044}[k]
        base = list(range(k-1))   # consec base 0..k-2 (bounded)
        Cfull = list(range(k))
        m_consec = float(meas_S7(Cfull))
        print(f"--- k={k}  cap_{k}={cap}  consec meas={m_consec:.5f} ---")
        print(f"  base (k-1 consec) = {base}")
        rows = []
        wmax = max(base)
        for w in range(wmax+1, 80):
            E = base + [w]
            if not primitive(E): continue
            m = float(meas_S7(E))
            rows.append((w, m, is_AP(E)))
        # print the trend
        mvals = [m for (w,m,ap) in rows]
        print(f"  w=k..79: meas range [{min(mvals):.5f}, {max(mvals):.5f}]")
        # show first few and any that exceed consec
        for (w,m,ap) in rows[:12]:
            mark = "  <- AP" if ap else ("  <-- > consec!" if m > m_consec+1e-12 else "")
            print(f"    w={w:3d}  meas={m:.5f}{mark}")
        over = [(w,m) for (w,m,ap) in rows if m > m_consec+1e-9 and not ap]
        if over:
            print(f"  NON-AP w with meas>consec: {[(w,round(m,5)) for w,m in over]}")
        else:
            print(f"  NO non-AP single-far-element beats consec for w in [{wmax+1},79]")
        # tail behaviour: average of last 20
        tail = mvals[-20:]
        print(f"  tail meas (w=60..79) avg={sum(tail)/len(tail):.5f}  vs consec {m_consec:.5f}")
        print()
