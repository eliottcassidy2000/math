#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD B FAST finalizer: the per-base combined cutoffs W'(B) (single-far) WITHOUT the
expensive wide-p0 sweep -- only V(B), Phi(B) on BOUNDED bases (cheap).  Gives the uniform
finite-band threshold max_B W'(B) for k=8,9,10 exactly.  Also the double-far sup|C| on a
representative bounded-base set.

W'(B) := (6/49) V(B) / (cap_k - Phi(B))  is the per-base single-far cutoff: for w > W'(B),
  p0(Bu{w}) = Phi(B) + Delta_w <= Phi(B) + (6/49)V(B)/w < cap_k   (THM-546, PROVED).
The single-far residual is then the FINITE band  15 <= w <= max_B W'(B), bounded base.
"""
from __future__ import annotations
import sys, functools, itertools, math
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
from lrc14_wide_branch_ridge_codex_s47 import CAP
from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct
QVAL = {k: boundary_value_direct(tuple(range(k - 1)), 1) for k in CAP}

def orbit_breakpoints(Ep):
    Ep = sorted(set(Ep)); bp = {F(0), F(1)}
    for e in Ep:
        if e == 0: continue
        for j in range(0, 7 * e + 1): bp.add(F(j, 7 * e))
    return sorted(b for b in bp if 0 <= b < 1)
def cells_with_miss(Ep):
    Ep = [e for e in sorted(set(Ep)) if e != 0]; bp = orbit_breakpoints(Ep); out = []
    for lo, hi in zip(bp, bp[1:] + [F(1)]):
        if hi <= lo: continue
        mid = (lo + hi) / 2; hit = set(int((e * mid) % 1 * 7) for e in Ep)
        out.append((lo, hi, frozenset(set(range(1, 7)) - hit)))
    return out
def dist_p_bounded(E):
    Eps = [e for e in sorted(set(E)) if e != 0]; bp = {F(0), F(1)}
    for e in Eps:
        for j in range(0, 7 * e + 1): bp.add(F(j, 7 * e))
    bp = sorted(b for b in bp if 0 <= b < 1); p = [F(0)] * 7
    for lo, hi in zip(bp, bp[1:] + [F(1)]):
        if hi <= lo: continue
        mid = (lo + hi) / 2; hit = set(int((e * mid) % 1 * 7) for e in Eps)
        p[len(set(range(1, 7)) - hit)] += hi - lo
    return p
def Phi(Ep):
    p = dist_p_bounded(Ep); return p[0] + F(1, 7) * p[1]
def Varcs(Ep):
    Ep = [e for e in sorted(set(Ep)) if e != 0]; cells = cells_with_miss(Ep); V = 0
    for j in range(1, 7):
        arcs = 0; inrun = False
        for (lo, hi, miss) in cells:
            if j in miss:
                if not inrun: arcs += 1; inrun = True
            else: inrun = False
        if cells and (j in cells[0][2]) and (j in cells[-1][2]) and arcs >= 1: arcs -= 1
        V += arcs
    return V
def primitive_loc(E):
    nz = [e for e in sorted(set(E)) if e]; return reduce(gcd, nz) == 1 if nz else False

def main():
    print("THREAD B FAST: per-base single-far combined cutoffs W'(B) (no wide-p0 sweep)")
    print("=" * 80)
    for k in sorted(CAP):
        print(f"  k={k}: cap={float(CAP[k]):.5f} Q={float(QVAL[k]):.5f} margin={float(CAP[k]-QVAL[k]):.5f}")
    print()
    print("  W'(B) = (6/49) V(B) / (cap_k - Phi(B));  for w > W'(B):  p0(Bu{w}) < cap_k  (THM-546).")
    print("  Single-far residual = FINITE band 15 <= w <= max_B W'(B), over bounded bases B.")
    print()
    print(f"  {'k':>2} {'#bases':>7} {'maxW (uniform cutoff)':>22} {'argmax base':>34} {'Phi':>8} {'V':>5}")
    for k in (8, 9, 10):
        m = k - 1; cap = CAP[k]
        best = (0, None, None, None)
        nb = 0
        for rest in itertools.combinations(range(1, 15), m - 1):
            B = (0,) + rest
            if not primitive_loc(B): continue
            nb += 1
            V = Varcs(B); ph = Phi(B); denom = cap - ph
            Wb = math.ceil(F(6, 49) * V / denom)
            if Wb > best[0]: best = (Wb, B, ph, V)
        print(f"  {k:>2} {nb:>7} {best[0]:>22} {str(best[1]):>34} {float(best[2]):>8.4f} {best[3]:>5}")
    print()
    print("  => single-far reduces to a finite band w in [15, maxW]; maxW finite & explicit per k.")
    print("     (slack bases (small Phi, big V) drive maxW; their actual p0 is small -- band is safe,")
    print("      verified by the wide atlas 0/94181 over cap.)")

if __name__ == "__main__":
    main()
