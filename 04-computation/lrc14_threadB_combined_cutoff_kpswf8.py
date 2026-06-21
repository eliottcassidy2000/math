#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD B (1), the CORRECT per-base combined cutoff (the actual single-far proof).

The right majorant is NOT Phi(B) <= Q then +error, but the per-base
    p0(B u {w}) = Phi(B) + Delta_w <= Phi(B) + (6/49) V(B) / w   [EXACT decomp + PROVED THM-546]
so  p0 < cap_k   whenever   w > W'(B) := (6/49) V(B) / (cap_k - Phi(B)).
This is tight per-base: high-V bases (slack regime) have SMALL Phi(B), so cap-Phi is large and
W'(B) is moderate; large-Phi bases (binding, near consec) have small V, so W'(B) is small too.

Proof of single-far wide:
  (TAIL)  for each bounded primitive base B and w > W'(B):  p0(Bu{w}) < cap_k.   PROVED (THM-546).
  (BAND)  finite exact check of all (B, 15<=w<=W'(B)).
We compute max W'(B) (the uniform cutoff), and run the BAND check.
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
MARGIN = {k: CAP[k] - QVAL[k] for k in CAP}

def orbit_breakpoints(Ep):
    Ep = sorted(set(Ep)); bp = {F(0), F(1)}
    for e in Ep:
        if e == 0: continue
        for j in range(0, 7 * e + 1): bp.add(F(j, 7 * e))
    return sorted(b for b in bp if 0 <= b < 1)

def cells_with_miss(Ep, bp=None):
    Ep = [e for e in sorted(set(Ep)) if e != 0]
    if bp is None: bp = orbit_breakpoints(Ep)
    out = []
    for lo, hi in zip(bp, bp[1:] + [F(1)]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        hit = set(int((e * mid) % 1 * 7) for e in Ep)
        out.append((lo, hi, frozenset(set(range(1, 7)) - hit)))
    return out

def dist_p(E):
    Eps = [e for e in sorted(set(E)) if e != 0]; bp = {F(0), F(1)}
    for e in Eps:
        for j in range(0, 7 * e + 1): bp.add(F(j, 7 * e))
    bp = sorted(b for b in bp if 0 <= b < 1); p = [F(0)] * 7
    for lo, hi in zip(bp, bp[1:] + [F(1)]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        hit = set(int((e * mid) % 1 * 7) for e in Eps)
        p[len(set(range(1, 7)) - hit)] += hi - lo
    return p

def p0(E):
    p = dist_p(E); return p[0]

def Phi(Ep):
    p = dist_p(Ep); return p[0] + F(1, 7) * p[1]

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
    print("THREAD B (1): per-base COMBINED cutoff W'(B)=(6/49)V(B)/(cap-Phi(B)) + band check")
    print("=" * 80)
    for k in sorted(CAP):
        print(f"  k={k}: cap={float(CAP[k]):.5f} Q(k-1)={float(QVAL[k]):.5f} margin={float(MARGIN[k]):.5f}")
    print()
    for k in (8, 9, 10):
        m = k - 1; cap = CAP[k]
        maxW = (0, None); n_bases = 0; n_cfg = 0
        worst_p0 = (F(0), None)
        # also track: does any base have cap-Phi(B) <= 0 (would break the cutoff)?  Phi<=Q<cap so no.
        for rest in itertools.combinations(range(1, 15), m - 1):
            B = (0,) + rest
            if not primitive_loc(B): continue
            n_bases += 1
            V = Varcs(B); ph = Phi(B)
            denom = cap - ph
            assert denom > 0, f"cap-Phi<=0 at B={B}"   # Phi(B)<=Q(k-1)<cap, always positive
            Wb = math.ceil(F(6, 49) * V / denom)        # combined per-base cutoff
            if Wb > maxW[0]: maxW = (Wb, B)
            for w in range(15, Wb + 1):
                E = tuple(sorted(set(B) | {w}))
                if not primitive_loc(E): continue
                n_cfg += 1
                pv = p0(E)
                if pv > worst_p0[0]: worst_p0 = (pv, (B, w))
        print(f"  k={k}: bounded bases={n_bases}, max combined cutoff W'={maxW[0]} at B={maxW[1]}")
        print(f"        band configs checked (15<=w<=W'(B)): {n_cfg}")
        print(f"        max p0 in band = {float(worst_p0[0]):.6f} at {worst_p0[1]}  cap={float(cap):.6f}  "
              f"margin={float(cap-worst_p0[0]):.6f}  "
              f"{'ALL<cap -> SINGLE-FAR CLOSED' if worst_p0[0]<cap else '*** VIOLATION ***'}")
    print()
    print("=" * 80)
    print("VERDICT (single-far): TAIL w>W'(B) PROVED by THM-546; BAND 15<=w<=W'(B) exact-checked.")
    print("If max p0 in band < cap for all k, single-far wide is PROVED < cap (k=8,9,10 here).")
    print("=" * 80)

if __name__ == "__main__":
    main()
