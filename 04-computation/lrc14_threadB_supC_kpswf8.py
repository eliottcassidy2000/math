#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD B (2) finalizer: sup|C| (the joint 2-far curvature) over bounded bases + certified RHS.
C(B;f1,f2) = I_B(f1,f2) - Phi_2(B), the joint 2D Erdos-Turan-Koksma gap (the open piece)."""
from __future__ import annotations
import sys, functools, itertools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
from lrc14_wide_branch_ridge_codex_s47 import CAP
from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct

def dist_p(E):
    Eps = [e for e in sorted(set(E)) if e != 0]; bp = {F(0), F(1)}
    for e in Eps:
        for j in range(0, 7 * e + 1): bp.add(F(j, 7 * e))
    bp = sorted(b for b in bp if 0 <= b < 1); p = [F(0)] * 7
    for lo, hi in zip(bp, bp[1:] + [F(1)]):
        if hi <= lo: continue
        mid = (lo + hi) / 2; hit = set(int((e * mid) % 1 * 7) for e in Eps)
        p[len(set(range(1, 7)) - hit)] += hi - lo
    return p
def p0(E): return dist_p(E)[0]
def Phi2(Ep):
    p = dist_p(Ep); return (2 * p[2] - p[1]) / F(49)
def prim(E):
    nz = [e for e in sorted(set(E)) if e]; return reduce(gcd, nz) == 1 if nz else False
def curv(B, f1, f2):
    pB = p0(B); pB1 = p0(tuple(sorted(set(B)|{f1}))); pB2 = p0(tuple(sorted(set(B)|{f2})))
    pB12 = p0(tuple(sorted(set(B)|{f1,f2})))
    return (pB12 - pB1 - pB2 + pB) - Phi2(B)
def Varcs(Ep):
    Ep = [e for e in sorted(set(Ep)) if e != 0]; bp = {F(0), F(1)}
    for e in Ep:
        for j in range(0, 7 * e + 1): bp.add(F(j, 7 * e))
    bp = sorted(b for b in bp if 0 <= b < 1); cells = []
    for lo, hi in zip(bp, bp[1:] + [F(1)]):
        if hi <= lo: continue
        mid = (lo + hi) / 2; hit = set(int((e * mid) % 1 * 7) for e in Ep)
        cells.append((lo, hi, frozenset(set(range(1, 7)) - hit)))
    V = 0
    for j in range(1, 7):
        arcs = 0; inrun = False
        for (lo, hi, miss) in cells:
            if j in miss:
                if not inrun: arcs += 1; inrun = True
            else: inrun = False
        if cells and (j in cells[0][2]) and (j in cells[-1][2]) and arcs >= 1: arcs -= 1
        V += arcs
    return V

def main():
    print("THREAD B (2): sup|C(B;f1,f2)| over bounded bases + certified double-far RHS vs cap")
    print("=" * 80)
    # sup|C| over bounded 7-bases (k=9): small exhaustive far + adjacent-saturated tail
    supC = (F(0), None); nb = 0
    for rest in itertools.combinations(range(1, 15), 6):
        B = (0,) + rest
        if not prim(B): continue
        nb += 1
        pairs = [(f, f + 1) for f in (15, 50)]            # adjacent (resonance), incl saturating
        for f1 in range(15, 22):
            for f2 in range(f1 + 1, 23): pairs.append((f1, f2))
        for (f1, f2) in pairs:
            C = curv(B, f1, f2)
            if abs(C) > supC[0]: supC = (abs(C), (B, f1, f2))
    print(f"sup|C| over {nb} bounded 7-bases, small+adjacent far pairs:")
    print(f"  sup|C| = {supC[0]} = {float(supC[0]):.6f}  at B={supC[1][0]} f=({supC[1][1]},{supC[1][2]})")
    print()
    # certified RHS for binding base consec_8: p0 <= P2 + (6/49)V(1/f1+1/f2) + sup|C|
    B = tuple(range(8)); V = Varcs(B); P2 = boundary_value_direct(B, 2); cap = CAP[9]
    print(f"binding base consec_8: V={V} P2={float(P2):.5f} cap9={float(cap):.5f} sup|C|(used)={float(supC[0]):.5f}")
    print(f"  {'far>=':>8} {'peel':>9} {'RHS':>9} {'<cap?':>20}")
    for fmin in (15, 30, 60, 100, 200):
        peel = F(6, 49) * V * (F(1, fmin) + F(1, fmin + 1))
        rhs = P2 + peel + supC[0]
        verdict = "<cap PROVED" if rhs < cap else ">=cap (finite check below fmin)"
        print(f"  {fmin:>8} {float(peel):>9.5f} {float(rhs):>9.5f} {verdict:>30}")
    print()
    print("Reading: with a UNIFORM sup|C| input (the open joint constant), the certified RHS drops")
    print("below cap once the far speeds exceed ~", end="")
    # find smallest fmin making RHS<cap
    fstar = None
    for fmin in range(15, 400):
        peel = F(6, 49) * V * (F(1, fmin) + F(1, fmin + 1))
        if P2 + peel + supC[0] < cap:
            fstar = fmin; break
    print(f"{fstar}.  Below it, a finite check (bounded base x far pairs <= {fstar}) closes the rest.")
    print("The SOLE open input is the uniform bound on sup|C| = sup|I_B - Phi_2| (OPEN-Q-108).")

if __name__ == "__main__":
    main()
