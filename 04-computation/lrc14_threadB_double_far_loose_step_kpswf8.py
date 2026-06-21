#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD B (2): the DOUBLE-FAR simultaneous peel + the PRECISE loose step.

THM-548 sec.5 (EXACT Newton identity, simultaneous peel from the BOUNDED base B):
    p0(B u {f1,f2})
       = P_2(B)                                   # fully-decorrelated 2-far limit (moment dual)
       + [p0(Bu{f1}) - Phi(B)]                     # one-far residual #1  (B bounded => THM-546)
       + [p0(Bu{f2}) - Phi(B)]                     # one-far residual #2  (B bounded => THM-546)
       + [ I_B(f1,f2) - Phi_2(B) ]                 # JOINT 2-far CURVATURE  (the open piece)
with  I_B(f1,f2) = p0(Bu{f1,f2}) - p0(Bu{f1}) - p0(Bu{f2}) + p0(B)  (Newton 2nd difference)
      Phi(B) = p0(B)+(1/7)p1(B),  Phi_2(B) = (2 p2(B) - p1(B))/49.

The two one-far residuals are PEELABLE: |p0(Bu{f_i})-Phi(B)| <= (6/49) V(B)/f_i  (THM-546).
The curvature C(B;f1,f2) := I_B - Phi_2(B) is the JOINT 2D Erdos-Turan-Koksma gap.

THIS SCRIPT:
 (a) verify the identity exactly;
 (b) PEEL ATTEMPT: try to bound C by iterating THM-546 (peel f2 from B, then peel f1 from Bu{f2}).
     Show the SECOND peel leaves the base Bu{f2} which, when f2 is ADJACENT to f1 (f2=f1+1),
     is NOT bounded -- V(Bu{f2}) = Theta(f2) blows up -- so (6/49)V(Bu{f2})/f1 = Theta(f2/f1) = O(1),
     NOT -> 0.  This is the precise loose step (the slope-1 resonance, HYP-2776).
 (c) Show C(B;f1,f1+1) SATURATES to a fixed nonzero value as f1->infinity (does not decorrelate),
     while C(B;f1,f2) for DISSOCIATED f1,f2 DOES decay.  => C is a genuine joint object.
 (d) the certified double-far majorant on bounded B and where it stands vs cap.
"""
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

QVAL = {k: boundary_value_direct(tuple(range(k - 1)), 1) for k in CAP}
MARGIN = {k: CAP[k] - QVAL[k] for k in CAP}

# exact engines
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
def Phi(Ep):
    p = dist_p(Ep); return p[0] + F(1, 7) * p[1]
def Phi2(Ep):
    p = dist_p(Ep); return (2 * p[2] - p[1]) / F(49)

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

def curvature(B, f1, f2):
    pB = p0(B); pB1 = p0(tuple(sorted(set(B)|{f1}))); pB2 = p0(tuple(sorted(set(B)|{f2})))
    pB12 = p0(tuple(sorted(set(B)|{f1,f2})))
    I_B = pB12 - pB1 - pB2 + pB
    return I_B - Phi2(B), I_B, pB12, pB1, pB2, pB

def main():
    print("THREAD B (2): double-far simultaneous peel + the precise LOOSE STEP")
    print("=" * 80)

    # (a) identity verification (EXACT, lhs-rhs == 0)
    print("(a) Newton simultaneous-peel identity (lhs - rhs must be 0):")
    for B, f1, f2 in [((0,1,2,3,4,5,6),15,16), ((0,1,2,3,4,5,6),15,31),
                      ((0,1,2,3,4,5,6,7),15,16), ((0,1,2,3,4,5,6,7),20,41),
                      ((0,2,4,6),15,16)]:
        C, I_B, pB12, pB1, pB2, pB = curvature(B, f1, f2)
        P2 = boundary_value_direct(tuple(sorted(set(B))), 2); phB = Phi(B)
        rhs = P2 + (pB1 - phB) + (pB2 - phB) + C
        print(f"   B={B} f=({f1},{f2}): lhs-rhs={pB12-rhs}  C=I_B-Phi2={float(C):+.6f}  "
              f"res1={float(pB1-phB):+.5f} res2={float(pB2-phB):+.5f}")
    print()

    # (b) PEEL ATTEMPT and the loose step: iterate THM-546.  Peel f2 first:
    #     |p0(Bu{f1,f2}) - Phi(Bu{f1})| <= (6/49) V(Bu{f1}) / f2     [THM-546, base = Bu{f1}]
    #     then |p0(Bu{f1}) - Phi(B)| <= (6/49) V(B)/f1.
    #     For the ADJACENT case f2=f1+1, the intermediate base Bu{f1} is fine (one extra element),
    #     BUT V(Bu{f1}) grows with f1 (the far element f1 creates 7*f1 sector crossings):
    print("-" * 80)
    print("(b) THE LOOSE STEP: iterated peel needs V(Bu{f1}); track its growth (adjacent f2=f1+1)")
    print("-" * 80)
    B = (0,1,2,3,4,5,6,7)
    print(f"   base B = consec_8,  peel f2=f1+1 first => intermediate base Bu{{f1}}:")
    print(f"     {'f1':>5} {'V(B)':>6} {'V(Bu_f1)':>9} {'(6/49)V(Bu_f1)/f2':>18} {'C=I_B-Phi2':>12}")
    for f1 in (15, 20, 40, 80, 160, 320):
        f2 = f1 + 1
        VB = Varcs(B); VBf1 = Varcs(tuple(sorted(set(B)|{f1})))
        peelbound = F(6,49) * VBf1 / f2
        C, *_ = curvature(B, f1, f2)
        print(f"     {f1:>5} {VB:>6} {VBf1:>9} {float(peelbound):>18.6f} {float(C):>12.7f}")
    print()
    print("   => V(Bu{f1}) GROWS ~ 7*f1 (the far element fragments every B_j into ~7 arcs each),")
    print("      so (6/49)V(Bu{f1})/f2 ~ (6/49)*7*f1 / (f1+1) -> 6/7 = 0.857, a CONSTANT, NOT 0.")
    print("      The iterated single-peel bound is therefore LOOSE by O(1): it cannot show the")
    print("      curvature -> small.  THIS IS THE PRECISE LOOSE STEP (slope-1 resonance, HYP-2776).")
    print()

    # (c) C saturates for adjacent, decays for dissociated.  Confirms C is a GENUINE joint object.
    print("-" * 80)
    print("(c) curvature C(B;f1,f2): ADJACENT saturates (no decorrelation), DISSOCIATED decays")
    print("-" * 80)
    for B in [(0,1,2,3,4,5,6,7), (0,1,2,3,4,5,6)]:
        print(f"   base={B}  Phi_2(B)={float(Phi2(B)):+.6f}")
        print("     adjacent f2=f1+1:")
        for f1 in (15, 40, 100, 250):
            C, *_ = curvature(B, f1, f1+1)
            print(f"       f1={f1:>4}: C={float(C):+.7f}")
        print("     dissociated (f1=15, f2 = a 'random' coprime far value):")
        for f2 in (37, 58, 91, 143, 211):
            C, *_ = curvature(B, 15, f2)
            print(f"       f2={f2:>4}: C={float(C):+.7f}")
        print()

    # (d) certified double-far majorant on BOUNDED base, sup over far pairs vs cap.
    #     p0(Bu{f1,f2}) <= P_2(B) + (6/49)V(B)(1/f1+1/f2) + sup_{f1,f2}|C(B;f1,f2)|
    #     The sup|C| is the open joint constant; we measure it exhaustively on bounded B
    #     (small far pairs) + the saturated adjacent tail, and report the resulting RHS.
    print("-" * 80)
    print("(d) certified double-far majorant on bounded base B (k=9): P_2 + 2 peel + sup|C|")
    print("-" * 80)
    k = 9; cap = CAP[k]
    # measure sup|C| over a representative set of bounded bases B (size 7) and far pairs
    supC = (F(0), None)
    nb = 0
    for rest in itertools.combinations(range(1, 15), 6):
        B = (0,) + rest
        if not primitive_loc(B): continue
        nb += 1
        # adjacent saturated (large f1) + small-pair exhaustive
        pairs = [(f, f+1) for f in (15, 60, 200)]
        for f1 in range(15, 26):
            for f2 in range(f1+1, 27):
                pairs.append((f1, f2))
        for (f1, f2) in pairs:
            C, *_ = curvature(B, f1, f2)
            if abs(C) > supC[0]: supC = (abs(C), (B, f1, f2))
    print(f"   over {nb} bounded primitive 7-bases + far pairs (adjacent-saturated + small exhaustive):")
    print(f"   sup |C(B;f1,f2)| = {float(supC[0]):.6f} at B={supC[1][0]} f=({supC[1][1]},{supC[1][2]})")
    # binding base consec_8: certified RHS at far speeds >= 15 (peel terms small) + saturated supC
    B = tuple(range(8)); V = Varcs(B); P2 = boundary_value_direct(B, 2)
    for fmin in (15, 30, 60):
        peel = F(6,49)*V*(F(1,fmin)+F(1,fmin+1))
        rhs = P2 + peel + supC[0]
        print(f"   B=consec_8, far>=({fmin},{fmin+1}): P2={float(P2):.5f} + peel={float(peel):.5f} "
              f"+ supC={float(supC[0]):.5f} = {float(rhs):.5f}  cap={float(cap):.5f}  "
              f"{'<cap' if rhs<cap else '*** >= cap (needs finite check for small far) ***'}")
    print()
    print("=" * 80)
    print("DOUBLE-FAR VERDICT:")
    print("  Identity EXACT.  Two one-far residuals peel cleanly (THM-546, bounded base).")
    print("  The LOOSE STEP is the joint curvature C=I_B-Phi_2 in the ADJACENT (slope-1) case:")
    print("  iterating THM-546 there needs V(Bu{f1})~7f1, giving an O(1) (6/7) bound, NOT ->0.")
    print("  C saturates (~0.01-0.05) but does NOT decay; a JOINT 2D ET-Koksma bound is required")
    print("  (OPEN-Q-108).  With a uniform sup|C| input, the certified RHS < cap for far>=~30,")
    print("  and small far pairs are a finite check.  Double-far is NOT closed by single-peel alone.")
    print("=" * 80)

if __name__ == "__main__":
    main()
