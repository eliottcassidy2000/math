#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD B (kps-2026-06-21-S26-wf8): the LOCALIZED err+ bound on the BINDING family.

The wide bound is `span(E) > 14  =>  p0(E) = measS7(E) < cap_k`.  Thread-1's regime split:
 - BINDING regime: p0 near cap  <=>  LARGE base B (~ consec_{k-1}) + far elements just past 14.
   Here Phi(B) is large (near the proven plateau max Q(k-1)) but V(B) is BOUNDED, so the
   THM-546 comb bound |Delta_w| <= (6/49) V(B) / w controls the error EXPLICITLY.
 - SLACK regime: large error <=> SMALL/clustered base => p0_decorr small => cap - p0 >> err.

This script does Thread B exactly:

(1) SINGLE FAR.  E = B u {w}, B a bounded near-consec (k-1)-set, w > 14.  Exactly
        p0(E) = Phi(B) + Delta_w,    Phi(B) := p0(B) + (1/7) p1(B)   [HYP-2642 recursion]
    with the PROVEN plateau-max  Phi(B) <= Q(k-1) = Phi(consec_{k-1})  (THM plateau lemma)
    and the PROVEN comb bound (THM-546 signed-Abel form)
        |Delta_w| <= (6/49) * V(B) / w,   V(B) := sum_{j=1..6} #arcs(B_j),  #arcs(B_j) <= 7*sum_{e in B} e... (bounded base).
    Wide => w >= 15.  The WORST (binding) base is B = consec_{k-1} (largest Phi AND largest V
    among span<=14 bases) at the SMALLEST wide speed w = 15.  Compute the explicit RHS
        Q(k-1) + (6/49) V(consec_{k-1}) / 15   vs cap_k.
    Does it PROVE p0 < cap for single-far?  (We also sweep ALL bounded bases B and the actual
    worst V to give the genuine RHS, and verify against the EXACT p0.)

(2) DOUBLE FAR.  E = B u {f1, f2}.  THM-548 sec.5 simultaneous (not iterated) peel:
        p0(B u {f1,f2}) = P_2(B) + [p0(Bu{f1}) - Phi(B)] + [p0(Bu{f2}) - Phi(B)] + [I_B(f1,f2) - Phi_2(B)]
    The two ONE-far residuals are over the BOUNDED base B => THM-546 gives <= (6/49)V(B)/f_i each.
    The CURVATURE term I_B(f1,f2) - Phi_2(B) is the joint 2D Erdos-Turan-Koksma gap.
    Attempt the bound and identify EXACTLY where the single-peel argument is loose
    (the adjacent f2 = f1 + 1 slope-1 resonance, HYP-2776).

EXACT rationals throughout.  Engine functions are imported VERBATIM from the resonance-direct
file (pure math, no stale CAP); CAP/Q come from the AUTHORITATIVE wide_branch_ridge engine.
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

# ---- authoritative CAP and the moment-dual boundary value (decorrelated baseline) ----
from lrc14_wide_branch_ridge_codex_s47 import CAP, primitive
from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct

QVAL = {k: boundary_value_direct(tuple(range(k - 1)), 1) for k in CAP}
MARGIN = {k: CAP[k] - QVAL[k] for k in CAP}

# ---------------------------------------------------------------------------
# EXACT ENGINES (ported verbatim from lrc14_ck_resonance-direct_kps-S17-wf.py)
# ---------------------------------------------------------------------------

def G0(y):
    """Centered antiderivative of (1_{[0,1/7)} - 1/7), periodised; |G0| <= 6/49.
    On [0,1/7): slope 6/7; on [1/7,1): slope -1/7."""
    y = y - int(y)
    if y < 0:
        y += 1
    return y * F(6, 7) if y < F(1, 7) else F(6, 49) - (y - F(1, 7)) * F(1, 7)

def orbit_breakpoints(Ep):
    Ep = sorted(set(Ep))
    bp = {F(0), F(1)}
    for e in Ep:
        if e == 0:
            continue
        for j in range(0, 7 * e + 1):
            bp.add(F(j, 7 * e))
    return sorted(b for b in bp if 0 <= b < 1)

def cells_with_miss(Ep, bp=None):
    Ep = [e for e in sorted(set(Ep)) if e != 0]
    if bp is None:
        bp = orbit_breakpoints(Ep)
    out = []
    for lo, hi in zip(bp, bp[1:] + [F(1)]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        hit = set(int((e * mid) % 1 * 7) for e in Ep)
        miss = set(range(1, 7)) - hit
        out.append((lo, hi, frozenset(miss)))
    return out

def wDelta_signed(Ep, w, bp=None):
    """w*Delta_w (SIGNED, exact). Sum over |miss|==1 cells of G0(w*hi - s/7) - G0(w*lo - s/7)."""
    Ep = [e for e in sorted(set(Ep)) if e != 0]
    if bp is None:
        bp = orbit_breakpoints(Ep)
    D = F(0)
    for lo, hi in zip(bp, bp[1:] + [F(1)]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        hit = set(int((e * mid) % 1 * 7) for e in Ep)
        miss = set(range(1, 7)) - hit
        if len(miss) == 1:
            s = next(iter(miss))
            D += G0(w * hi - F(s, 7)) - G0(w * lo - F(s, 7))
    return D

def p0(E):
    Eps = [e for e in sorted(set(E)) if e != 0]
    bp = {F(0), F(1)}
    for e in Eps:
        for j in range(0, 7 * e + 1):
            bp.add(F(j, 7 * e))
    bp = sorted(b for b in bp if 0 <= b < 1)
    tot = F(0)
    for lo, hi in zip(bp, bp[1:] + [F(1)]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        hit = set(int((e * mid) % 1 * 7) for e in Eps)
        if set(range(1, 7)) <= hit:
            tot += hi - lo
    return tot

def dist_p(E):
    Eps = [e for e in sorted(set(E)) if e != 0]
    bp = {F(0), F(1)}
    for e in Eps:
        for j in range(0, 7 * e + 1):
            bp.add(F(j, 7 * e))
    bp = sorted(b for b in bp if 0 <= b < 1)
    p = [F(0)] * 7
    for lo, hi in zip(bp, bp[1:] + [F(1)]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        hit = set(int((e * mid) % 1 * 7) for e in Eps)
        p[len(set(range(1, 7)) - hit)] += hi - lo
    return p

def Phi(Ep):
    """Plateau Phi(E') = p0(E') + (1/7) p1(E')."""
    p = dist_p(Ep)
    return p[0] + F(1, 7) * p[1]

def Phi2(Ep):
    """Two-far decorrelated limit Phi_2 = (2 p2 - p1) / 49 (THM-548)."""
    p = dist_p(Ep)
    return (2 * p[2] - p[1]) / F(49)

def Varcs(Ep):
    """V(E') := sum_{j=1..6} #arcs(B_j(E')) = number of MAXIMAL runs of cells whose miss-set
    contains sector j, summed over j=1..6.  This is the BV arc-complexity in THM-546's bound.
    (An 'arc' of B_j = a maximal interval of x on which E' misses sector j.)"""
    Ep = [e for e in sorted(set(Ep)) if e != 0]
    cells = cells_with_miss(Ep)  # (lo,hi,miss) in x-order, covering [0,1)
    V = 0
    for j in range(1, 7):
        inrun = False
        arcs = 0
        for (lo, hi, miss) in cells:
            if j in miss:
                if not inrun:
                    arcs += 1
                    inrun = True
            else:
                inrun = False
        # wrap-around: B_j is periodic, merge first & last cell if both contain j
        if cells and (j in cells[0][2]) and (j in cells[-1][2]) and arcs >= 1:
            # the run crossing 0 was double-counted as two arcs; merge
            # only merge if the first cell starts a run and last cell ends in a run
            arcs -= 1
        V += arcs
    return V

def primitive_loc(E):
    nz = [e for e in sorted(set(E)) if e]
    return reduce(gcd, nz) == 1 if nz else False

# ---------------------------------------------------------------------------
# CROSS-CHECKS
# ---------------------------------------------------------------------------

def cross_checks():
    print("=" * 80)
    print("(0) ENGINE CROSS-CHECKS (exact)")
    print("=" * 80)
    # (a) Delta_w = p0(B u {w}) - Phi(B), and wDelta_signed = w * Delta_w
    ok = True
    for B, w in [((0,1,2,3,4,5,6,7),30), ((0,1,2,3,4,5,6,7),15),
                 (tuple(range(8)),16), ((0,2,4,6,8,10,12,14),15)]:
        E = tuple(sorted(set(B) | {w}))
        Dw = p0(E) - Phi(B)
        wDw = wDelta_signed(B, w)
        ok &= (wDw == w * Dw)
        print(f"  B={B} w={w}: Delta_w={float(Dw):+.6f}  w*Delta_w(engine)={wDw}  "
              f"match={'OK' if wDw == w*Dw else 'MISMATCH'}")
    # (b) THM-546 comb bound: |Delta_w| <= (6/49) V(B) / w  ?  (verify on many B,w)
    print()
    print("  THM-546 comb bound check  |Delta_w| <= (6/49) V(B) / w :")
    bad = 0; tested = 0
    import random
    rng = random.Random(3)
    bases = [tuple(range(7)), tuple(range(8)), (0,2,4,6,8,10,12,14),
             (0,1,2,4,8,9,11), (0,1,3,5,7,9,11,13)]
    while len(bases) < 40:
        Bx = tuple(sorted(set([0] + rng.sample(range(1, 15), rng.choice([6,7,8])))))
        if Bx[0] == 0 and len(Bx) >= 6:
            bases.append(Bx)
    worst_ratio = F(0)
    for B in bases:
        V = Varcs(B)
        for w in range(15, 80):
            Dw = abs(p0(tuple(sorted(set(B)|{w}))) - Phi(B))
            rhs = F(6, 49) * V / w
            tested += 1
            if Dw > rhs:
                bad += 1
                if bad <= 5:
                    print(f"    VIOLATION B={B} w={w}: |Delta_w|={float(Dw):.6f} > (6/49)V/w={float(rhs):.6f} (V={V})")
            if rhs > 0 and Dw / rhs > worst_ratio:
                worst_ratio = Dw / rhs
    print(f"  tested {tested} (B,w) pairs (w=15..79); violations of (6/49)V/w bound: {bad}")
    print(f"  worst |Delta_w| / [(6/49)V/w] ratio = {float(worst_ratio):.4f}  "
          f"({'bound HOLDS (<=1)' if worst_ratio <= 1 else '*** bound FAILS ***'})")
    return ok and bad == 0

# ---------------------------------------------------------------------------
# (1) SINGLE FAR -- the explicit binding bound
# ---------------------------------------------------------------------------

def single_far():
    print()
    print("=" * 80)
    print("(1) SINGLE FAR:  p0(B u {w}) = Phi(B) + Delta_w,  binding RHS vs cap_k")
    print("=" * 80)
    print("  Decomposition (EXACT):  p0(E) = Phi(B) + Delta_w  [HYP-2642].")
    print("  Plateau-max (PROVED):   Phi(B) <= Q(k-1) = Phi(consec_{k-1}).")
    print("  Comb bound (THM-546):   |Delta_w| <= (6/49) V(B) / w,  w >= 15 (wide).")
    print("  => p0(E) <= Q(k-1) + (6/49) V(B) / 15.  Worst base = consec_{k-1} (max Phi & V).")
    print()
    print(f"  {'k':>2} {'cap_k':>9} {'Q(k-1)':>9} {'margin':>8} {'V(con)':>6} "
          f"{'(6/49)V/15':>10} {'RHS_binding':>11} {'RHS<cap?':>9}")
    closes_all_single = True
    for k in sorted(CAP):
        m = k - 1
        Bcon = tuple(range(m))
        V = Varcs(Bcon)
        Q = QVAL[k]
        worst_err = F(6, 49) * V / 15   # smallest wide w = 15
        rhs = Q + worst_err
        ok = rhs < CAP[k]
        closes_all_single &= ok
        print(f"  {k:>2} {float(CAP[k]):>9.5f} {float(Q):>9.5f} {float(MARGIN[k]):>8.5f} "
              f"{V:>6} {float(worst_err):>10.5f} {float(rhs):>11.5f} {str(ok):>9}")
    print()
    print(f"  ===> Crude binding RHS Q(k-1)+(6/49)V(consec)/15 < cap_k for all k=8..12 ?  {closes_all_single}")
    print()

    # Now the GENUINE worst over ALL bounded bases B (span<=14) and ACTUAL w-decay:
    # the binding case is not literally consec at w=15 -- err and Phi trade off.  Sweep
    # all bounded primitive bases, all wide w, take the EXACT sup of p0 and the sup of the
    # certified RHS = Phi(B) + (6/49)V(B)/w  (a rigorous per-(B,w) majorant of p0).
    print("-" * 80)
    print("  EXACT sweep: all bounded primitive bases B (0 in B, max<=14), wide w in [15, WMAX].")
    print("  Report (a) sup actual p0, (b) sup certified majorant Phi(B)+(6/49)V(B)/w.")
    print("-" * 80)
    WMAX = 60
    for k in (8, 9, 10):   # base size m = k-1 = 7,8,9
        m = k - 1
        best_p0 = (F(0), None)
        best_maj = (F(0), None)   # sup over (B,w) of certified majorant
        # iterate bounded bases: 0 plus (m-1) elements from 1..14
        nb = 0
        for rest in itertools.combinations(range(1, 15), m - 1):
            B = (0,) + rest
            if not primitive_loc(B):
                continue
            nb += 1
            V = Varcs(B); ph = Phi(B)
            # for the majorant, worst w is the smallest wide w (=15); but actual p0 sweep needs all w
            maj15 = ph + F(6, 49) * V / 15
            if maj15 > best_maj[0]:
                best_maj = (maj15, (B, 15))
            for w in range(15, WMAX + 1):
                E = tuple(sorted(set(B) | {w}))
                if not primitive_loc(E):
                    continue
                pv = p0(E)
                if pv > best_p0[0]:
                    best_p0 = (pv, (B, w))
        cap = CAP[k]
        print(f"  k={k} (base size {m}, {nb} bounded primitive bases):")
        print(f"    sup ACTUAL p0  = {float(best_p0[0]):.6f}  at B={best_p0[1][0]} w={best_p0[1][1]}"
              f"   margin to cap = {float(cap - best_p0[0]):.6f}  {'<cap OK' if best_p0[0]<cap else '*** EXCEEDS ***'}")
        print(f"    sup CERTIFIED majorant Phi(B)+(6/49)V(B)/15 = {float(best_maj[0]):.6f}  at B={best_maj[1][0]}"
              f"   margin to cap = {float(cap - best_maj[0]):.6f}  "
              f"{'PROVES <cap' if best_maj[0]<cap else '*** majorant >= cap (loose) ***'}")
    print()
    print("  NOTE: the certified majorant uses the PROVED THM-546 bound per (B,w); if its SUP < cap")
    print("  then single-far wide is PROVED < cap (no consec-extremality needed -- direct).")

# ---------------------------------------------------------------------------
# (2) DOUBLE FAR -- simultaneous peel + the loose step
# ---------------------------------------------------------------------------

def Phi_full_2far_identity_check(B, f1, f2):
    """Verify THM-548 sec.5 Newton identity (EXACT):
       p0(B u {f1,f2}) = P_2(B) + [p0(Bu{f1})-Phi(B)] + [p0(Bu{f2})-Phi(B)] + [I_B - Phi_2(B)]
    where P_2(B) = boundary_value_direct(B,2), I_B = p0(B u {f1,f2}) - p0(Bu{f1}) - p0(Bu{f2}) + p0(B)
    (the two-far curvature, = Delta_{f1,f2}).  Returns (lhs, rhs, curvature, dev_curv)."""
    pB   = p0(B)
    pB1  = p0(tuple(sorted(set(B)|{f1})))
    pB2  = p0(tuple(sorted(set(B)|{f2})))
    pB12 = p0(tuple(sorted(set(B)|{f1,f2})))
    P2   = boundary_value_direct(tuple(sorted(set(B))), 2)
    phB  = Phi(B)
    ph2  = Phi2(B)
    I_B  = pB12 - pB1 - pB2 + pB            # Delta_{f1,f2} (two-far Newton difference)
    dev1 = pB1 - phB
    dev2 = pB2 - phB
    devc = I_B - ph2
    rhs  = P2 + dev1 + dev2 + devc
    return pB12, rhs, I_B, devc, dev1, dev2, P2, ph2

def double_far():
    print()
    print("=" * 80)
    print("(2) DOUBLE FAR:  simultaneous peel (THM-548 sec.5) + the loose step")
    print("=" * 80)
    print("  Newton identity (EXACT):")
    print("    p0(Bu{f1,f2}) = P_2(B) + [p0(Bu{f1})-Phi(B)] + [p0(Bu{f2})-Phi(B)] + [I_B - Phi_2(B)]")
    print("  P_2(B) = fully-decorrelated 2-far limit (moment dual); the two ONE-far residuals")
    print("  are over the BOUNDED base B => THM-546:  |p0(Bu{f_i})-Phi(B)| <= (6/49)V(B)/f_i.")
    print("  The CURVATURE I_B - Phi_2(B) is the joint 2D ET-Koksma gap (the open piece).")
    print()
    # (a) verify the identity exactly on several configs
    print("  Identity check (lhs - rhs should be 0):")
    for B, f1, f2 in [((0,1,2,3,4,5,6),15,16), ((0,1,2,3,4,5,6),15,31),
                      ((0,1,2,3,4,5,6,7),15,16), ((0,1,2,3,4,5,6,7),20,41)]:
        lhs, rhs, I_B, devc, d1, d2, P2, ph2 = Phi_full_2far_identity_check(B, f1, f2)
        print(f"    B={B} f=({f1},{f2}): lhs-rhs={lhs-rhs}  curvature dev I_B-Phi2={float(devc):+.6f}")
    print()

    # (b) the curvature term: scan over bounded B and far pairs, find sup |I_B - Phi_2|,
    #     and show it PEAKS at the adjacent pair f2=f1+1 (slope-1 resonance, HYP-2776).
    print("-" * 80)
    print("  (b) curvature deviation  I_B - Phi_2(B)  as a function of the far pair")
    print("      (base B = consec_{k-1}); ADJACENT f2=f1+1 is the resonance peak.")
    print("-" * 80)
    for m in (7, 8):  # k = 8, 9
        B = tuple(range(m))
        print(f"   base = consec_{m}  (k={m+1}):  Phi_2(B) = {float(Phi2(B)):+.6f}")
        # adjacent pairs (f1, f1+1) as f1 grows -> saturating resonance value
        print("     adjacent f2=f1+1:   f1 |  I_B - Phi_2")
        sat = None
        for f1 in (15, 20, 40, 80, 160, 320):
            f2 = f1 + 1
            _, _, I_B, devc, *_ = Phi_full_2far_identity_check(B, f1, f2)
            print(f"       {f1:>5} | {float(devc):+.7f}")
            sat = devc
        # separated pairs (dissociated) -> curvature decays toward 0
        print("     separated (f1=20, f2 sweep, pick worst & a dissociated one):")
        worst = (F(0), None); diss = None
        for f2 in range(22, 90):
            f1 = 20
            E = tuple(sorted({f1, f2}))
            _, _, I_B, devc, *_ = Phi_full_2far_identity_check(B, f1, f2)
            if abs(devc) > worst[0]:
                worst = (abs(devc), (f1, f2, devc))
            if f2 == 81:  # f2/f1 = 81/20 dissociated-ish
                diss = devc
        print(f"       worst |I_B-Phi_2| over f2 in [22,89] (f1=20): {float(worst[0]):.7f} at {worst[1]}")
        print()

    # (c) THE LOOSE STEP: the single-peel (iterated) bound vs the actual curvature.
    #     Iterating THM-546 (peel f2, then f1) leaves E'=B u {f1} which, for ADJACENT
    #     f2=f1+1, is NOT decorrelated from f2 -- the V of the intermediate set is fine
    #     (B bounded, f1 just one element) BUT the curvature I_B-Phi_2 is a genuine
    #     joint 2D object that the per-element (6/49)V/w bound does NOT see.  Quantify.
    print("-" * 80)
    print("  (c) THE LOOSE STEP -- why iterated single-peel does not, by itself, bound curvature")
    print("-" * 80)
    print("  Iterated peel:  p0(Bu{f1,f2}) - Phi_2-type decorrelation is NOT just the sum of two")
    print("  one-far (6/49)V/w terms.  The simultaneous-peel identity isolates the EXTRA joint")
    print("  curvature  C := I_B - Phi_2(B).  For ADJACENT f2=f1+1, C saturates to a fixed")
    print("  nonzero value as f1->infinity (the relation f2-f1=1 never decorrelates), so it is")
    print("  NOT bounded by (6/49)V/f1 -> 0.  This is the precise loose step (HYP-2776).")
    print()
    # Assemble the FULL certified double-far majorant on the binding base and check vs cap.
    print("  Certified double-far majorant on binding base B=consec_{k-1}, worst far config:")
    print(f"    {'k':>2} {'cap':>8} {'P_2(B)':>8} {'2*(6/49)V/15':>12} {'sup|C|':>8} {'RHS':>9} {'<cap?':>7}")
    for k in (8, 9, 10):
        m = k - 1
        B = tuple(range(m))
        V = Varcs(B)
        P2 = boundary_value_direct(B, 2)
        twopeel = 2 * F(6, 49) * V / 15        # both far at smallest wide speed 15 (and 16)
        # sup curvature: take the adjacent-pair saturated value (worst), measured large f1
        supC = F(0)
        for f1 in (15, 50, 200):
            _, _, _, devc, *_ = Phi_full_2far_identity_check(B, f1, f1 + 1)
            if abs(devc) > supC:
                supC = abs(devc)
        # also small-f exhaustive worst over f1,f2 in [15,40]
        for f1 in range(15, 41):
            for f2 in range(f1 + 1, 42):
                _, _, _, devc, *_ = Phi_full_2far_identity_check(B, f1, f2)
                if abs(devc) > supC:
                    supC = abs(devc)
        rhs = P2 + twopeel + supC
        ok = rhs < CAP[k]
        print(f"    {k:>2} {float(CAP[k]):>8.5f} {float(P2):>8.5f} {float(twopeel):>12.5f} "
              f"{float(supC):>8.5f} {float(rhs):>9.5f} {str(ok):>7}")
    print()
    print("  Interpretation: when the two one-far residuals are SMALL (w>=15 not too small) the")
    print("  RHS is P_2(B) + small + sup|C|.  P_2(B) << cap (margin grows in r), so the only")
    print("  thing standing between this and a PROOF is a uniform bound on sup|C| = sup|I_B-Phi_2|")
    print("  over ALL bounded B and ALL far pairs -- the joint 2D ET-Koksma constant (OPEN-Q-108).")

def main():
    print("THREAD B (kps-S26-wf8): localized err+ bound on the BINDING family via THM-546")
    print()
    for k in sorted(CAP):
        print(f"  k={k}: cap={float(CAP[k]):.5f} Q(k-1)={float(QVAL[k]):.5f} margin={float(MARGIN[k]):.5f}")
    ok = cross_checks()
    single_far()
    double_far()
    print()
    print("=" * 80)
    print(f"engine cross-checks all pass: {ok}")
    print("=" * 80)

if __name__ == "__main__":
    main()
