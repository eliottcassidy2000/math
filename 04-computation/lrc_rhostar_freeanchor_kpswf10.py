#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_rhostar_freeanchor_kpswf10.py  (kind-pasteur 2026-06-21, THM-527 Thread A)

THE FIX: free-anchor rho*.  rho*_Vmax can be 0 on admissible covering sets
(established).  But criterion C(S) takes the BEST removed runner v, not just
Vmax.  We test whether the FREE-ANCHOR density

   Rho*_free(S) := max_{v in S} meas{ tau : ||u tau|| >= 1/14 for all u != v
                                            AND ||v tau|| >= 1/14
                                            AND the v-period certifies via C }

stays positive on the rho*=0 shapes.  More directly and rigorously, we work with
the actual criterion-C quantity per anchor:

   For each v in S, the "good v-periods" are the level-1/14 safe gaps of v that
   contain a point safe for ALL of S.  The COUNT of such good periods, divided by
   v, is the v-anchored good density rho*_v.  C(S) via v  <=>  W(S\{v})>1/(7v),
   which (slow-fast) corresponds to rho*_v > 0 in the w0->inf limit.

We DIRECTLY compute, for each anchor v, the EXACT good-period density
   rho*_v(S) = meas{ tau in [0,1) : ||u tau|| >= 1/14  for all u in S }  RESTRICTED
              to the structure "tau in a v-safe gap with a full-S-safe subarc",
which at the measure level is simply the GLOBAL safe set measure

   L(S) = meas{ tau : ||u tau|| >= 1/14  for all u in S }   (the lonely measure)

split by which v-gap it lies in.  KEY POINT: L(S) > 0  <=>  some global witness
exists at the MEASURE level; but a CRITERION-C certificate via v needs a safe
subarc INSIDE A SINGLE v-gap wide enough.  We compute BOTH:
  (i) L(S) (global lonely measure) -- is it > 0?  (M(S)>=1/14 needs the CLOSED
      safe set nonempty, slightly weaker than L>0, but L>0 => M>1/14 strictly).
  (ii) per-anchor criterion margin W(S\{v})*7v -- which v rescue?

This tells us: on the rho*_Vmax=0 admissible sets, (a) is L(S)>0 (a measure
witness)?  (b) which anchors give C?  => the free-anchor route's health.
"""
import itertools
import sys
from fractions import Fraction as Fr
from math import gcd
from functools import reduce


def safe_set_measure(S, thr=Fr(1, 14)):
    """L(S) = meas{ tau : ||v tau|| >= thr for all v in S }, exact."""
    bps = {Fr(0), Fr(1)}
    for v in S:
        for m in range(0, v):
            for r in (1, 13):
                x = Fr(14 * m + r, 14 * v)
                if 0 < x < 1:
                    bps.add(x)
    pts = sorted(bps)
    tot = Fr(0)
    for x0, x1 in zip(pts, pts[1:]):
        mid = (x0 + x1) / 2
        ok = True
        for v in S:
            f = (Fr(v) * mid) % 1
            dd = f if f <= Fr(1, 2) else 1 - f
            if dd < thr:
                ok = False
                break
        if ok:
            tot += (x1 - x0)
    return tot


def W_widest_arc(A, thr=Fr(1, 14)):
    bps = {Fr(0), Fr(1)}
    for a in A:
        for m in range(0, a):
            for r in (1, 13):
                v = Fr(14 * m + r, 14 * a)
                if 0 < v < 1:
                    bps.add(v)
    pts = sorted(bps)
    arcs = []
    for x0, x1 in zip(pts, pts[1:]):
        mid = (x0 + x1) / 2
        ok = True
        for a in A:
            f = (Fr(a) * mid) % 1
            dd = f if f <= Fr(1, 2) else 1 - f
            if dd < thr:
                ok = False
                break
        if ok:
            if arcs and arcs[-1][1] == x0:
                arcs[-1] = (arcs[-1][0], x1)
            else:
                arcs.append((x0, x1))
    bestw = Fr(0)
    bestarc = None
    for a, b in arcs:
        if b - a > bestw:
            bestw = b - a
            bestarc = (a, b)
    return bestw, bestarc


def per_anchor_C(S):
    """For each v: margin = 7v*W(S\{v}). Return sorted list (margin, v)."""
    S = sorted(set(S))
    res = []
    for v in S:
        A = [u for u in S if u != v]
        W, _ = W_widest_arc(A)
        res.append((W * 7 * v, v))
    res.sort(reverse=True)
    return res


def main():
    print("=" * 78)
    print("THM-527 Thread A: FREE-ANCHOR rescue of the rho*_Vmax=0 admissible sets")
    print("=" * 78)
    sys.stdout.flush()

    # admissible covering+primitive S3 sets with rho*_Vmax = 0 (from the
    # admissibility probe). Use the concrete Vmax instances found.
    sets = [
        [1, 2, 3, 12, 20, 21, 22, 23, 24, 26, 27, 28, 31],   # E1 Vmax=31, M=6/43
        [1, 2, 3, 12, 26, 27, 28, 29, 30, 32, 33, 34, 37],   # E1 Vmax=37, M=4/31
        [1, 2, 3, 13, 20, 22, 23, 24, 25, 26, 27, 28, 29],   # E2 Vmax=29, M=5/47
        [1, 2, 3, 12, 17, 20, 22, 23, 24, 25, 26, 27, 28],   # E3 Vmax=28, M=3/25
        [1, 2, 3, 13, 22, 23, 24, 26, 27, 28, 29, 30, 31],   # E4 Vmax=31, M=8/57
    ]
    for S in sets:
        Vmax = max(S)
        L = safe_set_measure(S)
        anchors = per_anchor_C(S)
        rescuers = [(float(m), v) for (m, v) in anchors if m > 1]
        print(f"\nS = {S}")
        print(f"  Vmax = {Vmax};  L(S) = lonely measure = {L} = {float(L):.6f}  "
              f"({'>0 measure witness' if L > 0 else 'L=0 (closed-set only)'})")
        print(f"  per-anchor C margin 7v*W(S\\v), top 5:")
        for (m, v) in anchors[:5]:
            tag = "  <-- Vmax" if v == Vmax else ""
            tag2 = "  RESCUES (C holds)" if m > 1 else ""
            print(f"      v={v:3d}: margin={float(m):.4f}{tag}{tag2}")
        vmax_margin = next(m for (m, v) in anchors if v == Vmax)
        print(f"  => via-Vmax margin = {float(vmax_margin):.4f} "
              f"({'C via Vmax' if vmax_margin>1 else 'Vmax FAILS (consistent with rho*_Vmax=0)'})")
        print(f"  => # anchors that rescue (C holds): {len(rescuers)}  "
              f"best at v={anchors[0][1]} margin={float(anchors[0][0]):.4f}")
        sys.stdout.flush()

    # ---- census: over MANY rho*_Vmax=0 admissible sets, does SOME anchor always
    #      rescue (C holds via free v)?  And is L(S) always > 0? ----
    print("\n" + "=" * 78)
    print("CENSUS: build admissible covering+primitive S3 with rho*_Vmax=0,")
    print("check free-anchor C-rescue + L(S)>0 rate")
    print("=" * 78)
    # regenerate from the 4 zero-shapes over a Vmax sweep
    zeros = [
        ([0, 3, 4, 5, 7, 8, 9, 10, 11], [1, 2, 3, 12]),
        ([0, 1, 2, 3, 4, 5, 6, 7, 9], [1, 2, 3, 13]),
        ([0, 1, 2, 3, 4, 5, 6, 8, 11], [1, 2, 3, 12]),
        ([0, 1, 2, 3, 4, 5, 7, 8, 9], [1, 2, 3, 13]),
    ]
    tot = 0
    no_rescue = 0
    L_zero = 0
    worstL = None
    worst_minmargin = None
    for (E, P) in zeros:
        sp = max(E)
        for Vmax in range(14 + sp, 14 + sp + 120):
            L = [Vmax - e for e in E]
            if len(set(L)) != len(L) or min(L) <= 13:
                continue
            S = sorted(set(P) | set(L))
            if len(S) != 13 or reduce(gcd, S) != 1:
                continue
            # covering?
            if not all(any(v % q == 0 for v in S) for q in range(2, 15)):
                continue
            tot += 1
            anchors = per_anchor_C(S)
            bestm = anchors[0][0]
            if bestm <= 1:
                no_rescue += 1
            Lm = safe_set_measure(S)
            if Lm == 0:
                L_zero += 1
            if worstL is None or Lm < worstL:
                worstL = Lm
            if worst_minmargin is None or bestm < worst_minmargin:
                worst_minmargin = bestm
    print(f"  admissible rho*_Vmax=0 sets built: {tot}")
    print(f"  sets where NO anchor rescues C (max margin <= 1): {no_rescue}")
    print(f"  sets with L(S)=0 (no measure witness): {L_zero}")
    print(f"  worst (smallest) best-anchor margin = {float(worst_minmargin):.4f} "
          f"({'all rescued, margin>1' if worst_minmargin>1 else 'SOME NOT rescued'})")
    print(f"  worst (smallest) L(S) = {float(worstL):.6f}")

    print("\n" + "=" * 78)
    print("READING:")
    print("  If every rho*_Vmax=0 admissible set is C-rescued by a FREE anchor")
    print("  (max margin > 1) and L(S) > 0, then the CORRECT object is the")
    print("  FREE-ANCHOR density (max over v), and the via-Vmax pin was the bug.")
    print("  The compactness floor should be posed for max_v rho*_v, not rho*_Vmax.")
    print("=" * 78)
    print("\nDONE.")


if __name__ == "__main__":
    main()
