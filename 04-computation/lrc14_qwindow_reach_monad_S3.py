#!/usr/bin/env python3
"""
monad-explorer-2026-07-07-S3 -- the q-WINDOW EXPLICIT REACH LEMMA and the kill-set
classification of the true 1/14 residual (HYP-4847).

CONTEXT (the domination audit, same session): kps-S28's spread13_lonely (GREEN,
LRCSpread13.lean) kills every ratio<=13 family at the 1/14 bar, so the k=13/P=empty
floor leg -- the target of the week's mu/E[maxgap]/E[U] burst including my HYP-4827 --
is dominated for Vmax >= 82.  The TRUE residual is the MIXED shapes: small part P != empty
(ratio > 13 automatically) with a top cluster.  For those, Part A's O(#arcs/Vmax) sketch
is replaced here by a fully explicit construction:

LEMMA (q-window reach; to verify + calibrate).  S = P u C, P subset {1..13}, C the top
cluster with co-offsets e = Vmax - u in [0, D_c].  Let q in {2..6} with q NOT dividing any
p' in P ("q survives P"), gcd(p,q)=1.  Witness ansatz tau = p/q + theta/Vmax:
  * P-side:   ||p' tau|| >= 1/q - 13(1+theta)/Vmax >= 1/14      (margin 1/q - 1/14 >= 2/21)
  * cluster:  u tau = u p/q + (1 - e/Vmax)(theta) mod 1: teeth sit on the q-grid
              (u p/q mod 1, <= q values) shifted by <= D_c(1+theta)/Vmax:
              the theta-line misses all 1/14-neighborhoods on a free set of measure
              >= 1 - q/7 - C*q*D_c/Vmax.
  Sufficient condition (to calibrate): D_c/Vmax < (7-q)/(14q) * safety, Vmax >= V0(q).

CLASSIFICATION (kill-sets).  q in {2..6} survives P iff P has no multiple of q.
P kills ALL of {2..6} iff P meets {6,12}, {4,8,12}, {3,6,9,12}, {5,10} simultaneously
(evenness follows).  Minimal kill-sets: {12,5}, {12,10}, {6,4,5}, {6,4,10}, {6,8,5},
{6,8,10}.  Per skeleton leg: k=12 (|P|=1): NO kill-all P -- every k=12 shape has a
surviving q-window; k=11 (|P|=2): exactly {12,5},{12,10}; k=10,9,8: supersets counted here.

This script: (1) verify the lemma numerically over adversarial mixed shapes -- find the
TRUE D_c/Vmax boundary per q vs the predicted (7-q)/(14q); (2) enumerate kill-all P per
|P|; (3) probe the genuine hard core (kill-all P + cluster): is M >= 1/14 still holding
via other mechanisms (it must, if LRC14 is true -- measure the margin); (4) sanity: the
q-window predictions on THM-530's pathology shapes.

Tournament Analysis declaration:
  vertices: shape classes (surviving-q vs kill-all-P vs wide-cluster);
  pairwise observable: does the explicit witness ansatz certify 1/14 on the class;
  switch/gauge: orient toward the class with no certificate = the honest residual;
  tie path: lemma calibration -> kill enumeration -> hard-core probe -> residual map.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import random

import numpy as np


def reach_at(v, tau_num, tau_den):
    """min_i ||v_i * tau|| at exact rational tau."""
    t = F(tau_num, tau_den)
    return min(abs(vi * t - round(vi * t)) for vi in v)


def qwindow_witness(v, Vmax, q, thetas=4096):
    """Best min-dist over tau = p/q + theta/Vmax, gcd(p,q)=1, theta on a grid (float)."""
    best = -1.0
    va = np.array(v, dtype=np.float64)
    for p in range(1, q):
        if gcd(p, q) != 1:
            continue
        th = (np.arange(thetas) + 0.5) / thetas
        taus = p / q + th / Vmax
        ph = np.outer(taus, va)
        d = np.abs(ph - np.rint(ph)).min(axis=1)
        m = d.max()
        if m > best:
            best = m
    return best


def make_shape(P, offsets, Vmax):
    """Family = P u {Vmax - e for e in offsets} (0 in offsets => Vmax included)."""
    C = sorted(set(Vmax - e for e in offsets))
    v = sorted(set(list(P) + C))
    return v


def survives(q, P):
    return all(p % q != 0 for p in P)


if __name__ == "__main__":
    random.seed(4847)
    print("=" * 78)
    print("PART 1 -- calibrate the q-window reach boundary D_c/Vmax per q")
    print("=" * 78)
    print("  predicted sufficient: D_c/Vmax < (7-q)/(14q); measure the TRUE boundary")
    print("  (worst over random cluster shapes at each D_c; P chosen q-surviving)")
    Vmax = 100003  # large prime, clean
    for q, P in [(2, (1, 3, 5)), (3, (1, 2, 5)), (5, (1, 2, 3)), (6, (1, 5, 7))]:
        pred = (7 - q) / (14 * q)
        line = f"  q={q} P={P} pred<{pred:.4f}: "
        for frac_dc in [0.25 * pred, 0.5 * pred, 0.9 * pred, 1.5 * pred, 3 * pred, 8 * pred]:
            Dc = int(frac_dc * Vmax)
            worst = 1.0
            for trial in range(24):
                k = 13 - len(P)
                offs = [0] + sorted(random.sample(range(1, max(Dc, k + 1) + 1), k - 1))
                v = make_shape(P, offs, Vmax)
                if len(v) != 13:
                    continue
                w = qwindow_witness(v, Vmax, q, 2048)
                if w < worst:
                    worst = w
            line += f"{frac_dc/pred:.2f}x:{worst:.4f}{'Y' if worst >= 1/14 else 'n'} "
        print(line + "  [Y = witness >= 1/14 certified on all 24 shapes]")

    print()
    print("=" * 78)
    print("PART 2 -- kill-set classification: P with NO surviving q in {2..6}")
    print("=" * 78)
    for psz in range(1, 6):
        kills = [P for P in combinations(range(1, 14), psz)
                 if not any(survives(q, P) for q in range(2, 7))]
        k = 13 - psz
        ex = f"; e.g. {kills[:4]}" if kills else ""
        print(f"  |P|={psz} (k={k}): {len(kills)}/{len(list(combinations(range(1,14),psz)))} "
              f"kill-all P's{ex}")

    print()
    print("=" * 78)
    print("PART 3 -- the genuine hard core: kill-all P + tight cluster; direct M probe")
    print("=" * 78)
    print("  (these shapes have NO q<=6 window; LRC14 says they are still lonely --")
    print("   measure the actual witness structure at moderate Vmax)")
    hard_cases = [
        ((5, 12), "k=11 kill-all #1"),
        ((10, 12), "k=11 kill-all #2"),
        ((4, 5, 6), "k=10 kill-all minimal"),
        ((1, 2, 3, 5, 7, 8, 9, 11, 12, 13), "m_P-attaining P* (THM-530)"),
    ]
    for P, tag in hard_cases:
        k = 13 - len(P)
        for Vmax_t, Dc in [(1009, 12), (1009, 60), (10007, 12), (10007, 400)]:
            offs = [0] + sorted(random.sample(range(1, max(Dc, k + 1) + 1), k - 1))
            v = make_shape(P, offs, Vmax_t)
            if len(v) != 13:
                continue
            # global witness search on a fine tau grid + local refinement
            res = 400000
            taus = (np.arange(res) + 0.5) / res
            va = np.array(v, dtype=np.float64)
            # process in chunks to bound memory
            best, bt = -1.0, 0.0
            for s in range(0, res, 100000):
                tt = taus[s:s + 100000]
                ph = np.outer(tt, va)
                d = np.abs(ph - np.rint(ph)).min(axis=1)
                i = int(d.argmax())
                if d[i] > best:
                    best, bt = float(d[i]), float(tt[i])
            print(f"  P={P} ({tag}) Vmax={Vmax_t} Dc={Dc}: M >= {best:.5f} "
                  f"({'OK' if best >= 1/14 else '!! BELOW 1/14 !!'}) at tau~{bt:.6f} "
                  f"~ {F(bt).limit_denominator(50)}")

    print()
    print("=" * 78)
    print("PART 4 -- verdict inputs")
    print("=" * 78)
    print("  If Part 1 confirms the (7-q)/(14q) boundary (with a safety factor), the")
    print("  q-window lemma gives EXPLICIT unconditional reach for every mixed shape")
    print("  with a surviving q and tight-enough cluster -- replacing Part A's sketch")
    print("  on that domain.  The honest 1/14 residual = kill-all-P shapes (Part 2 count)")
    print("  + wide clusters (D_c/Vmax above the boundary; recurse via peel/THM-608).")
