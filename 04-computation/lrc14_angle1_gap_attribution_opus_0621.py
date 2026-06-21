#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_angle1_gap_attribution_opus_0621.py  (opus-2026-06-21, ANGLE 1, part 5)

EXACT GAP ATTRIBUTION: where does consec's surplus over its nearest rivals live?

From part-4: consec is the global max but NOT a unique local max; the swap-landscape
is rugged.  So a hill-climb proof is dead.  We instead dissect the SURPLUS
   Delta(E) = sum_a W_a(consec) - sum_a W_a(E)  > 0
into per-resonance and per-sector pieces, for the NEAREST rivals (largest W < consec),
to find the conserved quantity consec optimizes.

KEY MODEL FACT (parts 1-3): for full-residue E, at each resonance a:
  - same y=0 start config {0,0,1,2,3,4,5,6} (Z/7, residue 0 doubled);
  - W_a = total covered y-measure in [-1/14,1/14];
  - the contiguous WINDOW through center has length T+(a)+T-(a) where T+,T- are the
    first-death times right/left, each = 1/(7 * v_bind) for the BINDING clock speed.

NEW HYPOTHESES (data-driven):
 (P1) BINDING-SPEED LAW: T+(a) and T-(a) are each 1/(7*b) for an integer b = the
      speed (|e|) of the binding clock on that side.  consec's window is governed by
      the TWO LARGEST speeds present (6 and 7).  Conjecture: sum of windows is
      controlled by sum_a 1/(7 b_a^+) + 1/(7 b_a^-) and consec MINIMIZES the binding
      speeds b -> MAXIMIZES the windows.  Test the binding-speed law exactly.
 (P2) DISCONNECTED-ARC BONUS: W_a - window = extra covered arcs away from center.
      Where do these come from?  consec has a BIG disconnected bonus at a=3,4,6.
      Conjecture: the bonus is consec's secret weapon -- it loses some windows but
      gains via disconnected arcs.  Quantify per rival.
 (P3) THE CONSERVED SUM: is sum_a (T+(a)+T-(a)) [windows only] ALSO maximized by
      consec, or is it only the TOTAL (windows + disconnected) that consec wins?
      If consec loses on windows-only but wins on total, the disconnected arcs are
      essential -> the survival-window picture is INCOMPLETE (a sharp structural
      statement; tells us what a proof must capture).
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def sector_of_point(e, a, y):
    pos = F(e*a) + F(7*e)*y
    return (pos.numerator // pos.denominator) % 7

def covered_all_at(E, a, y):
    return len({sector_of_point(e, a, y) for e in E}) == 7

def breakpoints(E, a):
    half = F(1, 14); bps = {F(0), -half, half}
    for e in E:
        if e == 0: continue
        lo_val = F(7*e)*(-half) + F(e*a); hi_val = F(7*e)*(half) + F(e*a)
        lo_i = min(lo_val, hi_val); hi_i = max(lo_val, hi_val)
        m = lo_i.numerator // lo_i.denominator
        while m <= hi_i.numerator // hi_i.denominator + 1:
            y = F(m - e*a, 7*e)
            if -half <= y <= half: bps.add(y)
            m += 1
    return sorted(bps)

def W_a_total(E, a):
    bps = breakpoints(E, a); tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        if covered_all_at(E, a, (lo+hi)/2): tot += hi - lo
    return tot

def window_TpTm(E, a):
    half = F(1, 14); bps = breakpoints(E, a); ivals = list(zip(bps, bps[1:]))
    Tp = F(0)
    for lo, hi in ivals:
        if hi <= 0: continue
        lo2 = max(lo, F(0))
        if covered_all_at(E, a, (lo2+hi)/2): Tp = hi
        else:
            if lo2 == F(0): Tp = F(0)
            break
    Tp = min(Tp, half)
    Tm = F(0)
    for lo, hi in reversed(ivals):
        if lo >= 0: continue
        hi2 = min(hi, F(0))
        if covered_all_at(E, a, (lo+hi2)/2): Tm = -lo
        else:
            if hi2 == F(0): Tm = F(0)
            break
    Tm = min(Tm, half)
    return Tp, Tm

def is_full_residue(E): return frozenset(e % 7 for e in E) == frozenset(range(7))
def consec(k): return list(range(k))
def measS7(E): return sum(W_a_total(E, a) for a in range(1, 7))
def window_sum(E): return sum(sum(window_TpTm(E, a)) for a in range(1, 7))
def disconn_sum(E): return measS7(E) - window_sum(E)

if __name__ == "__main__":
    print("="*80)
    print("GAP ATTRIBUTION: consec surplus dissected (windows vs disconnected arcs)")
    print("="*80)
    k = 8; C = consec(k)
    msC = measS7(C); winC = window_sum(C); disC = disconn_sum(C)
    print(f"\nconsec={C}")
    print(f"  total W = {float(msC):.6f}   windows = {float(winC):.6f}   disconnected = {float(disC):.6f}")

    # nearest rivals
    span = 14
    shapes = [(0,)+c for c in itertools.combinations(range(1, span+1), k-1)]
    shapes = [E for E in shapes if is_full_residue(E) and E != tuple(C)]
    ranked = sorted(shapes, key=lambda E: -measS7(list(E)))[:6]

    # (P3) does consec win on WINDOWS-ONLY too, or only on TOTAL?
    print(f"\n[P3] windows-only vs total for top rivals (does consec need disconnected arcs?):")
    print(f"  {'shape':32s} {'total':>9} {'windows':>9} {'disconn':>9}")
    print(f"  {'consec':32s} {float(msC):>9.6f} {float(winC):>9.6f} {float(disC):>9.6f}")
    win_beats = 0
    for E in ranked:
        E = list(E)
        print(f"  {str(E):32s} {float(measS7(E)):>9.6f} {float(window_sum(E)):>9.6f} {float(disconn_sum(E)):>9.6f}")
    # check across whole stratum: is consec max on windows-only?
    best_win = max(shapes+[tuple(C)], key=lambda E: window_sum(list(E)))
    nwin_beat = sum(1 for E in shapes if window_sum(list(E)) > winC)
    print(f"\n  consec window_sum max over stratum? argmax={'consec' if best_win==tuple(C) else best_win}; "
          f"#shapes beating consec on windows-only = {nwin_beat}")

    # (P1) binding-speed law: per a, is T+ = 1/(7 b) for integer b?
    print(f"\n[P1] BINDING-SPEED LAW for consec: T+, T- as 1/(7*b):")
    for a in range(1, 7):
        Tp, Tm = window_TpTm(C, a)
        bp = (1/(7*Tp)) if Tp>0 else None
        bm = (1/(7*Tm)) if Tm>0 else None
        print(f"  a={a}: T+={Tp}=1/(7*{bp})  T-={Tm}=1/(7*{bm})")

    # (P2) disconnected-arc bonus per a
    print(f"\n[P2] DISCONNECTED-arc bonus per a (W_a - window_a) for consec:")
    for a in range(1, 7):
        Wt = W_a_total(C, a); Tp, Tm = window_TpTm(C, a); win = Tp+Tm
        print(f"  a={a}: W_a={float(Wt):.6f}  window={float(win):.6f}  bonus={float(Wt-win):.6f}")

    # summary: which component dominates the surplus over the nearest rival?
    R = list(ranked[0])
    print(f"\nNEAREST RIVAL {R}:")
    print(f"  total surplus  Delta = {float(msC-measS7(R)):.6f}")
    print(f"  windows surplus      = {float(winC-window_sum(R)):.6f}")
    print(f"  disconnect surplus   = {float(disC-disconn_sum(R)):.6f}")
    print(f"\n  => consec's edge is "
          + ("mostly WINDOWS" if (winC-window_sum(R))>(disC-disconn_sum(R)) else "mostly DISCONNECTED ARCS"))
