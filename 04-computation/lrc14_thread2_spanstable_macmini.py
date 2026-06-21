#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_thread2_spanstable_macmini.py  (mac-mini-2026-06-21, THREAD 2)

The bounded-span scan (lrc14_thread2_crude_bound) showed, at k=12:
   max measS7 = 0.64518   < B_A=0.84094  < B_C=0.74143  < cap_12 = 6/7 = 0.85714.
But a TRUE bound must hold for ALL spans (LRC realizable k<=13, but the abstract
problem is unbounded).  This script answers the load-bearing question:

  (1) Does the crude bound B_C (pair-intersection)  -- and B_A (single-sector) --
      keep GROWING with span at k=12, and if so does it ever CROSS cap_12 = 6/7?
      Sweep span s = 11..S, exact, primitive E={0}+subset of [1,s], |E|=12.

  (2) The SPAN-FREE single-element fact:  for a single nonzero integer e and uniform
      x in [0,1], frac(e x) is EXACTLY uniform, so each sector j is hit by e with
      probability EXACTLY 1/7.  We use this to derive a span-free crude bound:
         measS7(E) <= P(sector j hit by E)  and  P(sector j missed by E) >= ??? .
      We MEASURE P(sector j missed) for the worst j over the stratum to see whether
      it stays >= 1/7  (which would give measS7 <= 6/7 span-free).

  (3) The decorrelated/iid LIMIT:  as the set spreads (dissociated), the elements act
      as independent uniform Z/7 samples; measS7 -> 7! S(k,7)/7^k (the iid surjection
      prob), which is TINY (<< 6/7).  So large span DRIVES measS7 DOWN, not up.
      We confirm the iid limit and that the bound's growth saturates BELOW cap.

Exact Fractions throughout; reuses the VERIFIED measS7.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd, comb, factorial
sys.stdout.reconfigure(line_buffering=True)

P = 7
CAP12 = F(6, 7)

def breakpoints(E):
    bps = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * abs(e) + 1):
            bps.add(F(m, 7 * abs(e)))
    return sorted(b for b in bps if 0 <= b <= 1)

def cells_with_sectors(E):
    bps = breakpoints(E); out = []
    for a in range(len(bps) - 1):
        lo, hi = bps[a], bps[a + 1]; mid = (lo + hi) / 2
        secs = frozenset(int(((e * mid) % 1) * 7) for e in E)
        out.append((hi - lo, secs))
    return out

def measS7_from_cells(cells):
    return sum(L for (L, s) in cells if len(s) == 7)

def P_sector_hit(cells, j):  return sum(L for (L, s) in cells if j in s)
def P_sector_miss(cells, j): return sum(L for (L, s) in cells if j not in s)
def P_pair_both_hit(cells, i, j): return sum(L for (L, s) in cells if i in s and j in s)

def bound_A(cells): return min(P_sector_hit(cells, j) for j in range(7))
def bound_C(cells): return min(P_pair_both_hit(cells, i, j) for i in range(7) for j in range(i+1, 7))
def worst_miss(cells): return max(P_sector_miss(cells, j) for j in range(7))  # largest single-sector miss prob

def is_primitive(E):
    g = 0
    for e in E: g = gcd(g, abs(e))
    return g == 1

def stratum_max(k, span):
    mxm = (F(0), None); mxA = F(0); mxC = F(0); min_worstmiss = F(2); argwm = None
    for combo in itertools.combinations(range(1, span + 1), k - 1):
        E = (0,) + combo
        if not is_primitive(E): continue
        cells = cells_with_sectors(E)
        m = measS7_from_cells(cells)
        if m > mxm[0]: mxm = (m, E)
        a = bound_A(cells); c = bound_C(cells); wm = worst_miss(cells)
        if a > mxA: mxA = a
        if c > mxC: mxC = c
        # the SPAN-FREE missing-sector check wants the SMALLEST worst-miss over E
        # (i.e. the E that comes CLOSEST to having every sector almost-surely hit).
        if wm < min_worstmiss: min_worstmiss = wm; argwm = E
    return mxm, mxA, mxC, min_worstmiss, argwm

def stirling2(n, k): return sum((-1)**(k-j)*comb(k,j)*j**n for j in range(k+1))//factorial(k)

if __name__ == "__main__":
    k = 12
    print("=" * 92)
    print(f"THREAD 2 span-stability at k={k}; cap_12 = 6/7 = {float(CAP12):.5f}")
    print("=" * 92)
    print("iid surjection limit 7! S(12,7)/7^12 =",
          F(factorial(7)*stirling2(12,7), 7**12), "=",
          float(F(factorial(7)*stirling2(12,7), 7**12)),
          "(spread -> this, << cap)")
    print()
    print("%-6s %-10s %-10s %-10s %-14s %-10s" %
          ("span", "maxMeas", "B_A", "B_C", "min_worstmiss", "B_C<=6/7?"))
    for span in range(11, 21):
        (m, argm), A, C, wm, argwm = stratum_max(k, span)
        # only show progress; min_worstmiss>=1/7 would give the span-free measS7<=6/7
        print("%-6d %-10.5f %-10.5f %-10.5f %-14.5f %-10s" %
              (span, float(m), float(A), float(C), float(wm), str(C <= CAP12)))
    print()
    print("min_worstmiss >= 1/7 = 0.142857 would PROVE measS7 <= 6/7 span-free (single-sector).")
    print("DONE.")
