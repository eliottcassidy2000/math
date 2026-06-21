#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_thread2_crude_bound_macmini.py  (mac-mini-2026-06-21, THREAD 2)

GOAL (corrected architecture, HYP-2778):  find a SIMPLE, NON-aggregate upper bound
   B(k)   with    max_E measS7(E)  <=  B(k)  <=  cap_k     for all k >= 12.
It need NOT be tight; it must just stay <= cap_k for k>=12 and NOT rely on consec-max.

caps (exact, from the 4-thread audit ALL-PASS):
   cap_8..cap_13 = 2243/5880, 1979/4004, 55/91, 66/91, 6/7, 1.

measS7 (the VERIFIED exact code): Z/7 vertex coloring color(e,x)=floor(7 frac(e x));
   measS7(E) = Lebesgue meas { x in [0,1] : { floor(7 frac(e x)) : e in E } = all of Z/7 }.

This script:
 (1) reproduces the EXACT measS7 (Fraction Lebesgue) verbatim from the landscape script;
 (2) computes the exact stratum-max of measS7 over k-subsets {0}∪combo, combo in [1,span], for k=12..16;
 (3) tests several CRUDE upper bounds, each computed EXACTLY per E, and reports whether
     B(E) >= measS7(E) (validity) and whether max_E B(E) <= cap_k (clears cap):
        (A) single-sector intersection:   measS7 <= min_j P(sector j hit)
        (B) Bonferroni-2 on the complement:
              measS7 = 1 - P(∪_j miss_j)
                     <= 1 - [ Σ_j P(miss_j)  -  Σ_{i<j} P(miss_i ∩ miss_j) ]   (Bonferroni lower bd on the union)
        (C) pair-intersection:   measS7 <= min_{i<j} P(sector i hit AND sector j hit)
        (D) triple-of-complement / a constant decorrelated value + resonance R<=(20/7)/p.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd, comb, factorial
sys.stdout.reconfigure(line_buffering=True)

P = 7
CAPS = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91),
        11: F(66, 91), 12: F(6, 7), 13: F(1, 1)}

# ----------------------------------------------------------------------------
# EXACT measS7 (verbatim mechanics from lrc14_measS7_extremal_landscape S6)
# ----------------------------------------------------------------------------
def breakpoints(E):
    bps = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * abs(e) + 1):
            bps.add(F(m, 7 * abs(e)))
    return sorted(b for b in bps if 0 <= b <= 1)

def sector_at(e, mid):
    return int(((e * mid) % 1) * 7)

def measS7(E):
    bps = breakpoints(E); tot = F(0)
    for a in range(len(bps) - 1):
        lo, hi = bps[a], bps[a + 1]; mid = (lo + hi) / 2
        if len(set(sector_at(e, mid) for e in E)) == 7:
            tot += hi - lo
    return tot

# Per-cell sector profile: returns list of (length, frozenset_of_sectors_present)
def cell_profile(E):
    bps = breakpoints(E); cells = []
    for a in range(len(bps) - 1):
        lo, hi = bps[a], bps[a + 1]; mid = (lo + hi) / 2
        secs = frozenset(sector_at(e, mid) for e in E)
        cells.append((hi - lo, secs))
    return cells

# ----------------------------------------------------------------------------
# CRUDE bound candidates, each EXACT, each computed from the cell profile.
# ----------------------------------------------------------------------------
def P_sector_hit(cells, j):
    return sum(L for (L, secs) in cells if j in secs)

def P_sector_miss(cells, j):
    return sum(L for (L, secs) in cells if j not in secs)

def P_pair_both_miss(cells, i, j):
    return sum(L for (L, secs) in cells if i not in secs and j not in secs)

def P_pair_both_hit(cells, i, j):
    return sum(L for (L, secs) in cells if i in secs and j in secs)

def bound_A_single(cells):
    # measS7 <= min_j P(sector j hit)
    return min(P_sector_hit(cells, j) for j in range(7))

def bound_B_bonf2(cells):
    # measS7 = 1 - P(union miss). Bonferroni-2 LOWER bound on the union:
    #   P(union) >= S1 - S2 .  =>  measS7 <= 1 - (S1 - S2) = 1 - S1 + S2.
    S1 = sum(P_sector_miss(cells, j) for j in range(7))
    S2 = sum(P_pair_both_miss(cells, i, j) for i in range(7) for j in range(i + 1, 7))
    return F(1) - S1 + S2

def bound_C_pair(cells):
    # measS7 <= min_{i<j} P(sector i hit AND sector j hit)
    return min(P_pair_both_hit(cells, i, j) for i in range(7) for j in range(i + 1, 7))

def is_primitive(E):
    g = 0
    for e in E: g = gcd(g, abs(e))
    return g == 1

# ----------------------------------------------------------------------------
# Stratum scan: E = {0} ∪ combo, combo subset of [1, span], |E| = k.
# ----------------------------------------------------------------------------
def stratum_scan(k, span, primitive_only=True):
    best_meas = (F(0), None)
    # track, over all E in stratum, the per-E bound values and validity
    maxB = {'A': F(0), 'B': F(0), 'C': F(0)}
    valid = {'A': True, 'B': True, 'C': True}
    invalid_witness = {'A': None, 'B': None, 'C': None}
    nE = 0
    for combo in itertools.combinations(range(1, span + 1), k - 1):
        E = (0,) + combo
        if primitive_only and not is_primitive(E):
            continue
        nE += 1
        cells = cell_profile(E)
        m = measS7(E)
        if m > best_meas[0]:
            best_meas = (m, E)
        bA = bound_A_single(cells)
        bB = bound_B_bonf2(cells)
        bC = bound_C_pair(cells)
        for name, b in (('A', bA), ('B', bB), ('C', bC)):
            if b > maxB[name]:
                maxB[name] = b
            if b < m:  # bound must be >= measS7
                if valid[name]:
                    valid[name] = False
                    invalid_witness[name] = (E, m, b)
    return best_meas, maxB, valid, invalid_witness, nE

if __name__ == "__main__":
    print("=" * 96)
    print("THREAD 2: crude bound B(k) with max measS7 <= B(k) <= cap_k for k>=12")
    print("=" * 96)
    # quick sanity: reproduce known binding value
    print("\nSANITY: measS7(consec_10) =", measS7(list(range(10))), "=", float(measS7(list(range(10)))),
          " (expect 8899/17640 = 0.5045)")

    # the stratum for thread-2 is bounded span; LRC(14) => relevant span is bounded.
    # We scan span up to a feasible cap per k. The cap is meaningful for k>=12.
    SCAN = {12: 16, 13: 16, 14: 16, 15: 17, 16: 18}
    print("\n%-3s %-6s %-9s %-10s | %-9s %-9s %-9s | %-9s %-9s %-9s" %
          ("k", "span", "maxMeas", "cap_k", "maxB_A", "maxB_B", "maxB_C",
           "A<=cap?", "B<=cap?", "C<=cap?"))
    for k in [12, 13, 14, 15, 16]:
        span = SCAN[k]
        cap = CAPS.get(k, F(1))
        (mx, argmx), maxB, valid, inv, nE = stratum_scan(k, span)
        print("%-3d %-6d %-9.5f %-10s | %-9.5f %-9.5f %-9.5f | %-9s %-9s %-9s  (nE=%d)" %
              (k, span, float(mx), str(cap),
               float(maxB['A']), float(maxB['B']), float(maxB['C']),
               str(maxB['A'] <= cap), str(maxB['B'] <= cap), str(maxB['C'] <= cap), nE))
        for name in ('A', 'B', 'C'):
            if not valid[name]:
                E, m, b = inv[name]
                print(f"    !!! bound {name} INVALID (b<measS7) at E={E}: measS7={float(m):.5f} > b={float(b):.5f}")
        print(f"    argmax measS7 = {argmx}  measS7={mx}={float(mx):.5f}")
    print("\nDONE.")
