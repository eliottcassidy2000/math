#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_angle1_exchange_stability_opus_0621.py  (opus-2026-06-21, ANGLE 1, part 4)

THE EXCHANGE-STABILITY CERTIFICATE for LAYER 3 (consec maximizes sum_a W_a).

Part-3 found: consec is a STRICT LOCAL MAX of sum_a W_a under all residue-preserving
single-clock swaps (32/32 swaps strictly drop, 0 ties, 0 rises) for k=8.

This script tests whether exchange-stability is a faithful CERTIFICATE of global
maximality, by checking:

 (E1) UNIQUENESS: is consec the ONLY exchange-stable full-residue shape (within a
      span box)?  A swap-move connects shapes differing in one clock (same residue).
      If consec is the unique local max and the swap-graph is connected with no other
      local maxima, then local=global (a discrete Morse / hill-climbing certificate).
      We compute ALL local maxima of sum_a W_a under residue-preserving swaps.

 (E2) HILL-CLIMB GLOBALITY: from EVERY full-residue shape, does greedy
      residue-preserving swap-ascent reach consec (up to scaling)?  If yes for all
      starts, consec is the unique attractor => global max on the stratum.

 (E3) k=9, k=10 generalization: repeat E1/E2 to confirm the certificate is not a
      k=8 accident.

 (E4) MONOTONE-ALONG-PATH: along the ascent path, is sum_a W_a strictly increasing
      each step?  (needed for hill-climb to be valid; report any plateau.)

All EXACT Fractions.  Honest: a single non-consec local max REFUTES the clean
certificate (forces a global argument); zero => strong VERIFIED certificate.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce, lru_cache
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def sector_of_point(e, a, y):
    pos = F(e*a) + F(7*e)*y
    return (pos.numerator // pos.denominator) % 7

def covered_all_at(E, a, y):
    return len({sector_of_point(e, a, y) for e in E}) == 7

def W_a_total(E, a):
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
    bps = sorted(bps); tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        if covered_all_at(E, a, (lo+hi)/2): tot += hi - lo
    return tot

@lru_cache(maxsize=None)
def measS7(E):  # E a sorted tuple
    return sum(W_a_total(list(E), a) for a in range(1, 7))

def is_full_residue(E): return frozenset(e % 7 for e in E) == frozenset(range(7))
def consec(k): return tuple(range(k))
def primitive(E):
    g = reduce(gcd, [abs(e) for e in E if e != 0], 0)
    return g == 1

def swap_neighbors(E, span):
    """all full-residue shapes reachable by changing ONE clock to another rep of
       its residue, within [0,span], distinct, size preserved."""
    E = tuple(sorted(E)); Eset = set(E); k = len(E); out = []
    for e in E:
        r = e % 7
        for ep in range(0, span+1):
            if ep % 7 != r or ep == e or ep in Eset: continue
            E2 = tuple(sorted((Eset - {e}) | {ep}))
            if len(E2) == k and is_full_residue(E2):
                out.append(E2)
    return out

def local_maxima(k, span):
    shapes = [(0,)+c for c in itertools.combinations(range(1, span+1), k-1)]
    shapes = [E for E in shapes if is_full_residue(E)]
    maxima = []
    for E in shapes:
        m = measS7(E)
        nbrs = swap_neighbors(E, span)
        if all(measS7(nb) <= m for nb in nbrs):
            # strict local max: no neighbor >= (allow ties = scalings)
            strictly = all(measS7(nb) < m for nb in nbrs if nb != E)
            maxima.append((E, m, strictly))
    return shapes, maxima

def hillclimb(E, span):
    """greedy ascent via residue-preserving swaps; return endpoint + monotone flag."""
    E = tuple(sorted(E)); path = [E]; mono = True
    while True:
        m = measS7(E)
        nbrs = swap_neighbors(E, span)
        best = E; bm = m
        for nb in nbrs:
            mn = measS7(nb)
            if mn > bm: bm = mn; best = nb
        if best == E: break
        if bm <= m: mono = False
        E = best; path.append(E)
    return E, path, mono

if __name__ == "__main__":
    print("="*80)
    print("EXCHANGE-STABILITY CERTIFICATE for LAYER 3")
    print("="*80)

    for k, span in [(8, 14), (9, 15), (10, 16)]:
        print(f"\n{'#'*70}\n### k={k}, span<= {span}\n{'#'*70}")
        C = consec(k); msC = measS7(C)
        print(f"  consec={C}  sum W_a={float(msC):.6f}={msC}")

        # (E1) all local maxima
        shapes, maxima = local_maxima(k, span)
        print(f"\n  [E1] full-residue shapes={len(shapes)}; local maxima (no swap-nbr strictly higher):")
        # group maxima by primitive class
        prim_maxima = [(E,m,st) for (E,m,st) in maxima]
        print(f"       #local maxima (incl scalings) = {len(maxima)}")
        # which are scalings of consec?
        def is_scaling_of_consec(E):
            g = reduce(gcd, [e for e in E if e>0], 0)
            return tuple(e//g for e in E) == C if g else False
        nonconsec = [(E,m,st) for (E,m,st) in maxima if not is_scaling_of_consec(E)]
        consec_scalings = [(E,m,st) for (E,m,st) in maxima if is_scaling_of_consec(E)]
        print(f"       consec & its scalings among maxima: {len(consec_scalings)}  "
              f"e.g. {[E for E,_,_ in consec_scalings][:4]}")
        print(f"       NON-consec local maxima: {len(nonconsec)}")
        if nonconsec:
            print(f"       *** consec is NOT the unique local max ***")
            for E,m,st in sorted(nonconsec, key=lambda t:-t[1])[:5]:
                print(f"          {E}: W={float(m):.6f} (strict={st})  ratio to consec={float(m/msC):.4f}")
        else:
            print(f"       => EVERY local max is a scaling of consec.")

        # (E2)/(E4) hill-climb globality from all starts
        reached = 0; not_reached = 0; nonmono = 0; bad_ends = set()
        for E in shapes:
            end, path, mono = hillclimb(E, span)
            if is_scaling_of_consec(end): reached += 1
            else:
                not_reached += 1; bad_ends.add(end)
            if not mono: nonmono += 1
        print(f"\n  [E2] hill-climb (residue-preserving swap-ascent) from all {len(shapes)} starts:")
        print(f"       reach consec (up to scaling): {reached}/{len(shapes)}")
        print(f"       stuck elsewhere: {not_reached}" + (f"  ends={list(bad_ends)[:5]}" if bad_ends else ""))
        print(f"       non-monotone ascents: {nonmono}")
        if not_reached == 0 and len(nonconsec) == 0:
            print(f"       => CERTIFICATE HOLDS at k={k}: consec is the UNIQUE swap-local-max")
            print(f"          and the unique hill-climb attractor => global max on stratum.")
