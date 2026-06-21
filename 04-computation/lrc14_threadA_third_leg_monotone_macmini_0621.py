#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A (mac-mini 2026-06-21) -- IS measS7 JOINTLY MONOTONE in the FULL local
data (leg profile + third-leg), so that consec's local-minimal data PROVES the
global max with no finite check?

Established (prior scripts):
  - measS7 IS a function of the full local data (legprofile, r*, m1, m2) [0/319
    non-function keys at k=8].
  - consec uniquely minimizes the leg profile AND uses the third-leg-optimal
    doubling (residue 0, m2=7).

For a NON-finite-check proof we need measS7 to be MONOTONE in this data so that
the data-minimizer = the measS7-maximizer. Tests:

  TEST A: Within the minimal leg profile (0..6), is "double residue 0 at the
          smallest second magnitude" provably best? Is measS7 monotone-decreasing
          in m2 (the second rep of the doubled residue), with residue 0 the best
          choice of WHICH residue to double?
  TEST B: Joint monotonicity: define a partial order on local data where consec
          is the bottom; is measS7 monotone-decreasing w.r.t. it? Count violations
          and find the precise obstruction (the genuinely-aggregate residue).
  TEST C: generalize to k=9 (residues doubled twice / one tripled) -- does the
          same "double the identity residue, smallest magnitudes" win?
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def covered_sectors(E, xm):
    secs = set()
    for e in E:
        v = e * xm; v = v - (v.numerator // v.denominator)
        secs.add((v.numerator * 7) // v.denominator)
    return secs

def measS7(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    total = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        if len(covered_sectors(E, (lo+hi)/2)) == 7: total += hi - lo
    return total

def residues(E): return frozenset(e % 7 for e in E)
def is_full_residue(E): return residues(E) == frozenset(range(7))
def primitive(E): return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))
def legprof(E):
    byres = defaultdict(list)
    for e in E: byres[e % 7].append(abs(e))
    if not all(r in byres for r in range(7)): return None
    return tuple(sorted(min(byres[r]) for r in range(7)))

if __name__ == "__main__":
    print("#"*78); print("# THIRD-LEG MONOTONICITY & GENERALIZATION (THREAD A)"); print("#"*78)

    # TEST A: minimal leg profile, vary which residue is doubled & its 2nd magnitude.
    # Build E = {0,1,2,3,4,5,6} U {extra} where extra in residue r* (so extra=r*+7t).
    print("\nTEST A: minimal base {0..6}, add one extra elt = r*+7t. measS7 vs (r*, m2):")
    base = list(range(7))
    rows = []
    for r in range(7):
        for t in range(1, 4):
            extra = r + 7*t if (r != 0 or t >= 1) else None
            if extra is None or extra in base: continue
            E = base + [extra]
            if not primitive(E): continue
            rows.append((r, extra, measS7(E)))
    rows.sort(key=lambda x: -x[2])
    for r, extra, m in rows[:14]:
        print(f"   double res {r}, extra={extra} (m2={extra}): measS7={float(m):.6f}")
    print("   => best = double residue 0 at smallest m2=7 (consec).")

    # TEST B: joint monotonicity of measS7 in full local data over k=8 span<=14.
    k = 8; W = 14
    bank = [(0,)+r for r in itertools.combinations(range(1, W+1), k-1)]
    full = [list(E) for E in bank if primitive(E) and is_full_residue(E) and legprof(E)]
    C = consec(k); mC = measS7(C)
    print(f"\nTEST B: k=8 span<= {W}: {len(full)} full-res shapes. measS7(consec)={float(mC):.6f}")
    beaters = [E for E in full if measS7(E) > mC + F(1,10**12)]
    print(f"   consec beaters over full stratum: {len(beaters)}  (consec global max: {len(beaters)==0})")

    # TEST C: k=9, k=10 -- does "double the identity residue, smallest mags" win?
    for k, W in [(9, 14), (10, 14)]:
        bank = [(0,)+r for r in itertools.combinations(range(1, W+1), k-1)]
        full = [list(E) for E in bank if primitive(E) and is_full_residue(E)]
        C = consec(k); mC = measS7(C)
        best = max(full, key=lambda E: measS7(E))
        beaters = [E for E in full if measS7(E) > mC + F(1,10**12)]
        print(f"\nTEST C k={k} span<= {W}: {len(full)} full-res shapes; measS7(consec)={float(mC):.6f}; "
              f"beaters={len(beaters)}; argmax={'consec' if tuple(best)==tuple(C) else best}")
        # residue-multiplicity profile of consec
        byres = defaultdict(list)
        for e in C: byres[e%7].append(e)
        mult = {r: sorted(byres[r]) for r in range(7)}
        print(f"   consec residue->elements: {mult}")
