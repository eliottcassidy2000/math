#!/usr/bin/env python3
"""
lrc14_spread_core_floor_boxeph_S19.py  (HYP-6132 follow-on, boxeph-2026-07-12-S19)

Pin the SPREAD-CORE transport floors exactly by exhaustive census at small Vmax:

  F12(Vmax) = min { M(H) : H a 12-subset of [1..Vmax], longest-AP(H) <= 7 }
  F11(Vmax) = same over 11-subsets
  ...and the TRANSPORTABLE variants (2H covers every d in {2..14}), which are the cores
  that actually ship as v_c = 2c*H u {deltas} (primitive + DC + spread at every scale).

Why: MISTAKE-140 / death-star-S14. The compressed (non-spread) floor is 1/13 sharp (V_L).
The spread stratum's floor was only bounded (<= 1/11 via H*). The transport turns any
spread core floor into a scale-free bound on min M(spread DC, diam >= D):
  j=1 (12-core + 1 detuned):  M(v_c) in [1/13, F12] (detuned dispatch g=2, monad)
  j=2 (11-core + 2 detuned):  M(v_c) in [d=2 dispatch floor, F11] (THM-678 territory)

Method: exact breakpoint evaluator with a CAP -- abort a core as soon as any candidate
time certifies margin >= cap (it cannot improve the running minimum). Small moduli first.
Pure Python, exact integers.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import sys

def longest_ap_le(vv, cap=7):
    s = set(vv)
    for i in range(len(vv)):
        for j in range(i + 1, len(vv)):
            d = vv[j] - vv[i]
            L = 2
            x = vv[j] + d
            while x in s: L += 1; x += d
            x = vv[i] - d
            while x in s: L += 1; x -= d
            if L > cap: return False
    return True

def dc_supporting(vv):
    # 2H covers every d in 2..14  <=>  for each d, some h with d | 2h
    for d in range(2, 15):
        ok = False
        for h in vv:
            if (2 * h) % d == 0: ok = True; break
        if not ok: return False
    return True

def M_capped(v, capn, capd):
    """Exact M(v) as (num, den) if M < cap, else None (aborted: M >= cap).
    Moduli: q <= 60 first, then pair-sums (i<=j) and |diffs| ascending."""
    qs = sorted(set(range(2, 61)))
    big = set()
    n = len(v)
    for i in range(n):
        for j in range(i, n):
            big.add(v[i] + v[j])
            if v[i] != v[j]: big.add(abs(v[i] - v[j]))
    qs = qs + sorted(q for q in big if q > 60)
    mb, qb = 0, 1
    for q in qs:
        if q < 2: continue
        for p in range(1, q // 2 + 1):
            m = q
            for x in v:
                r = (x * p) % q
                if r > q - r: r = q - r
                if r < m:
                    m = r
                    if m * qb <= mb * q: break  # cannot beat current best
            if m * qb > mb * q:
                mb, qb = m, q
                if mb * capd >= capn * qb:  # margin >= cap => M >= cap, abort
                    return None
    return (mb, qb)

def census(k, vmax, need_dc, cap):
    """Min M over spread k-subsets of [1..vmax]; returns (min_frac, argmin, n_evaluated)."""
    capn, capd = cap
    best = None; arg = None; nev = 0; nsp = 0
    for c in combinations(range(1, vmax + 1), k):
        if not longest_ap_le(c, 7): continue
        nsp += 1
        if need_dc and not dc_supporting(c): continue
        nev += 1
        r = M_capped(list(c), capn, capd)
        if r is not None:
            m, q = r
            if best is None or m * best[1] < best[0] * q:
                best = (m, q); arg = c
                capn, capd = m, q  # tighten cap to current record
    return best, arg, nev, nsp

def main():
    vmax = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    print(f"spread-core floor census, Vmax={vmax} (exact, capped evaluator)")
    print(f"reference: 1/11 = {float(F(1,11)):.5f} (H*), 1/12 = {float(F(1,12)):.5f}, 1/13 = {float(F(1,13)):.5f}")
    for k, cap in ((12, (1, 10)), (11, (1, 10))):
        for need_dc in (True, False):
            best, arg, nev, nsp = census(k, vmax, need_dc, cap)
            tag = "transportable (2H DC-covering)" if need_dc else "any spread core"
            if best:
                m, q = best
                print(f"k={k} {tag:32s}: min M = {m}/{q} = {m/q:.6f}  at {arg}   [{nev} evaluated / {nsp} spread]")
            else:
                print(f"k={k} {tag:32s}: NO core below cap {cap[0]}/{cap[1]}   [{nev} evaluated / {nsp} spread]")
    print("\ntransport consequence: min M(spread primitive DC, diam >= D) <= F12 for every D;")
    print("with two detuned slots (j=2), <= F11 for every D (THM-678 d=2 dispatch supplies the floor).")

if __name__ == '__main__':
    main()
