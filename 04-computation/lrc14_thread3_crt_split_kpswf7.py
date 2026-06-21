#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD 3 -- ANGLE 1: CRT / composite-14 2-adic x 7-adic PRODUCT split of measS7.

LRC(14): a config E of speeds; the lonely point fails (is covered) if every x in (0,1)
hits all residues. The standard "all 6 nonzero residues mod 7 covered" measure measS7(E)
(meas of x in [0,1] with {e*x mod 1 -> sector in Z/7} hitting all of {1..6}) is the 7-adic
APEX of the staircase. 14 = 2 * 7. The polynomial method (Rosenfeld/Fermat) needs k+1 PRIME;
14 composite is exactly why LRC(14) is open.

QUESTION: does measS7 (the 7-part) FACTOR through a 2-adic (doubling / even-odd) part and a
7-adic (apex) part, so that consec-max separates into two coprime sub-problems recombined by CRT?

Concretely the cover is "hit all of Z/14 \ {0} = (Z/2 x Z/7) minus 0". By CRT a residue mod 14
is a pair (parity, residue mod 7). Probe whether:
  (P1) the covered-x set for mod-14 = (covered for mod-2 part) AND (covered for mod-7 part)?
       i.e. does meas14(E) = measParity(E) * measS7(E)  as a PRODUCT (independence)?
  (P2) is consec the maximiser of EACH factor separately?
  (P3) the 2-adic part: doubling map x->2x; even/odd residue. Is the parity-cover trivial
       (always full) so the whole content is the 7-adic apex? (=> CRT split is trivial-but-clean)

5-min probe: exact rational meas of the parity-cover and the mod-7 cover, and the joint mod-14
cover, over many configs; test the product law and per-factor consec-max.
"""
import sys, itertools, random
from fractions import Fraction as F
from functools import reduce
from math import gcd
sys.stdout.reconfigure(line_buffering=True)

def breakpoints(E, mod):
    """all x in [0,1] where some e*x crosses a k/mod boundary; mod=7 or 2 or 14."""
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(int(e))
        if ae == 0: continue
        for m in range(mod * ae + 1):
            bps.add(F(m, mod * ae))
    return sorted(b for b in bps if 0 <= b <= 1)

def sectors(E, x, mod):
    s = set()
    for e in E:
        v = F(e) * x
        v = v - F(v.numerator // v.denominator)   # frac part
        s.add((v.numerator * mod) // v.denominator)
    return s

def meas_cover(E, mod, need):
    """meas of x where the sector-set (mod 'mod') contains all of 'need' (a set)."""
    bps = breakpoints(E, mod)
    tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        if need <= sectors(E, mid, mod):
            tot += hi - lo
    return tot

def measS7(E):   # need all of {1..6} mod 7
    return meas_cover(E, 7, set(range(1, 7)))

def measParity(E):  # need both 0 and 1 mod 2 (i.e. hit an even AND an odd multiple) -- nonzero parity classes of Z/2 is just {1}; but cover-all-of Z/14\{0} parity-projection = need {0,1} mod 2? careful.
    # Z/14 \ {0} projects onto Z/2 as {0,1} (both parities appear among 1..13). So "all residues covered" requires both parities present.
    return meas_cover(E, 2, {0, 1})

def meas14(E):  # need all of {1,...,13} mod 14
    return meas_cover(E, 14, set(range(1, 14)))

def rand_config(k, hi, seed):
    random.seed(seed)
    while True:
        E = sorted(set([0] + random.sample(range(1, hi), k - 1)))
        if len(E) == k and reduce(gcd, E) == 1:
            return E

def main():
    print("=== ANGLE 1: CRT 2-adic x 7-adic split of the LRC(14) cover ===\n")
    # (P3) parity cover: is it trivially full?
    print("(P3) PARITY cover measParity(E) (meas x hitting both even+odd sectors mod 2):")
    for k in (8, 9, 10):
        cnt_full = 0; vals = []
        for s in range(40):
            E = rand_config(k, 15, 1000 + s)
            mp = measParity(E); vals.append(mp)
            if mp == 1: cnt_full += 1
        print(f"  k={k}: {cnt_full}/40 configs have FULL parity cover; "
              f"min={min(vals)} ({float(min(vals)):.4f})  max={max(vals)} ({float(max(vals)):.4f})")
    print()

    # (P1) product law: meas14 ?= measParity * measS7  (independence of the two CRT factors)
    print("(P1) PRODUCT LAW test: meas14 vs measParity*measS7 (CRT independence):")
    nviol = 0; ntot = 0; worst = F(0)
    for k in (7, 8):
        for s in range(25):
            E = rand_config(k, 13, 2000 + s)
            m14 = meas14(E); m7 = measS7(E); m2 = measParity(E)
            prod = m2 * m7
            ntot += 1
            d = abs(m14 - prod)
            if d > 0: nviol += 1
            if d > worst: worst = d
            if s < 3 and k == 7:
                print(f"  E={E}: meas14={float(m14):.4f} measS7={float(m7):.4f} "
                      f"measPar={float(m2):.4f} prod={float(prod):.4f} diff={float(d):.4f}")
    print(f"  k=7,8: {nviol}/{ntot} configs VIOLATE the product law; worst |diff|={float(worst):.5f}")
    print("  => if product law holds (worst~0), the cover FACTORS by CRT (2-adic indep of 7-adic).")
    print("  => if it fails, the 2 and 7 parts are CORRELATED (resonance lives in the cross term).\n")

    # (P2) per-factor consec-max: is consec the max of measS7 AND of measParity (within bounded stratum)?
    print("(P2) per-factor consec-max (bounded stratum, k=8): consec=(0..7)")
    for k in (8, 9):
        consec = list(range(k))
        cs7 = measS7(consec); cs2 = measParity(consec); cs14 = meas14(consec)
        beat7 = beat2 = beat14 = 0; n = 0
        for s in range(120):
            E = rand_config(k, 14, 3000 + s)
            n += 1
            if measS7(E) > cs7: beat7 += 1
            if measParity(E) > cs2: beat2 += 1
            if meas14(E) > cs14: beat14 += 1
        print(f"  k={k}: consec measS7={float(cs7):.4f} measPar={float(cs2):.4f} meas14={float(cs14):.4f} | "
              f"random beats: S7={beat7}/{n} Par={beat2}/{n} 14={beat14}/{n}")
    print("\nVERDICT depends on output above.")

if __name__ == "__main__":
    main()
