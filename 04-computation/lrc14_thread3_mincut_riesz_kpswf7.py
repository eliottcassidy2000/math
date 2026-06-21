#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD 3 -- ANGLE 6 (min-cut/Foster electrical) + ANGLE 5 (Riesz-product certificate).

Setup recap (verified canon):
  measS7(E) = meas{x in [0,1] : {floor(7*frac(e*x))}_{e in E} covers all of {1..6}}.
  Conductance invariant HYP-2760: measS7 is an EXACT function of the residue-leg structure;
  c_r = sum_{e: e mod 7 == r} 1/|e| (parallel-clock conductance per sector r).
  Foster: sum_r g(r) = 112 conserved (g = R_eff on C_7). consec maximizes total conductance.

ANGLE 6 (min-cut): the cover survives at x iff every sector 1..6 is "hit". Model the 6 nonzero
sectors as a flow network where each speed e contributes capacity to the sectors it can reach
quickly (binding speed = slowest clock per sector). The lonely-region (uncovered) is the
MIN-CUT: the cheapest sector to leave uncovered. measS7 = the survival = product/min over a
bottleneck. PROBE: is measS7 EXACTLY min_r (something), i.e. governed by the single weakest
sector (a min-cut), or is it a genuine product (all sectors)? If min-cut, consec-max <=> consec
maximizes the MIN sector conductance (a clean max-min = LP-dual = certificate).

ANGLE 5 (Riesz product): the decorrelated product sum_t P_t p_t is a finite Riesz product
(arXiv:2511.16636). A Riesz product mu = prod_e (1 + cos(2pi e x)) is a positive measure; its
total mass over the cover region is a DUAL CERTIFICATE for an upper bound on measS7. PROBE:
build the Walsh/Fourier expansion measS7 = (1/64) sum_S qhat_S (canon HYP-2758) and check whether
the consec config maximizes the leading Riesz coefficient and whether a truncated Riesz product
(positive!) gives a valid UPPER bound hugging Q(k-1).
"""
import sys, itertools, random
from fractions import Fraction as F
from functools import reduce
from math import gcd
sys.stdout.reconfigure(line_buffering=True)

def covered_sectors(E, xm):
    s = set()
    for e in E:
        v = F(e) * xm; v = v - F(v.numerator // v.denominator)
        s.add((v.numerator * 7) // v.denominator)
    return s

def measS7(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        if covered_sectors(E, (lo + hi) / 2) >= set(range(1, 7)): tot += hi - lo
    return tot

def meas_miss_sector(E, r):
    """meas of x where sector r (in 1..6) is the/an UNCOVERED sector (sector r not hit)."""
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        if r not in covered_sectors(E, (lo + hi) / 2): tot += hi - lo
    return tot

def rand_config(k, hi, seed):
    random.seed(seed)
    while True:
        E = sorted(set([0] + random.sample(range(1, hi), k - 1)))
        if len(E) == k and reduce(gcd, E) == 1:
            return E

def main():
    print("=== ANGLE 6: is measS7 a MIN-CUT (single weakest sector) or a genuine product? ===\n")
    # measS7 = meas all covered. uncovered = union over r of {sector r missing}.
    # By incl-excl, the uncovered measure U = meas(some sector missing). If sectors were
    # disjoint-failing (min-cut/dominated by ONE sector), U ~ max_r meas_miss(r). Else overlaps.
    print("Compare 1-measS7 (true uncovered) vs sum_r miss_r (union bound) vs max_r miss_r (min-cut):")
    for k in (8, 9, 10):
        consec = list(range(k))
        m = measS7(consec); unc = 1 - m
        miss = [meas_miss_sector(consec, r) for r in range(1, 7)]
        smiss = sum(miss); mmiss = max(miss)
        print(f"  k={k} consec: measS7={float(m):.4f} uncovered={float(unc):.4f} | "
              f"sum_r miss={float(smiss):.4f}  max_r miss={float(mmiss):.4f}  "
              f"| union/true={float(smiss/unc) if unc>0 else 0:.3f} mincut/true={float(mmiss/unc) if unc>0 else 0:.3f}")
    print("  => mincut/true ~1 means ONE sector dominates failure (min-cut model exact);")
    print("     <<1 means many sectors jointly fail (genuine product, no single bottleneck).\n")

    # Does consec MAXIMIZE the min-cut margin = minimize max_r miss_r ? (max-min sector survival)
    print("ANGLE 6b: does consec minimize max_r miss_r (= maximize the weakest sector's survival)?")
    for k in (8, 9):
        consec = list(range(k))
        c_maxmiss = max(meas_miss_sector(consec, r) for r in range(1, 7))
        beat = 0; n = 0
        for s in range(60):
            E = rand_config(k, 14, 7000 + s); n += 1
            mm = max(meas_miss_sector(E, r) for r in range(1, 7))
            if mm < c_maxmiss: beat += 1
        print(f"  k={k}: consec max_r miss={float(c_maxmiss):.4f} | {beat}/{n} random configs have SMALLER max-miss (beat consec)")
    print("  => 0 beats: consec IS the max-min (min-cut LP optimum). >0: min-cut model not the extremal certificate.\n")

    print("=== ANGLE 5: Riesz / Walsh leading coefficient & positive-measure upper bound ===\n")
    # measS7 indicator expanded over the 7-adic Walsh-like basis. Cheap surrogate: the
    # 'energy' sum_r c_r^2 (Riesz product leading term) and whether consec maxes the
    # smallest c_r (the binding Riesz factor). c_r = sum_{e mod7==r} 1/|e|.
    def cvec(E):
        c = [F(0)] * 7
        for e in E:
            e = int(e)
            if e == 0: continue
            c[e % 7] += F(1, abs(e))
        return c
    print("c_r vector (parallel conductance per residue) and min nonzero c_r (Riesz binding factor):")
    for k in (8, 9):
        consec = list(range(k))
        c = cvec(consec)
        nz = [c[r] for r in range(1, 7)]
        print(f"  k={k} consec c_r(r=1..6)={[float(x) for x in nz]} min={float(min(nz)):.4f} sum_sq={float(sum(x*x for x in nz)):.4f}")
        # does consec maximize min_r c_r over bounded bases?
        beat = 0; n = 0
        cmin = min(nz)
        for s in range(80):
            E = rand_config(k, 14, 8000 + s); n += 1
            cc = cvec(E); mn = min(cc[r] for r in range(1, 7))
            if mn > cmin: beat += 1
        print(f"     {beat}/{n} random configs have LARGER min_r c_r (beat consec on min-conductance)")
    print("\nVERDICT depends on output.")

if __name__ == "__main__":
    main()
