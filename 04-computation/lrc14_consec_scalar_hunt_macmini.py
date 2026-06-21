#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_consec_scalar_hunt_macmini.py  (mac-mini 2026-06-21, THREAD D NEW LEAD)

GOAL: find a CLEANER single governing scalar than additive energy (HYP-2735) that
ORDERS measS7 monotonically among primitive k-sets, giving a non-finite-check proof
of consec-max.

The honest gap: consec [1..8] and AP-d2 [2,4,...,16] have EQUAL additive energy 344
yet measS7 = 0.345 vs 0.309. So AE is necessary-not-sufficient. What distinguishes them?

CANDIDATE SCALARS (each a single number computed from E):
  AE   = sum_s r(s)^2                       (additive energy; the current lead)
  AE3  = sum_s r3(s)^2 with r3 = 3-fold rep counts  (higher additive energy)
  MED  = multiplicative energy sum_s m(s)^2, m=#{(a,b):ab=s}
  diam = max-min (span)
  invE = sum_{a in E} 1/a                    (the "inverse mass" — sectors are frac(e v))
  sumlog = sum log(e)                        (geometric span proxy)
  sqrtE = sum sqrt(e)
  recipE = #distinct ratios e_i/e_j
  COVAR = the COVERAGE second moment: E_v[ (#sectors hit)^2 ] -- directly tied to p0
  CSPEC = E_v[ #sectors hit ] (first moment) and var.
  D1    = #distinct DIFFERENCES |a-b|
  D1mult= sum_d (mult of difference d)^2  (difference energy)

KEY IDEA being tested: measS7 = P(all 7 sectors hit). The "sector(e v)" for fixed v are
NOT independent; the dependence is governed by how the e's RESONATE = small linear
relations sum n_i e_i = 0 with small |n|. consec has the DENSEST such relations
(e_{i+1}-e_i=1, all gaps equal). The right scalar might be a RELATION-WEIGHTED count.

This script computes ALL candidate scalars + measS7 (exact) over primitive k-sets and
reports the SPEARMAN-style monotonicity (fraction of concordant pairs) of each scalar
with measS7, and specifically whether any scalar SEPARATES consec from AP-d2 in the
right direction.
"""
import itertools
from fractions import Fraction as Fr
from math import gcd, log, sqrt
from collections import Counter

P = 7
def sector(yf): return int(P * yf)

def breakpoints(E):
    bp = {Fr(0), Fr(1)}
    for e in E:
        if e == 0: continue
        for t in range(0, P * e):
            bp.add(Fr(t, P * e))
    return sorted(bp)

def measS7(E):
    E = [int(e) for e in E if int(e) != 0]
    xs = breakpoints(E); tot = Fr(0)
    for a, b in zip(xs, xs[1:]):
        mid = (a + b) / 2
        if len(set(sector((e * mid) % 1) for e in E)) == P:
            tot += (b - a)
    return tot

def coverage_moments(E):
    """E_v[#sectors hit] and E_v[(#sectors hit)^2] over v, exact."""
    E = [int(e) for e in E if int(e) != 0]
    xs = breakpoints(E); m1 = Fr(0); m2 = Fr(0)
    for a, b in zip(xs, xs[1:]):
        mid = (a + b) / 2
        h = len(set(sector((e * mid) % 1) for e in E))
        L = b - a
        m1 += h * L; m2 += h * h * L
    return m1, m2

def add_energy(E):
    r = Counter()
    for a in E:
        for b in E: r[a + b] += 1
    return sum(v * v for v in r.values())

def diff_energy(E):
    r = Counter()
    for a in E:
        for b in E: r[a - b] += 1
    return sum(v * v for v in r.values())

def mult_energy(E):
    r = Counter()
    for a in E:
        for b in E: r[a * b] += 1
    return sum(v * v for v in r.values())

def relation_weight(E, Nmax=3):
    """sum over small integer relations sum n_i e_i = 0, |n_i|<=Nmax, n != 0, of 1.
    Counts the small linear dependencies — the true source of sector correlation."""
    E = list(E); k = len(E)
    cnt = 0
    # too big for k=8 full; restrict to relations supported on <= 3 elements (the dominant ones)
    for sup in itertools.combinations(range(k), 3):
        for coeffs in itertools.product(range(-Nmax, Nmax+1), repeat=3):
            if all(c == 0 for c in coeffs): continue
            if sum(coeffs[t]*E[sup[t]] for t in range(3)) == 0:
                cnt += 1
    return cnt

def scalars(E):
    E = list(E)
    m1, m2 = coverage_moments(E)
    return {
        "AE": add_energy(E),
        "DiffE": diff_energy(E),
        "MultE": mult_energy(E),
        "span": max(E) - min(E),
        "invE": float(sum(Fr(1, e) for e in E)),
        "sumlog": sum(log(e) for e in E),
        "Ndiff": len(set(a - b for a in E for b in E if a != b)),
        "cov_m1": float(m1),
        "cov_m2": float(m2),     # second moment of coverage -- closest to p0
        "rel3": relation_weight(E, 3),
    }

def concordance(rows, key, target="p0", direction=+1):
    """fraction of pairs (A,B) with target strictly ordered that the scalar orders the same way."""
    n = len(rows); conc = 0; disc = 0
    for i in range(n):
        for j in range(i+1, n):
            ti = rows[i][target]; tj = rows[j][target]
            if ti == tj: continue
            si = rows[i][key]; sj = rows[j][key]
            if si == sj: continue
            # does scalar (in given direction) agree with target order?
            agree = (si - sj) * (ti - tj) * direction > 0
            if agree: conc += 1
            else: disc += 1
    return conc / (conc + disc) if (conc+disc) else float('nan')

def main():
    print("#"*80)
    print("# CLEAN-SCALAR HUNT for consec-max (THREAD D)")
    print("#"*80)

    # The two test families: consec vs AP-d2 (the honest-gap pair)
    consec = list(range(1, 9))
    apd2 = [2*i for i in range(1, 9)]    # [2,4,...,16] -- but primitive? gcd=2. dilation-equiv to [1..8]??
    print("\n=== the honest-gap pair ===")
    for name, E in [("consec[1..8]", consec), ("AP-d2[2..16]", apd2),
                    ("AP-d2 prim[1,3,5,7,9,11,13,15]", [1,3,5,7,9,11,13,15])]:
        v = measS7(E)
        print(f"  {name:32s} p0={float(v):.5f} AE={add_energy(E)} "
              f"DiffE={diff_energy(E)} span={max(E)-min(E)}")
    # NOTE: [2,4,...,16] = 2*[1..8] so by dilation invariance measS7 MUST equal consec.
    # The real AP-d2 comparison in the ledger must be the PRIMITIVE odd AP [1,3,..,15].

    print("\n  (dilation check) measS7([2,4..16]) == measS7([1..8])? ",
          measS7(apd2) == measS7(consec))

    # build the population: primitive 8-sets, min=1, max<=15 (manageable, ~6435)
    print("\n=== building population (primitive 8-sets, min=1, max<=15) ===")
    rows = []
    for combo in itertools.combinations(range(2, 16), 7):
        E = (1,) + combo
        g = 0
        for e in E: g = gcd(g, e)
        if g != 1: continue
        sc = scalars(E)
        sc["p0"] = float(measS7(E))
        sc["E"] = E
        rows.append(sc)
    print(f"  population size = {len(rows)}")

    # concordance of each scalar with p0
    print("\n=== concordance (Kendall-style) of each scalar with measS7 (1.0=perfect order) ===")
    keys = ["AE","DiffE","MultE","span","invE","sumlog","Ndiff","cov_m1","cov_m2","rel3"]
    dirs = {"AE":+1,"DiffE":+1,"MultE":+1,"span":-1,"invE":+1,"sumlog":-1,"Ndiff":-1,
            "cov_m1":+1,"cov_m2":+1,"rel3":+1}
    res = []
    for key in keys:
        c = concordance(rows, key, "p0", dirs[key])
        res.append((c, key))
    for c, key in sorted(res, reverse=True):
        print(f"  {key:8s} dir={dirs[key]:+d}  concordance={c:.4f}")

    # Does cov_m2 (coverage 2nd moment) order p0 better than AE? Inclusion-exclusion:
    # p0 = P(cov=7).  cov_m2 = E[cov^2].  These should be tightly related.
    print("\n=== top of population by measS7 (is consec #1? what scalar predicts it?) ===")
    top = sorted(rows, key=lambda r: -r["p0"])[:6]
    for r in top:
        print(f"  p0={r['p0']:.5f} AE={r['AE']} DiffE={r['DiffE']} cov_m2={r['cov_m2']:.4f} "
              f"rel3={r['rel3']} span={r['span']} E={r['E']}")

    print("\nDONE.")

if __name__ == "__main__":
    main()
