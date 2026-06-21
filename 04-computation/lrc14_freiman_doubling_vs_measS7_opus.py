#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE 3 -- FREIMAN INVERSE / additive-rigidity (opus, 2026-06-21).

CLAIM (Freiman-rigidity hypothesis): the measS7-maximizer is the MINIMAL-DOUBLING
(most additively-structured) set. By an inverse sumset theorem (Freiman/Green-Ruzsa)
small doubling => AP up to Freiman isomorphism, and among the relevant family the
AP/consec set is the doubling-minimizer. So:
    argmax measS7  ==  argmin |E+E|/|E|  ==  AP (consec, gcd-normalized).

We TEST this falsifiable correlation over the full-residue stratum (k=8, span <= W):
  for each shape E:
    measS7(E)                      (canonical occupancy: P(N=0) = pi[7])
    doubling D(E)   = |E+E| / |E|
    sumset size     = |E+E|
    diffset size    = |E-E|
    additive energy Ea(E) = #{(a,b,c,d) in E^4 : a+b=c+d}
  Then:
   (T1) Is argmax measS7 == argmin doubling? (and == argmax additive energy?)
   (T2) Spearman/rank correlation of measS7 vs doubling (should be NEGATIVE,
        ideally monotone), and vs additive energy (POSITIVE).
   (T3) Top-10 measS7 shapes vs top-10 minimal-doubling shapes: overlap?
   (T4) Is the maximizer FORCED to be an AP? (does min-doubling among full-residue
        shapes uniquely select consec / an AP?)

HONESTY: a correlation alone is NOT a proof. The decisive test is whether the
RANKING agrees exactly (argmax==argmin) AND whether the maximizer is provably AP.
We report exactly where the correlation BREAKS (it likely does: measS7 is a torus
covering measure, not literally a sumset functional). The break tells us whether
the Freiman route is alive or a clean refutation.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

# ----- canonical measS7 (occupancy, P(N=0) = pi[7]) -----
def occupancy_pi(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7*ae+1): bps.add(F(m, 7*ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    pi = [F(0)]*8
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        xm = (lo+hi)/2
        hit = set()
        for e in E:
            v = e*xm; v = v - (v.numerator // v.denominator)
            hit.add((v.numerator*7)//v.denominator)
        pi[len(hit)] += hi - lo
    return pi

def measS7(E):
    return occupancy_pi(E)[7]

# ----- additive functionals -----
def sumset(E):
    return {a+b for a in E for b in E}
def diffset(E):
    return {a-b for a in E for b in E}
def doubling(E):
    return F(len(sumset(E)), len(E))
def additive_energy(E):
    cnt = defaultdict(int)
    for a in E:
        for b in E:
            cnt[a+b] += 1
    return sum(v*v for v in cnt.values())   # #{(a,b,c,d): a+b=c+d}

# ----- stratum -----
def residues(E): return frozenset(e % 7 for e in E)
def is_full_residue(E): return residues(E) == frozenset(range(7))
def primitive(E): return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

def is_AP(E):
    """Is E (as a set of integers) an arithmetic progression?"""
    s = sorted(set(E))
    if len(s) < 2: return True
    d = s[1]-s[0]
    return all(s[i+1]-s[i] == d for i in range(len(s)-1))

def spearman(xs, ys):
    """Spearman rank correlation (exact-ish via float on ranks). Ties averaged."""
    def ranks(vals):
        idx = sorted(range(len(vals)), key=lambda i: vals[i])
        r = [0.0]*len(vals)
        i = 0
        while i < len(idx):
            j = i
            while j+1 < len(idx) and vals[idx[j+1]] == vals[idx[i]]:
                j += 1
            avg = (i + j)/2.0 + 1
            for k in range(i, j+1): r[idx[k]] = avg
            i = j+1
        return r
    rx, ry = ranks(xs), ranks(ys)
    n = len(xs)
    mx = sum(rx)/n; my = sum(ry)/n
    cov = sum((rx[i]-mx)*(ry[i]-my) for i in range(n))
    sx = sum((rx[i]-mx)**2 for i in range(n))**0.5
    sy = sum((ry[i]-my)**2 for i in range(n))**0.5
    if sx == 0 or sy == 0: return float('nan')
    return cov/(sx*sy)

if __name__ == "__main__":
    print("#"*78)
    print("# FREIMAN DOUBLING vs measS7  (ANGLE 3 -- additive rigidity)")
    print("#"*78)
    for k, W in [(8, 11), (8, 12), (8, 13)]:
        C = consec(k)
        bank = [(0,)+r for r in itertools.combinations(range(1, W+1), k-1)]
        full = [list(E) for E in bank if primitive(E) and is_full_residue(E)]
        print(f"\n{'='*70}\nk={k}, span<= {W}: {len(full)} full-residue primitive shapes")
        rows = []
        for E in full:
            rows.append((E, measS7(E), doubling(E), len(sumset(E)),
                         len(diffset(E)), additive_energy(E), is_AP(E)))
        # T1: argmax measS7, argmin doubling, argmax energy
        best_meas = max(rows, key=lambda r: r[1])
        min_doub  = min(rows, key=lambda r: r[2])
        max_en    = max(rows, key=lambda r: r[5])
        print(f"  argmax measS7 : {best_meas[0]}  measS7={float(best_meas[1]):.6f}  "
              f"doub={best_meas[2]}={float(best_meas[2]):.3f}  energy={best_meas[5]}  AP={best_meas[6]}")
        print(f"  argmin doublng: {min_doub[0]}  measS7={float(min_doub[1]):.6f}  "
              f"doub={min_doub[2]}={float(min_doub[2]):.3f}  energy={min_doub[5]}  AP={min_doub[6]}")
        print(f"  argmax energy : {max_en[0]}  measS7={float(max_en[1]):.6f}  "
              f"doub={max_en[2]}={float(max_en[2]):.3f}  energy={max_en[5]}  AP={max_en[6]}")
        consec_is_max = (tuple(best_meas[0]) == tuple(C))
        consec_is_mindoub = (tuple(min_doub[0]) == tuple(C))
        print(f"  consec={C}: is measS7-max? {consec_is_max}  is doubling-min? {consec_is_mindoub}")
        print(f"  >>> argmax measS7 == argmin doubling ? {tuple(best_meas[0])==tuple(min_doub[0])}")
        print(f"  >>> argmax measS7 == argmax energy   ? {tuple(best_meas[0])==tuple(max_en[0])}")
        # T2: rank correlations
        ms = [float(r[1]) for r in rows]
        db = [float(r[2]) for r in rows]
        en = [float(r[5]) for r in rows]
        ss = [float(r[3]) for r in rows]
        print(f"  Spearman(measS7, doubling)  = {spearman(ms, db):+.4f}  (want NEGATIVE)")
        print(f"  Spearman(measS7, energy)    = {spearman(ms, en):+.4f}  (want POSITIVE)")
        print(f"  Spearman(measS7, |E+E|)     = {spearman(ms, ss):+.4f}  (want NEGATIVE)")
        # T3: top-10 overlap
        top_meas = sorted(rows, key=lambda r: -r[1])[:10]
        top_doub = sorted(rows, key=lambda r: r[2])[:10]
        set_meas = {tuple(r[0]) for r in top_meas}
        set_doub = {tuple(r[0]) for r in top_doub}
        print(f"  top-10 measS7 ∩ top-10 min-doubling = {len(set_meas & set_doub)}/10")
        # T4: is the maximizer forced to be AP? How many distinct doubling values
        # hit the minimum, and are they all APs?
        mind = min_doub[2]
        at_min = [r for r in rows if r[2] == mind]
        print(f"  shapes achieving MIN doubling ({mind}): {len(at_min)}; all AP? "
              f"{all(r[6] for r in at_min)}")
        for r in at_min[:8]:
            print(f"      {r[0]}  measS7={float(r[1]):.6f}  AP={r[6]}")
        # Among the top-5 measS7, list doubling+AP
        print("  TOP-5 measS7 shapes:")
        for r in top_meas[:5]:
            print(f"      {r[0]}  measS7={float(r[1]):.6f}  doub={float(r[2]):.3f}  "
                  f"|E+E|={r[3]} |E-E|={r[4]} energy={r[5]} AP={r[6]}")
        print("  TOP-5 min-doubling shapes:")
        for r in sorted(rows, key=lambda r:(r[2], -r[1]))[:5]:
            print(f"      {r[0]}  measS7={float(r[1]):.6f}  doub={float(r[2]):.3f}  "
                  f"|E+E|={r[3]} |E-E|={r[4]} energy={r[5]} AP={r[6]}")
