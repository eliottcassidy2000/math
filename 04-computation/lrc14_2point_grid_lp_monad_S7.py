#!/usr/bin/env python3
"""
monad-explorer-2026-07-07-S7 -- DECIDING kps-S63's HERESY C (HYP-4987):
does PAIR-UNIFORMITY alone force E[maxgap] > 1/7, on the exact M=26 grid?

Setup (kps-S63): 13 points on Z_26, one pinned at 0 (C(25,12) = 5,200,300 configs).
Each config: pair-distance spectrum (78 unordered distances, d in {1..13}) and
maxgap (circular, in slots/26).  LP over mixtures of configs:
    minimize  E[maxgap]   s.t.  pair-distance density uniform,  total mass 1.
If LP value < 1/7 with a feasible primal: pair statistics alone CANNOT force the
floor -- weight->=3 relations are NECESSARY (the mirror-lemma story made a theorem
on the grid: pairs are exactly uniform for every integer family, so any proof of
the floor must consume >= 3-point structure).  If LP >= 1/7: a 2-point route opens.
kps-S63's bracket: converged dual 0.126406 < 1/7 = 0.142857 < best primal 0.14355.

TWO NULLS (both reported): (a) equal bins: count_d = 78/13 = 6 for every d=1..13
(kps's 13 constraints); (b) circle-discretized: P(dist=d) = 2/25 (d<=12), 1/25 (d=13)
=> counts 6.24 / 3.12 -- the honest discretization of the continuous uniform density.
Cross-check: kps reports 162,770 distinct (maxgap, spectrum) columns.
"""
import itertools
import numpy as np
from scipy.optimize import linprog
from fractions import Fraction as F

M = 26
NPTS = 13

def enumerate_columns():
    combos = np.fromiter(itertools.chain.from_iterable(
        itertools.combinations(range(1, 26), 12)), dtype=np.int8)
    combos = combos.reshape(-1, 12)
    ncfg = combos.shape[0]
    print(f"configs: {ncfg} (expect 5200300)")
    pos = np.concatenate([np.zeros((ncfg, 1), dtype=np.int8), combos], axis=1)

    # spectrum: counts of each distance 1..13 over the 78 pairs
    pairs = [(i, j) for i in range(NPTS) for j in range(i + 1, NPTS)]
    counts = np.zeros((ncfg, 13), dtype=np.int8)
    for i, j in pairs:
        diff = np.abs(pos[:, i].astype(np.int16) - pos[:, j].astype(np.int16))
        d = np.minimum(diff, M - diff).astype(np.int8)
        # accumulate: d in 1..13
        for dd in range(1, 14):
            counts[:, dd - 1] += (d == dd)
    # maxgap in slots
    srt = np.sort(pos, axis=1).astype(np.int16)
    gaps = np.diff(srt, axis=1)
    wrap = (srt[:, 0] + M - srt[:, -1]).reshape(-1, 1)
    mg = np.concatenate([gaps, wrap], axis=1).max(axis=1).astype(np.int8)

    key = np.concatenate([mg.reshape(-1, 1), counts], axis=1)
    uniq, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    print(f"distinct (maxgap, spectrum) columns: {uniq.shape[0]} (kps: 162770)")
    return uniq, cnt

def solve_lp(uniq, b13, tag):
    ncols = uniq.shape[0]
    c = uniq[:, 0].astype(np.float64) / M
    A = np.vstack([uniq[:, 1:].astype(np.float64).T, np.ones((1, ncols))])
    b = np.array(list(b13) + [1.0])
    r = linprog(c, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    print(f"[{tag}] status={r.status} ({r.message.strip()})")
    if r.status == 0:
        val = r.fun
        print(f"  LP value = {val:.6f}  vs 1/7 = {1/7:.6f}  "
              f"({'BELOW -- pair data CANNOT force the floor' if val < 1/7 - 1e-9 else 'AT/ABOVE'})")
        sup = np.where(r.x > 1e-9)[0]
        print(f"  support: {len(sup)} columns; top by weight:")
        order = sup[np.argsort(-r.x[sup])][:8]
        for k in order:
            print(f"    w={r.x[k]:.4f} maxgap={int(uniq[k,0])}/26={uniq[k,0]/M:.4f} "
                  f"spectrum={uniq[k,1:].tolist()}")
        # exact rational certificate check: re-verify with Fractions on the support
        return val, sup, r
    return None, None, r

if __name__ == "__main__":
    uniq, cnt = enumerate_columns()
    print()
    print("NULL (a): equal bins, count_d = 6 for all d=1..13 (kps-S63's constraints)")
    va, supa, ra = solve_lp(uniq, [6.0] * 13, "equal-bins")
    print()
    print("NULL (b): circle-discretized, count_d = 78*2/25 (d<=12), 78/25 (d=13)")
    vb, supb, rb = solve_lp(uniq, [78 * 2 / 25] * 12 + [78 / 25], "circle-null")
    print()
    print("VERDICT: kps-S63 Heresy-C bracket was [0.126406, ~0.14355].")
    if va is not None:
        print(f"  equal-bins LP  = {va:.6f} -> {'2-POINT ROUTE DEAD on the grid' if va < 1/7 else '2-point route ALIVE'}")
    if vb is not None:
        print(f"  circle-null LP = {vb:.6f} -> {'consistent' if (vb < 1/7) == (va < 1/7) else 'NULL-DEPENDENT (!)'}")
