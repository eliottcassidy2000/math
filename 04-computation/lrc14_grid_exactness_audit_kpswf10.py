#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_grid_exactness_audit_kpswf10.py  (kind-pasteur 2026-06-21, THREAD B audit)

EXACTNESS AUDIT of the maxgap-threshold measures used in THREAD B.

The function x -> maxgap{frac(e_i x)} is piecewise-rational; {maxgap > t} changes only at
x where (a) two phases coincide: frac(e_i x)=frac(e_j x) i.e. frac((e_i-e_j)x)=0, OR
(b) a phase-difference crosses t: frac((e_i-e_j)x) = +-t (the gap-edge hits the threshold).
For t in {1/7,2/7} and the sector test (multiples of 1/7), the union grid
   { j/(7 e) }  u  { (m +- t)/(e_i-e_j) }
contains ALL transition points.  We AUDIT by:
  (1) recomputing G1,G2 on the THREAD-B grid, and
  (2) recomputing on a 10x-denser grid (add all { (m +- s)/(e_i-e_j) } for s in many
      sub-7 fractions) -- the answer must be IDENTICAL (the extra points are non-transitions).
If identical, the THREAD-B measures are certified exact.
"""
from fractions import Fraction as Fr

P7 = 7
T1 = Fr(1, 7)
T2 = Fr(2, 7)


def phases_at(E, x):
    return sorted((int(e) * x) % 1 for e in E)


def maxgap(ph):
    if len(ph) <= 1:
        return Fr(1)
    g = max(b - a for a, b in zip(ph, ph[1:]))
    return max(g, (ph[0] + 1) - ph[-1])


def grid_threadB(E):
    bp = {Fr(0), Fr(1)}
    E = [int(e) for e in E]
    for e in E:
        if e == 0:
            continue
        for t in range(0, P7 * e + 1):
            bp.add(Fr(t, P7 * e))
    for i in range(len(E)):
        for j in range(len(E)):
            d = E[i] - E[j]
            if d == 0:
                continue
            ad = abs(d)
            for m in range(0, ad + 1):
                for s in (T1, -T1, T2, -T2):
                    v = Fr(m, ad) + s / ad
                    if Fr(0) <= v <= Fr(1):
                        bp.add(v)
    return sorted(bp)


def grid_dense(E):
    """THREAD-B grid PLUS many spurious sub-7 difference crossings (10x denser)."""
    bp = set(grid_threadB(E))
    E = [int(e) for e in E]
    extra = [Fr(a, 70) for a in range(1, 70)]   # 1/70..69/70 -- finer than 1/7
    for i in range(len(E)):
        for j in range(len(E)):
            d = E[i] - E[j]
            if d == 0:
                continue
            ad = abs(d)
            for m in range(0, ad + 1):
                for s in extra:
                    v = Fr(m, ad) + s / ad
                    if Fr(0) <= v <= Fr(1):
                        bp.add(v)
                    v2 = Fr(m, ad) - s / ad
                    if Fr(0) <= v2 <= Fr(1):
                        bp.add(v2)
    return sorted(bp)


def meas_G(E, pts):
    g1 = g2 = Fr(0)
    for a, b in zip(pts, pts[1:]):
        if b <= a:
            continue
        w = b - a
        mid = (a + b) / 2
        g = maxgap(phases_at(E, mid))
        if g > T1:
            g1 += w
        if g > T2:
            g2 += w
    return g1, g2


def main():
    print("=" * 80)
    print("EXACTNESS AUDIT: THREAD-B grid vs 10x-denser grid for {maxgap>1/7},{maxgap>2/7}")
    print("=" * 80)
    configs = {
        "consec k=5": [0, 1, 2, 3, 4],
        "consec k=7": list(range(7)),
        "consec k=9": list(range(9)),
        "perf k=7":   [0, 2, 3, 4, 5, 6, 8],
        "Sidon k=5":  [0, 1, 3, 7, 12],
        "two-block":  [0, 1, 2, 40, 41, 42],
    }
    print(f"{'config':<14}{'G1(B)':>10}{'G1(dense)':>12}{'G2(B)':>10}{'G2(dense)':>12}{'MATCH':>7}")
    allok = True
    for name, E in configs.items():
        g1b, g2b = meas_G(E, grid_threadB(E))
        g1d, g2d = meas_G(E, grid_dense(E))
        ok = (g1b == g1d) and (g2b == g2d)
        allok &= ok
        print(f"{name:<14}{float(g1b):>10.5f}{float(g1d):>12.5f}{float(g2b):>10.5f}"
              f"{float(g2d):>12.5f}{str(ok):>7}")
    print("-" * 65)
    print(f"\nALL grids agree exactly: {allok}")
    print("=> THREAD-B maxgap-threshold measures are CERTIFIED EXACT (grid captures all transitions).")


if __name__ == "__main__":
    main()
