#!/usr/bin/env python3
"""S92 (THM-5134 (E), exploratory, not pre-registered): does the abundancy index predict the multiplicity of F1?

r_sigma(N) = number of ordered triples (h1, h2, s) of the sheet sigma with value N.  THM-5134 (E): the local
densities of the sheet at N multiply to S(N) = (8/pi^2) sigma(n)/n, n the odd part of N (sheet s >= 1) or of N + 1
(sheet s = 0).  The squarefree-level product of THM-5128 (D3) is S_sf(N) = (8/pi^2) prod_{p | n} (1 + p/(p^2 - 1)).
For the odd N of a window [2^KLO, 2^KHI) this script computes r_0, r_1 exactly, the two candidate predictors, their
correlations with r, and the mean of r over the N grouped by v_3(n) and by v_5(n) against the limit ratios.
Usage: python3 series_vs_counts.py KHI KLO      (default 25 24)
"""
import math
import sys

import numpy as np
from numba import njit, prange

NT = 8


@njit(cache=True)
def build_q(size):
    q = np.zeros(size, np.int32)
    for h in range(1, size):
        half = h >> 1
        q[h] = q[half] + (((half + 1) >> 1) if (h & 1) else 0)
    return q


@njit(parallel=True, cache=True)
def counts(q, A, B, K, r):
    """r[t, sheet, (N - A) // 2] += multiplicity, N odd in [A, B)."""
    H = q.shape[0]
    for t in prange(NT):
        for s in range(K + 1):
            sheet = 0 if s == 0 else 1
            d0 = 1 if s == 0 else 0
            h1 = np.int64(t)
            while h1 < H:
                y1 = 2 * h1 + 1
                Y = y1 << s
                x1 = y1 - 2 * np.int64(q[h1])
                if (x1 + 2 * Y) * y1 > 2 * (B + d0):
                    if (y1 + 2 * Y) * y1 > 4 * (B + d0):
                        break
                    h1 += NT
                    continue
                y2hi = (2 * (B + d0)) // (x1 + 2 * Y) + 1
                y2lo = (A + d0) // (x1 + Y)
                h2lo = (y2lo - 1) // 2
                if h2lo < h1:
                    h2lo = h1
                h2hi = y2hi // 2 + 1
                if h2hi >= H:
                    h2hi = H - 1
                for h2 in range(h2lo, h2hi + 1):
                    y2 = 2 * h2 + 1
                    x2 = y2 - 2 * np.int64(q[h2])
                    N = x1 * x2 + Y * y2 - d0
                    if N < A or N >= B:
                        continue
                    r[t, sheet, (N - A) >> 1] += 1 if h2 == h1 else 2
                h1 += NT


@njit(cache=True)
def abundancy_odd(limit):
    """a[n] = sigma(n)/n and sf[n] = prod_{p | n} (1 + p/(p^2 - 1)) for odd n < limit (index n)."""
    a = np.zeros(limit, np.float64)
    sf = np.ones(limit, np.float64)
    isp = np.ones(limit, np.uint8)
    for d in range(1, limit, 2):
        inv = 1.0 / d
        for m in range(d, limit, 2 * d):
            a[m] += inv
    for p in range(3, limit, 2):
        if isp[p]:
            f = 1.0 + p / (p * p - 1.0)
            for m in range(p, limit, 2 * p):
                sf[m] *= f
                if m > p:
                    isp[m] = 0
    return a, sf


def main():
    KHI = int(sys.argv[1]) if len(sys.argv) > 1 else 25
    KLO = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    A, B = (1 << KLO) + 1, (1 << KHI) + 1          # odd N in [2^KLO + 1, 2^KHI)
    q = build_q((B >> 1) + 4)
    half = (B - A) >> 1
    r = np.zeros((NT, 2, half), np.int32)
    counts(q, A, B, KHI, r)
    r = r.sum(axis=0).astype(np.float64)
    N = A + 2 * np.arange(half, dtype=np.int64)
    a, sf = abundancy_odd(B + 2)
    n1 = N
    n0 = N + 1
    while True:
        even = (n0 & 1) == 0
        if not even.any():
            break
        n0 = np.where(even, n0 >> 1, n0)
    print(f"odd N in [2^{KLO}, 2^{KHI}): {half} values; mean r_0 = {r[0].mean():.4f}, mean r_1 = {r[1].mean():.4f}; "
          f"mean of sigma(n)/n over the window: sheet 0 {a[n0].mean():.4f}, sheet 1 {a[n1].mean():.4f} (pi^2/8 = {math.pi ** 2 / 8:.4f})")
    ok = True
    for sheet, n in ((0, n0), (1, n1)):
        x_full, x_sf = a[n], sf[n]
        c_full = np.corrcoef(r[sheet], x_full)[0, 1]
        c_sf = np.corrcoef(r[sheet], x_sf)[0, 1]
        print(f"sheet {sheet}: correlation of r with sigma(n)/n: {c_full:.4f}; with the squarefree-level product: {c_sf:.4f}")
        ok &= c_full > 0.1
        # groups by 3-adic and 5-adic valuation of n
        for p in (3, 5):
            v = np.zeros(half, np.int64)
            m = n.copy()
            for _ in range(12):
                div = (m % p) == 0
                v += div
                m = np.where(div, m // p, m)
            base = r[sheet][v == 0].mean()
            cells = []
            for vv in range(0, 5):
                sel = v == vv
                if sel.sum() < 200:
                    continue
                lim = (1 + 1.0 / p) * (1 - p ** (-(vv + 1.0))) / (1 - 1.0 / (p * p))     # L(v)/L(0)
                cells.append(f"v={vv}: {r[sheet][sel].mean() / base:.3f} (limit {lim:.3f}, {int(sel.sum())} N)")
            print(f"    mean of r by v_{p}(n), relative to v = 0: " + "; ".join(cells))
        # quintiles of the abundancy index
        order = np.argsort(x_full, kind="stable")
        qs = np.array_split(order, 5)
        print("    quintiles of sigma(n)/n: mean abundancy -> mean r: "
              + "; ".join(f"{x_full[i].mean():.3f} -> {r[sheet][i].mean():.3f}" for i in qs))
        ok &= r[sheet][qs[-1]].mean() > r[sheet][qs[0]].mean()
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
