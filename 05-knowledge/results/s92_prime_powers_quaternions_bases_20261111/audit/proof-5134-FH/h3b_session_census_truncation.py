#!/usr/bin/env python3
"""Auditor proof-5134-FH.  Why the triple counts of the session's base-b censuses differ from the complete census.

every_base.py (4), sheet s >= 1, N_b <= 2^24:
   (i)  the table of q is cut at h < 2^21 (q32 = q[: 1 << 21]);
   (ii) a smaller word h1 is skipped ("continue") when its own diagonal value N_b(h1, h1, s) exceeds X, although
        N_b(h1, h2, s) can be <= X for h2 > h1 (psi is not monotone).
bridges_exact.py S10, sheet s = 0, N_b <= 2^22:
   (iii) the table of psi is cut at h < 2^12.
This script counts the ordered triples with N_b <= X
   - completely (valid lower bounds only),
   - with a cap on the words, with and without rule (ii),
and compares with the totals printed in every_base.out and bridges_exact.out.
"""
import numpy as np
from numba import njit


def psi_table(nbits):
    q = np.zeros(1 << nbits, dtype=np.int64)
    for k in range(nbits):
        lo, hi = 1 << k, 1 << (k + 1)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        q[lo:hi] = q[half] + np.where((idx & 1) == 1, (half + 1) >> 1, 0)
    return np.arange(1 << nbits, dtype=np.int64) - q


@njit(cache=True)
def count(psi, b, X, sheet, Hcap, diag_skip):
    """ordered triples with N_b <= X, words below Hcap; diag_skip: skip h1 when N_b(h1, h1, s) > X.
    returns (total, number lost to words >= Hcap is obtained by difference outside)"""
    total = 0
    s0 = 0 if sheet == 0 else 1
    smax = 0 if sheet == 0 else 62
    for s in range(s0, smax + 1):
        if sheet == 1 and (b << (s - 1)) > X:
            break
        for h1 in range(Hcap):
            y1 = 2 * h1 + 1
            lb1 = 1 + b * ((h1 + 1) // 2)
            if sheet == 0:
                seam_d = b * (y1 * y1 - 1) // 2
            else:
                seam_d = (b << (s - 1)) * y1 * y1
            if lb1 * lb1 + seam_d > X:
                break
            X1 = 1 + b * psi[h1]
            if diag_skip and X1 * X1 + seam_d > X:
                continue
            for h2 in range(h1, Hcap):
                y2 = 2 * h2 + 1
                if sheet == 0:
                    seam = b * (y1 * y2 - 1) // 2
                else:
                    seam = (b << (s - 1)) * y1 * y2
                if X1 * (1 + b * ((h2 + 1) // 2)) + seam > X:
                    break
                N = X1 * (1 + b * psi[h2]) + seam
                if N <= X:
                    total += 1 if h2 == h1 else 2
    return total


def main():
    psi = psi_table(24)
    print("A. every_base.py (4): sheet s >= 1, N_b <= 2^24")
    printed = {2: 55892593, 3: 35516685, 4: 24945520, 5: 18836356, 6: 14916908, 10: 7632149}
    X = 1 << 24
    print("   b   complete   cap 2^21   cap 2^21 + diagonal skip   printed in every_base.out   missing   share")
    for b in (2, 3, 4, 5, 6, 10):
        full = count(psi, b, X, 1, 1 << 24, False)
        cap = count(psi, b, X, 1, 1 << 21, False)
        both = count(psi, b, X, 1, 1 << 21, True)
        print(f"  {b:2d}  {full:9d}  {cap:9d}  {both:9d}                 {printed[b]:9d}             {full - printed[b]:8d}  {(full - printed[b]) / full:.4%}"
              f"   reproduces the printed total: {both == printed[b]}")
    print("B. bridges_exact.py S10: sheet s = 0, N_b <= 2^22")
    printed0 = {2: 2845759, 3: 1832488, 5: 985666, 10: 379472}
    X = 1 << 22
    print("   b   complete   cap 2^12   printed in bridges_exact.out   share of the complete census")
    for b in (2, 3, 5, 10):
        full = count(psi, b, X, 0, 1 << 24, False)
        cap = count(psi, b, X, 0, 1 << 12, False)
        print(f"  {b:2d}  {full:9d}  {cap:9d}  {printed0[b]:9d}                      {printed0[b] / full:.4f}"
              f"   reproduces the printed total: {cap == printed0[b]}")
    # an explicit instance of (ii)
    def psi1(h):
        return int(psi[h])
    h1, h2, s, b = 16, 17, 1, 2
    X1, X2 = 1 + b * psi1(h1), 1 + b * psi1(h2)
    y1, y2 = 2 * h1 + 1, 2 * h2 + 1
    print(f"C. instance of (ii) at b = 2, s = 1: N(16,16) = {X1 * X1 + (b << (s - 1)) * y1 * y1}, N(16,17) = {X1 * X2 + (b << (s - 1)) * y1 * y2}:"
          f" for 3201 <= X < 3267 the triple (16, 17, 1) is in T(X) and the rule skips h1 = 16")


if __name__ == "__main__":
    main()
