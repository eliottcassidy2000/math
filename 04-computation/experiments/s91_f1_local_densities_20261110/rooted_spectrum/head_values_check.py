#!/usr/bin/env python3
"""S91 (HYP-9661 (C)): head values D - E + U over the rooted census of order <= 9.

Input: rooted_le9.npz (census.py), rows (x = H(B), s, e = E, Yplus, Yminus = U, d, strong, out-degree of the root).
The head value of (B, y) is H(B + z) = D - E + U (z -> y, every other vertex -> z).  Printed: the odd N < 9153 that
are not head values (i) at any root, (ii) at a root that is not a sink (E < D), (iii) with D <= N/3
(U - E >= 2D), (iv) with D <= N/3 at a root of out-degree >= 1.  9153 is the first odd number that the bases of
order <= 9 do not reach as a head.
Checks: (i) misses exactly 7, 21; (ii) exactly 1, 7, 21; (iii) and (iv) exactly 1, 7, 11, 13, 21, 47, 67.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LIM = 9153


def main():
    z = np.load(os.path.join(HERE, "rooted_le9.npz"))
    rows = np.concatenate([z[f"m{m}"] for m in range(1, 10)], axis=0)
    x, e, ym, outdeg = rows[:, 0], rows[:, 2], rows[:, 4], rows[:, 7]
    N = x - e + ym
    odd = np.arange(1, LIM, 2)
    ok = True

    def missing(mask, name, expected):
        v = np.unique(N[mask])
        miss = np.setdiff1d(odd, v[v < LIM]).tolist()
        print(f"{name}: odd N < {LIM} missed: {miss}")
        return miss == expected

    ok &= missing(np.ones(len(N), bool), "(i) all heads", [7, 21])
    ok &= missing(e < x, "(ii) heads at a root that is not a sink", [1, 7, 21])
    ok &= missing((ym - e) >= 2 * x, "(iii) heads with D <= N/3", [1, 7, 11, 13, 21, 47, 67])
    ok &= missing(((ym - e) >= 2 * x) & (outdeg >= 1), "(iv) heads with D <= N/3, root of out-degree >= 1",
                  [1, 7, 11, 13, 21, 47, 67])
    v = np.unique(N)
    first = int(np.setdiff1d(np.arange(23, 1 << 15, 2), v)[0])
    print(f"first odd N > 21 that is not a head value of a base of order <= 9: {first}")
    ok &= first == LIM
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
