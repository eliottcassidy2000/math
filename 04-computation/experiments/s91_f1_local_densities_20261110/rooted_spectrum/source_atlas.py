#!/usr/bin/env python3
"""S91: atlas of singleton-marked offset sources through order 9 (the "supply" question of Session261).

A singleton-marked source is a rooted tournament (B, y) with offset D = H(B) and slope U = Yminus(B, y), the number
of ordered covers (P1, P2) of B by two vertex-disjoint paths with P2 starting at y (P1 may be empty); in the census
rows (x, s, e, Yplus, Yminus, d, strong, outdeg) these are columns 0 and 4 (by converse also 0 and 3).  The inherited
tail law (THM-5126 (1'), Session261 constructors/NINE_OFFSET_EXTENSION.md) gives tournaments of count D + 2^t U for
every t >= 1.  Check here: the sources of Session261, (5, 9), (9, 9), (17, 18), (33, 36), (65, 288), are in the census.
The grid point
        N(delta, r, s) = 1 + 2^s + delta 2^(r+s)        (delta odd, r, s >= 1)
is supplied for every r >= 1 by a source with D = 2^s + 1 and U = delta 2^j, j <= s (tail length t = r + s - j).
The pendant increment H(B + v) - H(B) = Yminus(B, y) - e(B, y) (v beats only y) satisfies the domination-defect
bound >= 2^outdeg(y) (Session261 EXTENSION4 with the star interface K_(1, outdeg)).

Input: rooted_le9.npz (the complete rooted census of orders 1..9 written by census.py; not kept in the repository:
regenerate with `python3 census.py`, about 30 s).
Output (FINITE-EXACT, all rooted tournaments of order <= 9):
  * for each level s = 1..9 (D = 3, 5, ..., 513): number of slopes U, the odd parts delta supplied for all r
    (j <= s), the least odd delta not supplied, and the least delta not supplied even with a finite head (any j);
  * the row of the 9-grid and of the 1-, 3-, 5-, 7-grids: least j at each level;
  * the lower bound U >= 2^outdeg(y) (domination-defect inequality for the pendant vertex) and its equality cases.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    z = np.load(os.path.join(HERE, "rooted_le9.npz"))
    rows = np.concatenate([z[f"m{m}"] for m in range(1, 10)], axis=0)
    x, s, e, yp, ym, d, strong, outdeg = (rows[:, i] for i in range(8))
    ok = True
    # converse symmetry of the row set: (x, Yminus - e) and (x, Yplus - s) have the same set of pairs
    P1 = set(zip(x.tolist(), ym.tolist()))
    P2 = set(zip(x.tolist(), yp.tolist()))
    ok &= P1 == P2
    print("pairs (D, U) from Yminus and from Yplus agree as sets (converse):", P1 == P2, " #pairs:", len(P1))
    # domination-defect lower bound for the pendant vertex: U >= 2^k, k = out-degree of the root (k >= 1); U = 0 iff k = 0
    known = [(5, 9), (9, 9), (17, 18), (33, 36), (65, 288)]
    ok &= all(pr in P1 for pr in known)
    print("the five sources of Session261 are in the census:", all(pr in P1 for pr in known))
    U = ym - e
    k = outdeg
    viol = int(((k >= 1) & (U < (1 << np.minimum(k, 60)))).sum()) + int(((k == 0) & (U != 0)).sum())
    eq = (k >= 1) & (U == (1 << np.minimum(k, 60)))
    print(f"Yminus - e >= 2^outdeg(y) for every rooted tournament of order <= 9 (0 at a sink): violations {viol}; "
          f"equality in {int(eq.sum())} of {len(U)} rows; values of D at equality: {sorted(set(x[eq].tolist()))[:12]} ...")
    ok &= viol == 0
    byD = {}
    for D, u in P1:
        byD.setdefault(D, set()).add(u)
    print("\nlevel s, D = 2^s + 1: #slopes U > 0 | odd parts delta < 200 supplied for all r (j <= s): count, least missing | "
          "least delta missing with any j")
    table = {}
    for lev in range(1, 10):
        D = (1 << lev) + 1
        us = sorted(u for u in byD.get(D, ()) if u > 0)
        jmin = {}
        for u in us:
            j = (u & -u).bit_length() - 1
            dl = u >> j
            if dl not in jmin or j < jmin[dl]:
                jmin[dl] = j
        table[lev] = jmin
        full = sorted(dl for dl, j in jmin.items() if j <= lev and dl < 200)
        miss_full = next(dl for dl in range(1, 10 ** 6, 2) if not (dl in jmin and jmin[dl] <= lev))
        miss_any = next(dl for dl in range(1, 10 ** 6, 2) if dl not in jmin)
        print(f"  s = {lev}, D = {D:4d}: {len(us):6d} slopes | {len(full):3d} of 100, least missing {miss_full:4d} | {miss_any:4d}")
    print("\nleast j with (2^s + 1, delta 2^j) a source of order <= 9 ('.' = none; '*' = j > s, finite head only):")
    print("  delta | " + " ".join(f"s={lev}" for lev in range(1, 10)))
    for dl in list(range(1, 42, 2)) + [43, 45, 63, 65, 67, 85, 127, 129, 163]:
        cells = []
        for lev in range(1, 10):
            j = table[lev].get(dl)
            cells.append("  ." if j is None else (f"{j:3d}" if j <= lev else f"{j:2d}*"))
        print(f"  {dl:5d} | " + " ".join(cells))
    # the 9-grid of Session261: levels 2..6 exist with j = 0, 0, 1, 2, 5 (their sources); what about 7, 8, 9?
    nine = {lev: table[lev].get(9) for lev in range(1, 10)}
    print("\n9-grid: least j by level:", nine)
    ok &= all(nine[lev] is not None and nine[lev] <= lev for lev in (2, 3, 4, 5, 6))
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
