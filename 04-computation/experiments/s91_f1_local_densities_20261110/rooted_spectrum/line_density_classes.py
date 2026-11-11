#!/usr/bin/env python3
"""S91 (HYP-9661 (A'), exploratory, not pre-registered): the line counts r(N) by residue class of N.

r_10(N) = number of one-vertex sources (D, U) with D + 2U = N from rooted tournaments of order <= 10 (with appended
sinks), as in line_density_order10.py.  The source set is closed under (D, U) -> (m D, m U) for every attained m
(put a tournament with m Hamiltonian paths before the base), so a line with m | N inherits the sources of the line
N/m.  Printed: for p = 3, 5, 7, 9, 11, 13, 15 the mean of r(N) / (N/6) over odd N in [2^11, 2^13) by class of N
modulo p (class 0 first), the ratio class 0 / other classes, and the share of sources on a line with p | D
(then p | U as well when p | N).  Check: the class 0 is the favoured class at p = 3, 5, 9, 15 (attained moduli).
"""
import os
import sys

import numpy as np

import line_density_order10 as ld

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    NMAX = ld.NMAX
    M = np.zeros((NMAX >> 1, (NMAX // 3 >> 1) + 2), dtype=np.uint8)
    U = 1
    while 1 + 2 * U < NMAX:
        M[(1 + 2 * U) >> 1, 0] = 1
        U *= 2
    z = np.load(os.path.join(HERE, "reps_le8.npz"))
    for n in range(1, 9):
        ld.scan(n, z["n%d" % n].astype(np.int64), M, NMAX)
    ld.scan(9, np.load(os.path.join(HERE, "codes9.npy")).astype(np.int64), M, NMAX)
    r = M.sum(axis=1).astype(np.float64)
    lo, hi = 1 << 11, 1 << 13
    odd = np.arange(lo + 1, hi, 2)
    share = r[odd >> 1] / (odd / 6.0)
    ok = True
    Dodd = 2 * np.arange(M.shape[1]) + 1
    print(f"odd N in [{lo}, {hi}): mean share of lattice points that are sources = {share.mean():.4f}")
    print(" p | share at N = 0 (mod p) | share at the other classes (min, max of class means) | ratio | "
          "share of the sources of a line p | N that have p | D")
    for p in (3, 5, 7, 9, 11, 13, 15):
        cm = np.array([share[odd % p == c].mean() for c in range(p)])
        rest = share[odd % p != 0].mean()
        sel = odd[odd % p == 0]
        num = 0.0
        den = 0.0
        for N in sel:
            row = M[N >> 1]
            den += row.sum()
            num += row[(Dodd % p) == 0].sum()
        print(f"{p:2d} | {cm[0]:.4f} | {rest:.4f} ({cm[1:].min():.4f}, {cm[1:].max():.4f}) | {cm[0] / rest:.3f} | "
              f"{num / den:.3f}  (1/p = {1 / p:.3f})")
        if p in (3, 5, 9, 15):
            ok &= cm[0] > cm[1:].max()
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
