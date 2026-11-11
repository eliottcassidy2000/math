#!/usr/bin/env python3
"""exploration: structure of the invariant subspaces of T_k (numerical, high precision)."""
import sys

import mpmath as mp
import numpy as np

mp.mp.dps = 40


def twisted_matrix(m, k):
    inv2 = pow(2, -1, m)
    zeta = np.exp(2j * np.pi / m)
    T = np.zeros((2 * m, 2 * m), dtype=np.complex128)
    for eta in range(m):
        for beta in range(2):
            i = 2 * eta + beta
            T[i, 2 * ((2 * eta) % m) + 0] += 1.0
            e = (k * ((eta + beta) * inv2 % m)) % m
            T[i, 2 * ((2 * eta + 1) % m) + 1] += zeta ** e
    return T


def main():
    m = int(sys.argv[1])
    for k in range(1, m):
        T = twisted_matrix(m, k)
        ev, R = np.linalg.eig(T)
        order = np.argsort(-np.abs(ev))
        ev = ev[order]
        print(f"m = {m} k = {k}: |ev| =", np.round(np.abs(ev), 6))
        print("   arg/2pi =", np.round(np.angle(ev) / (2 * np.pi), 6))
        print("   sum ev =", np.round(ev.sum(), 8), " prod =", np.round(np.prod(ev), 8), " |prod| =", round(abs(np.prod(ev)), 8))
        # power sums
        for j in range(1, 8):
            print(f"   tr T^{j} = {np.round(np.trace(np.linalg.matrix_power(T, j)), 8)}")


if __name__ == "__main__":
    main()
