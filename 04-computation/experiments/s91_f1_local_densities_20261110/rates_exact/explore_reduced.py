#!/usr/bin/env python3
"""exploration: the complement symmetry C and the reduced m x m operators T_k^+, T_k^-; timing of eigvals."""
import sys
import time

import numpy as np

from explore_invariant import twisted_matrix


def reduced(m, k, sign):
    """(T f)(eta) = f(2 eta) + sign * w^(-(eta+1)) f(-2 eta - 2),  w = zeta^(k/2)."""
    inv2 = pow(2, -1, m)
    kk = (k * inv2) % m
    eta = np.arange(m)
    M = np.zeros((m, m), dtype=np.complex128)
    np.add.at(M, (eta, (2 * eta) % m), 1.0)
    np.add.at(M, (eta, (-2 * eta - 2) % m), sign * np.exp(-2j * np.pi * ((kk * (eta + 1)) % m) / m))
    return M


def main():
    worst = 0.0
    for m in (3, 5, 7, 9, 11, 13, 15, 21, 31, 45):
        for k in range(1, m):
            T = twisted_matrix(m, k)
            # C u (eta, beta) = w^(-(eta - beta)) u(-1 - eta, 1 - beta)
            inv2 = pow(2, -1, m)
            C = np.zeros((2 * m, 2 * m), dtype=np.complex128)
            for eta in range(m):
                for beta in range(2):
                    e = (k * inv2 * (eta - beta)) % m
                    C[2 * eta + beta, 2 * ((-1 - eta) % m) + (1 - beta)] = np.exp(-2j * np.pi * e / m)
            comm = np.abs(T @ C - C @ T).max()
            sq = np.abs(C @ C - np.exp(2j * np.pi * k / m) * np.eye(2 * m)).max()
            ev = np.sort_complex(np.linalg.eigvals(T))
            ev2 = np.sort_complex(np.concatenate([np.linalg.eigvals(reduced(m, k, +1)), np.linalg.eigvals(reduced(m, k, -1))]))
            # match multisets greedily
            d = 0.0
            rest = list(ev2)
            for e in ev:
                j = int(np.argmin([abs(e - r) for r in rest]))
                d = max(d, abs(e - rest[j]))
                rest.pop(j)
            worst = max(worst, comm, sq, d)
    print("max |TC - CT|, |C^2 - zeta^k|, spectrum mismatch over the tested (m, k):", worst)
    for m in (301, 601):
        t0 = time.time()
        for k in (1, 2, 3):
            np.linalg.eigvals(reduced(m, k, +1))
        print(f"eigvals of the reduced {m} x {m}: {(time.time() - t0) / 3:.3f} s each")
    ok = worst < 1e-9
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
