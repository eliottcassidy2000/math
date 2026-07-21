#!/usr/bin/env python3
"""
snf_skew_energy_klein_S399.py  (klein-2026-07-21-S399, ZOO ATLAS gap-fill V.2)

Two of the "6 graph invariants with no tournament analog" (invariant-zoo gap
map II.c), computed on the skew-adjacency matrix S = A - A^T:

  * Smith normal form (integer elementary divisors of S)
  * skew-energy  = sum of |eigenvalues| of S (graph-energy analog)

Finding: the SNF is DEGENERATE at odd n (odd-dimensional skew-symmetric =>
det 0, all-unit elementary divisors) and SUBSUMED BY THE PFAFFIAN at even n
(product of divisors = det = Pf^2; n=4 gives only {(1,1,1,1) Pf=1,
(1,1,3,3) Pf=3}).  So "compute SNF for tournaments" closes as a negative -- it
adds nothing beyond the Pfaffian THM-1475.  Skew-energy adds marginal spectral
resolution but is weak.  See ZOO ATLAS Part V.2.
"""
import numpy as np

exec(open('04-computation/directed_wowii_klein_S397.py').read().split('# ---- collect invariants')[0])


def skew(om, n):
    S = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                S[i][j] = 1 if beats(om, i, j) else -1
    return S


def snf(M):
    """Integer Smith normal form -> tuple of |elementary divisors|."""
    A = [row[:] for row in M]; R = len(A); C = len(A[0]) if R else 0
    divs = []; t = 0

    def swaprows(i, j): A[i], A[j] = A[j], A[i]

    def swapcols(i, j):
        for r in range(R):
            A[r][i], A[r][j] = A[r][j], A[r][i]

    while t < min(R, C):
        piv = None; best = None
        for i in range(t, R):
            for j in range(t, C):
                if A[i][j] != 0 and (best is None or abs(A[i][j]) < best):
                    best = abs(A[i][j]); piv = (i, j)
        if piv is None:
            break
        swaprows(t, piv[0]); swapcols(t, piv[1])
        done = False
        while not done:
            done = True
            for i in range(t + 1, R):
                if A[i][t] != 0:
                    q = A[i][t] // A[t][t]
                    for j in range(t, C):
                        A[i][j] -= q * A[t][j]
                    if A[i][t] != 0:
                        swaprows(t, i); done = False
            for j in range(t + 1, C):
                if A[t][j] != 0:
                    q = A[t][j] // A[t][t]
                    for i in range(t, R):
                        A[i][j] -= q * A[i][t]
                    if A[t][j] != 0:
                        swapcols(t, j); done = False
        ok = True
        for i in range(t + 1, R):
            for j in range(t + 1, C):
                if A[i][j] % A[t][t] != 0:
                    for k in range(t, C):
                        A[t][k] += A[i][k]
                    ok = False; break
            if not ok:
                break
        if ok:
            divs.append(abs(A[t][t])); t += 1
    return tuple(divs)


def energy(S):
    ev = np.linalg.eigvals(np.array(S, dtype=float))
    return round(float(np.sum(np.abs(ev))), 4)


if __name__ == '__main__':
    print("SNF elementary divisors + skew-energy of S=A-A^T over iso classes")
    print("=" * 72)
    for n in (3, 4, 5, 6, 7):
        cls = classes(n)
        fps = {}; snf_only = {}
        for om in cls:
            S = skew(om, n); sf = snf(S); en = energy(S)
            fps[(sf, en)] = fps.get((sf, en), 0) + 1
            snf_only[sf] = snf_only.get(sf, 0) + 1
        print(f" n={n}: {len(cls)} iso classes | distinct SNF {len(snf_only)} | "
              f"distinct (SNF,energy) {len(fps)} | biggest merged cell {max(fps.values())}")
        if n <= 5:
            for sf in sorted(snf_only):
                print(f"     SNF {sf}: {snf_only[sf]} classes")
    print("\n READING: SNF degenerate at odd n; subsumed by Pfaffian at even n. Closed negative.")
