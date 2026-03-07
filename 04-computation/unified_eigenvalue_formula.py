#!/usr/bin/env python3
"""
THEOREM: Unified eigenvalue formula for M(T_full_n) at ALL n.

For the transitive tournament on n vertices, the transfer matrix M has
tridiagonal structure with M[i,i] = (-1)^i, M[i,i±1] = (-1)^i.

EIGENVALUES: Let K = (n+1)/2. Then:

  mu_j = lambda_j^2 = 3 + 2*cos(j*pi/K)   for j = 1, ..., floor(n/2)

  Eigenvalues: ±sqrt(mu_j) for each j (paired).

  Plus lambda_0 = 1 when n is odd (unpaired eigenvalue).

CASES:
  - Odd n = 2k+1: K = k+1 (integer), j = 1..k, plus lambda_0 = 1
  - Even n = 2k:   K = k+1/2 (half-integer), j = 1..k, full ±pairing

PROOF: The period-2 unit cell transfer matrix T = [[1,1],[1,-1]] gives
Bloch dispersion mu = 3 + 2cos(theta). Boundary conditions:
  - Odd n (free boundaries): theta quantizes as j*pi/(k+1) [integer K]
  - Even n (anti-periodic):  theta quantizes as j*pi/(k+1/2) [half-integer K]

COROLLARIES:
1. rho(M) = sqrt(3 + 2*cos(pi/K)) -> sqrt(5) as n -> infinity.
   Gap: sqrt(5) - rho ~ pi^2/(2n) for large n.
2. Full ±pairing at even n; unpaired lambda=1 at odd n.
3. Spectral bandwidth: mu in [3-2, 3+2] = [1, 5], so lambda in [-sqrt(5), sqrt(5)].
4. Band spectrum in limit: [-sqrt(5), -1] ∪ [1, sqrt(5)] (with optional {1} at odd n).

VERIFICATION: Exact match for n = 2, 3, 4, ..., 25.

opus-2026-03-06-S11b (continued)
"""
import numpy as np

def make_M_transitive(n):
    M = np.zeros((n, n))
    for i in range(n):
        M[i,i] = (-1)**i
        if i < n-1:
            M[i,i+1] = (-1)**i
            M[i+1,i] = (-1)**i
    return M

if __name__ == '__main__':
    print("Unified eigenvalue formula: mu_j = 3 + 2cos(j*pi/K), K=(n+1)/2")
    print("=" * 70)

    all_match = True
    for n in range(2, 26):
        K = (n + 1) / 2
        k = n // 2

        # Predicted eigenvalues
        mu_pred = sorted([3 + 2*np.cos(j*np.pi/K) for j in range(1, k+1)])
        if n % 2 == 1:
            lambda_pred = sorted([1.0] + [np.sqrt(m) for m in mu_pred]
                                 + [-np.sqrt(m) for m in mu_pred])
        else:
            lambda_pred = sorted([np.sqrt(m) for m in mu_pred]
                                 + [-np.sqrt(m) for m in mu_pred])

        # Actual eigenvalues
        M = make_M_transitive(n)
        lambda_actual = sorted(np.linalg.eigvalsh(M))

        match = np.allclose(lambda_pred, lambda_actual, atol=1e-10)
        if not match:
            all_match = False
        rho = max(abs(e) for e in lambda_actual)
        print(f"n={n:2d}: {'MATCH' if match else 'FAIL'}  "
              f"rho={rho:.6f}  gap_to_sqrt5={np.sqrt(5)-rho:.6f}")

    print(f"\n{'ALL MATCH' if all_match else 'SOME FAIL'}")

    # Determinant formula
    print(f"\nDeterminant from eigenvalue product:")
    for n in range(2, 16):
        K = (n + 1) / 2
        k = n // 2
        M = make_M_transitive(n)
        det_actual = int(round(np.linalg.det(M)))

        # From eigenvalues: product of ±sqrt(mu_j) pairs gives -mu_j each
        det_formula = 1
        if n % 2 == 1:
            det_formula = 1  # lambda=1 contributes
        for j in range(1, k+1):
            mu = 3 + 2*np.cos(j*np.pi/K)
            det_formula *= -mu  # each ±pair contributes -mu
        det_formula = int(round(det_formula))

        print(f"  n={n:2d}: det={det_actual:8d}, formula={det_formula:8d}, "
              f"match={det_actual == det_formula}")
