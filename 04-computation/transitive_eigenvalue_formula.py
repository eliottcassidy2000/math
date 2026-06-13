#!/usr/bin/env python3
"""
THEOREM: Exact eigenvalue formula for M(T_full_{2k+1}).

For the transitive tournament on 2k+1 vertices, the transfer matrix
M(T_full_n) has eigenvalues:

  lambda_j = ±sqrt(3 + 2*cos(j*pi/(k+1)))   for j = 1, ..., k
  lambda_0 = 1

PROOF: M is tridiagonal with M[i,i] = (-1)^i, M[i,i+1] = M[i+1,i] = (-1)^i.
This has period-2 structure. The period-2 "unit cell" transfer matrix is
T = [[1,1],[1,-1]] with eigenvalues lambda satisfying the Bloch condition:
  |tr(T)/(2*det(T)^{1/2})| = |lambda^2 - 3| / 2 <= 1
  => lambda^2 in [1, 5]

For finite n = 2k+1 (odd), the boundary conditions quantize the Bloch angle:
  theta_j = j*pi/(k+1),  j = 1, ..., k
  lambda^2 = 3 + 2*cos(theta_j)

This gives the band spectrum [-sqrt(5), -1] ∪ {1} ∪ [1, sqrt(5)] in the
limit k -> infinity.

COROLLARIES:
1. ρ(M) = sqrt(3 + 2*cos(pi/(k+1))) -> sqrt(5) as k -> infinity.
2. Spectral gap: min |lambda| = 1 (the unpaired eigenvalue).
3. tr(M^{2k+1}) = 1 for all k >= 0 (eigenvalue pairing cancellation).
4. det(M) = product of eigenvalues = 1 * prod (3+2cos(j*pi/(k+1)))
           = prod_{j=1}^{k} (3+2cos(j*pi/(k+1)))  [sign: (-1)^{k*(k-1)/2}]

VERIFICATION: Exact for n = 3, 5, 7, 9, 11, 13, 15, 17, 19.

NOTE: For EVEN n = 2k, the boundary conditions change, giving a different
quantization: the eigenvalues are NOT simply 3 + 2*cos(j*pi/(k+1)).
The even-n case requires separate analysis.

opus-2026-03-06-S11b
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
    print("Exact eigenvalue formula for M(T_full_{2k+1})")
    print("=" * 60)

    for n in range(3, 20, 2):
        k = (n - 1) // 2

        # Predicted eigenvalues
        mu_pred = [3 + 2*np.cos(j*np.pi/(k+1)) for j in range(1, k+1)]
        lambda_pred = sorted([1.0] + [np.sqrt(mu) for mu in mu_pred] + [-np.sqrt(mu) for mu in mu_pred])

        # Actual eigenvalues
        M = make_M_transitive(n)
        lambda_actual = sorted(np.linalg.eigvalsh(M))

        match = np.allclose(lambda_pred, lambda_actual, atol=1e-10)
        print(f"  n={n:2d} (k={k}): match={match}")
        if not match:
            print(f"    predicted: {[round(e, 6) for e in lambda_pred]}")
            print(f"    actual:    {[round(e, 6) for e in lambda_actual]}")

    # Fibonacci determinant verification
    print(f"\nFibonacci determinant: det(M) = (-1)^{{n(n-1)/2}} F(n+1)")
    fib = [1, 1]
    for i in range(20):
        fib.append(fib[-1] + fib[-2])

    for n in range(3, 16, 2):
        k = (n - 1) // 2
        M = make_M_transitive(n)
        det_actual = int(round(np.linalg.det(M)))
        det_fib = (-1)**(n*(n-1)//2) * fib[n+1]

        # Also compute from formula
        det_formula = 1  # lambda=1 contributes
        for j in range(1, k+1):
            mu = 3 + 2*np.cos(j*np.pi/(k+1))
            det_formula *= -mu  # lambda pair ±sqrt(mu) contributes -mu
        det_formula = int(round(det_formula))

        print(f"  n={n:2d}: det={det_actual:8d}, Fib={det_fib:8d}, formula={det_formula:8d}, "
              f"match={det_actual == det_fib == det_formula}")
