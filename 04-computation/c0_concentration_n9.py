#!/usr/bin/env python3
"""
c_0 CONCENTRATION AT n=9 — opus-2026-03-06-S11b (continued³)

EXTENSION OF c_0 CONCENTRATION THEOREM TO n=9:

At n=7 regular: c_2=9I, c_4=180I, c_6=720I ALL universal. Only c_0 varies.
At n=9 regular: c_4=720I, c_6=10080I, c_8=40320I universal.
  BUT c_2 IS NOT UNIVERSAL (varies between regular tournaments).
  Both c_0 and c_2 vary.

UNIVERSAL COEFFICIENT TABLE (per-vertex diagonal values):
  n=5: c_4 = 4! = 24
  n=7: c_6 = 6! = 720,  c_4 = 180 = 6!/4,  c_2 = 9
  n=9: c_8 = 8! = 40320, c_6 = 10080 = 8!/4, c_4 = 720 = 6!

PATTERN:
  c_{n-1} = (n-1)!  (always)
  c_{n-3} = (n-1)!/4  (verified n=7,9)

The BOUNDARY between universal and non-universal coefficients:
  n=5: c_4 universal (top only)
  n=7: c_2, c_4, c_6 universal (top 3)
  n=9: c_4, c_6, c_8 universal (top 3)

The number of non-universal coefficients grows with n:
  n=5: 1 (c_0)
  n=7: 1 (c_0)
  n=9: 2 (c_0, c_2)

CONJECTURE: At n=2k+1, the non-universal coefficients are c_0, c_2, ..., c_{2j}
where the boundary j depends on n. The top ⌊(n-1)/2⌋ - j coefficients are universal.
"""

import numpy as np
import random

def random_regular_tournament(n):
    k = (n-1) // 2
    while True:
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
        if all(sum(row) == k for row in A):
            return A

def signed_position_diagonal_weighted(A, n, r_val):
    """Compute diag(M(r)) = sum_P (-1)^{pos(v,P)} * prod(weights)."""
    full = (1 << n) - 1
    dp_fwd = [[0.0]*n for _ in range(1 << n)]
    for v in range(n):
        dp_fwd[1 << v][v] = 1.0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp_fwd[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                w = r_val + (A[v][u] - 0.5)
                dp_fwd[mask | (1 << u)][u] += dp_fwd[mask][v] * w

    dp_bwd = [[0.0]*n for _ in range(1 << n)]
    for v in range(n):
        dp_bwd[1 << v][v] = 1.0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp_bwd[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                w = r_val + (A[u][v] - 0.5)
                dp_bwd[mask | (1 << u)][u] += dp_bwd[mask][v] * w

    diag = [0.0] * n
    for v in range(n):
        for mask_before in range(1 << n):
            if mask_before & (1 << v): continue
            mask_with_v = mask_before | (1 << v)
            if dp_fwd[mask_with_v][v] == 0: continue
            mask_after = full ^ mask_before
            if not (mask_after & (1 << v)): continue
            if dp_bwd[mask_after][v] == 0: continue
            k = bin(mask_before).count('1')
            diag[v] += ((-1)**k) * dp_fwd[mask_with_v][v] * dp_bwd[mask_after][v]
    return diag

if __name__ == "__main__":
    n = 9
    random.seed(42)

    print("=" * 70)
    print(f"n={n}: c_k COEFFICIENT UNIVERSALITY")
    print("=" * 70)

    r_values = [0.0, 0.3, 0.5, 0.7, 1.0]
    r_powers = np.array([[r**(2*k) for k in range(5)] for r in r_values])

    for trial in range(2):
        A = random_regular_tournament(n)
        diag_at_r = [signed_position_diagonal_weighted(A, n, r) for r in r_values]
        H = round(sum(diag_at_r[2]))

        vertex_c = np.zeros((n, 5))
        for v in range(n):
            y = [diag_at_r[i][v] for i in range(5)]
            vertex_c[v, :] = np.linalg.solve(r_powers, y)

        print(f"\nTrial {trial+1}, H={H}:")
        for k_idx, k_name in enumerate(['c_0', 'c_2', 'c_4', 'c_6', 'c_8']):
            vals = vertex_c[:, k_idx]
            print(f"  {k_name}: range=[{np.min(vals):.1f}, {np.max(vals):.1f}], "
                  f"scalar={'YES' if np.max(vals)-np.min(vals)<0.5 else 'NO'}")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)
