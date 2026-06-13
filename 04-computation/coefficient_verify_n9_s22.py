#!/usr/bin/env python3
"""
coefficient_verify_n9_s22.py — opus-2026-04-04-S22

VERIFY the coefficient formula a(n) = (n-2)!/2^{n-4} at n=8→9.
Prediction: a(9) = 7!/2^5 = 5040/32 = 157.5.

Also verify the GROWTH RATIO formula: a(n)/a(n-1) = (n-2)/2.
At n=9: growth ratio = 7/2 = 3.5. So a(9) = 45 × 3.5 = 157.5.

This is the deep structure: the simple coefficient in ΔH ≈ a(n)·Δc₃
is EXACTLY (n-2)!/2^{n-4}, arising from the completeness of tournaments.
"""

import numpy as np
from collections import defaultdict
import random
import time
import sys

random.seed(42)

def compute_H(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << n) - 1])

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def extend(A_sub, n_sub, sigma):
    n = n_sub + 1
    A = [[0]*n for _ in range(n)]
    for i in range(n_sub):
        for j in range(n_sub):
            A[i][j] = A_sub[i][j]
    for i in range(n_sub):
        if sigma[i]:
            A[n_sub][i] = 1
        else:
            A[i][n_sub] = 1
    return A

def delta_c3(A_sub, n_sub, sigma):
    return sum(A_sub[i][j] for i in range(n_sub) if sigma[i] == 1
                           for j in range(n_sub) if sigma[j] == 0)


def verify_coefficient(n_target, n_parents=300, n_ext=8):
    """Measure the simple coefficient at level n-1 → n."""
    n_sub = n_target - 1
    predicted = 1
    for k in range(4, n_target + 1):
        predicted = predicted * (k - 2) / 2
    # Also: predicted = factorial(n_target - 2) / 2**(n_target - 4)
    from math import factorial
    predicted_formula = factorial(n_target - 2) / (2 ** (n_target - 4))

    print(f"\n{'='*70}")
    print(f"  COEFFICIENT VERIFICATION at n={n_sub} → n={n_target}")
    print(f"  Predicted: a = (n-2)!/2^{{n-4}} = {n_target-2}!/2^{n_target-4} = {predicted_formula:.2f}")
    print(f"  Growth ratio: (n-2)/2 = {(n_target-2)/2:.2f}")
    print(f"  Sampling: {n_parents} parents × {n_ext} extensions = {n_parents*n_ext} data points")
    print(f"{'='*70}")

    data = []
    t0 = time.time()
    for p in range(n_parents):
        A_sub = random_tournament(n_sub)
        H_sub = compute_H(A_sub, n_sub)
        for _ in range(n_ext):
            sigma = tuple(random.randint(0,1) for _ in range(n_sub))
            A_n = extend(A_sub, n_sub, sigma)
            H_n = compute_H(A_n, n_target)
            dc3 = delta_c3(A_sub, n_sub, sigma)
            data.append((H_sub, H_n - H_sub, dc3))
        if (p+1) % 100 == 0:
            print(f"    ... {p+1}/{n_parents} parents ({time.time()-t0:.1f}s)", file=sys.stderr)

    H_sub_arr = np.array([d[0] for d in data], dtype=float)
    dH_arr = np.array([d[1] for d in data], dtype=float)
    dc3_arr = np.array([d[2] for d in data], dtype=float)
    N = len(data)

    # Simple regression
    X = np.column_stack([dc3_arr, np.ones(N)])
    c, _, _, _ = np.linalg.lstsq(X, dH_arr, rcond=None)
    R2 = 1 - np.sum((dH_arr - X@c)**2) / np.sum((dH_arr - dH_arr.mean())**2)

    print(f"\n  Measured: a = {c[0]:.4f} (intercept = {c[1]:.4f}), R²={R2:.4f}")
    print(f"  Predicted: a = {predicted_formula:.4f}")
    print(f"  Ratio measured/predicted = {c[0]/predicted_formula:.6f}")

    # Interaction model
    X3 = np.column_stack([dc3_arr * H_sub_arr, dc3_arr, H_sub_arr, np.ones(N)])
    c3, _, _, _ = np.linalg.lstsq(X3, dH_arr, rcond=None)
    R3 = 1 - np.sum((dH_arr - X3@c3)**2) / np.sum((dH_arr - dH_arr.mean())**2)
    print(f"  Interaction: a_inter = {c3[0]:.6f}, R²={R3:.4f}")

    return c[0], predicted_formula, R2


if __name__ == "__main__":
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  COEFFICIENT FORMULA VERIFICATION: a(n) = (n-2)!/2^{n-4}      ║")
    print("║  opus-2026-04-04-S22                                           ║")
    print("╚══════════════════════════════════════════════════════════════════╝")

    t0 = time.time()

    # Verify at n=8→9
    a_meas, a_pred, R2 = verify_coefficient(9, n_parents=200, n_ext=8)

    # Also spot-check n=7→8 and n=6→7 with the formula
    for n_target in [7, 8]:
        from math import factorial
        pred = factorial(n_target - 2) / (2 ** (n_target - 4))
        print(f"\n  Spot-check n={n_target}: predicted a = {pred:.4f}")

    # Print the full table
    print(f"\n{'='*70}")
    print(f"  COMPLETE COEFFICIENT TABLE")
    print(f"{'='*70}")
    from math import factorial
    print(f"  {'n':>3}  {'(n-2)!':>10}  {'2^(n-4)':>10}  {'a = ratio':>12}  {'growth':>8}")
    for n in range(4, 12):
        f = factorial(n-2)
        p = 2**(n-4)
        a = f / p
        if n > 4:
            growth = f"×{(n-2)/2:.1f}"
        else:
            growth = ""
        print(f"  {n:3d}  {f:10d}  {p:10d}  {a:12.2f}  {growth:>8}")

    print(f"\n  Asymptotic: a(n) ~ sqrt(2π(n-2)) · ((n-2)/(2e))^(n-2)")
    print(f"  Since (n-2)/(2e) > 1 for n ≥ 8, the coefficient grows SUPER-EXPONENTIALLY.")

    t1 = time.time()
    print(f"\n  Total time: {t1-t0:.1f}s")
