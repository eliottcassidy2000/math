#!/usr/bin/env python3
"""
walsh_omega_fast.py — opus-2026-03-13-S70

Fast Walsh-Hadamard transform of Omega_m for n=5,6.
Uses the FFT-like WHT algorithm: O(N log N) instead of O(N^2).

Key finding to verify: Omega_m has nonzero Walsh coefficients at specific degrees:
  n=5: Omega_2 → {0,2}, Omega_3 → {0,4}, Omega_4 → {0,2,4,6}
  n=6: ???
"""

import numpy as np
from math import comb
from itertools import permutations
from collections import Counter

def adj_matrix(n, T_bits):
    A = np.zeros((n,n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (T_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_regular_m_paths(A, m):
    """Count regular m-paths."""
    n = A.shape[0]
    count = 0
    def dfs(path, depth):
        nonlocal count
        if depth == m:
            count += 1
            return
        last = path[-1]
        for v in range(n):
            if v not in path and A[last][v] == 1:
                if depth >= 1 and A[path[-2]][v] != 1:
                    continue
                path.append(v)
                dfs(path, depth + 1)
                path.pop()
    for start in range(n):
        dfs([start], 0)
    return count

def hamilton_count(A):
    n = A.shape[0]
    count = 0
    for perm in permutations(range(n)):
        if all(A[perm[i]][perm[i+1]] == 1 for i in range(n-1)):
            count += 1
    return count

def fast_wht(f_array):
    """Fast Walsh-Hadamard Transform.
    f_array[i] = f(i) for i in {0,...,2^n-1}.
    Returns hat_f[S] = (1/2^n) sum_i f(i) (-1)^{popcount(i&S)}.
    """
    N = len(f_array)
    n = N.bit_length() - 1
    a = f_array.astype(float).copy()

    # In-place Walsh-Hadamard transform
    h = 1
    while h < N:
        for i in range(0, N, h * 2):
            for j in range(i, i + h):
                x = a[j]
                y = a[j + h]
                a[j] = x + y
                a[j + h] = x - y
        h *= 2

    return a / N

def popcount(x):
    return bin(x).count('1')

def analyze_walsh_degrees(hat_f, num_pairs, label="f"):
    """Analyze which Walsh degrees are nonzero."""
    N = len(hat_f)
    degree_stats = {}
    for deg in range(num_pairs + 1):
        nonzero = 0
        max_abs = 0
        total = 0
        for S in range(N):
            if popcount(S) == deg:
                total += 1
                if abs(hat_f[S]) > 1e-10:
                    nonzero += 1
                    max_abs = max(max_abs, abs(hat_f[S]))
        if total > 0:
            degree_stats[deg] = (nonzero, total, max_abs)

    degrees_present = [d for d, (nz, _, _) in degree_stats.items() if nz > 0]
    print(f"  {label}: nonzero at degrees {degrees_present}")
    for d in degrees_present:
        nz, tot, mx = degree_stats[d]
        print(f"    deg {d:2d}: {nz:5d}/{tot:5d} nonzero, max |coeff| = {mx:.8f}")
    return degrees_present

# ============================================================
# n = 5
# ============================================================
print("="*70)
print("n = 5: FAST WALSH ANALYSIS")
print("="*70)

n = 5
num_pairs = n*(n-1)//2
N = 2**num_pairs

# Compute all function values
print(f"Computing Omega_m and H for all {N} tournaments...")
omega_arrays = {m: np.zeros(N) for m in range(n)}
H_array = np.zeros(N)

for bits in range(N):
    A = adj_matrix(n, bits)
    for m in range(n):
        omega_arrays[m][bits] = count_regular_m_paths(A, m)
    H_array[bits] = hamilton_count(A)

print("Computing Walsh transforms...")
for m in range(n):
    hat_omega = fast_wht(omega_arrays[m])
    analyze_walsh_degrees(hat_omega, num_pairs, f"Omega_{m}")

hat_H = fast_wht(H_array)
analyze_walsh_degrees(hat_H, num_pairs, "H")

# Check: H = Omega_{n-1}?
diff = H_array - omega_arrays[n-1]
print(f"\n  H - Omega_{n-1}: min={diff.min():.0f}, max={diff.max():.0f}, mean={diff.mean():.4f}")
if diff.max() == 0 and diff.min() == 0:
    print(f"  H = Omega_{n-1} at n={n} ✓")
else:
    print(f"  H ≠ Omega_{n-1} at n={n} ✗")
    hat_diff = fast_wht(diff)
    analyze_walsh_degrees(hat_diff, num_pairs, "H-Omega_{n-1}")

# Chi
chi_array = sum((-1)**m * omega_arrays[m] for m in range(n))
hat_chi = fast_wht(chi_array)
analyze_walsh_degrees(hat_chi, num_pairs, "chi")

# ============================================================
# n = 6
# ============================================================
print(f"\n{'='*70}")
print("n = 6: FAST WALSH ANALYSIS")
print("="*70)

n = 6
num_pairs = n*(n-1)//2
N = 2**num_pairs

print(f"Computing Omega_m and H for all {N} tournaments...")
omega_arrays = {m: np.zeros(N) for m in range(n)}
H_array = np.zeros(N)

for bits in range(N):
    A = adj_matrix(n, bits)
    for m in range(n):
        omega_arrays[m][bits] = count_regular_m_paths(A, m)
    H_array[bits] = hamilton_count(A)
    if bits % 5000 == 0 and bits > 0:
        print(f"  {bits}/{N} done...")

print("Computing Walsh transforms...")
for m in range(n):
    hat_omega = fast_wht(omega_arrays[m])
    analyze_walsh_degrees(hat_omega, num_pairs, f"Omega_{m}")

hat_H = fast_wht(H_array)
analyze_walsh_degrees(hat_H, num_pairs, "H")

# H vs Omega_{n-1}
diff = H_array - omega_arrays[n-1]
print(f"\n  H - Omega_{n-1}: min={diff.min():.0f}, max={diff.max():.0f}, mean={diff.mean():.4f}")
hat_diff = fast_wht(diff)
analyze_walsh_degrees(hat_diff, num_pairs, "H-Omega_{n-1}")

# Chi
chi_array = sum((-1)**m * omega_arrays[m] for m in range(n))
hat_chi = fast_wht(chi_array)
analyze_walsh_degrees(hat_chi, num_pairs, "chi")

# Cross-correlations
print(f"\n  Cross-correlations:")
for m1 in range(2, n):
    for m2 in range(m1+1, n):
        corr = np.corrcoef(omega_arrays[m1], omega_arrays[m2])[0,1]
        print(f"    corr(Omega_{m1}, Omega_{m2}) = {corr:.6f}")

# Omega_m distribution
print(f"\n  Omega_m distributions:")
for m in range(n):
    vals = omega_arrays[m]
    print(f"    Omega_{m}: unique values = {sorted(set(vals.astype(int)))}")

# Verify Omega_2 = C(n,3) - t_3
print(f"\n  Omega_2 formula check:")
from collections import defaultdict
t3_array = np.zeros(N)
for bits in range(min(1000, N)):
    A = adj_matrix(n, bits)
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                    t3 += 1
    t3_array[bits] = t3

mismatches = sum(1 for bits in range(1000) if omega_arrays[2][bits] != comb(n,3) - t3_array[bits])
print(f"    Omega_2 = C({n},3) - t_3: {mismatches}/1000 mismatches")

print("\nDONE.")
