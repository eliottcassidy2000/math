#!/usr/bin/env python3
"""
walsh_omega_degree.py — opus-2026-03-13-S70

Test the Walsh-degree stratification of Omega_m:
  Omega_m has nonzero Walsh coefficients ONLY at degree 2*(m-1) (for m >= 2).

At n=5: confirmed Omega_2 → deg 2, Omega_3 → deg 4.
Now test at n=6: expect Omega_2 → deg 2, Omega_3 → deg 4, Omega_4 → deg 6(?).

Also analyze what GRAPH STRUCTURES the nonzero Walsh modes correspond to.
"""

import numpy as np
from math import comb
from itertools import combinations, permutations
from collections import Counter, defaultdict

def adj_matrix(n, T_bits):
    """Build adjacency matrix from bit encoding of tournament."""
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
    """Count regular m-paths in tournament A.
    Regular m-path: (v_0,...,v_m) with v_i→v_{i+1} AND v_{i-1}→v_{i+1} for 0<i<m.
    """
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
                # Check regularity: if depth >= 1, need A[path[-2]][v] == 1
                if depth >= 1 and A[path[-2]][v] != 1:
                    continue
                path.append(v)
                dfs(path, depth + 1)
                path.pop()

    for start in range(n):
        dfs([start], 0)

    return count

def count_3cycles(A):
    """Count 3-cycles in tournament."""
    n = A.shape[0]
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                    count += 1
    return count

def walsh_coefficient(f_vals, S_set, n):
    """Compute Walsh coefficient hat_f[S] for a specific set S of edge indices."""
    num_pairs = n*(n-1)//2
    val = 0
    for T_bits, f_val in f_vals.items():
        sign = 1
        for idx in S_set:
            if (T_bits >> idx) & 1:
                sign *= -1
        val += sign * f_val
    return val / (2**num_pairs)

def edge_index(n, i, j):
    """Get the edge index for pair (i,j) with i<j."""
    idx = 0
    for a in range(n):
        for b in range(a+1, n):
            if a == i and b == j:
                return idx
            idx += 1
    return -1

# ============================================================
# n = 5: Verify known results
# ============================================================
print("="*70)
print("n = 5: WALSH DEGREE OF Omega_m")
print("="*70)

n = 5
num_pairs = n*(n-1)//2

omega_vals = {m: {} for m in range(n)}
for bits in range(2**num_pairs):
    A = adj_matrix(n, bits)
    for m in range(n):
        omega_vals[m][bits] = count_regular_m_paths(A, m)

# Check which Walsh degrees are nonzero for each Omega_m
for m in range(n):
    degrees_present = set()
    for deg in range(num_pairs + 1):
        # Sample a few subsets of this degree
        found_nonzero = False
        for S in combinations(range(num_pairs), deg):
            coeff = walsh_coefficient(omega_vals[m], set(S), n)
            if abs(coeff) > 1e-10:
                found_nonzero = True
                break
        if found_nonzero:
            degrees_present.add(deg)
    print(f"  Omega_{m}: nonzero Walsh degrees = {sorted(degrees_present)}")

# ============================================================
# n = 6: Full analysis
# ============================================================
print(f"\n{'='*70}")
print("n = 6: WALSH DEGREE OF Omega_m")
print("="*70)

n = 6
num_pairs = n*(n-1)//2
total = 2**num_pairs
print(f"  Total tournaments: {total}")

# Compute Omega for all tournaments
print("  Computing Omega_m for all tournaments...")
omega_vals_6 = {m: {} for m in range(n)}

for bits in range(total):
    A = adj_matrix(n, bits)
    for m in range(n):
        omega_vals_6[m][bits] = count_regular_m_paths(A, m)

    if bits % 10000 == 0 and bits > 0:
        print(f"    {bits}/{total} done...")

print("  Done.")

# Statistics
for m in range(n):
    vals = list(omega_vals_6[m].values())
    print(f"  Omega_{m}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.2f}, #unique={len(set(vals))}")

# Walsh degree analysis: check each degree
print(f"\n  Walsh degree analysis:")
for m in range(n):
    degrees_present = set()
    for deg in range(num_pairs + 1):
        # Check if ANY coefficient at this degree is nonzero
        # For efficiency, sample randomly if C(num_pairs, deg) is large
        total_combos = comb(num_pairs, deg)
        if total_combos <= 500:
            check_all = True
            combos_to_check = list(combinations(range(num_pairs), deg))
        else:
            check_all = False
            # Random sample of 200 subsets
            import random
            random.seed(42 + deg + m*100)
            all_combos = list(combinations(range(num_pairs), deg))
            combos_to_check = random.sample(all_combos, min(200, len(all_combos)))

        found_nonzero = False
        for S in combos_to_check:
            coeff = walsh_coefficient(omega_vals_6[m], set(S), n)
            if abs(coeff) > 1e-10:
                found_nonzero = True
                break

        if found_nonzero:
            degrees_present.add(deg)
            if not check_all:
                print(f"    Omega_{m}: degree {deg} PRESENT (sampled)")
        elif not check_all:
            # Might have missed it — flag as uncertain
            pass

    print(f"  Omega_{m}: nonzero Walsh degrees (checked) = {sorted(degrees_present)}")

# Cross-check: chi = sum (-1)^m Omega_m
print(f"\n  Chi verification:")
chi_6 = {}
for bits in range(total):
    chi_6[bits] = sum((-1)**m * omega_vals_6[m][bits] for m in range(n))

chi_vals = sorted(set(chi_6.values()))
print(f"    chi values: {chi_vals}")

# Verify Omega_2 = C(n,3) - t_3
print(f"\n  Omega_2 = C({n},3) - t_3 check:")
mismatches = 0
for bits in range(min(1000, total)):
    A = adj_matrix(n, bits)
    t3 = count_3cycles(A)
    expected = comb(n, 3) - t3
    if omega_vals_6[2][bits] != expected:
        mismatches += 1
print(f"    {mismatches}/1000 mismatches")

# Omega_3 correlation analysis
print(f"\n  Correlation analysis:")
o3 = [omega_vals_6[3][b] for b in range(total)]
o4 = [omega_vals_6[4][b] for b in range(total)]
H = [omega_vals_6[5][b] for b in range(total)]

# t3 for all
t3_all = []
for bits in range(total):
    A = adj_matrix(n, bits)
    t3_all.append(count_3cycles(A))

print(f"    corr(t3, Omega_3) = {np.corrcoef(t3_all, o3)[0,1]:.6f}")
print(f"    corr(t3, Omega_4) = {np.corrcoef(t3_all, o4)[0,1]:.6f}")
print(f"    corr(Omega_3, Omega_4) = {np.corrcoef(o3, o4)[0,1]:.6f}")
print(f"    corr(H, Omega_3) = {np.corrcoef(H, o3)[0,1]:.6f}")
print(f"    corr(H, Omega_4) = {np.corrcoef(H, o4)[0,1]:.6f}")

# Check if Omega_3 is pure degree 4
print(f"\n  Omega_3 detailed degree check:")
for deg in [0, 2, 4, 6, 8]:
    combos = list(combinations(range(num_pairs), deg))
    nonzero_count = 0
    max_abs = 0
    for S in combos:
        coeff = walsh_coefficient(omega_vals_6[3], set(S), n)
        if abs(coeff) > 1e-10:
            nonzero_count += 1
            max_abs = max(max_abs, abs(coeff))
    print(f"    degree {deg}: {nonzero_count}/{len(combos)} nonzero (max |coeff| = {max_abs:.6f})")

# Check if Omega_4 is pure degree 6
print(f"\n  Omega_4 detailed degree check:")
for deg in [0, 2, 4, 6, 8]:
    combos = list(combinations(range(num_pairs), deg))
    nonzero_count = 0
    max_abs = 0
    for S in combos:
        coeff = walsh_coefficient(omega_vals_6[4], set(S), n)
        if abs(coeff) > 1e-10:
            nonzero_count += 1
            max_abs = max(max_abs, abs(coeff))
    print(f"    degree {deg}: {nonzero_count}/{len(combos)} nonzero (max |coeff| = {max_abs:.6f})")

print("\nDONE.")
