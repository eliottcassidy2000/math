#!/usr/bin/env python3
"""
omega_vs_H_walsh.py — opus-2026-03-13-S70

Investigate the difference between H (unrestricted HP count) and Omega_{n-1}
(regular HP count) in Walsh terms.

KEY QUESTION: Omega_{n-1} ≤ H, and they differ by the regularity condition.
The Walsh degree of H is {0, 2, ..., n-1} for odd n.
The Walsh degree of Omega_{n-1} might include HIGHER degrees.

Also: verify chi = sum (-1)^m Omega_m and compare with known Walsh structure.
"""

import numpy as np
from math import comb
from itertools import combinations, permutations
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

def hamilton_count(A):
    n = A.shape[0]
    count = 0
    for perm in permutations(range(n)):
        if all(A[perm[i]][perm[i+1]] == 1 for i in range(n-1)):
            count += 1
    return count

def regular_hp_count(A):
    """Count regular Hamiltonian paths = Omega_{n-1}."""
    n = A.shape[0]
    count = 0
    for perm in permutations(range(n)):
        # Check path condition
        if not all(A[perm[i]][perm[i+1]] == 1 for i in range(n-1)):
            continue
        # Check regularity: v_{i-1} → v_{i+1} for all 0 < i < n-1
        regular = True
        for i in range(1, n-1):
            if A[perm[i-1]][perm[i+1]] != 1:
                regular = False
                break
        if regular:
            count += 1
    return count

def walsh_coefficient(f_vals, S_set, num_pairs):
    val = 0
    for T_bits, f_val in f_vals.items():
        sign = 1
        for idx in S_set:
            if (T_bits >> idx) & 1:
                sign *= -1
        val += sign * f_val
    return val / (2**num_pairs)

# ============================================================
# n = 5 analysis
# ============================================================
n = 5
num_pairs = n*(n-1)//2
total = 2**num_pairs

print("="*70)
print(f"OMEGA_{n-1} vs H at n={n}")
print("="*70)

H_vals = {}
omega_last_vals = {}
diff_vals = {}

for bits in range(total):
    A = adj_matrix(n, bits)
    H = hamilton_count(A)
    omega_last = regular_hp_count(A)
    H_vals[bits] = H
    omega_last_vals[bits] = omega_last
    diff_vals[bits] = H - omega_last

print(f"  H range: [{min(H_vals.values())}, {max(H_vals.values())}], mean={np.mean(list(H_vals.values())):.4f}")
print(f"  Omega_{n-1} range: [{min(omega_last_vals.values())}, {max(omega_last_vals.values())}], mean={np.mean(list(omega_last_vals.values())):.4f}")
print(f"  Diff H-Omega range: [{min(diff_vals.values())}, {max(diff_vals.values())}]")
print(f"  Diff distribution: {dict(sorted(Counter(diff_vals.values()).items()))}")

# Walsh degree analysis of the DIFFERENCE H - Omega_{n-1}
print(f"\n  Walsh degree of H - Omega_{n-1}:")
for deg in range(num_pairs + 1):
    combos = list(combinations(range(num_pairs), deg))
    nonzero = 0
    max_coeff = 0
    for S in combos:
        coeff = walsh_coefficient(diff_vals, set(S), num_pairs)
        if abs(coeff) > 1e-10:
            nonzero += 1
            max_coeff = max(max_coeff, abs(coeff))
    if nonzero > 0:
        print(f"    degree {deg}: {nonzero}/{len(combos)} nonzero, max |coeff| = {max_coeff:.6f}")

# Compare Walsh degrees of H and Omega_{n-1}
print(f"\n  Walsh degree comparison:")
for deg in range(min(9, num_pairs+1)):
    combos = list(combinations(range(num_pairs), deg))
    H_nonzero = 0
    O_nonzero = 0
    for S in combos:
        h_coeff = walsh_coefficient(H_vals, set(S), num_pairs)
        o_coeff = walsh_coefficient(omega_last_vals, set(S), num_pairs)
        if abs(h_coeff) > 1e-10:
            H_nonzero += 1
        if abs(o_coeff) > 1e-10:
            O_nonzero += 1
    print(f"    degree {deg}: H has {H_nonzero}, Omega_{n-1} has {O_nonzero}")

# What is the regularity condition in Walsh terms?
# H path (v_0,...,v_{n-1}) is regular iff v_{i-1}→v_{i+1} for 0<i<n-1.
# These are skip-one edges. The regularity involves edges (v_{i-1},v_{i+1}).
# For n=5: skip edges are (v_0,v_2), (v_1,v_3), (v_2,v_4) — 3 additional edges.
# A Hamiltonian path uses 4 edges, and regularity adds 3 skip edges = 7 edges total.
# So a regular HP involves 7 of the 10 edges, giving potential Walsh degree up to 7.
# But we see degree 6, not 7. So some cancellation occurs.

print(f"\n  NOTE: H involves 4 edges per HP (degree 4 max)")
print(f"  Omega_{n-1} involves 4 + 3 skip edges = 7 edges per regular HP")
print(f"  But Walsh degree of Omega_{n-1} is max 6 (one degree of cancellation)")

# Examine the degree-6 Walsh coefficients of Omega_{n-1}
print(f"\n  Degree-6 Walsh coefficients of Omega_{n-1}:")
for S in combinations(range(num_pairs), 6):
    coeff = walsh_coefficient(omega_last_vals, set(S), num_pairs)
    if abs(coeff) > 1e-10:
        # Decode edges
        edges = []
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if idx in S:
                    edges.append((i,j))
                idx += 1
        print(f"    S={set(S)}, edges={edges}, coeff={coeff:.6f}")

# Check the graph structure of these degree-6 modes
print(f"\n  Structural analysis of degree-6 modes:")
for S in combinations(range(num_pairs), 6):
    coeff = walsh_coefficient(omega_last_vals, set(S), num_pairs)
    if abs(coeff) > 1e-10:
        edges = []
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if idx in S:
                    edges.append((i,j))
                idx += 1
        # Vertex degrees
        deg_count = {}
        for i, j in edges:
            deg_count[i] = deg_count.get(i, 0) + 1
            deg_count[j] = deg_count.get(j, 0) + 1
        vertices_used = sorted(deg_count.keys())
        complement_edges = set(range(num_pairs)) - set(S)
        # 10 total edges, 6 used, so complement has 4 edges
        comp_edges = []
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if idx in complement_edges:
                    comp_edges.append((i,j))
                idx += 1
        print(f"    edges={edges}, complement={comp_edges}")
        print(f"      degrees={deg_count}, vertices={vertices_used}")

print("\nDONE.")
