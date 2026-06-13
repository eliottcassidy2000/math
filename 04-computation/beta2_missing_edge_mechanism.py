#!/usr/bin/env python3
"""WHY does removing an edge create β_2?
Focus on the simplest case: 1 missing edge from a tournament at n=5.

For each of the 40 oriented graphs with β_2>0 and 1 missing pair,
find the explicit 2-cycle and understand what tournament structure kills it."""
import numpy as np
from itertools import combinations
import sys
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
    path_betti_numbers
)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
n_pairs = len(pairs)

print("=" * 70)
print("MECHANISM: HOW MISSING EDGES CREATE β_2")
print("=" * 70)

# Find all 1-missing oriented graphs with β_2 > 0
examples = []
for mask in range(3**n_pairs):
    A = [[0]*n for _ in range(n)]
    m = mask
    for i, (a, b) in enumerate(pairs):
        choice = m % 3
        m //= 3
        if choice == 1:
            A[a][b] = 1
        elif choice == 2:
            A[b][a] = 1

    missing = [(i,j) for i in range(n) for j in range(i+1,n)
               if A[i][j]==0 and A[j][i]==0]
    if len(missing) != 1:
        continue

    betti = path_betti_numbers(A, n, max_dim=3)
    if betti[2] > 0:
        examples.append((A, betti, missing[0]))

print(f"Found {len(examples)} examples with 1 missing pair and β_2 > 0")

# Analyze first few in detail
for idx, (A, betti, (miss_a, miss_b)) in enumerate(examples[:5]):
    print(f"\n{'─'*60}")
    print(f"Example {idx+1}: missing pair ({miss_a},{miss_b}), β = {betti}")
    for i in range(n):
        out = [j for j in range(n) if A[i][j]]
        print(f"  {i} → {out}")

    # Find the 2-cycle
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a1_list = [tuple(p) for p in a1]
    a2_list = [tuple(p) for p in a2]
    a3_list = [tuple(p) for p in a3]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0

    bd2 = build_full_boundary_matrix(a2_list, a1_list)
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim = dim_om2 - rank2

    U, S, Vt = np.linalg.svd(bd2_om, full_matrices=True)
    ker_basis = (om2 @ Vt[rank2:, :].T).T

    print(f"\n  dim(Ω_2) = {dim_om2}, rank(∂_2) = {rank2}, ker(∂_2|Ω_2) = {ker_dim}")
    print(f"  |A_2| = {len(a2_list)}, |A_3| = {len(a3_list)}")

    for ci in range(ker_dim):
        cyc = ker_basis[ci]
        print(f"\n  2-cycle #{ci+1}:")
        for j in range(len(a2_list)):
            if abs(cyc[j]) > 1e-8:
                p = a2_list[j]
                # Classify
                is_tt = A[p[0]][p[2]] == 1
                tag = "TT" if is_tt else "non-TT"
                # Check: does this path use the missing pair?
                uses_missing = (set([miss_a, miss_b]) <= set(p))
                miss_tag = " ★USES MISSING★" if uses_missing else ""
                print(f"    {cyc[j]:+.4f} * {p} [{tag}]{miss_tag}")

    # What Ω_3 elements exist?
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0
    print(f"\n  dim(Ω_3) = {dim_om3}")

    # Compare with both tournament completions
    for direction in [0, 1]:
        B = [row[:] for row in A]
        if direction == 0:
            B[miss_a][miss_b] = 1
        else:
            B[miss_b][miss_a] = 1
        b2 = enumerate_allowed_paths(B, n, 2)
        b3 = enumerate_allowed_paths(B, n, 3)
        b2_list = [tuple(p) for p in b2]
        b3_list = [tuple(p) for p in b3]

        # What new 2-paths appear?
        new_2paths = set(b2_list) - set(a2_list)
        # What new 3-paths appear?
        new_3paths = set(b3_list) - set(a3_list)

        edge_str = f"{miss_a}→{miss_b}" if direction==0 else f"{miss_b}→{miss_a}"
        betti_B = path_betti_numbers(B, n, max_dim=3)
        print(f"\n  Adding edge {edge_str}: β={betti_B}")
        print(f"    New 2-paths: {len(new_2paths)} — {sorted(new_2paths)[:6]}{'...' if len(new_2paths)>6 else ''}")
        print(f"    New 3-paths: {len(new_3paths)}")

        # The new 3-paths should provide the boundary that kills the 2-cycle
        if new_3paths:
            # Check: do the new 3-paths span the missing boundary?
            b_om2 = compute_omega_basis(B, n, 2, b2, enumerate_allowed_paths(B, n, 1))
            b_om3 = compute_omega_basis(B, n, 3, b3, b2)
            print(f"    dim(Ω_3|tournament) = {b_om3.shape[1] if b_om3.ndim==2 else 0} vs dim(Ω_3|oriented) = {dim_om3}")

# Summary statistics
print(f"\n\n{'='*70}")
print("STRUCTURE OF ALL 40 EXAMPLES")
print("="*70)

out_deg_types = Counter()
in_deg_types = Counter()
relationship_types = Counter()

for A, betti, (miss_a, miss_b) in examples:
    # What are the out-degrees of miss_a and miss_b?
    out_a = sum(A[miss_a])
    out_b = sum(A[miss_b])
    in_a = sum(A[j][miss_a] for j in range(n))
    in_b = sum(A[j][miss_b] for j in range(n))

    # What vertices are they connected to?
    shared_out = [v for v in range(n) if v != miss_a and v != miss_b
                  and A[miss_a][v] and A[miss_b][v]]
    shared_in = [v for v in range(n) if v != miss_a and v != miss_b
                 and A[v][miss_a] and A[v][miss_b]]

    key = (out_a, out_b, in_a, in_b, len(shared_out), len(shared_in))
    relationship_types[key] += 1

print(f"  (out_a, out_b, in_a, in_b, shared_out, shared_in): count")
for k in sorted(relationship_types):
    print(f"    {k}: {relationship_types[k]}")

print("\nDone.")
