#!/usr/bin/env python3
"""
n=8 tournament path homology and β_2=0 verification at n=6,7.

Also: complement duality at n=6.
"""
import numpy as np
from itertools import combinations
from collections import Counter
import random
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    path_betti_numbers, enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix, count_3cycles, ham_path_count
)

random.seed(42)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

# ===== β_2=0 verification at n=6,7 =====
print("=" * 70)
print("β_2=0 VERIFICATION: ker(∂_2) = im(∂_3)?")
print("=" * 70)

for n in [6, 7]:
    print(f"\nn={n}: sampling 100 tournaments")
    all_match = True
    for trial in range(100):
        A = random_tournament(n)

        a2 = enumerate_allowed_paths(A, n, 2)
        a1 = enumerate_allowed_paths(A, n, 1)
        a3 = enumerate_allowed_paths(A, n, 3)

        om2 = compute_omega_basis(A, n, 2, a2, a1)
        om3 = compute_omega_basis(A, n, 3, a3, a2)

        dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
        dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

        if dim_om2 > 0:
            bd2 = build_full_boundary_matrix(a2, a1)
            bd2_om = bd2 @ om2
            U, S, Vt = np.linalg.svd(bd2_om, full_matrices=True)
            rank2 = sum(s > 1e-8 for s in S)
            ker2 = dim_om2 - rank2

            if dim_om3 > 0:
                bd3 = build_full_boundary_matrix(a3, a2)
                bd3_om = bd3 @ om3
                rank3 = np.linalg.matrix_rank(bd3_om, tol=1e-8)
            else:
                rank3 = 0

            if ker2 != rank3:
                print(f"  MISMATCH at trial {trial}: ker(∂_2)={ker2}, im(∂_3)={rank3}")
                all_match = False

        if trial % 20 == 19:
            print(f"  ... {trial+1} done", flush=True)

    if all_match:
        print(f"  ker(∂_2) = im(∂_3) for ALL {100 if n < 8 else 30} trials ✓")

# ===== Complement duality at n=6 =====
print("\n" + "=" * 70)
print("COMPLEMENT DUALITY n=6")
print("=" * 70)

n = 6
matches = 0
total = 0
for trial in range(100):
    A = random_tournament(n)
    A_op = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                A_op[i][j] = 1 - A[i][j]
    b1 = path_betti_numbers(A, n, max_dim=5)
    b2 = path_betti_numbers(A_op, n, max_dim=5)
    total += 1
    if b1 == b2:
        matches += 1
    elif total <= 3:
        print(f"  MISMATCH: β(T)={b1}, β(T^op)={b2}")
    if trial % 20 == 19:
        print(f"  ... {trial+1} done", flush=True)

print(f"  β(T) = β(T^op): {matches}/{total}")

# ===== n=8: smaller sample =====
print("\n" + "=" * 70)
print("n=8 TOURNAMENTS (10 samples)")
print("=" * 70)

n = 8
betti_dist = Counter()
for trial in range(10):
    A = random_tournament(n)
    t3 = count_3cycles(A, n)
    betti = path_betti_numbers(A, n, max_dim=7)
    bt = tuple(betti)
    betti_dist[bt] += 1
    print(f"  Trial {trial}: t3={t3}, β={betti}", flush=True)

print(f"\n  Betti distribution:")
for bt in sorted(betti_dist.keys()):
    print(f"    β={list(bt)}: {betti_dist[bt]}")

print("\nDone.")
