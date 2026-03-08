#!/usr/bin/env python3
"""
beta2_min_surplus.py - Find minimum surplus across all tournaments

Known:
  n=5: min surplus = 0 (280/1024 = 27.3%)
  n=6: min surplus = 1 (80/32768 = 0.2%)
  n=7: min surplus = 9 (sampled)

Questions:
1. What is the exact minimum surplus at n=7? (exhaustive if feasible, else heavy sampling)
2. What is the surplus for the TRANSITIVE tournament? (base case for induction)
3. Is there a formula for min surplus?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os, random
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def transitive_tournament(n):
    """The transitive (acyclic) tournament: i->j iff i < j."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    return A

def compute_surplus(A, n):
    allowed = {}
    for p in range(5):
        allowed[p] = enumerate_allowed_paths(A, n, p)
        if not allowed[p]:
            break

    omega_basis = {}
    for p in range(4):
        if p not in allowed or not allowed[p]:
            omega_basis[p] = np.zeros((0, 0))
            continue
        basis = compute_omega_basis(A, n, p, allowed[p],
                                     allowed[p-1] if p >= 1 and p-1 in allowed else [])
        omega_basis[p] = basis

    dim_O2 = omega_basis[2].shape[1] if omega_basis[2].ndim == 2 else 0
    dim_O3 = omega_basis[3].shape[1] if omega_basis[3].ndim == 2 else 0

    if dim_O2 == 0:
        Z2 = 0
    else:
        bd2 = build_full_boundary_matrix(allowed[2], allowed.get(1, []))
        bd2_om = bd2 @ omega_basis[2]
        S_v = np.linalg.svd(bd2_om, compute_uv=False)
        rk2 = int(np.sum(np.abs(S_v) > 1e-8))
        Z2 = dim_O2 - rk2

    surplus = dim_O3 - Z2
    return surplus, dim_O2, dim_O3, Z2

# Transitive tournament surplus at each n
print("=" * 70)
print("TRANSITIVE TOURNAMENT SURPLUS")
print("=" * 70)

for n in range(3, 10):
    A = transitive_tournament(n)
    surplus, O2, O3, Z2 = compute_surplus(A, n)
    print(f"  n={n}: O2={O2}, O3={O3}, Z2={Z2}, surplus={surplus}")

# Exhaustive at n=5, n=6
print(f"\n{'='*70}")
print("EXHAUSTIVE SURPLUS DISTRIBUTION")
print(f"{'='*70}")

for n in [5, 6]:
    n_arcs = n*(n-1)//2
    total = 1 << n_arcs
    t0 = time.time()

    surplus_dist = Counter()
    for bits in range(total):
        A = build_adj(n, bits)
        surplus, _, _, _ = compute_surplus(A, n)
        surplus_dist[surplus] += 1

    dt = time.time() - t0
    print(f"\n  n={n}: ({dt:.0f}s)")
    for s in sorted(surplus_dist.keys()):
        print(f"    surplus={s}: {surplus_dist[s]} ({100*surplus_dist[s]/total:.1f}%)")
    print(f"    min surplus = {min(surplus_dist.keys())}")

# Heavy sampling at n=7
print(f"\n{'='*70}")
print("HEAVY SAMPLING AT n=7")
print(f"{'='*70}")

n = 7
N_SAMPLES = 10000
t0 = time.time()
min_surplus = float('inf')
surplus_dist = Counter()

for trial in range(N_SAMPLES):
    if trial % 2000 == 0 and trial > 0:
        dt = time.time() - t0
        print(f"  ... {trial}/{N_SAMPLES} ({dt:.0f}s), min_surplus={min_surplus}")

    A = random_tournament(n)
    surplus, O2, O3, Z2 = compute_surplus(A, n)
    surplus_dist[surplus] += 1

    if surplus < min_surplus:
        min_surplus = surplus
        score = tuple(sorted(sum(row) for row in A))
        print(f"  New min surplus={surplus}: O2={O2}, O3={O3}, Z2={Z2}, score={score}")

dt = time.time() - t0
print(f"\n  n=7 ({dt:.0f}s):")
for s in sorted(surplus_dist.keys())[:15]:
    print(f"    surplus={s}: {surplus_dist[s]} ({100*surplus_dist[s]/N_SAMPLES:.1f}%)")
print(f"    min surplus = {min(surplus_dist.keys())}")

# Heavy sampling at n=8
print(f"\n{'='*70}")
print("HEAVY SAMPLING AT n=8")
print(f"{'='*70}")

n = 8
N_SAMPLES = 3000
t0 = time.time()
min_surplus = float('inf')
surplus_dist = Counter()

for trial in range(N_SAMPLES):
    if trial % 1000 == 0 and trial > 0:
        dt = time.time() - t0
        print(f"  ... {trial}/{N_SAMPLES} ({dt:.0f}s), min_surplus={min_surplus}")

    A = random_tournament(n)
    surplus, O2, O3, Z2 = compute_surplus(A, n)
    surplus_dist[surplus] += 1

    if surplus < min_surplus:
        min_surplus = surplus
        score = tuple(sorted(sum(row) for row in A))
        print(f"  New min surplus={surplus}: O2={O2}, O3={O3}, Z2={Z2}, score={score}")

dt = time.time() - t0
print(f"\n  n=8 ({dt:.0f}s):")
for s in sorted(surplus_dist.keys())[:15]:
    print(f"    surplus={s}: {surplus_dist[s]} ({100*surplus_dist[s]/N_SAMPLES:.1f}%)")
print(f"    min surplus = {min(surplus_dist.keys())}")

# Summary
print(f"\n{'='*70}")
print("SUMMARY: MINIMUM SURPLUS BY n")
print(f"{'='*70}")

print("\nDone.")
