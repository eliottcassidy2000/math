#!/usr/bin/env python3
"""
beta2_transitive_perturbation.py - Analyze surplus change from transitive under arc flips

For the transitive tournament T_n (total order 0 < 1 < ... < n-1):
  - O_2 = C(n,3) (every 3-vertex subset gives a TT 2-path)
  - O_3 = C(n,4) (every 4-vertex subset gives a DT 3-path)
  - Z_2 = C(n-1,3) (how? let's verify)
  - surplus = C(n-1,4) (verified above)

Flipping arc (u->v) where u < v in the transitive order:
  - This creates a "backward" arc, breaking transitivity at exactly one pair
  - The effect on surplus depends on how many paths use this arc

Questions:
1. What's the exact surplus formula for T_n with one arc flipped?
2. Is there a formula for the minimum surplus after one flip?
3. Can we prove surplus >= 0 for ALL tournaments by induction on flips from T_n?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def transitive_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    return A

def compute_surplus_full(A, n):
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


print("=" * 70)
print("TRANSITIVE TOURNAMENT PERTURBATION ANALYSIS")
print("=" * 70)

for n in [5, 6, 7, 8]:
    print(f"\n{'='*70}")
    print(f"n = {n}")
    print(f"{'='*70}")

    A = transitive_tournament(n)
    surplus0, O2_0, O3_0, Z2_0 = compute_surplus_full(A, n)
    print(f"  Transitive: O2={O2_0}, O3={O3_0}, Z2={Z2_0}, surplus={surplus0}")

    # Try flipping each arc
    print(f"\n  Single arc flips from transitive:")
    flip_data = []
    for u in range(n):
        for v in range(u+1, n):
            # Flip u->v to v->u
            B = [row[:] for row in A]
            B[u][v] = 0
            B[v][u] = 1

            surplus, O2, O3, Z2 = compute_surplus_full(B, n)
            delta = surplus - surplus0
            flip_data.append({
                'u': u, 'v': v, 'gap': v-u,
                'surplus': surplus, 'delta': delta,
                'O2': O2, 'O3': O3, 'Z2': Z2
            })

    # Group by gap = v - u
    by_gap = defaultdict(list)
    for d in flip_data:
        by_gap[d['gap']].append(d)

    for gap in sorted(by_gap.keys()):
        surpluses = sorted(set(d['surplus'] for d in by_gap[gap]))
        deltas = sorted(set(d['delta'] for d in by_gap[gap]))
        count = len(by_gap[gap])
        # Are all flips with same gap equivalent?
        unique_surplus = len(set(d['surplus'] for d in by_gap[gap]))
        print(f"    gap={gap}: {count} flips, surplus={surpluses}, delta={deltas}")

    min_s = min(d['surplus'] for d in flip_data)
    print(f"  Min surplus after 1 flip: {min_s}")

    # For n<=7, try 2 flips from transitive
    if n <= 7:
        print(f"\n  Two consecutive arc flips from transitive:")
        min_s2 = float('inf')
        all_s2 = Counter()
        for fd1 in flip_data:
            # Start from 1-flip result
            B = [row[:] for row in A]
            B[fd1['u']][fd1['v']] = 0
            B[fd1['v']][fd1['u']] = 1

            for u2 in range(n):
                for v2 in range(n):
                    if u2 == v2 or B[u2][v2] == 0:
                        continue

                    C = [row[:] for row in B]
                    C[u2][v2] = 0
                    C[v2][u2] = 1

                    surplus, O2, O3, Z2 = compute_surplus_full(C, n)
                    all_s2[surplus] += 1
                    if surplus < min_s2:
                        min_s2 = surplus
                        print(f"    New min surplus={surplus}: flip1=({fd1['u']},{fd1['v']}), flip2=({u2},{v2})")

        print(f"  Min surplus after 2 flips: {min_s2}")
        print(f"  Surplus distribution (bottom 5):")
        for s in sorted(all_s2.keys())[:5]:
            print(f"    surplus={s}: {all_s2[s]}")

print("\nDone.")
