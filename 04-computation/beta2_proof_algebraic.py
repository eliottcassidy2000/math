#!/usr/bin/env python3
"""β_2 = 0 ALGEBRAIC PROOF EXPLORATION

Key insight from oriented graph analysis:
β_2 > 0 requires "parallel paths" through vertices with identical neighborhoods.
Tournaments prevent this because every pair has an edge.

APPROACH: For a 2-cycle z ∈ ker(∂_2) ∩ Ω_2, study how it decomposes
over vertex pairs {u,v}. The tournament edge u→v (or v→u) provides
"connecting" 3-paths that fill the cycle.

Specifically, for paths (a,b,c) and (a,b',c) in z where b≠b',
the edge between b and b' creates either:
- A DT path (b,b',c,?) or (?,a,b,b'), or
- A cancellation chain through b and b'.

TEST: For each 2-cycle, check if it can be decomposed into
"edge-specific" contributions, one per tournament edge.
"""
import numpy as np
from itertools import combinations
import sys, time
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix
)

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

# ===== APPROACH 1: Rank counting =====
# β_2 = 0 ⟺ rank(∂_2|Ω_2) + rank(∂_3|Ω_3) = dim(Ω_2)
# We know rank(∂_2|Ω_2) = dim(Ω_2) - dim(ker ∂_2|Ω_2)
# We know rank(∂_1|Ω_1) = n-1 (tournament is strongly connected? No, just connected)
# β_1 = dim(ker ∂_1|Ω_1) - rank(∂_2|Ω_2)
# So rank(∂_2|Ω_2) = dim(Ω_1) - (n-1) - β_1 = C(n,2) - n + 1 - β_1

# For β_2 = 0: rank(∂_3|Ω_3) = dim(Ω_2) - C(n,2) + n - 1 + β_1

print("=" * 70)
print("RANK COUNTING FOR β_2 = 0")
print("=" * 70)

for n in [4, 5, 6]:
    print(f"\n--- n={n} ---")
    t0 = time.time()

    data = defaultdict(list)
    count = 0
    for A in all_tournaments_gen(n):
        count += 1
        if n == 6 and count % 5000 == 0:
            print(f"  ... {count}/32768 ({time.time()-t0:.0f}s)", flush=True)

        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        a3 = enumerate_allowed_paths(A, n, 3)
        a1_list = [tuple(p) for p in a1]
        a2_list = [tuple(p) for p in a2]
        a3_list = [tuple(p) for p in a3]

        om2 = compute_omega_basis(A, n, 2, a2, a1)
        dim_om2 = om2.shape[1] if om2.ndim == 2 else 0

        om3 = compute_omega_basis(A, n, 3, a3, a2)
        dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

        # rank(∂_2)
        if dim_om2 > 0:
            bd2 = build_full_boundary_matrix(a2_list, a1_list)
            bd2_om = bd2 @ om2
            rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
        else:
            rank2 = 0

        # rank(∂_3)
        if dim_om3 > 0 and dim_om2 > 0:
            bd3 = build_full_boundary_matrix(a3_list, a2_list)
            im3 = bd3 @ om3
            im3_in_om2, _, _, _ = np.linalg.lstsq(om2, im3, rcond=None)
            rank3 = np.linalg.matrix_rank(im3_in_om2, tol=1e-8)
        else:
            rank3 = 0

        ker2 = dim_om2 - rank2
        beta1 = len(a1_list) - (n - 1) - rank2  # assuming β_0 = 1
        beta2 = ker2 - rank3

        t3 = sum(1 for a, b, c in combinations(range(n), 3)
                 if (A[a][b] and A[b][c] and A[c][a]) or
                    (A[b][a] and A[a][c] and A[c][b]))

        key = (t3, dim_om2, dim_om3, rank2, rank3, beta1, beta2)
        data[key] = data.get(key, 0) + 1

    t1 = time.time()
    print(f"  Done in {t1-t0:.1f}s")

    print(f"\n  (t3, dim_Ω2, dim_Ω3, rk_∂2, rk_∂3, β1, β2): count")
    for key in sorted(data.keys()):
        t3, d2, d3, r2, r3, b1, b2 = key
        # Only show if interesting or n is small
        if n == 6 and data[key] < 200:
            continue
        check = "✓" if b2 == 0 else "✗"
        # rank formula check: r2 should equal C(n,2) - n + 1 - b1
        pred_r2 = n*(n-1)//2 - n + 1 - b1
        r2_ok = "✓" if pred_r2 == r2 else "✗"
        # β_2 = 0 requires: r3 = d2 - r2 = d2 - C(n,2) + n - 1 + b1
        need_r3 = d2 - pred_r2
        r3_ok = "✓" if r3 == need_r3 else "✗"
        print(f"    t3={t3}, Ω2={d2}, Ω3={d3}, rk∂2={r2}({r2_ok}), rk∂3={r3}(need {need_r3},{r3_ok}), β1={b1}, β2={b2} {check}: {data[key]}")

# ===== APPROACH 2: Is there a formula for rank(∂_3)? =====
print(f"\n\n{'='*70}")
print("IS rank(∂_3) = dim(Ω_3) - dim(ker ∂_3)?")
print("What determines ker(∂_3)?")
print("="*70)

for n in [4, 5]:
    print(f"\n--- n={n} ---")
    ker3_data = defaultdict(list)

    for A in all_tournaments_gen(n):
        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        a3 = enumerate_allowed_paths(A, n, 3)
        a4 = enumerate_allowed_paths(A, n, 4) if n >= 5 else []
        a2_list = [tuple(p) for p in a2]
        a3_list = [tuple(p) for p in a3]
        a4_list = [tuple(p) for p in a4]

        om3 = compute_omega_basis(A, n, 3, a3, a2)
        dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

        if n >= 5 and dim_om3 > 0:
            om4 = compute_omega_basis(A, n, 4, a4, a3)
            dim_om4 = om4.shape[1] if om4.ndim == 2 else 0
        else:
            dim_om4 = 0

        om2 = compute_omega_basis(A, n, 2, a2, a1)
        dim_om2 = om2.shape[1] if om2.ndim == 2 else 0

        if dim_om3 > 0 and dim_om2 > 0:
            bd3 = build_full_boundary_matrix(a3_list, a2_list)
            im3 = bd3 @ om3
            im3_in_om2, _, _, _ = np.linalg.lstsq(om2, im3, rcond=None)
            rank3 = np.linalg.matrix_rank(im3_in_om2, tol=1e-8)
        else:
            rank3 = 0

        ker3 = dim_om3 - rank3

        # rank(∂_4)
        if dim_om4 > 0 and dim_om3 > 0:
            bd4 = build_full_boundary_matrix(a4_list, a3_list)
            im4 = bd4 @ om4
            im4_in_om3, _, _, _ = np.linalg.lstsq(om3, im4, rcond=None)
            rank4 = np.linalg.matrix_rank(im4_in_om3, tol=1e-8)
        else:
            rank4 = 0

        beta3 = ker3 - rank4

        t3 = sum(1 for a, b, c in combinations(range(n), 3)
                 if (A[a][b] and A[b][c] and A[c][a]) or
                    (A[b][a] and A[a][c] and A[c][b]))

        key = (t3, dim_om3, rank3, ker3, dim_om4, rank4, beta3)
        ker3_data[key] = ker3_data.get(key, 0) + 1

    print(f"  (t3, Ω3, rk∂3, ker∂3, Ω4, rk∂4, β3): count")
    for key in sorted(ker3_data.keys()):
        t3, d3, r3, k3, d4, r4, b3 = key
        print(f"    t3={t3}, Ω3={d3}, rk={r3}, ker={k3}, Ω4={d4}, rk∂4={r4}, β3={b3}: {ker3_data[key]}")

print("\nDone.")
