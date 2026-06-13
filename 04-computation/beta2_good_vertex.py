#!/usr/bin/env python3
"""
beta2_good_vertex.py - Find a "good" vertex for LES induction

For the inductive proof of beta2=0:
  LES: 0 -> H_2(T) -> H_2(R(v)) -> H_1(T\v) -> H_1(T) -> ...
  H_2(T)=0 iff the connecting map is injective.

Question: For every tournament T, is there a vertex v such that
either beta1(T\v) = 0 OR the map H_1(T\v) -> H_1(T) is injective?

If so, we can always find a "good" vertex for the induction step.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix, path_betti_numbers
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


def check_h1_map_injective(A, n, v):
    """Check if H_1(T\\v) -> H_1(T) is injective."""
    verts = [i for i in range(n) if i != v]
    n_sub = n - 1
    A_sub = [[A[verts[i]][verts[j]] for j in range(n_sub)] for i in range(n_sub)]

    betti_sub = path_betti_numbers(A_sub, n_sub, max_dim=1)
    if betti_sub[1] == 0:
        return True  # Trivially injective

    # Compute Z_1(T\v) and B_1(T\v)
    a1_s = enumerate_allowed_paths(A_sub, n_sub, 1)
    a2_s = enumerate_allowed_paths(A_sub, n_sub, 2)
    om1_s = compute_omega_basis(A_sub, n_sub, 1, a1_s, [(i,) for i in range(n_sub)])
    om2_s = compute_omega_basis(A_sub, n_sub, 2, a2_s, a1_s)
    bd1_s = build_full_boundary_matrix(a1_s, [(i,) for i in range(n_sub)])
    bd2_s = build_full_boundary_matrix(a2_s, a1_s)

    d_om1_s = om1_s.shape[1] if om1_s.ndim == 2 else 0
    d_om2_s = om2_s.shape[1] if om2_s.ndim == 2 else 0

    if d_om1_s == 0:
        return True

    # Z_1 basis
    bd1_om_s = bd1_s @ om1_s
    U_s, S_s, Vt_s = np.linalg.svd(bd1_om_s, full_matrices=True)
    rk1 = sum(s > 1e-8 for s in S_s)
    z1_om_s = Vt_s[rk1:].T
    z1_a1_s = om1_s @ z1_om_s

    # B_1 basis
    if d_om2_s > 0:
        b1_a1_s = bd2_s @ om2_s
        rk_b1_s = np.linalg.matrix_rank(b1_a1_s, tol=1e-8)
    else:
        b1_a1_s = np.zeros((len(a1_s), 0))
        rk_b1_s = 0

    dim_h1_s = z1_a1_s.shape[1] - rk_b1_s
    if dim_h1_s == 0:
        return True

    # Map to T
    a1_t = enumerate_allowed_paths(A, n, 1)
    a1_t_idx = {path: i for i, path in enumerate(a1_t)}

    inc = np.zeros((len(a1_t), len(a1_s)))
    for j, path_s in enumerate(a1_s):
        path_t = tuple(verts[k] for k in path_s)
        if path_t in a1_t_idx:
            inc[a1_t_idx[path_t], j] = 1

    # B_1(T)
    a2_t = enumerate_allowed_paths(A, n, 2)
    om2_t = compute_omega_basis(A, n, 2, a2_t, a1_t)
    bd2_t = build_full_boundary_matrix(a2_t, a1_t)
    d_om2_t = om2_t.shape[1] if om2_t.ndim == 2 else 0

    if d_om2_t > 0:
        b1_t = bd2_t @ om2_t
    else:
        b1_t = np.zeros((len(a1_t), 0))

    # Find H_1(T\v) representatives not in B_1(T\v)
    # Then check if they land in B_1(T)
    for col in range(z1_a1_s.shape[1]):
        z = z1_a1_s[:, col]
        # Check if z is in B_1(T\v)
        if rk_b1_s > 0:
            res = np.linalg.lstsq(b1_a1_s, z, rcond=None)
            residual = np.linalg.norm(b1_a1_s @ res[0] - z)
            if residual < 1e-6:
                continue  # z is a boundary, skip

        # z is a non-boundary cycle; check if i(z) is in B_1(T)
        z_t = inc @ z
        if b1_t.shape[1] > 0:
            res = np.linalg.lstsq(b1_t, z_t, rcond=None)
            residual = np.linalg.norm(b1_t @ res[0] - z_t)
            if residual < 1e-6:
                return False  # i(z) is a boundary in T -> not injective

    return True


# ============================================================
# For EVERY tournament, does there exist a "good" vertex?
# ============================================================
print("=" * 70)
print("GOOD VERTEX EXISTENCE CHECK")
print("=" * 70)

for n in [5, 6]:
    total = 1 << (n*(n-1)//2)
    has_good_vertex = 0
    no_good_vertex = 0

    if n <= 5:
        rng = range(total)
    else:
        import random
        random.seed(42)
        rng = [random.randint(0, total-1) for _ in range(500)]

    for bits in rng:
        A = build_adj(n, bits)

        found_good = False
        for v in range(n):
            if check_h1_map_injective(A, n, v):
                found_good = True
                break

        if found_good:
            has_good_vertex += 1
        else:
            no_good_vertex += 1
            scores = tuple(sorted([sum(row) for row in A]))
            if no_good_vertex <= 5:
                print(f"  NO good vertex: bits={bits}, scores={scores}")

    cnt = len(list(rng)) if n > 5 else total
    print(f"\nn={n}: has good vertex in {has_good_vertex}/{cnt}, no good vertex: {no_good_vertex}")


# ============================================================
# Alternative: which vertices have beta1(T\v) = 0?
# ============================================================
print(f"\n{'='*70}")
print("VERTICES WITH BETA1(T\\v) = 0")
print("=" * 70)

n = 5
total = 1 << (n*(n-1)//2)

# For each tournament, which vertices have beta1(T\v) = 0?
for bits in [0, 341, 10]:
    A = build_adj(n, bits)
    scores_full = tuple(sorted([sum(row) for row in A]))

    for v in range(n):
        verts = [i for i in range(n) if i != v]
        A_sub = [[A[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]
        betti_sub = path_betti_numbers(A_sub, n-1, max_dim=1)
        scores_sub = tuple(sorted([sum(row) for row in A_sub]))
        du_v = sum(A[v])
        print(f"  T(bits={bits}), delete v={v} (d_v={du_v}): scores_sub={scores_sub}, beta1={betti_sub[1]}")


# ============================================================
# KEY: Is the "source" (max out-degree) vertex always good?
# Or the "sink" (min out-degree) vertex?
# ============================================================
print(f"\n{'='*70}")
print("SOURCE/SINK VERTEX ANALYSIS")
print("=" * 70)

n = 6
import random
random.seed(42)

source_good = 0
sink_good = 0
total_tests = 0

for _ in range(500):
    bits = random.randint(0, (1 << (n*(n-1)//2)) - 1)
    A = build_adj(n, bits)
    degrees = [sum(A[i]) for i in range(n)]

    # Source = max out-degree
    source = degrees.index(max(degrees))
    # Sink = min out-degree
    sink = degrees.index(min(degrees))

    total_tests += 1
    if check_h1_map_injective(A, n, source):
        source_good += 1
    if check_h1_map_injective(A, n, sink):
        sink_good += 1

print(f"n={n}, {total_tests} tests:")
print(f"  Source (max d) is good: {source_good}/{total_tests} ({100*source_good/total_tests:.1f}%)")
print(f"  Sink (min d) is good: {sink_good}/{total_tests} ({100*sink_good/total_tests:.1f}%)")


# ============================================================
# What about n=7?
# ============================================================
print(f"\n{'='*70}")
print("N=7 SPOT CHECK")
print("=" * 70)

n = 7
random.seed(42)

for _ in range(100):
    bits = random.randint(0, (1 << (n*(n-1)//2)) - 1)
    A = build_adj(n, bits)

    found_good = False
    for v in range(n):
        verts = [i for i in range(n) if i != v]
        A_sub = [[A[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]
        betti_sub = path_betti_numbers(A_sub, n-1, max_dim=1)
        if betti_sub[1] == 0:
            found_good = True
            break

    if not found_good:
        # Need to check injectivity
        for v in range(n):
            if check_h1_map_injective(A, n, v):
                found_good = True
                break

    if not found_good:
        scores = tuple(sorted([sum(row) for row in A]))
        print(f"  NO good vertex at n=7: bits={bits}, scores={scores}")

print("n=7 check complete (100 samples)")


print("\n\nDone.")
