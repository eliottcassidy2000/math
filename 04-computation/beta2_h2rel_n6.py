#!/usr/bin/env python3
"""
beta2_h2rel_n6.py — Verify H₂^rel = dim(ker i_*) at n=6

At n=5 we proved H₂(T,T\v) = max(0, β₁(T\v) - β₁(T)) for all 5120 pairs.
This script checks whether this holds at n=6 (32768 tournaments × 6 vertices).

The LES gives β₂(T) = H₂^rel - dim(ker i_*). If H₂^rel = dim(ker i_*),
then β₂(T) = 0.

Author: opus-2026-03-08-S45
"""
import sys
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved


def delete_vertex(A, n, v):
    B = []
    for i in range(n):
        if i == v: continue
        row = []
        for j in range(n):
            if j == v: continue
            row.append(A[i][j])
        B.append(row)
    return B


def local_to_global_map(n, v):
    return [i for i in range(n) if i != v]


def embed_Tv_paths_in_T(allowed_T, allowed_Tv, omega_Tv, n, v):
    if omega_Tv.ndim != 2 or omega_Tv.shape[1] == 0:
        return np.zeros((len(allowed_T), 0))
    loc2glob = local_to_global_map(n, v)
    T_list = [tuple(x) for x in allowed_T]
    Tv_list = [tuple(x) for x in allowed_Tv]
    T_idx = {p: i for i, p in enumerate(T_list)}
    incl = np.zeros((len(T_list), len(Tv_list)))
    for j, path_local in enumerate(Tv_list):
        path_global = tuple(loc2glob[k] for k in path_local)
        if path_global in T_idx:
            incl[T_idx[path_global], j] = 1
    return incl @ omega_Tv


def compute_h2_rel(A, n, v):
    B = delete_vertex(A, n, v)
    om_T = {}; ap_T = {}; om_Tv = {}; ap_Tv = {}

    for p in range(5):
        ap_T[p] = enumerate_allowed_paths(A, n, p)
        if p == 0: om_T[p] = np.eye(n)
        elif ap_T[p]: om_T[p] = compute_omega_basis(A, n, p, ap_T[p], ap_T[p-1])
        else: om_T[p] = np.zeros((0, 0))

    for p in range(4):
        ap_Tv[p] = enumerate_allowed_paths(B, n-1, p)
        if p == 0: om_Tv[p] = np.eye(n-1)
        elif ap_Tv[p]: om_Tv[p] = compute_omega_basis(B, n-1, p, ap_Tv[p], ap_Tv[p-1])
        else: om_Tv[p] = np.zeros((0, 0))

    def dim_om(om):
        return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

    d2T = dim_om(om_T[2])
    d3T = dim_om(om_T[3])

    if d2T == 0: return 0

    emb1 = embed_Tv_paths_in_T(ap_T[1], ap_Tv[1], om_Tv[1], n, v)
    emb2 = embed_Tv_paths_in_T(ap_T[2], ap_Tv[2], om_Tv[2], n, v)

    bd2_A = build_full_boundary_matrix(ap_T[2], ap_T[1])
    bd2_om = bd2_A @ om_T[2]

    if emb1.shape[1] > 0:
        M = np.hstack([bd2_om, -emb1])
    else:
        M = bd2_om

    U, S, Vt = np.linalg.svd(M, full_matrices=True)
    rk_M = int(sum(s > 1e-8 for s in S))
    null_space = Vt[rk_M:].T if rk_M < Vt.shape[0] else np.zeros((Vt.shape[1], 0))
    c_part = null_space[:d2T, :]
    preimage_A = om_T[2] @ c_part

    if emb2.shape[1] > 0 and preimage_A.shape[1] > 0:
        combined = np.hstack([emb2, preimage_A])
    elif preimage_A.shape[1] > 0:
        combined = preimage_A
    else:
        combined = emb2
    rk_combined = np.linalg.matrix_rank(combined, tol=1e-8) if combined.shape[1] > 0 else 0
    rk_emb2 = np.linalg.matrix_rank(emb2, tol=1e-8) if emb2.shape[1] > 0 else 0
    dim_ker_R2 = rk_combined - rk_emb2

    if d3T > 0:
        bd3_A = build_full_boundary_matrix(ap_T[3], ap_T[2])
        bd3_om = bd3_A @ om_T[3]
        if emb2.shape[1] > 0:
            combined_3 = np.hstack([emb2, bd3_om])
        else:
            combined_3 = bd3_om
        rk_combined_3 = np.linalg.matrix_rank(combined_3, tol=1e-8)
        dim_im_R3 = rk_combined_3 - rk_emb2
    else:
        dim_im_R3 = 0

    return dim_ker_R2 - dim_im_R3


def compute_betti_1(A, n):
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    if not ap1: return 0
    ap0 = enumerate_allowed_paths(A, n, 0)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    d1 = om1.shape[1] if om1.ndim == 2 else 0

    bd1 = build_full_boundary_matrix(ap1, ap0)
    bd1_om = bd1 @ om1
    S1 = np.linalg.svd(np.linalg.lstsq(np.eye(n), bd1_om, rcond=None)[0], compute_uv=False)
    rk1 = int(sum(s > 1e-8 for s in S1))

    if ap2:
        om2 = compute_omega_basis(A, n, 2, ap2, ap1)
        d2 = om2.shape[1] if om2.ndim == 2 else 0
        if d2 > 0:
            bd2 = build_full_boundary_matrix(ap2, ap1)
            bd2_om = bd2 @ om2
            S2 = np.linalg.svd(np.linalg.lstsq(om1, bd2_om, rcond=None)[0], compute_uv=False)
            rk2 = int(sum(s > 1e-8 for s in S2))
        else: rk2 = 0
    else: rk2 = 0

    return (d1 - rk1) - rk2


# ===== MAIN =====
import time
print("=" * 70)
print("H₂^rel CHECK AT n=6")
print("=" * 70)

n = 6
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(pairs)
total_T = 1 << m
print(f"n={n}, {total_T} tournaments, {total_T * n} (T,v) pairs")

h2_dist = Counter()
breakdown = Counter()
bugs = 0
total = 0
t0 = time.time()

for bits in range(total_T):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    b1_T = compute_betti_1(A, n)

    for v in range(n):
        B = delete_vertex(A, n, v)
        b1_Tv = compute_betti_1(B, n-1)
        h2r = compute_h2_rel(A, n, v)
        h2_dist[h2r] += 1
        total += 1

        les_min = max(0, b1_Tv - b1_T)
        if h2r < les_min:
            bugs += 1
            if bugs <= 5:
                print(f"  BUG: bits={bits}, v={v}, h2_rel={h2r}, b1_T={b1_T}, b1_Tv={b1_Tv}")

        expected = max(0, b1_Tv - b1_T)
        if h2r != expected:
            breakdown[('mismatch', b1_T, b1_Tv, h2r)] += 1
        else:
            breakdown[('match',)] += 1

    if bits % 1000 == 0 and bits > 0:
        elapsed = time.time() - t0
        rate = bits / elapsed
        eta = (total_T - bits) / rate
        print(f"  ... {bits}/{total_T} ({elapsed:.0f}s, ETA {eta:.0f}s), bugs={bugs}")

elapsed = time.time() - t0
print(f"\nCompleted in {elapsed:.1f}s")
print(f"\nH₂^rel distribution ({total} pairs):")
for val, count in sorted(h2_dist.items()):
    print(f"  H₂^rel = {val}: {count}")

print(f"\nLES violations: {bugs}")
if bugs == 0:
    print("ALL CONSISTENT WITH LES! ✓")

match_count = breakdown.get(('match',), 0)
mismatch_entries = {k: v for k, v in breakdown.items() if k[0] == 'mismatch'}
print(f"\nH₂^rel = max(0, β₁(T\\v) - β₁(T)):")
print(f"  Match: {match_count}/{total}")
if mismatch_entries:
    print(f"  Mismatches:")
    for k, cnt in sorted(mismatch_entries.items()):
        print(f"    b1_T={k[1]}, b1_Tv={k[2]}, h2_rel={k[3]}: {cnt}")

print("\nDone.")
