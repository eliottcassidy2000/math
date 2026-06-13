#!/usr/bin/env python3
"""
Investigate n=6 tournaments where DT boundaries DON'T span ker(∂_2|Ω_2).
These 960 cases need cancellation 3-chains from Ω_3 \ DT.
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

n = 6
print(f"{'='*70}")
print(f"INVESTIGATING DT-FAILURE CASES AT n={n}")
print(f"{'='*70}")

t0 = time.time()
count = 0
failures = []

for A in all_tournaments_gen(n):
    count += 1
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a2_list = [tuple(p) for p in a2]
    a3_list = [tuple(p) for p in a3]
    a2_set = set(a2_list)

    n_dt = 0
    dt_paths = []
    for p in a3_list:
        faces = [(p[1],p[2],p[3]), (p[0],p[2],p[3]), (p[0],p[1],p[3]), (p[0],p[1],p[2])]
        if all(f in a2_set for f in faces):
            n_dt += 1
            dt_paths.append(p)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        continue

    bd2 = build_full_boundary_matrix(a2_list, [tuple(p) for p in a1])
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker2 = dim_om2 - rank2

    if ker2 > 0 and n_dt > 0:
        a3_idx = {p: i for i, p in enumerate(a3_list)}
        bd3 = build_full_boundary_matrix(a3_list, a2_list)
        dt_bds = np.zeros((len(a2_list), len(dt_paths)))
        for j, p in enumerate(dt_paths):
            dt_bds[:, j] = bd3[:, a3_idx[p]]
        dt_bds_om2, _, _, _ = np.linalg.lstsq(om2, dt_bds, rcond=None)
        im_dt_rank = np.linalg.matrix_rank(dt_bds_om2, tol=1e-8)

        if im_dt_rank < ker2:
            om3 = compute_omega_basis(A, n, 3, a3, a2)
            dim_om3 = om3.shape[1] if om3.ndim == 2 else 0
            t3 = sum(1 for a, b, c in combinations(range(n), 3)
                     if (A[a][b] and A[b][c] and A[c][a]) or
                        (A[b][a] and A[a][c] and A[c][b]))
            scores = tuple(sorted([sum(A[i]) for i in range(n)]))

            failures.append({
                'idx': count,
                'ker': ker2,
                'im_dt': im_dt_rank,
                'gap': ker2 - im_dt_rank,
                'n_dt': n_dt,
                'dim_om2': dim_om2,
                'dim_om3': dim_om3,
                'n_a3': len(a3_list),
                't3': t3,
                'scores': scores,
                'A': [row[:] for row in A],
            })

    elif ker2 > 0 and n_dt == 0:
        om3 = compute_omega_basis(A, n, 3, a3, a2)
        dim_om3 = om3.shape[1] if om3.ndim == 2 else 0
        t3 = sum(1 for a, b, c in combinations(range(n), 3)
                 if (A[a][b] and A[b][c] and A[c][a]) or
                    (A[b][a] and A[a][c] and A[c][b]))
        scores = tuple(sorted([sum(A[i]) for i in range(n)]))
        failures.append({
            'idx': count, 'ker': ker2, 'im_dt': 0, 'gap': ker2,
            'n_dt': 0, 'dim_om2': dim_om2, 'dim_om3': dim_om3,
            'n_a3': len(a3_list), 't3': t3, 'scores': scores,
            'A': [row[:] for row in A],
        })

t1 = time.time()
print(f"\nFound {len(failures)} failures in {t1-t0:.1f}s")

# ===== Analyze failure characteristics =====
print(f"\n--- Failure characteristics ---")
gap_dist = Counter(f['gap'] for f in failures)
print(f"  Gap (ker - im_DT) distribution: {dict(gap_dist)}")

score_dist = Counter(f['scores'] for f in failures)
print(f"\n  Score sequences:")
for s in sorted(score_dist.keys()):
    print(f"    {s}: {score_dist[s]}")

t3_dist = Counter(f['t3'] for f in failures)
print(f"\n  t3 distribution:")
for t in sorted(t3_dist.keys()):
    print(f"    t3={t}: {t3_dist[t]}")

# Key dimensions
print(f"\n  (ker, im_DT, dim_Om3, n_dt, dim_Om2) patterns:")
pattern_dist = Counter((f['ker'], f['im_dt'], f['dim_om3'], f['n_dt'], f['dim_om2'])
                       for f in failures)
for key in sorted(pattern_dist.keys()):
    print(f"    ker={key[0]}, im_DT={key[1]}, Om3={key[2]}, |DT|={key[3]}, Om2={key[4]}: {pattern_dist[key]}")

# ===== For ONE failure case, show what the extra Ω_3 elements look like =====
print(f"\n\n{'='*70}")
print("DETAILED ANALYSIS OF ONE FAILURE")
print("="*70)

f = failures[0]
A = f['A']
print(f"Tournament #{f['idx']}: scores={f['scores']}, t3={f['t3']}")
for i in range(n):
    nbrs = [j for j in range(n) if A[i][j] == 1]
    print(f"  {i} -> {nbrs}")

a1 = enumerate_allowed_paths(A, n, 1)
a2 = enumerate_allowed_paths(A, n, 2)
a3 = enumerate_allowed_paths(A, n, 3)
a2_list = [tuple(p) for p in a2]
a3_list = [tuple(p) for p in a3]
a2_set = set(a2_list)

om3 = compute_omega_basis(A, n, 3, a3, a2)
dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

# Classify Ω_3 basis elements
dt_in_om3 = 0
cancel_in_om3 = 0
for j in range(dim_om3):
    col = om3[:, j]
    nz = np.abs(col) > 1e-8
    n_nonzero = np.sum(nz)
    if n_nonzero == 1:
        # Single path — should be DT
        idx_nz = np.where(nz)[0][0]
        p = a3_list[idx_nz]
        is_dt = all(f in a2_set for f in [(p[1],p[2],p[3]),
                    (p[0],p[2],p[3]), (p[0],p[1],p[3]), (p[0],p[1],p[2])])
        if is_dt:
            dt_in_om3 += 1
        else:
            print(f"  WARNING: single-path Ω_3 element NOT DT: {p}")
    else:
        cancel_in_om3 += 1
        if cancel_in_om3 <= 3:
            # Show what it looks like
            paths = [(a3_list[k], col[k]) for k in range(len(a3_list)) if abs(col[k]) > 1e-8]
            print(f"\n  Cancellation 3-chain #{cancel_in_om3} ({len(paths)} terms):")
            for path, coeff in paths:
                # Check which faces are not in A_2
                faces = [(path[1],path[2],path[3]), (path[0],path[2],path[3]),
                         (path[0],path[1],path[3]), (path[0],path[1],path[2])]
                bad = [f for f in faces if f not in a2_set]
                bad_str = f" [bad faces: {bad}]" if bad else ""
                print(f"    {coeff:+.4f} * {path}{bad_str}")

print(f"\n  Ω_3 breakdown: {dt_in_om3} DT (individual) + {cancel_in_om3} cancellation chains")
print(f"  ker(∂_2)={f['ker']}, im(∂_3|DT)={f['im_dt']}, dim(Ω_3)={dim_om3}")

# ===== Verify full β_2 = 0 for this case =====
om2 = compute_omega_basis(A, n, 2, a2, a1)
bd3 = build_full_boundary_matrix(a3_list, a2_list)
bd3_om3 = bd3 @ om3
bd3_in_om2, _, _, _ = np.linalg.lstsq(om2, bd3_om3, rcond=None)
full_rank = np.linalg.matrix_rank(bd3_in_om2, tol=1e-8)
print(f"  im(∂_3|Ω_3) rank = {full_rank} (need {f['ker']} for β_2=0)")
print(f"  β_2 = {f['ker'] - full_rank}")

print("\nDone.")
