#!/usr/bin/env python3
"""
β_2 = 0 analysis at n=6 (all 32768 tournaments).
Verify: ker(∂_2) ≤ |DT| and im(∂_3|DT) = ker(∂_2|Ω_2).
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
print(f"β_2 = 0 ANALYSIS AT n={n} (ALL {2**((n*(n-1))//2)} TOURNAMENTS)")
print(f"{'='*70}")

t0 = time.time()
count = 0
ker_dt_dist = Counter()
gap_dist = Counter()
gap_formula_match = 0
dt_eq_indiv = 0
beta2_nonzero = 0

for A in all_tournaments_gen(n):
    count += 1
    if count % 5000 == 0:
        print(f"  ... {count}/32768 ({time.time()-t0:.0f}s)", flush=True)

    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a2_list = [tuple(p) for p in a2]
    a3_list = [tuple(p) for p in a3]
    a2_set = set(a2_list)

    # TT count
    n_tt = sum(1 for p in a2_list if A[p[0]][p[2]] == 1)

    # DT count (= paths with all faces in A_2)
    n_dt = 0
    dt_paths = []
    for p in a3_list:
        faces = [(p[1],p[2],p[3]), (p[0],p[2],p[3]), (p[0],p[1],p[3]), (p[0],p[1],p[2])]
        if all(f in a2_set for f in faces):
            n_dt += 1
            dt_paths.append(p)

    # Gap formula: Σ(mult-1)
    face_mult = defaultdict(int)
    for p in a2_list:
        a, b, c = p
        if A[c][a] == 1:
            face_mult[(a, c)] += 1
    predicted_gap = sum(max(0, m - 1) for m in face_mult.values())

    # Ω_2 and Ω_3
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0

    actual_gap = dim_om2 - n_tt
    if predicted_gap == actual_gap:
        gap_formula_match += 1

    if dim_om2 == 0:
        ker_dt_dist[(0, n_dt)] += 1
        continue

    # ker(∂_2|Ω_2)
    bd2 = build_full_boundary_matrix(a2_list, [tuple(p) for p in a1])
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker2 = dim_om2 - rank2

    ker_dt_dist[(ker2, n_dt)] += 1

    # Verify β_2 = 0 via DT boundaries
    if ker2 > 0 and n_dt > 0:
        a3_idx = {p: i for i, p in enumerate(a3_list)}
        bd3 = build_full_boundary_matrix(a3_list, a2_list)

        dt_bds = np.zeros((len(a2_list), len(dt_paths)))
        for j, p in enumerate(dt_paths):
            dt_bds[:, j] = bd3[:, a3_idx[p]]

        # Project to Ω_2
        dt_bds_om2, _, _, _ = np.linalg.lstsq(om2, dt_bds, rcond=None)
        im_dt_rank = np.linalg.matrix_rank(dt_bds_om2, tol=1e-8)

        if im_dt_rank < ker2:
            beta2_nonzero += 1
            print(f"  FAILURE: T#{count}, ker={ker2}, im_DT={im_dt_rank}")

    elif ker2 > 0 and n_dt == 0:
        beta2_nonzero += 1
        print(f"  FAILURE: T#{count}, ker={ker2}, no DT paths!")

t1 = time.time()
print(f"\nDone in {t1-t0:.1f}s")
print(f"\nRESULTS:")
print(f"  Gap formula matches: {gap_formula_match}/{count}")
print(f"  β_2 ≠ 0: {beta2_nonzero}/{count}")
print(f"\n  (ker, |DT|): count")
for key in sorted(ker_dt_dist.keys()):
    ker2, ndt = key
    print(f"    ker={ker2}, |DT|={ndt}: {ker_dt_dist[key]}")
always_le = all(k <= d for (k, d) in ker_dt_dist.keys())
print(f"\n  ker ≤ |DT| always? {always_le}")

print(f"\nAll done.")
