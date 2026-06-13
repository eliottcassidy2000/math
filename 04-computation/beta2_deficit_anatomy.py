#!/usr/bin/env python3
"""
beta2_deficit_anatomy.py - Understand WHY cancellation covers the DT deficit

DT deficit occurs when the boundary vectors of DT paths don't span all of Z_2.
Cancellation pairs (differences of 3-paths sharing a bad face) provide extra
dimensions. WHY do they always suffice?

Key question: is there a STRUCTURAL reason that the deficit is bounded
and cancellation always covers it?

APPROACH:
1. At n=6, find all tournaments with DT deficit > 0 (960 of them)
2. Analyze what Z_2 vectors are NOT in im(DT boundary)
3. Show that cancellation pairs target EXACTLY the missing vectors
4. Find the common structure

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


print("=" * 70)
print("DEFICIT ANATOMY AT n=6")
print("=" * 70)

n = 6
n_arcs = n*(n-1)//2
total = 1 << n_arcs

# Step 1: Find all deficit tournaments
deficit_tours = []

t0 = time.time()
for bits in range(total):
    if bits % 10000 == 0 and bits > 0:
        dt = time.time() - t0
        print(f"  ... {bits}/{total} ({dt:.0f}s), found {len(deficit_tours)}")

    A = build_adj(n, bits)
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    d_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if d_om2 == 0:
        continue

    om1 = compute_omega_basis(A, n, 1, a1, enumerate_allowed_paths(A, n, 0))
    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    coords2, _, _, _ = np.linalg.lstsq(om1, bd2_om, rcond=None)
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    dim_Z2 = d_om2 - rk2
    if dim_Z2 == 0:
        continue

    # DT
    bd3 = build_full_boundary_matrix(a3, a2)
    dt_idx = [i for i, p in enumerate(a3) if A[p[0]][p[2]] and A[p[1]][p[3]]]
    if not dt_idx:
        continue

    V_dt = np.zeros((len(a3), len(dt_idx)))
    for j, idx in enumerate(dt_idx):
        V_dt[idx, j] = 1
    bd3_dt = bd3 @ V_dt
    coords_dt, _, _, _ = np.linalg.lstsq(om2, bd3_dt, rcond=None)
    rk_dt = np.linalg.matrix_rank(coords_dt, tol=1e-8)

    deficit = dim_Z2 - rk_dt
    if deficit > 0:
        scores = tuple(sorted([sum(row) for row in A]))
        c3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                 if max(A[i][j]+A[i][k], A[j][i]+A[j][k], A[k][i]+A[k][j]) < 2)
        deficit_tours.append({
            'bits': bits, 'scores': scores, 'c3': c3,
            'deficit': deficit, 'dim_Z2': dim_Z2,
            'n_dt': len(dt_idx), 'rk_dt': rk_dt,
        })

dt = time.time() - t0
print(f"\nDone scanning in {dt:.0f}s")
print(f"Deficit tournaments: {len(deficit_tours)}")

# Step 2: Characterize deficit tournaments
score_dist = Counter()
c3_dist = Counter()
deficit_dist = Counter()

for t in deficit_tours:
    score_dist[t['scores']] += 1
    c3_dist[t['c3']] += 1
    deficit_dist[t['deficit']] += 1

print(f"\nDeficit distribution: {dict(sorted(deficit_dist.items()))}")
print(f"\nScore distribution of deficit tournaments:")
for scores in sorted(score_dist.keys()):
    print(f"  {scores}: {score_dist[scores]}")
print(f"\nc_3 distribution: {dict(sorted(c3_dist.items()))}")


# Step 3: For a few deficit tournaments, analyze what Z_2 vector is missing
print(f"\n{'='*70}")
print("DETAILED ANALYSIS: What Z_2 vectors are not covered by DT?")
print("=" * 70)

for t in deficit_tours[:5]:
    bits = t['bits']
    A = build_adj(n, bits)

    a0 = enumerate_allowed_paths(A, n, 0)
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om1 = compute_omega_basis(A, n, 1, a1, a0)
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    d_om2 = om2.shape[1]

    # Z_2 basis
    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    coords2, _, _, _ = np.linalg.lstsq(om1, bd2_om, rcond=None)
    U2, S2, Vt2 = np.linalg.svd(coords2, full_matrices=True)
    rk2 = int(np.sum(np.abs(S2) > 1e-8))
    Z2_basis = Vt2[rk2:].T  # d_om2 x dim_Z2

    # DT image
    bd3 = build_full_boundary_matrix(a3, a2)
    dt_idx = [i for i, p in enumerate(a3) if A[p[0]][p[2]] and A[p[1]][p[3]]]
    V_dt = np.zeros((len(a3), len(dt_idx)))
    for j, idx in enumerate(dt_idx):
        V_dt[idx, j] = 1
    bd3_dt = bd3 @ V_dt
    coords_dt, _, _, _ = np.linalg.lstsq(om2, bd3_dt, rcond=None)
    rk_dt = np.linalg.matrix_rank(coords_dt, tol=1e-8)

    # Find Z_2 vectors NOT in im(DT)
    # Project Z_2 basis onto complement of im(DT)
    U_dt, S_dt, Vt_dt = np.linalg.svd(coords_dt, full_matrices=True)
    cols_dt = U_dt[:, :rk_dt]  # basis for im(DT) in Omega_2 coords
    proj_dt = cols_dt @ cols_dt.T  # projection onto im(DT)
    complement = np.eye(d_om2) - proj_dt  # projection onto complement

    # Project Z_2 basis onto complement
    Z2_comp = complement @ Z2_basis
    U_comp, S_comp, _ = np.linalg.svd(Z2_comp, full_matrices=False)
    missing_dim = int(np.sum(np.abs(S_comp) > 1e-8))

    print(f"\nT#{bits}, scores={t['scores']}, c3={t['c3']}")
    print(f"  Z_2 dim={t['dim_Z2']}, DT rk={t['rk_dt']}, deficit={t['deficit']}")
    print(f"  Missing from DT: {missing_dim} dimensions")

    # Show the missing Z_2 vector(s)
    for j in range(min(missing_dim, 2)):
        if S_comp[j] < 1e-8:
            continue
        missing_vec = U_comp[:, j] * S_comp[j]
        # Convert to path coords
        z_paths = om2 @ missing_vec
        nonzero = [(i, z_paths[i]) for i in range(len(z_paths)) if abs(z_paths[i]) > 1e-8]
        print(f"  Missing vector {j} ({len(nonzero)} nonzero paths):")
        for idx, coeff in sorted(nonzero, key=lambda x: abs(x[1]), reverse=True)[:6]:
            a_, b_, c_ = a2[idx]
            tt = "TT" if A[a_][c_] else "NT"
            print(f"    {coeff:+.4f} * ({a_},{b_},{c_}) [{tt}]")

    # Now check: which cancellation pair fills it?
    bad02 = defaultdict(list)
    bad13 = defaultdict(list)
    for i, p in enumerate(a3):
        a_, b_, c_, d_ = p
        if not A[a_][c_]:
            bad02[(a_, c_, d_)].append(i)
        if not A[b_][d_]:
            bad13[(a_, b_, d_)].append(i)

    cancel_vecs = []
    cancel_info = []
    for group_key, indices in bad02.items():
        if len(indices) >= 2:
            for j in range(1, len(indices)):
                v = np.zeros(len(a3))
                v[indices[0]] = 1
                v[indices[j]] = -1
                cancel_vecs.append(v)
                cancel_info.append(('02', group_key, a3[indices[0]], a3[indices[j]]))
    for group_key, indices in bad13.items():
        if len(indices) >= 2:
            for j in range(1, len(indices)):
                v = np.zeros(len(a3))
                v[indices[0]] = 1
                v[indices[j]] = -1
                cancel_vecs.append(v)
                cancel_info.append(('13', group_key, a3[indices[0]], a3[indices[j]]))

    # Check which cancellation pairs fill the missing direction
    for k, (cv, ci) in enumerate(zip(cancel_vecs, cancel_info)):
        bd_cv = bd3 @ cv
        cv_om, _, _, _ = np.linalg.lstsq(om2, bd_cv.reshape(-1, 1), rcond=None)
        cv_comp = complement @ cv_om.flatten()
        if np.linalg.norm(cv_comp) > 1e-8:
            print(f"  FILLING cancel pair: type={ci[0]}, bad_face={ci[1]}")
            print(f"    path1={ci[2]}, path2={ci[3]}")
            print(f"    boundary projection norm in complement: {np.linalg.norm(cv_comp):.6f}")
            break


# Step 4: Is the deficit related to specific vertex/arc structure?
print(f"\n{'='*70}")
print("STRUCTURAL ANALYSIS OF DEFICIT")
print("=" * 70)

# Check: is deficit correlated with beta_1?
print("\nDeficit vs beta_1:")
b1_data = Counter()
for t in deficit_tours[:200]:
    bits = t['bits']
    A = build_adj(n, bits)
    a0 = enumerate_allowed_paths(A, n, 0)
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    om1 = compute_omega_basis(A, n, 1, a1, a0)
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    bd1 = build_full_boundary_matrix(a1, a0)
    rk1 = np.linalg.matrix_rank(bd1 @ om1, tol=1e-8)
    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    coords2, _, _, _ = np.linalg.lstsq(om1, bd2_om, rcond=None)
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z1 = om1.shape[1] - rk1
    b1 = z1 - rk2
    b1_data[(t['deficit'], b1)] += 1

for key in sorted(b1_data.keys()):
    print(f"  deficit={key[0]}, beta_1={key[1]}: {b1_data[key]}")


# Step 5: What TYPES of 3-cycles create the deficit?
print("\nTypes of 3-cycles in deficit tournaments:")
# The score sequence (1,2,2,3,3,4) has 960 tournaments and ALL have deficit
# (720 with deficit=1, 240 with deficit=2)
# What's special about this score sequence?

# Count: how many 3-cycles, and what's the "3-cycle structure"?
print(f"\nScore (1,2,2,3,3,4): {score_dist.get((1,2,2,3,3,4), 0)} deficit tournaments")
print(f"Score (2,2,2,3,3,3): {score_dist.get((2,2,2,3,3,3), 0)} deficit tournaments")

# Check which scores have NO deficit
all_scores_n6 = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    all_scores_n6[scores] += 1

print(f"\nAll score sequences at n=6 (total tournaments):")
for scores in sorted(all_scores_n6.keys()):
    total_count = all_scores_n6[scores]
    deficit_count = score_dist.get(scores, 0)
    pct = 100 * deficit_count / total_count
    if deficit_count > 0:
        print(f"  {scores}: {deficit_count}/{total_count} ({pct:.1f}%) have deficit")
    else:
        print(f"  {scores}: 0/{total_count} (0%) have deficit")


# Step 6: What is Omega_3 doing in deficit cases?
print(f"\n{'='*70}")
print("OMEGA_3 STRUCTURE IN DEFICIT CASES")
print("=" * 70)

# Compare dim_Om3 for deficit vs non-deficit
om3_data = {'deficit': Counter(), 'no_deficit': Counter()}

for bits in range(min(total, 5000)):
    A = build_adj(n, bits)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a1 = enumerate_allowed_paths(A, n, 1)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om2 = om2.shape[1] if om2.ndim == 2 else 0
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    # Quick check if deficit
    scores = tuple(sorted([sum(row) for row in A]))
    has_deficit = scores in [(1,2,2,3,3,4), (2,2,2,3,3,3)]

    key = (d_om2, d_om3, len(a3))
    if has_deficit:
        om3_data['deficit'][key] += 1
    else:
        om3_data['no_deficit'][key] += 1

print("Deficit cases:")
for key in sorted(om3_data['deficit'].keys()):
    print(f"  Om2={key[0]}, Om3={key[1]}, |A_3|={key[2]}: {om3_data['deficit'][key]}")

print("\nNon-deficit (first 5):")
for key in sorted(om3_data['no_deficit'].keys())[:5]:
    print(f"  Om2={key[0]}, Om3={key[1]}, |A_3|={key[2]}: {om3_data['no_deficit'][key]}")


print("\nDone.")
