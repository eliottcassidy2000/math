#!/usr/bin/env python3
"""
beta2_filling_algebraic.py - Algebraic structure of DT+cancellation filling

For a tournament T, the 3-paths that fill Z_2 come from two sources:
1. DT (Diamond Trick): paths (a,b,c,d) with a->c AND b->d
   These are automatically in Omega_3.
2. Cancellation: DIFFERENCES of paths sharing a bad face
   If (a,b,c,d) has bad face (a,c) [because c->a], and there are
   k >= 2 paths with the same bad face, we get k-1 elements of Omega_3.

QUESTION 1: What determines whether DT alone suffices?
QUESTION 2: Can we prove the combined space always fills Z_2?
QUESTION 3: Is there a direct formula for rk(im d_3) in terms of
            combinatorial invariants?

APPROACH: For each tournament, decompose the Omega_3 space into:
- DT subspace
- Cancellation subspace
- Remainder (if any)
And track their contributions to rk(d_3).

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


def analyze_filling(A, n):
    """Decompose the Omega_3 filling of Z_2.

    Returns dict with:
    - dim_Z2: dimension of Z_2
    - rk_d3: rank of d_3 (should equal dim_Z2 if beta_2 = 0)
    - n_dt: number of DT paths
    - n_cancel: number of cancellation pairs
    - rk_dt_fill: rank of d_3 restricted to DT subspace
    - rk_cancel_fill: rank of d_3 restricted to DT+cancellation
    - dim_Om3: dimension of Omega_3
    - n_A3: number of allowed 3-paths
    - n_bad02: number of 3-paths with bad face at positions (0,2)
    - n_bad13: number of 3-paths with bad face at positions (1,3)
    """
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    d_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if d_om2 == 0:
        return {'dim_Z2': 0, 'rk_d3': 0, 'n_dt': 0, 'dim_Om3': 0}

    # Z_2
    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    om1 = compute_omega_basis(A, n, 1, a1, enumerate_allowed_paths(A, n, 0))
    coords2, _, _, _ = np.linalg.lstsq(om1, bd2_om, rcond=None)
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    dim_Z2 = d_om2 - rk2

    if dim_Z2 == 0 or not a3:
        return {'dim_Z2': dim_Z2, 'rk_d3': 0, 'n_dt': 0, 'dim_Om3': 0,
                'rk_dt_fill': 0, 'rk_cancel_fill': 0}

    # Omega_3
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    # rk(d_3)
    if d_om3 > 0:
        bd3 = build_full_boundary_matrix(a3, a2)
        bd3_om = bd3 @ om3
        coords3, _, _, _ = np.linalg.lstsq(om2, bd3_om, rcond=None)
        rk_d3 = np.linalg.matrix_rank(coords3, tol=1e-8)
    else:
        bd3 = build_full_boundary_matrix(a3, a2)
        rk_d3 = 0

    # DT classification
    dt_idx = []
    bad02_groups = defaultdict(list)  # (a,c,d) -> list of path indices
    bad13_groups = defaultdict(list)  # (a,b,d) -> list of path indices
    n_bad02 = 0
    n_bad13 = 0

    for i, p in enumerate(a3):
        a_, b_, c_, d_ = p
        is_02 = A[a_][c_]  # a->c
        is_13 = A[b_][d_]  # b->d

        if is_02 and is_13:
            dt_idx.append(i)
        else:
            if not is_02:
                bad02_groups[(a_, c_, d_)].append(i)
                n_bad02 += 1
            if not is_13:
                bad13_groups[(a_, b_, d_)].append(i)
                n_bad13 += 1

    # DT filling
    if dt_idx:
        V_dt = np.zeros((len(a3), len(dt_idx)))
        for j, idx in enumerate(dt_idx):
            V_dt[idx, j] = 1
        bd3_dt = bd3 @ V_dt
        coords_dt, _, _, _ = np.linalg.lstsq(om2, bd3_dt, rcond=None)
        rk_dt_fill = np.linalg.matrix_rank(coords_dt, tol=1e-8)
    else:
        rk_dt_fill = 0

    # Cancellation pairs
    cancel_vecs = []
    for group_key, indices in bad02_groups.items():
        if len(indices) >= 2:
            for j in range(1, len(indices)):
                v = np.zeros(len(a3))
                v[indices[0]] = 1
                v[indices[j]] = -1
                cancel_vecs.append(v)
    for group_key, indices in bad13_groups.items():
        if len(indices) >= 2:
            for j in range(1, len(indices)):
                v = np.zeros(len(a3))
                v[indices[0]] = 1
                v[indices[j]] = -1
                cancel_vecs.append(v)

    n_cancel = len(cancel_vecs)

    # DT + cancellation filling
    all_vecs = []
    for idx in dt_idx:
        v = np.zeros(len(a3))
        v[idx] = 1
        all_vecs.append(v)
    all_vecs.extend(cancel_vecs)

    if all_vecs:
        V_all = np.column_stack(all_vecs)
        bd3_all = bd3 @ V_all
        coords_all, _, _, _ = np.linalg.lstsq(om2, bd3_all, rcond=None)
        rk_cancel_fill = np.linalg.matrix_rank(coords_all, tol=1e-8)
    else:
        rk_cancel_fill = 0

    return {
        'dim_Z2': dim_Z2, 'rk_d3': rk_d3,
        'n_dt': len(dt_idx), 'dim_Om3': d_om3,
        'rk_dt_fill': rk_dt_fill, 'rk_cancel_fill': rk_cancel_fill,
        'n_cancel': n_cancel, 'n_A3': len(a3),
        'n_bad02': n_bad02, 'n_bad13': n_bad13,
    }


# ============================================================
# ANALYSIS 1: n=5 exhaustive
# ============================================================
print("=" * 70)
print("ANALYSIS 1: Filling decomposition at n=5")
print("=" * 70)

n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

dt_suffices = 0
cancel_needed = 0
cancel_suffices = 0
cancel_fails = 0
no_cycle = 0

filling_stats = Counter()

for bits in range(total):
    A = build_adj(n, bits)
    res = analyze_filling(A, n)

    if res['dim_Z2'] == 0:
        no_cycle += 1
        continue

    if res['rk_dt_fill'] >= res['dim_Z2']:
        dt_suffices += 1
    elif res['rk_cancel_fill'] >= res['dim_Z2']:
        cancel_needed += 1
        cancel_suffices += 1
    else:
        cancel_fails += 1

    filling_stats[(res['dim_Z2'], res['n_dt'], res['rk_dt_fill'],
                   res['n_cancel'], res['rk_cancel_fill'])] += 1

print(f"No Z_2 cycle: {no_cycle}")
print(f"DT alone fills Z_2: {dt_suffices}")
print(f"Cancel needed and works: {cancel_suffices}")
print(f"Cancel fails: {cancel_fails}")

print(f"\nFilling patterns (Z2, n_DT, rk_DT_fill, n_cancel, rk_cancel_fill):")
for key in sorted(filling_stats.keys()):
    gap = key[0] - key[2]  # Z2 - rk_DT
    print(f"  Z2={key[0]}, DT={key[1]}, rk_DT={key[2]}, cancel={key[3]}, "
          f"rk_DT+cancel={key[4]}, gap={gap}: {filling_stats[key]}")


# ============================================================
# ANALYSIS 2: n=6 exhaustive — DT deficit structure
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: Filling decomposition at n=6")
print("=" * 70)

n = 6
n_arcs = n*(n-1)//2
total = 1 << n_arcs

dt_suffices_6 = 0
cancel_needed_6 = 0
cancel_suffices_6 = 0
cancel_fails_6 = 0
no_cycle_6 = 0

deficit_dist = Counter()  # Z2 - rk_DT distribution

t0 = time.time()
for bits in range(total):
    if bits % 10000 == 0 and bits > 0:
        dt = time.time() - t0
        print(f"  ... {bits}/{total} ({dt:.0f}s)")
    A = build_adj(n, bits)
    res = analyze_filling(A, n)

    if res['dim_Z2'] == 0:
        no_cycle_6 += 1
        continue

    deficit = res['dim_Z2'] - res['rk_dt_fill']
    deficit_dist[deficit] += 1

    if deficit == 0:
        dt_suffices_6 += 1
    elif res['rk_cancel_fill'] >= res['dim_Z2']:
        cancel_needed_6 += 1
        cancel_suffices_6 += 1
    else:
        cancel_fails_6 += 1

dt = time.time() - t0
print(f"\nCompleted in {dt:.0f}s")
print(f"No Z_2 cycle: {no_cycle_6}")
print(f"DT alone fills Z_2: {dt_suffices_6}")
print(f"Cancel needed and works: {cancel_suffices_6}")
print(f"Cancel fails: {cancel_fails_6}")
print(f"DT deficit distribution: {dict(sorted(deficit_dist.items()))}")


# ============================================================
# ANALYSIS 3: WHY does DT + cancellation work?
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: Structure of cancellation pairs")
print("=" * 70)

# For n=5, focus on tournaments where DT alone fails
n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

print(f"\nn=5: Cases where DT alone fails — what do cancellation pairs look like?")
example_count = 0
for bits in range(total):
    A = build_adj(n, bits)
    res = analyze_filling(A, n)
    if res['dim_Z2'] == 0 or res['rk_dt_fill'] >= res['dim_Z2']:
        continue
    if example_count >= 3:
        continue
    example_count += 1

    a3 = enumerate_allowed_paths(A, n, 3)
    a2 = enumerate_allowed_paths(A, n, 2)

    scores = [sum(row) for row in A]
    c3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
             if max(A[i][j]+A[i][k], A[j][i]+A[j][k], A[k][i]+A[k][j]) < 2)

    print(f"\n  T#{bits}, scores={scores}, c_3={c3}")
    print(f"    |A_3|={len(a3)}, dim_Om3={res.get('dim_Om3', '?')}")
    print(f"    Z_2={res['dim_Z2']}, rk_DT={res['rk_dt_fill']}, "
          f"rk_DT+cancel={res['rk_cancel_fill']}")

    # Show the bad face groups
    bad02 = defaultdict(list)
    bad13 = defaultdict(list)
    for i, p in enumerate(a3):
        a_, b_, c_, d_ = p
        if not A[a_][c_]:
            bad02[(a_, c_, d_)].append(p)
        if not A[b_][d_]:
            bad13[(a_, b_, d_)].append(p)

    # Only show groups with >= 2 members (these give cancellation elements)
    n_useful_02 = sum(1 for g in bad02.values() if len(g) >= 2)
    n_useful_13 = sum(1 for g in bad13.values() if len(g) >= 2)
    print(f"    Bad face groups with cancel: {n_useful_02} (type 02) + {n_useful_13} (type 13)")

    for key, paths in sorted(bad02.items()):
        if len(paths) >= 2:
            print(f"      Face ({key[0]},{key[1]},{key[2]}): {len(paths)} paths = {paths[:4]}")
    for key, paths in sorted(bad13.items()):
        if len(paths) >= 2:
            print(f"      Face ({key[0]},{key[1]},{key[2]}): {len(paths)} paths = {paths[:4]}")


# ============================================================
# ANALYSIS 4: DT paths as function of 4-vertex sub-tournaments
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: DT path counts and structure")
print("=" * 70)

# A DT path (a,b,c,d) requires: a->b, b->c, c->d, a->c, b->d
# That's 5 arcs out of 6 for the 4-vertex set {a,b,c,d}.
# The 6th arc is between a and d (can go either way).
# So DT paths correspond to "almost-transitive" 4-vertex sub-tournaments.

# Count: for a transitive 4-tournament a>b>c>d (meaning a->b, a->c, a->d, b->c, b->d, c->d):
# DT paths: (a,b,c,d) has a->c YES, b->d YES. DT!
# (a,b,d,c)? needs b->d YES, a->d YES, d->c?? NO: c->d. Not allowed.
# (a,c,d,_)? etc.
# Actually, in a FULLY transitive 4-tournament, how many DT 3-paths are there?

print("\nDT paths in transitive 4-tournament {0>1>2>3}:")
A4 = [[0,1,1,1],[0,0,1,1],[0,0,0,1],[0,0,0,0]]
for a in range(4):
    for b in range(4):
        if b == a or not A4[a][b]:
            continue
        for c in range(4):
            if c == a or c == b or not A4[b][c]:
                continue
            for d in range(4):
                if d == a or d == b or d == c or not A4[c][d]:
                    continue
                is_dt = A4[a][c] and A4[b][d]
                if is_dt:
                    print(f"  DT: ({a},{b},{c},{d})")

# Count DT paths per 4-vertex set type
print("\n\nDT count by 4-vertex sub-tournament type:")
from itertools import combinations

# At n=5, each 4-vertex subset has a specific tournament type
# Count DT for each sub-tournament and relate to the score sequence

for bits in range(min(total, 50)):
    A = build_adj(5, bits)
    scores = [sum(row) for row in A]
    total_dt = 0
    for quad in combinations(range(5), 4):
        a, b, c, d = quad
        sub_scores = []
        for i in quad:
            s = sum(A[i][j] for j in quad if j != i)
            sub_scores.append(s)
        sub_scores_sorted = tuple(sorted(sub_scores))

        # Count DT paths in this 4-set
        dt_count = 0
        for perm in [(a,b,c,d), (a,b,d,c), (a,c,b,d), (a,c,d,b), (a,d,b,c), (a,d,c,b),
                     (b,a,c,d), (b,a,d,c), (b,c,a,d), (b,c,d,a), (b,d,a,c), (b,d,c,a),
                     (c,a,b,d), (c,a,d,b), (c,b,a,d), (c,b,d,a), (c,d,a,b), (c,d,b,a),
                     (d,a,b,c), (d,a,c,b), (d,b,a,c), (d,b,c,a), (d,c,a,b), (d,c,b,a)]:
            p0, p1, p2, p3 = perm
            if A[p0][p1] and A[p1][p2] and A[p2][p3] and A[p0][p2] and A[p1][p3]:
                dt_count += 1
        total_dt += dt_count
    break  # Just one example

print(f"  bits=0: total DT = {total_dt}")


# ============================================================
# ANALYSIS 5: Is rk(DT boundary) = |DT| - something?
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: DT rank vs count at n=5")
print("=" * 70)

dt_rank_data = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    d_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if d_om2 == 0 or not a3:
        continue

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

    dt_rank_data[(len(dt_idx), rk_dt)] += 1

print("(|DT|, rk(DT boundary)) distribution:")
for key in sorted(dt_rank_data.keys()):
    print(f"  |DT|={key[0]}, rk={key[1]}: {dt_rank_data[key]}")


# ============================================================
# ANALYSIS 6: Is there a formula for n_DT?
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 6: Formula for number of DT paths")
print("=" * 70)

# A DT path (a,b,c,d) requires: a->b, b->c, c->d, a->c, b->d
# This is: a->b->c->d with additional a->c and b->d
# The number of such paths relates to the number of "almost transitive" 4-sets

dt_by_score = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    a3 = enumerate_allowed_paths(A, n, 3)
    n_dt = sum(1 for p in a3 if A[p[0]][p[2]] and A[p[1]][p[3]])
    dt_by_score[scores] += n_dt

n_by_score = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    n_by_score[scores] += 1

print("Average |DT| by score sequence:")
for scores in sorted(n_by_score.keys()):
    avg = dt_by_score[scores] / n_by_score[scores]
    print(f"  {scores}: avg |DT| = {avg:.1f}")


print("\nDone.")
