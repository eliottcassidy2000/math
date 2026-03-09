#!/usr/bin/env python3
"""beta2_cone_debug.py - Debug why d3_err is nonzero at n=7 Paley

The B*alpha=z system solves perfectly, but the reconstructed 3-chain
may contain 3-paths that are NOT allowed in A_3 (not in the path complex).

Check: for each T' 2-path (a,b,c), are (v,a,b,c) and (a,b,c,v) both
in A_3 (i.e., allowed 3-paths)?

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, random
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved


def paley_tournament(p):
    qr = set()
    for a in range(1, p):
        qr.add((a * a) % p)
    A = [[0] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                A[i][j] = 1
    return A


def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


def is_allowed(path, A):
    """Check if a path is allowed (all consecutive pairs are arcs)."""
    for i in range(len(path) - 1):
        if A[path[i]][path[i + 1]] != 1:
            return False
    return True


print("=" * 70)
print("DEBUG: Cone path allowedness at n=7 Paley")
print("=" * 70)

n = 7
A = paley_tournament(n)
paths3 = enumerate_allowed_paths(A, n, 3)
path3_set = set(tuple(p) for p in paths3)
path3_idx = {tuple(p): i for i, p in enumerate(paths3)}

v = 0
others = [x for x in range(n) if x != v]
B_sub, vlist = get_induced(A, n, others)
n_prime = len(others)
paths2_Tp = enumerate_allowed_paths(B_sub, n_prime, 2)
Tp_paths_orig = [tuple(vlist[x] for x in p) for p in paths2_Tp]

print(f"\nv={v}, T' vertices = {vlist}")
print(f"T' has {len(Tp_paths_orig)} allowed 2-paths")

front_ok = 0
back_ok = 0
both_ok = 0
neither = 0
front_only = 0
back_only = 0

for (a, b, c) in Tp_paths_orig:
    front = (v, a, b, c)
    back = (a, b, c, v)

    f_ok = front in path3_set and is_allowed(front, A)
    b_ok = back in path3_set and is_allowed(back, A)

    # More detailed: check arc existence
    f_arc = A[v][a] == 1  # v -> a required for front
    b_arc = A[c][v] == 1  # c -> v required for back

    if f_ok:
        front_ok += 1
    if b_ok:
        back_ok += 1
    if f_ok and b_ok:
        both_ok += 1
    if f_ok and not b_ok:
        front_only += 1
    if b_ok and not f_ok:
        back_only += 1
    if not f_ok and not b_ok:
        neither += 1

print(f"\nFront (v,a,b,c) allowed: {front_ok}/{len(Tp_paths_orig)}")
print(f"Back (a,b,c,v) allowed: {back_ok}/{len(Tp_paths_orig)}")
print(f"Both allowed: {both_ok}/{len(Tp_paths_orig)}")
print(f"Front only: {front_only}")
print(f"Back only: {back_only}")
print(f"Neither: {neither}")

# Show examples where paths fail
print("\nExamples of T' 2-paths where cone fails:")
count = 0
for (a, b, c) in Tp_paths_orig:
    front = (v, a, b, c)
    back = (a, b, c, v)
    f_ok = is_allowed(front, A)
    b_ok = is_allowed(back, A)

    if not f_ok or not b_ok:
        count += 1
        if count <= 5:
            f_arc = A[v][a]
            b_arc = A[c][v]
            print(f"  ({a},{b},{c}): v->{a}={f_arc}, {c}->v={b_arc}, "
                  f"front_ok={f_ok}, back_ok={b_ok}")

print(f"  Total failing: {count}")


# ============================================================
# REVISED APPROACH: Only use (a,b,c) where BOTH are allowed
# ============================================================
print(f"\n{'=' * 70}")
print("REVISED: Use only T' paths where BOTH cone paths are allowed")
print("=" * 70)

paths2 = enumerate_allowed_paths(A, n, 2)
path2_idx = {tuple(p): i for i, p in enumerate(paths2)}

# Filter T' paths
valid_Tp = [(a, b, c) for (a, b, c) in Tp_paths_orig
            if is_allowed((v, a, b, c), A) and is_allowed((a, b, c, v), A)]

print(f"Valid T' paths (both cone allowed): {len(valid_Tp)}/{len(Tp_paths_orig)}")

# Build filtered B matrix
B_fill = np.zeros((len(paths2), len(valid_Tp)))
for j, (a, b, c) in enumerate(valid_Tp):
    terms = [
        ((v, b, c), -1), ((v, a, c), +1), ((v, a, b), -1),
        ((b, c, v), +1), ((a, c, v), -1), ((a, b, v), +1),
    ]
    for path, coeff in terms:
        if path in path2_idx:
            B_fill[path2_idx[path], j] += coeff

rank_B = np.linalg.matrix_rank(B_fill, tol=1e-8)
print(f"rank(B_filtered) = {rank_B}")

# Test swap cycles with filtered B
P = [a for a in range(n) if a != v and A[v][a] == 1]
Q = [b for b in range(n) if b != v and A[b][v] == 1]
arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]

m = len(arcs_PQ)
rows = []
for a in P:
    row = [0] * m
    for j, (a2, b2) in enumerate(arcs_PQ):
        if a2 == a:
            row[j] = 1
    if any(r != 0 for r in row):
        rows.append(row)
for b in Q:
    row = [0] * m
    for j, (a2, b2) in enumerate(arcs_PQ):
        if b2 == b:
            row[j] = 1
    if any(r != 0 for r in row):
        rows.append(row)

C_mat = np.array(rows, dtype=float)
Sc = np.linalg.svd(C_mat, compute_uv=False)
rank_C = sum(s > 1e-8 for s in Sc)
ker_dim = m - rank_C
print(f"Swap cycle dim = {ker_dim}")

if ker_dim > 0:
    _, _, Vt = np.linalg.svd(C_mat, full_matrices=True)
    ker_basis = Vt[rank_C:]

    D3 = build_full_boundary_matrix(
        [tuple(p) for p in paths3], [tuple(p) for p in paths2]
    )

    for ki in range(min(ker_basis.shape[0], 3)):
        M_vec = ker_basis[ki]
        z = np.zeros(len(paths2))
        for j, (a, b) in enumerate(arcs_PQ):
            coeff = M_vec[j]
            if abs(coeff) < 1e-12:
                continue
            if (a, b, v) in path2_idx and (v, a, b) in path2_idx:
                z[path2_idx[(a, b, v)]] += coeff
                z[path2_idx[(v, a, b)]] -= coeff

        if np.max(np.abs(z)) < 1e-12:
            continue

        # Solve filtered system
        alpha, _, _, _ = np.linalg.lstsq(B_fill, z, rcond=None)
        err = np.max(np.abs(B_fill @ alpha - z))

        # Build actual 3-chain and verify
        w_full = np.zeros(len(paths3))
        for j, (a, b, c) in enumerate(valid_Tp):
            if abs(alpha[j]) < 1e-12:
                continue
            front = (v, a, b, c)
            back = (a, b, c, v)
            if front in path3_idx:
                w_full[path3_idx[front]] += alpha[j]
            if back in path3_idx:
                w_full[path3_idx[back]] += alpha[j]

        d3w = D3 @ w_full
        d3_err = np.max(np.abs(d3w - z))

        print(f"\n  Swap cycle {ki}: B_err={err:.2e}, d3_err={d3_err:.2e}")

        # Check omega membership
        omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
        if omega3.ndim == 2:
            proj, _, _, _ = np.linalg.lstsq(omega3, w_full, rcond=None)
            om_err = np.max(np.abs(w_full - omega3 @ proj))
            print(f"  w in Omega_3? err={om_err:.2e}")

            # Try solving in Omega_3 directly
            D3_om = D3 @ omega3
            c_om, _, _, _ = np.linalg.lstsq(D3_om, z, rcond=None)
            om_fill_err = np.max(np.abs(D3_om @ c_om - z))
            print(f"  Omega_3 direct fill: err={om_fill_err:.2e}")


# ============================================================
# PART 3: Does B*alpha=z work because z is ALWAYS in col(B)?
# ============================================================
print(f"\n{'=' * 70}")
print("PART 3: Full B (unfiltered) column space at n=5")
print("=" * 70)
print("Check: do the B columns span all swap cycles AND are in Omega_2?")

n = 5
for tidx, A in enumerate(all_tournaments_gen(n)):
    if tidx > 50:
        break

    paths1 = enumerate_allowed_paths(A, n, 1)
    paths2 = enumerate_allowed_paths(A, n, 2)
    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    path2_idx = {tuple(p): i for i, p in enumerate(paths2)}

    for v in range(n):
        P = [a for a in range(n) if a != v and A[v][a] == 1]
        Q = [b for b in range(n) if b != v and A[b][v] == 1]
        arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]
        if len(arcs_PQ) < 2:
            continue

        m = len(arcs_PQ)
        rows = []
        for a in P:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if a2 == a:
                    row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)
        for b in Q:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if b2 == b:
                    row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)

        if not rows:
            continue
        C_mat = np.array(rows, dtype=float)
        Sc = np.linalg.svd(C_mat, compute_uv=False)
        rank_C = sum(s > 1e-8 for s in Sc)
        ker_dim = m - rank_C
        if ker_dim == 0:
            continue

        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
        Tp_paths_orig = [tuple(vlist[x] for x in p) for p in paths2_Tp]

        B_fill = np.zeros((len(paths2), len(Tp_paths_orig)))
        for j, (a, b, c) in enumerate(Tp_paths_orig):
            terms = [
                ((v, b, c), -1), ((v, a, c), +1), ((v, a, b), -1),
                ((b, c, v), +1), ((a, c, v), -1), ((a, b, v), +1),
            ]
            for path, coeff in terms:
                if path in path2_idx:
                    B_fill[path2_idx[path], j] += coeff

        # Check each B column is in Omega_2
        dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
        if dim_O2 > 0:
            cols_in_omega = 0
            for j in range(B_fill.shape[1]):
                col = B_fill[:, j]
                if np.max(np.abs(col)) < 1e-12:
                    cols_in_omega += 1
                    continue
                proj, _, _, _ = np.linalg.lstsq(omega2, col, rcond=None)
                err = np.max(np.abs(col - omega2 @ proj))
                if err < 1e-8:
                    cols_in_omega += 1

            print(f"  T#{tidx} v={v}: B cols in Omega_2: "
                  f"{cols_in_omega}/{B_fill.shape[1]}")


print("\n\nDone.")
