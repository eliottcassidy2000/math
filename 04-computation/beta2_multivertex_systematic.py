#!/usr/bin/env python3
"""beta2_multivertex_systematic.py - Systematic multi-vertex cone test

The single-vertex filtered cone fails at n=8 (0.2%).
The multi-vertex cone (combining ALL vertices) appears to always work.

Test the multi-vertex cone systematically at n=7,8,9 and analyze
the rank structure.

ALSO: test an intermediate approach - for EACH swap cycle at vertex v,
does there exist SOME vertex u (possibly != v) whose single-vertex
cone fills it?

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time, random
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved

random.seed(42)


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


def compute_swap_cycles(A, n, v, paths2, path2_idx):
    """Compute swap cycle space for vertex v. Returns list of z vectors."""
    P = [a for a in range(n) if a != v and A[v][a] == 1]
    Q = [b for b in range(n) if b != v and A[b][v] == 1]
    arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]
    if len(arcs_PQ) < 2:
        return []

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
        return []

    C_mat = np.array(rows, dtype=float)
    Sc = np.linalg.svd(C_mat, compute_uv=False)
    rank_C = sum(s > 1e-8 for s in Sc)
    ker_dim = m - rank_C
    if ker_dim == 0:
        return []

    _, _, Vt = np.linalg.svd(C_mat, full_matrices=True)
    ker_basis = Vt[rank_C:]

    swap_vecs = []
    for ki in range(ker_basis.shape[0]):
        M_vec = ker_basis[ki]
        z = np.zeros(len(paths2))
        for j, (a, b) in enumerate(arcs_PQ):
            coeff = M_vec[j]
            if abs(coeff) < 1e-12:
                continue
            if (a, b, v) in path2_idx and (v, a, b) in path2_idx:
                z[path2_idx[(a, b, v)]] += coeff
                z[path2_idx[(v, a, b)]] -= coeff
        if np.max(np.abs(z)) > 1e-12:
            swap_vecs.append(z)

    return swap_vecs


def build_cone_B(A, n, v, paths2, path2_idx, filtered=True):
    """Build cone B matrix for vertex v.
    If filtered, only use T' paths where both cone paths are allowed."""
    others = [x for x in range(n) if x != v]
    B_sub, vlist = get_induced(A, n, others)
    paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
    Tp_all = [tuple(vlist[x] for x in p) for p in paths2_Tp]

    if filtered:
        valid = [(a, b, c) for (a, b, c) in Tp_all
                 if A[v][a] == 1 and A[c][v] == 1]
    else:
        valid = Tp_all

    if not valid:
        return np.zeros((len(paths2), 0))

    B = np.zeros((len(paths2), len(valid)))
    for j, (a, b, c) in enumerate(valid):
        terms = [
            ((v, b, c), -1), ((v, a, c), +1), ((v, a, b), -1),
            ((b, c, v), +1), ((a, c, v), -1), ((a, b, v), +1),
        ]
        for path, coeff in terms:
            if path in path2_idx:
                B[path2_idx[path], j] += coeff
    return B


# ============================================================
# PART 1: For each swap cycle, find ANY vertex that fills it
# ============================================================
print("=" * 70)
print("PART 1: CAN SOME OTHER VERTEX FILL EACH SWAP CYCLE?")
print("=" * 70)

for n in [7, 8]:
    random.seed(42)
    num_trials = 500 if n <= 7 else 500
    total_swaps = 0
    filled_by_own = 0
    filled_by_other = 0
    not_filled = 0

    t0 = time.time()
    for trial in range(num_trials):
        A = random_tournament(n)
        paths2 = enumerate_allowed_paths(A, n, 2)
        path2_idx = {tuple(p): i for i, p in enumerate(paths2)}

        # Pre-compute cone B for each vertex
        cone_Bs = {}
        for u in range(n):
            cone_Bs[u] = build_cone_B(A, n, u, paths2, path2_idx,
                                      filtered=True)

        for v in range(n):
            swaps = compute_swap_cycles(A, n, v, paths2, path2_idx)
            for z in swaps:
                total_swaps += 1

                # Try own vertex first
                B_v = cone_Bs[v]
                if B_v.shape[1] > 0:
                    alpha, _, _, _ = np.linalg.lstsq(B_v, z, rcond=None)
                    err = np.max(np.abs(B_v @ alpha - z))
                    if err < 1e-6:
                        filled_by_own += 1
                        continue

                # Try every other vertex
                found = False
                for u in range(n):
                    if u == v:
                        continue
                    B_u = cone_Bs[u]
                    if B_u.shape[1] == 0:
                        continue
                    alpha, _, _, _ = np.linalg.lstsq(B_u, z, rcond=None)
                    err = np.max(np.abs(B_u @ alpha - z))
                    if err < 1e-6:
                        filled_by_other += 1
                        found = True
                        break

                if not found:
                    not_filled += 1

        if (trial + 1) % 100 == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {trial + 1}/{num_trials} ({elapsed:.0f}s) "
                  f"swaps={total_swaps} own={filled_by_own} "
                  f"other={filled_by_other} none={not_filled}")

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): {total_swaps} swap cycles")
    print(f"  Filled by own vertex: {filled_by_own} "
          f"({100*filled_by_own/max(1,total_swaps):.1f}%)")
    print(f"  Filled by other vertex: {filled_by_other} "
          f"({100*filled_by_other/max(1,total_swaps):.1f}%)")
    print(f"  NOT filled by any single vertex: {not_filled} "
          f"({100*not_filled/max(1,total_swaps):.1f}%)")


# ============================================================
# PART 2: Multi-vertex cone rank analysis
# ============================================================
print(f"\n{'=' * 70}")
print("PART 2: MULTI-VERTEX CONE RANK")
print("=" * 70)

for n in [7, 8, 9]:
    random.seed(42)
    num_trials = 200 if n <= 8 else 50
    total_swaps = 0
    multi_ok = 0
    multi_fail = 0

    t0 = time.time()
    for trial in range(num_trials):
        A = random_tournament(n)
        paths2 = enumerate_allowed_paths(A, n, 2)
        path2_idx = {tuple(p): i for i, p in enumerate(paths2)}

        # Build multi-vertex B
        all_cols = []
        for u in range(n):
            B_u = build_cone_B(A, n, u, paths2, path2_idx, filtered=True)
            for j in range(B_u.shape[1]):
                all_cols.append(B_u[:, j])

        if not all_cols:
            continue

        B_multi = np.column_stack(all_cols)
        rank_multi = np.linalg.matrix_rank(B_multi, tol=1e-8)

        # Get all swap cycles
        all_swaps = []
        for v in range(n):
            swaps = compute_swap_cycles(A, n, v, paths2, path2_idx)
            all_swaps.extend(swaps)

        if not all_swaps:
            continue

        total_swaps += len(all_swaps)
        swap_mat = np.column_stack(all_swaps)
        swap_rank = np.linalg.matrix_rank(swap_mat, tol=1e-8)

        # Check swap space in col(B_multi)
        combined = np.hstack([B_multi, swap_mat])
        combined_rank = np.linalg.matrix_rank(combined, tol=1e-8)

        if combined_rank == rank_multi:
            multi_ok += 1
        else:
            multi_fail += 1
            if multi_fail <= 3:
                scores = tuple(sorted([sum(row) for row in A]))
                print(f"  n={n} FAIL trial {trial}: scores={scores}")
                print(f"    rank_multi={rank_multi}, swap_rank={swap_rank}, "
                      f"combined_rank={combined_rank}")

        if (trial + 1) % (num_trials // 4) == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {trial + 1}/{num_trials} ({elapsed:.0f}s) "
                  f"ok={multi_ok} fail={multi_fail}")

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): multi-vertex cone "
          f"OK={multi_ok}/{num_trials}, FAIL={multi_fail}")


# ============================================================
# PART 3: Unfiltered single-vertex cone (using ALL T' paths)
# ============================================================
print(f"\n{'=' * 70}")
print("PART 3: UNFILTERED CONE (ALL T' paths, not just v->a, c->v)")
print("=" * 70)
print("Even non-allowed cone paths may contribute correctly to B")

for n in [7, 8]:
    random.seed(42)
    num_trials = 500
    total_tested = 0
    unf_ok = 0
    unf_fail = 0

    t0 = time.time()
    for trial in range(num_trials):
        A = random_tournament(n)
        paths2 = enumerate_allowed_paths(A, n, 2)
        path2_idx = {tuple(p): i for i, p in enumerate(paths2)}

        any_fail = False
        for v in range(n):
            swaps = compute_swap_cycles(A, n, v, paths2, path2_idx)
            if not swaps:
                continue

            B_unf = build_cone_B(A, n, v, paths2, path2_idx, filtered=False)
            if B_unf.shape[1] == 0:
                any_fail = True
                break

            for z in swaps:
                total_tested += 1
                alpha, _, _, _ = np.linalg.lstsq(B_unf, z, rcond=None)
                err = np.max(np.abs(B_unf @ alpha - z))
                if err > 1e-6:
                    any_fail = True
                    break
            if any_fail:
                break

        if any_fail:
            unf_fail += 1
        else:
            unf_ok += 1

        if (trial + 1) % 100 == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {trial + 1}/{num_trials} ({elapsed:.0f}s) "
                  f"ok={unf_ok} fail={unf_fail}")

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): unfiltered single-vertex cone "
          f"OK={unf_ok}/{num_trials}, FAIL={unf_fail}")


print("\n\nDone.")
