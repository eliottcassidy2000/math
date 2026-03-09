#!/usr/bin/env python3
"""beta2_filtered_cone.py - Filtered cone construction for beta_2 = 0

CORRECTED from beta2_cone_proof.py: only use T' 2-paths (a,b,c) where
BOTH (v,a,b,c) and (a,b,c,v) are allowed 3-paths in T.

This requires:
- v -> a (for front cone)
- c -> v (for back cone)

These are exactly the T' 2-paths (a,b,c) where a in P(v) and c in Q(v),
i.e., the "doubly-reachable" paths starting from out-neighborhood and
ending in in-neighborhood of v.

Part 1: n=5 exhaustive
Part 2: n=6 exhaustive
Part 3: n=7 sampled + Paley
Part 4: Rank analysis and surplus structure

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


def all_tournaments_gen(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0] * n for _ in range(n)]
        for idx, (i, j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


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
    for i in range(len(path) - 1):
        if A[path[i]][path[i + 1]] != 1:
            return False
    return True


def test_filtered_cone(A, n, verbose=False):
    """Test filtered cone-from-T' for all vertices v.

    For each v, only use T' 2-paths (a,b,c) where both
    (v,a,b,c) and (a,b,c,v) are allowed 3-paths.

    Returns (all_ok, num_tested, details).
    """
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)
    path2_idx = {tuple(p): i for i, p in enumerate(paths2)}
    path3_idx = {tuple(p): i for i, p in enumerate(paths3)}

    D3 = build_full_boundary_matrix(
        [tuple(p) for p in paths3], [tuple(p) for p in paths2]
    )

    all_ok = True
    num_tested = 0
    details = []

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

        _, _, Vt = np.linalg.svd(C_mat, full_matrices=True)
        ker_basis = Vt[rank_C:]

        # Get T' 2-paths, filtered to those where both cone paths are allowed
        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
        Tp_all = [tuple(vlist[x] for x in p) for p in paths2_Tp]

        # Filter: need v->a AND c->v
        valid_Tp = [(a, b, c) for (a, b, c) in Tp_all
                    if A[v][a] == 1 and A[c][v] == 1]

        if len(valid_Tp) == 0:
            # No valid cone paths - check if there are swap cycles
            has_nonzero = False
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
                    has_nonzero = True
                    break
            if has_nonzero:
                all_ok = False
                details.append({
                    'v': v, 'status': 'NO_VALID_TP',
                    'num_Tp_all': len(Tp_all), 'num_valid': 0,
                })
            continue

        # Build B matrix using only valid T' paths
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

            if np.max(np.abs(z)) < 1e-12:
                continue

            num_tested += 1

            alpha, _, _, _ = np.linalg.lstsq(B_fill, z, rcond=None)
            B_err = np.max(np.abs(B_fill @ alpha - z))

            if B_err > 1e-6:
                all_ok = False
                details.append({
                    'v': v, 'status': 'FAIL',
                    'B_err': B_err, 'rank_B': rank_B,
                    'num_valid': len(valid_Tp),
                })
                continue

            # Verify: build 3-chain and check d_3(w) = z
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

            if d3_err > 1e-6:
                all_ok = False
                details.append({
                    'v': v, 'status': 'D3_FAIL',
                    'd3_err': d3_err, 'rank_B': rank_B,
                })
            elif verbose:
                details.append({
                    'v': v, 'status': 'OK',
                    'B_err': B_err, 'd3_err': d3_err,
                    'rank_B': rank_B, 'num_valid': len(valid_Tp),
                    'swap_dim': ker_basis.shape[0],
                })

    return all_ok, num_tested, details


# ============================================================
# PART 1: n=5 exhaustive
# ============================================================
print("=" * 70)
print("PART 1: FILTERED CONE at n=5 (EXHAUSTIVE)")
print("=" * 70)

n = 5
t0 = time.time()
total_tested = 0
total_ok = 0
total_fail = 0

for A in all_tournaments_gen(n):
    ok, tested, details = test_filtered_cone(A, n)
    total_tested += tested
    if ok:
        total_ok += 1
    else:
        total_fail += 1
        scores = tuple(sorted([sum(row) for row in A]))
        if total_fail <= 3:
            print(f"  FAIL scores={scores}: {details}")

elapsed = time.time() - t0
print(f"\nn=5 ({elapsed:.0f}s): {total_tested} swap cycles")
print(f"  OK: {total_ok}/1024, Failures: {total_fail}")


# ============================================================
# PART 2: n=6 exhaustive
# ============================================================
print(f"\n{'=' * 70}")
print("PART 2: FILTERED CONE at n=6 (EXHAUSTIVE)")
print("=" * 70)

n = 6
t0 = time.time()
total_tested = 0
total_ok = 0
total_fail = 0

for tidx, A in enumerate(all_tournaments_gen(n)):
    ok, tested, details = test_filtered_cone(A, n)
    total_tested += tested
    if ok:
        total_ok += 1
    else:
        total_fail += 1
        if total_fail <= 5:
            scores = tuple(sorted([sum(row) for row in A]))
            print(f"  FAIL T#{tidx} scores={scores}:")
            for d in details:
                print(f"    {d}")

    if (tidx + 1) % 5000 == 0:
        elapsed = time.time() - t0
        print(f"  {tidx + 1}/32768 ({elapsed:.0f}s) tested={total_tested} "
              f"ok={total_ok} fail={total_fail}")

elapsed = time.time() - t0
print(f"\nn=6 ({elapsed:.0f}s): {total_tested} swap cycles")
print(f"  OK: {total_ok}/32768, Failures: {total_fail}")


# ============================================================
# PART 3: n=7 sampled
# ============================================================
print(f"\n{'=' * 70}")
print("PART 3: FILTERED CONE at n=7 (500 random + Paley)")
print("=" * 70)

n = 7
t0 = time.time()
total_tested = 0
total_ok = 0
total_fail = 0

for trial in range(500):
    A = random_tournament(n)
    ok, tested, details = test_filtered_cone(A, n)
    total_tested += tested
    if ok:
        total_ok += 1
    else:
        total_fail += 1
        if total_fail <= 3:
            scores = tuple(sorted([sum(row) for row in A]))
            print(f"  FAIL trial {trial}: scores={scores}")
            for d in details:
                print(f"    {d}")

    if (trial + 1) % 100 == 0:
        elapsed = time.time() - t0
        print(f"  {trial + 1}/500 ({elapsed:.0f}s) tested={total_tested} "
              f"ok={total_ok} fail={total_fail}")

# Paley
A_paley = paley_tournament(7)
ok_paley, tested_paley, det_paley = test_filtered_cone(
    A_paley, n, verbose=True
)
total_tested += tested_paley

elapsed = time.time() - t0
print(f"\nn=7 ({elapsed:.0f}s): {total_tested} swap cycles")
print(f"  OK: {total_ok}/500 random, Paley={'OK' if ok_paley else 'FAIL'}")
print(f"  Failures: {total_fail}")

if ok_paley and det_paley:
    print(f"\n  Paley details:")
    for d in det_paley[:3]:
        print(f"    v={d['v']}: valid={d['num_valid']}, "
              f"rank={d['rank_B']}, B_err={d['B_err']:.2e}, "
              f"d3_err={d['d3_err']:.2e}")


# ============================================================
# PART 4: Rank surplus analysis
# ============================================================
print(f"\n{'=' * 70}")
print("PART 4: RANK SURPLUS at n=5,6")
print("=" * 70)

for n in [5, 6]:
    rank_data = []
    gen = all_tournaments_gen(n)

    for A in gen:
        paths2 = enumerate_allowed_paths(A, n, 2)
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

            _, _, Vt = np.linalg.svd(C_mat, full_matrices=True)
            ker_basis = Vt[rank_C:]

            others = [x for x in range(n) if x != v]
            B_sub, vlist = get_induced(A, n, others)
            paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
            Tp_all = [tuple(vlist[x] for x in p) for p in paths2_Tp]
            valid_Tp = [(a, b, c) for (a, b, c) in Tp_all
                        if A[v][a] == 1 and A[c][v] == 1]

            if not valid_Tp:
                continue

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

            # Compute actual swap dim
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

            if not swap_vecs:
                continue

            swap_mat = np.column_stack(swap_vecs)
            swap_dim = np.linalg.matrix_rank(swap_mat, tol=1e-8)

            rank_data.append({
                'swap_dim': swap_dim, 'num_valid': len(valid_Tp),
                'num_all': len(Tp_all), 'rank_B': rank_B,
                'd_out': len(P), 'd_in': len(Q),
            })

    print(f"\nn={n}: {len(rank_data)} cases with swap cycles")
    if rank_data:
        surpluses = [d['rank_B'] - d['swap_dim'] for d in rank_data]
        print(f"  rank(B) - swap_dim: min={min(surpluses)}, "
              f"max={max(surpluses)}, mean={np.mean(surpluses):.1f}")
        print(f"  ALL rank(B) >= swap_dim? "
              f"{'YES' if min(surpluses) >= 0 else 'NO'}")

        patterns = Counter()
        for d in rank_data:
            key = (d['swap_dim'], d['num_valid'], d['rank_B'])
            patterns[key] += 1

        print(f"\n  Top patterns (swap_dim, #valid_Tp, rank_B):")
        for (sd, nv, rb), count in sorted(
                patterns.items(), key=lambda x: -x[1])[:10]:
            print(f"    ({sd}, {nv}, {rb}): {count}  surplus={rb - sd}")


print("\n\nDone.")
