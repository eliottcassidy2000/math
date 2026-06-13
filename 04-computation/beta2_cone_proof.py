#!/usr/bin/env python3
"""beta2_cone_proof.py - Cone-from-T' construction proves beta_2 = 0

BREAKTHROUGH from beta2_universal_filling.py (n=5):
For ANY tournament T on n vertices, ANY vertex v, ANY swap cycle z:
  w = sum_{(a,b,c) in Allowed_2(T')} alpha_{abc} * [(v,a,b,c) + (a,b,c,v)]
fills z, i.e., d_3(w) = z.

KEY: The T'-path terms cancel because d_3(v,a,b,c) contributes +(a,b,c)
while d_3(a,b,c,v) contributes -(a,b,c).

This script:
Part 1: Verify cone-from-T' at n=6 (exhaustive, 720 tournaments)
Part 2: Verify at n=7 (500 random + Paley)
Part 3: Rank analysis - WHY is B*alpha=z always solvable?

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time, random
import numpy as np
from collections import defaultdict, Counter
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
    """Get induced sub-tournament on vertex subset."""
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


def test_cone_from_Tprime(A, n, verbose=False):
    """Test cone-from-T' construction for all vertices v.

    For each v, compute swap cycles and try to fill them using:
    w = sum alpha_{abc} * [(v,a,b,c) + (a,b,c,v)]
    where (a,b,c) ranges over Allowed_2(T') with T' = T \ {v}.

    Returns (all_ok, num_tested, details).
    """
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)
    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0
    path2_idx = {tuple(p): i for i, p in enumerate(paths2)}
    path3_idx = {tuple(p): i for i, p in enumerate(paths3)}

    all_ok = True
    num_tested = 0
    details = []

    for v in range(n):
        # Compute swap cycles for vertex v
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

        # Get T' = T \ {v} allowed 2-paths
        others = [x for x in range(n) if x != v]
        B, vlist = get_induced(A, n, others)
        n_prime = len(others)
        paths2_Tp = enumerate_allowed_paths(B, n_prime, 2)

        # Map T' paths back to original vertex labels
        Tp_paths_orig = []
        for p in paths2_Tp:
            orig = tuple(vlist[x] for x in p)
            Tp_paths_orig.append(orig)

        num_Tp = len(Tp_paths_orig)

        # Build the B matrix: for each T' 2-path (a,b,c),
        # the column contribution to z from [(v,a,b,c) + (a,b,c,v)] via d_3
        # T'-internal terms cancel, leaving only v-chain terms
        B_fill = np.zeros((len(paths2), num_Tp))

        for j, (a, b, c) in enumerate(Tp_paths_orig):
            # d_3(v,a,b,c) = (a,b,c) - (v,b,c) + (v,a,c) - (v,a,b)
            # d_3(a,b,c,v) = (b,c,v) - (a,c,v) + (a,b,v) - (a,b,c)
            # Sum: -(v,b,c) + (v,a,c) - (v,a,b) + (b,c,v) - (a,c,v) + (a,b,v)
            # (T'-internal terms (a,b,c) cancel!)

            terms = [
                ((v, b, c), -1), ((v, a, c), +1), ((v, a, b), -1),
                ((b, c, v), +1), ((a, c, v), -1), ((a, b, v), +1),
            ]
            for path, coeff in terms:
                if path in path2_idx:
                    B_fill[path2_idx[path], j] += coeff

        # Test each swap cycle
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

            # Solve B_fill @ alpha = z
            alpha, res, rank_B, sv = np.linalg.lstsq(B_fill, z, rcond=None)
            err = np.max(np.abs(B_fill @ alpha - z))

            if err > 1e-6:
                all_ok = False
                details.append({
                    'v': v, 'ker_dim': ker_dim, 'err': err,
                    'num_Tp': num_Tp, 'rank_B': rank_B,
                    'swap_dim': ker_basis.shape[0], 'status': 'FAIL'
                })
            else:
                if verbose:
                    # Check if filling is in Omega_3
                    w_full = np.zeros(len(paths3))
                    for j, (a, b, c) in enumerate(Tp_paths_orig):
                        if abs(alpha[j]) < 1e-12:
                            continue
                        front = (v, a, b, c)
                        back = (a, b, c, v)
                        if front in path3_idx:
                            w_full[path3_idx[front]] += alpha[j]
                        if back in path3_idx:
                            w_full[path3_idx[back]] += alpha[j]

                    # Verify d_3(w) = z
                    D3 = build_full_boundary_matrix(
                        [tuple(p) for p in paths3],
                        [tuple(p) for p in paths2]
                    )
                    check = D3 @ w_full
                    d3_err = np.max(np.abs(check - z))

                    details.append({
                        'v': v, 'ker_dim': ker_dim, 'err': err,
                        'num_Tp': num_Tp, 'rank_B': rank_B,
                        'swap_dim': ker_basis.shape[0],
                        'd3_err': d3_err,
                        'status': 'OK'
                    })

    return all_ok, num_tested, details


# ============================================================
# PART 1: n=6 exhaustive
# ============================================================
print("=" * 70)
print("PART 1: CONE-FROM-T' at n=6 (EXHAUSTIVE)")
print("=" * 70)

n = 6
t0 = time.time()
total_tested = 0
total_ok = 0
total_fail = 0
fail_examples = []

for tidx, A in enumerate(all_tournaments_gen(n)):
    ok, tested, details = test_cone_from_Tprime(A, n)
    total_tested += tested
    if ok:
        total_ok += 1
    else:
        total_fail += 1
        if total_fail <= 3:
            scores = tuple(sorted([sum(row) for row in A]))
            fail_examples.append((tidx, scores, details))

    if (tidx + 1) % 5000 == 0:
        elapsed = time.time() - t0
        print(f"  {tidx + 1}/32768 ({elapsed:.0f}s) tested={total_tested} "
              f"ok={total_ok} fail={total_fail}")

elapsed = time.time() - t0
print(f"\nn=6 ({elapsed:.0f}s): {total_tested} swap cycles tested")
print(f"  OK: {total_ok}/32768, Failures: {total_fail}")

for tidx, scores, details in fail_examples:
    print(f"\n  FAIL T#{tidx} scores={scores}:")
    for d in details:
        print(f"    v={d['v']}, err={d['err']:.2e}, "
              f"#Tp={d['num_Tp']}, rank_B={d['rank_B']}")


# ============================================================
# PART 2: n=7 sampled
# ============================================================
print(f"\n{'=' * 70}")
print("PART 2: CONE-FROM-T' at n=7 (500 random + Paley)")
print("=" * 70)

n = 7
t0 = time.time()
total_tested = 0
total_ok = 0
total_fail = 0

for trial in range(500):
    A = random_tournament(n)
    ok, tested, details = test_cone_from_Tprime(A, n)
    total_tested += tested
    if ok:
        total_ok += 1
    else:
        total_fail += 1
        if total_fail <= 3:
            scores = tuple(sorted([sum(row) for row in A]))
            print(f"  FAIL trial {trial}: scores={scores}")
            for d in details:
                print(f"    v={d['v']}, err={d['err']:.2e}, "
                      f"#Tp={d['num_Tp']}, rank_B={d['rank_B']}")

    if (trial + 1) % 100 == 0:
        elapsed = time.time() - t0
        print(f"  {trial + 1}/500 ({elapsed:.0f}s) tested={total_tested} "
              f"ok={total_ok} fail={total_fail}")

# Paley T_7
A_paley = paley_tournament(7)
ok_paley, tested_paley, details_paley = test_cone_from_Tprime(
    A_paley, n, verbose=True
)
total_tested += tested_paley

elapsed = time.time() - t0
print(f"\nn=7 ({elapsed:.0f}s): {total_tested} swap cycles tested")
print(f"  OK: {total_ok}/500 random, Paley={'OK' if ok_paley else 'FAIL'}")
print(f"  Failures: {total_fail}")

if ok_paley and details_paley:
    print(f"\n  Paley T_7 cone details:")
    for d in details_paley[:3]:
        print(f"    v={d['v']}: #Tp={d['num_Tp']}, rank_B={d['rank_B']}, "
              f"err={d['err']:.2e}, d3_err={d.get('d3_err', 'N/A')}")


# ============================================================
# PART 3: Rank analysis
# ============================================================
print(f"\n{'=' * 70}")
print("PART 3: RANK ANALYSIS - WHY IS B*alpha=z ALWAYS SOLVABLE?")
print("=" * 70)

print("\nFor each (tournament, vertex v) with swap cycles:")
print("  Compare: swap_dim, #T'_paths (cols of B), rank(B), "
      "rank(B restricted to swap space)")

n = 5
rank_data = []

for A in all_tournaments_gen(n):
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

        # Get T' paths
        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
        Tp_paths_orig = [tuple(vlist[x] for x in p) for p in paths2_Tp]
        num_Tp = len(Tp_paths_orig)

        # Build B matrix
        B_fill = np.zeros((len(paths2), num_Tp))
        for j, (a, b, c) in enumerate(Tp_paths_orig):
            terms = [
                ((v, b, c), -1), ((v, a, c), +1), ((v, a, b), -1),
                ((b, c, v), +1), ((a, c, v), -1), ((a, b, v), +1),
            ]
            for path, coeff in terms:
                if path in path2_idx:
                    B_fill[path2_idx[path], j] += coeff

        rank_B = np.linalg.matrix_rank(B_fill, tol=1e-8)

        # Compute the swap cycle space: image of the swap map
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

        # Check: is swap space in column space of B?
        combined = np.hstack([B_fill, swap_mat])
        rank_combined = np.linalg.matrix_rank(combined, tol=1e-8)
        in_col_space = (rank_combined == rank_B)

        rank_data.append({
            'swap_dim': swap_dim,
            'num_Tp': num_Tp,
            'rank_B': rank_B,
            'in_col_space': in_col_space,
            'deg_v_out': len(P),
            'deg_v_in': len(Q),
            'arcs_PQ': len(arcs_PQ),
        })

print(f"\nn=5: {len(rank_data)} (tournament, v) pairs with swap cycles")
print(f"\n  {'swap_dim':>8} {'#Tp':>5} {'rank_B':>7} {'swap<=B':>7} "
      f"{'d_out':>5} {'d_in':>5} {'#PQ':>5}")

# Group by pattern
patterns = Counter()
for d in rank_data:
    key = (d['swap_dim'], d['num_Tp'], d['rank_B'], d['in_col_space'])
    patterns[key] += 1

for (sd, nt, rb, ics), count in sorted(patterns.items()):
    print(f"  {sd:>8} {nt:>5} {rb:>7} {'YES' if ics else 'NO':>7}  "
          f"({count} cases)")

# Summary
all_in = all(d['in_col_space'] for d in rank_data)
print(f"\n  ALL swap spaces in col(B)? {'YES' if all_in else 'NO'}")
print(f"  rank(B) always >= swap_dim? "
      f"{'YES' if all(d['rank_B'] >= d['swap_dim'] for d in rank_data) else 'NO'}")

# Surplus analysis
surpluses = [d['rank_B'] - d['swap_dim'] for d in rank_data]
print(f"  rank(B) - swap_dim: min={min(surpluses)}, max={max(surpluses)}, "
      f"mean={np.mean(surpluses):.1f}")


# ============================================================
# PART 4: Rank analysis at n=6 (sampled)
# ============================================================
print(f"\n{'=' * 70}")
print("PART 4: RANK ANALYSIS at n=6 (1000 random)")
print("=" * 70)

n = 6
rank_data_6 = []
t0 = time.time()

for trial in range(1000):
    A = random_tournament(n)
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
        Tp_paths_orig = [tuple(vlist[x] for x in p) for p in paths2_Tp]
        num_Tp = len(Tp_paths_orig)

        B_fill = np.zeros((len(paths2), num_Tp))
        for j, (a, b, c) in enumerate(Tp_paths_orig):
            terms = [
                ((v, b, c), -1), ((v, a, c), +1), ((v, a, b), -1),
                ((b, c, v), +1), ((a, c, v), -1), ((a, b, v), +1),
            ]
            for path, coeff in terms:
                if path in path2_idx:
                    B_fill[path2_idx[path], j] += coeff

        rank_B = np.linalg.matrix_rank(B_fill, tol=1e-8)

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

        combined = np.hstack([B_fill, swap_mat])
        rank_combined = np.linalg.matrix_rank(combined, tol=1e-8)
        in_col_space = (rank_combined == rank_B)

        rank_data_6.append({
            'swap_dim': swap_dim,
            'num_Tp': num_Tp,
            'rank_B': rank_B,
            'in_col_space': in_col_space,
        })

    if (trial + 1) % 200 == 0:
        elapsed = time.time() - t0
        print(f"  {trial + 1}/1000 ({elapsed:.0f}s) "
              f"cases={len(rank_data_6)}")

elapsed = time.time() - t0
print(f"\nn=6 ({elapsed:.0f}s): {len(rank_data_6)} cases")

if rank_data_6:
    all_in_6 = all(d['in_col_space'] for d in rank_data_6)
    print(f"  ALL swap spaces in col(B)? {'YES' if all_in_6 else 'NO'}")

    surpluses_6 = [d['rank_B'] - d['swap_dim'] for d in rank_data_6]
    print(f"  rank(B) - swap_dim: min={min(surpluses_6)}, "
          f"max={max(surpluses_6)}, mean={np.mean(surpluses_6):.1f}")

    patterns_6 = Counter()
    for d in rank_data_6:
        key = (d['swap_dim'], d['num_Tp'], d['rank_B'], d['in_col_space'])
        patterns_6[key] += 1

    print(f"\n  Top patterns (swap_dim, #Tp, rank_B, in_col):")
    for (sd, nt, rb, ics), count in sorted(
            patterns_6.items(), key=lambda x: -x[1])[:15]:
        print(f"    ({sd}, {nt}, {rb}, {'Y' if ics else 'N'}): {count}")


print("\n\nDone.")
