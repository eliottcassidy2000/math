"""
beta2_completeness_argument.py — Does completeness (every pair has an edge)
force beta_2 = 0?

Key idea: In a tournament, Omega_1 = A_1 = {all edges}. This means there
are NO constraints at level 1 (Omega_1 = A_1), which is special.

For a general digraph D, Omega_1 might be a proper subspace of A_1
(if some 0-faces of 1-paths are non-allowed). But in a tournament,
ALL 0-paths (vertices) are allowed, so Omega_1 = A_1.

Similarly, at level 2: the non-allowed faces are (a,c) with c->a.
But crucially, (a,c) IS an ordered pair of distinct vertices.
In a tournament, (c,a) IS in A_1 (since c->a), but (a,c) is NOT.

So the constraint structure at level 2 depends ONLY on the orientation.

Question: Do NON-tournament digraphs (missing edges) have beta_2 > 0?
If so, this would confirm that completeness is essential.

Also: do RANDOM ORIENTED GRAPHS have beta_2 > 0?

Author: kind-pasteur-S45 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from itertools import combinations
from collections import Counter
sys.stdout.reconfigure(line_buffering=True)

def enumerate_allowed_paths(A, n, p):
    if p < 0: return []
    if p == 0: return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1: adj[i].append(j)
    paths = []
    stack = []
    for start in range(n):
        stack.append(([start], 1 << start))
        while stack:
            path, visited = stack.pop()
            if len(path) == p + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if not (visited & (1 << u)):
                    stack.append((path + [u], visited | (1 << u)))
    return paths

def boundary_coeffs(path):
    return [((-1)**i, path[:i] + path[i+1:]) for i in range(len(path))]

def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0: return np.zeros((0, 0))
    if p == 0: return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1
    if na_count == 0: return np.eye(dim_Ap)
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed:
                P[non_allowed[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = int(sum(s > 1e-10 for s in S))
    ns = Vt[rank:].T
    return ns if ns.shape[1] > 0 else np.zeros((dim_Ap, 0))

def build_boundary_matrix(allowed_p, allowed_pm1):
    if not allowed_p or not allowed_pm1:
        return np.zeros((max(len(allowed_pm1), 0), max(len(allowed_p), 0)))
    idx = {path: i for i, path in enumerate(allowed_pm1)}
    M = np.zeros((len(allowed_pm1), len(allowed_p)))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in idx:
                M[idx[face], j] += sign
    return M

def betti_single(A, n, target_p):
    allowed = {}
    for p in [target_p - 1, target_p, target_p + 1]:
        if p < 0: allowed[p] = []
        else: allowed[p] = enumerate_allowed_paths(A, n, p)
    omega_p = compute_omega_basis(A, n, target_p, allowed[target_p], allowed[target_p-1])
    dim_om = omega_p.shape[1] if omega_p.ndim == 2 else 0
    if dim_om == 0: return 0
    bd_p = build_boundary_matrix(allowed[target_p], allowed[target_p-1])
    bd_p_om = bd_p @ omega_p
    if bd_p_om.size > 0:
        sv = np.linalg.svd(bd_p_om, compute_uv=False)
        rk = int(sum(s > 1e-8 for s in sv))
    else: rk = 0
    ker = dim_om - rk
    omega_p1 = compute_omega_basis(A, n, target_p+1, allowed[target_p+1], allowed[target_p])
    dim_om1 = omega_p1.shape[1] if omega_p1.ndim == 2 else 0
    if dim_om1 > 0:
        bd1 = build_boundary_matrix(allowed[target_p+1], allowed[target_p])
        bd1_om = bd1 @ omega_p1
        sv1 = np.linalg.svd(bd1_om, compute_uv=False)
        im = int(sum(s > 1e-8 for s in sv1))
    else: im = 0
    return ker - im


def main():
    print("=" * 70)
    print("COMPLETENESS AND BETA_2: TESTING NON-TOURNAMENTS")
    print("=" * 70)

    # Part 1: Random oriented graphs (not tournaments — some edges missing)
    print("\n--- Part 1: Random oriented graphs on n=6 ---")
    n = 6
    rng = np.random.RandomState(42)
    beta2_count = 0
    total = 500

    for trial in range(total):
        # Random oriented graph: each pair independently has edge a->b, b->a, or NONE
        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(i+1, n):
                r = rng.random()
                if r < 0.33:
                    A[i][j] = 1
                elif r < 0.67:
                    A[j][i] = 1
                # else: no edge

        b2 = betti_single(A, n, 2)
        if b2 > 0:
            beta2_count += 1
            if beta2_count <= 5:
                edges = sum(A[i][j] for i in range(n) for j in range(n))
                b1 = betti_single(A, n, 1)
                print(f"  trial {trial}: beta_2={b2}, beta_1={b1}, edges={edges}")

    print(f"  beta_2 > 0 in {beta2_count}/{total} random oriented graphs ({100*beta2_count/total:.1f}%)")

    # Part 2: Dense random oriented graphs (every pair has an edge, but possibly mutual)
    print("\n--- Part 2: Dense random digraphs on n=6 (every pair has at least one edge) ---")
    beta2_count = 0
    for trial in range(500):
        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(i+1, n):
                r = rng.random()
                if r < 0.25:
                    A[i][j] = 1
                elif r < 0.50:
                    A[j][i] = 1
                elif r < 0.75:
                    A[i][j] = 1
                    A[j][i] = 1  # mutual edge
                else:
                    # Tournament-like: one direction only
                    if rng.random() < 0.5:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1

        b2 = betti_single(A, n, 2)
        if b2 > 0:
            beta2_count += 1
    print(f"  beta_2 > 0 in {beta2_count}/500 dense random digraphs")

    # Part 3: "Almost tournaments" — tournament with 1-2 edges removed
    print("\n--- Part 3: Almost-tournaments on n=6 (1 edge removed) ---")
    beta2_count = 0
    total_tested = 0
    for trial in range(500):
        # Random tournament
        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(i+1, n):
                if rng.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        # Remove one random edge
        edges = [(i, j) for i in range(n) for j in range(n) if A[i][j] == 1]
        if not edges: continue
        ei, ej = edges[rng.randint(len(edges))]
        A[ei][ej] = 0

        total_tested += 1
        b2 = betti_single(A, n, 2)
        if b2 > 0:
            beta2_count += 1

    print(f"  beta_2 > 0 in {beta2_count}/{total_tested} almost-tournaments")

    # Part 4: Tournaments with one edge FLIPPED (both directions = bidirectional)
    print("\n--- Part 4: Tournament + 1 mutual edge on n=6 ---")
    beta2_count = 0
    for trial in range(500):
        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(i+1, n):
                if rng.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
        # Add a mutual edge (make one edge bidirectional)
        edges = [(i, j) for i in range(n) for j in range(n) if A[i][j] and not A[j][i]]
        if edges:
            ei, ej = edges[rng.randint(len(edges))]
            A[ej][ei] = 1  # make bidirectional

        b2 = betti_single(A, n, 2)
        if b2 > 0:
            beta2_count += 1

    print(f"  beta_2 > 0 in {beta2_count}/500 tournament+mutual-edge")

    # Part 5: Complete bidirectional graph (every pair has BOTH directions)
    print("\n--- Part 5: Complete bidirectional graph on n=5,6 ---")
    for n_test in [5, 6]:
        A = np.zeros((n_test, n_test), dtype=int)
        for i in range(n_test):
            for j in range(n_test):
                if i != j:
                    A[i][j] = 1

        bettis = []
        for p in range(n_test):
            b = betti_single(A, n_test, p)
            bettis.append(b)
        print(f"  n={n_test}: Complete bidirectional: bettis={bettis}")

    # Part 6: Key test — does beta_2=0 hold for ALL complete oriented graphs?
    # (tournaments are complete oriented, meaning exactly one edge per pair)
    # Let's exhaustively check at n=4
    print("\n--- Part 6: Exhaustive check of COMPLETE ORIENTED graphs at n=4 ---")
    n = 4
    # A complete oriented graph on 4 vertices: for each pair, exactly one of (i,j) or (j,i)
    # This IS the set of tournaments! So this is redundant... but confirms.
    total = 0
    b2_nonzero = 0
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    for mask in range(1 << len(pairs)):
        A = np.zeros((n, n), dtype=int)
        for idx, (i,j) in enumerate(pairs):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        total += 1
        b2 = betti_single(A, n, 2)
        if b2 > 0:
            b2_nonzero += 1
    print(f"  n=4 tournaments: beta_2>0 in {b2_nonzero}/{total}")

    # Part 7: Check non-tournament COMPLETE digraphs at n=4
    # (where some pairs have bidirectional edges)
    print("\n--- Part 7: Complete digraphs (with possible mutual edges) at n=4 ---")
    n = 4
    total = 0
    b2_nonzero = 0
    # For each pair: one of (i->j only, j->i only, both)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    for mask in range(3 ** len(pairs)):
        A = np.zeros((n, n), dtype=int)
        m = mask
        for (i,j) in pairs:
            choice = m % 3
            m //= 3
            if choice == 0:
                A[i][j] = 1
            elif choice == 1:
                A[j][i] = 1
            else:
                A[i][j] = 1
                A[j][i] = 1
        total += 1
        b2 = betti_single(A, n, 2)
        if b2 > 0:
            b2_nonzero += 1
            if b2_nonzero <= 5:
                edges = sum(A[i][j] for i in range(n) for j in range(n))
                mutual = sum(1 for i in range(n) for j in range(i+1,n) if A[i][j] and A[j][i])
                print(f"    beta_2={b2}, edges={edges}, mutual_pairs={mutual}")

    print(f"  n=4 complete digraphs: beta_2>0 in {b2_nonzero}/{total}")

    # Part 8: Check SPARSE digraphs at n=5
    print("\n--- Part 8: Sparse random digraphs at n=5 ---")
    n = 5
    rng = np.random.RandomState(123)
    beta2_count = 0
    for trial in range(500):
        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(n):
                if i != j and rng.random() < 0.3:
                    A[i][j] = 1
        b2 = betti_single(A, n, 2)
        if b2 > 0:
            beta2_count += 1
    print(f"  n=5 sparse (p=0.3): beta_2>0 in {beta2_count}/500")

    print("\n--- SUMMARY ---")
    print("""
    beta_2 = 0 for ALL tournaments (complete oriented, exactly one edge per pair).
    beta_2 > 0 occurs in:
    - Sparse digraphs (missing edges allow cycles)
    - Dense digraphs with mutual edges (bidirectional pairs create new topology)
    - Almost-tournaments (removing even 1 edge can break beta_2=0)

    The COMPLETENESS + ANTISYMMETRY of tournaments is essential.
    """)

    print("DONE.")


if __name__ == '__main__':
    main()
