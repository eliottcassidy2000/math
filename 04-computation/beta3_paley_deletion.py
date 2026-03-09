"""
Check if n=6 beta_3=1 tournaments with H=45 are Paley T_7 deletions.

Known: At n=6, beta_3=1 has two score classes:
- Score (2,2,2,3,3,3): 240 tournaments with H=45
- Score (1,1,1,4,4,4): 80 tournaments with H=9

Paley T_7 - v gives a tournament on 6 vertices with score (2,2,3,3,3,x)
where x depends on the removed vertex. By vertex transitivity, all
deletions give isomorphic tournaments.
"""
import numpy as np
from itertools import permutations, combinations
from collections import defaultdict

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

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
    rank = sum(s > 1e-10 for s in S)
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

def path_betti_numbers(A, n, max_dim=None):
    if max_dim is None: max_dim = n - 1
    allowed = {}
    for p in range(-1, max_dim + 2):
        allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)
    omega = {}
    omega_dims = []
    for p in range(max_dim + 2):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
        omega_dims.append(omega[p].shape[1] if omega[p].ndim == 2 else 0)
    betti = []
    for p in range(max_dim + 1):
        dim_om = omega_dims[p]
        if dim_om == 0:
            betti.append(0)
            continue
        bd_p = build_boundary_matrix(allowed[p], allowed[p-1])
        bd_p_om = bd_p @ omega[p]
        if bd_p_om.size > 0:
            sv = np.linalg.svd(bd_p_om, compute_uv=False)
            rk = sum(s > 1e-8 for s in sv)
        else:
            rk = 0
        ker = dim_om - rk
        dim_om1 = omega_dims[p+1]
        if dim_om1 > 0:
            bd1 = build_boundary_matrix(allowed[p+1], allowed[p])
            bd1_om = bd1 @ omega[p+1]
            sv1 = np.linalg.svd(bd1_om, compute_uv=False)
            im = sum(s > 1e-8 for s in sv1)
        else:
            im = 0
        betti.append(ker - im)
    return betti, omega_dims[:max_dim+1]

def ham_path_count(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        bits = bin(mask).count('1')
        if bits == 1: continue
        for v in range(n):
            if not (mask & (1 << v)): continue
            prev_mask = mask ^ (1 << v)
            for u in range(n):
                if (prev_mask & (1 << u)) and A[u][v]:
                    if (mask, v) not in dp:
                        dp[(mask, v)] = 0
                    dp[(mask, v)] += dp.get((prev_mask, u), 0)
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_directed_cycles(A, n, k):
    count = 0
    for subset in combinations(range(n), k):
        for perm in permutations(subset):
            if all(A[perm[i]][perm[(i+1) % k]] for i in range(k)):
                count += 1
    return count // k

def delete_vertex(A, n, v):
    keep = [i for i in range(n) if i != v]
    B = np.zeros((n-1, n-1), dtype=int)
    for i, ki in enumerate(keep):
        for j, kj in enumerate(keep):
            B[i][j] = A[ki][kj]
    return B

def main():
    # Part 1: Paley T_7 deletion
    print("=" * 70)
    print("PART 1: Paley T_7 vertex deletion")
    print("=" * 70)
    n = 7
    A = np.zeros((n, n), dtype=int)
    qr = {1, 2, 4}
    for i in range(n):
        for j in range(n):
            if i != j and ((j - i) % n) in qr:
                A[i][j] = 1

    # Delete vertex 0
    B = delete_vertex(A, n, 0)
    scores = tuple(sorted(sum(B[i][j] for j in range(n-1)) for i in range(n-1)))
    H_del = ham_path_count(B, n-1)
    betti_del, odims_del = path_betti_numbers(B, n-1)
    c3_del = count_directed_cycles(B, n-1, 3)
    c5_del = count_directed_cycles(B, n-1, 5)
    print(f"T_7 - v0: scores={scores}, H={H_del}, c3={c3_del}, c5={c5_del}")
    print(f"  betti = {betti_del}")
    print(f"  omega_dims = {odims_del}")

    # Part 2: Find ALL n=6 tournaments with score (2,2,2,3,3,3) and H=45
    print("\n" + "=" * 70)
    print("PART 2: n=6 H=45 score (2,2,2,3,3,3) — Betti analysis")
    print("=" * 70)

    n6 = 6
    total = 2**(n6*(n6-1)//2)
    h45_count = 0
    h45_betti_dist = defaultdict(int)
    h45_odims_dist = defaultdict(int)
    h45_c3_dist = defaultdict(int)

    for bits in range(total):
        A6 = bits_to_adj(bits, n6)
        scores = tuple(sorted(sum(A6[i][j] for j in range(n6)) for i in range(n6)))
        if scores != (2,2,2,3,3,3):
            continue
        H = ham_path_count(A6, n6)
        if H != 45:
            continue
        h45_count += 1
        betti, odims = path_betti_numbers(A6, n6)
        h45_betti_dist[tuple(betti)] += 1
        h45_odims_dist[tuple(odims)] += 1
        c3 = count_directed_cycles(A6, n6, 3)
        h45_c3_dist[c3] += 1

    print(f"n=6 H=45 score (2,2,2,3,3,3): {h45_count} tournaments")
    print(f"Betti distribution: {dict(h45_betti_dist)}")
    print(f"Omega dims distribution: {dict(h45_odims_dist)}")
    print(f"c3 distribution: {dict(h45_c3_dist)}")

    # Part 3: n=6 H=45 — are they ALL Paley deletions?
    print("\n" + "=" * 70)
    print("PART 3: Are H=45 tournaments isomorphic to Paley deletion?")
    print("=" * 70)

    # Score of Paley deletion
    paley_del_scores = tuple(sorted(sum(B[i][j] for j in range(n-1)) for i in range(n-1)))
    print(f"Paley deletion scores: {paley_del_scores}")
    paley_del_c3 = count_directed_cycles(B, n-1, 3)
    paley_del_c5 = count_directed_cycles(B, n-1, 5)
    print(f"Paley deletion c3={paley_del_c3}, c5={paley_del_c5}, H={H_del}")

    # Part 4: n=6 score (1,1,1,4,4,4) H=9 — what are these?
    print("\n" + "=" * 70)
    print("PART 4: n=6 score (1,1,1,4,4,4) H=9 — structure")
    print("=" * 70)

    h9_count = 0
    h9_betti_dist = defaultdict(int)
    h9_odims_dist = defaultdict(int)

    for bits in range(total):
        A6 = bits_to_adj(bits, n6)
        scores = tuple(sorted(sum(A6[i][j] for j in range(n6)) for i in range(n6)))
        if scores != (1,1,1,4,4,4):
            continue
        H = ham_path_count(A6, n6)
        if H != 9:
            continue
        h9_count += 1
        betti, odims = path_betti_numbers(A6, n6)
        h9_betti_dist[tuple(betti)] += 1
        h9_odims_dist[tuple(odims)] += 1

        if h9_count <= 3:
            c3 = count_directed_cycles(A6, n6, 3)
            c5 = count_directed_cycles(A6, n6, 5)
            print(f"  bits={bits}: c3={c3}, c5={c5}, betti={list(betti)}, omega={list(odims)}")

    print(f"\nn=6 H=9 score (1,1,1,4,4,4): {h9_count} tournaments")
    print(f"Betti distribution: {dict(h9_betti_dist)}")
    print(f"Omega dims distribution: {dict(h9_odims_dist)}")

    # Part 5: Global summary — max H per score class and beta_3
    print("\n" + "=" * 70)
    print("PART 5: Max H per score class at n=6 and beta_3 status")
    print("=" * 70)

    max_H_by_score = defaultdict(int)
    beta3_by_score_H = defaultdict(set)

    for bits in range(total):
        A6 = bits_to_adj(bits, n6)
        scores = tuple(sorted(sum(A6[i][j] for j in range(n6)) for i in range(n6)))
        H = ham_path_count(A6, n6)
        max_H_by_score[scores] = max(max_H_by_score[scores], H)
        betti, _ = path_betti_numbers(A6, n6)
        if betti[3] > 0:
            beta3_by_score_H[scores].add(H)

    for scores in sorted(max_H_by_score.keys()):
        max_H = max_H_by_score[scores]
        b3_H = sorted(beta3_by_score_H[scores]) if scores in beta3_by_score_H else []
        marker = " <-- beta_3>0 = max_H!" if b3_H and max_H in b3_H else ""
        print(f"  {scores}: max_H={max_H}, beta_3>0 at H={b3_H}{marker}")

if __name__ == '__main__':
    main()
