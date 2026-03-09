"""
Sample n=8 tournaments for Betti numbers.

At n=8, the H-maximizer has H=661 (OEIS A038375).
Question: Does it have nonzero beta_5 (= beta_{n-3})?

We sample n=8 tournaments, focusing on H-maximizers.
n=8 is expensive, so we sample carefully.
"""
import numpy as np
from itertools import combinations
from collections import defaultdict
import sys

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
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

def main():
    rng = np.random.RandomState(42)
    n = 8

    # Part 1: The known n=8 H-maximizer
    # SC maximizer at n=8 has H=661, score (3,3,3,3,4,4,4,4)
    # We need to find/construct it
    print("=" * 70)
    print("PART 1: Searching for n=8 H-maximizer (H=661)")
    print("=" * 70)

    max_H = 0
    max_A = None
    max_scores = None

    # Sample a lot to find high-H tournaments
    for trial in range(5000):
        A = random_tournament(n, rng)
        H = ham_path_count(A, n)
        if H > max_H:
            max_H = H
            max_A = A.copy()
            max_scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
            print(f"  trial {trial}: H={H}, scores={max_scores}")

        if trial % 1000 == 999:
            print(f"  ... {trial+1} trials done, max H so far: {max_H}")

    print(f"\nBest found: H={max_H}, scores={max_scores}")

    if max_H >= 600:
        print(f"\nComputing Betti numbers for best tournament (H={max_H})...")
        # Only compute up to beta_5 to save time
        betti, odims = path_betti_numbers(max_A, n, max_dim=min(7, n-1))
        print(f"  betti = {betti}")
        print(f"  omega_dims = {odims}")
        chi = sum((-1)**p * betti[p] for p in range(len(betti)))
        print(f"  chi = {chi}")

    # Part 2: Paley T_7 deletion at n=6, then compute Betti at n=8 via
    # a different approach: use deletion of Paley.
    # Actually, at n=8, Paley T_p doesn't exist (8 is not prime).
    # The H-maximizer at n=8 is known: H=661.

    # Part 3: Random n=8 tournaments — just check a few for beta structure
    print("\n" + "=" * 70)
    print("PART 3: Random n=8 — beta structure (5 samples)")
    print("=" * 70)

    for trial in range(5):
        A = random_tournament(n, rng)
        H = ham_path_count(A, n)
        scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        print(f"\n  trial {trial}: H={H}, scores={scores}")
        betti, odims = path_betti_numbers(A, n, max_dim=5)
        print(f"  betti (up to p=5) = {betti}")
        print(f"  omega_dims (up to p=5) = {odims}")

if __name__ == '__main__':
    main()
