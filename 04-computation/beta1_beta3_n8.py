"""
Check beta_1 and beta_3 mutual exclusivity at n=8.

At n=6,7: beta_1 and beta_3 are mutually exclusive.
At n=8: chi can be 2 (from beta_4=1), so theoretically beta_1+beta_3 > 1
is possible. Check if it happens.

Also: can beta_1 or beta_3 be > 1 at n=8?
"""
import numpy as np
from collections import defaultdict

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

def betti_single(A, n, target_p):
    allowed = {}
    for p in [target_p - 1, target_p, target_p + 1]:
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)
    omega_p = compute_omega_basis(A, n, target_p, allowed[target_p], allowed[target_p-1])
    dim_om = omega_p.shape[1] if omega_p.ndim == 2 else 0
    if dim_om == 0:
        return 0
    bd_p = build_boundary_matrix(allowed[target_p], allowed[target_p-1])
    bd_p_om = bd_p @ omega_p
    if bd_p_om.size > 0:
        sv = np.linalg.svd(bd_p_om, compute_uv=False)
        rk = sum(s > 1e-8 for s in sv)
    else:
        rk = 0
    ker = dim_om - rk
    omega_p1 = compute_omega_basis(A, n, target_p+1, allowed[target_p+1], allowed[target_p])
    dim_om1 = omega_p1.shape[1] if omega_p1.ndim == 2 else 0
    if dim_om1 > 0:
        bd1 = build_boundary_matrix(allowed[target_p+1], allowed[target_p])
        bd1_om = bd1 @ omega_p1
        sv1 = np.linalg.svd(bd1_om, compute_uv=False)
        im = sum(s > 1e-8 for s in sv1)
    else:
        im = 0
    return ker - im

def main():
    n = 8
    rng = np.random.RandomState(42)

    print("=" * 70)
    print(f"beta_1, beta_3 at n={n}: mutual exclusivity check")
    print("=" * 70)

    b1_vals = []
    b3_vals = []
    b1_and_b3 = 0
    total = 0

    for trial in range(300):
        A = random_tournament(n, rng)
        b1 = betti_single(A, n, 1)
        b3 = betti_single(A, n, 3)
        b1_vals.append(b1)
        b3_vals.append(b3)
        total += 1

        if b1 > 0 and b3 > 0:
            b1_and_b3 += 1
            scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
            print(f"  COOCCURRENCE at trial {trial}: b1={b1}, b3={b3}, scores={scores}")

        if trial % 50 == 49:
            print(f"  ... {trial+1} trials done")

    b1_pos = sum(1 for v in b1_vals if v > 0)
    b3_pos = sum(1 for v in b3_vals if v > 0)
    b1_max = max(b1_vals)
    b3_max = max(b3_vals)

    print(f"\nTotal: {total}")
    print(f"beta_1 > 0: {b1_pos}/{total} ({100*b1_pos/total:.1f}%)")
    print(f"beta_3 > 0: {b3_pos}/{total} ({100*b3_pos/total:.1f}%)")
    print(f"beta_1 AND beta_3 > 0: {b1_and_b3}/{total}")
    print(f"max beta_1: {b1_max}")
    print(f"max beta_3: {b3_max}")

    # Distribution
    b1_dist = defaultdict(int)
    b3_dist = defaultdict(int)
    for v in b1_vals:
        b1_dist[v] += 1
    for v in b3_vals:
        b3_dist[v] += 1
    print(f"\nbeta_1 distribution: {dict(sorted(b1_dist.items()))}")
    print(f"beta_3 distribution: {dict(sorted(b3_dist.items()))}")

if __name__ == '__main__':
    main()
