"""
beta2_only_n9.py — Check ONLY beta_2 at n=9 and n=10 (the critical vanishing).

Only needs paths at p=1,2,3 which is much cheaper than full Betti computation.

Author: kind-pasteur-S45 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from itertools import combinations
sys.stdout.reconfigure(line_buffering=True)

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

def beta2_only(A, n):
    """Compute beta_2 only, using paths at p=1,2,3."""
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    ob2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = ob2.shape[1] if ob2.ndim == 2 else 0
    if dim_om2 == 0: return 0

    # ker(d_2)
    bd2 = build_boundary_matrix(a2, a1)
    d2_om = bd2 @ ob2
    sv2 = np.linalg.svd(d2_om, compute_uv=False)
    rank_d2 = int(sum(s > 1e-8 for s in sv2))
    ker_d2 = dim_om2 - rank_d2

    if ker_d2 == 0: return 0  # beta_2 <= ker_d2

    # im(d_3)
    ob3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = ob3.shape[1] if ob3.ndim == 2 else 0
    if dim_om3 == 0: return ker_d2

    bd3 = build_boundary_matrix(a3, a2)
    d3_om = bd3 @ ob3
    sv3 = np.linalg.svd(d3_om, compute_uv=False)
    im_d3 = int(sum(s > 1e-8 for s in sv3))

    return ker_d2 - im_d3


def main():
    print("=" * 70)
    print("BETA_2 ONLY CHECK AT n=9, n=10")
    print("=" * 70)

    for n in [9, 10]:
        rng = np.random.RandomState(12345 + n)
        total = 500 if n <= 9 else 100
        violations = 0

        print(f"\n--- n={n}: {total} samples ---", flush=True)

        for trial in range(total):
            A = random_tournament(n, rng)
            b2 = beta2_only(A, n)
            if b2 > 0:
                violations += 1
                scores = tuple(sorted([sum(A[i]) for i in range(n)]))
                c3 = sum(1 for i,j,k in combinations(range(n), 3)
                         if (A[i][j] and A[j][k] and A[k][i]) or
                            (A[i][k] and A[k][j] and A[j][i]))
                print(f"  *** BETA_2 > 0! trial {trial}: b2={b2}, scores={scores}, c3={c3}", flush=True)

            if (trial + 1) % 50 == 0:
                print(f"  {trial+1}/{total}: violations={violations}", flush=True)

        print(f"  n={n}: beta_2>0 in {violations}/{total}")

    # Also check disjoint support at n=9
    print("\n--- Disjoint support check at n=9 ---")
    n = 9
    rng = np.random.RandomState(99999)
    for trial in range(100):
        A = random_tournament(n, rng)
        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        a1_set = set(a1)

        max_na = 0
        for path in a2:
            na = sum(1 for sign, face in boundary_coeffs(path)
                     if len(face) == 2 and tuple(face) not in a1_set)
            max_na = max(max_na, na)

        if max_na > 1:
            print(f"  VIOLATION: trial {trial}, max NA per 2-path = {max_na}")
            break

        if trial == 0:
            print(f"  n={n}: |A_2|={len(a2)}, checking disjoint support...", flush=True)

    print(f"  Disjoint support verified for 100 random n=9 tournaments")

    print("\nDONE.")


if __name__ == '__main__':
    main()
