"""
Even Betti number vanishing: beta_{2k} = 0 for all tournaments?

This is a KEY structural claim underlying THM-095 (seesaw proof).
If all even Betti numbers vanish, adjacent odd Betti numbers are
mutually exclusive through the seesaw mechanism.

Check: beta_0, beta_2, beta_4 for all tournaments at various n.

beta_0 = 1 for connected digraphs. Every tournament has a Hamiltonian
path (Redei), so the undirected underlying graph is connected.
Actually for GLMY path homology, beta_0 = # weakly connected components.
Tournaments are always weakly connected. So beta_0 = 1 always.

beta_2 = ? (main question)
beta_4 = ? (at n >= 7)
"""

import numpy as np
from math import comb

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
    print("=" * 70)
    print("EVEN BETTI NUMBER VANISHING CHECK")
    print("=" * 70)

    for n in range(3, 10):
        print(f"\n--- n={n} ---")
        if n <= 6:
            N = 2**(n*(n-1)//2)
            exhaustive = True
        elif n == 7:
            N = 2000
            exhaustive = False
        elif n == 8:
            N = 500
            exhaustive = False
        else:
            N = 200
            exhaustive = False

        rng = np.random.RandomState(42 + n)

        even_bettis = {}  # p -> list of nonzero values
        max_even_p = min(n-1, 6)  # check even p up to 6

        cnt = 0
        for trial in range(N):
            if exhaustive:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            for p in range(0, max_even_p + 1, 2):
                if p >= n:
                    break
                b = betti_single(A, n, p)
                if b != 0:
                    if p not in even_bettis:
                        even_bettis[p] = []
                    even_bettis[p].append((trial, b))

            cnt += 1
            if not exhaustive and cnt % 500 == 0:
                print(f"  {cnt}/{N} checked")

        if not even_bettis:
            print(f"  ALL even Betti numbers = 0 ({cnt} {'exhaustive' if exhaustive else 'sampled'})")
        else:
            for p in sorted(even_bettis.keys()):
                examples = even_bettis[p]
                print(f"  beta_{p} NONZERO: {len(examples)} cases")
                for trial, b in examples[:3]:
                    print(f"    trial {trial}: beta_{p} = {b}")

        # Also check odd Betti numbers for reference
        print(f"  Checking odd Betti numbers...")
        odd_counts = {}
        for trial in range(min(N, 500 if n >= 8 else N)):
            if exhaustive:
                A = bits_to_adj(trial, n)
            else:
                rng2 = np.random.RandomState(42 + n)
                for _ in range(trial + 1):
                    A = random_tournament(n, rng2)

            for p in range(1, min(n, 8), 2):
                b = betti_single(A, n, p)
                if b > 0:
                    odd_counts[p] = odd_counts.get(p, 0) + 1

        for p in sorted(odd_counts.keys()):
            print(f"    beta_{p} > 0: {odd_counts[p]} times")

    print("\n" + "=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    print("If all even Betti numbers vanish, then adjacent odd Betti numbers")
    print("are mutually exclusive via the seesaw mechanism (THM-095).")

if __name__ == '__main__':
    main()
