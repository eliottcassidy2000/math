"""
Euler characteristic analysis for GLMY path homology of tournaments.

Key finding: Paley T_7 has chi = 1 + 0 - 0 + 0 + 6 - 0 + 0 = 7 = n.

Questions:
1. Is chi = 1 always? (Euler char of a "contractible" space)
2. Is chi = n for Paley/VT?
3. Distribution of chi across all tournaments?
4. chi = sum(-1)^p * dim(Omega_p)? (Euler char of chain complex)
"""
import numpy as np
from itertools import combinations
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

def main():
    # Part 1: Euler characteristic distribution at n=5
    print("=" * 70)
    print("PART 1: Euler characteristic at n=5 (exhaustive)")
    print("=" * 70)
    n = 5
    total = 2**(n*(n-1)//2)
    chi_dist = defaultdict(int)
    chi_omega_dist = defaultdict(int)  # chi from Omega dims

    for bits in range(total):
        A = bits_to_adj(bits, n)
        betti, odims = path_betti_numbers(A, n)
        chi = sum((-1)**p * betti[p] for p in range(n))
        chi_omega = sum((-1)**p * odims[p] for p in range(n))
        chi_dist[chi] += 1
        chi_omega_dist[chi_omega] += 1

    print(f"n={n}: chi (from Betti) distribution: {dict(sorted(chi_dist.items()))}")
    print(f"n={n}: chi (from Omega dims) distribution: {dict(sorted(chi_omega_dist.items()))}")

    # Part 2: Euler characteristic at n=6
    print("\n" + "=" * 70)
    print("PART 2: Euler characteristic at n=6 (exhaustive)")
    print("=" * 70)
    n = 6
    total = 2**(n*(n-1)//2)
    chi_dist = defaultdict(int)
    chi_omega_dist = defaultdict(int)

    for bits in range(total):
        A = bits_to_adj(bits, n)
        betti, odims = path_betti_numbers(A, n)
        chi = sum((-1)**p * betti[p] for p in range(n))
        chi_omega = sum((-1)**p * odims[p] for p in range(n))
        chi_dist[chi] += 1
        chi_omega_dist[chi_omega] += 1

    print(f"n={n}: chi (Betti) distribution: {dict(sorted(chi_dist.items()))}")
    print(f"n={n}: chi (Omega) distribution: {dict(sorted(chi_omega_dist.items()))}")

    # Part 3: Omega Euler characteristic = what?
    print("\n" + "=" * 70)
    print("PART 3: Omega Euler char interpretation")
    print("=" * 70)
    print("For exact sequences: chi(chain complex) = alternating sum of dims")
    print("Since d^2 = 0, chi(Omega) = sum(-1)^p dim(Omega_p) = sum(-1)^p beta_p")
    print("Wait — they should be equal! That's the rank-nullity theorem.")
    print("Let me verify...")

    n = 5
    total = 2**(n*(n-1)//2)
    all_equal = True
    for bits in range(total):
        A = bits_to_adj(bits, n)
        betti, odims = path_betti_numbers(A, n)
        chi_b = sum((-1)**p * betti[p] for p in range(n))
        chi_o = sum((-1)**p * odims[p] for p in range(n))
        if chi_b != chi_o:
            all_equal = False
            print(f"  MISMATCH at bits={bits}: chi_betti={chi_b}, chi_omega={chi_o}")
            break
    if all_equal:
        print(f"  n=5: chi(Betti) = chi(Omega) for ALL {total} tournaments (as expected)")

    n = 6
    total = 2**(n*(n-1)//2)
    all_equal = True
    for bits in range(total):
        A = bits_to_adj(bits, n)
        betti, odims = path_betti_numbers(A, n)
        chi_b = sum((-1)**p * betti[p] for p in range(n))
        chi_o = sum((-1)**p * odims[p] for p in range(n))
        if chi_b != chi_o:
            all_equal = False
            print(f"  MISMATCH at bits={bits}: chi_betti={chi_b}, chi_omega={chi_o}")
            break
    if all_equal:
        print(f"  n=6: chi(Betti) = chi(Omega) for ALL {total} tournaments (as expected)")

    # Part 4: Can chi be expressed in terms of tournament invariants?
    print("\n" + "=" * 70)
    print("PART 4: Chi in terms of n, H, c3, etc.")
    print("=" * 70)

    for n in [3, 4, 5]:
        total = 2**(n*(n-1)//2)
        print(f"\nn={n}:")
        for bits in range(total):
            A = bits_to_adj(bits, n)
            betti, odims = path_betti_numbers(A, n)
            chi = sum((-1)**p * betti[p] for p in range(n))
            # Count directed 3-cycles
            c3 = 0
            for a, b, c in combinations(range(n), 3):
                if A[a][b] and A[b][c] and A[c][a]: c3 += 1
                if A[a][c] and A[c][b] and A[b][a]: c3 += 1
            scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
            print(f"  scores={scores}, c3={c3}, chi={chi}, betti={list(betti)}")

    # At n=5: chi = 1 - beta_1 (since beta_2=beta_3=beta_4=0)
    # So chi in {0, 1}: chi=1 when beta_1=0, chi=0 when beta_1=1
    print("\nKey observation: chi = 1 - beta_1 at n=5 (since higher betti all 0)")
    print("At n=6: chi = 1 - beta_1 + beta_3 (since beta_2=beta_4=beta_5=0)")

if __name__ == '__main__':
    main()
