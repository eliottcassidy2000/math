"""
Characterize n=6 tournaments with beta_3 = 1.

From exhaustive computation: exactly 320/32768 (1.0%) of n=6 tournaments have beta_3 = 1.

Questions:
1. What are the scores/structure of these tournaments?
2. Are they related to H-maximizers?
3. What are their cycle counts?
4. Does beta_1 and beta_3 mutual exclusivity hold at n=6?
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

def main():
    n = 6
    total = 2**(n*(n-1)//2)

    print("=" * 70)
    print(f"PART 1: Characterizing n={n} tournaments with beta_3 = 1")
    print("=" * 70)

    beta3_pos = []
    beta1_pos = []
    score_dist_b3 = defaultdict(int)
    H_dist_b3 = defaultdict(int)
    both_pos = 0

    for bits in range(total):
        A = bits_to_adj(bits, n)
        betti, odims = path_betti_numbers(A, n)

        b1 = betti[1]
        b3 = betti[3]

        if b1 > 0 and b3 > 0:
            both_pos += 1

        if b3 > 0:
            scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
            H = ham_path_count(A, n)
            c3 = count_directed_cycles(A, n, 3)
            c5 = count_directed_cycles(A, n, 5)
            score_dist_b3[scores] += 1
            H_dist_b3[H] += 1
            beta3_pos.append((bits, scores, H, c3, c5, list(betti), list(odims)))

        if b1 > 0:
            beta1_pos.append(bits)

    print(f"\nbeta_3 > 0: {len(beta3_pos)}/{total}")
    print(f"beta_1 > 0: {len(beta1_pos)}/{total}")
    print(f"BOTH beta_1 AND beta_3 > 0: {both_pos}/{total}")

    print(f"\nScore distribution of beta_3 = 1 tournaments:")
    for scores, count in sorted(score_dist_b3.items(), key=lambda x: -x[1]):
        print(f"  {scores}: {count}")

    print(f"\nH distribution of beta_3 = 1 tournaments:")
    for H_val, count in sorted(H_dist_b3.items()):
        print(f"  H={H_val}: {count}")

    # Show a few examples
    print(f"\nFirst 5 examples:")
    for i, (bits, scores, H, c3, c5, betti, odims) in enumerate(beta3_pos[:5]):
        print(f"  bits={bits}: scores={scores}, H={H}, c3={c3}, c5={c5}")
        print(f"    betti={betti}, omega_dims={odims}")

    # Part 2: Check if beta_3 = 1 tournaments are H-maximizers within their score class
    print("\n" + "=" * 70)
    print("PART 2: beta_3 vs max H within score class")
    print("=" * 70)

    # Build max H per score class
    max_H_by_score = defaultdict(int)
    for bits in range(total):
        A = bits_to_adj(bits, n)
        scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        H = ham_path_count(A, n)
        max_H_by_score[scores] = max(max_H_by_score[scores], H)

    for scores, count in sorted(score_dist_b3.items()):
        max_H = max_H_by_score[scores]
        h_vals = sorted(set(H for _, s, H, _, _, _, _ in beta3_pos if s == scores))
        print(f"  Score {scores}: max_H_class={max_H}, beta_3=1 H values: {h_vals}")

    # Part 3: Omega dimensions palindromicity check
    print("\n" + "=" * 70)
    print("PART 3: Omega dims palindromicity")
    print("=" * 70)

    n_palindromic = 0
    n_not = 0
    palindromic_H = defaultdict(int)

    for bits in range(total):
        A = bits_to_adj(bits, n)
        _, odims = path_betti_numbers(A, n)
        is_pal = all(odims[p] == odims[n-1-p] for p in range(n))
        if is_pal:
            n_palindromic += 1
            H = ham_path_count(A, n)
            palindromic_H[H] += 1
        else:
            n_not += 1

    print(f"Palindromic omega dims: {n_palindromic}/{total} ({100*n_palindromic/total:.1f}%)")
    print(f"Not palindromic: {n_not}/{total}")
    if palindromic_H:
        print(f"H values for palindromic: {dict(sorted(palindromic_H.items()))}")

if __name__ == '__main__':
    main()
