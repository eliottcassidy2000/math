"""
Full Betti vector analysis at n=7.

Questions:
1. What Betti vectors occur? Distribution?
2. Are beta_1 and beta_3 mutually exclusive (as at n=6)?
3. Are beta_1 and beta_4 mutually exclusive?
4. Chi distribution: {0,1} or larger values?
5. Does beta_3 ever occur at n=7?
6. For non-regular tournaments: what Betti do H-maximizers have?
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
    n = 7
    rng = np.random.RandomState(42)

    print("=" * 70)
    print(f"Full Betti vector analysis at n={n} (sampling)")
    print("=" * 70)

    betti_dist = defaultdict(int)
    chi_dist = defaultdict(int)
    betti_by_score = defaultdict(lambda: defaultdict(int))
    max_H_by_score = defaultdict(int)
    betti_of_hmax = defaultdict(lambda: defaultdict(int))
    total = 0

    # Track specific combinations
    b1_pos = 0
    b3_pos = 0
    b4_pos = 0
    b1_and_b3 = 0
    b1_and_b4 = 0
    b3_and_b4 = 0

    for trial in range(2000):
        A = random_tournament(n, rng)
        H = ham_path_count(A, n)
        scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))

        # Only compute full Betti for first 500
        if trial < 500:
            betti, odims = path_betti_numbers(A, n)
            bv = tuple(betti)
            betti_dist[bv] += 1
            chi = sum((-1)**p * betti[p] for p in range(len(betti)))
            chi_dist[chi] += 1
            betti_by_score[scores][bv] += 1

            if betti[1] > 0: b1_pos += 1
            if len(betti) > 3 and betti[3] > 0: b3_pos += 1
            if len(betti) > 4 and betti[4] > 0: b4_pos += 1
            if betti[1] > 0 and len(betti) > 3 and betti[3] > 0: b1_and_b3 += 1
            if betti[1] > 0 and len(betti) > 4 and betti[4] > 0: b1_and_b4 += 1
            if len(betti) > 3 and betti[3] > 0 and len(betti) > 4 and betti[4] > 0: b3_and_b4 += 1
            total += 1

        max_H_by_score[scores] = max(max_H_by_score[scores], H)

        if trial % 200 == 199:
            print(f"  ... {trial+1} trials done")

    print(f"\nTotal with full Betti: {total}")
    print(f"\nBetti vector distribution:")
    for bv in sorted(betti_dist.keys(), key=lambda x: -betti_dist[x]):
        count = betti_dist[bv]
        pct = 100*count/total
        chi = sum((-1)**p * bv[p] for p in range(len(bv)))
        print(f"  {list(bv)}: {count} ({pct:.1f}%), chi={chi}")

    print(f"\nChi distribution:")
    for chi in sorted(chi_dist.keys()):
        print(f"  chi={chi}: {chi_dist[chi]} ({100*chi_dist[chi]/total:.1f}%)")

    print(f"\nMutual exclusivity:")
    print(f"  beta_1 > 0: {b1_pos}/{total} ({100*b1_pos/total:.1f}%)")
    print(f"  beta_3 > 0: {b3_pos}/{total} ({100*b3_pos/total:.1f}%)")
    print(f"  beta_4 > 0: {b4_pos}/{total} ({100*b4_pos/total:.1f}%)")
    print(f"  beta_1 AND beta_3 > 0: {b1_and_b3}")
    print(f"  beta_1 AND beta_4 > 0: {b1_and_b4}")
    print(f"  beta_3 AND beta_4 > 0: {b3_and_b4}")

    # For the regular score class, use targeted search
    print("\n" + "=" * 70)
    print("Regular tournament (3,3,3,3,3,3,3) Betti analysis")
    print("=" * 70)

    reg_betti_dist = defaultdict(int)
    reg_count = 0
    rng2 = np.random.RandomState(100)
    for trial in range(50000):
        A = random_tournament(n, rng2)
        scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        if scores != (3,3,3,3,3,3,3):
            continue
        H = ham_path_count(A, n)
        betti, odims = path_betti_numbers(A, n)
        bv = tuple(betti)
        reg_betti_dist[(H, bv)] += 1
        reg_count += 1
        if reg_count % 10 == 0:
            print(f"  ... {reg_count} regular tournaments processed")
        if reg_count >= 30:
            break

    print(f"\nRegular tournament Betti by H:")
    for (H, bv), count in sorted(reg_betti_dist.items()):
        chi = sum((-1)**p * bv[p] for p in range(len(bv)))
        print(f"  H={H}: betti={list(bv)}, chi={chi}, count={count}")

if __name__ == '__main__':
    main()
