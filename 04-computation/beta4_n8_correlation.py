"""
Check whether beta_4 > 0 correlates with H=661 at n=8.

Finding from beta5_n8_multi.py:
- beta_5 = 0 for ALL 20 H=661 tournaments (pattern breaks)
- But tournament #3 had beta_4 = 1 with H=661
- 2/3 H=661 had beta_4=0, 1/3 had beta_4=1

Question: Is beta_4 > 0 exclusive to H-maximizers at n=8?
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
    """Compute just beta_{target_p}."""
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
    target_score = (3,3,3,3,4,4,4,4)

    print("=" * 70)
    print(f"beta_4 vs H correlation for n={n} SC-score tournaments")
    print("=" * 70)

    rng = np.random.RandomState(789)

    # Collect SC-score tournaments, compute H and beta_4
    beta4_by_H = defaultdict(list)
    count = 0
    target = 200  # aim for 200 SC-score tournaments

    trials = 0
    while count < target and trials < 500000:
        A = random_tournament(n, rng)
        trials += 1
        scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        if scores != target_score:
            continue

        count += 1
        H = ham_path_count(A, n)
        b4 = betti_single(A, n, 4)
        beta4_by_H[H].append(b4)

        if count % 25 == 0:
            print(f"  ... {count} SC-score tournaments processed ({trials} total trials)")

    print(f"\nTotal: {count} SC-score tournaments from {trials} trials")
    print(f"\nbeta_4 by H value:")
    for H in sorted(beta4_by_H.keys()):
        vals = beta4_by_H[H]
        n_pos = sum(1 for v in vals if v > 0)
        n_tot = len(vals)
        max_b4 = max(vals) if vals else 0
        print(f"  H={H}: {n_tot} tournaments, beta_4>0 in {n_pos}/{n_tot}"
              f" ({100*n_pos/n_tot:.1f}%), max_beta_4={max_b4}")

    # Summary
    all_h661 = beta4_by_H.get(661, [])
    all_non661 = []
    for H, vals in beta4_by_H.items():
        if H != 661:
            all_non661.extend(vals)

    if all_h661:
        n_pos_661 = sum(1 for v in all_h661 if v > 0)
        print(f"\nH=661: beta_4>0 in {n_pos_661}/{len(all_h661)}")
    if all_non661:
        n_pos_other = sum(1 for v in all_non661 if v > 0)
        print(f"H!=661: beta_4>0 in {n_pos_other}/{len(all_non661)}")

if __name__ == '__main__':
    main()
