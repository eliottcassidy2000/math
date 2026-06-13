"""
Does beta_3 > 0 correlate with H-maximization within score class at n=7?

At n=6: beta_3>0 iff H = max within score class (PROVED exhaustive).
At n=7: beta_3>0 occurs for ~7% of tournaments. Is it score-class max?

Also: what characterizes beta_1>0 vs beta_3>0?
Hypothesis: beta_1>0 = "low H" within score class, beta_3>0 = "high H".
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
    n = 7
    rng = np.random.RandomState(42)

    print("=" * 70)
    print(f"beta_3 and beta_1 vs H-maximization at n={n}")
    print("=" * 70)

    # Collect data: score, H, beta_1, beta_3
    data_by_score = defaultdict(list)
    count = 0

    for trial in range(5000):
        A = random_tournament(n, rng)
        scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        H = ham_path_count(A, n)
        b1 = betti_single(A, n, 1)
        b3 = betti_single(A, n, 3)
        data_by_score[scores].append((H, b1, b3))
        count += 1

        if trial % 500 == 499:
            print(f"  ... {trial+1} trials done")

    print(f"\nTotal: {count} tournaments across {len(data_by_score)} score classes")

    # For each score class: what H values have beta_3>0? What H has beta_1>0?
    print("\n" + "=" * 70)
    print("Score class analysis: beta_3>0 vs H rank")
    print("=" * 70)

    for scores in sorted(data_by_score.keys()):
        entries = data_by_score[scores]
        n_entries = len(entries)
        if n_entries < 3:
            continue

        max_H = max(e[0] for e in entries)
        min_H = min(e[0] for e in entries)

        b3_pos = [(H, b1, b3) for H, b1, b3 in entries if b3 > 0]
        b1_pos = [(H, b1, b3) for H, b1, b3 in entries if b1 > 0]

        if not b3_pos and not b1_pos:
            continue  # Skip all-trivial score classes

        b3_H = sorted(set(H for H, _, _ in b3_pos)) if b3_pos else []
        b1_H = sorted(set(H for H, _, _ in b1_pos)) if b1_pos else []

        all_H = sorted(set(H for H, _, _ in entries))

        print(f"\n  {scores}: {n_entries} samples, H in {all_H[:3]}...{all_H[-3:] if len(all_H)>3 else ''}")
        print(f"    max_H={max_H}, min_H={min_H}")
        if b3_pos:
            print(f"    beta_3>0: {len(b3_pos)}/{n_entries} at H={b3_H}")
            is_max = all(H == max_H for H in b3_H)
            print(f"    beta_3>0 only at max_H? {'YES' if is_max else 'NO'}")
        if b1_pos:
            print(f"    beta_1>0: {len(b1_pos)}/{n_entries} at H={b1_H}")

    # Global summary
    all_b3 = []
    all_b1 = []
    all_zero = []
    for scores, entries in data_by_score.items():
        max_H = max(e[0] for e in entries)
        for H, b1, b3 in entries:
            if b3 > 0:
                all_b3.append((scores, H, H == max_H))
            elif b1 > 0:
                all_b1.append((scores, H, H == max_H))
            else:
                all_zero.append((scores, H, H == max_H))

    print("\n" + "=" * 70)
    print("Global summary")
    print("=" * 70)
    if all_b3:
        n_at_max = sum(1 for _, _, is_max in all_b3 if is_max)
        print(f"beta_3>0: {len(all_b3)} total, {n_at_max}/{len(all_b3)} at score-class max H"
              f" ({100*n_at_max/len(all_b3):.1f}%)")
    if all_b1:
        n_at_max = sum(1 for _, _, is_max in all_b1 if is_max)
        print(f"beta_1>0: {len(all_b1)} total, {n_at_max}/{len(all_b1)} at score-class max H"
              f" ({100*n_at_max/len(all_b1):.1f}%)")

if __name__ == '__main__':
    main()
