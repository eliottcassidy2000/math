"""Quick even Betti check at n=7,8 and beta_3*beta_5 test."""
import numpy as np
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

def betti_single(A, n, target_p):
    allowed = {}
    for p in [target_p - 1, target_p, target_p + 1]:
        if p < 0: allowed[p] = []
        else: allowed[p] = enumerate_allowed_paths(A, n, p)
    omega_p = compute_omega_basis(A, n, target_p, allowed[target_p], allowed[target_p-1])
    dim_om = omega_p.shape[1] if omega_p.ndim == 2 else 0
    if dim_om == 0: return 0
    bd_p = build_boundary_matrix(allowed[target_p], allowed[target_p-1])
    bd_p_om = bd_p @ omega_p
    if bd_p_om.size > 0:
        sv = np.linalg.svd(bd_p_om, compute_uv=False)
        rk = sum(s > 1e-8 for s in sv)
    else: rk = 0
    ker = dim_om - rk
    omega_p1 = compute_omega_basis(A, n, target_p+1, allowed[target_p+1], allowed[target_p])
    dim_om1 = omega_p1.shape[1] if omega_p1.ndim == 2 else 0
    if dim_om1 > 0:
        bd1 = build_boundary_matrix(allowed[target_p+1], allowed[target_p])
        bd1_om = bd1 @ omega_p1
        sv1 = np.linalg.svd(bd1_om, compute_uv=False)
        im = sum(s > 1e-8 for s in sv1)
    else: im = 0
    return ker - im

print("EVEN BETTI VANISHING + BETA_3*BETA_5 CHECK", flush=True)

# Part 1: Even betti at n=7
print("\n--- n=7: beta_2, beta_4 (500 samples) ---", flush=True)
rng = np.random.RandomState(42)
violations = 0
for trial in range(500):
    A = random_tournament(7, rng)
    for p in [2, 4]:
        b = betti_single(A, 7, p)
        if b != 0:
            violations += 1
            print(f"  VIOLATION: trial {trial}, beta_{p} = {b}", flush=True)
    if (trial+1) % 100 == 0:
        print(f"  {trial+1}/500, violations: {violations}", flush=True)
print(f"  Result: {violations} violations", flush=True)

# Part 2: Even betti at n=8
print("\n--- n=8: beta_2, beta_4, beta_6 (200 samples) ---", flush=True)
rng = np.random.RandomState(43)
violations = 0
for trial in range(200):
    A = random_tournament(8, rng)
    for p in [2, 4, 6]:
        b = betti_single(A, 8, p)
        if b != 0:
            violations += 1
            print(f"  VIOLATION: trial {trial}, beta_{p} = {b}", flush=True)
    if (trial+1) % 50 == 0:
        print(f"  {trial+1}/200, violations: {violations}", flush=True)
print(f"  Result: {violations} violations", flush=True)

# Part 3: beta_3*beta_5 at n=8
print("\n--- n=8: beta_3*beta_5 mutual exclusivity (300 samples) ---", flush=True)
rng = np.random.RandomState(777)
b3_cnt, b5_cnt, both_cnt = 0, 0, 0
for trial in range(300):
    A = random_tournament(8, rng)
    b3 = betti_single(A, 8, 3)
    b5 = betti_single(A, 8, 5)
    if b3 > 0: b3_cnt += 1
    if b5 > 0: b5_cnt += 1
    if b3 > 0 and b5 > 0:
        both_cnt += 1
        print(f"  BOTH NONZERO: trial {trial}, beta_3={b3}, beta_5={b5}", flush=True)
    if (trial+1) % 100 == 0:
        print(f"  {trial+1}/300: b3>0={b3_cnt}, b5>0={b5_cnt}, both={both_cnt}", flush=True)

print(f"\n  RESULT: beta_3>0={b3_cnt}, beta_5>0={b5_cnt}, both>0={both_cnt}", flush=True)
if both_cnt == 0:
    print("  beta_3*beta_5 = 0 CONFIRMED (extends seesaw to next pair)", flush=True)
