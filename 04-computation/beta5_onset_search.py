"""
Beta_5 onset search: When does beta_5 first appear?

Onset pattern so far:
- beta_1 first at n=3 (d=2)
- beta_3 first at n=6 (d=5)
- beta_5 first at n=? (d=?)

If arithmetic progression with diff 3: onset at d=8 (n=9) or d=11 (n=12)?
But beta_5=0 at n=9 (from betti_rate_scaling.py). So either:
1. Onset at n=10 (d=9)?
2. Onset at n>10?

The pattern beta_{2k+1} first at n=3k could also fit:
- beta_1 (k=0): n=3*0+3 = 3 YES
- beta_3 (k=1): n=3*1+3 = 6 YES
- beta_5 (k=2): n=3*2+3 = 9? But beta_5=0 at n=9...

Alternative: beta_{2k+1} first at n=2(2k+1)+1 = 4k+3:
- beta_1: n=3 YES
- beta_3: n=7? NO (first at n=6)

Alternative: beta_p first at n=p+3:
- beta_1: n=4? NO (first at n=3)
- beta_3: n=6 YES

Hmm, let me just check computationally at n=10,11 with targeted search.
We need tournaments that maximize cyclic content to find beta_5>0.
Strategy: try near-regular tournaments, Paley-like structures.
"""
import numpy as np
from math import comb
from itertools import combinations

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
    print("BETA_5 ONSET SEARCH")
    print("=" * 70)

    # Part 1: Check Paley tournaments
    # Paley T_p has high cyclic content — best candidate for beta_5
    print("\n--- Part 1: Paley tournaments ---")

    for p in [7, 11]:
        n = p
        A = np.zeros((n, n), dtype=int)
        qr = {k*k % n for k in range(1, n)}
        for i in range(n):
            for j in range(n):
                if i != j and ((j-i) % n) in qr:
                    A[i][j] = 1

        print(f"\n  Paley T_{p}:")
        for target_p in [1, 3, 5, 7]:
            if target_p >= n: break
            b = betti_single(A, n, target_p)
            print(f"    beta_{target_p} = {b}")

    # Part 2: Random search at n=9 with many samples
    print("\n--- Part 2: Random n=9, checking beta_5 ---")
    n = 9
    rng = np.random.RandomState(42)
    N = 200
    b5_pos = 0
    b3_pos = 0
    cnt = 0

    for trial in range(N):
        A = random_tournament(n, rng)
        b5 = betti_single(A, n, 5)
        if b5 > 0:
            b5_pos += 1
            print(f"    FOUND beta_5={b5} at trial {trial}!")

        b3 = betti_single(A, n, 3)
        if b3 > 0: b3_pos += 1
        cnt += 1

        if cnt % 50 == 0:
            print(f"    {cnt}/{N}: beta_3>0 = {100*b3_pos/cnt:.1f}%, beta_5>0 = {100*b5_pos/cnt:.1f}%")

    print(f"  n=9 result: beta_5>0 = {b5_pos}/{cnt} ({100*b5_pos/cnt:.1f}%)")
    print(f"              beta_3>0 = {b3_pos}/{cnt} ({100*b3_pos/cnt:.1f}%)")

    # Part 3: Exhaustive search for beta_5 among regular n=9 tournaments
    # Regular n=9: scores all 4. These have maximal cyclic content.
    print("\n--- Part 3: Regular n=9 (score 4) tournaments ---")
    n = 9
    rng = np.random.RandomState(100)
    N = 100  # Sample regular tournaments
    b5_pos = 0
    b5_max = 0
    b3_pos = 0
    cnt = 0

    for trial in range(N):
        # Generate near-regular tournament by repeated random flips
        # Start with a random tournament and reject if not regular
        for _ in range(10000):
            A = random_tournament(n, rng)
            scores = [sum(A[i]) for i in range(n)]
            if all(s == 4 for s in scores):
                break
        else:
            # Very hard to hit exactly regular by random sampling
            # Use a different approach: start from Paley and flip
            continue

        scores = [sum(A[i]) for i in range(n)]
        if not all(s == 4 for s in scores):
            continue

        b5 = betti_single(A, n, 5)
        if b5 > 0:
            b5_pos += 1
            b5_max = max(b5_max, b5)
            print(f"    FOUND beta_5={b5} at trial {trial}!")

        b3 = betti_single(A, n, 3)
        if b3 > 0: b3_pos += 1
        cnt += 1

    if cnt > 0:
        print(f"  Regular n=9 result: {cnt} samples")
        print(f"    beta_5>0 = {b5_pos}/{cnt}, max = {b5_max}")
        print(f"    beta_3>0 = {b3_pos}/{cnt}")
    else:
        print(f"  Could not generate enough regular tournaments by random sampling")

    # Part 4: Check beta_5 on Paley T_11 specifically
    print("\n--- Part 4: Paley T_11 detailed ---")
    n = 11
    A = np.zeros((n, n), dtype=int)
    qr = {k*k % n for k in range(1, n)}
    for i in range(n):
        for j in range(n):
            if i != j and ((j-i) % n) in qr:
                A[i][j] = 1

    print(f"  Paley T_11 QR = {qr}")
    for target_p in [1, 3, 5, 7, 9]:
        if target_p >= n:
            break
        print(f"    Computing beta_{target_p}...", end=" ", flush=True)
        b = betti_single(A, n, target_p)
        print(f"beta_{target_p} = {b}")

if __name__ == '__main__':
    main()
