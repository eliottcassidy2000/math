"""
Beta_4 > 0 at n=8: detailed investigation.

THIS IS A SURPRISE. Even Betti vanishing fails at n=8 with beta_4 > 0.
Need to understand:
1. How common is beta_4 > 0?
2. What do these tournaments look like?
3. What are the full Betti vectors?
4. Does beta_2 = 0 still hold? (critical for THM-095)
5. Is chi still in {0,1}?
"""
import numpy as np
from math import comb

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

def main():
    print("=" * 70)
    print("BETA_4 > 0 AT n=8: DETAILED INVESTIGATION")
    print("=" * 70)

    n = 8

    # Part 1: Find beta_4 > 0 cases with full Betti vector
    print("\n--- Part 1: Full Betti vectors for beta_4 > 0 ---", flush=True)
    rng = np.random.RandomState(43)  # same seed as the quick check
    b4_cases = []
    all_betti_profiles = {}

    for trial in range(200):
        A = random_tournament(n, rng)

        # Compute full Betti vector
        bettis = []
        for p in range(n):
            b = betti_single(A, n, p)
            bettis.append(b)

        profile = tuple(bettis)
        all_betti_profiles[profile] = all_betti_profiles.get(profile, 0) + 1

        if bettis[4] > 0:
            scores = tuple(sorted([sum(A[i]) for i in range(n)]))
            c3 = sum(1 for i in range(n) for j in range(n) for k in range(n)
                     if i < j < k and A[i][j] and A[j][k] and A[k][i])
            c3 += sum(1 for i in range(n) for j in range(n) for k in range(n)
                      if i < j < k and A[j][i] and A[i][k] and A[k][j])
            print(f"\n  Trial {trial}: betti = {bettis}", flush=True)
            print(f"    scores = {scores}, c3 = {c3}", flush=True)
            chi = sum((-1)**p * bettis[p] for p in range(n))
            print(f"    chi = {chi}", flush=True)
            b4_cases.append((trial, bettis, scores, chi))

        if (trial+1) % 50 == 0:
            print(f"  {trial+1}/200 done", flush=True)

    print(f"\n  All Betti profiles (200 samples):", flush=True)
    for profile, cnt in sorted(all_betti_profiles.items(), key=lambda x: -x[1]):
        chi = sum((-1)**p * profile[p] for p in range(n))
        nonzero = [(p, b) for p, b in enumerate(profile) if b > 0 and p > 0]
        print(f"    {list(profile)}: {cnt} ({100*cnt/200:.1f}%), chi={chi}, nonzero={nonzero}", flush=True)

    # Part 2: Larger search for prevalence
    print(f"\n--- Part 2: Prevalence of beta_4>0 at n=8 (1000 samples) ---", flush=True)
    rng = np.random.RandomState(12345)
    b4_pos = 0
    b2_pos = 0
    b6_pos = 0
    for trial in range(1000):
        A = random_tournament(n, rng)
        b2 = betti_single(A, n, 2)
        b4 = betti_single(A, n, 4)
        b6 = betti_single(A, n, 6)
        if b2 > 0: b2_pos += 1
        if b4 > 0: b4_pos += 1
        if b6 > 0: b6_pos += 1
        if (trial+1) % 200 == 0:
            print(f"  {trial+1}/1000: b2>0={b2_pos}, b4>0={b4_pos}, b6>0={b6_pos}", flush=True)

    print(f"\n  n=8, 1000 samples:", flush=True)
    print(f"    beta_2>0: {b2_pos} ({100*b2_pos/1000:.1f}%)", flush=True)
    print(f"    beta_4>0: {b4_pos} ({100*b4_pos/1000:.1f}%)", flush=True)
    print(f"    beta_6>0: {b6_pos} ({100*b6_pos/1000:.1f}%)", flush=True)

    # Part 3: Does beta_4>0 coexist with beta_3>0 or beta_5>0?
    print(f"\n--- Part 3: Coexistence of beta_4 with odd bettis ---", flush=True)
    rng = np.random.RandomState(99999)
    coexist = {(3,4): 0, (4,5): 0, (3,5): 0, (1,4): 0}
    counts = {1: 0, 3: 0, 4: 0, 5: 0}
    for trial in range(500):
        A = random_tournament(n, rng)
        bs = {}
        for p in [1, 3, 4, 5]:
            bs[p] = betti_single(A, n, p)
            if bs[p] > 0: counts[p] = counts.get(p, 0) + 1
        for a, b in coexist:
            if bs[a] > 0 and bs[b] > 0:
                coexist[(a,b)] += 1
        if (trial+1) % 100 == 0:
            print(f"  {trial+1}/500: counts={counts}, coexist={coexist}", flush=True)

    print(f"\n  Final: counts={counts}, coexistence={coexist}", flush=True)

    # Part 4: Check if chi can be negative now
    print(f"\n--- Part 4: Chi distribution ---", flush=True)
    chi_vals = set()
    rng = np.random.RandomState(11111)
    for trial in range(500):
        A = random_tournament(n, rng)
        bettis = []
        for p in range(n):
            bettis.append(betti_single(A, n, p))
        chi = sum((-1)**p * bettis[p] for p in range(n))
        chi_vals.add(chi)
        if chi < 0 or chi > 1:
            print(f"  NEW CHI: trial {trial}, chi={chi}, betti={bettis}", flush=True)
        if (trial+1) % 100 == 0:
            print(f"  {trial+1}/500: chi values seen: {sorted(chi_vals)}", flush=True)

    print(f"\n  All chi values: {sorted(chi_vals)}", flush=True)

if __name__ == '__main__':
    main()
