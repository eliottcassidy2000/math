"""
n9_exotic_profiles.py — Exotic Betti profiles at n=9.

Questions:
1. Does beta_4 become more common at n=9?
2. Does beta_6 first appear at n=9?
3. Is beta_2 still 0?
4. What are the Betti profile types and their frequencies?

Author: kind-pasteur-S45 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from itertools import combinations
from collections import Counter
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
        rk = int(sum(s > 1e-8 for s in sv))
    else: rk = 0
    ker = dim_om - rk

    omega_p1 = compute_omega_basis(A, n, target_p+1, allowed[target_p+1], allowed[target_p])
    dim_om1 = omega_p1.shape[1] if omega_p1.ndim == 2 else 0
    if dim_om1 > 0:
        bd1 = build_boundary_matrix(allowed[target_p+1], allowed[target_p])
        bd1_om = bd1 @ omega_p1
        sv1 = np.linalg.svd(bd1_om, compute_uv=False)
        im = int(sum(s > 1e-8 for s in sv1))
    else: im = 0
    return ker - im

def count_directed_3cycles(A, n):
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: c3 += 1
        if A[i][k] and A[k][j] and A[j][i]: c3 += 1
    return c3

def is_strongly_connected(A, n):
    visited_fwd = set()
    queue = [0]
    visited_fwd.add(0)
    while queue:
        v = queue.pop(0)
        for u in range(n):
            if A[v][u] and u not in visited_fwd:
                visited_fwd.add(u)
                queue.append(u)
    if len(visited_fwd) != n: return False
    visited_bwd = set()
    queue = [0]
    visited_bwd.add(0)
    while queue:
        v = queue.pop(0)
        for u in range(n):
            if A[u][v] and u not in visited_bwd:
                visited_bwd.add(u)
                queue.append(u)
    return len(visited_bwd) == n


def main():
    print("=" * 70)
    print("N=9 EXOTIC BETTI PROFILES")
    print("=" * 70)

    n = 9
    rng = np.random.RandomState(77777)
    total = 1000
    exotic = []
    profile_counts = Counter()
    beta2_violations = 0

    for trial in range(total):
        A = random_tournament(n, rng)

        bettis = [1]  # beta_0 = 1 always
        for p in range(1, n):
            b = betti_single(A, n, p)
            bettis.append(b)

        profile = tuple(bettis)
        profile_counts[profile] += 1

        if bettis[2] > 0:
            beta2_violations += 1
            print(f"  *** BETA_2 > 0 at trial {trial}! betti={bettis}", flush=True)

        nonzero = [p for p in range(1, n) if bettis[p] > 0]
        if len(nonzero) > 0:
            scores = tuple(sorted([sum(A[i]) for i in range(n)]))
            c3 = count_directed_3cycles(A, n)
            sc = is_strongly_connected(A, n)
            chi = sum((-1)**p * bettis[p] for p in range(n))
            exotic.append({
                'trial': trial, 'bettis': bettis, 'scores': scores,
                'c3': c3, 'sc': sc, 'chi': chi
            })

        if (trial + 1) % 200 == 0:
            print(f"  {trial+1}/{total} done, beta2_violations={beta2_violations}", flush=True)

    print(f"\n--- RESULTS ---")
    print(f"  Total: {total}, beta_2 violations: {beta2_violations}")
    print(f"\n  Profile frequencies:")
    for profile, cnt in profile_counts.most_common():
        nonzero = [p for p in range(n) if profile[p] > 0]
        pct = 100 * cnt / total
        if cnt >= 2 or any(profile[p] > 0 for p in [2, 4, 5, 6, 7]):
            print(f"    {list(profile)}: {cnt} ({pct:.1f}%) nonzero={nonzero}")

    # Analysis by type
    print(f"\n  Exotic type analysis:")
    type_counts = Counter()
    for e in exotic:
        nonzero = tuple(p for p in range(1, n) if e['bettis'][p] > 0)
        type_counts[nonzero] += 1

    for nonzero, cnt in type_counts.most_common():
        cases = [e for e in exotic if tuple(p for p in range(1, n) if e['bettis'][p] > 0) == nonzero]
        chi_vals = sorted(set(e['chi'] for e in cases))
        c3_range = (min(e['c3'] for e in cases), max(e['c3'] for e in cases))
        print(f"    nonzero at p={nonzero}: {cnt} cases, chi={chi_vals}, c3={c3_range}")

    # Special focus on even Betti numbers
    for p_even in [2, 4, 6]:
        cases = [e for e in exotic if e['bettis'][p_even] > 0]
        if cases:
            print(f"\n  beta_{p_even} > 0 cases ({len(cases)}):")
            for e in cases[:10]:
                print(f"    trial {e['trial']}: betti={e['bettis']}, scores={e['scores']}, c3={e['c3']}, chi={e['chi']}")
            if len(cases) > 10:
                vals = [e['bettis'][p_even] for e in cases]
                print(f"    ... {len(cases)} total, beta_{p_even} values: {sorted(Counter(vals).items())}")

    # Coexistence check
    print(f"\n  Coexistence checks:")
    for p1, p2 in [(1,3), (1,5), (3,4), (3,5), (4,5), (1,4), (4,6)]:
        coexist = sum(1 for e in exotic if e['bettis'][p1] > 0 and e['bettis'][p2] > 0)
        print(f"    beta_{p1} + beta_{p2}: {coexist} cases")


if __name__ == '__main__':
    main()
