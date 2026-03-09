"""
n9_fast_betti.py — Fast check of key Betti numbers at n=9.

Strategy: Only compute beta_2 and beta_4 (the key even ones).
Skip beta_1, beta_3 etc. unless needed.
Use path enumeration caching.

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

def compute_omega_dim_fast(A, n, p, allowed_p, allowed_pm1):
    """Return dim(Omega_p) without computing the full basis."""
    dim_Ap = len(allowed_p)
    if dim_Ap == 0: return 0
    if p == 0: return dim_Ap
    allowed_pm1_set = set(allowed_pm1)
    non_allowed = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1
    if na_count == 0: return dim_Ap
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed:
                P[non_allowed[face], j] += sign
    sv = np.linalg.svd(P, compute_uv=False)
    rank = int(sum(s > 1e-10 for s in sv))
    return dim_Ap - rank

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

def betti_single(A, n, target_p, allowed_cache=None):
    """Compute beta_{target_p} using cached path enumerations if available."""
    if allowed_cache is None:
        allowed_cache = {}

    for p in [target_p - 1, target_p, target_p + 1]:
        if p not in allowed_cache:
            if p < 0:
                allowed_cache[p] = []
            else:
                allowed_cache[p] = enumerate_allowed_paths(A, n, p)

    omega_p = compute_omega_basis(A, n, target_p, allowed_cache[target_p], allowed_cache[target_p-1])
    dim_om = omega_p.shape[1] if omega_p.ndim == 2 else 0
    if dim_om == 0: return 0

    bd_p = build_boundary_matrix(allowed_cache[target_p], allowed_cache[target_p-1])
    bd_p_om = bd_p @ omega_p
    if bd_p_om.size > 0:
        sv = np.linalg.svd(bd_p_om, compute_uv=False)
        rk = int(sum(s > 1e-8 for s in sv))
    else: rk = 0
    ker = dim_om - rk

    omega_p1 = compute_omega_basis(A, n, target_p+1, allowed_cache[target_p+1], allowed_cache[target_p])
    dim_om1 = omega_p1.shape[1] if omega_p1.ndim == 2 else 0
    if dim_om1 > 0:
        bd1 = build_boundary_matrix(allowed_cache[target_p+1], allowed_cache[target_p])
        bd1_om = bd1 @ omega_p1
        sv1 = np.linalg.svd(bd1_om, compute_uv=False)
        im = int(sum(s > 1e-8 for s in sv1))
    else: im = 0
    return ker - im


def main():
    print("=" * 70)
    print("FAST BETTI CHECK AT N=9")
    print("=" * 70)

    n = 9
    total = 300
    rng = np.random.RandomState(54321)

    beta2_violations = 0
    beta4_count = 0
    beta6_count = 0
    exotic_profiles = []

    for trial in range(total):
        A = random_tournament(n, rng)

        # Cache path enumerations
        allowed_cache = {}
        for p in range(-1, n+1):
            if p < 0:
                allowed_cache[p] = []
            elif p <= 5:
                # Only enumerate up to p=5 (for beta_4 we need p=3,4,5)
                allowed_cache[p] = enumerate_allowed_paths(A, n, p)
                if trial == 0:
                    print(f"  p={p}: |A_p|={len(allowed_cache[p])}", flush=True)
            else:
                allowed_cache[p] = enumerate_allowed_paths(A, n, p)
                if trial == 0:
                    print(f"  p={p}: |A_p|={len(allowed_cache[p])}", flush=True)

        # Compute beta_2 (needs p=1,2,3)
        b2 = betti_single(A, n, 2, allowed_cache)
        if b2 > 0:
            beta2_violations += 1
            print(f"  *** BETA_2 > 0 at trial {trial}! b2={b2}", flush=True)

        # Compute beta_4 (needs p=3,4,5)
        b4 = betti_single(A, n, 4, allowed_cache)
        if b4 > 0:
            beta4_count += 1

        # Compute beta_6 only if we can (needs p=5,6,7 — might be expensive)
        # Actually at n=9, p=7 paths have 8 vertices = very few
        b6 = betti_single(A, n, 6, allowed_cache)
        if b6 > 0:
            beta6_count += 1

        # If anything exotic, compute full profile
        if b2 > 0 or b4 > 0 or b6 > 0:
            bettis = [1]
            for p in range(1, n):
                if p == 2:
                    bettis.append(b2)
                elif p == 4:
                    bettis.append(b4)
                elif p == 6:
                    bettis.append(b6)
                else:
                    bettis.append(betti_single(A, n, p, allowed_cache))
            scores = tuple(sorted([sum(A[i]) for i in range(n)]))
            c3 = sum(1 for i,j,k in combinations(range(n), 3)
                     if (A[i][j] and A[j][k] and A[k][i]) or
                        (A[i][k] and A[k][j] and A[j][i]))
            chi = sum((-1)**p * bettis[p] for p in range(n))
            exotic_profiles.append({
                'trial': trial, 'bettis': bettis, 'scores': scores, 'c3': c3, 'chi': chi
            })
            print(f"  trial {trial}: betti={bettis}, scores={scores}, c3={c3}, chi={chi}", flush=True)

        if (trial + 1) % 50 == 0:
            print(f"  {trial+1}/{total}: b2_viol={beta2_violations}, b4={beta4_count}, b6={beta6_count}", flush=True)

    print(f"\n--- SUMMARY ---")
    print(f"  n=9, {total} samples")
    print(f"  beta_2 > 0: {beta2_violations}")
    print(f"  beta_4 > 0: {beta4_count} ({100*beta4_count/total:.1f}%)")
    print(f"  beta_6 > 0: {beta6_count} ({100*beta6_count/total:.1f}%)")

    if exotic_profiles:
        print(f"\n  Exotic profiles:")
        type_counts = Counter()
        for e in exotic_profiles:
            nonzero = tuple(p for p in range(1, n) if e['bettis'][p] > 0)
            type_counts[nonzero] += 1
        for nz, cnt in type_counts.most_common():
            print(f"    nonzero at p={nz}: {cnt} cases")

    print("\nDONE.")


if __name__ == '__main__':
    main()
