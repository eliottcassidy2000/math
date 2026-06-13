#!/usr/bin/env python3
"""
consecutive_seesaw.py — opus-2026-03-09-S54

HYPOTHESIS: beta_k * beta_{k+1} = 0 for all tournaments and all k >= 1.
("Consecutive seesaw" — adjacent Betti numbers can't both be nonzero.)

This would extend the known:
  - beta_2 = 0 always (THM-108)
  - beta_{2k-1} * beta_{2k+1} = 0 (adjacent odd seesaw)

If true, combined with beta_2 = 0, the only possible patterns are:
  - Single nonzero odd Betti: (...,0,b_{2k+1},0,...) — the seesaw pattern
  - Single nonzero even Betti: (...,0,b_{2k},0,...) — e.g., Paley T_7 with b_4=6

RELEVANCE TO PROOF: If beta_3(T) * beta_4(T) = 0, then when beta_3(T) = 1,
we have beta_4(T) = 0. This means H_4(T) = 0, which constrains the LES
connecting map δ: H_4(T,T\v) → H_3(T\v).
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter

def tournament_from_bits(n, bits):
    T = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                T[i][j] = 1
            else:
                T[j][i] = 1
            idx += 1
    return T

def is_allowed_path(T, path):
    for i in range(len(path)-1):
        if not T[path[i]][path[i+1]]:
            return False
    return len(path) == len(set(path))

def compute_betti(T, n, max_p_limit=None):
    tol = 1e-8
    if max_p_limit is None:
        max_p_limit = n - 1

    all_paths = {}
    for p in range(min(n, max_p_limit + 1)):
        paths = []
        for verts in combinations(range(n), p+1):
            for perm in permutations(verts):
                if is_allowed_path(T, perm):
                    paths.append(perm)
        all_paths[p] = paths

    max_p = max(all_paths.keys())

    # Build Omega
    omega = {}
    for p in range(max_p + 1):
        a_p = all_paths.get(p, [])
        if not a_p:
            omega[p] = np.zeros((0, 0))
            continue
        a_pm1_set = set(all_paths.get(p-1, [])) if p > 0 else set()
        na = {}
        for sigma in a_p:
            for i in range(1, len(sigma)-1):
                face = sigma[:i] + sigma[i+1:]
                if p > 0 and face not in a_pm1_set:
                    if face not in na:
                        na[face] = len(na)
        if not na:
            omega[p] = np.eye(len(a_p))
        else:
            mat = np.zeros((len(na), len(a_p)))
            for j, sigma in enumerate(a_p):
                for i in range(1, len(sigma)-1):
                    face = sigma[:i] + sigma[i+1:]
                    if face in na:
                        mat[na[face], j] += (-1)**i
            try:
                U, S, Vt = np.linalg.svd(mat, full_matrices=True)
                rank = int(np.sum(S > tol))
                if len(a_p) - rank == 0:
                    omega[p] = np.zeros((len(a_p), 0))
                else:
                    omega[p] = Vt[rank:].T
            except:
                omega[p] = np.eye(len(a_p))

    # Boundary maps
    boundary = {}
    for p in range(1, max_p + 1):
        a_p = all_paths.get(p, [])
        a_pm1 = all_paths.get(p-1, [])
        if not a_p or not a_pm1:
            boundary[p] = np.zeros((0, 0))
            continue
        idx = {path: i for i, path in enumerate(a_pm1)}
        mat = np.zeros((len(a_pm1), len(a_p)))
        for j, sigma in enumerate(a_p):
            for k in range(len(sigma)):
                face = sigma[:k] + sigma[k+1:]
                if face in idx:
                    mat[idx[face], j] += (-1)**k
        boundary[p] = mat

    # Projected boundary ranks
    dims = {p: omega[p].shape[1] if omega[p].ndim == 2 else 0 for p in range(max_p+1)}
    ranks = {}
    for p in range(1, max_p + 1):
        Om_p = omega.get(p, np.zeros((0,0)))
        Om_pm1 = omega.get(p-1, np.zeros((0,0)))
        if Om_p.ndim < 2 or Om_p.shape[1] == 0 or Om_pm1.ndim < 2 or Om_pm1.shape[1] == 0:
            ranks[p] = 0
            continue
        dp = Om_pm1.T @ boundary[p] @ Om_p
        try:
            S = np.linalg.svd(dp, compute_uv=False)
            ranks[p] = int(np.sum(S > tol))
        except:
            ranks[p] = 0

    betti = {}
    for p in range(max_p + 1):
        ker = dims[p] - ranks.get(p, 0)
        im_next = ranks.get(p+1, 0)
        betti[p] = ker - im_next

    return betti

def main():
    print("=" * 70)
    print("CONSECUTIVE SEESAW: beta_k * beta_{k+1} = 0 ?")
    print("=" * 70)

    # Part 1: Exhaustive at n=6
    n = 6
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs

    violations_6 = []
    betti_profiles_6 = Counter()

    print(f"\nPart 1: Exhaustive check at n={n}")
    for bits in range(n_total):
        if bits % 5000 == 0:
            print(f"  ... {bits}/{n_total}", flush=True)
        T = tournament_from_bits(n, bits)
        try:
            betti = compute_betti(T, n)
        except:
            continue

        profile = tuple(betti.get(p, 0) for p in range(n))
        betti_profiles_6[profile] += 1

        # Check consecutive seesaw
        for k in range(1, n-1):
            if betti.get(k, 0) > 0 and betti.get(k+1, 0) > 0:
                violations_6.append((bits, k, betti.get(k,0), betti.get(k+1,0)))

    print(f"\n  Betti profiles at n=6:")
    for profile, cnt in betti_profiles_6.most_common():
        print(f"    {list(profile)}: {cnt}")

    print(f"\n  Consecutive seesaw violations at n=6: {len(violations_6)}")
    if violations_6:
        for bits, k, bk, bk1 in violations_6[:5]:
            print(f"    bits={bits}: beta_{k}={bk}, beta_{k+1}={bk1}")

    # Part 2: Sampled at n=7
    n = 7
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs
    rng = np.random.RandomState(42)

    violations_7 = []
    betti_profiles_7 = Counter()
    n_samples = 3000

    print(f"\nPart 2: Sampled check at n={n} ({n_samples} samples)")
    for trial in range(n_samples):
        if trial % 500 == 0:
            print(f"  ... {trial}/{n_samples}", flush=True)
        bits = rng.randint(0, n_total)
        T = tournament_from_bits(n, bits)
        try:
            betti = compute_betti(T, n)
        except:
            continue

        profile = tuple(betti.get(p, 0) for p in range(n))
        betti_profiles_7[profile] += 1

        for k in range(1, n-1):
            if betti.get(k, 0) > 0 and betti.get(k+1, 0) > 0:
                violations_7.append((bits, k, betti.get(k,0), betti.get(k+1,0), profile))

    print(f"\n  Betti profiles at n=7:")
    for profile, cnt in betti_profiles_7.most_common(15):
        print(f"    {list(profile)}: {cnt}")

    print(f"\n  Consecutive seesaw violations at n=7: {len(violations_7)}")
    if violations_7:
        for bits, k, bk, bk1, prof in violations_7[:5]:
            print(f"    bits={bits}: beta_{k}={bk}, beta_{k+1}={bk1}, profile={list(prof)}")

    # Part 3: Focus on beta_3 * beta_4
    print(f"\n{'='*70}")
    print("Part 3: beta_3 * beta_4 specifically")
    print("=" * 70)

    b3_b4_dist = Counter()
    for profile, cnt in betti_profiles_7.items():
        b3_b4_dist[(profile[3], profile[4])] += cnt

    print(f"\n  (beta_3, beta_4) distribution at n=7:")
    for (b3, b4), cnt in sorted(b3_b4_dist.items()):
        product = b3 * b4
        marker = " ***VIOLATION***" if product > 0 else ""
        print(f"    ({b3}, {b4}): {cnt}{marker}")

    # Part 4: Also check adjacent even-odd
    print(f"\n{'='*70}")
    print("Part 4: ALL adjacency products at n=7")
    print("=" * 70)

    for k in range(1, n-1):
        products = Counter()
        for profile, cnt in betti_profiles_7.items():
            products[profile[k] * profile[k+1]] += cnt
        nonzero = sum(c for p, c in products.items() if p > 0)
        print(f"  beta_{k} * beta_{k+1}: products = {dict(sorted(products.items()))}, violations = {nonzero}")

    # Part 5: Check at n=8 with small sample
    n = 8
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs
    rng8 = np.random.RandomState(99)
    n_samples_8 = 200

    violations_8 = []
    betti_profiles_8 = Counter()

    print(f"\nPart 5: Sampled check at n={n} ({n_samples_8} samples)")
    for trial in range(n_samples_8):
        if trial % 50 == 0:
            print(f"  ... {trial}/{n_samples_8}", flush=True)
        bits = rng8.randint(0, n_total)
        T = tournament_from_bits(n, bits)
        try:
            betti = compute_betti(T, n)
        except:
            continue

        profile = tuple(betti.get(p, 0) for p in range(n))
        betti_profiles_8[profile] += 1

        for k in range(1, n-1):
            if betti.get(k, 0) > 0 and betti.get(k+1, 0) > 0:
                violations_8.append((bits, k, betti.get(k,0), betti.get(k+1,0), profile))

    print(f"\n  Betti profiles at n=8:")
    for profile, cnt in betti_profiles_8.most_common(15):
        print(f"    {list(profile)}: {cnt}")

    print(f"\n  Consecutive seesaw violations at n=8: {len(violations_8)}")
    if violations_8:
        for bits, k, bk, bk1, prof in violations_8[:5]:
            print(f"    bits={bits}: beta_{k}={bk}, beta_{k+1}={bk1}, profile={list(prof)}")

    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print("=" * 70)
    print(f"  n=6: {len(violations_6)} consecutive seesaw violations (exhaustive)")
    print(f"  n=7: {len(violations_7)} consecutive seesaw violations ({n_samples} samples)")
    print(f"  n=8: {len(violations_8)} consecutive seesaw violations ({n_samples_8} samples)")

if __name__ == '__main__':
    main()
