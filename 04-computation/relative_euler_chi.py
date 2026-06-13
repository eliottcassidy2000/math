#!/usr/bin/env python3
"""
relative_euler_chi.py — opus-2026-03-09-S52

Investigate the Euler characteristic of the relative complex.

Key equation: chi(T,T\v) = sum (-1)^p H_p(T,T\v) = chi(T) - chi(T\v)

If chi is fixed/predictable, it constrains the relative Betti numbers.
For tournaments: chi = 1 always? Or does it vary?

Also investigate: sum of relative Betti = sum H_p(T,T\v).
If most H_p = 0 and at most one odd H_p = 1, what does this imply?
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

def compute_full_betti(T, n, max_p=None):
    """Compute all Betti numbers."""
    if max_p is None:
        max_p = n - 1
    tol = 1e-8
    all_paths = {}
    for p in range(0, max_p + 2):
        paths = []
        if p + 1 <= n:
            for verts in combinations(range(n), p+1):
                for perm in permutations(verts):
                    if is_allowed_path(T, perm):
                        paths.append(perm)
        all_paths[p] = paths

    def build_omega(pd, mp):
        omega = {}
        for p in range(0, mp + 1):
            a_p = pd[p]
            if not a_p:
                omega[p] = np.zeros((0, 0))
                continue
            a_pm1_set = set(pd[p-1]) if p > 0 else set()
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
                U, S, Vt = np.linalg.svd(mat, full_matrices=True)
                rank = int(np.sum(S > tol))
                null_dim = len(a_p) - rank
                if null_dim == 0:
                    omega[p] = np.zeros((len(a_p), 0))
                else:
                    omega[p] = Vt[rank:].T
        return omega

    omega = build_omega(all_paths, max_p)

    boundary = {}
    for p in range(1, max_p + 1):
        a_p = all_paths[p]
        a_pm1 = all_paths[p-1]
        if not a_p or not a_pm1:
            boundary[p] = np.zeros((0, 0))
            continue
        idx = {path: i for i, path in enumerate(a_pm1)}
        mat = np.zeros((len(a_pm1), len(a_p)))
        for j, sigma in enumerate(a_p):
            for i in range(len(sigma)):
                face = sigma[:i] + sigma[i+1:]
                if face in idx:
                    mat[idx[face], j] += (-1)**i
        boundary[p] = mat

    ranks = {}
    dims = {}
    for p in range(max_p + 1):
        dims[p] = omega[p].shape[1] if omega[p].ndim == 2 else 0
    for p in range(1, max_p + 1):
        Om_p = omega[p]
        Om_pm1 = omega[p-1]
        if Om_p.ndim < 2 or Om_p.shape[1] == 0 or Om_pm1.ndim < 2 or Om_pm1.shape[1] == 0:
            ranks[p] = 0
            continue
        dp = Om_pm1.T @ boundary[p] @ Om_p
        S = np.linalg.svd(dp, compute_uv=False)
        ranks[p] = int(np.sum(S > tol))

    betti = {}
    for p in range(max_p + 1):
        ker = dims[p] - ranks.get(p, 0)
        im_next = ranks.get(p+1, 0)
        betti[p] = ker - im_next

    return betti, dims

def main():
    print("=" * 70)
    print("EULER CHARACTERISTIC OF TOURNAMENTS AND RELATIVE COMPLEXES")
    print("=" * 70)

    for n in [4, 5, 6]:
        print(f"\n{'='*70}")
        print(f"n = {n}")
        print(f"{'='*70}")

        num_arcs = n*(n-1)//2
        n_total = 1 << num_arcs

        chi_dist = Counter()
        omega_chi_dist = Counter()
        betti_profiles = Counter()

        limit = n_total if n <= 6 else 200

        for bits in range(limit):
            if bits % 5000 == 0 and bits > 0:
                print(f"  ... {bits}/{limit}")
            T = tournament_from_bits(n, bits)
            betti, dims = compute_full_betti(T, n)

            # Euler characteristic of the chain complex
            chi = sum((-1)**p * betti.get(p, 0) for p in range(n))
            chi_dist[chi] += 1

            # Omega Euler characteristic
            omega_chi = sum((-1)**p * dims.get(p, 0) for p in range(n))
            omega_chi_dist[omega_chi] += 1

            # Betti profile
            profile = tuple(betti.get(p, 0) for p in range(n))
            betti_profiles[profile] += 1

        print(f"\n  Euler characteristic chi = sum (-1)^p beta_p:")
        print(f"    Distribution: {dict(sorted(chi_dist.items()))}")

        print(f"\n  Omega Euler characteristic = sum (-1)^p dim(Omega_p):")
        print(f"    Distribution: {dict(sorted(omega_chi_dist.items()))}")

        print(f"\n  Betti profiles (beta_0, beta_1, ..., beta_{n-1}):")
        for profile, count in sorted(betti_profiles.items(), key=lambda x: -x[1])[:15]:
            chi_val = sum((-1)**p * profile[p] for p in range(len(profile)))
            print(f"    {profile}: {count} (chi={chi_val})")

    # ===========================================
    # Part 2: Relative Euler characteristic
    # ===========================================
    print(f"\n{'='*70}")
    print("RELATIVE EULER CHARACTERISTIC: chi(T) - chi(T\\v)")
    print("=" * 70)

    n = 6
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs

    rel_chi_dist = Counter()
    chi_T_dist = Counter()
    chi_Tv_dist = Counter()

    for bits in range(n_total):
        if bits % 5000 == 0 and bits > 0:
            print(f"  ... {bits}/{n_total}")
        T = tournament_from_bits(n, bits)
        betti_T, _ = compute_full_betti(T, n)
        chi_T = sum((-1)**p * betti_T.get(p, 0) for p in range(n))

        for v in range(1):  # just v=0 for speed
            remaining = [i for i in range(n) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
            betti_sub, _ = compute_full_betti(T_sub, n-1)
            chi_Tv = sum((-1)**p * betti_sub.get(p, 0) for p in range(n-1))

            rel_chi = chi_T - chi_Tv
            rel_chi_dist[rel_chi] += 1
            chi_T_dist[chi_T] += 1
            chi_Tv_dist[chi_Tv] += 1

    print(f"\n  chi(T) distribution: {dict(sorted(chi_T_dist.items()))}")
    print(f"  chi(T\\v) distribution: {dict(sorted(chi_Tv_dist.items()))}")
    print(f"  chi(T) - chi(T\\v) distribution: {dict(sorted(rel_chi_dist.items()))}")

    # ===========================================
    # Part 3: If beta_even = 0, what does chi tell us?
    # ===========================================
    print(f"\n{'='*70}")
    print("SIMPLIFICATION: IF ALL beta_even = 0")
    print("=" * 70)
    print("""
    beta_0 = 1 always.
    beta_2 = 0 (THM-108).
    beta_4 can be > 0 (Paley T_7 has beta_4 = 6).

    But at n <= 6: beta_4 = 0 (no 4-cycles in Omega).
    So chi(T) = 1 - beta_1 + 0 - beta_3 + 0 - beta_5 = 1 - beta_1 - beta_3 - beta_5.

    By seesaw: beta_1 * beta_3 = 0, so at most one of beta_1, beta_3 is nonzero.

    chi(T) = 1 - max(beta_1, beta_3) - beta_5 (at n <= 6).

    If chi is constant (= 1):
      max(beta_1, beta_3) + beta_5 = 0, so all = 0.
    If chi varies but is ≥ 0:
      max(beta_1, beta_3) + beta_5 ≤ 1.
      Since beta_5 ≥ 0 and max(beta_1, beta_3) ≤ 1:
      This would imply max(beta_1, beta_3) ∈ {0, 1} and beta_5 = 0 when nonzero.
    """)

    # Verify: is chi always 1 for tournaments?
    print(f"  Checking if chi = 1 always:")
    for n_check in [3, 4, 5, 6]:
        num_arcs = n_check*(n_check-1)//2
        n_total = 1 << num_arcs
        chis = set()
        for bits in range(n_total):
            T = tournament_from_bits(n_check, bits)
            betti, _ = compute_full_betti(T, n_check)
            chi = sum((-1)**p * betti.get(p, 0) for p in range(n_check))
            chis.add(chi)
        print(f"    n={n_check}: chi values = {sorted(chis)}")

if __name__ == '__main__':
    main()
