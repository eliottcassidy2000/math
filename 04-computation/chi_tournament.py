#!/usr/bin/env python3
"""
chi_tournament.py — opus-2026-03-09-S52

CRITICAL DISCOVERY: chi(T) = sum (-1)^p beta_p ∈ {0, 1} at n ≤ 6!

If chi ≥ 0 always, and beta_2 = 0, then at n where beta_4 = 0:
  1 - beta_1 - beta_3 ≥ 0
  beta_1 + beta_3 ≤ 1
  Boolean odd Betti follows (with seesaw: beta_1 * beta_3 = 0)!

Check: does chi ≥ 0 always? Does chi ∈ {0, 1} extend?
What happens at n=7 where beta_4 can be nonzero?
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

def compute_omega_dims(T, n):
    """Compute dim(Omega_p) for all p — enough for Euler characteristic."""
    tol = 1e-8
    all_paths = {}
    for p in range(n):
        paths = []
        if p + 1 <= n:
            for verts in combinations(range(n), p+1):
                for perm in permutations(verts):
                    if is_allowed_path(T, perm):
                        paths.append(perm)
        all_paths[p] = paths

    omega_dims = {}
    for p in range(n):
        a_p = all_paths[p]
        if not a_p:
            omega_dims[p] = 0
            continue
        a_pm1_set = set(all_paths[p-1]) if p > 0 else set()
        na = {}
        for sigma in a_p:
            for i in range(1, len(sigma)-1):
                face = sigma[:i] + sigma[i+1:]
                if p > 0 and face not in a_pm1_set:
                    if face not in na:
                        na[face] = len(na)
        if not na:
            omega_dims[p] = len(a_p)
        else:
            mat = np.zeros((len(na), len(a_p)))
            for j, sigma in enumerate(a_p):
                for i in range(1, len(sigma)-1):
                    face = sigma[:i] + sigma[i+1:]
                    if face in na:
                        mat[na[face], j] += (-1)**i
            try:
                S = np.linalg.svd(mat, compute_uv=False)
                rank = int(np.sum(S > tol))
            except:
                U, S, Vt = np.linalg.svd(mat, full_matrices=False)
                rank = int(np.sum(S > tol))
            omega_dims[p] = len(a_p) - rank

    return omega_dims

def main():
    print("=" * 70)
    print("EULER CHARACTERISTIC OF TOURNAMENTS: chi(T) = sum (-1)^p dim(Omega_p)")
    print("=" * 70)

    # Exhaustive for n ≤ 6, sampled for n=7
    for n in [3, 4, 5, 6]:
        num_arcs = n*(n-1)//2
        n_total = 1 << num_arcs

        chi_dist = Counter()
        for bits in range(n_total):
            T = tournament_from_bits(n, bits)
            dims = compute_omega_dims(T, n)
            chi = sum((-1)**p * dims.get(p, 0) for p in range(n))
            chi_dist[chi] += 1

        print(f"\n  n={n}: chi distribution = {dict(sorted(chi_dist.items()))}")
        print(f"    Total: {sum(chi_dist.values())}, chi ≥ 0: {all(k >= 0 for k in chi_dist)}")

    # n=7 sampled
    print(f"\n  n=7 (sampled):")
    n = 7
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs
    rng = np.random.RandomState(42)
    sample_size = 500

    chi_dist = Counter()
    dim_profiles = Counter()
    for trial in range(sample_size):
        if trial % 100 == 0:
            print(f"    ... {trial}/{sample_size}")
        bits = rng.randint(0, n_total)
        T = tournament_from_bits(n, bits)
        dims = compute_omega_dims(T, n)
        chi = sum((-1)**p * dims.get(p, 0) for p in range(n))
        chi_dist[chi] += 1
        dim_profiles[tuple(dims.get(p, 0) for p in range(n))] += 1

    print(f"  n=7: chi distribution = {dict(sorted(chi_dist.items()))}")
    print(f"    chi ≥ 0: {all(k >= 0 for k in chi_dist)}")
    print(f"    min chi: {min(chi_dist.keys())}, max chi: {max(chi_dist.keys())}")

    # Check Paley T_7
    print(f"\n  Paley T_7:")
    qr = {1, 2, 4}
    T_paley = [[int((j-i) % 7 in qr) if i != j else 0 for j in range(7)] for i in range(7)]
    dims = compute_omega_dims(T_paley, 7)
    chi = sum((-1)**p * dims.get(p, 0) for p in range(7))
    print(f"    Omega dims: {[dims.get(p,0) for p in range(7)]}")
    print(f"    chi = {chi}")

    # Also check: Omega Euler chi formula
    # chi = sum (-1)^p dim(Omega_p)
    # For tournaments, dim(Omega_0) = n always.
    # Question: is there a formula for chi in terms of tournament invariants?
    print(f"\n{'='*70}")
    print("CHI vs TOURNAMENT INVARIANTS")
    print("=" * 70)

    # Check if chi correlates with number of 3-cycles
    n = 6
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs

    chi_by_t3 = {}
    for bits in range(n_total):
        T = tournament_from_bits(n, bits)
        # Count 3-cycles
        t3 = 0
        for a in range(n):
            for b in range(a+1, n):
                for c in range(b+1, n):
                    if T[a][b] + T[b][c] + T[c][a] in [0, 3]:
                        t3 += 1
        # Also count Hamiltonian paths (but that's expensive, skip)
        dims = compute_omega_dims(T, n)
        chi = sum((-1)**p * dims.get(p, 0) for p in range(n))
        if t3 not in chi_by_t3:
            chi_by_t3[t3] = Counter()
        chi_by_t3[t3][chi] += 1

    print(f"\n  n=6: chi by number of 3-cycles (t3):")
    for t3 in sorted(chi_by_t3.keys()):
        print(f"    t3={t3}: chi distribution = {dict(sorted(chi_by_t3[t3].items()))}")

    # KEY QUESTION: is chi = 1 iff beta_1 = 0?
    # chi = 1 - beta_1 - beta_3 (at n=6)
    # chi = 1 iff beta_1 = beta_3 = 0
    # chi = 0 iff exactly one of beta_1, beta_3 is 1
    print(f"\n  chi = 1 iff acyclic (beta_1 = beta_3 = 0)")
    print(f"  chi = 0 iff exactly one of beta_1, beta_3 = 1")

    # Final check: dim(Omega_0) for tournaments
    print(f"\n{'='*70}")
    print("OMEGA DIMENSIONS — IS THERE A PATTERN?")
    print("=" * 70)

    for n_check in [3, 4, 5, 6]:
        num_arcs = n_check*(n_check-1)//2
        n_total = 1 << num_arcs
        dim_sets = {p: set() for p in range(n_check)}
        for bits in range(n_total):
            T = tournament_from_bits(n_check, bits)
            dims = compute_omega_dims(T, n_check)
            for p in range(n_check):
                dim_sets[p].add(dims.get(p, 0))

        print(f"\n  n={n_check}:")
        for p in range(n_check):
            vals = sorted(dim_sets[p])
            if len(vals) <= 10:
                print(f"    dim(Omega_{p}): {vals}")
            else:
                print(f"    dim(Omega_{p}): min={min(vals)}, max={max(vals)}, |range|={len(vals)}")

if __name__ == '__main__':
    main()
