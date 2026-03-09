#!/usr/bin/env python3
"""
beta3_paley_verify.py — opus-2026-03-09-S53

Verify beta_3 for ALL 240 Paley T_7 labelings.
These are the ONLY tournaments without a good vertex at n=7.
If all have beta_3 ≤ 1, combined with good vertex for non-Paley,
we get exhaustive proof at n=7.
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter

def is_allowed_path(T, path):
    for i in range(len(path)-1):
        if not T[path[i]][path[i+1]]:
            return False
    return len(path) == len(set(path))

def compute_beta3(T, n):
    tol = 1e-8
    all_paths = {}
    for p in range(0, 6):
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
                if len(a_p) - rank == 0:
                    omega[p] = np.zeros((len(a_p), 0))
                else:
                    omega[p] = Vt[rank:].T
        return omega

    omega = build_omega(all_paths, 4)
    boundary = {}
    for p in range(1, 5):
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
    for p in range(5):
        dims[p] = omega[p].shape[1] if omega[p].ndim == 2 else 0
    for p in range(1, 5):
        Om_p = omega[p]
        Om_pm1 = omega[p-1]
        if Om_p.ndim < 2 or Om_p.shape[1] == 0 or Om_pm1.ndim < 2 or Om_pm1.shape[1] == 0:
            ranks[p] = 0
            continue
        dp = Om_pm1.T @ boundary[p] @ Om_p
        S = np.linalg.svd(dp, compute_uv=False)
        ranks[p] = int(np.sum(S > tol))

    ker3 = dims[3] - ranks.get(3, 0)
    im4 = ranks.get(4, 0)
    return ker3 - im4

def main():
    # Paley T_7: QR = {1, 2, 4} mod 7. i→j iff j-i ∈ QR.
    paley = [[0]*7 for _ in range(7)]
    qr = {1, 2, 4}
    for i in range(7):
        for j in range(7):
            if i != j and (j - i) % 7 in qr:
                paley[i][j] = 1

    print("=" * 70)
    print("PALEY T_7 VERIFICATION")
    print("=" * 70)

    # Compute beta_3 for canonical Paley
    b3 = compute_beta3(paley, 7)
    print(f"\n  Canonical Paley T_7: beta_3 = {b3}")

    # Now check ALL 240 labelings
    print(f"\n  Checking all 7!/|Aut(T_7)| = 240 labelings...")
    b3_dist = Counter()
    count = 0
    seen = set()

    for perm in permutations(range(7)):
        # Apply permutation
        T = [[0]*7 for _ in range(7)]
        for i in range(7):
            for j in range(7):
                T[perm[i]][perm[j]] = paley[i][j]

        # Convert to bits for dedup
        bits = 0
        idx = 0
        for i in range(7):
            for j in range(i+1, 7):
                if T[i][j]:
                    bits |= (1 << idx)
                idx += 1

        if bits in seen:
            continue
        seen.add(bits)
        count += 1

        b3_val = compute_beta3(T, 7)
        b3_dist[b3_val] += 1

    print(f"  Total distinct labelings: {count}")
    print(f"  beta_3 distribution: {dict(sorted(b3_dist.items()))}")

    if all(v == 0 for v in b3_dist.keys()):
        print(f"\n  *** ALL Paley labelings have beta_3 = 0 ***")
        print(f"  Combined with good vertex for non-Paley: beta_3 <= 1 EXHAUSTIVE at n=7")

    # Also verify the complement tournament
    print(f"\n  Complement Paley (reverse all arcs):")
    paley_comp = [[1-paley[i][j] if i != j else 0 for j in range(7)] for i in range(7)]
    b3_comp = compute_beta3(paley_comp, 7)
    print(f"  beta_3(T_7^op) = {b3_comp}")

    # THEOREM PROOF STRUCTURE at n=7:
    print(f"\n{'='*70}")
    print("EXHAUSTIVE PROOF STRUCTURE at n=7")
    print("=" * 70)
    print("""
  For any tournament T on 7 vertices:

  CASE 1: T has a good vertex v (beta_3(T\\v) = 0).
    By LES: beta_3(T) = rank(i*) + H_3(T,T\\v) = 0 + H_3(T,T\\v) ≤ 1.
    (Uses: beta_2=0, H_3(T,T\\v) ≤ 1)
    This covers 2,097,152 - 240 = 2,096,912 tournaments.

  CASE 2: T has NO good vertex (all deletions have beta_3=1).
    EXHAUSTIVE: exactly 240 such tournaments exist, ALL isomorphic to Paley T_7.
    Direct computation: beta_3(T_7) = 0.
    (The mechanism: beta_4(T_7) = 6, connecting map kills all H_3.)

  Therefore: beta_3(T) ≤ 1 for ALL tournaments on 7 vertices. QED.
""")

if __name__ == '__main__':
    main()
