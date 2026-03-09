#!/usr/bin/env python3
"""
beta3_regular7_deep.py — opus-2026-03-09-S53

The 4 regular tournaments on 7 vertices where ALL 7 deletions have beta_3=1
yet beta_3(T)=0 are deeply interesting.

This means:
- beta_3(T\v) = 1 for all v
- H_3(T,T\v) = beta_3(T) - rank(i_*) = 0 - rank(i_*) ≤ 0
  But H_3 ≥ 0, so H_3(T,T\v) = 0 and rank(i_*) = 0 for all v.

Wait - that's impossible if beta_3(T)=0. The map i_*: H_3(T\v) → H_3(T)
goes FROM a 1-dim space TO a 0-dim space, so rank(i_*)=0 always.
Then H_3(T,T\v) = 0 - 0 = 0.

So these tournaments have:
  beta_3(T) = 0, beta_3(T\v) = 1 for all v, H_3(T,T\v) = 0 for all v.

The LES gives: ... → H_3(T\v)=Z → H_3(T)=0 → H_3(T,T\v)=0 → H_2(T\v)=0 → ...
So the map H_3(T\v) → H_3(T) has im=0 (trivially, target is 0).
And H_3(T,T\v)=0 means ker(H_3(T,T\v) → H_2(T\v))=0.

This is consistent. The 3-cycle in T\v gets "killed" when we embed in T.

Let me find these tournaments and study their structure.
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

def compute_beta3(T, n):
    tol = 1e-8
    all_paths = {}
    max_p = min(5, n)
    for p in range(0, max_p + 1):
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

def compute_full_betti(T, n):
    """Compute all Betti numbers."""
    tol = 1e-8
    all_paths = {}
    max_p = min(n, 7)
    for p in range(0, max_p):
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
            a_p = pd.get(p, [])
            if not a_p:
                omega[p] = np.zeros((0, 0))
                continue
            a_pm1_set = set(pd.get(p-1, [])) if p > 0 else set()
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

    mp = max_p - 1
    omega = build_omega(all_paths, mp)
    boundary = {}
    for p in range(1, mp + 1):
        a_p = all_paths.get(p, [])
        a_pm1 = all_paths.get(p-1, [])
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
    for p in range(mp + 1):
        dims[p] = omega[p].shape[1] if omega[p].ndim == 2 else 0
    for p in range(1, mp + 1):
        Om_p = omega[p]
        Om_pm1 = omega[p-1]
        if Om_p.ndim < 2 or Om_p.shape[1] == 0 or Om_pm1.ndim < 2 or Om_pm1.shape[1] == 0:
            ranks[p] = 0
            continue
        dp = Om_pm1.T @ boundary[p] @ Om_p
        S = np.linalg.svd(dp, compute_uv=False)
        ranks[p] = int(np.sum(S > tol))

    betti = []
    for p in range(mp + 1):
        ker = dims[p] - ranks.get(p, 0)
        im_next = ranks.get(p+1, 0)
        betti.append(ker - im_next)
    return betti

def main():
    n = 7
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs

    print("=" * 70)
    print("DEEP ANALYSIS: Regular tournaments on 7 vertices")
    print("Finding those with ALL deletions having beta_3=1")
    print("=" * 70)

    rng = np.random.RandomState(42)

    # Find the special tournaments
    special = []
    total_regular = 0

    # Need more samples to find them
    for trial in range(200000):
        bits = rng.randint(0, n_total)
        T = tournament_from_bits(n, bits)
        scores = tuple(sorted([sum(T[i][j] for j in range(n) if j != i) for i in range(n)]))
        if scores != (3,3,3,3,3,3,3):
            continue

        total_regular += 1

        # Check all deletions
        all_b3_1 = True
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
            try:
                b3_sub = compute_beta3(T_sub, n-1)
                if b3_sub != 1:
                    all_b3_1 = False
                    break
            except:
                all_b3_1 = False
                break

        if all_b3_1:
            special.append((bits, T))
            # Compute full Betti for T
            betti = compute_full_betti(T, n)
            print(f"\n  Found special tournament #{len(special)} (bits={bits})")
            print(f"    Betti numbers: {betti}")

            # Check adjacency structure
            print(f"    Adjacency matrix:")
            for row in T:
                print(f"      {row}")

            # Check: is this the Paley tournament?
            # Paley T_7: QR = {1,2,4}, i→j if j-i in QR (mod 7)
            paley = [[0]*7 for _ in range(7)]
            qr = {1, 2, 4}
            for i in range(7):
                for j in range(7):
                    if i != j and (j - i) % 7 in qr:
                        paley[i][j] = 1

            # Check if T is isomorphic to Paley (by checking all 7! relabelings)
            is_paley = False
            for perm in permutations(range(7)):
                match = True
                for i in range(7):
                    for j in range(7):
                        if T[i][j] != paley[perm[i]][perm[j]]:
                            match = False
                            break
                    if not match:
                        break
                if match:
                    is_paley = True
                    break
            print(f"    Isomorphic to Paley T_7? {is_paley}")

            # Check deletion beta profiles
            for v in range(n):
                remaining = [i for i in range(n) if i != v]
                T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
                del_betti = compute_full_betti(T_sub, n-1)
                print(f"    T\\{v}: betti = {del_betti}")

            if len(special) >= 6:
                break

    print(f"\n  Total regular tournaments sampled: {total_regular}")
    print(f"  Special tournaments (all deletions beta_3=1): {len(special)}")

    # Key analysis: for these tournaments
    # beta_3(T) = 0, beta_3(T\v) = 1 for all v
    # By LES: H_3(T\v) → H_3(T) → H_3(T,T\v) → H_2(T\v)
    #          Z        → 0       → ?             → 0
    # So H_3(T,T\v) = 0 (exact seq: 0 → H_3(T,T\v) → 0)
    # And the map Z → 0 has kernel Z = im(connecting map from H_4(T,T\v))
    # So H_4(T,T\v) → H_3(T\v) is surjective!
    # This means rank(delta_4: H_4(T,T\v) → H_3(T\v)) = 1

    print(f"\n{'='*70}")
    print("LES ANALYSIS for special tournaments")
    print("='*70")
    print("""
  For these T: beta_3(T)=0, beta_3(T\\v)=1 for all v.

  LES: ... → H_4(T,T\\v) →δ H_3(T\\v) →i* H_3(T) → H_3(T,T\\v) → H_2(T\\v) → ...
                            Z=1           0            ?             0

  Since H_3(T)=0, the map i*: H_3(T\\v)→H_3(T) is zero.
  Exactness at H_3(T\\v): im(δ) = ker(i*) = H_3(T\\v) = Z.
  So δ: H_4(T,T\\v) → H_3(T\\v) is SURJECTIVE.

  Since H_3(T)=0 and H_2(T\\v)=0:
  H_3(T,T\\v) = 0 (sandwiched between zero maps).

  This means the relative 3-homology vanishes even though
  the absolute 3-homology of the subcomplex is nontrivial!
  The connecting map from relative H_4 does all the work.
""")

if __name__ == '__main__':
    main()
