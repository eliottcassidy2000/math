#!/usr/bin/env python3
"""
relative_odd_mechanism.py — opus-2026-03-09-S52

Investigate WHY H_{2k+1}(T, T\v) <= 1.

The relative complex R_p = Omega_p(T) / Omega_p(T\v) consists of
"v-dependent" elements — those involving vertex v.

Key observations to check:
1. R_p dimensions: how big is the relative complex?
2. Are the relative boundary maps "almost surjective"?
3. Is there a structural reason the cokernel at odd levels is ≤ 1?

Specific focus: H_1(T, T\v) and H_3(T, T\v) at n=6.
Both are ≤ 1. What's the mechanism?

For H_1: this is related to how v's cocycles interact with T\v.
  H_1(T, T\v) = 1 iff there's a new 1-cycle using v that's not
  killed by the relative boundary from R_2.

For H_3: this is related to how v's 3-cycles interact with T\v.
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

def compute_relative_details(T, n, v, max_p=5):
    """Compute detailed relative complex structure."""
    tol = 1e-8

    all_T = {}
    all_Tv = {}
    for p in range(0, max_p + 2):
        paths_T = []
        if p + 1 <= n:
            for verts in combinations(range(n), p+1):
                for perm in permutations(verts):
                    if is_allowed_path(T, perm):
                        paths_T.append(perm)
        all_T[p] = paths_T
        all_Tv[p] = [path for path in paths_T if v not in path]

    def build_omega(all_allowed, max_p_val):
        omega_proj = {}
        for p in range(0, max_p_val + 1):
            a_p = all_allowed[p]
            if not a_p:
                omega_proj[p] = np.zeros((0, 0))
                continue
            a_pm1_set = set(all_allowed[p-1]) if p > 0 else set()
            na_faces_map = {}
            for sigma in a_p:
                for i in range(1, len(sigma)-1):
                    face = sigma[:i] + sigma[i+1:]
                    if p > 0 and face not in a_pm1_set:
                        if face not in na_faces_map:
                            na_faces_map[face] = len(na_faces_map)
            if not na_faces_map:
                omega_proj[p] = np.eye(len(a_p))
            else:
                na_mat = np.zeros((len(na_faces_map), len(a_p)))
                for j, sigma in enumerate(a_p):
                    for i in range(1, len(sigma)-1):
                        face = sigma[:i] + sigma[i+1:]
                        if face in na_faces_map:
                            na_mat[na_faces_map[face], j] += (-1)**i
                U, S, Vt = np.linalg.svd(na_mat, full_matrices=True)
                rank = int(np.sum(S > tol))
                null_dim = len(a_p) - rank
                if null_dim == 0:
                    omega_proj[p] = np.zeros((len(a_p), 0))
                else:
                    omega_proj[p] = Vt[rank:].T
        return omega_proj

    omega_T = build_omega(all_T, max_p)
    omega_Tv = build_omega(all_Tv, max_p)

    boundary_T = {}
    for p in range(1, max_p + 1):
        a_p = all_T[p]
        a_pm1 = all_T[p-1]
        if not a_p or not a_pm1:
            boundary_T[p] = np.zeros((len(a_pm1) if a_pm1 else 0, len(a_p) if a_p else 0))
            continue
        idx_pm1 = {path: i for i, path in enumerate(a_pm1)}
        mat = np.zeros((len(a_pm1), len(a_p)))
        for j, sigma in enumerate(a_p):
            for i in range(len(sigma)):
                face = sigma[:i] + sigma[i+1:]
                if face in idx_pm1:
                    mat[idx_pm1[face], j] += (-1)**i
        boundary_T[p] = mat

    omega_Tv_in_T = {}
    for p in range(0, max_p + 1):
        a_p_T = all_T[p]
        a_p_Tv = all_Tv[p]
        om_Tv = omega_Tv[p]
        if om_Tv.ndim < 2 or om_Tv.shape[1] == 0 or not a_p_T:
            omega_Tv_in_T[p] = np.zeros((len(a_p_T) if a_p_T else 0, 0))
            continue
        idx_T = {path: i for i, path in enumerate(a_p_T)}
        mapping = np.zeros((len(a_p_T), len(a_p_Tv)))
        for j, path in enumerate(a_p_Tv):
            if path in idx_T:
                mapping[idx_T[path], j] = 1.0
        omega_Tv_in_T[p] = mapping @ om_Tv

    # Compute complement and relative complex
    complement = {}
    rel_dims = {}
    for p in range(0, max_p + 1):
        om_T_p = omega_T[p]
        om_Tv_p = omega_Tv_in_T[p]
        if om_T_p.ndim < 2 or om_T_p.shape[1] == 0:
            complement[p] = np.zeros((0, 0))
            rel_dims[p] = 0
            continue
        if om_Tv_p.shape[1] == 0:
            complement[p] = om_T_p
            rel_dims[p] = om_T_p.shape[1]
            continue
        proj = om_T_p.T @ om_Tv_p
        U, S, Vt = np.linalg.svd(proj, full_matrices=True)
        sub_rank = int(np.sum(S > tol))
        if sub_rank < om_T_p.shape[1]:
            comp_coords = U[:, sub_rank:]
            complement[p] = om_T_p @ comp_coords
            rel_dims[p] = comp_coords.shape[1]
        else:
            complement[p] = np.zeros((om_T_p.shape[0], 0))
            rel_dims[p] = 0

    rel_rank = {}
    for p in range(1, max_p + 1):
        comp_p = complement[p]
        om_T_pm1 = omega_T[p-1]
        comp_pm1 = complement[p-1]
        if comp_p.ndim < 2 or comp_p.shape[1] == 0 or comp_pm1.ndim < 2 or comp_pm1.shape[1] == 0:
            rel_rank[p] = 0
            continue
        if p not in boundary_T or boundary_T[p].shape[0] == 0:
            rel_rank[p] = 0
            continue
        dp_comp = boundary_T[p] @ comp_p
        dp_in_omega = om_T_pm1.T @ dp_comp
        comp_pm1_in_omega = om_T_pm1.T @ comp_pm1
        if comp_pm1_in_omega.shape[1] == 0:
            rel_rank[p] = 0
            continue
        Q, R = np.linalg.qr(comp_pm1_in_omega, mode='reduced')
        proj_dp = Q.T @ dp_in_omega
        S = np.linalg.svd(proj_dp, compute_uv=False)
        rel_rank[p] = int(np.sum(S > tol))

    rel_betti = {}
    for p in range(0, max_p + 1):
        ker = rel_dims[p] - rel_rank.get(p, 0)
        im_next = rel_rank.get(p+1, 0)
        rel_betti[p] = ker - im_next

    abs_dims = {p: omega_T[p].shape[1] if omega_T[p].ndim == 2 else 0 for p in range(max_p+1)}
    sub_dims = {p: omega_Tv[p].shape[1] if omega_Tv[p].ndim == 2 else 0 for p in range(max_p+1)}

    return {
        'rel_betti': rel_betti,
        'rel_dims': rel_dims,
        'rel_rank': rel_rank,
        'abs_dims': abs_dims,
        'sub_dims': sub_dims
    }

def main():
    print("=" * 70)
    print("RELATIVE COMPLEX STRUCTURE: WHY H_{odd}(T,T\\v) <= 1?")
    print("=" * 70)

    n = 6
    print(f"\nn = {n}: exhaustive analysis of relative complex dimensions")

    num_arcs = n*(n-1)//2
    total = 1 << num_arcs

    # Track detailed data
    data = {p: {'dims': [], 'ranks': [], 'ker': [], 'im': [], 'betti': []}
            for p in range(6)}

    for bits in range(total):
        if bits % 5000 == 0 and bits > 0:
            print(f"  ... {bits}/{total}")

        T = tournament_from_bits(n, bits)
        for v in range(n):
            res = compute_relative_details(T, n, v)
            for p in range(6):
                data[p]['dims'].append(res['rel_dims'].get(p, 0))
                data[p]['ranks'].append(res['rel_rank'].get(p, 0))
                ker = res['rel_dims'].get(p, 0) - res['rel_rank'].get(p, 0)
                im_next = res['rel_rank'].get(p+1, 0)
                data[p]['ker'].append(ker)
                data[p]['im'].append(im_next)
                data[p]['betti'].append(res['rel_betti'].get(p, 0))

    print(f"\n{'='*70}")
    print("RELATIVE COMPLEX STATISTICS")
    print("="*70)

    for p in range(6):
        dims = data[p]['dims']
        kers = data[p]['ker']
        ims = data[p]['im']
        bettis = data[p]['betti']

        if not dims or max(dims) == 0:
            continue

        print(f"\n  Level p={p}:")
        print(f"    dim(R_p): min={min(dims)}, max={max(dims)}, mean={np.mean(dims):.1f}")
        print(f"    ker(d_p^rel): min={min(kers)}, max={max(kers)}, mean={np.mean(kers):.1f}")
        print(f"    im(d_{p+1}^rel): min={min(ims)}, max={max(ims)}, mean={np.mean(ims):.1f}")
        print(f"    H_p(T,T\\v): {dict(Counter(bettis))}")

        # KEY: compute the "gap" = ker - im for each parity
        gap_counts = Counter(k - i for k, i in zip(kers, ims))
        print(f"    gap = ker - im: {dict(sorted(gap_counts.items()))}")

        # Deficiency: how far is rank(d_{p+1}^rel) from ker(d_p^rel)?
        if max(kers) > 0:
            deficiency = [k - i for k, i in zip(kers, ims) if k > 0]
            if deficiency:
                def_counts = Counter(deficiency)
                print(f"    deficiency (gap when ker>0): {dict(sorted(def_counts.items()))}")

    # CRITICAL: Focus on odd levels
    print(f"\n{'='*70}")
    print("ODD-LEVEL RELATIVE HOMOLOGY ANALYSIS")
    print("="*70)

    for p in [1, 3, 5]:
        bettis = data[p]['betti']
        dims = data[p]['dims']
        kers = data[p]['ker']
        ims = data[p]['im']

        print(f"\n  H_{p}(T, T\\v):")
        print(f"    Distribution: {dict(sorted(Counter(bettis).items()))}")

        # When H_p > 0: what are the ker and im values?
        pos_cases = [(d, k, i, b) for d, k, i, b in zip(dims, kers, ims, bettis) if b > 0]
        if pos_cases:
            print(f"    When H_{p} > 0 ({len(pos_cases)} cases):")
            ker_vals = Counter(c[1] for c in pos_cases)
            im_vals = Counter(c[2] for c in pos_cases)
            dim_vals = Counter(c[0] for c in pos_cases)
            print(f"      dim(R_{p}): {dict(sorted(dim_vals.items()))}")
            print(f"      ker(d_{p}^rel): {dict(sorted(ker_vals.items()))}")
            print(f"      im(d_{p+1}^rel): {dict(sorted(im_vals.items()))}")

    # EVEN-LEVEL for comparison
    print(f"\n{'='*70}")
    print("EVEN-LEVEL RELATIVE HOMOLOGY (for comparison)")
    print("="*70)

    for p in [2, 4]:
        bettis = data[p]['betti']
        if max(bettis) > 0:
            print(f"\n  H_{p}(T, T\\v):")
            print(f"    Distribution: {dict(sorted(Counter(bettis).items()))}")

            pos_cases = [(data[p]['dims'][i], data[p]['ker'][i], data[p]['im'][i], bettis[i])
                        for i in range(len(bettis)) if bettis[i] > 0]
            if pos_cases:
                print(f"    When H_{p} > 0 ({len(pos_cases)} cases):")
                ker_vals = Counter(c[1] for c in pos_cases)
                print(f"      ker(d_{p}^rel): {dict(sorted(ker_vals.items()))}")

if __name__ == '__main__':
    main()
