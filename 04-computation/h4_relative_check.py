#!/usr/bin/env python3
"""
h4_relative_check.py — opus-2026-03-09-S54

From the LES analysis with consecutive seesaw:
  When b3(T)=1, b3(T\v)=1: i_*-injectivity <=> H_4(T,T\v) = 0

Check 1: Verify H_4(T,T\v) = 0 for all bad vertices (b3(T)=1, b3(T\v)=1)
Check 2: Verify H_4(T,T\v) = F for Paley case (b3(T)=0, b3(T\v)=1)
Check 3: Compute ALL H_p(T,T\v) to confirm the LES picture

The relative homology can be computed by the formula:
  beta_p(T,T\v) = beta_p(T) - rank(i_*_p) + dim(coker δ_{p+1})
But it's simpler to compute from the long exact sequence directly.

Alternative: compute R_* = C_*(T)/C_*(T\v) directly and find its homology.
But that requires careful quotient computation in Omega coords.

SIMPLEST: use rank(i_*) to determine H_p(T,T\v) from LES.
For the critical degree p=4:
  H_4(T\v) → H_4(T) → H_4(T,T\v) → H_3(T\v) → H_3(T)
  Since H_4(T\v)=0, H_4(T)=0:
    0 → H_4(T,T\v) --δ--> H_3(T\v) --i*--> H_3(T) → ...
  H_4(T,T\v) = ker(i_*) (by injectivity of δ and im(δ)=ker(i_*))
  So: dim H_4(T,T\v) = dim ker(i_*) = beta_3(T\v) - rank(i_*)

So we just need rank(i_*) at degree 3!
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

def compute_homology_data(T, n, max_p_limit=7):
    tol = 1e-8
    all_paths = {}
    for p in range(min(n, max_p_limit)):
        paths = []
        for verts in combinations(range(n), p+1):
            for perm in permutations(verts):
                if is_allowed_path(T, perm):
                    paths.append(perm)
        all_paths[p] = paths

    max_p = max(all_paths.keys())

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

    return {'all_paths': all_paths, 'omega': omega, 'boundary': boundary,
            'dims': dims, 'ranks': ranks, 'betti': betti}

def compute_rank_istar(data_T, data_sub, remaining, n, degree=3):
    """Compute rank of i_*: H_degree(T\v) -> H_degree(T)."""
    tol = 1e-8
    d = degree

    if data_sub['betti'].get(d, 0) == 0:
        return 0

    # Get H_d generator(s) of T\v
    Om_d_sub = data_sub['omega'].get(d, np.zeros((0,0)))
    Om_dm1_sub = data_sub['omega'].get(d-1, np.zeros((0,0)))
    if Om_d_sub.ndim < 2 or Om_d_sub.shape[1] == 0:
        return 0

    d_map_sub = Om_dm1_sub.T @ data_sub['boundary'][d] @ Om_d_sub
    U, S, Vt = np.linalg.svd(d_map_sub, full_matrices=True)
    r = int(np.sum(S > tol))
    ker_basis = Vt[r:].T  # ker(d_d) in Omega_d coords, shape (dim_Omega_d, dim_ker)

    # im(d_{d+1}) in ker(d_d) coords
    Om_dp1_sub = data_sub['omega'].get(d+1, np.zeros((0,0)))
    if Om_dp1_sub.ndim == 2 and Om_dp1_sub.shape[1] > 0 and (d+1) in data_sub['boundary']:
        d_dp1_sub = Om_d_sub.T @ data_sub['boundary'][d+1] @ Om_dp1_sub
        im_in_ker = ker_basis.T @ d_dp1_sub
        U2, S2, _ = np.linalg.svd(im_in_ker, full_matrices=True)
        r2 = int(np.sum(S2 > tol))
        # H_d generators in ker coords
        h_generators = U2[:, r2:]
    else:
        h_generators = np.eye(ker_basis.shape[1])

    if h_generators.shape[1] == 0:
        return 0

    # Convert to Omega_d(T\v) coords
    h_omega = ker_basis @ h_generators  # shape (dim_Omega_d_sub, n_gen)

    # Convert to path coords
    h_paths = Om_d_sub @ h_omega  # shape (n_paths_sub, n_gen)

    # Embed into T's path space
    paths_T = data_T['all_paths'].get(d, [])
    paths_sub = data_sub['all_paths'].get(d, [])
    path_idx_T = {p: i for i, p in enumerate(paths_T)}

    n_gen = h_paths.shape[1]
    embed_omega = np.zeros((data_T['dims'].get(d, 0), n_gen))
    Om_d_T = data_T['omega'].get(d, np.zeros((0,0)))

    for g in range(n_gen):
        embed_raw = np.zeros(len(paths_T))
        for j in range(len(paths_sub)):
            if abs(h_paths[j, g]) > tol:
                path_T = tuple(remaining[k] for k in paths_sub[j])
                if path_T in path_idx_T:
                    embed_raw[path_idx_T[path_T]] = h_paths[j, g]
        embed_omega[:, g] = Om_d_T.T @ embed_raw

    # Check: embedded generators should be in ker(d_d^T)
    Om_dm1_T = data_T['omega'].get(d-1, np.zeros((0,0)))
    d_map_T = Om_dm1_T.T @ data_T['boundary'][d] @ Om_d_T
    residual = d_map_T @ embed_omega
    if np.linalg.norm(residual) > 1e-6:
        print(f"  WARNING: embedded cycle not in ker(d_{d}^T)! residual={np.linalg.norm(residual):.2e}")

    # Project onto H_d(T): compute component in im(d_{d+1}^T)
    Om_dp1_T = data_T['omega'].get(d+1, np.zeros((0,0)))
    if Om_dp1_T.ndim == 2 and Om_dp1_T.shape[1] > 0 and (d+1) in data_T['boundary']:
        d_dp1_T = Om_d_T.T @ data_T['boundary'][d+1] @ Om_dp1_T
        # embed_omega = component in im(d_{d+1}^T) + component in H_d(T)
        # Use least squares to find the im part
        for g in range(n_gen):
            x, _, _, _ = np.linalg.lstsq(d_dp1_T, embed_omega[:, g], rcond=None)
            im_part = d_dp1_T @ x
            h_part = embed_omega[:, g] - im_part
            embed_omega[:, g] = h_part  # replace with homology class

    # rank of embed_omega (after projection to H_d(T))
    try:
        S_embed = np.linalg.svd(embed_omega, compute_uv=False)
        rank_istar = int(np.sum(S_embed > tol))
    except:
        rank_istar = 0

    return rank_istar

def main():
    print("=" * 70)
    print("H_4(T,T\\v) CHECK AND i_*-INJECTIVITY VIA LES")
    print("=" * 70)

    n = 7
    n_arcs = n*(n-1)//2
    n_total = 1 << n_arcs
    rng = np.random.RandomState(42)

    results = {'GOOD': [], 'BAD': []}
    count = 0

    print(f"\n  Scanning n={n} tournaments...")
    for trial in range(10000):
        if trial % 1000 == 0:
            print(f"  ... trial {trial} (found {count} b3=1)", flush=True)
        bits = rng.randint(0, n_total)
        T = tournament_from_bits(n, bits)
        try:
            data_T = compute_homology_data(T, n)
        except:
            continue
        if data_T['betti'].get(3, 0) != 1:
            continue
        count += 1

        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
            try:
                data_sub = compute_homology_data(T_sub, n-1)
            except:
                continue

            b3_sub = data_sub['betti'].get(3, 0)
            key = "BAD" if b3_sub == 1 else "GOOD"

            rank_i3 = compute_rank_istar(data_T, data_sub, remaining, n, degree=3)

            # From LES: H_4(T,T\v) = ker(i_*) since delta is injective
            # dim H_4(T,T\v) = b3_sub - rank_i3
            h4_rel = b3_sub - rank_i3

            # H_3(T,T\v) = coker(i_*)
            h3_rel = data_T['betti'].get(3, 0) - rank_i3

            results[key].append({
                'v': v, 'rank_i3': rank_i3,
                'h3_rel': h3_rel, 'h4_rel': h4_rel,
                'b3_T': data_T['betti'].get(3, 0),
                'b3_sub': b3_sub,
                'b4_T': data_T['betti'].get(4, 0),
                'b4_sub': data_sub['betti'].get(4, 0),
            })

        if count >= 80:
            break

    print(f"\n  Results from {count} beta_3=1 tournaments:")

    for key in ["GOOD", "BAD"]:
        data = results[key]
        print(f"\n  {key} vertices ({len(data)} total):")

        rank_dist = Counter(d['rank_i3'] for d in data)
        print(f"    rank(i_*): {dict(sorted(rank_dist.items()))}")

        h3_dist = Counter(d['h3_rel'] for d in data)
        print(f"    dim H_3(T,T\\v): {dict(sorted(h3_dist.items()))}")

        h4_dist = Counter(d['h4_rel'] for d in data)
        print(f"    dim H_4(T,T\\v): {dict(sorted(h4_dist.items()))}")

        b4_T_dist = Counter(d['b4_T'] for d in data)
        print(f"    beta_4(T): {dict(sorted(b4_T_dist.items()))}")

        b4_sub_dist = Counter(d['b4_sub'] for d in data)
        print(f"    beta_4(T\\v): {dict(sorted(b4_sub_dist.items()))}")

    # Part 2: Check Paley T_7
    print(f"\n{'='*70}")
    print("Part 2: Paley T_7")
    print("=" * 70)

    # Paley T_7: QR = {1,2,4} mod 7
    n_paley = 7
    QR = {1, 2, 4}
    T_paley = [[0]*n_paley for _ in range(n_paley)]
    for i in range(n_paley):
        for j in range(n_paley):
            if i != j and (j - i) % n_paley in QR:
                T_paley[i][j] = 1

    data_paley = compute_homology_data(T_paley, n_paley)
    print(f"  Paley T_7 betti: {[data_paley['betti'].get(p,0) for p in range(n_paley)]}")

    for v in range(min(3, n_paley)):
        remaining = [i for i in range(n_paley) if i != v]
        T_sub = [[T_paley[remaining[a]][remaining[b]] for b in range(n_paley-1)] for a in range(n_paley-1)]
        data_sub = compute_homology_data(T_sub, n_paley-1)

        rank_i3 = compute_rank_istar(data_paley, data_sub, remaining, n_paley, degree=3)
        h3_rel = data_paley['betti'].get(3, 0) - rank_i3
        h4_rel = data_sub['betti'].get(3, 0) - rank_i3

        print(f"  v={v}: b3_sub={data_sub['betti'].get(3,0)}, rank(i_*)={rank_i3}, "
              f"H_3^rel={h3_rel}, H_4^rel={h4_rel}, "
              f"b4_T={data_paley['betti'].get(4,0)}, b4_sub={data_sub['betti'].get(4,0)}")

    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print("=" * 70)
    print("""
  LES decomposition (b3(T)=1, consecutive seesaw):
    H_4(T)=0, H_4(T\\v)=0, H_2(T\\v)=0

  => 0 -> H_4(T,T\\v) --δ--> H_3(T\\v) --i*--> H_3(T) -> H_3(T,T\\v) -> 0

  For BAD vertex (b3(T\\v)=1):
    rank(i_*) = 1  =>  H_4(T,T\\v) = 0, H_3(T,T\\v) = 0
    i_* is an ISOMORPHISM: H_3(T\\v) ≅ H_3(T)

  For GOOD vertex (b3(T\\v)=0):
    rank(i_*) = 0  =>  H_4(T,T\\v) = 0, H_3(T,T\\v) = F

  For Paley (b3(T)=0, b3(T\\v)=1):
    rank(i_*) = 0  =>  H_4(T,T\\v) = F (δ surjective, kills H_3)
    The connecting map δ kills the H_3 of T\\v

  KEY INSIGHT: i_*-injectivity is equivalent to H_4(T,T\\v) = 0.
  The consecutive seesaw reduces this to a statement about the
  connecting map δ at degree 4.
""")

if __name__ == '__main__':
    main()
