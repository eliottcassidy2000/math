"""
codim1_mechanism.py — WHY does codim(π_old(im d_4), π_old(ker d_3)) = 1?

From HYP-413: this is equivalent to the Ghost Cycle Theorem.

Key insight: codim_old = beta_3 iff π_old preserves the quotient structure.
Specifically: the natural map ker(d_3)/im(d_4) → π_old(ker)/π_old(im) is injective.

This map sends [z] ↦ [π_old(z)]. It's well-defined since π_old(im) ⊂ π_old(ker).
It's injective iff: π_old(z) ∈ π_old(im) implies z ∈ im(d_4).

Equivalently: if z ∈ ker(d_3) and π_old(z) = π_old(b) for some b ∈ im(d_4),
then z - b ∈ ker(π_old) ∩ ker(d_3), meaning z - b is a tv-only cycle.
So the question is: is z - b ∈ im(d_4)?

This is exactly the Ghost Cycle Theorem: tv-only cycles are boundaries.

NEW APPROACH: Study the map d_4 restricted to through-v 4-chains,
projected to old 3-coordinates. Call this map φ: Ω_4^tv → A_3^old.

φ = π_old ∘ d_4|_{Ω_4^tv}

The single-face theorem says each through-v 4-path contributes
exactly one old 3-face (the v-deletion face). So φ is essentially
the "v-deletion face" map.

If we can show: im(φ) + im(d_4^{old→old}) = π_old(ker(d_3)),
that would give codim_old = 1.

Actually: π_old(ker(d_3)) has dimension dim(ker) - dim(K_tv) = dim(ker) - dim(K_tv).
And π_old(im(d_4)) has dimension dim(im) - dim(B_tv).
We need: dim(ker) - dim(K_tv) - [dim(im) - dim(B_tv)] = 1.
i.e., (dim(ker) - dim(im)) - (dim(K_tv) - dim(B_tv)) = 1.
i.e., beta_3 = 1 + dim(K_tv) - dim(B_tv).

So codim_old = 1 iff K_tv = B_tv iff Ghost Cycle.

CIRCULAR. We need a DIRECT argument for why codim_old = 1.

New idea: relate π_old(ker(d_3)) to ker(d_3^{T\v}).
The non-v paths in T are EXACTLY the paths in T\v (no orphans, from S57).
So old 3-cycles in T are related to 3-cycles in T\v.

Specifically: if z ∈ ker(d_3^T) and π_old(z) ≠ 0, then π_old(z) might be
in ker(d_3^{T\v}) after relabeling. But π_old(z) is NOT necessarily a cycle
in T\v because d_3^T(z) = 0 involves BOTH old and tv 2-faces.

The old-part of d_3(z) = 0 means: the old 2-faces sum to minus the tv 2-faces.
But both are zero (d_3(z)=0 means ALL 2-face components sum to zero).

Wait: d_3(z) = 0 means the full boundary is zero. The boundary has both old
and tv 2-path components. For a MIXED cycle z (both old and tv 3-paths),
the old-to-old and tv-to-old boundary components must cancel, and similarly
for the old-to-tv and tv-to-tv components.

Actually: d_3 maps 3-paths to 2-paths. A 3-path has 4 faces (alternating sign).
An old 3-path (a,b,c,d) with v ∉ {a,b,c,d} has faces:
  (b,c,d), (a,c,d), (a,b,d), (a,b,c) — all old.

A tv 3-path with v at position k has faces:
  positions 0,1,2,3, deleting one vertex at a time.
  The face with v deleted is old; the other 3 faces are tv.

So d_3(old) → old only. d_3(tv) → tv + old.

For z ∈ ker(d_3): d_3(z) = 0.
Split z = z_old + z_tv. Then:
  d_3(z_old) + d_3(z_tv) = 0
  [d_3(z_old)]_old + [d_3(z_tv)]_old = 0  ... (in old 2-coords)
  [d_3(z_tv)]_tv = 0                       ... (in tv 2-coords)

Since d_3(old) is purely old: [d_3(z_old)]_old = d_3(z_old).
So: d_3(z_old) + [d_3(z_tv)]_old = 0, meaning d_3(z_old) = -[d_3(z_tv)]_old.

And: [d_3(z_tv)]_tv = 0 (the tv-to-tv part of d_3 applied to z_tv vanishes).

Now, π_old(z) = z_old. When is z_old a cycle in T\v?
z_old is a cycle in T\v iff d_3^{T\v}(z_old) = 0.
But d_3^{T\v} is NOT the same as d_3^T restricted to old paths!
In GLMY homology, the boundary map depends on the Omega constraints,
which differ between T and T\v.

Actually wait: at the A-path level, d_3 is the same map (alternating face).
The Omega constraints determine WHICH chains are valid, but the boundary
map on valid chains is the same combinatorial map.

So d_3^T(z_old) = d_3^{T\v}(relabel(z_old)) (modulo relabeling).
And d_3(z_old) = -[d_3(z_tv)]_old ≠ 0 in general.

So z_old is NOT a cycle in T\v (unless z is tv-only, in which case z_old=0).

This is the key obstruction: the old part of a T-cycle is NOT a T\v-cycle.

Let's compute: what is d_3^{T\v}(π_old(z)) for z ∈ ker(d_3^T)?
It equals d_3(z_old) = -[d_3(z_tv)]_old.

This residual boundary tells us HOW FAR z_old is from being a T\v-cycle.

Author: opus-2026-03-10-S59
"""
import sys
import time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, deletion_adj,
    enumerate_all_allowed,
    _build_constraint_matrix, _gauss_rank_np, _gauss_nullbasis_modp,
    full_chain_complex_modp, boundary_faces,
    RANK_PRIME
)

PRIME = RANK_PRIME


def analyze_old_projection_boundary(A, n, v):
    """For each z in ker(d_3^T), compute d_3^{T\v}(π_old(z))."""
    max_p = min(n - 1, 6)
    ap_T = enumerate_all_allowed(A, n, max_p)

    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]
    remaining_inv = {remaining[i]: i for i in range(n1)}
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    paths_3_T = ap_T.get(3, [])
    paths_4_T = ap_T.get(4, [])
    paths_2_T = ap_T.get(2, [])
    paths_3_Tv = ap_Tv.get(3, [])
    paths_2_Tv = ap_Tv.get(2, [])

    if not paths_3_T or not paths_4_T:
        return None

    tv3 = [i for i, p in enumerate(paths_3_T) if v in p]
    old3 = [i for i, p in enumerate(paths_3_T) if v not in p]
    tv2 = [i for i, p in enumerate(paths_2_T) if v in p]
    old2 = [i for i, p in enumerate(paths_2_T) if v not in p]

    if not tv3 or not old3:
        return None

    def get_omega(ap, deg):
        paths = ap.get(deg, [])
        if not paths:
            return np.zeros((0, 0), dtype=np.int64)
        P, nr, nc = _build_constraint_matrix(ap, deg, PRIME)
        if P is not None:
            _, nb = _gauss_nullbasis_modp(P, nr, nc, PRIME)
            return np.array(nb, dtype=np.int64) if nb else np.zeros((0, nc), dtype=np.int64)
        return np.eye(len(paths), dtype=np.int64)

    ob3_T = get_omega(ap_T, 3)
    ob4_T = get_omega(ap_T, 4)
    ob3_Tv = get_omega(ap_Tv, 3)

    # d_3^T in A-path coords
    idx2_T = {p: i for i, p in enumerate(paths_2_T)}
    bd3_T = np.zeros((len(paths_2_T), len(paths_3_T)), dtype=np.int64)
    for j, path in enumerate(paths_3_T):
        for sign, face in boundary_faces(path):
            if face in idx2_T:
                bd3_T[idx2_T[face], j] = (bd3_T[idx2_T[face], j] + sign) % PRIME

    # d_3^{T\v} in A-path coords
    idx2_Tv = {p: i for i, p in enumerate(paths_2_Tv)}
    bd3_Tv = np.zeros((len(paths_2_Tv), len(paths_3_Tv)), dtype=np.int64)
    for j, path in enumerate(paths_3_Tv):
        for sign, face in boundary_faces(path):
            if face in idx2_Tv:
                bd3_Tv[idx2_Tv[face], j] = (bd3_Tv[idx2_Tv[face], j] + sign) % PRIME

    # ker(d_3^T) in Omega coords
    d3o = bd3_T @ ob3_T.T % PRIME
    if d3o.size > 0:
        _, kv = _gauss_nullbasis_modp(d3o.astype(np.int32), d3o.shape[0], d3o.shape[1], PRIME)
        ker_d3_T = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob3_T.shape[0]), dtype=np.int64)
    else:
        ker_d3_T = np.eye(ob3_T.shape[0], dtype=np.int64)

    K = ker_d3_T @ ob3_T % PRIME  # ker(d_3^T) in A_3(T) coords
    dim_K = K.shape[0]

    if dim_K == 0:
        return None

    # im(d_4^T) in A_3(T) coords
    idx3_T = {p: i for i, p in enumerate(paths_3_T)}
    bd4_T = np.zeros((len(paths_3_T), len(paths_4_T)), dtype=np.int64)
    for j, path in enumerate(paths_4_T):
        for sign, face in boundary_faces(path):
            if face in idx3_T:
                bd4_T[idx3_T[face], j] = (bd4_T[idx3_T[face], j] + sign) % PRIME

    im_d4 = (bd4_T @ ob4_T.T % PRIME).T if ob4_T.shape[0] > 0 else np.zeros((0, len(paths_3_T)), dtype=np.int64)
    rk_d4 = int(_gauss_rank_np(im_d4.copy(), PRIME)) if im_d4.shape[0] > 0 else 0
    beta3 = dim_K - rk_d4

    if beta3 != 1:
        return None

    # Map old 3-paths of T to 3-paths of T\v
    old3_T_to_Tv = {}
    idx3_Tv = {p: i for i, p in enumerate(paths_3_Tv)}
    for i in old3:
        path_T = paths_3_T[i]
        path_Tv = tuple(remaining_inv[x] for x in path_T)
        if path_Tv in idx3_Tv:
            old3_T_to_Tv[i] = idx3_Tv[path_Tv]

    # π_old(K): restrict ker(d_3^T) to old coordinates, map to T\v coords
    # K_old_Tv[k, j] = K[k, old3[i]] where old3_T_to_Tv[old3[i]] = j
    K_old_Tv = np.zeros((dim_K, len(paths_3_Tv)), dtype=np.int64)
    for i_idx in old3:
        if i_idx in old3_T_to_Tv:
            j_Tv = old3_T_to_Tv[i_idx]
            K_old_Tv[:, j_Tv] = (K_old_Tv[:, j_Tv] + K[:, i_idx]) % PRIME

    # ker(d_3^{T\v}) in Omega coords
    d3o_Tv = bd3_Tv @ ob3_Tv.T % PRIME
    if d3o_Tv.size > 0:
        _, kv_Tv = _gauss_nullbasis_modp(d3o_Tv.astype(np.int32), d3o_Tv.shape[0], d3o_Tv.shape[1], PRIME)
        ker_d3_Tv = np.array(kv_Tv, dtype=np.int64) if kv_Tv else np.zeros((0, ob3_Tv.shape[0]), dtype=np.int64)
    else:
        ker_d3_Tv = np.eye(ob3_Tv.shape[0], dtype=np.int64)

    ker_d3_Tv_A = ker_d3_Tv @ ob3_Tv % PRIME

    # Compute d_3^{T\v}(π_old(z)) for each z in ker(d_3^T)
    # This is bd3_Tv @ K_old_Tv^T
    residual = bd3_Tv @ K_old_Tv.T % PRIME  # shape: (len(paths_2_Tv), dim_K)

    # How many ker(d_3^T) basis vectors have zero residual?
    n_zero_residual = 0
    for k in range(dim_K):
        col = residual[:, k] % PRIME
        if np.all(col == 0):
            n_zero_residual += 1

    rk_residual = int(_gauss_rank_np(residual.copy(), PRIME))

    # π_old(K) lies in what relation to ker(d_3^{T\v})?
    # π_old(K) ∩ ker(d_3^{T\v}): vectors in K_old_Tv that are T\v-cycles
    if ker_d3_Tv_A.shape[0] > 0 and K_old_Tv.shape[0] > 0:
        combined_ker = np.vstack([ker_d3_Tv_A, K_old_Tv]) % PRIME
        rk_comb_ker = int(_gauss_rank_np(combined_ker.copy(), PRIME))
        rk_ker_Tv = int(_gauss_rank_np(ker_d3_Tv_A.copy(), PRIME))
        rk_K_old = int(_gauss_rank_np(K_old_Tv.copy(), PRIME))
        dim_intersection = rk_ker_Tv + rk_K_old - rk_comb_ker
    else:
        dim_intersection = 0
        rk_ker_Tv = 0
        rk_K_old = 0

    # im(d_4^{T\v})
    paths_4_Tv = ap_Tv.get(4, [])
    ob4_Tv = get_omega(ap_Tv, 4)
    bd4_Tv = np.zeros((len(paths_3_Tv), len(paths_4_Tv)), dtype=np.int64) if paths_4_Tv else np.zeros((len(paths_3_Tv), 0), dtype=np.int64)
    for j, path in enumerate(paths_4_Tv):
        for sign, face in boundary_faces(path):
            if face in idx3_Tv:
                bd4_Tv[idx3_Tv[face], j] = (bd4_Tv[idx3_Tv[face], j] + sign) % PRIME

    im_d4_Tv = (bd4_Tv @ ob4_Tv.T % PRIME).T if ob4_Tv.shape[0] > 0 else np.zeros((0, len(paths_3_Tv)), dtype=np.int64)
    rk_d4_Tv = int(_gauss_rank_np(im_d4_Tv.copy(), PRIME)) if im_d4_Tv.shape[0] > 0 else 0
    beta3_Tv = rk_ker_Tv - rk_d4_Tv

    return {
        'dim_K': dim_K,
        'rk_d4': rk_d4,
        'beta3': beta3,
        'rk_residual': rk_residual,  # rank of d_3^{T\v}(π_old(K))
        'n_zero_residual': n_zero_residual,  # how many K basis vecs map to T\v-cycles
        'rk_K_old': rk_K_old,  # rank of π_old(K) in T\v coords
        'rk_ker_Tv': rk_ker_Tv,  # dim ker(d_3^{T\v})
        'rk_d4_Tv': rk_d4_Tv,  # dim im(d_4^{T\v})
        'beta3_Tv': beta3_Tv,
        'dim_intersection': dim_intersection,  # dim(π_old(K) ∩ ker(d_3^{T\v}))
        'n_old3': len(old3),
        'n_tv3': len(tv3),
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"OLD PROJECTION BOUNDARY RESIDUAL AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        t0 = time.time()
        target = 300 if n <= 7 else 200

        for trial in range(50000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) != 1:
                continue

            for v_cand in range(n):
                # Only BAD vertices
                remaining = [i for i in range(n) if i != v_cand]
                n1 = n - 1
                A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]
                cc_Tv = full_chain_complex_modp(A_sub, n1, min(n1 - 1, 6))
                if cc_Tv['bettis'].get(3, 0) != 1:
                    continue  # only BAD vertices

                r = analyze_old_projection_boundary(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)

        elapsed = time.time() - t0
        print(f"  {len(results)} BAD (T,v) pairs, {elapsed:.1f}s")

        for key in ['dim_K', 'rk_d4', 'rk_residual', 'n_zero_residual',
                     'rk_K_old', 'rk_ker_Tv', 'rk_d4_Tv', 'beta3_Tv',
                     'dim_intersection']:
            vals = [r[key] for r in results]
            print(f"  {key}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")

        # KEY: rk_K_old vs rk_ker_Tv
        # π_old(ker(d_3^T)) rank vs ker(d_3^{T\v}) rank
        print(f"\n  rk(π_old(K)) - rk(ker(d_3^Tv)):")
        diff = [r['rk_K_old'] - r['rk_ker_Tv'] for r in results]
        print(f"    {dict(sorted(Counter(diff).items()))}")

        # dim(π_old(K) ∩ ker(d_3^{T\v}))
        print(f"\n  dim(π_old(K) ∩ ker(d_3^Tv)):")
        print(f"    {dict(sorted(Counter(r['dim_intersection'] for r in results).items()))}")

        # KEY: does im(d_4^{T\v}) ⊂ π_old(K)?
        # The embedded im(d_4^{T\v}) should map into π_old(im(d_4^T)) ⊂ π_old(K)

        # Residual structure: rk_residual tells us how many directions in π_old(K)
        # are NOT T\v-cycles
        print(f"\n  rk(residual) = rk(d_3^Tv(π_old(K))):")
        print(f"    {dict(sorted(Counter(r['rk_residual'] for r in results).items()))}")

        # Decomposition: π_old(K) = [part in ker(d_3^Tv)] + [part with nonzero d_3^Tv-boundary]
        # dim(π_old(K)) = dim_intersection + (rk_K_old - dim_intersection)
        # The nonzero-boundary part has dimension = rk_residual (rank-nullity)
        # So: rk_K_old = dim_intersection + rk_residual? Not exactly...
        # Actually: rk_K_old - dim_intersection = number of independent directions
        # in π_old(K) that are NOT T\v-cycles
        non_cycle_dim = [r['rk_K_old'] - r['dim_intersection'] for r in results]
        print(f"\n  Non-T\\v-cycle directions in π_old(K):")
        print(f"    {dict(sorted(Counter(non_cycle_dim).items()))}")
        print(f"    rk_residual: {dict(sorted(Counter(r['rk_residual'] for r in results).items()))}")
        print(f"    (These should be equal)")


if __name__ == '__main__':
    main()
    print("\nDONE.")
