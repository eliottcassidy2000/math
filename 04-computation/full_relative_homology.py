"""
full_relative_homology.py — Complete H_*(T, T\v) for b3=1 tournaments

Compute ALL relative Betti numbers for (T,v) pairs.
If the relative complex is acyclic (all β_p^rel = 0),
then i: T\v → T is a homology equivalence.

This would be a much stronger result than just H_4^rel = 0.

Method: Compute using the SES of chain complexes and the LES.
Since we know β_k(T) and β_k(T\v) and the ranks of i_*,
we can determine β_k^rel from the LES.

Actually: compute i_* at ALL degrees and use exactness.

Author: opus-2026-03-10-S59
"""
import sys
import time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament,
    enumerate_all_allowed,
    _build_constraint_matrix, _gauss_rank_np, _gauss_nullbasis_modp,
    full_chain_complex_modp, boundary_faces,
    RANK_PRIME
)

PRIME = RANK_PRIME


def compute_istar_all_degrees(A, n, v, max_p=None):
    """Compute rank(i_*: H_p(T\v) → H_p(T)) for all p."""
    if max_p is None:
        max_p = min(n - 1, 6)

    cc_T = full_chain_complex_modp(A, n, max_p)
    bettis_T = cc_T['bettis']

    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]
    cc_Tv = full_chain_complex_modp(A_sub, n1, min(max_p, n1-1))
    bettis_Tv = cc_Tv['bettis']

    if bettis_T.get(3, 0) != 1 or bettis_Tv.get(3, 0) != 1:
        return None

    # For computing i_* at a specific degree p:
    # i_*: H_p(T\v) → H_p(T) = ker(d_p^{T\v})/im(d_{p+1}^{T\v}) → ker(d_p^T)/im(d_{p+1}^T)
    # The map embeds a T\v cycle (in A_p(T\v) coords) into A_p(T) coords,
    # then checks its class in H_p(T).
    #
    # rank(i_*) = β_p(T\v) iff i_* is injective.
    # rank(i_*) = dim(im(i_*))
    #
    # For computing: embed cycle generators of H_p(T\v) into A_p(T),
    # then compute how many are independent mod im(d_{p+1}^T).

    ap_T = enumerate_all_allowed(A, n, max_p)
    remaining_inv = {remaining[i]: i for i in range(n1)}
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    paths_T = {p: ap_T.get(p, []) for p in range(2, max_p + 1)}
    paths_Tv = {p: ap_Tv.get(p, []) for p in range(2, min(max_p, n1-1) + 1)}

    def get_omega(ap, deg):
        ps = ap.get(deg, [])
        if not ps:
            return np.zeros((0, 0), dtype=np.int64)
        P, nr, nc = _build_constraint_matrix(ap, deg, PRIME)
        if P is not None:
            _, nb = _gauss_nullbasis_modp(P, nr, nc, PRIME)
            return np.array(nb, dtype=np.int64) if nb else np.zeros((0, nc), dtype=np.int64)
        return np.eye(len(ps), dtype=np.int64)

    istar_ranks = {}

    for p in range(2, max_p + 1):
        bp_Tv = bettis_Tv.get(p, 0)
        bp_T = bettis_T.get(p, 0)

        if bp_Tv == 0:
            istar_ranks[p] = 0
            continue
        if bp_T == 0:
            istar_ranks[p] = 0
            continue

        # Need to compute i_* at degree p
        # Get cycle generators for T\v
        ob_Tv_p = get_omega(ap_Tv, p)
        ob_T_p = get_omega(ap_T, p)
        ob_Tv_pp1 = get_omega(ap_Tv, p+1)
        ob_T_pp1 = get_omega(ap_T, p+1)

        idx_T_pm1 = {path: i for i, path in enumerate(paths_T.get(p-1, []))}
        idx_Tv_pm1 = {path: i for i, path in enumerate(paths_Tv.get(p-1, []))}

        # d_p^{T\v} in Omega coords
        bd_Tv_p = np.zeros((len(paths_Tv.get(p-1, [])), len(paths_Tv.get(p, []))), dtype=np.int64)
        for j, path in enumerate(paths_Tv.get(p, [])):
            for sign, face in boundary_faces(path):
                if face in idx_Tv_pm1:
                    bd_Tv_p[idx_Tv_pm1[face], j] = (bd_Tv_p[idx_Tv_pm1[face], j] + sign) % PRIME

        d_Tv_p_omega = bd_Tv_p @ ob_Tv_p.T % PRIME if ob_Tv_p.shape[0] > 0 else np.zeros((len(paths_Tv.get(p-1, [])), 0), dtype=np.int64)

        # ker(d_p^{T\v})
        if d_Tv_p_omega.size > 0 and d_Tv_p_omega.shape[1] > 0:
            _, kv = _gauss_nullbasis_modp(
                d_Tv_p_omega.astype(np.int32), d_Tv_p_omega.shape[0], d_Tv_p_omega.shape[1], PRIME)
            ker_Tv = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob_Tv_p.shape[0]), dtype=np.int64)
        else:
            ker_Tv = np.eye(ob_Tv_p.shape[0], dtype=np.int64) if ob_Tv_p.shape[0] > 0 else np.zeros((0, 0), dtype=np.int64)

        # im(d_{p+1}^{T\v}): cycles modulo boundaries
        idx_Tv_p = {path: i for i, path in enumerate(paths_Tv.get(p, []))}
        bd_Tv_pp1 = np.zeros((len(paths_Tv.get(p, [])), len(paths_Tv.get(p+1, []))), dtype=np.int64)
        for j, path in enumerate(paths_Tv.get(p+1, [])):
            for sign, face in boundary_faces(path):
                if face in idx_Tv_p:
                    bd_Tv_pp1[idx_Tv_p[face], j] = (bd_Tv_pp1[idx_Tv_p[face], j] + sign) % PRIME

        im_d_Tv_pp1 = (bd_Tv_pp1 @ ob_Tv_pp1.T % PRIME).T if ob_Tv_pp1.shape[0] > 0 else np.zeros((0, len(paths_Tv.get(p, []))), dtype=np.int64)

        # H_p(T\v) generators: ker_Tv mod im_d_Tv_pp1
        # ker_Tv in A_p(T\v) coords
        ker_Tv_A = ker_Tv @ ob_Tv_p % PRIME if ker_Tv.shape[0] > 0 else np.zeros((0, len(paths_Tv.get(p, []))), dtype=np.int64)

        # Embed into A_p(T) coords
        embed_p = np.zeros((len(paths_Tv.get(p, [])), len(paths_T.get(p, []))), dtype=np.int64)
        idx_T_p = {path: i for i, path in enumerate(paths_T.get(p, []))}
        for jv, pv in enumerate(paths_Tv.get(p, [])):
            mapped = tuple(remaining[x] for x in pv)
            if mapped in idx_T_p:
                embed_p[jv, idx_T_p[mapped]] = 1

        ker_Tv_embedded = ker_Tv_A @ embed_p % PRIME if ker_Tv_A.shape[0] > 0 else np.zeros((0, len(paths_T.get(p, []))), dtype=np.int64)

        # im(d_{p+1}^T)
        idx_T_p_map = {path: i for i, path in enumerate(paths_T.get(p, []))}
        bd_T_pp1 = np.zeros((len(paths_T.get(p, [])), len(paths_T.get(p+1, []))), dtype=np.int64)
        for j, path in enumerate(paths_T.get(p+1, [])):
            for sign, face in boundary_faces(path):
                if face in idx_T_p_map:
                    bd_T_pp1[idx_T_p_map[face], j] = (bd_T_pp1[idx_T_p_map[face], j] + sign) % PRIME

        im_d_T_pp1 = (bd_T_pp1 @ ob_T_pp1.T % PRIME).T if ob_T_pp1.shape[0] > 0 else np.zeros((0, len(paths_T.get(p, []))), dtype=np.int64)

        # rank(i_*) = dim(ker_Tv_embedded + im_d_T_pp1) - dim(im_d_T_pp1)
        rk_im_T = int(_gauss_rank_np(im_d_T_pp1.copy(), PRIME)) if im_d_T_pp1.shape[0] > 0 else 0

        if ker_Tv_embedded.shape[0] > 0:
            if im_d_T_pp1.shape[0] > 0:
                combined = np.vstack([im_d_T_pp1, ker_Tv_embedded]) % PRIME
                rk_comb = int(_gauss_rank_np(combined.copy(), PRIME))
                istar_ranks[p] = rk_comb - rk_im_T
            else:
                istar_ranks[p] = int(_gauss_rank_np(ker_Tv_embedded.copy(), PRIME))
        else:
            istar_ranks[p] = 0

    # Compute relative Betti from LES
    # LES: ... → H_p(T) →f_p H_p^rel →δ_p H_{p-1}(T\v) →i_{p-1} H_{p-1}(T) → ...
    # From high to low:
    # At p = max_p+1: H_{max_p+1}^rel = 0 (assumption)
    # For each p: 0 → coker(f_{p+1}) → H_p^rel → ker(i_{p-1}) → 0
    # where coker(f_{p+1}) = H_p(T) / im(f_{p+1})
    # and ker(i_{p-1}) = β_{p-1}(T\v) - rk(i_{p-1})

    # Actually from exactness:
    # β_p^rel = dim(ker δ_p) + dim(im δ_p)
    # ker δ_p = im f_p, dim = β_p(T) - dim(ker f_p) = β_p(T) - dim(im δ_{p+1})
    # im δ_p = ker i_{p-1}, dim = β_{p-1}(T\v) - rk(i_{p-1})

    # So β_p^rel = [β_p(T) - dim(im δ_{p+1})] + [β_{p-1}(T\v) - rk(i_{p-1})]

    # Need dim(im δ_p) at each level. From im δ_p = ker i_{p-1}:
    # dim(im δ_p) = β_{p-1}(T\v) - rk(i_{p-1})

    rel_bettis = {}
    for p in range(2, max_p + 1):
        dim_im_delta_pp1 = bettis_T.get(p, 0) - istar_ranks.get(p, 0) if p <= max_p else 0
        # Wait: im(δ_{p+1}) = ker(i_p), dim = β_p(T\v) - rk(i_p)
        dim_im_delta_pp1 = bettis_Tv.get(p, 0) - istar_ranks.get(p, 0)

        dim_im_delta_p = bettis_Tv.get(p-1, 0) - istar_ranks.get(p-1, 0)

        bp_rel = (bettis_T.get(p, 0) - dim_im_delta_pp1) + dim_im_delta_p
        # Hmm, this doesn't look right. Let me reconsider.

        # From the LES: ... H_p(T) → H_p^rel → H_{p-1}(T\v) → H_{p-1}(T) → ...
        # β_p^rel = dim(ker δ_p) + dim(im δ_p)
        # dim(ker δ_p) = dim(im f_p) where f_p: H_p(T) → H_p^rel
        # dim(im δ_p) = dim(ker i_{p-1})

        # From exactness at H_p(T): ker f_p = im i_p (where i_p comes from the previous degree in the LES)
        # Wait, the LES is:
        # δ_{p+1}: H_{p+1}^rel → H_p(T\v) →i_p H_p(T) →f_p H_p^rel →δ_p H_{p-1}(T\v) →i_{p-1} ...
        #
        # Exactness at H_p(T): ker f_p = im i_p. So dim(ker f_p) = rk(i_p).
        # dim(im f_p) = β_p(T) - rk(i_p).
        #
        # Exactness at H_p^rel: ker δ_p = im f_p. So dim(ker δ_p) = β_p(T) - rk(i_p).
        #
        # Exactness at H_{p-1}(T\v): ker i_{p-1} = im δ_p. So dim(im δ_p) = β_{p-1}(T\v) - rk(i_{p-1}).
        #
        # Therefore: β_p^rel = [β_p(T) - rk(i_p)] + [β_{p-1}(T\v) - rk(i_{p-1})]

        bp_rel = (bettis_T.get(p, 0) - istar_ranks.get(p, 0)) + (bettis_Tv.get(p-1, 0) - istar_ranks.get(p-1, 0))
        rel_bettis[p] = bp_rel

    return {
        'bettis_T': bettis_T,
        'bettis_Tv': bettis_Tv,
        'istar_ranks': istar_ranks,
        'rel_bettis': rel_bettis,
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"FULL RELATIVE HOMOLOGY AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        t0 = time.time()
        target = 100 if n <= 7 else 200

        for trial in range(80000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)

            for v in range(n):
                r = compute_istar_all_degrees(A, n, v)
                if r is None:
                    continue
                results.append(r)
                if len(results) >= target:
                    break

        elapsed = time.time() - t0
        print(f"  {len(results)} (T,v) pairs, {elapsed:.1f}s")

        # rank(i_*) at each degree
        for p in range(2, min(n, 7)):
            vals = [r['istar_ranks'].get(p, 0) for r in results]
            bTv = [r['bettis_Tv'].get(p, 0) for r in results]
            if any(v != 0 for v in vals) or any(v != 0 for v in bTv):
                inj = sum(1 for i, v in enumerate(vals) if v == bTv[i])
                print(f"  i_* at degree {p}: rank dist={dict(Counter(vals))}, β_p(T\\v) dist={dict(Counter(bTv))}, injective: {inj}/{len(results)}")

        # Relative Betti numbers
        print(f"\n  Relative Betti numbers:")
        for p in range(2, min(n, 7)):
            vals = [r['rel_bettis'].get(p, 0) for r in results]
            if any(v != 0 for v in vals):
                print(f"  β_{p}^rel: {dict(Counter(vals))}")

        # Is relative complex acyclic?
        acyclic = sum(1 for r in results if all(r['rel_bettis'].get(p, 0) == 0 for p in range(2, min(n, 7))))
        print(f"\n  Relative complex acyclic: {acyclic}/{len(results)}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
