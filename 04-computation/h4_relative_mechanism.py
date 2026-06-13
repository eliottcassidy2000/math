"""
h4_relative_mechanism.py — WHY is H_4(T,T\v)=0 at n=7 for BAD vertices?

From the LES:  0 → H_4(T,T\v) → H_3(T\v) →[i_*] H_3(T)
(since b4(T)=b4(T\v)=0 at n=7)

So i_*-injectivity ⟺ H_4(T,T\v) = 0.

The relative chain complex:
  ... → Ω_5^rel →[d_5^rel]→ Ω_4^rel →[d_4^rel]→ Ω_3^rel → ...
where Ω_p^rel = Ω_p(T) / Ω_p(T\v)

H_4^rel = ker(d_4^rel) / im(d_5^rel)

Strategy: compute dim(Ω_p^rel) and ranks of d_p^rel for BAD vertices
at n=7 (where H_4^rel=0 always) vs n=8 (where H_4^rel=1 sometimes).

Author: opus-2026-03-09-S57
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, deletion_adj,
    enumerate_all_allowed, build_adj_lists,
    _build_constraint_matrix, _gauss_rank_np, _gauss_nullbasis_modp,
    full_chain_complex_modp, boundary_faces,
    RANK_PRIME
)

PRIME = RANK_PRIME


def compute_relative_h4(A, n, v, verbose=False):
    """Compute H_4(T, T\\v) via relative chain complex.

    Ω_p^rel = Ω_p(T) / Ω_p(T\\v)
    d_p^rel is the induced boundary map on quotients.

    H_4^rel = ker(d_4^rel) / im(d_5^rel)
    """
    max_p = min(n - 1, 6)

    # Full chain complex of T
    cc_T = full_chain_complex_modp(A, n, max_p)

    # Full chain complex of T\v
    Av, nv = deletion_adj(A, n, v)
    cc_Tv = full_chain_complex_modp(Av, nv, max_p)

    # Get Betti numbers
    b3_T = cc_T['bettis'].get(3, 0)
    b3_Tv = cc_Tv['bettis'].get(3, 0)
    b4_T = cc_T['bettis'].get(4, 0)
    b4_Tv = cc_Tv['bettis'].get(4, 0)

    if b3_Tv == 0:
        return None  # Not a BAD vertex

    # Now compute relative dimensions
    # Ω_p^rel has dim = dim(Ω_p(T)) - dim(Ω_p(T\v))
    # (This follows from T\v → T being an inclusion of simplicial complexes)

    dim_rel = {}
    for p in range(max_p + 1):
        dim_T = cc_T.get('omega_dims', {}).get(p, 0)
        dim_Tv = cc_Tv.get('omega_dims', {}).get(p, 0)
        dim_rel[p] = dim_T - dim_Tv

    # For H_4^rel, we need:
    # rank(d_4^rel) and rank(d_5^rel)
    #
    # From the short exact sequence of chain complexes:
    # 0 → Ω_*(T\v) → Ω_*(T) → Ω_*^rel → 0
    # we get:
    # rank(d_p^T) = rank(d_p^{T\v}) + rank(d_p^rel) - correction
    # Actually, it's more subtle. Let me compute directly.

    # The relative ranks can be computed from the LES:
    # ... → H_5(T) → H_5^rel → H_4(T\v) → H_4(T) → H_4^rel → H_3(T\v) → H_3(T) → ...

    # Since we know all Betti numbers, we can compute H_p^rel from the LES!
    # Exactness at each step.

    # Actually, the LES gives us CONSTRAINTS but not exact values without
    # knowing the ranks of the connecting maps.

    # The SIMPLEST approach: use the fact that
    # chi(relative) = chi(T) - chi(T\v)  (Euler characteristic is additive)

    chi_T = sum((-1)**p * cc_T['bettis'].get(p, 0) for p in range(n))
    chi_Tv = sum((-1)**p * cc_Tv['bettis'].get(p, 0) for p in range(n - 1))
    chi_rel = chi_T - chi_Tv

    # Also compute chi from Omega dimensions
    chi_omega_T = sum((-1)**p * cc_T.get('omega_dims', {}).get(p, 0) for p in range(max_p + 1))
    chi_omega_Tv = sum((-1)**p * cc_Tv.get('omega_dims', {}).get(p, 0) for p in range(max_p + 1))
    chi_omega_rel = chi_omega_T - chi_omega_Tv

    # For the concrete computation of H_4^rel, use the LES segment:
    # H_4(T) → H_4^rel → H_3(T\v) →[i_*] H_3(T) → H_3^rel → H_2(T\v)
    #
    # With b4(T) = 0: the map H_4(T) → H_4^rel is zero, so
    #   H_4^rel injects into H_3(T\v)
    # With b2(T\v) = 0: the map H_3^rel → H_2(T\v) is zero, so
    #   H_3(T) → H_3^rel is surjective
    #
    # Exactness: ker(i_*) = im(H_4^rel → H_3(T\v)) = H_4^rel (injective!)
    # So: dim H_4^rel = dim ker(i_*) = b3(T\v) - rank(i_*)
    #
    # If b3(T\v) = 1 and b3(T) = 1:
    #   rank(i_*) ∈ {0, 1}
    #   H_4^rel = 1 - rank(i_*)
    #
    # So H_4^rel = 0 ⟺ i_* injective ⟺ rank(i_*) = 1

    # Direct computation of rank(i_*): embed H_3(T\v) into H_3(T)
    # This requires computing the actual cycle spaces

    # Approach: compute Ω_3 basis for T and T\v, find ker(d_3) for each,
    # embed T\v's kernel into T's A_3, project to T's Ω_3, check if in ker(d_3^T)

    # Build allowed paths
    ap_T = enumerate_all_allowed(A, n, max_p)
    ap_Tv = enumerate_all_allowed(Av, nv, max_p)

    # Map T\v vertices to T vertices
    vmap = list(range(v)) + list(range(v + 1, n))  # T\v vertex i -> T vertex vmap[i]

    # Get Ω_3 basis for T\v (in T\v's A_3 coordinates)
    P3v, na3v, nc3v = _build_constraint_matrix(ap_Tv, 3, PRIME)
    paths_3_Tv = ap_Tv.get(3, [])
    if P3v is not None:
        r3v, nb3v = _gauss_nullbasis_modp(P3v, na3v, nc3v, PRIME)
        omega3_Tv = np.array(nb3v, dtype=np.int64) if nb3v else np.zeros((0, nc3v), dtype=np.int64)
    else:
        omega3_Tv = np.eye(len(paths_3_Tv), dtype=np.int64) if paths_3_Tv else np.zeros((0, 0), dtype=np.int64)

    dim_omega3_Tv = omega3_Tv.shape[0]

    # Get Ω_3 basis for T
    P3T, na3T, nc3T = _build_constraint_matrix(ap_T, 3, PRIME)
    paths_3_T = ap_T.get(3, [])
    if P3T is not None:
        r3T, nb3T = _gauss_nullbasis_modp(P3T, na3T, nc3T, PRIME)
        omega3_T = np.array(nb3T, dtype=np.int64) if nb3T else np.zeros((0, nc3T), dtype=np.int64)
    else:
        omega3_T = np.eye(len(paths_3_T), dtype=np.int64) if paths_3_T else np.zeros((0, 0), dtype=np.int64)

    dim_omega3_T = omega3_T.shape[0]

    # Build d_3 for T: Ω_3 → A_2
    paths_2_T = ap_T.get(2, [])
    idx2T = {p: i for i, p in enumerate(paths_2_T)}
    bd3_T = np.zeros((len(paths_2_T), len(paths_3_T)), dtype=np.int64)
    for j, path in enumerate(paths_3_T):
        for sign, face in boundary_faces(path):
            if face in idx2T:
                bd3_T[idx2T[face], j] = (bd3_T[idx2T[face], j] + sign) % PRIME

    # d_3 restricted to Ω_3: bd3_T @ omega3_T^T
    d3_omega_T = bd3_T @ omega3_T.T % PRIME

    # ker(d_3) in Ω_3(T)
    r_d3T, ker_d3T_vecs = _gauss_nullbasis_modp(
        d3_omega_T.astype(np.int32) if d3_omega_T.size > 0 else np.zeros((1, 1), dtype=np.int32),
        d3_omega_T.shape[0], d3_omega_T.shape[1], PRIME
    )
    ker_d3T = np.array(ker_d3T_vecs, dtype=np.int64) if ker_d3T_vecs else np.zeros((0, dim_omega3_T), dtype=np.int64)
    dim_ker_d3T = ker_d3T.shape[0]

    # ker(d_3) in A_3 coordinates (for T)
    if dim_ker_d3T > 0:
        ker_d3T_A3 = ker_d3T @ omega3_T % PRIME
    else:
        ker_d3T_A3 = np.zeros((0, len(paths_3_T)), dtype=np.int64)

    # Now embed T\v's ker(d_3) into T's A_3
    # Build d_3 for T\v
    paths_2_Tv = ap_Tv.get(2, [])
    idx2Tv = {p: i for i, p in enumerate(paths_2_Tv)}
    bd3_Tv = np.zeros((len(paths_2_Tv), len(paths_3_Tv)), dtype=np.int64)
    for j, path in enumerate(paths_3_Tv):
        for sign, face in boundary_faces(path):
            if face in idx2Tv:
                bd3_Tv[idx2Tv[face], j] = (bd3_Tv[idx2Tv[face], j] + sign) % PRIME

    d3_omega_Tv = bd3_Tv @ omega3_Tv.T % PRIME
    r_d3Tv, ker_d3Tv_vecs = _gauss_nullbasis_modp(
        d3_omega_Tv.astype(np.int32) if d3_omega_Tv.size > 0 else np.zeros((1, 1), dtype=np.int32),
        d3_omega_Tv.shape[0], d3_omega_Tv.shape[1], PRIME
    )
    ker_d3Tv = np.array(ker_d3Tv_vecs, dtype=np.int64) if ker_d3Tv_vecs else np.zeros((0, dim_omega3_Tv), dtype=np.int64)
    dim_ker_d3Tv = ker_d3Tv.shape[0]

    # ker(d_3^{T\v}) in T\v's A_3 coordinates
    if dim_ker_d3Tv > 0:
        ker_d3Tv_A3v = ker_d3Tv @ omega3_Tv % PRIME
    else:
        ker_d3Tv_A3v = np.zeros((0, len(paths_3_Tv)), dtype=np.int64)

    # Embed into T's A_3 coordinates
    # Map path (a,b,c,d) in T\v to (vmap[a], vmap[b], vmap[c], vmap[d]) in T
    idx3T = {p: i for i, p in enumerate(paths_3_T)}
    embed = np.zeros((len(paths_3_Tv), len(paths_3_T)), dtype=np.int64)
    for jv, path_v in enumerate(paths_3_Tv):
        mapped = tuple(vmap[x] for x in path_v)
        if mapped in idx3T:
            embed[jv, idx3T[mapped]] = 1

    # Embedded ker(d_3^{T\v}) in T's A_3
    if dim_ker_d3Tv > 0:
        embedded_ker = ker_d3Tv_A3v @ embed % PRIME
    else:
        embedded_ker = np.zeros((0, len(paths_3_T)), dtype=np.int64)

    # Now check rank of i_*: embedded_ker projected to H_3(T)
    # i.e., compute rank of embedded_ker in ker_d3T_A3 modulo im(d_4^T)

    # Build im(d_4^T) in A_3 coordinates, restricted to Ω_4
    paths_4_T = ap_T.get(4, [])
    P4T, na4T, nc4T = _build_constraint_matrix(ap_T, 4, PRIME)
    if P4T is not None and len(paths_4_T) > 0:
        r4T, nb4T = _gauss_nullbasis_modp(P4T, na4T, nc4T, PRIME)
        omega4_T = np.array(nb4T, dtype=np.int64) if nb4T else np.zeros((0, nc4T), dtype=np.int64)
    elif paths_4_T:
        omega4_T = np.eye(len(paths_4_T), dtype=np.int64)
    else:
        omega4_T = np.zeros((0, 0), dtype=np.int64)

    # d_4 boundary matrix
    idx3T_dict = {p: i for i, p in enumerate(paths_3_T)}
    bd4_T = np.zeros((len(paths_3_T), len(paths_4_T)), dtype=np.int64)
    for j, path in enumerate(paths_4_T):
        for sign, face in boundary_faces(path):
            if face in idx3T_dict:
                bd4_T[idx3T_dict[face], j] = (bd4_T[idx3T_dict[face], j] + sign) % PRIME

    # im(d_4) in A_3 = bd4_T @ omega4_T^T
    if omega4_T.shape[0] > 0:
        im_d4_T = (bd4_T @ omega4_T.T % PRIME).T  # rows = boundary vectors
    else:
        im_d4_T = np.zeros((0, len(paths_3_T)), dtype=np.int64)

    # Stack: im(d_4) + embedded_ker
    # rank(i_*) = rank([im_d4; embedded_ker]) - rank(im_d4) ... restricted to ker_d3
    # Actually simpler: project everything to the quotient ker(d_3)/im(d_4) = H_3(T)

    # Compute rank of [im_d4; embedded_ker] restricted to ker_d3
    # First get rank of im_d4 in ker_d3 (this is rank(d_4))
    # Then add embedded_ker and check if rank increases

    if im_d4_T.shape[0] > 0 and embedded_ker.shape[0] > 0:
        combined = np.vstack([im_d4_T, embedded_ker]) % PRIME
    elif embedded_ker.shape[0] > 0:
        combined = embedded_ker % PRIME
    else:
        combined = im_d4_T % PRIME if im_d4_T.shape[0] > 0 else np.zeros((0, len(paths_3_T)), dtype=np.int64)

    # We need everything in ker(d_3^T) space
    # But im(d_4) is automatically in ker(d_3) (d_3 ∘ d_4 = 0)
    # And embedded_ker is in ker(d_3^{T\v}) ⊂ ker(d_3^T)
    # (well, the embedding should preserve this, but let me verify)

    # Compute ranks directly in A_3 coordinates (the ker_d3 restriction is automatic)
    if im_d4_T.shape[0] > 0:
        rk_d4 = int(_gauss_rank_np(im_d4_T.copy().astype(np.int64), PRIME))
    else:
        rk_d4 = 0

    if combined.shape[0] > 0:
        rk_combined = int(_gauss_rank_np(combined.copy().astype(np.int64), PRIME))
    else:
        rk_combined = 0

    rank_istar = rk_combined - rk_d4
    h4_rel = dim_ker_d3Tv - rank_istar  # = b3(T\v) - rank(i_*)

    # Also compute rank(d_4) and rank(d_5) for relative
    rk_d4_T = cc_T.get('ranks', {}).get(4, 0)
    rk_d4_Tv = cc_Tv.get('ranks', {}).get(4, 0)
    rk_d5_T = cc_T.get('ranks', {}).get(5, 0)
    rk_d5_Tv = cc_Tv.get('ranks', {}).get(5, 0)

    result = {
        'b3_T': b3_T, 'b3_Tv': b3_Tv,
        'b4_T': b4_T, 'b4_Tv': b4_Tv,
        'dim_omega3_T': dim_omega3_T, 'dim_omega3_Tv': dim_omega3_Tv,
        'dim_omega4_T': omega4_T.shape[0],
        'dim_rel': dim_rel,
        'dim_ker_d3T': dim_ker_d3T,
        'dim_ker_d3Tv': dim_ker_d3Tv,
        'rk_d4_in_A3': rk_d4,
        'rk_combined': rk_combined,
        'rank_istar': rank_istar,
        'h4_rel': h4_rel,
        'chi_T': chi_T, 'chi_Tv': chi_Tv, 'chi_rel': chi_rel,
        'rk_d4_T': rk_d4_T, 'rk_d4_Tv': rk_d4_Tv,
        'rk_d5_T': rk_d5_T, 'rk_d5_Tv': rk_d5_Tv,
    }

    if verbose:
        print(f"    dim(Ω_3^rel) = {dim_rel.get(3, 0)}, dim(Ω_4^rel) = {dim_rel.get(4, 0)}, "
              f"dim(Ω_5^rel) = {dim_rel.get(5, 0)}")
        print(f"    dim_ker_d3(T)={dim_ker_d3T}, dim_ker_d3(T\\v)={dim_ker_d3Tv}")
        print(f"    rk(d_4)={rk_d4}, rk(combined)={rk_combined}, rank(i_*)={rank_istar}")
        print(f"    H_4^rel = {h4_rel}")
        print(f"    rk(d_4^T)={rk_d4_T}, rk(d_4^Tv)={rk_d4_Tv}")
        print(f"    rk(d_5^T)={rk_d5_T}, rk(d_5^Tv)={rk_d5_Tv}")

    return result


def main():
    print("=" * 70)
    print("H_4^rel MECHANISM: WHY IS H_4(T,T\\v)=0 AT n=7?")
    print("=" * 70)

    for n in [7, 8]:
        print(f"\n{'='*60}")
        print(f"n = {n}")
        print(f"{'='*60}")

        rng = np.random.RandomState(12345)

        h4_rel_dist = Counter()
        rank_istar_dist = Counter()
        dim_rel_4_vals = []
        dim_rel_5_vals = []
        found = 0
        target = 100 if n == 7 else 50
        t0 = time.time()

        for trial in range(20000):
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, n - 1)
            b3 = cc['bettis'].get(3, 0)

            if b3 != 1:
                continue

            # Check all vertices
            for v in range(n):
                Av, nv = deletion_adj(A, n, v)
                cc_v = full_chain_complex_modp(Av, nv, nv - 1)
                b3v = cc_v['bettis'].get(3, 0)

                if b3v == 0:
                    continue  # GOOD vertex, skip

                # BAD vertex
                verbose = (found < 3)
                if verbose:
                    scores = sorted([int(sum(A[i])) for i in range(n)])
                    print(f"\n  Trial {trial}, v={v}, scores={scores}")

                r = compute_relative_h4(A, n, v, verbose=verbose)
                if r is None:
                    continue

                h4_rel_dist[r['h4_rel']] += 1
                rank_istar_dist[r['rank_istar']] += 1
                dim_rel_4_vals.append(r['dim_rel'].get(4, 0))
                dim_rel_5_vals.append(r['dim_rel'].get(5, 0))

                found += 1

            if found >= target:
                break

        elapsed = time.time() - t0
        print(f"\n  n={n}: {found} BAD vertices analyzed, {elapsed:.1f}s")
        print(f"  H_4^rel distribution: {dict(sorted(h4_rel_dist.items()))}")
        print(f"  rank(i_*) distribution: {dict(sorted(rank_istar_dist.items()))}")

        if dim_rel_4_vals:
            print(f"  dim(Ω_4^rel): min={min(dim_rel_4_vals)}, max={max(dim_rel_4_vals)}, "
                  f"mean={np.mean(dim_rel_4_vals):.1f}")
            print(f"  dim(Ω_5^rel): min={min(dim_rel_5_vals)}, max={max(dim_rel_5_vals)}, "
                  f"mean={np.mean(dim_rel_5_vals):.1f}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
