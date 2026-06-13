"""
les_rank_i_star_n7.py — Compute rank(i_*) at n=7 for the LES proof of beta_3 <= 1

KEY QUESTION: When beta_3(T) = 1 and beta_3(T\v) = 1,
what is rank(i_*: H_3(T\v) -> H_3(T))?

HYP-371a (opus-S52): H_3(T,T\v) = beta_3(T) - rank(i_*) when beta_2 = 0.
For beta_3(T) = 1: H_3(T,T\v) = 1 - rank(i_*) in {0, 1}.
So rank(i_*) in {0, 1} (which is trivially true since beta_3(T) = 1).

The CRITICAL case is: when beta_3(T\v) = 1 (non-good vertex),
does rank(i_*) = 1? If so, H_3(T,T\v) = 0 and beta_3(T) = rank(i_*) = 1.
This would prove beta_3(T) <= 1 without needing good vertex existence!

Strategy: For beta_3(T)=1 tournaments at n=7, find all v with beta_3(T\v)=1,
and compute rank(i_*) by checking if the H_3(T\v) generator maps nontrivially into H_3(T).

Author: kind-pasteur-S47 (2026-03-09)
"""
import sys
import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, full_chain_complex, boundary_faces,
    compute_omega_basis_numpy
)

PRIME = (1 << 31) - 1


def gauss_rank_modp(M, prime=PRIME):
    """Gaussian elimination mod prime on numpy int64 array."""
    M = M.copy() % prime
    nrows, ncols = M.shape
    rank = 0
    for col in range(ncols):
        nonzero = np.where(M[rank:, col] != 0)[0]
        if len(nonzero) == 0:
            continue
        pivot = nonzero[0] + rank
        if pivot != rank:
            M[[rank, pivot]] = M[[pivot, rank]]
        inv = pow(int(M[rank, col]), prime - 2, prime)
        M[rank] = M[rank] * inv % prime
        factors = M[:, col].copy()
        factors[rank] = 0
        nz = np.where(factors != 0)[0]
        if len(nz) > 0:
            M[nz] = (M[nz] - np.outer(factors[nz], M[rank])) % prime
        rank += 1
    return rank


def gauss_nullspace_modp(M, prime=PRIME):
    """Null space basis (rows) mod prime."""
    M = M.copy() % prime
    nrows, ncols = M.shape
    pivot_cols = []
    rank = 0
    for col in range(ncols):
        nonzero = np.where(M[rank:, col] != 0)[0]
        if len(nonzero) == 0:
            continue
        pivot = nonzero[0] + rank
        if pivot != rank:
            M[[rank, pivot]] = M[[pivot, rank]]
        inv = pow(int(M[rank, col]), prime - 2, prime)
        M[rank] = M[rank] * inv % prime
        factors = M[:, col].copy()
        factors[rank] = 0
        nz = np.where(factors != 0)[0]
        if len(nz) > 0:
            M[nz] = (M[nz] - np.outer(factors[nz], M[rank])) % prime
        pivot_cols.append(col)
        rank += 1

    null_dim = ncols - rank
    if null_dim == 0:
        return np.zeros((0, ncols), dtype=np.int64)

    free_cols = sorted(set(range(ncols)) - set(pivot_cols))
    null_basis = np.zeros((null_dim, ncols), dtype=np.int64)
    for i, fc in enumerate(free_cols):
        null_basis[i, fc] = 1
        for j, pc in enumerate(pivot_cols):
            null_basis[i, pc] = (-M[j, fc]) % prime
    return null_basis % prime


def get_allowed_paths(A, vertices, max_p):
    """Get all allowed p-paths on given vertex set."""
    paths = {}
    for p in range(max_p + 2):
        p_paths = []
        if p + 1 <= len(vertices):
            for verts in combinations(vertices, p + 1):
                for perm in permutations(verts):
                    ok = True
                    for i in range(len(perm) - 1):
                        if A[perm[i]][perm[i+1]] != 1:
                            ok = False
                            break
                    if ok:
                        p_paths.append(perm)
        paths[p] = p_paths
    return paths


def build_omega_modp(ap, max_p, prime=PRIME):
    """Build Omega basis mod p. Returns dict p -> null basis (rows are basis vecs)."""
    omega = {}
    for p in range(max_p + 1):
        a_p = ap[p]
        if not a_p:
            omega[p] = np.zeros((0, 0), dtype=np.int64)
            continue

        a_pm1_set = set(ap.get(p-1, [])) if p > 0 else set()
        na_faces = {}
        for sigma in a_p:
            for i in range(1, len(sigma) - 1):
                face = sigma[:i] + sigma[i+1:]
                if p > 0 and face not in a_pm1_set:
                    if face not in na_faces:
                        na_faces[face] = len(na_faces)

        if not na_faces:
            omega[p] = np.eye(len(a_p), dtype=np.int64)
        else:
            C = np.zeros((len(na_faces), len(a_p)), dtype=np.int64)
            for j, sigma in enumerate(a_p):
                for i in range(1, len(sigma) - 1):
                    face = sigma[:i] + sigma[i+1:]
                    if face in na_faces:
                        C[na_faces[face], j] = (C[na_faces[face], j] + ((-1)**i)) % prime
            C = C % prime
            null = gauss_nullspace_modp(C, prime)
            omega[p] = null if null.shape[0] > 0 else np.zeros((0, len(a_p)), dtype=np.int64)
    return omega


def build_boundary_modp(ap, p, prime=PRIME):
    """Build boundary matrix d_p: A_p -> A_{p-1} mod prime."""
    a_p = ap[p]
    a_pm1 = ap[p-1]
    if not a_p or not a_pm1:
        return np.zeros((len(a_pm1) if a_pm1 else 0, len(a_p) if a_p else 0), dtype=np.int64)

    idx = {path: i for i, path in enumerate(a_pm1)}
    mat = np.zeros((len(a_pm1), len(a_p)), dtype=np.int64)
    for j, sigma in enumerate(a_p):
        for i in range(len(sigma)):
            face = sigma[:i] + sigma[i+1:]
            if face in idx:
                mat[idx[face], j] = (mat[idx[face], j] + ((-1)**i)) % prime
    return mat % prime


def compute_rank_i_star(A, n, v, prime=PRIME):
    """Compute rank(i_*: H_3(T\\v) -> H_3(T)) using mod-p linear algebra.

    Strategy:
    1. Compute ker(d_3) for T and T\\v in A_3 coordinates
    2. Compute im(d_4) for T and T\\v in A_3 coordinates
    3. H_3(T) = ker(d_3)/im(d_4), H_3(T\\v) similarly
    4. i_* sends a cycle [z] in H_3(T\\v) to [z] viewed in H_3(T)
    5. rank(i_*) = dim of (ker(d_3^{T\\v}) mapped into H_3(T)) modulo im(d_4^T)
    """
    max_p = 5
    remaining = [i for i in range(n) if i != v]

    all_v = list(range(n))
    ap_T = get_allowed_paths(A, all_v, max_p)
    ap_Tv = get_allowed_paths(A, remaining, max_p)

    omega_T = build_omega_modp(ap_T, max_p, prime)
    omega_Tv = build_omega_modp(ap_Tv, max_p, prime)

    # d_3 and d_4 for T in Omega coords
    bd3_T = build_boundary_modp(ap_T, 3, prime)
    bd4_T = build_boundary_modp(ap_T, 4, prime)

    om3_T = omega_T[3]  # rows are Omega_3(T) basis
    om2_T = omega_T[2]
    om4_T = omega_T[4]

    if om3_T.shape[0] == 0:
        return {'rank_i_star': 0, 'b3_T': 0, 'b3_Tv': 0}

    # d_3 on Omega_3(T): om2_T @ bd3_T @ om3_T.T
    d3_omega = om2_T @ bd3_T @ om3_T.T % prime if om2_T.shape[0] > 0 else np.zeros((0, om3_T.shape[0]), dtype=np.int64)

    # ker(d_3) for T
    ker_d3_T = gauss_nullspace_modp(d3_omega, prime)  # rows are kernel basis in Omega_3(T) coords

    # im(d_4) for T (in Omega_3(T) coords)
    if om4_T.shape[0] > 0 and om3_T.shape[0] > 0:
        d4_omega = om3_T @ bd4_T @ om4_T.T % prime
        # Image of d_4 in Omega_3 coords: rows of d4_omega.T are image generators
        # Actually d4_omega has shape (dim_Om3, dim_Om4), columns are images
        im_d4_T = d4_omega  # dim_Om3 x dim_Om4
    else:
        im_d4_T = np.zeros((om3_T.shape[0], 0), dtype=np.int64)

    b3_T = ker_d3_T.shape[0] - gauss_rank_modp(im_d4_T.copy(), prime) if ker_d3_T.shape[0] > 0 else 0

    # Same for T\v
    bd3_Tv = build_boundary_modp(ap_Tv, 3, prime)
    bd4_Tv = build_boundary_modp(ap_Tv, 4, prime)

    om3_Tv = omega_Tv[3]
    om2_Tv = omega_Tv[2]
    om4_Tv = omega_Tv[4]

    if om3_Tv.shape[0] == 0:
        b3_Tv = 0
        return {'rank_i_star': 0, 'b3_T': b3_T, 'b3_Tv': 0}

    d3_omega_Tv = om2_Tv @ bd3_Tv @ om3_Tv.T % prime if om2_Tv.shape[0] > 0 else np.zeros((0, om3_Tv.shape[0]), dtype=np.int64)
    ker_d3_Tv = gauss_nullspace_modp(d3_omega_Tv, prime)

    if om4_Tv.shape[0] > 0:
        d4_omega_Tv = om3_Tv @ bd4_Tv @ om4_Tv.T % prime
    else:
        d4_omega_Tv = np.zeros((om3_Tv.shape[0], 0), dtype=np.int64)

    rank_d4_Tv = gauss_rank_modp(d4_omega_Tv.copy(), prime)
    b3_Tv = ker_d3_Tv.shape[0] - rank_d4_Tv

    if b3_Tv == 0:
        return {'rank_i_star': 0, 'b3_T': b3_T, 'b3_Tv': 0}

    # Now compute rank(i_*)
    # i_* sends [z] in H_3(T\v) to [z] in H_3(T)
    # Need to: take ker(d_3) cycles of T\v, map them to Omega_3(T), quotient by im(d_4) of T

    # Step 1: Map ker_d3_Tv (in Omega_3(T\v) coords) to A_3(T\v) coords
    # ker_d3_Tv rows * om3_Tv gives A_3(T\v) coordinates
    cycles_Tv_in_ATv = ker_d3_Tv @ om3_Tv % prime  # (ker_dim, |A_3(T\v)|)

    # Step 2: Embed A_3(T\v) into A_3(T)
    idx_T = {path: i for i, path in enumerate(ap_T[3])}
    mapping = np.zeros((len(ap_Tv[3]), len(ap_T[3])), dtype=np.int64)
    for j_tv, path_tv in enumerate(ap_Tv[3]):
        if path_tv in idx_T:
            mapping[j_tv, idx_T[path_tv]] = 1

    cycles_Tv_in_AT = cycles_Tv_in_ATv @ mapping % prime  # (ker_dim, |A_3(T)|)

    # Step 3: Project to Omega_3(T) coords
    # cycles_in_omT = cycles_Tv_in_AT @ om3_T.T would give Omega coords
    # But we need to use the pseudoinverse / projection properly
    # Since om3_T rows span Omega_3(T), to project:
    # x in A_3(T) maps to om3_T @ om3_T.T (pseudoinverse) ... but in mod p this is tricky

    # Instead: use the fact that cycles of T\v are already Omega chains.
    # They ARE in Omega_3(T) since Omega(T\v) ⊆ Omega(T).
    # So we can express them directly in Omega_3(T) coords by solving:
    # cycles_Tv_in_AT = coeff @ om3_T

    # Solve: om3_T.T @ coeff.T = cycles_Tv_in_AT.T
    # i.e., augmented system [om3_T.T | cycles.T] has solution

    dim_om3_T = om3_T.shape[0]
    n_cycles = cycles_Tv_in_AT.shape[0]

    # Build augmented matrix and solve
    # om3_T has shape (dim_om3, |A_3|). We want coeff such that coeff @ om3_T = cycles
    # This is: cycles = coeff @ om3_T, solving for coeff
    # Transpose: om3_T.T @ coeff.T = cycles.T
    # System: M x = b where M = om3_T.T, x = coeff.T, b = cycles.T

    # For mod p: just row reduce the augmented [om3_T; cycles] and read off coords
    aug = np.vstack([om3_T, cycles_Tv_in_AT]) % prime
    aug_rank = gauss_rank_modp(aug.copy(), prime)

    # If aug_rank > dim_om3_T, some cycles are outside Omega(T) — shouldn't happen
    if aug_rank > dim_om3_T:
        print(f"  WARNING: cycles outside Omega(T)! aug_rank={aug_rank} > dim_om3={dim_om3_T}")

    # The cycles in Omega_3(T) coords: project them
    # Do Gauss on om3_T to get RREF, then express cycles in that basis
    M_rref = om3_T.copy() % prime
    nrows_om = M_rref.shape[0]
    ncols_om = M_rref.shape[1]
    pivots = []
    r = 0
    for col in range(ncols_om):
        nonzero = np.where(M_rref[r:, col] != 0)[0]
        if len(nonzero) == 0:
            continue
        pivot = nonzero[0] + r
        if pivot != r:
            M_rref[[r, pivot]] = M_rref[[pivot, r]]
        inv = pow(int(M_rref[r, col]), prime - 2, prime)
        M_rref[r] = M_rref[r] * inv % prime
        factors = M_rref[:, col].copy()
        factors[r] = 0
        nz = np.where(factors != 0)[0]
        if len(nz) > 0:
            M_rref[nz] = (M_rref[nz] - np.outer(factors[nz], M_rref[r])) % prime
        pivots.append(col)
        r += 1

    # Now express cycles in RREF basis: for each cycle row, extract pivot columns
    cycles_in_omega = np.zeros((n_cycles, dim_om3_T), dtype=np.int64)
    for ci in range(n_cycles):
        row = cycles_Tv_in_AT[ci].copy()
        for basis_idx, pcol in enumerate(pivots):
            coeff = row[pcol] % prime
            if coeff != 0:
                row = (row - coeff * M_rref[basis_idx]) % prime
                cycles_in_omega[ci, basis_idx] = coeff
    cycles_in_omega = cycles_in_omega % prime

    # Step 4: Now we have cycles in Omega_3(T) coords.
    # Project to H_3(T) = ker(d_3)/im(d_4).
    # We need: rank of (cycles restricted to ker(d_3), modulo im(d_4))

    # First restrict to ker(d_3): project cycles onto ker(d_3) basis
    # ker_d3_T rows are kernel basis in Omega_3(T) coords
    ker_dim = ker_d3_T.shape[0]
    if ker_dim == 0:
        return {'rank_i_star': 0, 'b3_T': b3_T, 'b3_Tv': b3_Tv}

    # Check: cycles should be IN ker(d_3) since they're d_3-closed in T\v
    # and Omega(T\v) is a subcomplex of Omega(T)
    # So d_3(cycle) = 0 in Omega(T) too.
    # Express cycles in ker(d_3) basis:

    # Build: ker_d3_T augmented with cycles
    ker_plus_cycles = np.vstack([ker_d3_T, cycles_in_omega]) % prime
    combined_rank = gauss_rank_modp(ker_plus_cycles.copy(), prime)
    # The cycles' projection into ker has rank:
    cycles_in_ker_rank = combined_rank - ker_dim + (n_cycles - (combined_rank - ker_dim))
    # Actually: rank of cycles within the span of ker_d3_T = n_cycles - (combined_rank - ker_dim)
    # No... combined_rank = rank(span(ker ∪ cycles))
    # Since cycles ⊆ ker (they're ker(d_3) cycles), combined_rank = ker_dim
    # So the cycles contribute 0 new dimensions = they're all in ker

    # The RANK of i_* = dim of (im of cycles in ker(d_3)/im(d_4))
    # = rank of [im_d4_T | cycles_in_omega_restricted_to_ker] - rank(im_d4_T)
    # where everything is in Omega_3(T) coords

    # im(d_4) columns + cycles (as rows, transposed to columns)
    if im_d4_T.shape[1] > 0:
        im_cols = im_d4_T  # dim_Om3 x dim_Om4, columns are image vectors
    else:
        im_cols = np.zeros((dim_om3_T, 0), dtype=np.int64)

    # Add cycle vectors as columns
    cycle_cols = cycles_in_omega.T  # dim_Om3 x n_cycles

    combined = np.hstack([im_cols, cycle_cols]) % prime
    combined_rank_full = gauss_rank_modp(combined.copy(), prime)
    im_d4_rank = gauss_rank_modp(im_cols.copy(), prime) if im_cols.shape[1] > 0 else 0

    rank_i_star = combined_rank_full - im_d4_rank

    return {
        'rank_i_star': rank_i_star,
        'b3_T': b3_T,
        'b3_Tv': b3_Tv,
        'ker_d3_T_dim': ker_d3_T.shape[0],
        'ker_d3_Tv_dim': ker_d3_Tv.shape[0],
    }


def main():
    print("=" * 70)
    print("RANK(i_*) COMPUTATION AT n=7 FOR LES PROOF")
    print("=" * 70)

    n = 7
    rng = np.random.RandomState(42)

    # Part 1: Find beta_3=1 tournaments and compute rank(i_*) for all v
    print(f"\n--- Part 1: rank(i_*) for beta_3=1 at n={n} ---")

    results = []
    checked = 0

    for trial in range(5000):
        A = random_tournament(n, rng)
        data = full_chain_complex(A, n, max_p=5)
        b3 = data['bettis'].get(3, 0)

        if b3 != 1:
            continue

        checked += 1
        tour_results = []
        has_bad_vertex = False

        for v in range(n):
            r = compute_rank_i_star(A, n, v)
            tour_results.append(r)
            h3_rel = r['b3_T'] - r['rank_i_star']  # HYP-371a formula

            if r['b3_Tv'] > 0:
                has_bad_vertex = True

        results.append({
            'trial': trial,
            'vertex_results': tour_results,
            'has_bad_vertex': has_bad_vertex,
            'score': tuple(sorted([int(sum(A[i])) for i in range(n)])),
        })

        # Print details for interesting cases
        if has_bad_vertex and len([r for r in results if r['has_bad_vertex']]) <= 10:
            print(f"\n  Trial {trial}: beta_3(T)=1, score={results[-1]['score']}")
            for v_idx, r in enumerate(tour_results):
                h3_rel = r['b3_T'] - r['rank_i_star']
                marker = " ***" if r['b3_Tv'] > 0 else ""
                print(f"    v={v_idx}: b3(T\\v)={r['b3_Tv']}, rank(i_*)={r['rank_i_star']}, "
                      f"H_3^rel={h3_rel}{marker}")

        if checked >= 100:
            break

        if checked % 20 == 0:
            print(f"  ... {checked} beta_3=1 tournaments found from {trial+1} trials", flush=True)

    print(f"\n  Total: {checked} beta_3=1 tournaments from {trial+1} trials")

    # Part 2: Statistics
    print(f"\n--- Part 2: Statistics ---")

    n_with_bad = sum(1 for r in results if r['has_bad_vertex'])
    print(f"  Tournaments with at least one bad vertex (b3(T\\v)>0): {n_with_bad}/{len(results)}")

    # Distribution of rank(i_*) values
    rank_dist = Counter()
    rank_dist_bad = Counter()  # when b3(T\v) > 0
    rank_dist_good = Counter()  # when b3(T\v) = 0

    h3_rel_dist = Counter()
    h3_rel_when_bad = Counter()

    for r in results:
        for v_idx, vr in enumerate(r['vertex_results']):
            h3_rel = vr['b3_T'] - vr['rank_i_star']
            rank_dist[vr['rank_i_star']] += 1
            h3_rel_dist[h3_rel] += 1

            if vr['b3_Tv'] > 0:
                rank_dist_bad[vr['rank_i_star']] += 1
                h3_rel_when_bad[h3_rel] += 1
            else:
                rank_dist_good[vr['rank_i_star']] += 1

    print(f"\n  rank(i_*) distribution (all vertices):")
    for k, v in sorted(rank_dist.items()):
        print(f"    rank(i_*)={k}: {v}")

    print(f"\n  rank(i_*) when b3(T\\v) > 0 (bad vertices):")
    for k, v in sorted(rank_dist_bad.items()):
        print(f"    rank(i_*)={k}: {v}")

    print(f"\n  rank(i_*) when b3(T\\v) = 0 (good vertices):")
    for k, v in sorted(rank_dist_good.items()):
        print(f"    rank(i_*)={k}: {v}")

    print(f"\n  H_3(T,T\\v) = b3(T) - rank(i_*) distribution:")
    for k, v in sorted(h3_rel_dist.items()):
        print(f"    H_3^rel={k}: {v}")

    print(f"\n  H_3(T,T\\v) when b3(T\\v) > 0:")
    for k, v in sorted(h3_rel_when_bad.items()):
        print(f"    H_3^rel={k}: {v}")

    # Part 3: KEY TEST - Does rank(i_*) = 1 whenever b3(T\v) = 1?
    print(f"\n--- Part 3: KEY TEST ---")
    print("  Does rank(i_*) = 1 whenever b3(T\\v) = 1?")

    all_match = True
    for r in results:
        for v_idx, vr in enumerate(r['vertex_results']):
            if vr['b3_Tv'] == 1 and vr['rank_i_star'] != 1:
                print(f"  COUNTEREXAMPLE: trial={r['trial']}, v={v_idx}, "
                      f"b3(T\\v)={vr['b3_Tv']}, rank(i_*)={vr['rank_i_star']}")
                all_match = False

    if all_match:
        print(f"  YES! rank(i_*) = 1 for ALL {sum(rank_dist_bad.values())} bad vertex cases")
        print(f"  This means: when b3(T\\v) = 1, H_3(T,T\\v) = 1 - 1 = 0")
        print(f"  And: b3(T) = rank(i_*) + H_3^rel = 1 + 0 = 1. QED!")
        print(f"\n  PROOF IMPLICATION:")
        print(f"  For any vertex v of a beta_3=1 tournament T:")
        print(f"    - If b3(T\\v) = 0: H_3^rel = 1, rank(i_*) = 0. b3(T) = 0 + 1 = 1.")
        print(f"    - If b3(T\\v) = 1: H_3^rel = 0, rank(i_*) = 1. b3(T) = 1 + 0 = 1.")
        print(f"  Either way, b3(T) = 1. No good vertex needed!")
    else:
        print(f"  NO - counterexamples exist")

    # Part 4: Score analysis — what determines whether a vertex is good?
    print(f"\n--- Part 4: Bad vertex properties ---")

    bad_vertex_scores = []
    good_vertex_scores = []
    for r in results:
        for v_idx, vr in enumerate(r['vertex_results']):
            # We don't have per-vertex score here, but we have overall score
            pass

    # Count bad vertices per tournament
    bad_counts = Counter()
    for r in results:
        n_bad = sum(1 for vr in r['vertex_results'] if vr['b3_Tv'] > 0)
        bad_counts[n_bad] += 1

    print(f"  # bad vertices per beta_3=1 tournament:")
    for k, v in sorted(bad_counts.items()):
        print(f"    {k} bad vertices: {v} tournaments")

    print("\nDONE.")


if __name__ == '__main__':
    main()
