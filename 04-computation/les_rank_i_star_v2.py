"""
les_rank_i_star_v2.py — FIXED computation of rank(i_*) for LES proof

Key fix: compute ranks from boundary maps in A-path coordinates,
not in Omega coordinates. The formula om2 @ bd3 @ om3.T does NOT give
the correct d_3 in Omega coords under mod-p arithmetic (om2 is not
an orthogonal projector mod p).

Correct approach:
- rank(d_p) = rank(bd_p @ omega_p_basis.T) since im(d_p) ⊆ Omega_{p-1}
- ker(d_p) computed as null space of bd_p @ omega_p_basis.T
- For i_*: compute Z = ker(d_3^{T\\v}) in A_3 coords, embed into A_3(T),
  then rank(i_*) = rank([im_d4_T | Z_embedded]) - rank(im_d4_T)

Author: kind-pasteur-S47 (2026-03-09)
"""
import sys
import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import random_tournament, full_chain_complex

PRIME = (1 << 31) - 1


def gauss_rank_modp(M, prime=PRIME):
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
    paths = {}
    vset = set(vertices)
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


def build_omega_basis_modp(ap, max_p, prime=PRIME):
    """Build Omega basis. Returns dict p -> matrix (rows are basis vectors in A_p coords)."""
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
    """Build boundary d_p: A_p -> A_{p-1} mod prime."""
    a_p = ap[p]
    a_pm1 = ap.get(p-1, [])
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


def compute_beta3_and_cycles(A, n, vertices, max_p=5, prime=PRIME):
    """Compute beta_3 and the cycle/boundary data for a tournament on given vertices.

    Returns dict with:
    - beta_3: the Betti number
    - ap3: allowed 3-paths
    - Z_A3: cycle basis in A_3 coordinates (rows)
    - im_d4_A3: image of d_4 in A_3 coordinates (columns)
    - rank_d3, rank_d4, dim_omega3
    """
    ap = get_allowed_paths(A, vertices, max_p)
    omega = build_omega_basis_modp(ap, max_p, prime)

    om3 = omega[3]  # rows = basis vectors in A_3
    om4 = omega.get(4, np.zeros((0, 0), dtype=np.int64))

    dim_om3 = om3.shape[0]
    dim_om4 = om4.shape[0]

    if dim_om3 == 0:
        return {
            'beta_3': 0, 'ap3': ap[3],
            'Z_A3': np.zeros((0, len(ap[3])), dtype=np.int64),
            'im_d4_A3': np.zeros((len(ap[3]), 0), dtype=np.int64),
            'rank_d3': 0, 'rank_d4': 0, 'dim_omega3': 0,
        }

    # d_3: Omega_3 -> A_2 via bd3 @ om3.T
    bd3 = build_boundary_modp(ap, 3, prime)
    d3_mat = bd3 @ om3.T % prime  # |A_2| x dim_Om3

    rank_d3 = gauss_rank_modp(d3_mat.copy(), prime)

    # ker(d_3) in Omega_3 coords = null space of d3_mat
    ker_d3 = gauss_nullspace_modp(d3_mat.copy(), prime)  # rows in Omega_3 coords
    ker_dim = ker_d3.shape[0]

    # d_4: Omega_4 -> A_3 via bd4 @ om4.T
    if dim_om4 > 0:
        bd4 = build_boundary_modp(ap, 4, prime)
        im_d4_A3 = bd4 @ om4.T % prime  # |A_3| x dim_Om4 (columns are image vectors)
        rank_d4 = gauss_rank_modp(im_d4_A3.copy(), prime)
    else:
        im_d4_A3 = np.zeros((len(ap[3]), 0), dtype=np.int64)
        rank_d4 = 0

    beta_3 = ker_dim - rank_d4

    # Cycle basis in A_3 coords: ker_d3 @ om3
    if ker_dim > 0:
        Z_A3 = ker_d3 @ om3 % prime  # rows are cycles in A_3 coords
    else:
        Z_A3 = np.zeros((0, len(ap[3])), dtype=np.int64)

    return {
        'beta_3': beta_3, 'ap3': ap[3],
        'Z_A3': Z_A3, 'im_d4_A3': im_d4_A3,
        'rank_d3': rank_d3, 'rank_d4': rank_d4,
        'dim_omega3': dim_om3, 'ker_dim': ker_dim,
    }


def compute_rank_i_star(A, n, v, prime=PRIME):
    """Compute rank(i_*: H_3(T\\v) -> H_3(T)).

    Strategy:
    rank(i_*) = rank([im_d4_T | Z_Tv_embedded]) - rank(im_d4_T)
    where Z_Tv_embedded = cycle basis of T\\v embedded in A_3(T).
    """
    all_v = list(range(n))
    remaining = [i for i in range(n) if i != v]

    data_T = compute_beta3_and_cycles(A, n, all_v)
    data_Tv = compute_beta3_and_cycles(A, n, remaining)

    b3_T = data_T['beta_3']
    b3_Tv = data_Tv['beta_3']

    if b3_Tv == 0:
        # rank(i_*) = 0 trivially
        return {
            'rank_i_star': 0,
            'b3_T': b3_T,
            'b3_Tv': 0,
            'h3_rel': b3_T,  # HYP-371a: H_3^rel = b3_T - 0 = b3_T
        }

    # Embed cycles of T\v into A_3(T) coordinates
    Z_Tv = data_Tv['Z_A3']  # rows in A_3(T\v) coords
    ap3_T = data_T['ap3']
    ap3_Tv = data_Tv['ap3']

    idx_T = {path: i for i, path in enumerate(ap3_T)}

    # Build embedding map: A_3(T\v) -> A_3(T)
    embed = np.zeros((len(ap3_Tv), len(ap3_T)), dtype=np.int64)
    for j, path in enumerate(ap3_Tv):
        if path in idx_T:
            embed[j, idx_T[path]] = 1

    Z_Tv_in_AT = Z_Tv @ embed % prime  # rows are cycles in A_3(T) coords

    # Compute rank(i_*) = rank([im_d4_T | Z_Tv_embedded]) - rank(im_d4_T)
    im_d4_T = data_T['im_d4_A3']  # columns in A_3(T)

    if im_d4_T.shape[1] > 0:
        combined = np.hstack([im_d4_T, Z_Tv_in_AT.T]) % prime
        rank_combined = gauss_rank_modp(combined.copy(), prime)
        rank_im_d4 = gauss_rank_modp(im_d4_T.copy(), prime)
    else:
        # No d_4 image in T — cycles of T\v are directly in H_3(T)
        rank_combined = gauss_rank_modp(Z_Tv_in_AT.T.copy(), prime) if Z_Tv_in_AT.shape[0] > 0 else 0
        rank_im_d4 = 0

    rank_i = rank_combined - rank_im_d4

    h3_rel = b3_T - rank_i  # HYP-371a formula

    return {
        'rank_i_star': rank_i,
        'b3_T': b3_T,
        'b3_Tv': b3_Tv,
        'h3_rel': h3_rel,
        'ker_d3_T': data_T['ker_dim'],
        'ker_d3_Tv': data_Tv['ker_dim'],
    }


def main():
    print("=" * 70)
    print("RANK(i_*) v2 — FIXED COMPUTATION")
    print("=" * 70)

    # Part 0: Verify at n=6 (known exhaustive results)
    print("\n--- Part 0: Verification at n=6 ---")
    from tournament_utils import bits_to_adj
    n = 6
    total = 2 ** (n*(n-1)//2)
    n6_checks = 0
    n6_ok = True

    for bits in range(total):
        A = bits_to_adj(bits, n)
        data = full_chain_complex(A, n, max_p=5)
        if data['bettis'].get(3, 0) != 1:
            continue

        for v in range(n):
            r = compute_rank_i_star(A, n, v)
            if r['b3_T'] != 1:
                print(f"  BUG: bits={bits}, v={v}, b3_T={r['b3_T']} (expected 1)")
                n6_ok = False
            if r['b3_Tv'] != 0:
                print(f"  BUG: bits={bits}, v={v}, b3_Tv={r['b3_Tv']} (expected 0)")
                n6_ok = False
            if r['rank_i_star'] != 0:
                print(f"  SURPRISE: bits={bits}, v={v}, rank(i_*)={r['rank_i_star']} (expected 0)")
            if r['h3_rel'] != 1:
                print(f"  BUG: bits={bits}, v={v}, H3_rel={r['h3_rel']} (expected 1)")
                n6_ok = False
            n6_checks += 1

        if bits % 5000 == 0 and bits > 0:
            print(f"  ... {bits}/{total}", flush=True)

    print(f"  Checked {n6_checks} (T,v) pairs at n=6: {'ALL OK' if n6_ok else 'BUGS FOUND'}")

    # Part 1: n=7 computation
    print("\n--- Part 1: rank(i_*) at n=7 ---")
    n = 7
    rng = np.random.RandomState(42)

    results = []
    checked = 0

    for trial in range(5000):
        A = random_tournament(n, rng)
        data = full_chain_complex(A, n, max_p=5)
        if data['bettis'].get(3, 0) != 1:
            continue

        checked += 1
        tour_results = []
        has_bad = False

        for v in range(n):
            r = compute_rank_i_star(A, n, v)
            tour_results.append(r)
            if r['b3_Tv'] > 0:
                has_bad = True

        results.append({
            'trial': trial,
            'vertex_results': tour_results,
            'has_bad': has_bad,
            'score': tuple(sorted([int(sum(A[i])) for i in range(n)])),
        })

        if has_bad and sum(1 for r in results if r['has_bad']) <= 5:
            print(f"\n  Trial {trial}: b3(T)=1, score={results[-1]['score']}")
            for vi, r in enumerate(tour_results):
                marker = " <-- BAD" if r['b3_Tv'] > 0 else ""
                print(f"    v={vi}: b3(T\\v)={r['b3_Tv']}, rank(i_*)={r['rank_i_star']}, "
                      f"H_3^rel={r['h3_rel']}{marker}")

        if checked >= 100:
            break

        if checked % 20 == 0:
            print(f"  ... {checked} beta_3=1 found from {trial+1} trials", flush=True)

    print(f"\n  Total: {checked} beta_3=1 from {trial+1} trials")

    # Part 2: Statistics
    print("\n--- Part 2: Statistics ---")

    n_bad = sum(1 for r in results if r['has_bad'])
    print(f"  Tours with bad vertex: {n_bad}/{len(results)}")

    rank_dist = Counter()
    rank_when_bad = Counter()
    h3_rel_dist = Counter()

    for r in results:
        for vr in r['vertex_results']:
            rank_dist[vr['rank_i_star']] += 1
            h3_rel_dist[vr['h3_rel']] += 1
            if vr['b3_Tv'] > 0:
                rank_when_bad[vr['rank_i_star']] += 1

    print(f"\n  rank(i_*) overall: {dict(sorted(rank_dist.items()))}")
    print(f"  rank(i_*) when b3(T\\v) > 0: {dict(sorted(rank_when_bad.items()))}")
    print(f"  H_3^rel distribution: {dict(sorted(h3_rel_dist.items()))}")

    # Part 3: KEY TEST
    print("\n--- Part 3: KEY TEST ---")
    print("  When b3(T\\v) = 1, is rank(i_*) = 1?")

    violations = []
    for r in results:
        for vi, vr in enumerate(r['vertex_results']):
            if vr['b3_Tv'] == 1 and vr['rank_i_star'] != 1:
                violations.append((r['trial'], vi, vr))

    if not violations:
        total_bad = sum(rank_when_bad.values())
        print(f"  YES! All {total_bad} bad-vertex cases have rank(i_*) = 1")
        print(f"  => When b3(T\\v) = 1: H_3^rel = 1 - 1 = 0")
        print(f"  => b3(T) = 1 = rank(i_*) + 0")
        print()
        print("  THEOREM (computational): For beta_3=1 tournaments at n=7:")
        print("  - If b3(T\\v) = 0: rank(i_*) = 0, H_3^rel = 1")
        print("  - If b3(T\\v) = 1: rank(i_*) = 1, H_3^rel = 0")
        print("  The LES always gives b3(T) = 1, without needing good vertex!")
    else:
        print(f"  NO! {len(violations)} counterexamples:")
        for trial, vi, vr in violations[:5]:
            print(f"    trial={trial}, v={vi}: b3_Tv={vr['b3_Tv']}, rank(i_*)={vr['rank_i_star']}")

    # Part 4: Check H_3^rel non-negativity (sanity)
    print("\n--- Part 4: Sanity check ---")
    neg_count = sum(1 for r in results for vr in r['vertex_results'] if vr['h3_rel'] < 0)
    print(f"  H_3^rel < 0 cases: {neg_count}")
    if neg_count > 0:
        print("  WARNING: negative H_3^rel indicates a computational bug!")
        for r in results:
            for vi, vr in enumerate(r['vertex_results']):
                if vr['h3_rel'] < 0:
                    print(f"    trial={r['trial']}, v={vi}: b3_T={vr['b3_T']}, "
                          f"rank_i={vr['rank_i_star']}, H_3^rel={vr['h3_rel']}")
                    break
            if neg_count > 0:
                break

    # Part 5: Additional check — extend to n=8
    print("\n--- Part 5: Quick n=8 check ---")
    n = 8
    rng8 = np.random.RandomState(123)
    n8_checked = 0
    n8_violations = 0

    for trial in range(3000):
        A = random_tournament(n, rng8)
        data = full_chain_complex(A, n, max_p=5)
        if data['bettis'].get(3, 0) != 1:
            continue

        n8_checked += 1
        for v in range(n):
            r = compute_rank_i_star(A, n, v)
            if r['b3_Tv'] == 1 and r['rank_i_star'] != 1:
                n8_violations += 1
                print(f"  n=8 VIOLATION: trial={trial}, v={v}, rank(i_*)={r['rank_i_star']}")
            if r['h3_rel'] < 0:
                print(f"  n=8 BUG: trial={trial}, v={v}, H_3^rel={r['h3_rel']}")

        if n8_checked >= 20:
            break

        if n8_checked % 5 == 0:
            print(f"  ... {n8_checked} beta_3=1 found at n=8", flush=True)

    print(f"  n=8: checked {n8_checked} tournaments, violations: {n8_violations}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
