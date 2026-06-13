"""
boundary_decomposition.py — Algebraic structure of im(d_4^new) projected to old space

KEY ALGEBRAIC OBSERVATION:
A new 5-path σ = (a₀,...,v,...,a₄) where v is at position i has boundary:
  d(σ) = Σⱼ (-1)ʲ (a₀,...,âⱼ,...,a₄)

Face j=i omits v → this is an OLD 4-path (doesn't use v).
Faces j≠i still contain v → these are NEW 4-paths.

So each new 5-path contributes exactly ONE old 4-path to the boundary.
If we project d_4^new onto the old-path subspace, each column has exactly
one nonzero entry (with sign (-1)^{position of v}).

QUESTION: Is this projection (when further projected to Omega) zero?
Or does it have rank, but the rank is already contained in im(d_4^old)?

This is the MECHANISM question for HYP-398.

Author: opus-2026-03-09-S56
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament,
    enumerate_all_allowed, boundary_faces,
    _gauss_rank_np, _gauss_nullbasis_modp,
    _build_constraint_matrix,
    full_chain_complex_modp,
    RANK_PRIME
)

PRIME = RANK_PRIME


def get_omega_basis(ap, deg, prime):
    paths = ap.get(deg, [])
    if not paths:
        return None
    P, nr, nc = _build_constraint_matrix(ap, deg, prime)
    if P is not None:
        r, nb = _gauss_nullbasis_modp(P, nr, nc, prime)
        return np.array(nb, dtype=np.int64) if nb else None
    else:
        return np.eye(len(paths), dtype=np.int64)


def analyze_decomposition(A, n, v, verbose=False):
    """Decompose d_4 image into old/new contributions in A-path coordinates."""
    max_p = n - 1
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]

    cc_T = full_chain_complex_modp(A, n, max_p)
    cc_Tv = full_chain_complex_modp(A_sub, n1, min(max_p, n1 - 1))

    if cc_T['bettis'].get(3, 0) < 1 or cc_Tv['bettis'].get(3, 0) < 1:
        return None

    ap_T = enumerate_all_allowed(A, n, max_p)

    paths_T_3 = ap_T.get(3, [])
    paths_T_4 = ap_T.get(4, [])

    if not paths_T_4:
        return None

    idx_T3 = {tuple(q): i for i, q in enumerate(paths_T_3)}

    # Classify 4-paths and 5-paths as old (no v) or new (uses v)
    old_3idx = [i for i, q in enumerate(paths_T_3) if v not in q]
    new_3idx = [i for i, q in enumerate(paths_T_3) if v in q]
    old_4idx = [i for i, q in enumerate(paths_T_4) if v not in q]
    new_4idx = [i for i, q in enumerate(paths_T_4) if v in q]

    # Full boundary matrix d_4: A_3 x A_4
    bd = np.zeros((len(paths_T_3), len(paths_T_4)), dtype=np.int64)
    for j, sigma in enumerate(paths_T_4):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx_T3:
                bd[idx_T3[ft], j] = (bd[idx_T3[ft], j] + sign) % PRIME

    # Extract blocks:
    # bd_new_to_old: old A_3 coords of boundaries of new A_4 paths
    # bd_new_to_new: new A_3 coords of boundaries of new A_4 paths
    bd_new_to_old = bd[np.ix_(old_3idx, new_4idx)] % PRIME
    bd_new_to_new = bd[np.ix_(new_3idx, new_4idx)] % PRIME

    # Each new 5-path should have exactly 1 old face
    nnz_per_col = np.sum(bd_new_to_old != 0, axis=0)

    # Rank of the old projection
    rank_old_proj = _gauss_rank_np(bd_new_to_old.copy(), PRIME)

    # Now project through Omega: does the Omega quotient kill some of this?
    ob_T3 = get_omega_basis(ap_T, 3, PRIME)
    ob_T4 = get_omega_basis(ap_T, 4, PRIME)

    if ob_T4 is None or ob_T3 is None:
        return None

    # Which Omega_4 basis vectors are new?
    new_omega4 = []
    old_omega4 = []
    for i in range(ob_T4.shape[0]):
        uses_v = any(ob_T4[i, j] % PRIME != 0 and v in paths_T_4[j]
                     for j in range(len(paths_T_4)))
        if uses_v:
            new_omega4.append(i)
        else:
            old_omega4.append(i)

    # d_4 in A-path coords, using Omega_4 basis
    # im(d_4^T) in A_3 coords = bd @ ob_T4.T
    im_total = bd @ ob_T4.T % PRIME

    # im(d_4^new) = bd @ ob_T4[new_omega4].T
    if new_omega4:
        im_new_omega = bd @ ob_T4[new_omega4].T % PRIME
        # Project to old 4-paths only
        im_new_omega_old = im_new_omega[old_3idx] % PRIME
        rank_omega_old_proj = _gauss_rank_np(im_new_omega_old.copy(), PRIME)
    else:
        rank_omega_old_proj = 0

    # Also: rank of old Omega_4 image in old A_3
    if old_omega4:
        im_old_omega = bd @ ob_T4[old_omega4].T % PRIME
        im_old_omega_old = im_old_omega[old_3idx] % PRIME
        rank_old_omega_old_proj = _gauss_rank_np(im_old_omega_old.copy(), PRIME)
    else:
        rank_old_omega_old_proj = 0

    # Combined rank of [old + new] Omega images in old A_3
    if old_omega4 and new_omega4:
        combined = np.hstack([im_old_omega[old_3idx], im_new_omega[old_3idx]]) % PRIME
        rank_combined_old = _gauss_rank_np(combined.copy(), PRIME)
    elif old_omega4:
        rank_combined_old = rank_old_omega_old_proj
    elif new_omega4:
        rank_combined_old = rank_omega_old_proj
    else:
        rank_combined_old = 0

    # KEY: does new add anything to old in the old-projection?
    new_adds_to_old = rank_combined_old - rank_old_omega_old_proj

    result = {
        'n_old_3': len(old_3idx),
        'n_new_3': len(new_3idx),
        'n_old_4': len(old_4idx),
        'n_new_4': len(new_4idx),
        'nnz_per_new5': list(nnz_per_col),
        'rank_raw_old_proj': rank_old_proj,
        'n_old_omega4': len(old_omega4),
        'n_new_omega4': len(new_omega4),
        'rank_omega_old_proj': rank_omega_old_proj,
        'rank_old_omega_old_proj': rank_old_omega_old_proj,
        'rank_combined_old': rank_combined_old,
        'new_adds_to_old': new_adds_to_old,
        'b3_T': cc_T['bettis'].get(3, 0),
        'b3_Tv': cc_Tv['bettis'].get(3, 0),
    }

    if verbose:
        print(f"    A_3: {len(old_3idx)} old + {len(new_3idx)} new = {len(paths_T_3)}")
        print(f"    A_4: {len(old_4idx)} old + {len(new_4idx)} new = {len(paths_T_4)}")
        print(f"    Omega_4: {len(old_omega4)} old + {len(new_omega4)} new")
        print(f"    Raw old proj of new A_4 boundaries: rank = {rank_old_proj}")
        nnz_vals = Counter(nnz_per_col)
        print(f"    #old faces per new 5-path: {dict(sorted(nnz_vals.items()))}")
        print(f"    Omega old proj of new Omega_4: rank = {rank_omega_old_proj}")
        print(f"    Omega old proj of old Omega_4: rank = {rank_old_omega_old_proj}")
        print(f"    Combined old proj: rank = {rank_combined_old}")
        print(f"    New adds to old projection: {new_adds_to_old}")

    return result


def main():
    print("=" * 70)
    print("BOUNDARY DECOMPOSITION — OLD/NEW FACE STRUCTURE")
    print("=" * 70)

    n = 7
    rng = np.random.RandomState(42)

    found = 0
    target = 40
    t0 = time.time()
    results = []

    while found < target:
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, n - 1)
        if cc['bettis'].get(3, 0) != 1:
            continue
        found += 1

        for vv in range(n):
            r = analyze_decomposition(A, n, vv, verbose=(found <= 3))
            if r is not None:
                results.append(r)

        if found % 10 == 0:
            elapsed = time.time() - t0
            print(f"  {found}/{target}, {len(results)} bad verts, {elapsed:.1f}s")

    elapsed = time.time() - t0
    print(f"\n  Done: {found} tours, {len(results)} bad verts, {elapsed:.1f}s")

    print(f"\n{'='*70}")
    print("RESULTS")
    print("=" * 70)

    # How many old faces does each new 5-path contribute?
    all_nnz = []
    for r in results:
        all_nnz.extend(r['nnz_per_new5'])
    nnz_dist = Counter(all_nnz)
    print(f"\n  #old faces per new A_4 path: {dict(sorted(nnz_dist.items()))}")

    # Raw old projection rank
    raw_dist = Counter(r['rank_raw_old_proj'] for r in results)
    print(f"\n  Raw old projection rank: {dict(sorted(raw_dist.items()))}")

    # Omega old projection of new Omega_4
    omega_new_dist = Counter(r['rank_omega_old_proj'] for r in results)
    print(f"  Omega old proj of new: {dict(sorted(omega_new_dist.items()))}")

    # Does new add to old?
    adds_dist = Counter(r['new_adds_to_old'] for r in results)
    print(f"  New adds to old (old proj): {dict(sorted(adds_dist.items()))}")

    if all(r['new_adds_to_old'] == 0 for r in results):
        print(f"\n  *** NEW Omega_4 boundaries add NOTHING to old-projection ***")
        print(f"  This means: im(d_4^new) projected to old A_3 is CONTAINED")
        print(f"  in im(d_4^old) projected to old A_3.")
        print(f"  STRONGER than HYP-398: not just at ker_d3 level, but at A_3 level!")
    else:
        # Check if it's contained when we consider ker_d3
        nonzero = sum(1 for r in results if r['new_adds_to_old'] != 0)
        print(f"\n  New adds to old in {nonzero}/{len(results)} cases")
        print(f"  (HYP-398 could still hold at ker_d3 level)")


if __name__ == '__main__':
    main()
    print("\nDONE.")
