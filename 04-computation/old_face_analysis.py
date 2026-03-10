"""
old_face_analysis.py — Why can't old faces of new boundaries reach H_3?

For BAD vertex v, each new 5-path through v contributes ONE old 4-path face.
When we take an Omega_4 linear combination of new 5-paths, the old-face
projection becomes a linear combination of old 4-paths.

QUESTION: In which subspace of old A_3 does this projection land?
Specifically: does it land in im(d_3^{T\v}) ⊕ im(d_4^{T\v})?

If the old projection of new Omega_4 boundaries lies in
  im(d_3^{T\v}) + im(d_4^{T\v})
then it can't affect ker(d_3) \ im(d_4^{T\v}), proving HYP-398.

Also: what is the GEOMETRIC meaning of the old face?
The old face (a₀,...,â_i,...,a₄) where v was at position i means:
we have a 5-path a₀→...→v→...→a₄ in T, and removing v gives
an old 4-path. But this 4-path might NOT be an allowed path in T\v
(it could have a "gap" where v was).

Actually, it IS an allowed 4-path in T since all edges are between
non-v vertices. But it might not be in Omega_3(T\v) — it could be
killed by the NA constraint.

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


def detailed_face_analysis(A, n, v, verbose=False):
    """Understand which old faces appear and their relationship to T\v chain complex."""
    max_p = n - 1
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]

    cc_T = full_chain_complex_modp(A, n, max_p)
    cc_Tv = full_chain_complex_modp(A_sub, n1, min(max_p, n1 - 1))

    if cc_T['bettis'].get(3, 0) != 1 or cc_Tv['bettis'].get(3, 0) != 1:
        return None

    ap_T = enumerate_all_allowed(A, n, max_p)
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    def embed_path(path_tv):
        return tuple(remaining[i] for i in path_tv)

    paths_T_3 = ap_T.get(3, [])
    paths_T_4 = ap_T.get(4, [])
    paths_Tv_3 = ap_Tv.get(3, [])

    idx_T3 = {tuple(q): i for i, q in enumerate(paths_T_3)}
    idx_Tv3 = {tuple(q): i for i, q in enumerate(paths_Tv_3)}

    # Embedded old 4-paths
    embedded_Tv3 = set()
    for q in paths_Tv_3:
        embedded_Tv3.add(embed_path(q))

    # For each new 5-path, find its old face and check if it's in A_3(T\v)
    new_5paths = [q for q in paths_T_4 if v in q]

    old_face_in_Tv = 0
    old_face_not_in_Tv = 0
    no_old_face = 0

    v_position_counts = Counter()

    for sigma in new_5paths:
        t = tuple(sigma)
        v_pos = list(t).index(v)
        v_position_counts[v_pos] += 1

        # Old face: remove v
        old_face = t[:v_pos] + t[v_pos+1:]
        if old_face in idx_T3 and v not in old_face:
            if old_face in embedded_Tv3:
                old_face_in_Tv += 1
            else:
                old_face_not_in_Tv += 1
        else:
            no_old_face += 1

    # KEY: For old faces that ARE in A_3(T\v), are they in Omega_3(T\v)?
    # An old face path is in A_3(T) automatically.
    # It's in A_3(T\v) iff all its vertices are in remaining.
    # But is it in Omega_3(T\v)?

    # Check: get NA relations for T\v at degree 3
    ob_Tv3 = get_omega_basis(ap_Tv, 3, PRIME)
    ob_T3 = get_omega_basis(ap_T, 3, PRIME)
    ob_T4 = get_omega_basis(ap_T, 4, PRIME)

    # Compute d_3(T\v) and ker(d_3(T\v))
    paths_Tv_2 = ap_Tv.get(2, [])
    idx_Tv2 = {tuple(q): i for i, q in enumerate(paths_Tv_2)}
    bd_Tv3 = np.zeros((len(paths_Tv_2), len(paths_Tv_3)), dtype=np.int64)
    for j, sigma in enumerate(paths_Tv_3):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx_Tv2:
                bd_Tv3[idx_Tv2[ft], j] = (bd_Tv3[idx_Tv2[ft], j] + sign) % PRIME

    d3_Tv_omega = bd_Tv3 @ ob_Tv3.T % PRIME
    M_d3Tv = [[int(x) % PRIME for x in d3_Tv_omega[r]] for r in range(d3_Tv_omega.shape[0])]
    _, ker_d3_Tv = _gauss_nullbasis_modp(M_d3Tv, d3_Tv_omega.shape[0], d3_Tv_omega.shape[1], PRIME)
    ker_d3_Tv = np.array(ker_d3_Tv, dtype=np.int64)

    # H_3 generator of T\v: find cycle not in im(d_4^Tv)
    paths_Tv_4 = ap_Tv.get(4, [])
    ob_Tv4 = get_omega_basis(ap_Tv, 4, PRIME)

    idx_Tv3_map = {tuple(q): i for i, q in enumerate(paths_Tv_3)}
    bd_Tv4 = np.zeros((len(paths_Tv_3), len(paths_Tv_4)), dtype=np.int64)
    for j, sigma in enumerate(paths_Tv_4):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx_Tv3_map:
                bd_Tv4[idx_Tv3_map[ft], j] = (bd_Tv4[idx_Tv3_map[ft], j] + sign) % PRIME

    if ob_Tv4 is not None and ob_Tv4.shape[0] > 0:
        im_d4_Tv = bd_Tv4 @ ob_Tv4.T % PRIME
    else:
        im_d4_Tv = np.zeros((len(paths_Tv_3), 1), dtype=np.int64)

    # im(d_4^Tv) in Omega_3(T\v) coords
    im_d4_Tv_omega = ob_Tv3 @ im_d4_Tv % PRIME
    rank_d4_Tv = _gauss_rank_np(im_d4_Tv_omega.copy(), PRIME)

    # Now: embed im(d_3(T\v)) and im(d_4(T\v)) into A_3(T) old-path space
    # and check if old-face-projection of new Omega_4 lands in their span

    # Build embedding matrix: A_3(T\v) → A_3(T)
    embed_3 = np.zeros((len(paths_T_3), len(paths_Tv_3)), dtype=np.int64)
    for j, ptv in enumerate(paths_Tv_3):
        emb = embed_path(ptv)
        if emb in idx_T3:
            embed_3[idx_T3[emb], j] = 1

    # im(d_3(T\v)) in A_3(T) = embed_3 @ bd_Tv3.T (actually we need the image in A_3)
    # d_3(T\v): A_3(T\v) → A_2(T\v)
    # Its image tells us nothing about A_3.
    # What we actually want: the KERNEL of the map A_3(T\v) → Omega_3(T\v)
    # (i.e., the NA relations) embedded into A_3(T).

    # Actually, let me think differently.
    # The old-face projection of d_4(new Omega_4) lands in A_3(T) restricted to old paths.
    # Old paths = paths not using v = embedded A_3(T\v) paths.
    # So the old-face projection is a vector in the span of embedded A_3(T\v).

    # Now HYP-398 says: when projected to ker(d_3^T) / im(d_4^T), this is zero.
    # Equivalently: the old-face projection of d_4(new) is in
    #   im(d_4^old) + (complement of ker(d_3^T) in A_3(T))

    # Let me compute: embed the full ker(d_3(T\v)) into A_3(T), then check
    # if old-face-projection of new d_4 lies in span of:
    #   (embedded im(d_4(T\v))) + (non-ker-d3-T part)

    # Actually the cleanest test: project old-face-proj to old Omega_3(T\v),
    # then to ker(d_3(T\v)), then check if result is in im(d_4(T\v)).

    # Step 1: New Omega_4 basis
    if ob_T4 is None:
        return None
    new_omega4 = []
    for i in range(ob_T4.shape[0]):
        uses_v = any(ob_T4[i, j] % PRIME != 0 and v in paths_T_4[j]
                     for j in range(len(paths_T_4)))
        if uses_v:
            new_omega4.append(i)
    if not new_omega4:
        return None

    # Step 2: d_4(new Omega_4) in A_3(T)
    bd_T4 = np.zeros((len(paths_T_3), len(paths_T_4)), dtype=np.int64)
    for j, sigma in enumerate(paths_T_4):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx_T3:
                bd_T4[idx_T3[ft], j] = (bd_T4[idx_T3[ft], j] + sign) % PRIME

    im_new = bd_T4 @ ob_T4[new_omega4].T % PRIME  # shape (|A_3(T)|, #new_omega4)

    # Step 3: Extract old-path components
    old_3idx = [i for i, q in enumerate(paths_T_3) if v not in q]
    im_new_old = im_new[old_3idx] % PRIME  # shape (#old_A3, #new_omega4)

    # Step 4: Map to T\v A_3 coordinates
    # old_3idx maps to embedded T\v paths
    # Build reverse map: old A_3(T) index → A_3(T\v) index
    old_to_Tv = {}
    for j, ptv in enumerate(paths_Tv_3):
        emb = embed_path(ptv)
        if emb in idx_T3:
            old_T_idx = idx_T3[emb]
            # Find position in old_3idx
            try:
                pos = old_3idx.index(old_T_idx)
                old_to_Tv[pos] = j
            except ValueError:
                pass

    # Convert im_new_old to T\v A_3 coordinates
    im_new_Tv = np.zeros((len(paths_Tv_3), len(new_omega4)), dtype=np.int64)
    for pos, Tv_idx in old_to_Tv.items():
        im_new_Tv[Tv_idx] = im_new_old[pos]

    # Step 5: Project to Omega_3(T\v) coordinates
    im_new_Tv_omega = ob_Tv3 @ im_new_Tv % PRIME  # shape (dim_Omega3_Tv, #new)

    # Step 6: Project to ker(d_3(T\v))
    im_new_in_ker = ker_d3_Tv @ im_new_Tv_omega % PRIME  # shape (ker_dim, #new)
    rank_in_ker = _gauss_rank_np(im_new_in_ker.copy(), PRIME)

    # Step 7: Check if this is already in im(d_4(T\v))
    # im_d4_Tv_omega is in Omega_3(T\v) coords
    # Project to ker(d_3) coords
    im_d4_in_ker = ker_d3_Tv @ im_d4_Tv_omega % PRIME
    rank_d4_in_ker = _gauss_rank_np(im_d4_in_ker.copy(), PRIME)

    # Combined
    combined = np.hstack([im_d4_in_ker, im_new_in_ker]) % PRIME
    rank_combined = _gauss_rank_np(combined.copy(), PRIME)

    new_adds_in_ker = rank_combined - rank_d4_in_ker

    result = {
        'old_face_in_Tv': old_face_in_Tv,
        'old_face_not_in_Tv': old_face_not_in_Tv,
        'no_old_face': no_old_face,
        'ker_dim_Tv': ker_d3_Tv.shape[0],
        'rank_d4_Tv': rank_d4_Tv,
        'rank_new_in_ker': rank_in_ker,
        'rank_d4_in_ker': rank_d4_in_ker,
        'rank_combined': rank_combined,
        'new_adds_in_ker': new_adds_in_ker,
        'v_positions': dict(sorted(v_position_counts.items())),
    }

    if verbose:
        print(f"    Old faces: {old_face_in_Tv} in A_3(T\\v), {old_face_not_in_Tv} NOT in A_3(T\\v)")
        print(f"    No old face: {no_old_face}")
        print(f"    v positions in new 5-paths: {result['v_positions']}")
        print(f"    ker_d3(T\\v) dim = {ker_d3_Tv.shape[0]}, rank(d_4^Tv) = {rank_d4_Tv}")
        print(f"    rank of new-old-proj in ker(d_3): {rank_in_ker}")
        print(f"    rank of d_4(T\\v) in ker(d_3): {rank_d4_in_ker}")
        print(f"    combined rank: {rank_combined}")
        print(f"    NEW KILLS in ker(d_3): {new_adds_in_ker}")

    return result


def main():
    print("=" * 70)
    print("OLD FACE ANALYSIS — DETAILED FACE STRUCTURE")
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
            r = detailed_face_analysis(A, n, vv, verbose=(found <= 3))
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

    # Are all old faces in A_3(T\v)?
    not_in_Tv = sum(r['old_face_not_in_Tv'] for r in results)
    print(f"\n  Old faces NOT in A_3(T\\v): {not_in_Tv}")
    if not_in_Tv == 0:
        print(f"  => ALL old faces of new 5-paths are in A_3(T\\v)!")
    else:
        print(f"  => Some old faces are NOT valid T\\v paths (exist in T but not T\\v)")

    # KEY: new additions to ker_d3
    kills_dist = Counter(r['new_adds_in_ker'] for r in results)
    print(f"\n  New kills in ker(d_3(T\\v)): {dict(sorted(kills_dist.items()))}")

    if all(r['new_adds_in_ker'] == 0 for r in results):
        print(f"\n  *** CONFIRMED: Old-face projection of new d_4, projected to")
        print(f"      ker(d_3(T\\v)), is ENTIRELY in im(d_4(T\\v))! ***")
        print(f"  This is HYP-398 from the T\\v perspective.")
    else:
        fail = sum(1 for r in results if r['new_adds_in_ker'] != 0)
        print(f"\n  HYP-398 fails via this route: {fail}/{len(results)}")

    # Rank comparisons
    print(f"\n  rank(new-old-proj in ker_d3): {Counter(r['rank_new_in_ker'] for r in results)}")
    print(f"  rank(d_4_Tv in ker_d3):       {Counter(r['rank_d4_in_ker'] for r in results)}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
