"""
single_face_structure.py — Verify the single-old-face structure

KEY INSIGHT: For a through-v 4-path (a,b,c,d,e) with v at position k,
the boundary d_4 has faces:
  (-1)^0 (b,c,d,e), (-1)^1 (a,c,d,e), ..., (-1)^4 (a,b,c,d)

The ONLY old face (not containing v) is the one where v is deleted.
All other faces still contain v, so they're "new" 3-paths.

In Omega basis, a single Omega_4 vector is a linear combination of
allowed 4-paths. Its boundary is a linear combination of faces.
The old-projected part is just the sum of the v-deletion faces.

This means: im(new_bdys_old) is spanned by vectors of a VERY specific form.
Each new Omega_4 vector contributes a vector that's the old-projection
of deleting v from each of its through-v 4-path constituents.

Does this structural constraint force the image to miss H_3(T\v)?

Author: opus-2026-03-09-S57
"""
import sys
import time
import numpy as np
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


def verify_single_face(A, n, v):
    """Verify that each through-v 4-path has exactly one old face in boundary."""
    ap_T = enumerate_all_allowed(A, n, min(n - 1, 6))
    paths_4_T = ap_T.get(4, [])
    paths_3_T = ap_T.get(3, [])

    n_through_v_4 = 0
    n_old_faces_per_path = []

    for path in paths_4_T:
        if v not in path:
            continue
        n_through_v_4 += 1

        old_face_count = 0
        for sign, face in boundary_faces(path):
            if v not in face:
                old_face_count += 1

        n_old_faces_per_path.append(old_face_count)

    return n_through_v_4, n_old_faces_per_path


def detailed_face_analysis(A, n, v):
    """Analyze what the old faces look like for through-v 4-paths."""
    max_p = min(n - 1, 6)
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]

    cc_T = full_chain_complex_modp(A, n, max_p)
    cc_Tv = full_chain_complex_modp(A_sub, n1, min(max_p, n1 - 1))

    if cc_T['bettis'].get(3, 0) != 1 or cc_Tv['bettis'].get(3, 0) != 1:
        return None

    ap_T = enumerate_all_allowed(A, n, max_p)
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    paths_3_T = ap_T.get(3, [])
    paths_4_T = ap_T.get(4, [])
    paths_3_Tv = ap_Tv.get(3, [])
    paths_2_Tv = ap_Tv.get(2, [])

    idx3T = {p: i for i, p in enumerate(paths_3_T)}
    idx3Tv = {p: i for i, p in enumerate(paths_3_Tv)}

    # For each through-v 4-path, find the old face and check if it's in A_3(T\v)
    old_3_idx = [i for i, p in enumerate(paths_3_T) if v not in p]
    old_3_set = set(old_3_idx)

    # Map from global old-3-path to T\v path index
    remaining_inv = {remaining[i]: i for i in range(n1)}
    def to_Tv_path(face):
        """Convert an old face (global indices) to T\v indices."""
        return tuple(remaining_inv[x] for x in face)

    n_old_face_in_Tv = 0
    n_old_face_total = 0
    n_through_v_4 = 0

    for path in paths_4_T:
        if v not in path:
            continue
        n_through_v_4 += 1

        for sign, face in boundary_faces(path):
            if v not in face:
                n_old_face_total += 1
                # Is this face an allowed 3-path in T\v?
                tv_face = to_Tv_path(face)
                if tv_face in idx3Tv:
                    n_old_face_in_Tv += 1

    # Now the key question: do the old faces of through-v 4-paths
    # span a subspace that intersects H_3(T\v) nontrivially?
    #
    # Build the "old face map": for each through-v 4-path,
    # output its old face as a vector in A_3(T\v)
    def get_omega(ap, deg):
        paths = ap.get(deg, [])
        if not paths:
            return np.zeros((0, 0), dtype=np.int64)
        P, nr, nc = _build_constraint_matrix(ap, deg, PRIME)
        if P is not None:
            _, nb = _gauss_nullbasis_modp(P, nr, nc, PRIME)
            return np.array(nb, dtype=np.int64) if nb else np.zeros((0, nc), dtype=np.int64)
        return np.eye(len(paths), dtype=np.int64)

    ob4_T = get_omega(ap_T, 4)
    ob3_Tv = get_omega(ap_Tv, 3)
    ob4_Tv = get_omega(ap_Tv, 4)

    # Build the "v-deletion map": 4-path through v -> 3-path in T\v
    # For each through-v 4-path index j, find its old face in T\v coords
    v_del_map = np.zeros((len(paths_4_T), len(paths_3_Tv)), dtype=np.int64)
    for j, path in enumerate(paths_4_T):
        if v not in path:
            continue
        for sign, face in boundary_faces(path):
            if v not in face:
                tv_face = to_Tv_path(face)
                if tv_face in idx3Tv:
                    v_del_map[j, idx3Tv[tv_face]] = (v_del_map[j, idx3Tv[tv_face]] + sign) % PRIME

    # Apply to Omega_4(T) basis vectors to get old-face images in A_3(T\v)
    if ob4_T.shape[0] > 0:
        omega4_old_face = (ob4_T @ v_del_map % PRIME)  # rows = old-face images per Omega_4 vec
    else:
        omega4_old_face = np.zeros((0, len(paths_3_Tv)), dtype=np.int64)

    # Restrict to through-v Omega_4 vectors
    old_4_path_set = set(i for i, p in enumerate(paths_4_T) if v not in p)
    new_omega4_rows = []
    for row_i in range(ob4_T.shape[0]):
        has_new = any(ob4_T[row_i, j] % PRIME != 0 for j in range(len(paths_4_T)) if j not in old_4_path_set)
        if has_new:
            new_omega4_rows.append(row_i)

    if new_omega4_rows:
        new_old_face_images = omega4_old_face[new_omega4_rows]
    else:
        new_old_face_images = np.zeros((0, len(paths_3_Tv)), dtype=np.int64)

    # Key: is the image of new_old_face_images in Omega_3(T\v) or in A_3(T\v)?
    # It's in A_3(T\v) coords. Project to Omega_3(T\v) using the constraint matrix.
    # Actually, we should check if these vectors are in Omega_3(T\v) or not.
    # But more importantly: are they in im(d_4^{T\v})?

    # im(d_4^{T\v}) in A_3(T\v)
    bd4_Tv = np.zeros((len(paths_3_Tv), ap_Tv.get(4, []).__len__()), dtype=np.int64)
    for j, path in enumerate(ap_Tv.get(4, [])):
        for sign, face in boundary_faces(path):
            if face in idx3Tv:
                bd4_Tv[idx3Tv[face], j] = (bd4_Tv[idx3Tv[face], j] + sign) % PRIME

    if ob4_Tv.shape[0] > 0:
        im_d4Tv = (bd4_Tv @ ob4_Tv.T % PRIME).T
    else:
        im_d4Tv = np.zeros((0, len(paths_3_Tv)), dtype=np.int64)

    rk_d4Tv = int(_gauss_rank_np(im_d4Tv.copy(), PRIME)) if im_d4Tv.shape[0] > 0 else 0

    # Are new_old_face_images in im(d_4^{T\v})?
    if im_d4Tv.shape[0] > 0 and new_old_face_images.shape[0] > 0:
        combined = np.vstack([im_d4Tv, new_old_face_images]) % PRIME
        rk_comb = int(_gauss_rank_np(combined.copy(), PRIME))
        new_in_im = (rk_comb == rk_d4Tv)
        extra_rank = rk_comb - rk_d4Tv
    else:
        new_in_im = True
        extra_rank = 0

    # Also check: are new_old_face_images even in ker(d_3^{T\v})?
    # d_3^{T\v} applied to each new_old_face_image
    bd3_Tv = np.zeros((len(ap_Tv.get(2, [])), len(paths_3_Tv)), dtype=np.int64)
    idx2Tv = {p: i for i, p in enumerate(ap_Tv.get(2, []))}
    for j, path in enumerate(paths_3_Tv):
        for sign, face in boundary_faces(path):
            if face in idx2Tv:
                bd3_Tv[idx2Tv[face], j] = (bd3_Tv[idx2Tv[face], j] + sign) % PRIME

    if new_old_face_images.shape[0] > 0:
        d3_applied = bd3_Tv @ new_old_face_images.T % PRIME
        n_in_ker = 0
        for col in range(d3_applied.shape[1]):
            if np.all(d3_applied[:, col] % PRIME == 0):
                n_in_ker += 1
    else:
        n_in_ker = 0

    out_deg = int(sum(A[v]))

    return {
        'n_through_v_4': n_through_v_4,
        'n_old_face_in_Tv': n_old_face_in_Tv,
        'n_old_face_total': n_old_face_total,
        'n_new_omega4': len(new_omega4_rows),
        'rk_new_old_face': int(_gauss_rank_np(new_old_face_images.copy(), PRIME)) if new_old_face_images.shape[0] > 0 else 0,
        'rk_d4Tv': rk_d4Tv,
        'new_in_im_d4Tv': new_in_im,
        'extra_rank_over_im': extra_rank,
        'n_in_ker_d3Tv': n_in_ker,
        'out_deg': out_deg,
    }


def main():
    # First: verify single-face property
    print("=" * 70)
    print("VERIFICATION: Through-v 4-paths have exactly 1 old face")
    print("=" * 70)

    rng = np.random.RandomState(42)
    all_single = True
    for _ in range(100):
        n = 7
        A = random_tournament(n, rng)
        for v in range(n):
            n4, counts = verify_single_face(A, n, v)
            if any(c != 1 for c in counts):
                print(f"  VIOLATION at n={n}, v={v}: counts={counts}")
                all_single = False
    for _ in range(50):
        n = 8
        A = random_tournament(n, rng)
        for v in range(n):
            n4, counts = verify_single_face(A, n, v)
            if any(c != 1 for c in counts):
                print(f"  VIOLATION at n={n}, v={v}: counts={counts}")
                all_single = False

    if all_single:
        print("  CONFIRMED: Every through-v 4-path has exactly 1 old face in boundary")

    # Detailed analysis
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"DETAILED OLD-FACE ANALYSIS AT n={n}")
        print(f"{'='*70}")

        rng2 = np.random.RandomState(42)
        results = []
        target = 150 if n == 7 else 100
        t0 = time.time()

        for trial in range(20000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng2)
            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) != 1:
                continue
            for v_cand in range(n):
                r = detailed_face_analysis(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)

            if len(results) % 50 == 0 and len(results) > 0:
                elapsed = time.time() - t0
                print(f"  Progress: {len(results)}, {elapsed:.1f}s")

        elapsed = time.time() - t0
        print(f"  {len(results)} BAD vertices, {elapsed:.1f}s")

        # Are old faces always allowed in T\v?
        always_in_Tv = all(r['n_old_face_in_Tv'] == r['n_old_face_total'] for r in results)
        print(f"\n  Old faces always in A_3(Tv): {always_in_Tv}")
        if not always_in_Tv:
            miss_pcts = [(1 - r['n_old_face_in_Tv']/max(r['n_old_face_total'],1)) for r in results]
            print(f"    miss rate: mean={np.mean(miss_pcts):.3f}")

        # New old-face images: in im(d_4^{T\v})?
        n_in_im = sum(1 for r in results if r['new_in_im_d4Tv'])
        print(f"\n  new_old_face ⊂ im(d_4^Tv): {n_in_im}/{len(results)}")
        extras = [r['extra_rank_over_im'] for r in results]
        print(f"  extra rank over im(d_4^Tv): min={min(extras)}, max={max(extras)}, mean={np.mean(extras):.1f}")

        # Are new old-face images in ker(d_3^{T\v})?
        in_ker = [r['n_in_ker_d3Tv'] for r in results]
        n_omega = [r['n_new_omega4'] for r in results]
        print(f"\n  New Omega_4 vectors: mean={np.mean(n_omega):.1f}")
        print(f"  Of those, old-face image in ker(d_3^Tv): mean={np.mean(in_ker):.1f}")
        in_ker_pct = [r['n_in_ker_d3Tv']/max(r['n_new_omega4'],1) for r in results]
        print(f"  Fraction in ker: mean={np.mean(in_ker_pct):.3f}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
