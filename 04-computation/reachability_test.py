"""
reachability_test.py — How much of ker(d_3^{T\v}) is "reachable" by old faces?

The old-face images of new Omega_4 vectors land in A_3(T\v) (mostly).
Their intersection with ker(d_3^{T\v}) has some dimension.
If this dimension < dim(ker(d_3^{T\v})) - rk(d_4^{T\v}),
then they can't reach H_3(T\v).

More precisely: project old-face images to ker(d_3^{T\v})/im(d_4^{T\v}).
If this projection is zero, they can't kill H_3.

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


def reachability_analysis(A, n, v):
    """Check if old-face images can reach H_3(T\v)."""
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

    paths_4_T = ap_T.get(4, [])
    paths_3_Tv = ap_Tv.get(3, [])
    paths_4_Tv = ap_Tv.get(4, [])
    paths_2_Tv = ap_Tv.get(2, [])

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

    idx3Tv = {p: i for i, p in enumerate(paths_3_Tv)}
    idx2Tv = {p: i for i, p in enumerate(paths_2_Tv)}

    remaining_inv = {remaining[i]: i for i in range(n1)}

    # Build v-deletion map: through-v 4-path -> T\v 3-path
    v_del_map = np.zeros((len(paths_4_T), len(paths_3_Tv)), dtype=np.int64)
    for j, path in enumerate(paths_4_T):
        if v not in path:
            continue
        for sign, face in boundary_faces(path):
            if v not in face:
                tv_face = tuple(remaining_inv.get(x, -1) for x in face)
                if -1 not in tv_face and tv_face in idx3Tv:
                    v_del_map[j, idx3Tv[tv_face]] = (v_del_map[j, idx3Tv[tv_face]] + sign) % PRIME

    # New Omega_4 vectors
    old_4_path_set = set(i for i, p in enumerate(paths_4_T) if v not in p)
    new_omega4_rows = [i for i in range(ob4_T.shape[0])
                       if any(ob4_T[i, j] % PRIME != 0 for j in range(len(paths_4_T)) if j not in old_4_path_set)]

    if new_omega4_rows:
        new_ob4 = ob4_T[new_omega4_rows]
        # Old-face images in A_3(T\v)
        old_face_images = new_ob4 @ v_del_map % PRIME
    else:
        old_face_images = np.zeros((0, len(paths_3_Tv)), dtype=np.int64)

    # ker(d_3^{T\v}) and im(d_4^{T\v}) in A_3(T\v)
    bd3_Tv = np.zeros((len(paths_2_Tv), len(paths_3_Tv)), dtype=np.int64)
    for j, path in enumerate(paths_3_Tv):
        for sign, face in boundary_faces(path):
            if face in idx2Tv:
                bd3_Tv[idx2Tv[face], j] = (bd3_Tv[idx2Tv[face], j] + sign) % PRIME

    d3o_Tv = bd3_Tv @ ob3_Tv.T % PRIME
    if d3o_Tv.size > 0:
        _, kv = _gauss_nullbasis_modp(d3o_Tv.astype(np.int32), d3o_Tv.shape[0], d3o_Tv.shape[1], PRIME)
        ker_d3Tv = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob3_Tv.shape[0]), dtype=np.int64)
    else:
        ker_d3Tv = np.eye(ob3_Tv.shape[0], dtype=np.int64)

    ker_d3Tv_A = ker_d3Tv @ ob3_Tv % PRIME if ker_d3Tv.shape[0] > 0 else np.zeros((0, len(paths_3_Tv)), dtype=np.int64)

    bd4_Tv_mat = np.zeros((len(paths_3_Tv), len(paths_4_Tv)), dtype=np.int64)
    for j, path in enumerate(paths_4_Tv):
        for sign, face in boundary_faces(path):
            if face in idx3Tv:
                bd4_Tv_mat[idx3Tv[face], j] = (bd4_Tv_mat[idx3Tv[face], j] + sign) % PRIME

    if ob4_Tv.shape[0] > 0:
        im_d4Tv = (bd4_Tv_mat @ ob4_Tv.T % PRIME).T
    else:
        im_d4Tv = np.zeros((0, len(paths_3_Tv)), dtype=np.int64)

    rk_d4Tv = int(_gauss_rank_np(im_d4Tv.copy(), PRIME)) if im_d4Tv.shape[0] > 0 else 0
    dim_ker = ker_d3Tv_A.shape[0]

    # Intersect old_face_images with ker(d_3^{T\v})
    # A vector is in ker iff d_3 applied to it is zero
    # For each old_face_image row, check if it's in ker
    if old_face_images.shape[0] > 0:
        d3_app = bd3_Tv @ old_face_images.T % PRIME
        in_ker_mask = [np.all(d3_app[:, col] % PRIME == 0) for col in range(d3_app.shape[1])]
        old_face_in_ker = old_face_images[in_ker_mask]
    else:
        old_face_in_ker = np.zeros((0, len(paths_3_Tv)), dtype=np.int64)

    rk_old_face_in_ker = int(_gauss_rank_np(old_face_in_ker.copy(), PRIME)) if old_face_in_ker.shape[0] > 0 else 0

    # How much of ker/im does it reach?
    # Combine with im(d_4^{T\v}) and check rank increase
    if im_d4Tv.shape[0] > 0 and old_face_in_ker.shape[0] > 0:
        combined = np.vstack([im_d4Tv, old_face_in_ker]) % PRIME
        rk_combined = int(_gauss_rank_np(combined.copy(), PRIME))
        reach_beyond_im = rk_combined - rk_d4Tv
    elif old_face_in_ker.shape[0] > 0:
        rk_combined = rk_old_face_in_ker
        reach_beyond_im = rk_old_face_in_ker
    else:
        rk_combined = rk_d4Tv
        reach_beyond_im = 0

    # The H_3 quotient has dim 1 (since b3=1).
    # reach_beyond_im = 1 means old faces CAN reach H_3 direction
    # reach_beyond_im = 0 means old faces stay in im(d_4) — CAN'T kill H_3

    out_deg = int(sum(A[v]))

    return {
        'out_deg': out_deg,
        'dim_ker_Tv': dim_ker,
        'rk_d4Tv': rk_d4Tv,
        'b3_Tv': dim_ker - rk_d4Tv,
        'n_new_omega4': len(new_omega4_rows),
        'rk_old_face_total': int(_gauss_rank_np(old_face_images.copy(), PRIME)) if old_face_images.shape[0] > 0 else 0,
        'n_old_face_in_ker': old_face_in_ker.shape[0],
        'rk_old_face_in_ker': rk_old_face_in_ker,
        'reach_beyond_im': reach_beyond_im,
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"REACHABILITY ANALYSIS AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        target = 200 if n == 7 else 120
        t0 = time.time()

        for trial in range(30000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) != 1:
                continue
            for v_cand in range(n):
                r = reachability_analysis(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)
            if len(results) % 50 == 0 and len(results) > 0:
                elapsed = time.time() - t0
                print(f"  Progress: {len(results)}, {elapsed:.1f}s")

        elapsed = time.time() - t0
        print(f"  {len(results)} BAD vertices, {elapsed:.1f}s")

        # KEY: reach_beyond_im distribution
        reach = [r['reach_beyond_im'] for r in results]
        from collections import Counter
        reach_dist = Counter(reach)
        print(f"\n  reach_beyond_im distribution: {dict(sorted(reach_dist.items()))}")
        print(f"  (0 = old faces stay in im(d4_Tv), can't kill H3)")
        print(f"  (1 = old faces reach H3 direction, CAN kill H3)")

        # Other stats
        print(f"\n  rk_old_face_in_ker: mean={np.mean([r['rk_old_face_in_ker'] for r in results]):.1f}")
        print(f"  rk_d4Tv:            mean={np.mean([r['rk_d4Tv'] for r in results]):.1f}")
        print(f"  dim_ker_Tv:         mean={np.mean([r['dim_ker_Tv'] for r in results]):.1f}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
