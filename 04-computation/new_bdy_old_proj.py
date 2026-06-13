"""
new_bdy_old_proj.py — Do through-v 4-chain boundaries project to im(d_4^{T\v})?

Key question: im(d_4^T)_old has rank ~19.3 at n=7 but im(d_4^{T\v}) has rank ~4.2.
The extra 15 directions come from through-v 4-chains whose boundaries have
nonzero old-path components.

If these extra old-projections all land in im(d_4^{T\v}), then im(d_4^T)_old
would equal im(d_4^{T\v}), contradicting our data. So they DON'T all land there.

But the question is: do they avoid the H_3(T\v) direction?

Formally: Let W = im(d_4^{T\v}) (in old coords), dim ~4.2.
Let V = im(d_4^T)_old (in old coords), dim ~19.3.
Both sit inside ker(d_3^T)_old, dim ~20.5.

The extra directions V \ W form a ~15-dim space.
The H_3(T\v) gen is a single direction in ker(d_3^T)_old \ V (at n=7, always).

What if: the extra directions V \ W are always ORTHOGONAL (in some sense)
to the H_3 direction? Or: there's an algebraic reason that through-v
boundaries can never produce the H_3 component.

Test: for each through-v 4-chain's boundary, project to old coords,
then check if it has a component along the H_3(T\v) direction.

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


def new_bdy_projection_analysis(A, n, v):
    """Check if through-v 4-boundaries projected to old coords avoid H_3(T\v)."""
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

    ob3_T = get_omega(ap_T, 3)
    ob4_T = get_omega(ap_T, 4)
    ob3_Tv = get_omega(ap_Tv, 3)
    ob4_Tv = get_omega(ap_Tv, 4)

    # Build d_4^T: A_4(T) -> A_3(T)
    idx3T = {p: i for i, p in enumerate(paths_3_T)}
    bd4_T = np.zeros((len(paths_3_T), len(paths_4_T)), dtype=np.int64)
    for j, path in enumerate(paths_4_T):
        for sign, face in boundary_faces(path):
            if face in idx3T:
                bd4_T[idx3T[face], j] = (bd4_T[idx3T[face], j] + sign) % PRIME

    # im(d_4^T) in A_3(T) = boundary images of Omega_4 basis vectors
    if ob4_T.shape[0] > 0:
        im_d4T_A = (bd4_T @ ob4_T.T % PRIME).T
    else:
        im_d4T_A = np.zeros((0, len(paths_3_T)), dtype=np.int64)

    # Classify Omega_4 basis vectors as old/new
    old_4_path_set = set(i for i, p in enumerate(paths_4_T) if v not in p)
    new_omega4_rows = []
    old_omega4_rows = []
    for row_i in range(ob4_T.shape[0]):
        has_new = any(ob4_T[row_i, j] % PRIME != 0 for j in range(len(paths_4_T)) if j not in old_4_path_set)
        if has_new:
            new_omega4_rows.append(row_i)
        else:
            old_omega4_rows.append(row_i)

    # Boundaries of new Omega_4 vectors
    if new_omega4_rows:
        new_bdys = im_d4T_A[new_omega4_rows]  # rows = boundary images
    else:
        new_bdys = np.zeros((0, len(paths_3_T)), dtype=np.int64)

    # Old 3-path indices
    old_3_idx = [i for i, p in enumerate(paths_3_T) if v not in p]

    # Project new boundaries to old coords
    if new_bdys.shape[0] > 0 and old_3_idx:
        new_bdys_old = new_bdys[:, old_3_idx] % PRIME
    else:
        new_bdys_old = np.zeros((0, len(old_3_idx)), dtype=np.int64)

    # Compute the H_3(T\v) generator in old A_3(T) coords
    # ker(d_3^{T\v}) in A_3(T\v)
    idx2Tv = {p: i for i, p in enumerate(paths_2_Tv)}
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

    # im(d_4^{T\v}) in A_3(T\v)
    idx3Tv = {p: i for i, p in enumerate(paths_3_Tv)}
    bd4_Tv = np.zeros((len(paths_3_Tv), len(paths_4_Tv)), dtype=np.int64)
    for j, path in enumerate(paths_4_Tv):
        for sign, face in boundary_faces(path):
            if face in idx3Tv:
                bd4_Tv[idx3Tv[face], j] = (bd4_Tv[idx3Tv[face], j] + sign) % PRIME

    if ob4_Tv.shape[0] > 0:
        im_d4Tv = (bd4_Tv @ ob4_Tv.T % PRIME).T
    else:
        im_d4Tv = np.zeros((0, len(paths_3_Tv)), dtype=np.int64)

    rk_d4Tv = int(_gauss_rank_np(im_d4Tv.copy(), PRIME)) if im_d4Tv.shape[0] > 0 else 0

    # Find H_3(T\v) generator: element of ker NOT in im
    h3Tv_gen = None
    for i in range(ker_d3Tv_A.shape[0]):
        if im_d4Tv.shape[0] > 0:
            test = np.vstack([im_d4Tv, ker_d3Tv_A[i:i+1]]) % PRIME
            rk = int(_gauss_rank_np(test.copy(), PRIME))
            if rk > rk_d4Tv:
                h3Tv_gen = ker_d3Tv_A[i]
                break
        else:
            h3Tv_gen = ker_d3Tv_A[i]
            break

    if h3Tv_gen is None:
        return None

    # Embed H_3(T\v) gen into A_3(T) old coords
    embed_old = np.zeros(len(old_3_idx), dtype=np.int64)
    old_3_set = set(old_3_idx)
    old_3_pos = {idx: pos for pos, idx in enumerate(old_3_idx)}
    for jv, path_v in enumerate(paths_3_Tv):
        mapped = tuple(remaining[x] for x in path_v)
        if mapped in idx3T:
            global_idx = idx3T[mapped]
            if global_idx in old_3_pos:
                embed_old[old_3_pos[global_idx]] = (embed_old[old_3_pos[global_idx]] + h3Tv_gen[jv]) % PRIME

    # Now: do the new_bdys_old span a space that contains embed_old?
    # Also: does span(im_d4Tv_old, new_bdys_old) contain embed_old?

    # Embed im(d_4^{T\v}) into old A_3(T) coords
    embed_map = np.zeros((len(paths_3_Tv), len(old_3_idx)), dtype=np.int64)
    for jv, path_v in enumerate(paths_3_Tv):
        mapped = tuple(remaining[x] for x in path_v)
        if mapped in idx3T:
            gi = idx3T[mapped]
            if gi in old_3_pos:
                embed_map[jv, old_3_pos[gi]] = 1

    emb_im_d4Tv_old = im_d4Tv @ embed_map % PRIME if im_d4Tv.shape[0] > 0 else np.zeros((0, len(old_3_idx)), dtype=np.int64)
    rk_emb_im = int(_gauss_rank_np(emb_im_d4Tv_old.copy(), PRIME)) if emb_im_d4Tv_old.shape[0] > 0 else 0

    # Test: is embed_old in span(new_bdys_old)?
    if new_bdys_old.shape[0] > 0:
        rk_new = int(_gauss_rank_np(new_bdys_old.copy(), PRIME))
        test_new = np.vstack([new_bdys_old, embed_old.reshape(1, -1)]) % PRIME
        rk_new_plus = int(_gauss_rank_np(test_new.copy(), PRIME))
        gen_in_new_bdy = (rk_new_plus == rk_new)
    else:
        rk_new = 0
        gen_in_new_bdy = False

    # Test: is embed_old in span(emb_im_d4Tv_old + new_bdys_old)?
    if emb_im_d4Tv_old.shape[0] > 0 and new_bdys_old.shape[0] > 0:
        combined_old = np.vstack([emb_im_d4Tv_old, new_bdys_old]) % PRIME
    elif new_bdys_old.shape[0] > 0:
        combined_old = new_bdys_old
    elif emb_im_d4Tv_old.shape[0] > 0:
        combined_old = emb_im_d4Tv_old
    else:
        combined_old = np.zeros((0, len(old_3_idx)), dtype=np.int64)

    if combined_old.shape[0] > 0:
        rk_comb = int(_gauss_rank_np(combined_old.copy(), PRIME))
        test_comb = np.vstack([combined_old, embed_old.reshape(1, -1)]) % PRIME
        rk_comb_plus = int(_gauss_rank_np(test_comb.copy(), PRIME))
        gen_in_combined = (rk_comb_plus == rk_comb)
    else:
        rk_comb = 0
        gen_in_combined = False

    # rank(i_*) for reference
    embed_full = np.zeros((len(paths_3_Tv), len(paths_3_T)), dtype=np.int64)
    for jv, path_v in enumerate(paths_3_Tv):
        mapped = tuple(remaining[x] for x in path_v)
        if mapped in idx3T:
            embed_full[jv, idx3T[mapped]] = 1

    embedded_ker = ker_d3Tv_A @ embed_full % PRIME
    if im_d4T_A.shape[0] > 0 and embedded_ker.shape[0] > 0:
        full_combined = np.vstack([im_d4T_A, embedded_ker]) % PRIME
    elif embedded_ker.shape[0] > 0:
        full_combined = embedded_ker
    else:
        full_combined = im_d4T_A
    rk_d4T = int(_gauss_rank_np(im_d4T_A.copy(), PRIME)) if im_d4T_A.shape[0] > 0 else 0
    rk_full = int(_gauss_rank_np(full_combined.copy(), PRIME)) if full_combined.shape[0] > 0 else 0
    rank_istar = rk_full - rk_d4T

    out_deg = int(sum(A[v]))

    return {
        'rank_istar': rank_istar,
        'out_deg': out_deg,
        'rk_new_bdys_old': rk_new if new_bdys_old.shape[0] > 0 else 0,
        'rk_emb_im_d4Tv': rk_emb_im,
        'rk_combined_old': rk_comb,
        'gen_in_new_bdy': gen_in_new_bdy,
        'gen_in_combined': gen_in_combined,
        'n_new_omega4': len(new_omega4_rows),
        'n_old_omega4': len(old_omega4_rows),
        'n_old_3': len(old_3_idx),
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"NEW BOUNDARY OLD-PROJECTION ANALYSIS AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        checked = 0
        target = 200 if n == 7 else 120
        t0 = time.time()

        for trial in range(30000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            checked += 1
            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) != 1:
                continue
            for v_cand in range(n):
                r = new_bdy_projection_analysis(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)
            if len(results) % 50 == 0 and len(results) > 0:
                elapsed = time.time() - t0
                print(f"  Progress: {len(results)}, {elapsed:.1f}s")

        elapsed = time.time() - t0
        fails = [r for r in results if r['rank_istar'] == 0]
        succs = [r for r in results if r['rank_istar'] == 1]
        print(f"  {len(results)} BAD vertices ({len(fails)} fail), {elapsed:.1f}s")

        # Is H_3(T\v) gen in span of new boundaries (old-projected)?
        in_new = sum(1 for r in results if r['gen_in_new_bdy'])
        in_comb = sum(1 for r in results if r['gen_in_combined'])
        print(f"\n  H3(Tv) gen IN span(new_bdys_old): {in_new}/{len(results)}")
        print(f"  H3(Tv) gen IN span(im_d4Tv_old + new_bdys_old): {in_comb}/{len(results)}")

        # Rank stats
        print(f"\n  rk(new_bdys_old): mean={np.mean([r['rk_new_bdys_old'] for r in results]):.1f}")
        print(f"  rk(emb_im_d4Tv):  mean={np.mean([r['rk_emb_im_d4Tv'] for r in results]):.1f}")
        print(f"  rk(combined_old): mean={np.mean([r['rk_combined_old'] for r in results]):.1f}")
        print(f"  n_old_3:          mean={np.mean([r['n_old_3'] for r in results]):.1f}")

        if fails:
            print(f"\n  FAILURES:")
            for r in fails[:5]:
                print(f"    out={r['out_deg']}: gen_in_new={r['gen_in_new_bdy']}, "
                      f"gen_in_comb={r['gen_in_combined']}, "
                      f"rk_new={r['rk_new_bdys_old']}, rk_emb_im={r['rk_emb_im_d4Tv']}")

        # Sample successes
        print(f"\n  SUCCESSES (sample):")
        for r in succs[:5]:
            print(f"    out={r['out_deg']}: gen_in_new={r['gen_in_new_bdy']}, "
                  f"gen_in_comb={r['gen_in_combined']}, "
                  f"rk_new={r['rk_new_bdys_old']}, rk_emb_im={r['rk_emb_im_d4Tv']}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
