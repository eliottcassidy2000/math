"""
orphan_face_mechanism.py — Test the orphan face hypothesis

KEY INSIGHT: The v-deletion face of a through-v 4-path is always an old
3-path in A_3(T). But it might NOT be in A_3(T\v) — it's an "orphan" face.

In A_3(T\v) coords: these orphans are invisible.
In A_3(T) old coords: they ARE visible and contribute to im(d_4^T)_old.

Hypothesis: rank(i_*)=0 failures happen when orphan faces span a direction
that captures the H_3(T\v) generator (embedded in A_3(T) old coords).

Test: compute im(d_4^T)_old using ONLY orphan face contributions.
If this already captures the gen, orphans are the mechanism.
If not, the gen gets captured only when orphan + non-orphan combine.

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


def orphan_analysis(A, n, v):
    """Analyze the role of orphan old faces in i_* failure."""
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

    idx3T = {p: i for i, p in enumerate(paths_3_T)}

    # Classify old 3-paths as "embedded" (in T\v) or "orphan" (in T only)
    remaining_inv = {remaining[i]: i for i in range(n1)}
    idx3Tv_set = set(ap_Tv.get(3, []))  # T\v 3-paths as tuples in T\v indexing

    embedded_old_3 = []  # indices in paths_3_T
    orphan_old_3 = []    # indices in paths_3_T

    for i, p in enumerate(paths_3_T):
        if v in p:
            continue  # new path, skip
        # Map to T\v indices
        tv_p = tuple(remaining_inv.get(x, -1) for x in p)
        if -1 not in tv_p and tv_p in idx3Tv_set:
            embedded_old_3.append(i)
        else:
            orphan_old_3.append(i)

    n_embedded = len(embedded_old_3)
    n_orphan = len(orphan_old_3)

    # Build d_4^T: full boundary matrix
    bd4_T = np.zeros((len(paths_3_T), len(paths_4_T)), dtype=np.int64)
    for j, path in enumerate(paths_4_T):
        for sign, face in boundary_faces(path):
            if face in idx3T:
                bd4_T[idx3T[face], j] = (bd4_T[idx3T[face], j] + sign) % PRIME

    # im(d_4^T) in A_3(T)
    if ob4_T.shape[0] > 0:
        im_d4T = (bd4_T @ ob4_T.T % PRIME).T
    else:
        im_d4T = np.zeros((0, len(paths_3_T)), dtype=np.int64)

    rk_d4T = int(_gauss_rank_np(im_d4T.copy(), PRIME)) if im_d4T.shape[0] > 0 else 0

    # All old indices
    old_3_idx = embedded_old_3 + orphan_old_3

    # Project im(d_4^T) to old coords
    if im_d4T.shape[0] > 0 and old_3_idx:
        im_old = im_d4T[:, old_3_idx] % PRIME
        rk_im_old = int(_gauss_rank_np(im_old.copy(), PRIME))
    else:
        im_old = np.zeros((0, len(old_3_idx)), dtype=np.int64)
        rk_im_old = 0

    # Project to embedded-only coords
    if im_d4T.shape[0] > 0 and embedded_old_3:
        im_embedded = im_d4T[:, embedded_old_3] % PRIME
        rk_im_embedded = int(_gauss_rank_np(im_embedded.copy(), PRIME))
    else:
        rk_im_embedded = 0

    # Project to orphan-only coords
    if im_d4T.shape[0] > 0 and orphan_old_3:
        im_orphan = im_d4T[:, orphan_old_3] % PRIME
        rk_im_orphan = int(_gauss_rank_np(im_orphan.copy(), PRIME))
    else:
        rk_im_orphan = 0

    # Embedded H_3(T\v) gen in old A_3(T) coords (embedded part only, orphan part = 0)
    # First get ker(d_3^{T\v}) and find H_3 gen
    ob3_Tv = get_omega(ap_Tv, 3)
    ob4_Tv = get_omega(ap_Tv, 4)
    paths_2_Tv = ap_Tv.get(2, [])
    paths_4_Tv = ap_Tv.get(4, [])
    idx2Tv = {p: i for i, p in enumerate(paths_2_Tv)}
    idx3Tv = {p: i for i, p in enumerate(paths_3_Tv)}

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

    # Embed ker and gen into A_3(T) (only on embedded old coords)
    embed = np.zeros((len(paths_3_Tv), len(paths_3_T)), dtype=np.int64)
    for jv, path_v in enumerate(paths_3_Tv):
        mapped = tuple(remaining[x] for x in path_v)
        if mapped in idx3T:
            embed[jv, idx3T[mapped]] = 1

    embedded_ker = ker_d3Tv_A @ embed % PRIME if ker_d3Tv_A.shape[0] > 0 else np.zeros((0, len(paths_3_T)), dtype=np.int64)

    # rank(i_*)
    if im_d4T.shape[0] > 0 and embedded_ker.shape[0] > 0:
        combined = np.vstack([im_d4T, embedded_ker]) % PRIME
    elif embedded_ker.shape[0] > 0:
        combined = embedded_ker
    else:
        combined = im_d4T
    rk_combined = int(_gauss_rank_np(combined.copy(), PRIME)) if combined.shape[0] > 0 else 0
    rank_istar = rk_combined - rk_d4T

    # Embedded gen projected to old coords only
    if embedded_ker.shape[0] > 0:
        emb_old = embedded_ker[:, old_3_idx] % PRIME
        # Check which old coords have nonzero: embedded vs orphan
        emb_embedded_part = embedded_ker[:, embedded_old_3] % PRIME if embedded_old_3 else np.zeros((embedded_ker.shape[0], 0), dtype=np.int64)
        emb_orphan_part = embedded_ker[:, orphan_old_3] % PRIME if orphan_old_3 else np.zeros((embedded_ker.shape[0], 0), dtype=np.int64)

        # The embedded gen should have zero on orphan coords (since embedding maps T\v -> T via relabeling)
        orphan_nonzero = 0
        if emb_orphan_part.shape[1] > 0:
            orphan_nonzero = int(np.sum(np.any(emb_orphan_part % PRIME != 0, axis=1)))
    else:
        orphan_nonzero = -1

    out_deg = int(sum(A[v]))

    return {
        'rank_istar': rank_istar,
        'out_deg': out_deg,
        'n_old_3': len(old_3_idx),
        'n_embedded_old': n_embedded,
        'n_orphan_old': n_orphan,
        'orphan_fraction': n_orphan / max(len(old_3_idx), 1),
        'rk_im_old': rk_im_old,
        'rk_im_embedded': rk_im_embedded,
        'rk_im_orphan': rk_im_orphan,
        'orphan_gen_nonzero': orphan_nonzero,
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"ORPHAN FACE ANALYSIS AT n={n}")
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
                r = orphan_analysis(A, n, v_cand)
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

        print(f"\n  Orphan old 3-paths:")
        print(f"    n_orphan: mean={np.mean([r['n_orphan_old'] for r in results]):.1f}")
        print(f"    n_embedded: mean={np.mean([r['n_embedded_old'] for r in results]):.1f}")
        print(f"    orphan fraction: mean={np.mean([r['orphan_fraction'] for r in results]):.3f}")

        print(f"\n  im(d_4^T) rank on different coord sets:")
        print(f"    old (all): mean={np.mean([r['rk_im_old'] for r in results]):.1f}")
        print(f"    embedded only: mean={np.mean([r['rk_im_embedded'] for r in results]):.1f}")
        print(f"    orphan only: mean={np.mean([r['rk_im_orphan'] for r in results]):.1f}")

        print(f"\n  Embedded gen has orphan-coord components: "
              f"{sum(1 for r in results if r['orphan_gen_nonzero'] > 0)}/{len(results)}")

        if fails:
            print(f"\n  FAILURES:")
            for r in fails[:5]:
                print(f"    out={r['out_deg']}: orphan_frac={r['orphan_fraction']:.3f}, "
                      f"rk_old={r['rk_im_old']}, rk_emb={r['rk_im_embedded']}, "
                      f"rk_orph={r['rk_im_orphan']}, gen_orphan={r['orphan_gen_nonzero']}")

        if succs:
            print(f"\n  SUCCESSES (sample):")
            for r in succs[:3]:
                print(f"    out={r['out_deg']}: orphan_frac={r['orphan_fraction']:.3f}, "
                      f"rk_old={r['rk_im_old']}, rk_emb={r['rk_im_embedded']}, "
                      f"rk_orph={r['rk_im_orphan']}, gen_orphan={r['orphan_gen_nonzero']}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
