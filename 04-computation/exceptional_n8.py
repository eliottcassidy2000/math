"""
exceptional_n8.py — Study the rare n=8 tournaments where rank(i_*)=0

At n=8 with β_3=1, about 0.6% of (T,v) pairs have rank(i_*)=0,
meaning H_4^rel(T,T\v) = 1. What makes these tournaments special?

Questions:
1. Are these always the SAME tournament (up to isomorphism)?
2. What is the score sequence of the exceptional tournaments?
3. What is the structure of the T\v deletion?
4. Do they have unusual Ω_4, Ω_5 dimensions?
5. What is the in-degree/out-degree of v in these cases?

Author: opus-2026-03-10-S59
"""
import sys
import time
import numpy as np
from collections import Counter
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


def score_sequence(A, n):
    """Sorted out-degrees."""
    return tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))


def check_istar(A, n, v):
    """Check rank(i_*) using the chain criterion."""
    max_p = min(n - 1, 6)
    ap_T = enumerate_all_allowed(A, n, max_p)

    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]
    remaining_inv = {remaining[i]: i for i in range(n1)}
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    cc_Tv = full_chain_complex_modp(A_sub, n1, min(n1-1, 5))
    b3_Tv = cc_Tv['bettis'].get(3, 0)
    if b3_Tv != 1:
        return None

    paths = {p: ap_T.get(p, []) for p in range(2, max_p + 1)}
    paths_Tv = {p: ap_Tv.get(p, []) for p in range(2, min(max_p, n1-1) + 1)}

    tv4 = [i for i, p in enumerate(paths.get(4, [])) if v in p]
    old3 = [i for i, p in enumerate(paths.get(3, [])) if v not in p]
    tv3 = [i for i, p in enumerate(paths.get(3, [])) if v in p]

    if not tv4:
        return None

    def get_omega(ap, deg):
        ps = ap.get(deg, [])
        if not ps:
            return np.zeros((0, 0), dtype=np.int64)
        P, nr, nc = _build_constraint_matrix(ap, deg, PRIME)
        if P is not None:
            _, nb = _gauss_nullbasis_modp(P, nr, nc, PRIME)
            return np.array(nb, dtype=np.int64) if nb else np.zeros((0, nc), dtype=np.int64)
        return np.eye(len(ps), dtype=np.int64)

    ob4_T = get_omega(ap_T, 4)
    ob4_Tv = get_omega(ap_Tv, 4)
    ob3_Tv = get_omega(ap_Tv, 3)

    idx3_T = {p: i for i, p in enumerate(paths.get(3, []))}
    idx3_Tv = {p: i for i, p in enumerate(paths_Tv.get(3, []))}
    idx2_Tv = {p: i for i, p in enumerate(paths_Tv.get(2, []))}

    bd4_T = np.zeros((len(paths.get(3, [])), len(paths.get(4, []))), dtype=np.int64)
    for j, path in enumerate(paths.get(4, [])):
        for sign, face in boundary_faces(path):
            if face in idx3_T:
                bd4_T[idx3_T[face], j] = (bd4_T[idx3_T[face], j] + sign) % PRIME

    d4_tvtv = bd4_T[np.ix_(tv3, tv4)] % PRIME
    d4_tvtv_omega = d4_tvtv @ ob4_T[:, tv4].T % PRIME if ob4_T.shape[0] > 0 else np.zeros((0, 0), dtype=np.int64)

    if d4_tvtv_omega.size > 0 and d4_tvtv_omega.shape[1] > 0:
        _, kv = _gauss_nullbasis_modp(
            d4_tvtv_omega.astype(np.int32), d4_tvtv_omega.shape[0], d4_tvtv_omega.shape[1], PRIME)
        ker_tvtv = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob4_T.shape[0]), dtype=np.int64)
    else:
        ker_tvtv = np.eye(ob4_T.shape[0], dtype=np.int64)

    ker_A4 = ker_tvtv @ ob4_T % PRIME if ker_tvtv.shape[0] > 0 else np.zeros((0, len(paths.get(4, []))), dtype=np.int64)

    # Full boundary projected to old3
    psi_full = (bd4_T[old3, :] @ ker_A4.T % PRIME).T if ker_A4.shape[0] > 0 else np.zeros((0, len(old3)), dtype=np.int64)

    # Map to T\v coords
    old3_to_Tv = {}
    for local_j, global_i in enumerate(old3):
        path_T = paths[3][global_i]
        path_Tv = tuple(remaining_inv[x] for x in path_T)
        if path_Tv in idx3_Tv:
            old3_to_Tv[local_j] = idx3_Tv[path_Tv]

    psi_Tv = np.zeros((psi_full.shape[0], len(paths_Tv.get(3, []))), dtype=np.int64)
    for local_j in range(len(old3)):
        if local_j in old3_to_Tv:
            j_Tv = old3_to_Tv[local_j]
            psi_Tv[:, j_Tv] = (psi_Tv[:, j_Tv] + psi_full[:, local_j]) % PRIME

    # im(d_4^{T\v})
    bd4_Tv = np.zeros((len(paths_Tv.get(3, [])), len(paths_Tv.get(4, []))), dtype=np.int64)
    for j, path in enumerate(paths_Tv.get(4, [])):
        for sign, face in boundary_faces(path):
            if face in idx3_Tv:
                bd4_Tv[idx3_Tv[face], j] = (bd4_Tv[idx3_Tv[face], j] + sign) % PRIME

    im_d4_Tv = (bd4_Tv @ ob4_Tv.T % PRIME).T if ob4_Tv.shape[0] > 0 else np.zeros((0, len(paths_Tv.get(3, []))), dtype=np.int64)
    rk_d4_Tv = int(_gauss_rank_np(im_d4_Tv.copy(), PRIME)) if im_d4_Tv.shape[0] > 0 else 0

    # Check containment
    if psi_Tv.shape[0] > 0 and im_d4_Tv.shape[0] > 0:
        combined = np.vstack([im_d4_Tv, psi_Tv]) % PRIME
        rk_combined = int(_gauss_rank_np(combined.copy(), PRIME))
        rank_istar = 0 if rk_combined > rk_d4_Tv else 1
    elif psi_Tv.shape[0] > 0:
        rk_psi = int(_gauss_rank_np(psi_Tv.copy(), PRIME))
        rank_istar = 0 if rk_psi > 0 else 1
    else:
        rank_istar = 1  # trivially injective

    # Collect structural data
    outdeg_v = sum(A[v][j] for j in range(n))
    indeg_v = sum(A[j][v] for j in range(n))

    return {
        'rank_istar': rank_istar,
        'score_seq': score_sequence(A, n),
        'score_seq_Tv': score_sequence(A_sub, n1),
        'outdeg_v': outdeg_v,
        'indeg_v': indeg_v,
        'dim_Omega4': ob4_T.shape[0],
        'dim_Omega5': get_omega(ap_T, 5).shape[0] if ap_T.get(5) else 0,
        'dim_tv4': len(tv4),
        'n_paths_3': len(paths.get(3, [])),
        'n_paths_4': len(paths.get(4, [])),
        'n_paths_5': len(paths.get(5, [])),
        'b3_Tv': b3_Tv,
    }


def main():
    print("FINDING EXCEPTIONAL n=8 TOURNAMENTS WITH rank(i_*)=0")
    print("="*70)

    rng = np.random.RandomState(42)
    exceptional = []
    normal = []
    t0 = time.time()

    for trial in range(200000):
        if len(exceptional) >= 20:
            break
        if trial % 1000 == 0 and trial > 0:
            print(f"  trial {trial}, found {len(exceptional)} exceptional, {len(normal)} normal")

        A = random_tournament(8, rng)
        cc = full_chain_complex_modp(A, 8, 7)
        if cc['bettis'].get(3, 0) != 1:
            continue

        for v in range(8):
            r = check_istar(A, 8, v)
            if r is None:
                continue

            if r['rank_istar'] == 0:
                r['A'] = [row[:] for row in A]
                r['v'] = v
                exceptional.append(r)
            elif len(normal) < 100:
                normal.append(r)

    elapsed = time.time() - t0
    print(f"\n  Found {len(exceptional)} exceptional cases in {elapsed:.1f}s")

    # Analyze exceptional vs normal
    print(f"\n{'='*70}")
    print(f"EXCEPTIONAL (rank i_* = 0):")
    print(f"{'='*70}")
    if exceptional:
        score_seqs = Counter(r['score_seq'] for r in exceptional)
        print(f"  Score sequences: {dict(score_seqs)}")

        score_seqs_Tv = Counter(r['score_seq_Tv'] for r in exceptional)
        print(f"  Score sequences T\\v: {dict(score_seqs_Tv)}")

        outdeg_dist = Counter(r['outdeg_v'] for r in exceptional)
        print(f"  Out-degree of v: {dict(sorted(outdeg_dist.items()))}")

        for key in ['dim_Omega4', 'dim_Omega5', 'dim_tv4', 'n_paths_3', 'n_paths_4', 'n_paths_5']:
            vals = [r[key] for r in exceptional]
            print(f"  {key}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")

    print(f"\n{'='*70}")
    print(f"NORMAL (rank i_* = 1):")
    print(f"{'='*70}")
    if normal:
        score_seqs = Counter(r['score_seq'] for r in normal)
        print(f"  Score sequences (top 5): {dict(score_seqs.most_common(5))}")

        outdeg_dist = Counter(r['outdeg_v'] for r in normal)
        print(f"  Out-degree of v: {dict(sorted(outdeg_dist.items()))}")

        for key in ['dim_Omega4', 'dim_Omega5', 'dim_tv4', 'n_paths_3', 'n_paths_4', 'n_paths_5']:
            vals = [r[key] for r in normal]
            print(f"  {key}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")

    # Are the exceptional tournaments isomorphic?
    if len(exceptional) >= 2:
        print(f"\n{'='*70}")
        print(f"ARE EXCEPTIONAL TOURNAMENTS ISOMORPHIC?")
        print(f"{'='*70}")
        for i, r in enumerate(exceptional[:10]):
            ss = r['score_seq']
            dv = r['outdeg_v']
            print(f"  #{i}: score_seq={ss}, v-outdeg={dv}, Ω4={r['dim_Omega4']}, Ω5={r['dim_Omega5']}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
