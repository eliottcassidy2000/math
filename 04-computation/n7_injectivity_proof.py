"""
n7_injectivity_proof.py — Investigate WHY i_* is always injective at n=7

Strategy: At n=7 with b3(T)=1, the key chain complex dimensions are:
  Ω_3(T) ~ 24-28, Ω_4(T) ~ 15-18
  ker(d_3^T) = rk(d_4) + 1 (since b3=1)

The question is: why can't the embedded H_3(T\v) generator fall in im(d_4)?

Approach 1: Check if ALL old 3-cycles (not through v) are in im(d_4^old)
  where d_4^old means boundaries of 4-chains not through v.
  If old ker(d_3) ⊂ im(d_4^old) + span(h3_gen), then any old cycle
  either IS the generator or is a boundary.

Approach 2: Check the dimension of "old" vs "new" parts of ker(d_3)
  and im(d_4).

Approach 3: Direct computation — for each BAD vertex, decompose
  ker(d_3^T) = im(d_4^old) ⊕ im(d_4^new) ⊕ span(h3_gen)
  and check where the embedded H_3(T\v) gen lives.

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
    enumerate_all_allowed,
    _build_constraint_matrix, _gauss_rank_np, _gauss_nullbasis_modp,
    full_chain_complex_modp, boundary_faces,
    RANK_PRIME
)

PRIME = RANK_PRIME


def decompose_ker_d3(A, n, v, verbose=False):
    """Decompose ker(d_3^T) into old/new parts relative to vertex v."""
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
    paths_2_T = ap_T.get(2, [])
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

    # Classify 3-paths and 4-paths as "old" (not through v) or "new" (through v)
    old_3_idx = []
    new_3_idx = []
    for i, p in enumerate(paths_3_T):
        if v in p:
            new_3_idx.append(i)
        else:
            old_3_idx.append(i)

    old_4_idx = []
    new_4_idx = []
    for i, p in enumerate(paths_4_T):
        if v in p:
            new_4_idx.append(i)
        else:
            old_4_idx.append(i)

    # d_3^T in A_3 -> A_2 coords
    idx2T = {p: i for i, p in enumerate(paths_2_T)}
    bd3_T = np.zeros((len(paths_2_T), len(paths_3_T)), dtype=np.int64)
    for j, path in enumerate(paths_3_T):
        for sign, face in boundary_faces(path):
            if face in idx2T:
                bd3_T[idx2T[face], j] = (bd3_T[idx2T[face], j] + sign) % PRIME

    # ker(d_3^T) in Omega_3 coords
    d3o_T = bd3_T @ ob3_T.T % PRIME
    if d3o_T.size > 0:
        _, ker_vecs = _gauss_nullbasis_modp(
            d3o_T.astype(np.int32), d3o_T.shape[0], d3o_T.shape[1], PRIME)
        ker_d3T_omega = np.array(ker_vecs, dtype=np.int64) if ker_vecs else np.zeros((0, ob3_T.shape[0]), dtype=np.int64)
    else:
        ker_d3T_omega = np.eye(ob3_T.shape[0], dtype=np.int64)

    # ker(d_3^T) in A_3 coords
    ker_d3T_A = ker_d3T_omega @ ob3_T % PRIME

    # d_4^T in A_4 -> A_3 coords
    idx3T = {p: i for i, p in enumerate(paths_3_T)}
    bd4_T = np.zeros((len(paths_3_T), len(paths_4_T)), dtype=np.int64)
    for j, path in enumerate(paths_4_T):
        for sign, face in boundary_faces(path):
            if face in idx3T:
                bd4_T[idx3T[face], j] = (bd4_T[idx3T[face], j] + sign) % PRIME

    # im(d_4^T) in A_3 coords
    if ob4_T.shape[0] > 0:
        im_d4_A = (bd4_T @ ob4_T.T % PRIME).T
    else:
        im_d4_A = np.zeros((0, len(paths_3_T)), dtype=np.int64)

    rk_d4 = int(_gauss_rank_np(im_d4_A.copy(), PRIME)) if im_d4_A.shape[0] > 0 else 0

    # === OLD part of im(d_4) ===
    # Boundaries of old 4-chains (not through v)
    # These are linear combos of ob4_T rows restricted to old-4 paths
    # Actually: old 4-chains in Omega_4(T) are Omega basis vectors
    # whose support is entirely on non-v paths.
    # But Omega basis vectors can mix old and new paths...
    # Instead: compute im(d_4) from old allowed 4-paths only

    # Simpler: compute d_4 applied to each Omega_4 basis vector,
    # classify as "old" if the Omega_4 vector has support only on old 4-paths
    old_omega4_rows = []
    new_omega4_rows = []
    for row_i in range(ob4_T.shape[0]):
        row = ob4_T[row_i]
        has_new = any(row[j] % PRIME != 0 for j in new_4_idx)
        if has_new:
            new_omega4_rows.append(row_i)
        else:
            old_omega4_rows.append(row_i)

    # im(d_4^old) = boundaries of old Omega_4 vectors
    if old_omega4_rows:
        old_ob4 = ob4_T[old_omega4_rows]
        im_d4_old = (bd4_T @ old_ob4.T % PRIME).T
    else:
        im_d4_old = np.zeros((0, len(paths_3_T)), dtype=np.int64)

    # im(d_4^new) = boundaries of new Omega_4 vectors
    if new_omega4_rows:
        new_ob4 = ob4_T[new_omega4_rows]
        im_d4_new = (bd4_T @ new_ob4.T % PRIME).T
    else:
        im_d4_new = np.zeros((0, len(paths_3_T)), dtype=np.int64)

    rk_d4_old = int(_gauss_rank_np(im_d4_old.copy(), PRIME)) if im_d4_old.shape[0] > 0 else 0
    rk_d4_new = int(_gauss_rank_np(im_d4_new.copy(), PRIME)) if im_d4_new.shape[0] > 0 else 0

    # === Embed ker(d_3^{T\v}) ===
    idx2Tv = {p: i for i, p in enumerate(paths_2_Tv)}
    bd3_Tv = np.zeros((len(paths_2_Tv), len(paths_3_Tv)), dtype=np.int64)
    for j, path in enumerate(paths_3_Tv):
        for sign, face in boundary_faces(path):
            if face in idx2Tv:
                bd3_Tv[idx2Tv[face], j] = (bd3_Tv[idx2Tv[face], j] + sign) % PRIME

    d3o_Tv = bd3_Tv @ ob3_Tv.T % PRIME
    if d3o_Tv.size > 0:
        _, ker_vecs_Tv = _gauss_nullbasis_modp(
            d3o_Tv.astype(np.int32), d3o_Tv.shape[0], d3o_Tv.shape[1], PRIME)
        ker_d3Tv = np.array(ker_vecs_Tv, dtype=np.int64) if ker_vecs_Tv else np.zeros((0, ob3_Tv.shape[0]), dtype=np.int64)
    else:
        ker_d3Tv = np.eye(ob3_Tv.shape[0], dtype=np.int64)

    ker_d3Tv_A = ker_d3Tv @ ob3_Tv % PRIME if ker_d3Tv.shape[0] > 0 else np.zeros((0, len(paths_3_Tv)), dtype=np.int64)

    embed = np.zeros((len(paths_3_Tv), len(paths_3_T)), dtype=np.int64)
    for jv, path_v in enumerate(paths_3_Tv):
        mapped = tuple(remaining[x] for x in path_v)
        if mapped in idx3T:
            embed[jv, idx3T[mapped]] = 1

    embedded_ker = ker_d3Tv_A @ embed % PRIME if ker_d3Tv_A.shape[0] > 0 else np.zeros((0, len(paths_3_T)), dtype=np.int64)

    # rank(i_*)
    if im_d4_A.shape[0] > 0 and embedded_ker.shape[0] > 0:
        combined = np.vstack([im_d4_A, embedded_ker]) % PRIME
    elif embedded_ker.shape[0] > 0:
        combined = embedded_ker
    else:
        combined = im_d4_A
    rk_combined = int(_gauss_rank_np(combined.copy(), PRIME)) if combined.shape[0] > 0 else 0
    rank_istar = rk_combined - rk_d4

    # === KEY TEST: Is the embedded ker(d_3^{T\v}) generator in im(d_4^old)? ===
    if im_d4_old.shape[0] > 0 and embedded_ker.shape[0] > 0:
        test_old = np.vstack([im_d4_old, embedded_ker]) % PRIME
        rk_test_old = int(_gauss_rank_np(test_old.copy(), PRIME))
        emb_in_old = (rk_test_old == rk_d4_old)
    else:
        emb_in_old = False

    # Is it in im(d_4^old) + im(d_4^new)? (= im(d_4) full)
    # This is just rank_istar == 0

    # What fraction of ker(d_3^T) is spanned by old cycles?
    # "Old cycles" = ker(d_3) vectors supported only on old 3-paths
    old_ker_rows = []
    for i in range(ker_d3T_A.shape[0]):
        row = ker_d3T_A[i]
        has_new = any(row[j] % PRIME != 0 for j in new_3_idx)
        if not has_new:
            old_ker_rows.append(i)

    n_old_ker = len(old_ker_rows)

    out_deg = int(sum(A[v]))

    return {
        'rank_istar': rank_istar,
        'out_deg': out_deg,
        'rk_d4': rk_d4,
        'rk_d4_old': rk_d4_old,
        'rk_d4_new': rk_d4_new,
        'dim_ker_d3T': ker_d3T_A.shape[0],
        'n_old_ker': n_old_ker,
        'n_old_3': len(old_3_idx),
        'n_new_3': len(new_3_idx),
        'n_old_4': len(old_4_idx),
        'n_new_4': len(new_4_idx),
        'n_old_omega4': len(old_omega4_rows),
        'n_new_omega4': len(new_omega4_rows),
        'emb_in_old_imd4': emb_in_old,
        'omega3_T': ob3_T.shape[0],
        'omega4_T': ob4_T.shape[0],
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"OLD/NEW DECOMPOSITION OF ker(d_3) AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        checked = 0
        target = 200 if n == 7 else 100
        t0 = time.time()

        for trial in range(20000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            checked += 1

            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) != 1:
                continue

            for v_cand in range(n):
                r = decompose_ker_d3(A, n, v_cand)
                if r is None:
                    continue
                r['trial'] = trial
                r['v'] = v_cand
                results.append(r)

            if len(results) % 50 == 0 and len(results) > 0:
                elapsed = time.time() - t0
                print(f"  Progress: {len(results)} BAD vertices, {elapsed:.1f}s")

        elapsed = time.time() - t0
        print(f"  {len(results)} BAD vertices collected from {checked} tournaments, {elapsed:.1f}s")

        fails = [r for r in results if r['rank_istar'] == 0]
        succs = [r for r in results if r['rank_istar'] == 1]
        print(f"  Failures: {len(fails)}, Successes: {len(succs)}")

        # Key decomposition stats
        print(f"\n  --- Dimension decomposition (all BAD vertices) ---")
        print(f"  dim ker(d3):   min={min(r['dim_ker_d3T'] for r in results)}, "
              f"max={max(r['dim_ker_d3T'] for r in results)}, "
              f"mean={np.mean([r['dim_ker_d3T'] for r in results]):.1f}")
        print(f"  rk(d4):        min={min(r['rk_d4'] for r in results)}, "
              f"max={max(r['rk_d4'] for r in results)}, "
              f"mean={np.mean([r['rk_d4'] for r in results]):.1f}")
        print(f"  rk(d4_old):    min={min(r['rk_d4_old'] for r in results)}, "
              f"max={max(r['rk_d4_old'] for r in results)}, "
              f"mean={np.mean([r['rk_d4_old'] for r in results]):.1f}")
        print(f"  rk(d4_new):    min={min(r['rk_d4_new'] for r in results)}, "
              f"max={max(r['rk_d4_new'] for r in results)}, "
              f"mean={np.mean([r['rk_d4_new'] for r in results]):.1f}")
        print(f"  rk(d4_old)+rk(d4_new) vs rk(d4):")
        diffs = [r['rk_d4_old'] + r['rk_d4_new'] - r['rk_d4'] for r in results]
        print(f"    difference: min={min(diffs)}, max={max(diffs)}, mean={np.mean(diffs):.2f}")
        print(f"    (>0 means old and new boundaries OVERLAP)")

        print(f"\n  n_old_ker (ker(d3) vectors with pure old support):")
        old_ker_vals = [r['n_old_ker'] for r in results]
        print(f"    min={min(old_ker_vals)}, max={max(old_ker_vals)}, mean={np.mean(old_ker_vals):.1f}")
        print(f"    total dim_ker: mean={np.mean([r['dim_ker_d3T'] for r in results]):.1f}")
        print(f"    fraction old: {np.mean([r['n_old_ker']/max(r['dim_ker_d3T'],1) for r in results]):.3f}")

        print(f"\n  Omega_4 old/new split:")
        print(f"    old omega4: mean={np.mean([r['n_old_omega4'] for r in results]):.1f}")
        print(f"    new omega4: mean={np.mean([r['n_new_omega4'] for r in results]):.1f}")
        print(f"    total: mean={np.mean([r['omega4_T'] for r in results]):.1f}")

        # KEY: embedded gen in old im(d4)?
        emb_in_old_count = sum(1 for r in results if r['emb_in_old_imd4'])
        print(f"\n  Embedded H3(Tv) gen IN im(d4_old): {emb_in_old_count}/{len(results)} "
              f"({100*emb_in_old_count/max(len(results),1):.1f}%)")

        if fails:
            print(f"\n  --- FAILURE-specific decomposition ---")
            for r in fails[:5]:
                print(f"    trial={r['trial']}, v={r['v']}, out_deg={r['out_deg']}")
                print(f"      dim_ker={r['dim_ker_d3T']}, rk_d4={r['rk_d4']}, "
                      f"rk_old={r['rk_d4_old']}, rk_new={r['rk_d4_new']}")
                print(f"      n_old_ker={r['n_old_ker']}, emb_in_old={r['emb_in_old_imd4']}")
                overlap = r['rk_d4_old'] + r['rk_d4_new'] - r['rk_d4']
                print(f"      old+new overlap = {overlap}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
