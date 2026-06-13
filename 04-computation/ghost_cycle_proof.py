"""
ghost_cycle_proof.py — PROOF: Ghost Cycle Theorem ⟺ HYP-408

THEOREM (Ghost Cycle Equivalence):
  Let T be a tournament with beta_3(T) = 1 and v a vertex with
  both through-v and non-through-v 3-paths. Let π_old be the old-projection
  (zero out through-v coordinates). Then the following are equivalent:

  (A) HYP-408: codim(π_old(im(d_4)), π_old(ker(d_3))) = 1
  (B) Ghost Cycle Theorem: every through-v-only cycle is a boundary
      (i.e., ker(π_old|_{ker(d_3)}) ⊂ im(d_4))

PROOF:
  Let K = ker(d_3), B = im(d_4). Then K/B = H_3(T), dim = beta_3 = 1.

  Define:
    K_tv = ker(π_old|_K) = {z ∈ K : π_old(z) = 0}  [tv-only cycles]
    B_tv = ker(π_old|_B) = {z ∈ B : π_old(z) = 0}  [tv-only boundaries]

  Since B ⊂ K, we have B_tv ⊂ K_tv (as subspaces).

  By rank-nullity:
    dim(K_tv) = dim(K) - dim(π_old(K))
    dim(B_tv) = dim(B) - dim(π_old(B))

  Therefore:
    dim(K_tv) - dim(B_tv) = [dim(K) - dim(B)] - [dim(π_old(K)) - dim(π_old(B))]
                           = beta_3 - codim(π_old(B), π_old(K))

  (A) ⟹ (B):
    If codim(π_old(B), π_old(K)) = beta_3 = 1, then dim(K_tv) - dim(B_tv) = 0.
    Since B_tv ⊂ K_tv and dim(B_tv) = dim(K_tv), we get B_tv = K_tv.
    Therefore K_tv ⊂ B = im(d_4). ∎

  (B) ⟹ (A):
    If K_tv ⊂ B, then K_tv = K_tv ∩ B = B_tv.
    So dim(K_tv) = dim(B_tv), hence codim(π_old(B), π_old(K)) = beta_3 = 1. ∎

  MOREOVER: when beta_3 = 1:
    codim(π_old(B), π_old(K)) ∈ {0, 1}  always
    (since 0 ≤ codim ≤ beta_3 = 1)

  codim = 0 means π_old(B) = π_old(K), i.e., old-projection doesn't
  distinguish K and B. Then dim(K_tv) - dim(B_tv) = 1, so there's a
  tv-only cycle that's NOT a boundary — it represents the H_3 generator.

  codim = 1 means the H_3 generator has a "genuinely old" component,
  and all tv-only cycles are boundaries.

COMPUTATIONAL VERIFICATION:
  This script verifies both directions of the equivalence.

Author: opus-2026-03-09-S58
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament,
    enumerate_all_allowed,
    _build_constraint_matrix, _gauss_rank_np, _gauss_nullbasis_modp,
    full_chain_complex_modp, boundary_faces,
    RANK_PRIME
)

PRIME = RANK_PRIME


def verify_equivalence(A, n, v):
    """Verify Ghost Cycle Theorem ⟺ HYP-408."""
    max_p = min(n - 1, 6)
    ap = enumerate_all_allowed(A, n, max_p)

    paths_3 = ap.get(3, [])
    paths_4 = ap.get(4, [])
    paths_2 = ap.get(2, [])

    if not paths_3 or not paths_4:
        return None

    tv3 = [i for i, p in enumerate(paths_3) if v in p]
    old3 = [i for i, p in enumerate(paths_3) if v not in p]

    if not tv3 or not old3:
        return None

    def get_omega(ap, deg):
        paths = ap.get(deg, [])
        if not paths:
            return np.zeros((0, 0), dtype=np.int64)
        P, nr, nc = _build_constraint_matrix(ap, deg, PRIME)
        if P is not None:
            _, nb = _gauss_nullbasis_modp(P, nr, nc, PRIME)
            return np.array(nb, dtype=np.int64) if nb else np.zeros((0, nc), dtype=np.int64)
        return np.eye(len(paths), dtype=np.int64)

    ob3 = get_omega(ap, 3)
    ob4 = get_omega(ap, 4)

    # Build d_3, d_4
    idx2 = {p: i for i, p in enumerate(paths_2)}
    idx3 = {p: i for i, p in enumerate(paths_3)}

    bd3 = np.zeros((len(paths_2), len(paths_3)), dtype=np.int64)
    for j, path in enumerate(paths_3):
        for sign, face in boundary_faces(path):
            if face in idx2:
                bd3[idx2[face], j] = (bd3[idx2[face], j] + sign) % PRIME

    bd4 = np.zeros((len(paths_3), len(paths_4)), dtype=np.int64)
    for j, path in enumerate(paths_4):
        for sign, face in boundary_faces(path):
            if face in idx3:
                bd4[idx3[face], j] = (bd4[idx3[face], j] + sign) % PRIME

    d3o = bd3 @ ob3.T % PRIME

    # K = ker(d_3) in A_3 coords
    if d3o.size > 0:
        _, kv = _gauss_nullbasis_modp(d3o.astype(np.int32), d3o.shape[0], d3o.shape[1], PRIME)
        ker_d3 = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    else:
        ker_d3 = np.eye(ob3.shape[0], dtype=np.int64)

    K = ker_d3 @ ob3 % PRIME  # ker(d_3) in A_3 coords

    # B = im(d_4) in A_3 coords (rows)
    im_d4 = (bd4 @ ob4.T % PRIME).T if ob4.shape[0] > 0 else np.zeros((0, len(paths_3)), dtype=np.int64)

    dim_K = K.shape[0]
    dim_B = int(_gauss_rank_np(im_d4.copy(), PRIME)) if im_d4.shape[0] > 0 else 0
    beta3 = dim_K - dim_B

    if beta3 != 1:
        return None

    # π_old: project to old coordinates
    old3_sorted = sorted(old3)

    pi_K = K[:, old3_sorted] % PRIME  # π_old(K)
    rk_pi_K = int(_gauss_rank_np(pi_K.copy(), PRIME))

    pi_B = im_d4[:, old3_sorted] % PRIME if im_d4.shape[0] > 0 else np.zeros((0, len(old3_sorted)), dtype=np.int64)
    rk_pi_B = int(_gauss_rank_np(pi_B.copy(), PRIME)) if pi_B.shape[0] > 0 else 0

    codim_old = rk_pi_K - rk_pi_B  # = codim(π_old(B), π_old(K))

    # K_tv = ker(π_old|_K)
    dim_K_tv = dim_K - rk_pi_K

    # B_tv = ker(π_old|_B)
    dim_B_tv = dim_B - rk_pi_B

    # Check: dim(K_tv) - dim(B_tv) = beta3 - codim_old
    predicted_gap = beta3 - codim_old
    actual_gap = dim_K_tv - dim_B_tv

    # Ghost Cycle Theorem: K_tv ⊂ B? (check actual containment)
    if dim_K_tv > 0:
        # Find basis for K_tv
        pi_K_T = pi_K.T % PRIME
        _, tv_basis = _gauss_nullbasis_modp(
            pi_K_T.astype(np.int32), pi_K_T.shape[0], pi_K_T.shape[1], PRIME)
        tv_coeffs = np.array(tv_basis, dtype=np.int64) if tv_basis else np.zeros((0, K.shape[0]), dtype=np.int64)
        K_tv_vecs = tv_coeffs @ K % PRIME if tv_coeffs.shape[0] > 0 else np.zeros((0, len(paths_3)), dtype=np.int64)

        if K_tv_vecs.shape[0] > 0 and im_d4.shape[0] > 0:
            combined = np.vstack([im_d4, K_tv_vecs]) % PRIME
            rk_comb = int(_gauss_rank_np(combined.copy(), PRIME))
            ghost_contained = (rk_comb == dim_B)
        else:
            ghost_contained = (K_tv_vecs.shape[0] == 0)
    else:
        ghost_contained = True

    return {
        'dim_K': dim_K,
        'dim_B': dim_B,
        'beta3': beta3,
        'rk_pi_K': rk_pi_K,
        'rk_pi_B': rk_pi_B,
        'codim_old': codim_old,
        'dim_K_tv': dim_K_tv,
        'dim_B_tv': dim_B_tv,
        'predicted_gap': predicted_gap,
        'actual_gap': actual_gap,
        'gap_matches': predicted_gap == actual_gap,
        'ghost_contained': ghost_contained,
        'hyp408': codim_old == 1,
        'equivalence_holds': ghost_contained == (codim_old == 1),
    }


def main():
    print("=" * 70)
    print("PROOF VERIFICATION: Ghost Cycle Theorem ⟺ HYP-408")
    print("=" * 70)
    print()
    print("THEOREM: For beta_3(T)=1, the following are equivalent:")
    print("  (A) codim(π_old(im d_4), π_old(ker d_3)) = 1   [HYP-408]")
    print("  (B) ker(π_old|_{ker d_3}) ⊂ im(d_4)            [Ghost Cycle Thm]")
    print()
    print("PROOF: dim(K_tv) - dim(B_tv) = beta_3 - codim_old")
    print("       Since B_tv ⊂ K_tv, equal dims ⟹ equal subspaces.")
    print("       codim_old = 1 ⟺ dim gap = 0 ⟺ B_tv = K_tv ⟺ K_tv ⊂ B. ∎")
    print()

    for n in [6, 7, 8]:
        print(f"\n{'='*60}")
        print(f"VERIFICATION AT n={n}")
        print(f"{'='*60}")

        rng = np.random.RandomState(42)
        results = []
        t0 = time.time()
        target = 500 if n <= 7 else 400

        for trial in range(80000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) != 1:
                continue

            for v_cand in range(n):
                r = verify_equivalence(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)

        elapsed = time.time() - t0
        print(f"  {len(results)} (T,v) pairs, {elapsed:.1f}s")

        # Check the algebraic identity
        all_match = sum(1 for r in results if r['gap_matches'])
        print(f"\n  dim(K_tv) - dim(B_tv) = beta_3 - codim_old: {all_match}/{len(results)}")

        # Check equivalence
        equiv = sum(1 for r in results if r['equivalence_holds'])
        print(f"  Ghost Cycle ⟺ HYP-408: {equiv}/{len(results)}")

        # Distribution of codim_old
        codim_dist = Counter(r['codim_old'] for r in results)
        print(f"\n  codim(π_old(B), π_old(K)) distribution: {dict(sorted(codim_dist.items()))}")

        # When codim_old = 1 (HYP-408 holds)
        hyp408 = [r for r in results if r['codim_old'] == 1]
        print(f"  HYP-408 holds: {len(hyp408)}/{len(results)}")
        if hyp408:
            gc = sum(1 for r in hyp408 if r['ghost_contained'])
            print(f"    Ghost Cycle holds: {gc}/{len(hyp408)}")

        # When codim_old = 0 (HYP-408 fails)
        hyp408_fail = [r for r in results if r['codim_old'] == 0]
        print(f"  HYP-408 fails (codim=0): {len(hyp408_fail)}/{len(results)}")
        if hyp408_fail:
            gc_fail = sum(1 for r in hyp408_fail if r['ghost_contained'])
            print(f"    Ghost Cycle holds: {gc_fail}/{len(hyp408_fail)}")
            print(f"    Ghost Cycle FAILS: {len(hyp408_fail) - gc_fail}/{len(hyp408_fail)}")

            # Dimension analysis of failures
            print(f"\n    dim(K_tv) when codim_old=0: {Counter(r['dim_K_tv'] for r in hyp408_fail)}")
            print(f"    dim(B_tv) when codim_old=0: {Counter(r['dim_B_tv'] for r in hyp408_fail)}")

        # Summary: does this match what we expect?
        print(f"\n  SUMMARY n={n}:")
        print(f"    HYP-408 rate: {len(hyp408)}/{len(results)} = {100*len(hyp408)/max(len(results),1):.1f}%")
        print(f"    Ghost Cycle always holds when HYP-408 holds: {'YES' if all(r['ghost_contained'] for r in hyp408) else 'NO'}")
        if hyp408_fail:
            gc_when_fail = sum(1 for r in hyp408_fail if r['ghost_contained'])
            print(f"    Ghost Cycle when HYP-408 FAILS: {gc_when_fail}/{len(hyp408_fail)}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
