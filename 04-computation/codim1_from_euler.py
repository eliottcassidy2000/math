"""
codim1_from_euler.py — Prove HYP-408 via Euler characteristic

HYP-408: codim(im(d_4^T)_old, ker(d_3^T)_old) = 1

Let π: A_3(T) → A_3(T)_old be the old projection (keep non-v paths, zero v paths).

Key: π is NOT a chain map (d ∘ π ≠ π ∘ d for through-v chains).
But the IMAGE of ker(d_3^T) under π has a specific dimension.

Define:
  K_p = ker(d_p^T) ∩ Ω_p(T)  [cycles in degree p]
  B_p = im(d_{p+1}^T) ∩ Ω_p(T)  [boundaries in degree p]

Then H_p(T) = K_p / B_p, dim = beta_p.

The old projection π maps K_3 → A_3(T)_old and B_3 → A_3(T)_old.
HYP-408 says: dim π(K_3) - dim π(B_3) = 1.

But we KNOW dim K_3 - dim B_3 = beta_3 = 1.
So HYP-408 says: π "preserves" the codimension.

This happens iff: ker(π|_{K_3}) ∩ (K_3 ∖ B_3) = ∅, i.e., no nonzero
element of K_3/B_3 is killed by π.

Equivalently: the H_3(T) generator has nonzero old-projection.
Equivalently: supp(h3_gen) ⊄ {through-v paths}.

We PROVED this computationally (0/1185 have vertex cover).

But can we prove it algebraically? The generator is a cycle —
a linear combination of 3-paths with zero boundary. If ALL paths
in its support pass through v, then v is in every support path.

A 3-path is (a,b,c,d) — 4 vertices. For v to be in it,
v ∈ {a,b,c,d}. So supp(h3_gen) ⊂ {through-v paths} means
every 3-path in the support uses vertex v.

Can we prove this is impossible for a nonzero cycle (∈ ker(d_3) ∖ im(d_4))?

Consider d_3 of a through-v 3-path (a,v,c,d) [v at position 1]:
  d_3 = (v,c,d) - (a,c,d) + (a,v,d) - (a,v,c)

The face (a,c,d) does NOT contain v — it's an old 2-path.
So d_3 of a through-v cycle has a nonzero old component (from faces where v is deleted).

For the cycle to be in ker(d_3), the sum of all these v-deletion faces must cancel.
So: if h3_gen is supported entirely on through-v paths,
then the sum of v-deletion 2-faces must be zero in A_2(T).

Let me check: is this sum-zero condition achievable?

Author: opus-2026-03-09-S58
"""
import sys
import time
import numpy as np
from collections import Counter, defaultdict
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


def through_v_cycle_test(A, n, v):
    """Check if any nonzero cycle in ker(d_3) is supported entirely on through-v paths."""
    max_p = min(n - 1, 6)
    ap = enumerate_all_allowed(A, n, max_p)

    paths_3 = ap.get(3, [])
    paths_2 = ap.get(2, [])

    if not paths_3:
        return None

    # Classify 3-paths
    through_v = [i for i, p in enumerate(paths_3) if v in p]
    not_through_v = [i for i, p in enumerate(paths_3) if v not in p]

    if not through_v:
        return {'has_through_v_cycle': False, 'reason': 'no_through_v_paths'}

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

    # d_3 in Omega coords
    idx2 = {p: i for i, p in enumerate(paths_2)}
    bd3 = np.zeros((len(paths_2), len(paths_3)), dtype=np.int64)
    for j, path in enumerate(paths_3):
        for sign, face in boundary_faces(path):
            if face in idx2:
                bd3[idx2[face], j] = (bd3[idx2[face], j] + sign) % PRIME

    d3o = bd3 @ ob3.T % PRIME

    # ker(d_3) in Omega coords
    if d3o.size > 0:
        _, kv = _gauss_nullbasis_modp(d3o.astype(np.int32), d3o.shape[0], d3o.shape[1], PRIME)
        ker_d3 = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    else:
        ker_d3 = np.eye(ob3.shape[0], dtype=np.int64)

    # ker(d_3) in A_3 coords
    ker_d3_A = ker_d3 @ ob3 % PRIME

    if ker_d3_A.shape[0] == 0:
        return {'has_through_v_cycle': False, 'reason': 'no_cycles'}

    # Check if any ker element is supported only on through-v paths
    through_v_set = set(through_v)
    not_through_v_set = set(not_through_v)

    # Restrict ker to non-through-v coordinates
    if not_through_v:
        ker_old = ker_d3_A[:, not_through_v] % PRIME
        rk_ker_old = int(_gauss_rank_np(ker_old.copy(), PRIME))
        ker_through_v_only = ker_d3_A.shape[0] - rk_ker_old  # dim of kernel of π|_{ker(d_3)}
    else:
        ker_through_v_only = ker_d3_A.shape[0]

    return {
        'has_through_v_cycle': ker_through_v_only > 0,
        'dim_ker_d3': ker_d3_A.shape[0],
        'dim_ker_d3_old_proj': rk_ker_old if not_through_v else 0,
        'dim_through_v_only': ker_through_v_only,
        'n_through_v': len(through_v),
        'n_not_through_v': len(not_through_v),
    }


def main():
    for n in [6, 7, 8]:
        print(f"\n{'='*70}")
        print(f"THROUGH-v CYCLE TEST AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        t0 = time.time()
        target = 1000 if n <= 7 else 300

        for trial in range(50000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) < 1:
                continue

            for v_cand in range(n):
                r = through_v_cycle_test(A, n, v_cand)
                if r is None:
                    continue
                r['trial'] = trial
                r['v'] = v_cand
                results.append(r)

        elapsed = time.time() - t0
        print(f"  {len(results)} (tournament, vertex) pairs, {elapsed:.1f}s")

        # How many have through-v-only cycles?
        has_tv = sum(1 for r in results if r.get('has_through_v_cycle', False))
        print(f"\n  Has through-v-only cycle: {has_tv}/{len(results)}")

        if has_tv > 0:
            print("  EXAMPLES:")
            for r in results:
                if r.get('has_through_v_cycle', False):
                    print(f"    trial={r['trial']}, v={r['v']}: "
                          f"dim_ker={r['dim_ker_d3']}, dim_old_proj={r['dim_ker_d3_old_proj']}, "
                          f"dim_tv_only={r['dim_through_v_only']}")
                    if sum(1 for r2 in results if r2.get('has_through_v_cycle', False)) > 5:
                        break

        # Distribution of dim_through_v_only
        tvonly_dist = Counter(r.get('dim_through_v_only', 0) for r in results)
        print(f"\n  dim(through-v-only cycles) distribution: {dict(sorted(tvonly_dist.items()))}")

        # KEY: dim(ker_d3) vs dim(ker_d3_old_proj)
        dim_diffs = [r['dim_ker_d3'] - r.get('dim_ker_d3_old_proj', 0)
                     for r in results if 'dim_ker_d3_old_proj' in r]
        if dim_diffs:
            diff_dist = Counter(dim_diffs)
            print(f"\n  dim(ker_d3) - dim(π(ker_d3)) distribution: {dict(sorted(diff_dist.items()))}")
            print(f"  (This is the number of ker(d_3) directions killed by old-projection)")


if __name__ == '__main__':
    main()
    print("\nDONE.")
