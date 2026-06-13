"""
h3_gen_support.py — Can the H_3(T) generator be supported entirely on through-v paths?

If the H_3(T) generator has nonzero old-coordinate projection for EVERY v,
then codim(im(d_4)_old, ker(d_3)_old) >= 1.
Combined with the fact that b3=1 (so the quotient has dim 1),
this gives codim = exactly 1 = HYP-408.

So HYP-408 reduces to: for every vertex v, the H_3(T) generator
has at least one non-v path with nonzero coefficient.

Equivalently: supp(h3_gen) is NOT contained in {through-v paths} for any v.
Equivalently: there's no v such that every path in supp(h3_gen) passes through v.

This is related to the "support vertex cover" — if some v is in every
support path, the generator is "concentrated at v".

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


def h3_generator_support(A, n):
    """Find the H_3(T) generator support and check vertex coverage."""
    max_p = min(n - 1, 6)

    cc = full_chain_complex_modp(A, n, max_p)
    if cc['bettis'].get(3, 0) != 1:
        return None

    ap = enumerate_all_allowed(A, n, max_p)
    paths_3 = ap.get(3, [])
    paths_4 = ap.get(4, [])
    paths_2 = ap.get(2, [])

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

    # ker(d_3)
    idx2 = {p: i for i, p in enumerate(paths_2)}
    bd3 = np.zeros((len(paths_2), len(paths_3)), dtype=np.int64)
    for j, path in enumerate(paths_3):
        for sign, face in boundary_faces(path):
            if face in idx2:
                bd3[idx2[face], j] = (bd3[idx2[face], j] + sign) % PRIME

    d3o = bd3 @ ob3.T % PRIME
    if d3o.size > 0:
        _, kv = _gauss_nullbasis_modp(d3o.astype(np.int32), d3o.shape[0], d3o.shape[1], PRIME)
        ker_d3 = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    else:
        ker_d3 = np.eye(ob3.shape[0], dtype=np.int64)

    ker_d3_A = ker_d3 @ ob3 % PRIME

    # im(d_4)
    idx3 = {p: i for i, p in enumerate(paths_3)}
    bd4 = np.zeros((len(paths_3), len(paths_4)), dtype=np.int64)
    for j, path in enumerate(paths_4):
        for sign, face in boundary_faces(path):
            if face in idx3:
                bd4[idx3[face], j] = (bd4[idx3[face], j] + sign) % PRIME

    if ob4.shape[0] > 0:
        im_d4 = (bd4 @ ob4.T % PRIME).T
    else:
        im_d4 = np.zeros((0, len(paths_3)), dtype=np.int64)

    rk_d4 = int(_gauss_rank_np(im_d4.copy(), PRIME)) if im_d4.shape[0] > 0 else 0

    # Find H_3 generator (ker element not in im)
    h3_gen = None
    for i in range(ker_d3_A.shape[0]):
        if im_d4.shape[0] > 0:
            test = np.vstack([im_d4, ker_d3_A[i:i+1]]) % PRIME
            rk = int(_gauss_rank_np(test.copy(), PRIME))
            if rk > rk_d4:
                h3_gen = ker_d3_A[i]
                break
        else:
            h3_gen = ker_d3_A[i]
            break

    if h3_gen is None:
        return None

    # Support analysis
    support_paths = [paths_3[j] for j in range(len(paths_3)) if h3_gen[j] % PRIME != 0]
    support_indices = [j for j in range(len(paths_3)) if h3_gen[j] % PRIME != 0]

    # For each vertex v, count support paths through v
    vertex_counts = {}
    for v in range(n):
        count = sum(1 for p in support_paths if v in p)
        vertex_counts[v] = count

    # Is there any v that's in ALL support paths?
    support_size = len(support_paths)
    vertex_cover = [v for v in range(n) if vertex_counts[v] == support_size]

    # For each v, what fraction of support goes through v?
    # The min over v of (support NOT through v) tells us the "old footprint"
    min_old_support = min(support_size - vertex_counts[v] for v in range(n))

    # Vertex with most support (most paths go through it)
    max_v = max(vertex_counts, key=vertex_counts.get)
    max_count = vertex_counts[max_v]

    return {
        'support_size': support_size,
        'total_paths_3': len(paths_3),
        'support_fraction': support_size / len(paths_3),
        'vertex_cover': vertex_cover,
        'has_cover': len(vertex_cover) > 0,
        'min_old_support': min_old_support,
        'max_through_v': max_count,
        'max_through_v_frac': max_count / support_size,
        'vertex_counts': vertex_counts,
        'scores': sorted([int(sum(A[i])) for i in range(n)]),
    }


def main():
    for n in [5, 6, 7, 8]:
        print(f"\n{'='*70}")
        print(f"H_3 GENERATOR SUPPORT ANALYSIS AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        target = 500 if n <= 7 else 200
        t0 = time.time()

        for trial in range(50000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            r = h3_generator_support(A, n)
            if r is None:
                continue
            results.append(r)

        elapsed = time.time() - t0
        print(f"  {len(results)} tournaments with b3=1, {elapsed:.1f}s")

        if not results:
            print("  No b3=1 tournaments found")
            continue

        # Has vertex cover?
        n_cover = sum(1 for r in results if r['has_cover'])
        print(f"\n  Generator has vertex cover (some v in ALL support paths): {n_cover}/{len(results)}")

        # Min old support
        min_olds = [r['min_old_support'] for r in results]
        print(f"  min_old_support (min over v of #support not through v):")
        print(f"    min={min(min_olds)}, max={max(min_olds)}, mean={np.mean(min_olds):.1f}")
        print(f"    distribution: {dict(sorted(Counter(min_olds).items()))}" if len(set(min_olds)) <= 20 else "")

        # Support size
        supp = [r['support_size'] for r in results]
        print(f"\n  Support size: min={min(supp)}, max={max(supp)}, mean={np.mean(supp):.1f}")
        print(f"  Total 3-paths: mean={np.mean([r['total_paths_3'] for r in results]):.1f}")
        print(f"  Support fraction: mean={np.mean([r['support_fraction'] for r in results]):.3f}")

        # Max through-v fraction
        max_thru = [r['max_through_v_frac'] for r in results]
        print(f"\n  Max through-v fraction: min={min(max_thru):.3f}, max={max(max_thru):.3f}, "
              f"mean={np.mean(max_thru):.3f}")

        # Per-vertex analysis: are out-degree extremes more likely to be "cover" vertices?
        if n >= 7:
            # Check which vertex has highest through-v count
            max_v_outdeg = []
            for r in results:
                max_v = max(r['vertex_counts'], key=r['vertex_counts'].get)
                out_deg = int(sum(A_dummy for A_dummy in [0]))  # can't access A here
                max_v_outdeg.append(max_v)  # just record vertex ID


if __name__ == '__main__':
    main()
    print("\nDONE.")
