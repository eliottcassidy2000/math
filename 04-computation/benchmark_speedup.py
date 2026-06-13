"""
benchmark_speedup.py — Compare old SVD vs new mod-p computation for path homology

Tests:
1. beta_3 computation: old (numpy SVD) vs new (mod-p Gauss)
2. beta_1 computation: old (full Omega_2) vs new (TT-boundary span)
3. Correctness verification: results must match
4. Scaling test: how timings grow with n

Author: kind-pasteur-S47 (2026-03-09)
"""
import sys
import time
import numpy as np
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

# Import both old-style and new utility functions
from tournament_utils import (
    bits_to_adj, random_tournament, enumerate_all_allowed,
    compute_omega_basis_numpy, boundary_faces,
    compute_betti_modp, compute_betti_hybrid, compute_beta1_fast,
    full_chain_complex
)


def compute_beta3_old(A, n):
    """Original SVD-based beta_3 (from beta3_good_vertex_characterization.py)."""
    ap = {}
    for p in range(min(6, n)):
        ap[p] = []
        if p == 0:
            ap[p] = [(v,) for v in range(n)]
            continue
        adj = [[] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if A[i][j] == 1:
                    adj[i].append(j)
        for start in range(n):
            stack = [([start], 1 << start)]
            while stack:
                path, visited = stack.pop()
                if len(path) == p + 1:
                    ap[p].append(tuple(path))
                    continue
                v = path[-1]
                for u in adj[v]:
                    if not (visited & (1 << u)):
                        stack.append((path + [u], visited | (1 << u)))

    omega_bases = {}
    omega_dims = {}
    for p in range(min(6, n)):
        if not ap.get(p, []):
            omega_bases[p] = np.zeros((0, 0))
            omega_dims[p] = 0
        else:
            omega_bases[p], omega_dims[p] = compute_omega_basis_numpy(ap, p)

    dim3 = omega_dims.get(3, 0)
    if dim3 == 0:
        return 0

    bd3 = np.zeros((len(ap.get(2, [])), len(ap[3])))
    idx2 = {path: i for i, path in enumerate(ap.get(2, []))}
    for j, path in enumerate(ap[3]):
        for sign, face in boundary_faces(path):
            if face in idx2:
                bd3[idx2[face], j] += sign

    O3 = omega_bases[3]
    d3_om = bd3 @ O3
    sv3 = np.linalg.svd(d3_om, compute_uv=False)
    rank_d3 = int(sum(s > 1e-8 for s in sv3))
    ker_d3 = dim3 - rank_d3

    if ker_d3 == 0:
        return 0

    dim4 = omega_dims.get(4, 0)
    if dim4 == 0:
        return ker_d3

    bd4 = np.zeros((len(ap[3]), len(ap.get(4, []))))
    idx3 = {path: i for i, path in enumerate(ap[3])}
    for j, path in enumerate(ap.get(4, [])):
        for sign, face in boundary_faces(path):
            if face in idx3:
                bd4[idx3[face], j] += sign

    O4 = omega_bases[4]
    d4_om = bd4 @ O4
    O3pinv = np.linalg.pinv(O3)
    d4_omega3 = O3pinv @ d4_om
    sv4 = np.linalg.svd(d4_omega3, compute_uv=False)
    rank_d4 = int(sum(s > 1e-8 for s in sv4))

    return ker_d3 - rank_d4


def compute_beta1_old(A, n):
    """Original SVD-based beta_1."""
    ap = enumerate_all_allowed(A, n, 3)
    omega_bases = {}
    omega_dims = {}
    for p in range(min(4, n)):
        if not ap.get(p, []):
            omega_bases[p] = np.zeros((0, 0))
            omega_dims[p] = 0
        else:
            omega_bases[p], omega_dims[p] = compute_omega_basis_numpy(ap, p)

    dim1 = omega_dims.get(1, 0)
    if dim1 == 0:
        return 0

    bd1 = np.zeros((len(ap.get(0, [])), len(ap[1])))
    idx0 = {path: i for i, path in enumerate(ap.get(0, []))}
    for j, path in enumerate(ap[1]):
        for sign, face in boundary_faces(path):
            if face in idx0:
                bd1[idx0[face], j] += sign

    O1 = omega_bases[1]
    d1_om = bd1 @ O1
    sv1 = np.linalg.svd(d1_om, compute_uv=False)
    rank_d1 = int(sum(s > 1e-8 for s in sv1))
    ker_d1 = dim1 - rank_d1

    if ker_d1 == 0:
        return 0

    dim2 = omega_dims.get(2, 0)
    if dim2 == 0:
        return ker_d1

    bd2 = np.zeros((len(ap[1]), len(ap[2])))
    idx1 = {path: i for i, path in enumerate(ap[1])}
    for j, path in enumerate(ap[2]):
        for sign, face in boundary_faces(path):
            if face in idx1:
                bd2[idx1[face], j] += sign

    O2 = omega_bases[2]
    d2_om = bd2 @ O2
    O1pinv = np.linalg.pinv(O1)
    d2_omega1 = O1pinv @ d2_om
    sv2 = np.linalg.svd(d2_omega1, compute_uv=False)
    rank_d2 = int(sum(s > 1e-8 for s in sv2))

    return ker_d1 - rank_d2


def main():
    print("=" * 70)
    print("BENCHMARK: Old SVD vs New Mod-p Path Homology Computation")
    print("=" * 70)

    # Part 1: Correctness verification at n=6 (exhaustive)
    print("\n--- Part 1: Correctness check at n=6 (all 32768 tournaments) ---")
    n = 6
    total = 2 ** (n*(n-1)//2)
    mismatches_b3 = 0
    mismatches_b1 = 0

    t0 = time.time()
    for bits in range(total):
        A = bits_to_adj(bits, n)
        old_b3 = compute_beta3_old(A, n)
        new_b3 = compute_betti_modp(A, n, 3, max_p=5)
        if old_b3 != new_b3:
            mismatches_b3 += 1
            print(f"  MISMATCH beta_3 at bits={bits}: old={old_b3} new={new_b3}")

        old_b1 = compute_beta1_old(A, n)
        new_b1 = compute_beta1_fast(A, n)
        if old_b1 != new_b1:
            mismatches_b1 += 1
            print(f"  MISMATCH beta_1 at bits={bits}: old={old_b1} new={new_b1}")

        if (bits + 1) % 5000 == 0:
            print(f"  {bits+1}/{total} checked...", flush=True)

    t1 = time.time()
    print(f"  Checked all {total}: beta_3 mismatches={mismatches_b3}, beta_1 mismatches={mismatches_b1}")
    print(f"  Total time: {t1-t0:.1f}s")

    # Part 2: Speed comparison at n=7
    print("\n--- Part 2: Speed comparison at n=7 (200 random tournaments) ---")
    n = 7
    N = 200
    rng = np.random.RandomState(42)

    # Old method - beta_3
    t0 = time.time()
    old_results_b3 = []
    for _ in range(N):
        A = random_tournament(n, rng)
        old_results_b3.append(compute_beta3_old(A, n))
    t_old_b3 = time.time() - t0

    # New method - beta_3
    rng2 = np.random.RandomState(42)  # same seed
    t0 = time.time()
    new_results_b3 = []
    for _ in range(N):
        A = random_tournament(n, rng2)
        new_results_b3.append(compute_betti_modp(A, n, 3, max_p=5))
    t_new_b3 = time.time() - t0

    # Hybrid method - beta_3
    rng2b = np.random.RandomState(42)
    t0 = time.time()
    hyb_results_b3 = []
    for _ in range(N):
        A = random_tournament(n, rng2b)
        hyb_results_b3.append(compute_betti_hybrid(A, n, 3, max_p=5))
    t_hyb_b3 = time.time() - t0

    match_b3 = sum(1 for a, b in zip(old_results_b3, new_results_b3) if a == b)
    match_hyb = sum(1 for a, b in zip(old_results_b3, hyb_results_b3) if a == b)
    print(f"  beta_3 old SVD:    {t_old_b3:.2f}s ({N} tours)")
    print(f"  beta_3 pure mod-p: {t_new_b3:.2f}s ({N} tours)")
    print(f"  beta_3 hybrid:     {t_hyb_b3:.2f}s ({N} tours)")
    print(f"  Speedup pure:   {t_old_b3/t_new_b3:.2f}x")
    print(f"  Speedup hybrid: {t_old_b3/t_hyb_b3:.2f}x")
    print(f"  Matches pure: {match_b3}/{N}, hybrid: {match_hyb}/{N}")

    # Old method - beta_1
    rng3 = np.random.RandomState(42)
    t0 = time.time()
    old_results_b1 = []
    for _ in range(N):
        A = random_tournament(n, rng3)
        old_results_b1.append(compute_beta1_old(A, n))
    t_old_b1 = time.time() - t0

    # New method - beta_1
    rng4 = np.random.RandomState(42)
    t0 = time.time()
    new_results_b1 = []
    for _ in range(N):
        A = random_tournament(n, rng4)
        new_results_b1.append(compute_beta1_fast(A, n))
    t_new_b1 = time.time() - t0

    match_b1 = sum(1 for a, b in zip(old_results_b1, new_results_b1) if a == b)
    print(f"\n  beta_1 old SVD:   {t_old_b1:.2f}s ({N} tours)")
    print(f"  beta_1 new TT:    {t_new_b1:.2f}s ({N} tours)")
    print(f"  Speedup: {t_old_b1/t_new_b1:.2f}x")
    print(f"  Matches: {match_b1}/{N}")

    # Part 3: full_chain_complex speed at n=7
    print("\n--- Part 3: full_chain_complex timing at n=7 ---")
    rng5 = np.random.RandomState(42)
    t0 = time.time()
    for _ in range(50):
        A = random_tournament(n, rng5)
        full_chain_complex(A, n, max_p=5)
    t_full = time.time() - t0
    print(f"  full_chain_complex (numpy SVD): {t_full:.2f}s for 50 tours")
    print(f"  Per tournament: {t_full/50*1000:.0f}ms")

    # New mod-p for same 50
    rng6 = np.random.RandomState(42)
    t0 = time.time()
    for _ in range(50):
        A = random_tournament(n, rng6)
        for p in range(6):
            compute_betti_modp(A, n, p, max_p=min(p+2, 5))
    t_modp_all = time.time() - t0
    print(f"  All bettis mod-p:               {t_modp_all:.2f}s for 50 tours")
    print(f"  Per tournament: {t_modp_all/50*1000:.0f}ms")

    # Part 4: Scaling with n (hybrid)
    print("\n--- Part 4: Per-tournament time scaling ---")
    for n_test in [5, 6, 7, 8]:
        rng_t = np.random.RandomState(42)
        count = 100 if n_test <= 7 else 20
        t0 = time.time()
        for _ in range(count):
            A = random_tournament(n_test, rng_t)
            compute_betti_hybrid(A, n_test, 3, max_p=5)
        t_hyb = time.time() - t0

        rng_t2 = np.random.RandomState(42)
        t0 = time.time()
        for _ in range(count):
            A = random_tournament(n_test, rng_t2)
            compute_beta3_old(A, n_test)
        t_old = time.time() - t0

        ratio = t_old / t_hyb if t_hyb > 0 else 0
        print(f"  n={n_test}: hybrid={t_hyb/count*1000:.1f}ms, SVD={t_old/count*1000:.1f}ms, ratio={ratio:.2f}x ({count} samples)")

    print("\nDONE.")


if __name__ == '__main__':
    main()
