#!/usr/bin/env python3
"""
beta5_n9_check.py -- Check whether beta_5 > 0 appears for n=9 tournaments.

Pattern from THM-099:
  n=6: beta_3 appears (S-phase)
  n=7: beta_4 appears (=6 for BIBD maximizers)
  n=8: beta_4 appears (=1 for some maximizers)
  n=9: does beta_5 appear?

We sample random n=9 tournaments, compute Betti numbers up to dim 7,
and track the highest-H tournaments found.

kind-pasteur-2026-03-08-S34
"""

import numpy as np
from itertools import permutations
from collections import defaultdict
import random
import time
import sys

# ============================================================
# PATH HOMOLOGY FUNCTIONS (copied from path_homology_v2.py to
# avoid running its module-level validation code on import)
# ============================================================

def enumerate_allowed_paths(A, n, p):
    """All sequences of p+1 distinct vertices following directed edges."""
    if p < 0:
        return []
    if p == 0:
        return [(v,) for v in range(n)]
    paths = []
    for perm in permutations(range(n), p + 1):
        ok = True
        for i in range(p):
            if A[perm[i]][perm[i+1]] != 1:
                ok = False
                break
        if ok:
            paths.append(perm)
    return paths

def boundary_coeffs(path):
    """Returns list of (sign, face_tuple) for boundary(path)."""
    p = len(path) - 1
    result = []
    for i in range(p + 1):
        face = path[:i] + path[i+1:]
        result.append(((-1)**i, face))
    return result

def build_full_boundary_matrix(allowed_p, allowed_pm1):
    if not allowed_p or not allowed_pm1:
        return np.zeros((max(len(allowed_pm1), 0), max(len(allowed_p), 0)))
    idx_pm1 = {path: i for i, path in enumerate(allowed_pm1)}
    M = np.zeros((len(allowed_pm1), len(allowed_p)))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in idx_pm1:
                M[idx_pm1[face], j] += sign
    return M

def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0:
        return np.zeros((0, 0))
    if p == 0:
        return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed_faces = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed_faces:
                    non_allowed_faces[face] = na_count
                    na_count += 1
    if na_count == 0:
        return np.eye(dim_Ap)
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed_faces:
                P[non_allowed_faces[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = sum(s > 1e-10 for s in S)
    null_space = Vt[rank:].T
    if null_space.shape[1] == 0:
        return np.zeros((dim_Ap, 0))
    return null_space

def path_betti_numbers(A, n, max_dim=None):
    """Compute GLMY path Betti numbers beta_0, ..., beta_{max_dim}."""
    if max_dim is None:
        max_dim = n - 1
    allowed = {}
    for p in range(-1, max_dim + 2):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)
    omega = {}
    for p in range(max_dim + 2):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
    betti = []
    for p in range(max_dim + 1):
        dim_omega_p = omega[p].shape[1] if omega[p].ndim == 2 else 0
        if dim_omega_p == 0:
            betti.append(0)
            continue
        bd_p = build_full_boundary_matrix(allowed[p], allowed[p-1])
        bd_p_omega = bd_p @ omega[p]
        if bd_p_omega.shape[0] > 0 and bd_p_omega.shape[1] > 0:
            S_p = np.linalg.svd(bd_p_omega, compute_uv=False)
            rank_p = sum(s > 1e-8 for s in S_p)
        else:
            rank_p = 0
        ker_dim = dim_omega_p - rank_p
        dim_omega_p1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 else 0
        if dim_omega_p1 > 0:
            bd_p1 = build_full_boundary_matrix(allowed[p+1], allowed[p])
            bd_p1_omega = bd_p1 @ omega[p+1]
            S_p1 = np.linalg.svd(bd_p1_omega, compute_uv=False)
            im_dim = sum(s > 1e-8 for s in S_p1)
        else:
            im_dim = 0
        beta_p = ker_dim - im_dim
        betti.append(max(0, beta_p))
    return betti


# ============================================================
# TOURNAMENT UTILITIES
# ============================================================

def random_tournament(n):
    """Generate a random tournament on n vertices as adjacency matrix."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def ham_path_count(A, n):
    """Count Hamiltonian paths using DP."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_3cycles(A, n):
    """Count directed 3-cycles."""
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[j][i] and A[k][j] and A[i][k]):
                    t3 += 1
    return t3

def score_sequence(A, n):
    """Compute sorted out-degree sequence."""
    return tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))

def tournament_to_bits(A, n):
    """Encode tournament as integer for deduplication."""
    bits = 0
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if A[i][j]:
                bits |= (1 << idx)
            idx += 1
    return bits


# ============================================================
# MAIN COMPUTATION
# ============================================================

def main():
    n = 9
    max_dim = 7  # Compute beta_0 through beta_7

    # How many random tournaments to sample
    NUM_SAMPLES = 500
    # Keep top-K by H value
    TOP_K = 20

    print("=" * 70)
    print(f"BETA_5 CHECK FOR n={n} TOURNAMENTS")
    print(f"Sampling {NUM_SAMPLES} random tournaments, Betti up to dim {max_dim}")
    print("=" * 70)

    # First, time a single computation to gauge feasibility
    print("\nTiming single Betti computation at n=9, max_dim=7...")
    A_test = random_tournament(n)
    t0 = time.time()
    betti_test = path_betti_numbers(A_test, n, max_dim=max_dim)
    t1 = time.time()
    dt = t1 - t0
    print(f"  Single computation: {dt:.2f}s, beta = {betti_test}")
    print(f"  Estimated total for {NUM_SAMPLES}: {dt * NUM_SAMPLES / 60:.1f} min")

    # If too slow, reduce max_dim or sample count
    if dt > 60:
        print("  WARNING: Very slow! Reducing max_dim to 5 and samples to 100.")
        max_dim = 5
        NUM_SAMPLES = 100
    elif dt > 20:
        print("  WARNING: Slow. Reducing samples to 200.")
        NUM_SAMPLES = 200
    elif dt > 5:
        print("  Moderate speed. Keeping 500 samples.")
    else:
        print("  Fast! Increasing to 1000 samples.")
        NUM_SAMPLES = 1000

    sys.stdout.flush()

    # Phase 1: Quick H-value survey (H is fast, Betti is slow)
    # Sample many tournaments, find the ones with highest H
    print(f"\n--- Phase 1: H-value survey (10000 random tournaments) ---")
    h_survey = []
    t_start = time.time()
    for trial in range(10000):
        A = random_tournament(n)
        H = ham_path_count(A, n)
        bits = tournament_to_bits(A, n)
        h_survey.append((H, bits, A))
    t_phase1 = time.time() - t_start
    print(f"  H survey done in {t_phase1:.1f}s")

    h_values = [x[0] for x in h_survey]
    print(f"  H range: [{min(h_values)}, {max(h_values)}]")
    print(f"  H mean: {np.mean(h_values):.1f}, H std: {np.std(h_values):.1f}")

    # Sort by H descending, keep unique top ones
    h_survey.sort(key=lambda x: -x[0])
    top_unique = []
    seen_bits = set()
    for H, bits, A in h_survey:
        if bits not in seen_bits:
            seen_bits.add(bits)
            top_unique.append((H, bits, A))
        if len(top_unique) >= TOP_K:
            break

    print(f"\n  Top {len(top_unique)} unique H values:")
    for i, (H, bits, A) in enumerate(top_unique):
        sc = score_sequence(A, n)
        t3 = count_3cycles(A, n)
        print(f"    #{i+1}: H={H}, t3={t3}, score={sc}")
    sys.stdout.flush()

    # Phase 2: Compute Betti numbers for random sample + top-H tournaments
    print(f"\n--- Phase 2: Betti computation (max_dim={max_dim}) ---")

    # Track results
    betti_dist = defaultdict(int)
    beta5_positive = []
    beta6_positive = []
    beta7_positive = []  # only if max_dim >= 7
    all_results = []

    # First: compute Betti for all top-H tournaments
    print(f"\n  Computing Betti for top {len(top_unique)} H-maximizers...")
    for i, (H, bits, A) in enumerate(top_unique):
        t0 = time.time()
        betti = path_betti_numbers(A, n, max_dim=max_dim)
        dt = time.time() - t0
        bt = tuple(betti)
        betti_dist[bt] += 1
        sc = score_sequence(A, n)
        t3 = count_3cycles(A, n)
        all_results.append((H, bt, t3, sc, bits))

        print(f"    #{i+1}: H={H}, t3={t3}, beta={betti} ({dt:.1f}s)")

        if len(betti) > 5 and betti[5] > 0:
            beta5_positive.append((H, bt, t3, sc, bits))
        if len(betti) > 6 and betti[6] > 0:
            beta6_positive.append((H, bt, t3, sc, bits))
        if len(betti) > 7 and betti[7] > 0:
            beta7_positive.append((H, bt, t3, sc, bits))
        sys.stdout.flush()

    # Then: random sample
    print(f"\n  Computing Betti for {NUM_SAMPLES} random tournaments...")
    t_start = time.time()
    for trial in range(NUM_SAMPLES):
        A = random_tournament(n)
        H = ham_path_count(A, n)
        bits = tournament_to_bits(A, n)

        betti = path_betti_numbers(A, n, max_dim=max_dim)
        bt = tuple(betti)
        betti_dist[bt] += 1
        sc = score_sequence(A, n)
        t3 = count_3cycles(A, n)
        all_results.append((H, bt, t3, sc, bits))

        if len(betti) > 5 and betti[5] > 0:
            beta5_positive.append((H, bt, t3, sc, bits))
        if len(betti) > 6 and betti[6] > 0:
            beta6_positive.append((H, bt, t3, sc, bits))
        if len(betti) > 7 and betti[7] > 0:
            beta7_positive.append((H, bt, t3, sc, bits))

        if (trial + 1) % 50 == 0:
            elapsed = time.time() - t_start
            rate = (trial + 1) / elapsed
            eta = (NUM_SAMPLES - trial - 1) / rate
            print(f"    {trial+1}/{NUM_SAMPLES} done ({rate:.1f}/s, ETA {eta:.0f}s), "
                  f"beta5>0: {len(beta5_positive)}, "
                  f"beta6>0: {len(beta6_positive)}")
            sys.stdout.flush()

    total_time = time.time() - t_start
    total_computed = len(top_unique) + NUM_SAMPLES
    print(f"\n  Random sample done in {total_time:.1f}s ({NUM_SAMPLES/total_time:.1f}/s)")

    # ============================================================
    # RESULTS
    # ============================================================
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)

    print(f"\nTotal tournaments computed: {total_computed}")

    print(f"\n--- Betti distribution (all {total_computed} tournaments) ---")
    for bt in sorted(betti_dist.keys()):
        print(f"  beta={list(bt)}: {betti_dist[bt]} tournaments")

    print(f"\n--- beta_5 > 0 instances: {len(beta5_positive)} ---")
    if beta5_positive:
        for H, bt, t3, sc, bits in beta5_positive[:20]:
            print(f"  H={H}, t3={t3}, beta={list(bt)}, score={sc}")
    else:
        print("  NONE FOUND")

    print(f"\n--- beta_6 > 0 instances: {len(beta6_positive)} ---")
    if beta6_positive:
        for H, bt, t3, sc, bits in beta6_positive[:20]:
            print(f"  H={H}, t3={t3}, beta={list(bt)}, score={sc}")
    else:
        print("  NONE FOUND")

    if max_dim >= 7:
        print(f"\n--- beta_7 > 0 instances: {len(beta7_positive)} ---")
        if beta7_positive:
            for H, bt, t3, sc, bits in beta7_positive[:20]:
                print(f"  H={H}, t3={t3}, beta={list(bt)}, score={sc}")
        else:
            print("  NONE FOUND")

    # Correlations with H
    print(f"\n--- Betti vs H correlation ---")
    h_by_betti = defaultdict(list)
    for H, bt, t3, sc, bits in all_results:
        h_by_betti[bt].append(H)

    print(f"  {'Betti vector':<40s} {'count':>6s} {'H_min':>8s} {'H_mean':>8s} {'H_max':>8s}")
    for bt in sorted(h_by_betti.keys()):
        vals = h_by_betti[bt]
        print(f"  {str(list(bt)):<40s} {len(vals):>6d} {min(vals):>8d} {np.mean(vals):>8.1f} {max(vals):>8d}")

    # Summary
    print(f"\n--- SUMMARY ---")
    max_nonzero_dim = 0
    for bt in betti_dist.keys():
        for d in range(len(bt)-1, -1, -1):
            if bt[d] > 0:
                max_nonzero_dim = max(max_nonzero_dim, d)
                break
    print(f"  Highest dimension with beta_d > 0: {max_nonzero_dim}")

    if beta5_positive:
        print(f"  beta_5 > 0: YES ({len(beta5_positive)} instances)")
        print(f"    Confirms pattern beta_{{n-4}} > 0 for n=9")
    else:
        # Check for beta_4, beta_3 etc.
        b4_count = sum(1 for bt in betti_dist if len(bt) > 4 and bt[4] > 0)
        b3_count = sum(1 for bt in betti_dist if len(bt) > 3 and bt[3] > 0)
        print(f"  beta_5 > 0: NOT FOUND in {total_computed} samples")
        print(f"  beta_4 > 0: {b4_count} distinct Betti vectors")
        print(f"  beta_3 > 0: {b3_count} distinct Betti vectors")

    # Highest H found overall
    all_results.sort(key=lambda x: -x[0])
    print(f"\n  Top 5 by H:")
    for H, bt, t3, sc, bits in all_results[:5]:
        print(f"    H={H}, t3={t3}, beta={list(bt)}, score={sc}")

    print("\nDone.")

if __name__ == "__main__":
    main()
