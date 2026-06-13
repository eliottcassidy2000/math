#!/usr/bin/env python3
"""
beta3_good_vertex_exhaustive_n7.py — opus-2026-03-09-S53

EXHAUSTIVE check at n=7: does every tournament have a good vertex?
(A vertex v where beta_3(T\v) = 0)

Strategy: for each of 2^21 = 2,097,152 tournaments T on 7 vertices,
check if there exists v with beta_3(T\v) = 0.

Optimization: beta_3(T\v) computation at n=6 is fast.
Short-circuit: stop checking deletions once we find beta_3=0.
Also: precompute which 6-vertex tournaments have beta_3=1 (only 320/32768).

The 320 beta_3=1 tournaments at n=6 are exactly:
  - 80 with score (1,1,1,4,4,4)
  - 240 with score (2,2,2,3,3,3)

So we can first check score sequences of deletions. If any deletion
has a score NOT in {(1,1,1,4,4,4), (2,2,2,3,3,3)}, it's automatically good.
This gives a very fast filter.
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter
import time

def tournament_from_bits(n, bits):
    T = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                T[i][j] = 1
            else:
                T[j][i] = 1
            idx += 1
    return T

def is_allowed_path(T, path):
    for i in range(len(path)-1):
        if not T[path[i]][path[i+1]]:
            return False
    return len(path) == len(set(path))

def compute_beta3_n6(T):
    """Optimized beta_3 for n=6 tournaments."""
    n = 6
    tol = 1e-8
    all_paths = {}
    for p in range(0, 6):
        paths = []
        if p + 1 <= n:
            for verts in combinations(range(n), p+1):
                for perm in permutations(verts):
                    if is_allowed_path(T, perm):
                        paths.append(perm)
        all_paths[p] = paths

    def build_omega(pd, mp):
        omega = {}
        for p in range(0, mp + 1):
            a_p = pd[p]
            if not a_p:
                omega[p] = np.zeros((0, 0))
                continue
            a_pm1_set = set(pd[p-1]) if p > 0 else set()
            na = {}
            for sigma in a_p:
                for i in range(1, len(sigma)-1):
                    face = sigma[:i] + sigma[i+1:]
                    if p > 0 and face not in a_pm1_set:
                        if face not in na:
                            na[face] = len(na)
            if not na:
                omega[p] = np.eye(len(a_p))
            else:
                mat = np.zeros((len(na), len(a_p)))
                for j, sigma in enumerate(a_p):
                    for i in range(1, len(sigma)-1):
                        face = sigma[:i] + sigma[i+1:]
                        if face in na:
                            mat[na[face], j] += (-1)**i
                U, S, Vt = np.linalg.svd(mat, full_matrices=True)
                rank = int(np.sum(S > tol))
                if len(a_p) - rank == 0:
                    omega[p] = np.zeros((len(a_p), 0))
                else:
                    omega[p] = Vt[rank:].T
        return omega

    omega = build_omega(all_paths, 4)
    boundary = {}
    for p in range(1, 5):
        a_p = all_paths[p]
        a_pm1 = all_paths[p-1]
        if not a_p or not a_pm1:
            boundary[p] = np.zeros((0, 0))
            continue
        idx = {path: i for i, path in enumerate(a_pm1)}
        mat = np.zeros((len(a_pm1), len(a_p)))
        for j, sigma in enumerate(a_p):
            for i in range(len(sigma)):
                face = sigma[:i] + sigma[i+1:]
                if face in idx:
                    mat[idx[face], j] += (-1)**i
        boundary[p] = mat

    ranks = {}
    dims = {}
    for p in range(5):
        dims[p] = omega[p].shape[1] if omega[p].ndim == 2 else 0
    for p in range(1, 5):
        Om_p = omega[p]
        Om_pm1 = omega[p-1]
        if Om_p.ndim < 2 or Om_p.shape[1] == 0 or Om_pm1.ndim < 2 or Om_pm1.shape[1] == 0:
            ranks[p] = 0
            continue
        dp = Om_pm1.T @ boundary[p] @ Om_p
        S = np.linalg.svd(dp, compute_uv=False)
        ranks[p] = int(np.sum(S > tol))

    ker3 = dims[3] - ranks.get(3, 0)
    im4 = ranks.get(4, 0)
    return ker3 - im4

def main():
    n = 7
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs  # 2,097,152

    beta3_1_scores = {(1,1,1,4,4,4), (2,2,2,3,3,3)}

    print("=" * 70)
    print("EXHAUSTIVE GOOD VERTEX CHECK at n=7")
    print(f"Total tournaments: {n_total:,}")
    print("=" * 70)

    # Step 1: Build a lookup table for 6-vertex tournaments
    # Map: bits -> beta_3 (for n=6 tournaments)
    # Actually 2^15 = 32768, feasible to precompute
    print("\n  Step 1: Precomputing beta_3 for all n=6 tournaments...")
    t0 = time.time()
    n6_arcs = 6*(6-1)//2  # 15
    n6_total = 1 << n6_arcs  # 32768
    beta3_lookup = {}
    count_b3_1 = 0

    for bits6 in range(n6_total):
        if bits6 % 5000 == 0:
            print(f"    ... {bits6}/{n6_total}", flush=True)
        T6 = tournament_from_bits(6, bits6)
        try:
            b3 = compute_beta3_n6(T6)
        except:
            b3 = -1
        beta3_lookup[bits6] = b3
        if b3 == 1:
            count_b3_1 += 1

    t1 = time.time()
    print(f"  Done in {t1-t0:.1f}s. beta_3=1 count: {count_b3_1}")

    # Step 2: For each n=7 tournament, check deletions
    # For each vertex v in {0,...,6}, deleting v gives a 6-vertex tournament.
    # We need to convert the 7-vertex bit encoding to a 6-vertex bit encoding.

    # Bit encoding: for n=7, the edges are ordered as:
    # (0,1), (0,2), (0,3), (0,4), (0,5), (0,6),
    # (1,2), (1,3), (1,4), (1,5), (1,6),
    # (2,3), (2,4), (2,5), (2,6),
    # (3,4), (3,5), (3,6),
    # (4,5), (4,6),
    # (5,6)
    # Total: 21 edges

    # When deleting vertex v, we need the bits for edges NOT involving v.
    # And remap the vertices {0,...,6}\{v} to {0,...,5}.

    # Precompute: for each v, which bit positions in the n=7 encoding
    # correspond to edges NOT involving v, and their mapping to the n=6 encoding.

    # Edge index for n=7
    edge_idx_7 = {}
    idx = 0
    for i in range(7):
        for j in range(i+1, 7):
            edge_idx_7[(i,j)] = idx
            idx += 1

    # For each v, compute the mapping
    deletion_maps = {}  # v -> list of (bit_pos_7, bit_pos_6)
    for v in range(7):
        remaining = [i for i in range(7) if i != v]
        remap = {remaining[k]: k for k in range(6)}
        mapping = []
        idx6 = 0
        for i in range(6):
            for j in range(i+1, 6):
                orig_i = remaining[i]
                orig_j = remaining[j]
                if orig_i > orig_j:
                    orig_i, orig_j = orig_j, orig_i
                bit7 = edge_idx_7[(orig_i, orig_j)]
                mapping.append((bit7, idx6))
                idx6 += 1
        deletion_maps[v] = mapping

    print(f"\n  Step 2: Checking all {n_total:,} tournaments for good vertex...")
    t2 = time.time()

    violations = 0
    checked = 0
    score_filter_successes = 0

    for bits7 in range(n_total):
        if bits7 % 200000 == 0:
            elapsed = time.time() - t2
            rate = bits7 / elapsed if elapsed > 0 else 0
            eta = (n_total - bits7) / rate if rate > 0 else 0
            print(f"    ... {bits7:,}/{n_total:,} ({elapsed:.0f}s, {rate:.0f}/s, ETA {eta:.0f}s)", flush=True)

        # For each deletion, first do a quick score check
        # If any deletion has score NOT in beta3_1_scores, it's automatically good
        found_good = False

        for v in range(7):
            # Extract the 6-vertex tournament bits
            bits6 = 0
            for bit7, bit6 in deletion_maps[v]:
                if bits7 & (1 << bit7):
                    bits6 |= (1 << bit6)

            # Score check first
            T6 = tournament_from_bits(6, bits6)
            scores = tuple(sorted([sum(T6[i][j] for j in range(6) if j != i) for i in range(6)]))
            if scores not in beta3_1_scores:
                found_good = True
                score_filter_successes += 1
                break

            # Score is compatible with beta_3=1, check actual beta_3
            b3 = beta3_lookup.get(bits6, -1)
            if b3 == 0:
                found_good = True
                break

        checked += 1
        if not found_good:
            violations += 1
            T7 = tournament_from_bits(7, bits7)
            scores7 = tuple(sorted([sum(T7[i][j] for j in range(7) if j != i) for i in range(7)]))
            print(f"    !!! VIOLATION at bits={bits7}, scores={scores7}")
            # Show deletion beta_3 values
            for v2 in range(7):
                bits6_2 = 0
                for bit7, bit6 in deletion_maps[v2]:
                    if bits7 & (1 << bit7):
                        bits6_2 |= (1 << bit6)
                b3_2 = beta3_lookup.get(bits6_2, -1)
                print(f"      T\\{v2}: beta_3={b3_2}")

    t3 = time.time()
    print(f"\n  Done in {t3-t2:.1f}s")
    print(f"  Checked: {checked:,}")
    print(f"  Score filter short-circuits: {score_filter_successes:,}")
    print(f"  VIOLATIONS (no good vertex): {violations}")

    if violations == 0:
        print(f"\n  *** THEOREM CONFIRMED EXHAUSTIVELY AT n=7 ***")
        print(f"  Every tournament on 7 vertices has a vertex v with beta_3(T\\v) = 0.")
        print(f"  Combined with H_3(T,T\\v) <= 1 (THM-111), this gives beta_3(T) <= 1.")
    else:
        print(f"\n  !!! {violations} VIOLATIONS FOUND !!!")
        print(f"  The good vertex property FAILS at n=7.")

if __name__ == '__main__':
    main()
