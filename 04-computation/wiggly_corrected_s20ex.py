#!/usr/bin/env python3
"""
wiggly_corrected_s20ex.py — CORRECTED wiggly metagraph (proper complement)
kind-pasteur-2026-03-24-S20ex

FIX: The complement of a tournament reverses ALL arcs (not just tiles).
In the tiling model, this means: build the adjacency, reverse all arcs,
then find the canonical form. NOT just flip tile bits.

This script:
1. Computes the CORRECT G_n/Z_2 from the tiling model
2. Builds the CORRECT W and A matrices
3. Re-runs the key analyses with the proper merging
4. Verifies V_merged matches (A000568 + SC)/2
"""

import sys
import numpy as np
from math import factorial
from itertools import permutations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  CORRECTED WIGGLY METAGRAPH (proper complement)")
print("  kind-pasteur-2026-03-24-S20ex")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in [4, 5, 6]:
    t0 = time.time()
    TILES = get_tiles(n)
    num_tiles = len(TILES)
    N = n
    VERTS = list(range(n, 0, -1))
    all_perms = list(permutations(range(N)))

    def bits_to_adj(bits):
        A = [[0]*N for _ in range(N)]
        for k in range(N-1): A[k][k+1] = 1
        for i in range(num_tiles):
            x, y = TILES[i]
            xi, yi = VERTS.index(x), VERTS.index(y)
            if bits[i] == 0: A[xi][yi] = 1
            else: A[yi][xi] = 1
        return A

    def adj_complement(A):
        """TRUE complement: reverse ALL arcs."""
        n = len(A)
        C = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j:
                    C[i][j] = 1 - A[i][j]
        return C

    def canonicalize(A):
        best = None
        for p in all_perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
            if best is None or s < best: best = s
        return best

    # Build all tilings and their canonical forms
    canon_map = {}  # mask -> canon string
    adj_cache = {}  # mask -> adjacency matrix
    for mask in range(1 << num_tiles):
        bits = [(mask >> k) & 1 for k in range(num_tiles)]
        A = bits_to_adj(bits)
        adj_cache[mask] = A
        canon_map[mask] = canonicalize(A)

    # TRUE complement: reverse ALL arcs in adjacency, then canonicalize
    comp_canon = {}  # canon -> comp_canon
    for cn in set(canon_map.values()):
        # Find a mask with this canon
        for mask, c in canon_map.items():
            if c == cn:
                A = adj_cache[mask]
                A_comp = adj_complement(A)
                comp_canon[cn] = canonicalize(A_comp)
                break

    # Merged classes (true complement merging)
    merged_map = {}
    for cn in set(canon_map.values()):
        merged_map[cn] = min(cn, comp_canon.get(cn, cn))

    merged_list = sorted(set(merged_map.values()))
    V = len(set(canon_map.values()))  # unmerged
    V_merged = len(merged_list)
    SC = sum(1 for cn in set(canon_map.values()) if comp_canon.get(cn, cn) == cn)
    mcn_idx = {mcn: i for i, mcn in enumerate(merged_list)}

    # Verify
    expected_V = {4: 4, 5: 12, 6: 56}
    expected_SC = {4: 2, 5: 8, 6: 12}  # self-complementary iso classes

    print(f"\n{'#'*60}")
    print(f"  n = {n}, tiles = {num_tiles}")
    print(f"  V_unmerged = {V} (expected {expected_V.get(n,'?')})")
    print(f"  SC = {SC} (expected {expected_SC.get(n,'?')})")
    print(f"  V_merged = {V_merged} (expected {(expected_V.get(n,0)+expected_SC.get(n,0))//2})")
    print(f"{'#'*60}")

    # Compute tiling counts per merged class
    tilings = np.zeros(V_merged)
    for mask in range(1 << num_tiles):
        mcn = merged_map[canon_map[mask]]
        tilings[mcn_idx[mcn]] += 1

    # Compute H per merged class (count Hamiltonian paths of representative)
    # H = number of labeled tournaments in class / |class| ... no
    # Actually need to count HP. For small n, use the tiling count identity:
    # #tilings(C) = H(C) if C is an unmerged class (each HP gives one tiling)
    # Wait: #tilings represents # labeled tournaments = n!/|Aut|
    # H = # Hamiltonian paths of any T in C
    # # labeled tournaments in C = n!/|Aut(C)|
    # These are different! H can vary independently of class size.
    # Need actual HP count.
    def count_ham_paths(A):
        """Count Hamiltonian paths in tournament A."""
        n = len(A)
        # DP: dp[mask][v] = # of paths visiting vertices in mask, ending at v
        dp = [[0]*n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)): continue
                if dp[mask][v] == 0: continue
                for u in range(n):
                    if mask & (1 << u): continue
                    if A[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]
        full = (1 << n) - 1
        return sum(dp[full])

    # Compute H for each merged class
    H_vals = {}
    for mcn in merged_list:
        for mask, cn in canon_map.items():
            if merged_map[cn] == mcn:
                A = adj_cache[mask]
                H_vals[mcn] = count_ham_paths(A)
                break

    m = num_tiles

    # Build W (wiggly weight matrix) on CORRECT merged classes
    W = np.zeros((V_merged, V_merged))
    for mask in range(1 << num_tiles):
        i = mcn_idx[merged_map[canon_map[mask]]]
        for wi in range(num_tiles):
            j = mcn_idx[merged_map[canon_map[mask ^ (1 << wi)]]]
            W[i, j] += 1

    W_off = W.copy()
    np.fill_diagonal(W_off, 0)
    A_unw = (W_off > 0).astype(float)
    E = int(A_unw.sum()) // 2
    sl = np.diag(W)

    print(f"\n  CORRECTED METAGRAPH:")
    print(f"    E = {E}")
    print(f"    Max degree = {int(A_unw.sum(axis=1).max())}")

    # Verify W@1 = tilings * m
    row_sums = W.sum(axis=1)
    print(f"    W@1 = tilings*m? {np.allclose(row_sums, tilings * m)}")

    # ================================================================
    # CORRECTED SPECTRAL ANALYSIS
    # ================================================================
    eig_A = np.sort(np.linalg.eigvalsh(A_unw))[::-1]
    eig_W = np.sort(np.linalg.eigvalsh(W_off))[::-1]

    print(f"\n  SPECTRAL (corrected):")
    print(f"    A eigenvalues (top 5): {[f'{x:.3f}' for x in eig_A[:5]]}")
    print(f"    W eigenvalues (top 5): {[f'{x:.3f}' for x in eig_W[:5]]}")

    # Normalized N
    D_t_si = np.diag(1.0 / np.sqrt(np.maximum(tilings, 1)))
    N_mat = D_t_si @ W_off @ D_t_si
    eig_N = np.sort(np.linalg.eigvalsh(N_mat))[::-1]

    if V_merged > 2:
        slope = np.polyfit(eig_A, eig_N, 1)
        resid = eig_N - np.polyval(slope, eig_A)
        r2 = 1 - np.var(resid) / np.var(eig_N)
        frob_resid = np.linalg.norm(N_mat - slope[0]*A_unw - slope[1]*np.eye(V_merged), 'fro')
        frob_N = np.linalg.norm(N_mat, 'fro')

        print(f"    N vs A: slope={slope[0]:.4f}, R^2={r2:.6f}")
        print(f"    Frobenius residual: {frob_resid:.4f} / {frob_N:.4f} = {frob_resid/frob_N*100:.1f}%")

    # ================================================================
    # CORRECTED H CORRELATION
    # ================================================================
    H_arr = np.array([H_vals[merged_list[i]] for i in range(V_merged)])
    t_arr = tilings

    evals, evecs = np.linalg.eigh(W_off)
    order = np.argsort(-evals)
    evals = evals[order]
    evecs = evecs[:, order]

    print(f"\n  H CORRELATION (corrected):")
    for k in range(min(3, V_merged)):
        v = evecs[:, k]
        cH = abs(np.corrcoef(v, H_arr)[0,1]) if np.std(H_arr) > 0 else 0
        cT = abs(np.corrcoef(v, t_arr)[0,1]) if np.std(t_arr) > 0 else 0
        print(f"    Eigvec {k} (eval={evals[k]:.1f}): |corr(H)|={cH:.3f}, |corr(tilings)|={cT:.3f}")

    # ================================================================
    # CORRECTED SELF-LOOP ANALYSIS
    # ================================================================
    sl_rate = sl / np.maximum(tilings * m, 1)
    print(f"\n  SELF-LOOP (corrected):")
    print(f"    SL rate: min={sl_rate.min():.4f}, max={sl_rate.max():.4f}, avg={sl_rate.mean():.4f}")

    # Per-class detail
    print(f"\n  PER-CLASS (sorted by H):")
    print(f"  {'idx':>4} {'H':>6} {'tilings':>8} {'deg':>4} {'W_deg':>8} {'SL':>8} {'SL%':>7}")
    for i in sorted(range(V_merged), key=lambda i: H_arr[i]):
        deg = int(A_unw[i,:].sum())
        w_deg = int(W_off[i,:].sum())
        print(f"  {i:4d} {H_arr[i]:6.0f} {tilings[i]:8.0f} {deg:4d} {w_deg:8.0f} {sl[i]:8.0f} {sl_rate[i]*100:6.2f}%")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
