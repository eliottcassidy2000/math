#!/usr/bin/env python3
"""
aw_precision_s20ew.py — Precisely what A and W encode, and joint tricks
kind-pasteur-2026-03-24-S20ew

A = unweighted adjacency (0/1 matrix, no self-loops)
W = wiggly-weighted adjacency (integer matrix, with self-loops on diagonal)

PRECISE ENCODING:
  A[i,j] = 1 iff classes i,j connected by at least one tile flip
  W[i,j] = #{(tiling T in class i, tile X) : T flip X is in class j}
  W[i,i] = #{(tiling T in class i, tile X) : T flip X is still in class i}

IDENTITIES:
  sum_j W[i,j] = tilings(i) * m    (each tiling has m tiles)
  A = sign(W_offdiag)               (same edge set)
  W symmetric                       (each line counted from both sides)

QUESTION: What does W know that A doesn't?
  W encodes tilings(i) via row sums.
  Does W also encode |Aut(i)| or H(i)?
  Is W[i,j]/(tilings(i)*tilings(j)) a function of the EDGE only?

TRICKS:
  1. Normalized matrix N = D_t^{-1/2} W D_t^{-1/2}: removes tiling variation
  2. Ratio R = W[i,j]/A[i,j]: weight per edge, is this a simple function?
  3. P = D_s^{-1} W: Markov matrix, stationary = tilings/sum
"""

import sys
import numpy as np
from math import factorial
from itertools import permutations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  A vs W: PRECISE INFORMATION CONTENT")
print("  kind-pasteur-2026-03-24-S20ew")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in [5, 6]:
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

    def canonicalize(A):
        best = None
        for p in all_perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
            if best is None or s < best: best = s
        return best

    canon_map = {}
    for mask in range(1 << num_tiles):
        bits = [(mask >> k) & 1 for k in range(num_tiles)]
        adj = bits_to_adj(bits)
        canon_map[mask] = canonicalize(adj)

    comp_map = {}
    for cn in set(canon_map.values()):
        for mask, c in canon_map.items():
            if c == cn:
                comp_map[cn] = canon_map[mask ^ ((1 << num_tiles) - 1)]
                break

    merged_map = {cn: min(cn, comp_map.get(cn, cn)) for cn in set(canon_map.values())}
    merged_list = sorted(set(merged_map.values()))
    V = len(merged_list)
    idx = {mcn: i for i, mcn in enumerate(merged_list)}

    # Class properties
    tilings = np.zeros(V)
    for mask in range(1 << num_tiles):
        mcn = merged_map[canon_map[mask]]
        tilings[idx[mcn]] += 1

    # Build W
    W = np.zeros((V, V))
    for mask in range(1 << num_tiles):
        i = idx[merged_map[canon_map[mask]]]
        for wi in range(num_tiles):
            j = idx[merged_map[canon_map[mask ^ (1 << wi)]]]
            W[i, j] += 1

    # A (unweighted, no self-loops)
    W_off = W.copy()
    np.fill_diagonal(W_off, 0)
    A = (W_off > 0).astype(float)
    E = int(A.sum()) // 2

    m = num_tiles
    sl = np.diag(W)

    print(f"\n{'#'*60}")
    print(f"  n = {n}, V = {V}, E = {E}, m = {m}")
    print(f"{'#'*60}")

    # ================================================================
    # IDENTITY 1: Row sum of W = tilings * m
    # ================================================================
    row_sums = W.sum(axis=1)
    expected = tilings * m
    print(f"\n  IDENTITY 1: W @ 1 = tilings * m")
    print(f"    Max |W@1 - t*m| = {np.max(np.abs(row_sums - expected)):.1f}")
    print(f"    EXACT: {np.allclose(row_sums, expected)}")

    # ================================================================
    # IDENTITY 2: W is symmetric
    # ================================================================
    print(f"\n  IDENTITY 2: W = W^T")
    print(f"    Max |W - W^T| = {np.max(np.abs(W - W.T)):.1f}")
    print(f"    EXACT: {np.allclose(W, W.T)}")

    # ================================================================
    # DECOMPOSITION: W = D_t * P (tiling diagonal times Markov)
    # ================================================================
    D_t = np.diag(tilings)
    D_t_inv = np.diag(1.0 / np.maximum(tilings, 1))
    P = D_t_inv @ W  # P[i,j] = W[i,j] / tilings(i)

    # P row sums = m (each tiling has m transitions)
    P_row_sums = P.sum(axis=1)
    print(f"\n  MARKOV: P = D_t^{{-1}} W, row sums = m = {m}")
    print(f"    Max |P@1 - m| = {np.max(np.abs(P_row_sums - m)):.6f}")

    # P is NOT symmetric in general (W is, but D_t * P is not)
    # But detailed balance: tilings(i) * P[i,j] = W[i,j] = W[j,i] = tilings(j) * P[j,i]
    # So P[i,j]/P[j,i] = tilings(j)/tilings(i) — reversibility with stationary = tilings/sum

    # ================================================================
    # KEY QUESTION: What does W know beyond A and tilings?
    # If W = f(A, tilings), then W carries no new info.
    # Test: W[i,j] = c * A[i,j] * sqrt(tilings(i) * tilings(j))?
    # ================================================================
    print(f"\n  WHAT DOES W KNOW BEYOND A AND TILINGS?")

    # Normalized matrix: N = D_t^{-1/2} W D_t^{-1/2}
    D_t_sqrt_inv = np.diag(1.0 / np.sqrt(np.maximum(tilings, 1)))
    N_mat = D_t_sqrt_inv @ W_off @ D_t_sqrt_inv

    # If N is proportional to A, then W = c * D_t^{1/2} A D_t^{1/2}
    # i.e., W[i,j] = c * sqrt(ti * tj) * A[i,j]
    n_vals = []
    for i in range(V):
        for j in range(i+1, V):
            if A[i,j] > 0:
                n_vals.append(N_mat[i,j])

    print(f"    N[i,j] = W[i,j] / sqrt(ti * tj) for edges:")
    print(f"    min = {min(n_vals):.4f}, max = {max(n_vals):.4f}, std/mean = {np.std(n_vals)/np.mean(n_vals):.4f}")
    print(f"    N constant? {np.std(n_vals)/np.mean(n_vals) < 0.01}")

    # Alternative: W[i,j] = c * A[i,j] * tilings(i) * tilings(j) / something?
    edge_ratio_product = []
    edge_ratio_sum = []
    edge_ratio_min = []
    for i in range(V):
        for j in range(i+1, V):
            if A[i,j] > 0:
                ti, tj = tilings[i], tilings[j]
                w = W_off[i,j]
                edge_ratio_product.append(w / (ti * tj))
                edge_ratio_sum.append(w / (ti + tj))
                edge_ratio_min.append(w / min(ti, tj))

    print(f"\n    W / (ti * tj): min={min(edge_ratio_product):.6f}, max={max(edge_ratio_product):.6f}, CV={np.std(edge_ratio_product)/np.mean(edge_ratio_product):.3f}")
    print(f"    W / (ti + tj): min={min(edge_ratio_sum):.4f}, max={max(edge_ratio_sum):.4f}, CV={np.std(edge_ratio_sum)/np.mean(edge_ratio_sum):.3f}")
    print(f"    W / min(ti,tj): min={min(edge_ratio_min):.4f}, max={max(edge_ratio_min):.4f}, CV={np.std(edge_ratio_min)/np.mean(edge_ratio_min):.3f}")

    # ================================================================
    # SPECTRAL TRICK: Compare eigenvalues of N, A, and P
    # ================================================================
    eig_A = np.sort(np.linalg.eigvalsh(A))[::-1]
    eig_N = np.sort(np.linalg.eigvalsh(N_mat))[::-1]
    eig_W = np.sort(np.linalg.eigvalsh(W_off))[::-1]

    # If N ~ c*A then eig_N ~ c*eig_A
    if V > 2:
        slope_NA = np.polyfit(eig_A, eig_N, 1)
        resid_NA = eig_N - np.polyval(slope_NA, eig_A)
        r2_NA = 1 - np.var(resid_NA) / np.var(eig_N)

        print(f"\n  SPECTRAL TRICK: N vs A eigenvalues")
        print(f"    Linear fit: eig_N = {slope_NA[0]:.4f} * eig_A + {slope_NA[1]:.4f}")
        print(f"    R^2 = {r2_NA:.6f}")
        print(f"    Max residual = {np.max(np.abs(resid_NA)):.4f}")
        print(f"    N ~ constant * A? {'YES' if r2_NA > 0.999 else 'NO'}")

    # ================================================================
    # THE TWO COMPUTATIONS TRICK
    # ================================================================
    print(f"\n  JOINT A+W TRICKS:")

    # Trick 1: tilings from W
    tilings_from_W = W.sum(axis=1) / m
    print(f"    1. tilings = W@1 / m (EXACT: {np.allclose(tilings_from_W, tilings)})")

    # Trick 2: degree from A
    degree = A.sum(axis=1)
    print(f"    2. degree = A@1 (trivial)")

    # Trick 3: Average weight per edge = W_off@1 / degree
    avg_weight = np.zeros(V)
    for i in range(V):
        if degree[i] > 0:
            avg_weight[i] = W_off[i,:].sum() / degree[i]
    print(f"    3. avg_weight = W_off_row / degree: min={avg_weight.min():.1f}, max={avg_weight.max():.1f}")

    # Trick 4: Self-loop rate from both
    sl_rate = sl / (tilings * m)  # = P[i,i] / m
    print(f"    4. SL_rate = diag(W) / (tilings*m): min={sl_rate.min():.4f}, max={sl_rate.max():.4f}")

    # Trick 5: "Excess weight" = how much heavier than uniform
    # If all edges had equal weight: W_off[i,j] = (tilings(i)*m - sl(i)) / degree(i)
    # The excess is W_off[i,j] / expected[i,j]
    print(f"    5. Weight excess (W actual / W if uniform):")
    excesses = []
    for i in range(V):
        for j in range(i+1, V):
            if A[i,j] > 0:
                expected_i = (tilings[i] * m - sl[i]) / degree[i]
                excess = W_off[i,j] / expected_i if expected_i > 0 else 0
                excesses.append(excess)
    print(f"       min={min(excesses):.4f}, max={max(excesses):.4f}, mean={np.mean(excesses):.4f}")

    # Trick 6: MUTUAL INFORMATION between A-structure and W-structure
    # The "structural" part of W (beyond tiling count) is N = D^{-1/2} W D^{-1/2}
    # Entropy of N eigenvalues (normalized)
    eig_N_pos = eig_N[eig_N > 1e-10]
    eig_N_norm = eig_N_pos / eig_N_pos.sum()
    entropy_N = -np.sum(eig_N_norm * np.log2(eig_N_norm))
    eig_A_pos = eig_A[eig_A > 1e-10]
    eig_A_norm = eig_A_pos / eig_A_pos.sum()
    entropy_A = -np.sum(eig_A_norm * np.log2(eig_A_norm))
    print(f"    6. Spectral entropy: A={entropy_A:.4f}, N={entropy_N:.4f}, ratio={entropy_N/entropy_A:.4f}")

    # ================================================================
    # RESIDUAL: What N has that A doesn't
    # ================================================================
    if V > 2 and r2_NA < 0.999:
        resid_matrix = N_mat - slope_NA[0] * A - slope_NA[1] * np.eye(V)
        eig_resid = np.sort(np.linalg.eigvalsh(resid_matrix))[::-1]
        print(f"\n  RESIDUAL (N - c*A - d*I):")
        print(f"    Top eigenvalues: {[f'{x:.4f}' for x in eig_resid[:5]]}")
        print(f"    Frobenius norm: {np.linalg.norm(resid_matrix, 'fro'):.4f}")
        print(f"    % of N norm: {np.linalg.norm(resid_matrix, 'fro')/np.linalg.norm(N_mat, 'fro')*100:.2f}%")

        # What does the residual correlate with?
        resid_diag = np.diag(resid_matrix)
        corr_tilings = np.corrcoef(resid_diag, tilings)[0,1]
        corr_sl = np.corrcoef(resid_diag, sl)[0,1]
        print(f"    Residual diagonal vs tilings: r={corr_tilings:.4f}")
        print(f"    Residual diagonal vs self-loops: r={corr_sl:.4f}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\n" + "=" * 70)
print("SUMMARY: A vs W INFORMATION CONTENT")
print("=" * 70)
print("""
A KNOWS:            W KNOWS:
  topology            topology (same edges)
  degree              degree (same)
  spectral gap        spectral gap (different value!)
  clustering          clustering (same topology)
  -                   tiling counts (from row sums)
  -                   self-loop structure
  -                   edge weights (thickness)

W = D_t^{1/2} @ N @ D_t^{1/2}  where N is the "structural" part.
If N ~ c * A: then W just scales A by tiling counts. No new topology.
If N != c * A: then W carries genuine new topological info beyond tilings.

THE TRICK: Compute A for topology, W for tilings and weights.
  W@1/m = tilings (EXACT). No enumeration needed if W is computable.
  A@1 = degree. Standard graph metric.
  diag(W)/(tilings*m) = SL rate per class. = D(n) contribution.
  N = D_t^{-1/2} W D_t^{-1/2} = "pure structure" after removing tilings.
  Compare eig(N) vs eig(A) to see if W adds topology.
""")

print("DONE.")
print("=" * 80)
