#!/usr/bin/env python3
"""
wiggly_metagraph_deep_s20ev.py — Complete analysis of the wiggly-weighted metagraph
kind-pasteur-2026-03-24-S20ev

The wiggly metagraph W_n has:
  - Same vertices as G_n/Z_2 (merged iso classes)
  - Same edges as G_n/Z_2
  - WEIGHTED edges: W[i,j] = total wiggly lines connecting classes i and j
  - WEIGHTED self-loops: W[i,i] = total wiggly self-loop lines at class i

Compute ALL standard metagraph invariants for both the unweighted A and weighted W:
  1. Degree sequences (unweighted vs weighted)
  2. Eigenvalues of A and W
  3. Laplacian spectra
  4. Spectral gap (Markov chain mixing)
  5. W/A ratio per edge (is W proportional to A?)
  6. Self-loop weight vs degree correlation
  7. H-gradient in weighted graph
  8. Spine/ribs/sea weighted decomposition
  9. Cheeger constant (weighted vs unweighted)
  10. Random walk stationary distribution
"""

import sys
import numpy as np
from math import factorial, comb
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  WIGGLY METAGRAPH: COMPLETE INVARIANT ANALYSIS")
print("  kind-pasteur-2026-03-24-S20ev")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in range(4, 7):
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

    # Build tilings and canonical forms
    canon_map = {}
    for mask in range(1 << num_tiles):
        bits = [(mask >> k) & 1 for k in range(num_tiles)]
        A = bits_to_adj(bits)
        canon_map[mask] = canonicalize(A)

    # Merged classes
    comp_map_cn = {}
    for cn in set(canon_map.values()):
        for mask, c in canon_map.items():
            if c == cn:
                comp_mask = mask ^ ((1 << num_tiles) - 1)
                comp_map_cn[cn] = canon_map[comp_mask]
                break

    merged_map = {cn: min(cn, comp_map_cn.get(cn, cn)) for cn in set(canon_map.values())}
    merged_list = sorted(set(merged_map.values()))
    V = len(merged_list)
    mcn_to_idx = {mcn: i for i, mcn in enumerate(merged_list)}

    # Class properties
    class_tilings = defaultdict(int)
    for mask in range(1 << num_tiles):
        mcn = merged_map[canon_map[mask]]
        class_tilings[mcn] += 1

    # Compute |Aut| and SC status
    class_info = {}
    for mcn in merged_list:
        for mask, cn in canon_map.items():
            if merged_map[cn] == mcn:
                size = sum(1 for m2 in range(1 << num_tiles) if merged_map[canon_map[m2]] == mcn)
                # Check if SC
                is_sc = (comp_map_cn.get(cn, cn) == cn)
                aut = factorial(N) * (1 if is_sc else 2) // size
                tilings = class_tilings[mcn]
                H = tilings * aut // (2 if not is_sc else 1)
                class_info[mcn] = {'sc': is_sc, 'aut': aut, 'tilings': tilings, 'H': H}
                break

    # Build weighted adjacency matrix
    W = np.zeros((V, V))  # wiggly weighted
    A_unw = np.zeros((V, V))  # unweighted (0/1)

    for mask in range(1 << num_tiles):
        mcn1 = merged_map[canon_map[mask]]
        i = mcn_to_idx[mcn1]
        for wi in range(num_tiles):
            flip_mask = mask ^ (1 << wi)
            mcn2 = merged_map[canon_map[flip_mask]]
            j = mcn_to_idx[mcn2]
            W[i, j] += 1

    # Make symmetric (we counted each direction separately)
    W = (W + W.T) / 2

    # Build unweighted
    A_unw = (W > 0).astype(float)
    np.fill_diagonal(A_unw, 0)  # no self-loops in unweighted

    # Self-loop weights (diagonal of W)
    self_loops = np.diag(W).copy()
    W_no_diag = W.copy()
    np.fill_diagonal(W_no_diag, 0)

    print(f"\n{'#'*70}")
    print(f"  n = {n}, V = {V}, tiles = {num_tiles}")
    print(f"{'#'*70}")

    # ================================================================
    # 1. DEGREE SEQUENCES
    # ================================================================
    deg_unw = A_unw.sum(axis=1)
    deg_w = W_no_diag.sum(axis=1)
    strength = W.sum(axis=1)  # total weight including self-loops

    print(f"\n  1. DEGREE SEQUENCES:")
    print(f"    Unweighted: min={deg_unw.min():.0f}, max={deg_unw.max():.0f}, avg={deg_unw.mean():.1f}")
    print(f"    Weighted (no SL): min={deg_w.min():.0f}, max={deg_w.max():.0f}, avg={deg_w.mean():.1f}")
    print(f"    Strength (with SL): min={strength.min():.0f}, max={strength.max():.0f}, avg={strength.mean():.1f}")
    print(f"    Self-loop: min={self_loops.min():.0f}, max={self_loops.max():.0f}, avg={self_loops.mean():.1f}")

    # ================================================================
    # 2. EIGENVALUES
    # ================================================================
    eig_A = sorted(np.linalg.eigvalsh(A_unw), reverse=True)
    eig_W = sorted(np.linalg.eigvalsh(W_no_diag), reverse=True)
    eig_W_full = sorted(np.linalg.eigvalsh(W), reverse=True)

    print(f"\n  2. EIGENVALUES:")
    print(f"    Unweighted A (top 5): {[f'{x:.3f}' for x in eig_A[:5]]}")
    print(f"    Weighted W (top 5):   {[f'{x:.3f}' for x in eig_W[:5]]}")
    print(f"    W+SL (top 5):         {[f'{x:.3f}' for x in eig_W_full[:5]]}")

    # Spectral gap
    gap_A = eig_A[0] - eig_A[1] if V > 1 else 0
    gap_W = eig_W[0] - eig_W[1] if V > 1 else 0
    print(f"    Spectral gap A: {gap_A:.4f}")
    print(f"    Spectral gap W: {gap_W:.4f}")

    # ================================================================
    # 3. LAPLACIAN
    # ================================================================
    D_unw = np.diag(deg_unw)
    L_unw = D_unw - A_unw
    D_w = np.diag(deg_w)
    L_w = D_w - W_no_diag

    eig_L_unw = sorted(np.linalg.eigvalsh(L_unw))
    eig_L_w = sorted(np.linalg.eigvalsh(L_w))

    print(f"\n  3. LAPLACIAN SPECTRUM:")
    print(f"    Unweighted (first 5): {[f'{x:.3f}' for x in eig_L_unw[:5]]}")
    print(f"    Weighted (first 5):   {[f'{x:.3f}' for x in eig_L_w[:5]]}")
    print(f"    Algebraic connectivity A: {eig_L_unw[1]:.4f}")
    print(f"    Algebraic connectivity W: {eig_L_w[1]:.4f}")

    # ================================================================
    # 4. W/A RATIO — is W proportional to A?
    # ================================================================
    ratios = []
    for i in range(V):
        for j in range(i+1, V):
            if A_unw[i,j] > 0:
                ratios.append(W_no_diag[i,j])

    if ratios:
        print(f"\n  4. W/A RATIO (edge weights):")
        print(f"    min={min(ratios):.1f}, max={max(ratios):.1f}, avg={sum(ratios)/len(ratios):.1f}")
        print(f"    W proportional to A? {max(ratios) - min(ratios) < 0.001}")
        # Distribution
        weight_hist = Counter(int(r) for r in ratios)
        print(f"    Weight distribution: {dict(sorted(weight_hist.items()))}")

    # ================================================================
    # 5. SELF-LOOP vs DEGREE CORRELATION
    # ================================================================
    if V > 2:
        sl_arr = self_loops
        dw_arr = deg_w
        corr = np.corrcoef(sl_arr, dw_arr)[0,1]
        print(f"\n  5. SELF-LOOP vs DEGREE CORRELATION: {corr:.4f}")

        # Self-loop as fraction of total strength
        sl_frac = sl_arr / np.maximum(strength, 1)
        print(f"    SL fraction: min={sl_frac.min():.4f}, max={sl_frac.max():.4f}, avg={sl_frac.mean():.4f}")

    # ================================================================
    # 6. SPINE/RIBS/SEA WEIGHTED
    # ================================================================
    sc_sc_w = 0; sc_ns_w = 0; ns_ns_w = 0
    sc_sc_c = 0; sc_ns_c = 0; ns_ns_c = 0

    for i in range(V):
        for j in range(i+1, V):
            if A_unw[i,j] > 0:
                mcn_i = merged_list[i]
                mcn_j = merged_list[j]
                sci = class_info[mcn_i]['sc']
                scj = class_info[mcn_j]['sc']
                w = W_no_diag[i,j]
                if sci and scj:
                    sc_sc_w += w; sc_sc_c += 1
                elif sci or scj:
                    sc_ns_w += w; sc_ns_c += 1
                else:
                    ns_ns_w += w; ns_ns_c += 1

    total_w = sc_sc_w + sc_ns_w + ns_ns_w
    total_c = sc_sc_c + sc_ns_c + ns_ns_c

    print(f"\n  6. SPINE/RIBS/SEA (weighted):")
    print(f"    {'Type':>8} {'Count':>6} {'Weight':>10} {'Avg weight':>12} {'% of total':>12}")
    for name, c, w in [('SC-SC', sc_sc_c, sc_sc_w), ('SC-NS', sc_ns_c, sc_ns_w), ('NS-NS', ns_ns_c, ns_ns_w)]:
        avg = w/c if c > 0 else 0
        pct = w/total_w*100 if total_w > 0 else 0
        print(f"    {name:>8} {c:6d} {w:10.0f} {avg:12.1f} {pct:11.1f}%")

    # ================================================================
    # 7. RANDOM WALK (Markov chain)
    # ================================================================
    # Transition matrix: P[i,j] = W[i,j] / strength[i]
    P = W / np.maximum(strength.reshape(-1,1), 1)
    # Stationary distribution (left eigenvector with eigenvalue 1)
    # = proportional to strength
    pi_stat = strength / strength.sum()

    print(f"\n  7. RANDOM WALK:")
    print(f"    Stationary dist (top 5): {sorted(pi_stat, reverse=True)[:5]}")
    print(f"    Entropy of stationary: {-sum(p*np.log2(p) for p in pi_stat if p > 0):.4f} (max={np.log2(V):.4f})")
    # Self-return probability
    self_return = sum(P[i,i] * pi_stat[i] for i in range(V))
    print(f"    Expected self-return: {self_return:.4f}")

    # Markov spectral gap (second eigenvalue of P)
    eig_P = sorted(np.linalg.eigvals(P).real, reverse=True)
    print(f"    Markov eigenvalues (top 5): {[f'{x:.4f}' for x in eig_P[:5]]}")
    if V > 1:
        markov_gap = 1 - eig_P[1]
        print(f"    Markov spectral gap: {markov_gap:.4f}")
        print(f"    Mixing time estimate: {1/markov_gap:.1f} steps")

    # ================================================================
    # 8. EIGENVECTOR CORRELATION WITH H
    # ================================================================
    H_vals = np.array([class_info[merged_list[i]]['H'] for i in range(V)])
    tilings_vals = np.array([class_info[merged_list[i]]['tilings'] for i in range(V)])

    # Eigenvectors of W
    evals_W, evecs_W = np.linalg.eigh(W_no_diag)
    idx = np.argsort(-evals_W)
    evals_W = evals_W[idx]
    evecs_W = evecs_W[:, idx]

    print(f"\n  8. EIGENVECTOR CORRELATION WITH H:")
    for k in range(min(5, V)):
        v = evecs_W[:, k]
        corr_H = abs(np.corrcoef(v, H_vals)[0,1]) if np.std(H_vals) > 0 else 0
        corr_T = abs(np.corrcoef(v, tilings_vals)[0,1]) if np.std(tilings_vals) > 0 else 0
        print(f"    Eigvec {k} (eval={evals_W[k]:.2f}): |corr(H)|={corr_H:.3f}, |corr(tilings)|={corr_T:.3f}")

    # ================================================================
    # 9. WEIGHTED DIAMETER & DISTANCES
    # ================================================================
    # Shortest path using 1/weight as distance
    from scipy.sparse.csgraph import shortest_path
    if np.any(W_no_diag > 0):
        dist_unw = shortest_path(A_unw, directed=False)
        # For weighted: use resistance distance (1/weight)
        W_resist = np.zeros_like(W_no_diag)
        W_resist[W_no_diag > 0] = 1.0 / W_no_diag[W_no_diag > 0]
        W_resist[W_no_diag == 0] = 0  # no edge = infinite distance
        # Use shortest path on resistance
        try:
            dist_w = shortest_path(W_resist, directed=False)
            diam_unw = int(dist_unw[dist_unw < np.inf].max())
            diam_w = dist_w[dist_w < np.inf].max()
            print(f"\n  9. DISTANCES:")
            print(f"    Unweighted diameter: {diam_unw}")
            print(f"    Resistance diameter: {diam_w:.4f}")
            print(f"    Avg unweighted distance: {dist_unw[dist_unw < np.inf].mean():.3f}")
            print(f"    Avg resistance distance: {dist_w[dist_w < np.inf].mean():.6f}")
        except:
            print(f"\n  9. DISTANCES: scipy not available for weighted shortest path")

    # ================================================================
    # 10. KEY COMPARISON: W vs A spectra
    # ================================================================
    print(f"\n  10. SPECTRAL COMPARISON W vs A:")
    # Check if W = c * A + d * I (affine relation)
    if V > 2 and np.std(eig_A[:V]) > 0:
        # Linear regression: eig_W = a * eig_A + b
        eig_A_arr = np.array(eig_A[:V])
        eig_W_arr = np.array(eig_W[:V])
        A_mat = np.vstack([eig_A_arr, np.ones(V)]).T
        result = np.linalg.lstsq(A_mat, eig_W_arr, rcond=None)
        slope, intercept = result[0]
        residual = eig_W_arr - (slope * eig_A_arr + intercept)
        r_squared = 1 - np.var(residual) / np.var(eig_W_arr) if np.var(eig_W_arr) > 0 else 0

        print(f"    Linear fit: eig_W ~ {slope:.2f} * eig_A + {intercept:.2f}")
        print(f"    R^2 = {r_squared:.6f}")
        print(f"    Max residual: {np.max(np.abs(residual)):.4f}")
        print(f"    W ~ {slope:.2f} * A? {'YES (perfect)' if r_squared > 0.9999 else 'NO (structure differs)'}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
