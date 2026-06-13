#!/usr/bin/env python3
"""
krawtchouk_h_n7_s306.py -- opus-2026-03-24-S306

KRAWTCHOUK POLYNOMIALS vs HAMILTONIAN PATH COUNT

The tiling space for n-tournaments is Q_m, the m-dimensional binary hypercube,
where m = C(n-1, 2). Each tiling is a binary string of length m. The Hamming
weight of a tiling is the number of 1-bits.

The Krawtchouk polynomials K_k(x; m) are the eigenvectors of the Hamming
association scheme on Q_m. They form an orthogonal basis for functions on
the hypercube.

Prior finding: K_1(x; m) = m - 2x correlates with H(T) at r = -0.94 (n=5)
and r = -0.88 (n=6).

This script computes:
1. For n=5,6,7: build all tilings, iso classes, H values
2. For each tiling: Hamming weight, K_1, K_2, K_3 values
3. Correlation of each K_k with H across all tilings
4. Also: correlation at the iso-class level (weighted by class size)
5. Regression: how much of H variance is explained by K_1, K_2, K_3?

Krawtchouk polynomials (binary Hamming scheme):
  K_k(x; m) = sum_{j=0}^{k} (-1)^j C(x,j) C(m-x, k-j)

First few:
  K_0(x; m) = 1
  K_1(x; m) = m - 2x
  K_2(x; m) = C(m,2) - 2(m-1)x + 2C(x,2)   [= 2x^2 - 2mx + m(m-1)/2]
  K_3(x; m) = C(m,3) - 3C(m-1,2)x + 3(m-2)C(x,2) - 2*2*C(x,3)
             [three-term recurrence is simpler]

Three-term recurrence: K_{k+1}(x) = (m - 2x - 2k) K_k(x) / (k+1) - (m - k + 1) K_{k-1}(x) / (k+1)
Wait -- standard recurrence:
  (k+1) K_{k+1}(x; m) = (m - 2x) K_k(x; m) - (m - k) K_{k-1}(x; m)

We use the explicit sum formula for correctness.

Author: opus-2026-03-24-S306
"""

import sys
import time
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter

sys.stdout.reconfigure(line_buffering=True)


def banner(title):
    print("\n" + "=" * 72)
    print("  " + title)
    print("=" * 72 + "\n")


# --- Krawtchouk polynomial ---
def krawtchouk(k, x, m):
    """Compute K_k(x; m) = sum_{j=0}^{k} (-1)^j C(x, j) C(m-x, k-j)."""
    val = 0
    for j in range(k + 1):
        val += ((-1) ** j) * comb(x, j) * comb(m - x, k - j)
    return val


# --- Canon (isomorphism class representative) ---
def canon(A, n):
    """Canonical form via score-group permutations."""
    sc = [sum(A[i]) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n):
        sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(sg.keys())]

    best = None

    def gp(gs):
        if not gs:
            yield []
            return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]):
                yield list(p) + r

    for p in gp(gs):
        f = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or f < best:
            best = f
    return best


# --- Hamiltonian path count via bitmask DP ---
def Hcount(A, n):
    """Count Hamiltonian paths in tournament with adjacency matrix A."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    full = (1 << n) - 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            c = dp.get((S, v), 0)
            if c == 0:
                continue
            for w in range(n):
                if S & (1 << w):
                    continue
                if A[v][w]:
                    dp[(S | (1 << w), w)] = dp.get((S | (1 << w), w), 0) + c
    return sum(dp.get((full, v), 0) for v in range(n))


# --- Build adjacency matrix from tiling bits ---
def bits_to_adj(bits, n, tile_pairs):
    """Convert tiling bits to adjacency matrix."""
    A = [[0] * n for _ in range(n)]
    # Base path: i -> i-1 for consecutive vertices (descending path n->n-1->...->1)
    # In 0-indexed: i+1 -> i, so A[i+1][i] = 1? No -- looking at existing code:
    # A[i][i-1] = 1 for i in range(1,n), meaning edge from i to i-1 (descending path)
    for i in range(1, n):
        A[i][i - 1] = 1
    for k, (lo, hi) in enumerate(tile_pairs):
        if (bits >> k) & 1:
            A[lo][hi] = 1
        else:
            A[hi][lo] = 1
    return A


def popcount(x):
    """Count number of set bits."""
    return bin(x).count('1')


# ============================================================
banner("KRAWTCHOUK POLYNOMIALS vs H(T) -- TILING HYPERCUBE ANALYSIS")
# ============================================================

print("Krawtchouk polynomials K_k(x; m) on the binary Hamming scheme Q_m")
print("K_k(x; m) = sum_{j=0}^{k} (-1)^j C(x,j) C(m-x, k-j)")
print()
print("Tiling: m = C(n-1, 2) bits, each bit = one non-base-path arc direction")
print("Hamming weight x = number of 1-bits in the tiling")
print()

# Quick sanity check on Krawtchouk
print("Sanity check: K_1(x; 6) for x=0..6:")
for x in range(7):
    print("  K_1(%d; 6) = %d  [expect 6 - 2*%d = %d]" %
          (x, krawtchouk(1, x, 6), x, 6 - 2 * x))
print()

for n in [5, 6, 7]:
    m = comb(n - 1, 2)
    total = 2 ** m
    banner("n = %d, m = %d, total tilings = %d" % (n, m, total))

    t0 = time.time()

    # Build tile pairs (non-base-path arcs)
    base_set = set((i, i + 1) for i in range(n - 1))
    all_p = [(i, j) for i in range(n) for j in range(i + 1, n)]
    tile_pairs = [(i, j) for i, j in all_p if (i, j) not in base_set]
    tile_pairs.sort()

    print("Tile pairs (%d): %s" % (len(tile_pairs), tile_pairs))
    print()

    # --- Phase 1: Compute iso classes for all tilings ---
    print("Phase 1: Computing iso classes for %d tilings..." % total)
    sys.stdout.flush()

    cmap = {}       # canon -> class_id
    tc = {}         # bits -> class_id
    members = defaultdict(list)  # class_id -> [bits]
    sizes = {}      # class_id -> count

    checkpoint = max(1, total // 10)
    for bits in range(total):
        if bits > 0 and bits % checkpoint == 0:
            elapsed = time.time() - t0
            print("  %d/%d tilings processed (%.1f s)" % (bits, total, elapsed))
            sys.stdout.flush()

        A = bits_to_adj(bits, n, tile_pairs)
        c = canon(A, n)
        if c not in cmap:
            cmap[c] = len(cmap)
            sizes[cmap[c]] = 0
        cid = cmap[c]
        sizes[cid] += 1
        tc[bits] = cid
        members[cid].append(bits)

    nc = len(cmap)
    elapsed1 = time.time() - t0
    print("  Done: %d iso classes found in %.1f s" % (nc, elapsed1))
    sys.stdout.flush()

    # --- Phase 2: Compute H for each iso class representative ---
    print("\nPhase 2: Computing H for %d iso class representatives..." % nc)
    sys.stdout.flush()

    H_class = {}  # class_id -> H value
    t1 = time.time()
    for cid in range(nc):
        b = members[cid][0]
        A = bits_to_adj(b, n, tile_pairs)
        H_class[cid] = Hcount(A, n)
    elapsed2 = time.time() - t1
    print("  Done in %.1f s" % elapsed2)
    print("  H range: [%d, %d]" % (min(H_class.values()), max(H_class.values())))
    print("  Distinct H values: %d" % len(set(H_class.values())))

    # --- Phase 3: For each tiling, compute Hamming weight, K_1, K_2, K_3, H ---
    print("\nPhase 3: Computing Krawtchouk values for all tilings...")
    sys.stdout.flush()

    # Arrays for all tilings
    weights = np.zeros(total, dtype=int)
    H_all = np.zeros(total, dtype=int)
    K1_all = np.zeros(total, dtype=float)
    K2_all = np.zeros(total, dtype=float)
    K3_all = np.zeros(total, dtype=float)

    # Precompute Krawtchouk values for each possible weight
    K1_table = {}
    K2_table = {}
    K3_table = {}
    for x in range(m + 1):
        K1_table[x] = krawtchouk(1, x, m)
        K2_table[x] = krawtchouk(2, x, m)
        K3_table[x] = krawtchouk(3, x, m)

    for bits in range(total):
        w = popcount(bits)
        weights[bits] = w
        H_all[bits] = H_class[tc[bits]]
        K1_all[bits] = K1_table[w]
        K2_all[bits] = K2_table[w]
        K3_all[bits] = K3_table[w]

    # --- Phase 4: Correlations across all tilings ---
    print("\nPhase 4: Correlations across ALL %d tilings" % total)
    print("-" * 50)

    # Pearson correlations
    corr_K1_H = np.corrcoef(K1_all, H_all)[0, 1]
    corr_K2_H = np.corrcoef(K2_all, H_all)[0, 1]
    corr_K3_H = np.corrcoef(K3_all, H_all)[0, 1]
    corr_w_H = np.corrcoef(weights.astype(float), H_all)[0, 1]

    print("  corr(K_1, H) = %.6f" % corr_K1_H)
    print("  corr(K_2, H) = %.6f" % corr_K2_H)
    print("  corr(K_3, H) = %.6f" % corr_K3_H)
    print("  corr(weight, H) = %.6f" % corr_w_H)
    print()
    print("  Note: K_1 = m - 2*weight, so corr(K_1, H) = -corr(weight, H)")
    print("  Verified: %.6f vs %.6f" % (corr_K1_H, -corr_w_H))

    # --- Phase 5: Multiple regression K_1, K_2, K_3 -> H ---
    print("\nPhase 5: Multiple regression K_1, K_2, K_3 -> H")
    print("-" * 50)

    # Design matrix: [1, K1, K2, K3]
    X = np.column_stack([np.ones(total), K1_all, K2_all, K3_all])
    y = H_all.astype(float)

    # Least squares
    beta, residuals, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    y_pred = X @ beta
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    print("  Regression: H ~ beta_0 + beta_1*K_1 + beta_2*K_2 + beta_3*K_3")
    print("  beta_0 (intercept) = %.6f" % beta[0])
    print("  beta_1 (K_1)       = %.6f" % beta[1])
    print("  beta_2 (K_2)       = %.6f" % beta[2])
    print("  beta_3 (K_3)       = %.6f" % beta[3])
    print("  R^2 = %.6f" % R2)
    print("  Residual std = %.4f" % np.sqrt(ss_res / total))
    print("  Mean H = %.4f, Std H = %.4f" % (np.mean(y), np.std(y)))

    # Also try K1 alone
    X1 = np.column_stack([np.ones(total), K1_all])
    beta1, _, _, _ = np.linalg.lstsq(X1, y, rcond=None)
    y_pred1 = X1 @ beta1
    ss_res1 = np.sum((y - y_pred1) ** 2)
    R2_1 = 1 - ss_res1 / ss_tot if ss_tot > 0 else 0
    print("\n  K_1-only: R^2 = %.6f" % R2_1)

    # K1 + K2
    X12 = np.column_stack([np.ones(total), K1_all, K2_all])
    beta12, _, _, _ = np.linalg.lstsq(X12, y, rcond=None)
    y_pred12 = X12 @ beta12
    ss_res12 = np.sum((y - y_pred12) ** 2)
    R2_12 = 1 - ss_res12 / ss_tot if ss_tot > 0 else 0
    print("  K_1+K_2: R^2 = %.6f" % R2_12)

    # --- Phase 6: Iso-class level analysis ---
    print("\nPhase 6: Iso-class level correlations")
    print("-" * 50)

    # For each iso class, compute mean Hamming weight of its tilings
    class_mean_w = {}
    class_std_w = {}
    for cid in range(nc):
        ws = [popcount(b) for b in members[cid]]
        class_mean_w[cid] = np.mean(ws)
        class_std_w[cid] = np.std(ws)

    # Arrays at class level
    H_cls = np.array([H_class[cid] for cid in range(nc)])
    w_cls = np.array([class_mean_w[cid] for cid in range(nc)])
    K1_cls = np.array([m - 2 * class_mean_w[cid] for cid in range(nc)])
    sz_cls = np.array([sizes[cid] for cid in range(nc)])

    # Unweighted correlation at class level
    corr_cls = np.corrcoef(K1_cls, H_cls)[0, 1]
    print("  Unweighted corr(K_1_mean, H) over %d classes = %.6f" % (nc, corr_cls))

    # Weighted correlation (by class size)
    def weighted_corr(x, y, w):
        sw = np.sum(w)
        mx = np.sum(w * x) / sw
        my = np.sum(w * y) / sw
        cov = np.sum(w * (x - mx) * (y - my)) / sw
        sx = np.sqrt(np.sum(w * (x - mx) ** 2) / sw)
        sy = np.sqrt(np.sum(w * (y - my) ** 2) / sw)
        return cov / (sx * sy) if sx > 0 and sy > 0 else 0

    corr_wcls = weighted_corr(K1_cls, H_cls, sz_cls)
    print("  Size-weighted corr(K_1_mean, H) = %.6f" % corr_wcls)

    # --- Phase 7: Weight distribution per H level ---
    print("\nPhase 7: Weight distribution by H level")
    print("-" * 50)

    H_levels = sorted(set(H_class.values()))
    print("  %-8s  %-8s  %-10s  %-10s  %-8s" %
          ("H", "#class", "mean_wt", "std_wt", "#tilings"))
    for h in H_levels:
        cids = [cid for cid in range(nc) if H_class[cid] == h]
        all_weights = []
        total_t = 0
        for cid in cids:
            for b in members[cid]:
                all_weights.append(popcount(b))
            total_t += sizes[cid]
        mw = np.mean(all_weights) if all_weights else 0
        sw = np.std(all_weights) if all_weights else 0
        print("  %-8d  %-8d  %-10.3f  %-10.3f  %-8d" %
              (h, len(cids), mw, sw, total_t))

    # --- Phase 8: Higher Krawtchouk at class level ---
    print("\nPhase 8: Higher Krawtchouk at iso-class level")
    print("-" * 50)

    # For each class, compute mean K_k over all member tilings
    for k in range(1, min(7, m + 1)):
        Kk_table = {x: krawtchouk(k, x, m) for x in range(m + 1)}
        Kk_cls = np.array([
            np.mean([Kk_table[popcount(b)] for b in members[cid]])
            for cid in range(nc)
        ])
        corr_k = np.corrcoef(Kk_cls, H_cls)[0, 1]
        corr_k_w = weighted_corr(Kk_cls, H_cls, sz_cls)
        print("  K_%d: unweighted corr = %+.6f, weighted corr = %+.6f" %
              (k, corr_k, corr_k_w))

    # --- Phase 9: Krawtchouk expansion of H ---
    print("\nPhase 9: Krawtchouk expansion coefficients of H")
    print("-" * 50)
    print("  H as function on Q_m. Expand H(bits) in Krawtchouk basis.")
    print("  hat{H}_k = (1/2^m) * sum_x H(x) * K_k(wt(x); m)")
    print("  (This is the k-th spectral coefficient in the Hamming scheme)")

    # Compute spectral coefficients
    max_k = min(7, m)
    for k in range(max_k + 1):
        Kk_table = {x: krawtchouk(k, x, m) for x in range(m + 1)}
        coeff = 0.0
        for bits in range(total):
            coeff += H_all[bits] * Kk_table[weights[bits]]
        coeff /= total  # normalize by 2^m
        # Also compute norm: <K_k, K_k> = sum K_k(wt)^2 * C(m,wt) / 2^m
        norm_k = sum(Kk_table[x] ** 2 * comb(m, x) for x in range(m + 1)) / total
        hat_norm = coeff / norm_k if norm_k > 0 else 0
        print("  hat{H}_%d = %.6f  (norm K_%d = %.4f, normalized = %.6f)" %
              (k, coeff, k, norm_k, hat_norm))

    # --- Phase 10: Summary statistics ---
    elapsed_total = time.time() - t0
    print("\nTotal time for n=%d: %.1f s" % (n, elapsed_total))

    print("\n" + "=" * 72)
    print("  SUMMARY for n=%d" % n)
    print("=" * 72)
    print("  m = %d, 2^m = %d, iso classes = %d" % (m, total, nc))
    print("  corr(K_1, H) = %.6f  [all tilings]" % corr_K1_H)
    print("  corr(K_2, H) = %.6f  [all tilings]" % corr_K2_H)
    print("  corr(K_3, H) = %.6f  [all tilings]" % corr_K3_H)
    print("  R^2(K_1 only)      = %.6f" % R2_1)
    print("  R^2(K_1 + K_2)     = %.6f" % R2_12)
    print("  R^2(K_1 + K_2 + K_3) = %.6f" % R2)
    print("  corr(K_1_mean, H) at class level = %.6f" % corr_cls)
    print()

banner("COMPUTATION COMPLETE")
