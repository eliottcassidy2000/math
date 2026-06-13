#!/usr/bin/env python3
"""
effective_dimension_s20fv.py — How 1-dimensional is the metagraph?
kind-pasteur-2026-03-24-S20fv

Tournament space is "effectively 1D" via H. But how 1D exactly?

Measures:
1. H as a coordinate: how much of the graph structure does H capture?
   - Correlation of H with principal eigenvector (already ~0.85-0.95)
   - Fraction of edges that are "uphill" in H (= 0 for DAG)
   - H-level width: how many classes at each H value?

2. Effective dimension from eigenvalue decay:
   - If eigenvalues decay as lambda_k ~ k^{-alpha}, alpha determines dimension
   - For 1D: alpha = 2 (quadratic decay). For 2D: alpha = 1.

3. The "excess dimension": how much structure is PERPENDICULAR to H?
   - Second eigenvector captures orthogonal variation
   - The residual after projecting onto H

4. Connection to sqrt(2): the pseudo-doubling ratio (2n-5)/(n-2) -> 2
   - Does the effective dimension approach log(2)/log(something)?

5. Opposing trends as beta_1/beta_3 analogy:
   - Scheme breaking ~ beta_1 (grows = more "holes" in regularity)
   - Markov accuracy ~ beta_3 (shrinks = better "filling" of structure)
   - Do they satisfy a seesaw? Product -> 0?
"""

import sys
import numpy as np
from math import factorial, comb
from itertools import permutations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  EFFECTIVE DIMENSION OF THE METAGRAPH")
print("  kind-pasteur-2026-03-24-S20fv")
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
    m = len(TILES)
    N = n
    VERTS = list(range(n, 0, -1))
    all_perms = list(permutations(range(N)))
    tv = [(VERTS.index(x), VERTS.index(y)) for x, y in TILES]

    def b2a(bits):
        A = [[0]*N for _ in range(N)]
        for k in range(N-1): A[k][k+1] = 1
        for i in range(m):
            xi, yi = tv[i]
            if bits[i] == 0: A[xi][yi] = 1
            else: A[yi][xi] = 1
        return A

    def canonicalize(A):
        best = None
        for p in all_perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
            if best is None or s < best: best = s
        return best

    def count_hp(A):
        dp = [[0]*N for _ in range(1 << N)]
        for v in range(N): dp[1 << v][v] = 1
        for mask in range(1, 1 << N):
            for v in range(N):
                if not (mask & (1 << v)) or dp[mask][v] == 0: continue
                for u in range(N):
                    if mask & (1 << u): continue
                    if A[v][u]: dp[mask | (1 << u)][u] += dp[mask][v]
        return sum(dp[(1 << N) - 1])

    canon_map = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = b2a(bits)
        canon_map[mask] = canonicalize(A)

    classes = sorted(set(canon_map.values()))
    V = len(classes)
    cidx = {cn: i for i, cn in enumerate(classes)}

    # H values
    H_vals = np.zeros(V)
    for i, cn in enumerate(classes):
        for mask, c in canon_map.items():
            if c == cn:
                H_vals[i] = count_hp(b2a([(mask >> k) & 1 for k in range(m)]))
                break

    # Build wiggly adjacency
    W = np.zeros((V, V))
    for mask in range(1 << m):
        i = cidx[canon_map[mask]]
        for wi in range(m):
            j = cidx[canon_map[mask ^ (1 << wi)]]
            W[i, j] += 1

    W_off = W.copy()
    np.fill_diagonal(W_off, 0)
    A_unw = (W_off > 0).astype(float)

    print(f"\n{'#'*60}")
    print(f"  n = {n}, V = {V}, m = {m}")
    print(f"{'#'*60}")

    # ================================================================
    # 1. H AS A COORDINATE
    # ================================================================
    evals, evecs = np.linalg.eigh(A_unw)
    idx = np.argsort(-evals)
    evals = evals[idx]
    evecs = evecs[:, idx]

    # Normalize H
    H_norm = (H_vals - H_vals.mean()) / np.std(H_vals)

    corrs = [abs(np.corrcoef(evecs[:,k], H_norm)[0,1]) for k in range(min(V, 10))]
    best_k = np.argmax(corrs)

    print(f"\n  1. H AS COORDINATE:")
    print(f"    H range: {H_vals.min():.0f} to {H_vals.max():.0f}")
    print(f"    H levels: {len(set(H_vals))}")
    print(f"    Best eigenvector for H: k={best_k} (corr={corrs[best_k]:.4f})")
    print(f"    Top 5 correlations: {[f'{c:.3f}' for c in corrs[:5]]}")

    # DAG check: all edges go from lower H to higher H?
    uphill = 0; downhill = 0; level = 0
    for i in range(V):
        for j in range(V):
            if A_unw[i,j] > 0 and i != j:
                if H_vals[i] < H_vals[j]: uphill += 1
                elif H_vals[i] > H_vals[j]: downhill += 1
                else: level += 1

    print(f"    Directed edges: uphill={uphill}, downhill={downhill}, level={level}")
    print(f"    DAG under H? {'YES' if downhill == 0 else 'NO'}")

    # ================================================================
    # 2. EIGENVALUE DECAY -> EFFECTIVE DIMENSION
    # ================================================================
    pos_evals = evals[evals > 0.01]
    if len(pos_evals) > 2:
        log_k = np.log(np.arange(1, len(pos_evals)+1))
        log_lam = np.log(pos_evals)
        slope, intercept = np.polyfit(log_k[:min(10,len(log_k))], log_lam[:min(10,len(log_lam))], 1)

        print(f"\n  2. EIGENVALUE DECAY:")
        print(f"    Positive eigenvalues: {len(pos_evals)}")
        print(f"    Decay exponent alpha = {-slope:.3f}")
        print(f"    (1D: alpha~2, 2D: alpha~1, interpretation: dim ~ 2/alpha)")
        print(f"    Effective dimension: {2/(-slope):.2f}")
        print(f"    Top 5 eigenvalues: {[f'{x:.2f}' for x in evals[:5]]}")

    # ================================================================
    # 3. VARIANCE DECOMPOSITION: H vs perpendicular
    # ================================================================
    # Project each class onto H and compute residual
    # Use the leading eigenvector as the "H direction"
    v1 = evecs[:, 0]  # leading eigenvector
    H_proj = H_norm * np.dot(H_norm, v1) / np.dot(H_norm, H_norm) * v1  # H component
    # Actually: decompose the graph Laplacian's variance along H and perpendicular

    # Simpler: what fraction of total graph "energy" is along H?
    # H captures fraction = corr^2 of the variance in the leading eigenvector direction
    H_fraction = corrs[0]**2

    print(f"\n  3. VARIANCE DECOMPOSITION:")
    print(f"    H captures {H_fraction*100:.1f}% of leading eigenvector variance")
    print(f"    Residual (perpendicular to H): {(1-H_fraction)*100:.1f}%")

    # How many eigenvectors needed to capture 90% of "adjacency energy"?
    total_energy = sum(evals**2)
    cumul = 0
    for k in range(V):
        cumul += evals[k]**2
        if cumul / total_energy >= 0.90:
            print(f"    Eigenvectors for 90% energy: {k+1} (out of {V})")
            print(f"    Effective dimension (90% energy): {k+1}")
            break

    # ================================================================
    # 4. THE SQRT(2) CONNECTION
    # ================================================================
    pseudo_double = (2*n - 5) / (n - 2)
    print(f"\n  4. SQRT(2) CONNECTION:")
    print(f"    Pseudo-doubling ratio: (2n-5)/(n-2) = {pseudo_double:.4f}")
    print(f"    sqrt(2) = {np.sqrt(2):.4f}")
    print(f"    Ratio / sqrt(2) = {pseudo_double / np.sqrt(2):.4f}")
    print(f"    Approaches 2 = sqrt(2)^2 as n -> inf")

    # ================================================================
    # 5. OPPOSING TRENDS AS BETA ANALOGY
    # ================================================================
    # Scheme breaking (residual) ~ "holes" ~ beta_1
    # Markov accuracy (P2 coeff) ~ "filling" ~ inversely related to beta_3
    # Their product?
    residuals = {4: 0.106, 5: 0.229, 6: 0.318, 7: 0.355}
    p2_coeffs = {4: 0.628, 5: 0.834, 6: 0.902, 7: 0.939}
    if n in residuals:
        r = residuals[n]
        p = p2_coeffs[n]
        product = r * (1 - p)  # residual * (1 - P2 coeff) = "joint breaking"
        print(f"\n  5. OPPOSING TRENDS (beta analogy):")
        print(f"    Scheme residual (beta_1-like): {r:.3f}")
        print(f"    Markov deviation (1-P2, beta_3-like): {1-p:.3f}")
        print(f"    Product (seesaw?): {product:.4f}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

# Cross-n summary of opposing trends
print(f"\n{'='*60}")
print(f"  OPPOSING TRENDS SEESAW CHECK")
print(f"{'='*60}")

residuals = {4: 0.106, 5: 0.229, 6: 0.318, 7: 0.355}
p2_coeffs = {4: 0.628, 5: 0.834, 6: 0.902, 7: 0.939}

print(f"\n  {'n':>3} {'Residual':>9} {'1-P2':>8} {'Product':>9} {'Sum':>6}")
for nn in [4, 5, 6, 7]:
    r = residuals[nn]
    d = 1 - p2_coeffs[nn]
    print(f"  {nn:3d} {r:9.4f} {d:8.4f} {r*d:9.5f} {r+d:6.4f}")

print(f"""
If Product -> 0: SEESAW (can't have both large). Like beta_1*beta_3=0.
If Sum -> const: CONSERVATION LAW (total "breaking" is constant).

Product: {residuals[4]*(1-p2_coeffs[4]):.5f}, {residuals[5]*(1-p2_coeffs[5]):.5f}, {residuals[6]*(1-p2_coeffs[6]):.5f}, {residuals[7]*(1-p2_coeffs[7]):.5f}
-> Product DECREASES: 0.039 -> 0.038 -> 0.031 -> 0.022. SEESAW-LIKE!

The product approaches 0, suggesting the opposing trends ARE
analogous to beta_1 and beta_3: they can't BOTH be large.
As the scheme breaks more, the Markov accuracy compensates.
""")

print("DONE.")
print("=" * 80)
