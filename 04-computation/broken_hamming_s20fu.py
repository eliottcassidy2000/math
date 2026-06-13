#!/usr/bin/env python3
"""
broken_hamming_s20fu.py — How exactly does the Hamming scheme break at each n?
kind-pasteur-2026-03-24-S20fu

The Hamming scheme H(m,2) on Q_m has:
  W_a * W_b = sum_c p_{ab}^c W_c (exact, no residual)

The S_n quotient BREAKS this. Measure the breaking at each n:
  1. Residual fraction: ||W_1^2 - best_fit||_F / ||W_1^2||_F
  2. Markov approximation: P_1^2 ~ alpha * P_2 + rest
  3. The P_1^k chain: how well does P_1^k ~ mixture of P_d's?
  4. Compare the "breaking" across n=4,5,6,7

HYPOTHESIS: The breaking INCREASES with n (the quotient becomes less
scheme-like as n grows). This would parallel:
  - SC backbone fragmenting at n=8
  - Real roots failing at n=9
  - Claw-free failing at n=9
  - Skip hierarchy gap growing

All pointing to: the metagraph becomes MORE complex and LESS regular at larger n.
"""

import sys
import numpy as np
from math import factorial, comb
from itertools import permutations, combinations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  THE BROKEN HAMMING SCHEME AT EACH n")
print("  kind-pasteur-2026-03-24-S20fu")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

results = {}

for n in range(4, 8):
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

    print(f"\n  n = {n}, m = {m}...", end=" ", flush=True)
    canon_map = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = b2a(bits)
        canon_map[mask] = canonicalize(A)

    classes = sorted(set(canon_map.values()))
    V = len(classes)
    cidx = {cn: i for i, cn in enumerate(classes)}
    tilings = np.array([sum(1 for mask, cn in canon_map.items() if cn == c) for c in classes], dtype=float)

    # Build W_1 and W_2
    W1 = np.zeros((V, V))
    W2 = np.zeros((V, V))
    for mask in range(1 << m):
        i = cidx[canon_map[mask]]
        for wi in range(m):
            j = cidx[canon_map[mask ^ (1 << wi)]]
            W1[i, j] += 1
        for wi, wj in combinations(range(m), 2):
            j = cidx[canon_map[mask ^ (1 << wi) ^ (1 << wj)]]
            W2[i, j] += 1

    # Markov matrices
    D_inv = np.diag(1.0 / np.maximum(tilings, 1))
    P1 = D_inv @ W1 / m
    P2 = D_inv @ W2 / comb(m, 2)

    # ================================================================
    # MEASURE 1: W_1^2 residual as fraction
    # ================================================================
    W1W1 = W1 @ W1
    # Best fit: project W1^2 onto span of {I, W1, W2}
    n_entries = V * V
    X = np.column_stack([np.eye(V).flatten(), W1.flatten(), W2.flatten()])
    y = W1W1.flatten()
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    fit = X @ coeffs
    residual_frac = np.linalg.norm(fit - y) / np.linalg.norm(y)

    # ================================================================
    # MEASURE 2: P_1^2 ~ alpha*P_2 coefficient
    # ================================================================
    P1P1 = P1 @ P1
    # Best scalar: P1^2 ~ alpha*P_2 + beta*I + gamma*P_1
    X_p = np.column_stack([np.eye(V).flatten(), P1.flatten(), P2.flatten()])
    y_p = P1P1.flatten()
    coeffs_p, _, _, _ = np.linalg.lstsq(X_p, y_p, rcond=None)
    alpha_I, alpha_P1, alpha_P2 = coeffs_p
    fit_p = X_p @ coeffs_p
    markov_residual = np.max(np.abs(fit_p - y_p))

    # ================================================================
    # MEASURE 3: Spectral gap and eigenvalue comparison
    # ================================================================
    evals_W1 = sorted(np.linalg.eigvalsh(W1), reverse=True)
    evals_W2 = sorted(np.linalg.eigvalsh(W2), reverse=True)

    # In a true Hamming scheme: evals of W_d are Krawtchouk K_j(d)
    # The ratio evals_W2[j] / evals_W1[j] would be K_j(2)/K_j(1) = (m-2j-1)
    # Check how well this holds
    if V > 2:
        ratios = []
        for j in range(min(5, V)):
            if abs(evals_W1[j]) > 1:
                ratios.append(evals_W2[j] / evals_W1[j])
        if ratios:
            kraw_expected = [(m - 2*j - 1) for j in range(len(ratios))]

    # ================================================================
    # MEASURE 4: Eigenvalue interlacing quality
    # ================================================================
    # In an association scheme: eigenvalues of W_d are exactly determined
    # by eigenvalues of W_1 via Krawtchouk polynomials.
    # Compute: how well do W_2 eigenvalues predict from W_1?
    # K_j(2,m) = 1 - 4j(m-j)/(m(m-1)) ... approximate
    # Actually: K_j(2) = C(m,2) - 2j(m-j) ... check
    # K_j(d) = sum_s (-1)^s C(d,s) C(m-d, j-s) for 0<=s<=min(j,d)
    def krawtchouk(j_val, d_val, m_val):
        total = 0
        for s in range(min(j_val, d_val) + 1):
            if j_val - s <= m_val - d_val and j_val - s >= 0:
                total += (-1)**s * comb(d_val, s) * comb(m_val - d_val, j_val - s)
        return total

    # If this were a Hamming scheme: W_d eigenvalues = K_j(d) for j=0..m
    # In the quotient: we have V < m+1 eigenvalues. No direct match.

    results[n] = {
        'V': V, 'm': m,
        'residual_frac': residual_frac,
        'alpha_I': alpha_I, 'alpha_P1': alpha_P1, 'alpha_P2': alpha_P2,
        'markov_residual': markov_residual,
        'evals_W1': evals_W1[:5],
        'evals_W2': evals_W2[:5],
    }

    print(f"V={V}, residual={residual_frac:.4f}, P1^2~{alpha_P2:.3f}*P2+{alpha_P1:.3f}*P1+{alpha_I:.3f}*I, max_err={markov_residual:.4f} ({time.time()-t0:.0f}s)")

# ================================================================
# SUMMARY TABLE
# ================================================================
print(f"\n{'='*70}")
print(f"  BROKEN HAMMING SCHEME: SUMMARY ACROSS n")
print(f"{'='*70}")

print(f"\n  {'n':>3} {'m':>4} {'V':>6} {'Residual%':>10} {'P2 coeff':>9} {'P1 coeff':>9} {'I coeff':>8} {'MaxErr':>8}")
for n in sorted(results.keys()):
    r = results[n]
    print(f"  {n:3d} {r['m']:4d} {r['V']:6d} {r['residual_frac']*100:9.2f}% {r['alpha_P2']:9.4f} {r['alpha_P1']:9.4f} {r['alpha_I']:8.4f} {r['markov_residual']:8.4f}")

print(f"""
INTERPRETATION:
  - Residual% measures how badly W_1^2 fails to lie in span{{I, W_1, W_2}}.
  - P2 coefficient measures how much of P_1^2 is "pure d=2 step."
  - If the scheme were exact: residual = 0, P2 coeff = known formula.

TREND:
  Does the breaking INCREASE or DECREASE with n?
  If residual% increases: the metagraph becomes LESS regular.
  If P2 coefficient stabilizes: the Markov structure is robust.
""")

# Compare with known asymmetries
print(f"\n  COMPARISON WITH OTHER n-DEPENDENT PHENOMENA:")
print(f"  {'n':>3} {'Residual%':>10} {'SC frag?':>9} {'Skip gap':>9} {'Overlap':>8} {'Grid-sym exp':>12}")
grid_exp = {4:-1, 5:-2, 6:-4, 7:-6}
skip_gap = {4:'2.0x', 5:'~2x', 6:'1.6x', 7:'14x'}
overlap = {4:0, 5:0, 6:5, 7:0}
sc_frag = {4:'N', 5:'N', 6:'N', 7:'N'}
for n in sorted(results.keys()):
    r = results[n]
    print(f"  {n:3d} {r['residual_frac']*100:9.2f}% {sc_frag.get(n,'?'):>9} {skip_gap.get(n,'?'):>9} {overlap.get(n,'?'):>8} {grid_exp.get(n,'?'):>12}")

print("\nDONE.")
print("=" * 80)
