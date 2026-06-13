#!/usr/bin/env python3
"""
waggly_scheme_s20ft.py — Waggly as an association scheme quotient
kind-pasteur-2026-03-24-S20ft

THE HAMMING SCHEME H(m, 2):
  Vertices: {0,1}^m (binary strings of length m)
  Relations: R_d = {(x,y) : Hamming(x,y) = d} for d = 0,1,...,m
  This is an ASSOCIATION SCHEME with m+1 classes.

  Eigenvalues: Krawtchouk polynomials K_j(d) = sum_s (-1)^s C(d,s) C(m-d,j-s)
  The adjacency matrix A_d has eigenvalues K_j(d) for j=0,...,m.

S_n ACTS on H(m,2) where m = C(n-1,2), preserving the Hamming distance.
The QUOTIENT Q_m / S_n inherits the association scheme structure.

QUESTION: Is the quotient Q_m / S_n itself an association scheme?
  If yes: its eigenvalues are determined by Krawtchouk polynomials
  restricted to the S_n-invariant subspace.
  This would give us the spectral gap, mixing time, etc. FOR FREE.

CONCRETE TEST: Compute the distance matrices at each d for the quotient
and check if they satisfy the association scheme axioms:
  A_d * A_e = sum_f p_{de}^f A_f (structure constants)

Also: the CODING THEORY interpretation.
  Each iso class C is a "codeword" in the quotient.
  The "minimum distance" d_min = minimum d such that d=d_min covers some edge.
  The "covering radius" = maximum d needed = diameter = k*.
"""

import sys
import numpy as np
from math import factorial, comb
from itertools import permutations, combinations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  WAGGLY AS ASSOCIATION SCHEME QUOTIENT")
print("  kind-pasteur-2026-03-24-S20ft")
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

    canon_map = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = b2a(bits)
        canon_map[mask] = canonicalize(A)

    classes = sorted(set(canon_map.values()))
    V = len(classes)
    cidx = {cn: i for i, cn in enumerate(classes)}
    class_masks = defaultdict(list)
    for mask, cn in canon_map.items():
        class_masks[cn].append(mask)

    tilings = np.array([len(class_masks[cn]) for cn in classes], dtype=float)

    print(f"\n{'#'*60}")
    print(f"  n = {n}, m = {m}, V = {V}")
    print(f"{'#'*60}")

    # Build DISTANCE-d WEIGHT MATRICES W_d[i,j] for the quotient
    max_d = min(m, 6)
    W_d = {}

    for d in range(1, max_d + 1):
        Wd = np.zeros((V, V))
        for mask in range(1 << m):
            i = cidx[canon_map[mask]]
            for combo in combinations(range(m), d):
                flip = mask
                for wi in combo:
                    flip ^= (1 << wi)
                j = cidx[canon_map[flip]]
                Wd[i, j] += 1
        W_d[d] = Wd
        print(f"  Built W_{d} ({time.time()-t0:.1f}s)")

    # ================================================================
    # CHECK: Is the quotient an association scheme?
    # Test: W_a @ W_b = sum_c p_{ab}^c W_c (up to diagonal correction)
    # ================================================================
    print(f"\n  ASSOCIATION SCHEME TEST:")
    print(f"  Testing W_1 @ W_1 = sum_d p_d * W_d + diag correction")

    W1W1 = W_d[1] @ W_d[1]

    # Try to express W1W1 as linear combination of W_d's and identity
    # W1W1 = sum_d alpha_d W_d + beta I
    # This is a linear system

    # Build the matrix equation
    n_mats = max_d + 1  # W_0=I, W_1, ..., W_{max_d}
    n_entries = V * V

    X = np.zeros((n_entries, n_mats))
    X[:, 0] = np.eye(V).flatten()  # W_0 = I
    for d in range(1, max_d + 1):
        X[:, d] = W_d[d].flatten()

    y = W1W1.flatten()

    # Least squares
    coeffs, residual, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    fit = X @ coeffs
    max_error = np.max(np.abs(fit - y))

    print(f"  W_1 @ W_1 = {' + '.join(f'{coeffs[d]:.2f}*W_{d}' for d in range(n_mats) if abs(coeffs[d]) > 0.01)}")
    print(f"  Max fit error: {max_error:.6f}")
    print(f"  Perfect fit (association scheme)? {max_error < 0.01}")

    if max_error > 0.01:
        print(f"  NOT an association scheme (W_1*W_1 not in span of W_d's)")
        # What's the residual structure?
        resid = (W1W1 - fit.reshape(V, V))
        print(f"  Residual Frobenius: {np.linalg.norm(resid):.4f}")
        print(f"  Residual as fraction of W1W1: {np.linalg.norm(resid)/np.linalg.norm(W1W1)*100:.2f}%")

    # ================================================================
    # NORMALIZED QUOTIENT: check if W_d / (tilings outer product) is cleaner
    # ================================================================
    print(f"\n  NORMALIZED TEST (remove tiling weights):")
    D_t_inv = np.diag(1.0 / np.maximum(tilings, 1))
    P1 = D_t_inv @ W_d[1] / m  # Markov transition at d=1

    P1P1 = P1 @ P1
    Pd = {}
    for d in range(1, max_d + 1):
        Pd[d] = D_t_inv @ W_d[d] / comb(m, d)  # normalized at each d

    # Try P1^2 as combo of P_d
    X_norm = np.zeros((n_entries, n_mats))
    X_norm[:, 0] = np.eye(V).flatten()
    for d in range(1, max_d + 1):
        X_norm[:, d] = Pd[d].flatten()

    y_norm = P1P1.flatten()
    coeffs_norm, _, _, _ = np.linalg.lstsq(X_norm, y_norm, rcond=None)
    fit_norm = X_norm @ coeffs_norm
    error_norm = np.max(np.abs(fit_norm - y_norm))

    print(f"  P_1^2 = {' + '.join(f'{coeffs_norm[d]:.4f}*P_{d}' for d in range(n_mats) if abs(coeffs_norm[d]) > 0.001)}")
    print(f"  Max error: {error_norm:.6f}")
    print(f"  Clean? {error_norm < 0.001}")

    # ================================================================
    # KRAWTCHOUK CHECK: eigenvalues of W_d
    # ================================================================
    print(f"\n  KRAWTCHOUK EIGENVALUE CHECK:")
    for d in range(1, min(4, max_d+1)):
        evals = sorted(np.linalg.eigvalsh(W_d[d].astype(float)), reverse=True)
        # Krawtchouk: K_j(d; m) = sum_s (-1)^s C(d,s) C(m-d, j-s)
        def krawtchouk(j, d_val, m_val):
            total = 0
            for s in range(min(j, d_val) + 1):
                total += (-1)**s * comb(d_val, s) * comb(m_val - d_val, j - s)
            return total

        kraw_evals = sorted([krawtchouk(j, d, m) * tilings.sum() / V for j in range(V)], reverse=True)
        # The actual eigenvalues of W_d in the quotient won't exactly match Krawtchouk
        # because the quotient breaks the Hamming scheme symmetry

        print(f"  W_{d} eigenvalues (top 5): {[f'{x:.1f}' for x in evals[:5]]}")
        if d == 1:
            print(f"  Krawtchouk K_j(1,{m}) * 2^m/V: {[f'{x:.1f}' for x in kraw_evals[:5]]}")

    # ================================================================
    # CODING INTERPRETATION
    # ================================================================
    print(f"\n  CODING THEORY INTERPRETATION:")
    print(f"  Each tiling is a 'codeword' in {{0,1}}^{m}")
    print(f"  ISO classes are 'code orbits' under S_{n}")
    print(f"  The 'code' C = set of all tilings has:")
    print(f"    Size: |C| = 2^{m} = {2**m}")
    print(f"    Dimension: {m} (it's the full space)")

    # Weight distribution: how many tilings at each Hamming weight?
    weight_dist = Counter()
    for mask in range(1 << m):
        hw = bin(mask).count('1')
        weight_dist[hw] += 1

    print(f"  Hamming weight distribution of tilings:")
    for w in sorted(weight_dist.keys()):
        print(f"    weight {w}: {weight_dist[w]} tilings = C({m},{w})")

    # Class-level: average Hamming weight per class
    print(f"\n  Average Hamming weight by class (sorted by H):")
    # (omit detail for brevity, just summary)
    class_hw = {}
    for cn in classes:
        total_hw = sum(bin(mask).count('1') for mask in class_masks[cn])
        class_hw[cn] = total_hw / len(class_masks[cn])

    hw_vals = list(class_hw.values())
    print(f"    Range: {min(hw_vals):.2f} to {max(hw_vals):.2f}")
    print(f"    Mean: {sum(hw_vals)/len(hw_vals):.2f} (expected: m/2 = {m/2})")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
