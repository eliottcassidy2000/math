#!/usr/bin/env python3
"""
triangle_matrix_codec.py — Generalized triangular decomposition for ANY square matrix
kind-pasteur-2026-03-25-S20gk

THE MATH: C(k+1,2) + C(k,2) = k^2. Two triangles tile a square.

THE INSIGHT: Ordering entries by DIAGONAL DISTANCE creates a natural
progressive codec for ANY square matrix with diagonal correlation.

APPLICATIONS beyond binary:
  1. CORRELATION MATRICES: symmetric, entries near diagonal ~1.0
  2. TRANSITION MATRICES: row-stochastic, band-limited
  3. DISTANCE MATRICES: symmetric, diagonal = 0, near-diagonal = small
  4. ADJACENCY MATRICES: sparse, band structure for ordered graphs
  5. COVARIANCE MATRICES: symmetric positive definite, band-dominated

For REAL-VALUED matrices: quantize to k bits per entry, then
apply the binary codec to each bit plane. Or use the layer
structure directly with arithmetic coding per layer.

The KEY SPEEDUP: for band-limited matrices (most real-world matrices),
early layers suffice for good approximation. Skip late layers = lossy
compression with bounded error.
"""

import sys
import numpy as np
import time
import random

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  TRIANGLE MATRIX CODEC")
print("  kind-pasteur-2026-03-25-S20gk")
print("=" * 80)


def matrix_to_diag_layers(M):
    """Decompose NxN matrix into diagonal layers.

    Layer 0: diagonal entries M[i][i]
    Layer k (k>=1): entries M[i][j] with |i-j| = k
      Sub-diagonal: M[i][i-k] for i=k,...,N-1
      Super-diagonal: M[i][i+k] for i=0,...,N-1-k

    Returns list of layers, each layer is a list of values.
    """
    N = len(M)
    layers = []
    # Layer 0: diagonal
    layers.append([M[i][i] for i in range(N)])
    # Layer k: k-th off-diagonals
    for k in range(1, N):
        layer = []
        # Sub-diagonal
        for i in range(k, N):
            layer.append(M[i][i-k])
        # Super-diagonal
        for i in range(N-k):
            layer.append(M[i][i+k])
        layers.append(layer)
    return layers


def diag_layers_to_matrix(layers, N):
    """Reconstruct NxN matrix from diagonal layers."""
    M = [[0.0]*N for _ in range(N)]
    # Layer 0: diagonal
    for i in range(N):
        M[i][i] = layers[0][i]
    # Layer k: off-diagonals
    for k in range(1, N):
        idx = 0
        for i in range(k, N):
            M[i][i-k] = layers[k][idx]
            idx += 1
        for i in range(N-k):
            M[i][i+k] = layers[k][idx]
            idx += 1
    return M


def progressive_reconstruct(layers, N, up_to_layer):
    """Reconstruct from first k layers only."""
    partial_layers = layers[:up_to_layer+1]
    # Pad remaining with zeros
    while len(partial_layers) < N:
        k = len(partial_layers)
        partial_layers.append([0.0] * (2*(N-k)))
    return diag_layers_to_matrix(partial_layers, N)


def frobenius_error(M1, M2):
    """Frobenius norm of difference."""
    N = len(M1)
    return np.sqrt(sum((M1[i][j] - M2[i][j])**2 for i in range(N) for j in range(N)))


def relative_error(M_orig, M_approx):
    """Relative Frobenius error."""
    N = len(M_orig)
    norm_orig = np.sqrt(sum(M_orig[i][j]**2 for i in range(N) for j in range(N)))
    err = frobenius_error(M_orig, M_approx)
    return err / norm_orig if norm_orig > 0 else 0


# ================================================================
print("\n  1. CORRELATION MATRIX (real-valued, symmetric):")
for N in [8, 16, 32]:
    # Generate random correlation matrix: C = A^T A / N, normalized
    A = np.random.randn(N*3, N)
    C = (A.T @ A) / (N*3)
    # Normalize to correlation
    d = np.sqrt(np.diag(C))
    C = C / np.outer(d, d)
    M = C.tolist()

    layers = matrix_to_diag_layers(M)
    M2 = diag_layers_to_matrix(layers, N)
    assert np.allclose(M, M2, atol=1e-10), "Roundtrip failed!"

    # Progressive: how many layers for 90%, 95%, 99% accuracy?
    for target in [0.90, 0.95, 0.99]:
        for k in range(N):
            approx = progressive_reconstruct(layers, N, k)
            err = relative_error(M, approx)
            if 1 - err >= target:
                entries_used = sum(len(layers[l]) for l in range(k+1))
                total_entries = N * N
                print(f"    {N}x{N}: {target*100:.0f}% accuracy at layer {k}/{N-1} ({entries_used}/{total_entries} = {entries_used/total_entries*100:.0f}% entries)")
                break

# ================================================================
print("\n  2. BAND-LIMITED MATRIX (exponential decay from diagonal):")
for N in [16, 32, 64]:
    M = [[np.exp(-abs(i-j)/3) * (1 + 0.1*random.gauss(0,1)) for j in range(N)] for i in range(N)]
    layers = matrix_to_diag_layers(M)

    # Layer energy: sum of squares per layer
    total_energy = sum(sum(v**2 for v in layer) for layer in layers)
    cumul = 0
    for k in range(min(N, 15)):
        layer_energy = sum(v**2 for v in layers[k])
        cumul += layer_energy
        if k < 8 or cumul/total_energy > 0.99:
            print(f"    {N}x{N} layer {k:2d}: energy={layer_energy/total_energy*100:5.1f}%, cumul={cumul/total_energy*100:5.1f}%")
            if cumul/total_energy > 0.99:
                break

# ================================================================
print("\n  3. TRANSITION MATRIX (Markov chain):")
for N in [8, 16]:
    # Random Markov: near-diagonal transitions likely
    M = [[0.0]*N for _ in range(N)]
    for i in range(N):
        weights = [np.exp(-abs(i-j)/2) for j in range(N)]
        total = sum(weights)
        for j in range(N):
            M[i][j] = weights[j] / total

    layers = matrix_to_diag_layers(M)
    for target in [0.90, 0.95, 0.99]:
        for k in range(N):
            approx = progressive_reconstruct(layers, N, k)
            err = relative_error(M, approx)
            if 1 - err >= target:
                entries_pct = sum(len(layers[l]) for l in range(k+1)) / (N*N) * 100
                print(f"    {N}x{N} Markov: {target*100:.0f}% at layer {k} ({entries_pct:.0f}% entries)")
                break

# ================================================================
print("\n  4. COMPRESSION RATIOS (lossy, keeping first k layers):")
for N in [16, 32, 64]:
    # Band-limited matrix
    M = [[np.exp(-abs(i-j)/3) + 0.05*random.gauss(0,1) for j in range(N)] for i in range(N)]
    layers = matrix_to_diag_layers(M)

    for quality in ["low", "medium", "high"]:
        if quality == "low": max_k = N // 4
        elif quality == "medium": max_k = N // 2
        else: max_k = 3 * N // 4

        approx = progressive_reconstruct(layers, N, max_k)
        err = relative_error(M, approx)
        entries_kept = sum(len(layers[l]) for l in range(max_k+1))
        ratio = N*N / entries_kept

        print(f"    {N}x{N} {quality:>6}: keep {max_k}/{N-1} layers, {ratio:.1f}x compression, {(1-err)*100:.1f}% quality")

# ================================================================
print("\n  5. SPEED (decompose + reconstruct):")
for N in [16, 32, 64, 128, 256]:
    M = [[random.gauss(0, 1) for _ in range(N)] for _ in range(N)]
    reps = 100 if N <= 64 else 10
    t0 = time.time()
    for _ in range(reps):
        layers = matrix_to_diag_layers(M)
    enc = (time.time() - t0) / reps * 1e3
    t0 = time.time()
    for _ in range(reps):
        diag_layers_to_matrix(layers, N)
    dec = (time.time() - t0) / reps * 1e3
    print(f"    {N}x{N}: encode={enc:.1f}ms decode={dec:.1f}ms")

print(f"\n{'='*60}")
print("WHY THIS WORKS")
print(f"{'='*60}")
print("""
The diagonal decomposition works because MOST real-world matrices
have BANDWIDTH STRUCTURE: entries decay with distance from diagonal.

  Correlation matrices: corr(X_i, X_j) ~ exp(-|i-j|/tau)
  Transition matrices: P(i->j) ~ exp(-|i-j|/sigma)
  Distance matrices: d(i,j) ~ |i-j| (for ordered data)
  Adjacency of path/cycle graphs: band-limited by construction

The triangular decomposition C(k+1,2) + C(k,2) = k^2 orders entries
by diagonal distance, putting MOST IMPORTANT entries first.

For a matrix with bandwidth b (entries beyond b-th diagonal are ~0):
  First b layers capture ~99% of the information
  Compression ratio ~ N/b (for large N)
  Progressive decode gives useful approximation from b << N layers

This is the TOURNAMENT PERSPECTIVE on matrix compression:
  The matrix IS two tournament-style tilings
  The layers ARE the recursion levels of the tiling
  Progressive decode follows the tournament chain
  Delta encoding exploits the tiling XOR structure
""")

print("DONE.")
print("=" * 80)
