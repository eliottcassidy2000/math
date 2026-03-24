#!/usr/bin/env python3
"""
tournament_media_codec.py — Binary matrix compression via tournament tiling
kind-pasteur-2026-03-25-S20gj

THE IDEA: Any N×N binary matrix = lower triangle + upper triangle + diagonal
= delta_N tiling + delta_{N-1} tiling + N diagonal bits.

The tournament chain structure gives:
  1. PROGRESSIVE: decode coarse structure first (early layers = low-skip)
  2. DELTA: consecutive frames share early layers -> trie compression
  3. LAYERED: natural multi-resolution (layer k = resolution level k)

For MEDIA:
  - Binary image (NxN): directly a square of two tilings
  - Grayscale: bit-plane decomposition -> stack of binary matrices
  - Video: sequence of binary matrices -> delta encoding between frames
  - Audio spectrogram: threshold to binary -> tournament compression

The key advantage: PROGRESSIVE DECODE with natural resolution ordering.
Early layers give the "coarse shape" (large-skip connections = long-range structure).
Late layers add "fine detail" (small-skip = local texture).

This inverts the usual image compression (DCT/wavelets give frequency bands).
Tournament compression gives SPATIAL HIERARCHY (core structure first).
"""

import sys
import numpy as np
import time
import random
from math import comb

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  TOURNAMENT MEDIA CODEC")
print("  kind-pasteur-2026-03-25-S20gj")
print("=" * 80)


def matrix_to_tilings(M):
    """Decompose NxN binary matrix into two tilings.

    Returns (lower_layers, upper_layers, diagonal, N) where:
      lower_layers: list of ints, layer k has k+1 bits (lower triangle)
      upper_layers: list of ints, layer k has k+1 bits (upper triangle)
      diagonal: list of N bits
    """
    N = len(M)
    assert all(len(row) == N for row in M), "Matrix must be square"

    # Lower triangle: M[i][j] for i > j
    lower_layers = []
    for i in range(1, N):
        layer = 0
        for j in range(i):
            if M[i][j]:
                layer |= (1 << j)
        lower_layers.append(layer)

    # Upper triangle: M[i][j] for i < j
    upper_layers = []
    for j in range(1, N):
        layer = 0
        for i in range(j):
            if M[i][j]:
                layer |= (1 << i)
        upper_layers.append(layer)

    # Diagonal
    diagonal = [M[i][i] for i in range(N)]

    return lower_layers, upper_layers, diagonal, N


def tilings_to_matrix(lower_layers, upper_layers, diagonal, N):
    """Reconstruct NxN binary matrix from two tilings."""
    M = [[0]*N for _ in range(N)]

    # Diagonal
    for i in range(N):
        M[i][i] = diagonal[i]

    # Lower triangle
    for i in range(1, N):
        for j in range(i):
            M[i][j] = (lower_layers[i-1] >> j) & 1

    # Upper triangle
    for j in range(1, N):
        for i in range(j):
            M[i][j] = (upper_layers[j-1] >> i) & 1

    return M


def encode_layers(layers):
    """Pack layers into bytes."""
    bits = 0
    offset = 0
    for k, layer in enumerate(layers):
        width = k + 1
        bits |= (layer & ((1 << width) - 1)) << offset
        offset += width
    n_bytes = (offset + 7) // 8
    return bits.to_bytes(n_bytes, 'little'), offset


def progressive_reconstruct(lower_layers, upper_layers, diagonal, N, up_to_layer):
    """Reconstruct partial matrix from first k layers."""
    k = min(up_to_layer, N-1)
    M = [[0]*N for _ in range(N)]
    for i in range(N):
        M[i][i] = diagonal[i]
    for i in range(1, k+1):
        for j in range(i):
            M[i][j] = (lower_layers[i-1] >> j) & 1
    for j in range(1, k+1):
        for i in range(j):
            M[i][j] = (upper_layers[j-1] >> i) & 1
    return M


def delta_layers(layers1, layers2):
    """XOR delta between two layer sequences."""
    return [l1 ^ l2 for l1, l2 in zip(layers1, layers2)]


def count_delta_bits(delta):
    """Count changed bits in delta."""
    return sum(bin(d).count('1') for d in delta)


# ================================================================
# DEMO
# ================================================================
print("\n  1. BASIC DECOMPOSITION:")
for N in [4, 8, 16, 32, 64]:
    # Random binary matrix
    M = [[random.randint(0, 1) for _ in range(N)] for _ in range(N)]
    lower, upper, diag, n = matrix_to_tilings(M)
    M2 = tilings_to_matrix(lower, upper, diag, n)
    assert M == M2, "Roundtrip failed!"

    lower_bytes, lower_bits = encode_layers(lower)
    upper_bytes, upper_bits = encode_layers(upper)
    raw_bits = N * N
    tiling_bits = lower_bits + upper_bits + N  # + diagonal
    tiling_bytes = len(lower_bytes) + len(upper_bytes) + (N + 7) // 8

    print(f"    {N}x{N}: raw={raw_bits} bits, tiling={tiling_bits} bits ({tiling_bits/raw_bits*100:.0f}%), {tiling_bytes} bytes")

# ================================================================
print("\n  2. PROGRESSIVE RECONSTRUCTION:")
N = 16
M = [[random.randint(0, 1) for _ in range(N)] for _ in range(N)]
lower, upper, diag, n = matrix_to_tilings(M)

for k in [2, 4, 8, 12, 15]:
    partial = progressive_reconstruct(lower, upper, diag, N, k)
    # Count correct pixels
    correct = sum(M[i][j] == partial[i][j] for i in range(N) for j in range(N))
    total = N * N
    # Bits used so far
    bits_used = sum(range(1, k+1)) * 2 + N  # lower + upper + diagonal
    print(f"    Layer {k:2d}/{N-1}: {correct}/{total} correct ({correct/total*100:.0f}%), {bits_used}/{N*N} bits ({bits_used/(N*N)*100:.0f}%)")

# ================================================================
print("\n  3. VIDEO DELTA COMPRESSION (consecutive frames):")
N = 32
frame1 = [[random.randint(0, 1) for _ in range(N)] for _ in range(N)]
# Frame 2: flip ~5% of pixels (small motion)
frame2 = [row[:] for row in frame1]
for _ in range(int(N*N*0.05)):
    i, j = random.randint(0, N-1), random.randint(0, N-1)
    frame2[i][j] = 1 - frame2[i][j]

l1, u1, d1, _ = matrix_to_tilings(frame1)
l2, u2, d2, _ = matrix_to_tilings(frame2)

dl = delta_layers(l1, l2)
du = delta_layers(u1, u2)
dd = [a ^ b for a, b in zip(d1, d2)]

delta_bits = count_delta_bits(dl) + count_delta_bits(du) + sum(dd)
full_bits = N * N
print(f"    Frame: {N}x{N} = {full_bits} bits")
print(f"    Delta (5% change): {delta_bits} bits ({delta_bits/full_bits*100:.1f}% of full)")
print(f"    Compression vs keyframe: {full_bits/max(delta_bits,1):.1f}x")

# More motion levels
for pct in [1, 5, 10, 25, 50]:
    frame_b = [row[:] for row in frame1]
    for _ in range(int(N*N*pct/100)):
        i, j = random.randint(0, N-1), random.randint(0, N-1)
        frame_b[i][j] = 1 - frame_b[i][j]
    lb, ub, db, _ = matrix_to_tilings(frame_b)
    dl = delta_layers(l1, lb)
    du = delta_layers(u1, ub)
    dd = [a ^ b for a, b in zip(d1, db)]
    delta_bits = count_delta_bits(dl) + count_delta_bits(du) + sum(dd)
    print(f"    {pct:2d}% change: delta = {delta_bits:4d} bits ({delta_bits/full_bits*100:.1f}%), {full_bits/max(delta_bits,1):.1f}x compression")

# ================================================================
print("\n  4. TRIE COMPRESSION (similar images):")
N = 16
base = [[random.randint(0, 1) for _ in range(N)] for _ in range(N)]
lb, ub, db, _ = matrix_to_tilings(base)

total_shared = 0
total_bits = 0
for trial in range(50):
    variant = [row[:] for row in base]
    # Change last few rows/cols (high layers)
    for _ in range(random.randint(5, 20)):
        i = random.randint(N//2, N-1)
        j = random.randint(0, N-1)
        variant[i][j] = 1 - variant[i][j]
    lv, uv, dv, _ = matrix_to_tilings(variant)

    # Shared prefix
    shared_l = 0
    for k in range(len(lb)):
        if lb[k] == lv[k]: shared_l += 1
        else: break
    shared_u = 0
    for k in range(len(ub)):
        if ub[k] == uv[k]: shared_u += 1
        else: break

    shared_bits = sum(range(1, shared_l+1)) + sum(range(1, shared_u+1))
    total_shared += shared_bits
    total_bits += N*N

avg_shared = total_shared / 50
print(f"    Base: {N}x{N}, 50 variants (bottom-half mutations)")
print(f"    Avg shared prefix bits: {avg_shared:.0f} / {N*N} = {avg_shared/(N*N)*100:.0f}%")
print(f"    Trie stores 50 images in ~{50*(N*N - avg_shared)/N/N:.1f} image equivalents")

# ================================================================
print("\n  5. SPEED:")
for N in [16, 32, 64, 128]:
    M = [[random.randint(0, 1) for _ in range(N)] for _ in range(N)]
    reps = 1000 if N <= 64 else 100
    t0 = time.time()
    for _ in range(reps):
        matrix_to_tilings(M)
    enc = (time.time() - t0) / reps * 1e6
    lower, upper, diag, n = matrix_to_tilings(M)
    t0 = time.time()
    for _ in range(reps):
        tilings_to_matrix(lower, upper, diag, n)
    dec = (time.time() - t0) / reps * 1e6
    print(f"    {N}x{N}: encode={enc:.0f}us decode={dec:.0f}us")

print("\nDONE.")
print("=" * 80)
