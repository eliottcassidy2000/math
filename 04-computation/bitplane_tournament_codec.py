#!/usr/bin/env python3
"""
bitplane_tournament_codec.py — Bit-plane decomposition + tournament binary codec
kind-pasteur-2026-03-25-S20gl

THE FIX: Don't compress real matrices directly. Instead:
  1. Quantize NxN matrix to k-bit integers (0..2^k-1)
  2. Decompose into k BIT PLANES (MSB to LSB)
  3. Apply the BINARY tournament codec to each plane
  4. Progressive decode: MSB planes first = coarse, add LSB for detail

Each bit plane is an NxN binary matrix = two tournament tilings.
The binary codec gives: 21x delta compression, 76% trie sharing,
progressive spatial decode within each plane.

COMBINED PROGRESSIVE:
  - Across planes: MSB first (value resolution)
  - Within planes: early layers first (spatial resolution)
  - DOUBLY PROGRESSIVE: both value and space improve monotonically
"""

import sys
import numpy as np
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  BIT-PLANE TOURNAMENT CODEC")
print("  kind-pasteur-2026-03-25-S20gl")
print("=" * 80)


def quantize(M, bits):
    """Quantize float matrix to k-bit integers."""
    M = np.array(M)
    vmin, vmax = M.min(), M.max()
    if vmax == vmin:
        return np.zeros_like(M, dtype=int), vmin, vmax
    levels = (1 << bits) - 1
    Q = np.round((M - vmin) / (vmax - vmin) * levels).astype(int)
    Q = np.clip(Q, 0, levels)
    return Q, vmin, vmax


def dequantize(Q, vmin, vmax, bits):
    """Dequantize k-bit integers to floats."""
    levels = (1 << bits) - 1
    if levels == 0:
        return np.full_like(Q, vmin, dtype=float)
    return vmin + Q.astype(float) / levels * (vmax - vmin)


def extract_bitplanes(Q, bits):
    """Extract k bit planes from quantized matrix. MSB first."""
    planes = []
    for b in range(bits-1, -1, -1):
        plane = ((Q >> b) & 1).tolist()
        planes.append(plane)
    return planes


def combine_bitplanes(planes, bits):
    """Reconstruct quantized matrix from bit planes."""
    N = len(planes[0])
    Q = [[0]*N for _ in range(N)]
    for b_idx, plane in enumerate(planes):
        bit_pos = bits - 1 - b_idx
        for i in range(N):
            for j in range(N):
                Q[i][j] |= (plane[i][j] << bit_pos)
    return Q


def binary_tiling_encode(plane):
    """Encode binary plane as tournament tiling layers."""
    N = len(plane)
    lower_layers = []
    for i in range(1, N):
        layer = 0
        for j in range(i):
            if plane[i][j]:
                layer |= (1 << j)
        lower_layers.append(layer)
    upper_layers = []
    for j in range(1, N):
        layer = 0
        for i in range(j):
            if plane[i][j]:
                layer |= (1 << i)
        upper_layers.append(layer)
    diag = [plane[i][i] for i in range(N)]
    return lower_layers, upper_layers, diag


def binary_tiling_decode(lower_layers, upper_layers, diag, N):
    """Decode binary plane from tiling layers."""
    plane = [[0]*N for _ in range(N)]
    for i in range(N):
        plane[i][i] = diag[i]
    for i in range(1, N):
        for j in range(i):
            plane[i][j] = (lower_layers[i-1] >> j) & 1
    for j in range(1, N):
        for i in range(j):
            plane[i][j] = (upper_layers[j-1] >> i) & 1
    return plane


def delta_encode_layers(layers1, layers2):
    """Delta between two layer sequences."""
    return [l1 ^ l2 for l1, l2 in zip(layers1, layers2)]


def count_delta_bits(delta):
    return sum(bin(d).count('1') for d in delta)


def encode_matrix(M, bits=8):
    """Full encode: quantize + bitplane + tiling."""
    Q, vmin, vmax = quantize(M, bits)
    planes = extract_bitplanes(Q, bits)
    encoded_planes = []
    for plane in planes:
        lower, upper, diag = binary_tiling_encode(plane)
        encoded_planes.append((lower, upper, diag))
    return encoded_planes, vmin, vmax, bits


def decode_matrix(encoded_planes, vmin, vmax, bits, N, up_to_plane=None):
    """Decode with optional progressive (first k planes only)."""
    if up_to_plane is None:
        up_to_plane = len(encoded_planes)
    planes = []
    for p_idx in range(up_to_plane):
        lower, upper, diag = encoded_planes[p_idx]
        plane = binary_tiling_decode(lower, upper, diag, N)
        planes.append(plane)
    # Pad remaining planes with zeros (progressive lossy)
    while len(planes) < bits:
        planes.append([[0]*N for _ in range(N)])
    Q = combine_bitplanes(planes, bits)
    return dequantize(np.array(Q), vmin, vmax, bits)


# ================================================================
# DEMO
# ================================================================
print("\n  1. PROGRESSIVE VALUE DECODE:")
for N in [16, 32]:
    np.random.seed(42)
    M = np.array([[np.exp(-abs(i-j)/3.0) * (1 + 0.1*np.random.randn()) for j in range(N)] for i in range(N)])
    bits = 8

    encoded, vmin, vmax, b = encode_matrix(M.tolist(), bits)
    full_decode = decode_matrix(encoded, vmin, vmax, b, N)
    full_err = np.linalg.norm(M - full_decode) / np.linalg.norm(M)

    print(f"\n  {N}x{N}, {bits}-bit quantization (full error: {full_err*100:.2f}%):")
    print(f"  {'Planes':>6} {'Quality':>8} {'Bits':>10} {'Compress':>9}")
    for k in range(1, bits+1):
        partial = decode_matrix(encoded, vmin, vmax, b, N, up_to_plane=k)
        err = np.linalg.norm(M - partial) / np.linalg.norm(M)
        total_bits = k * N * N
        total_possible = bits * N * N
        compress = total_possible / total_bits
        print(f"  {k:6d} {(1-err)*100:7.2f}% {total_bits:10d} {compress:8.1f}x")

# ================================================================
print("\n  2. VIDEO DELTA (bitplane level):")
N = 32
bits = 8
np.random.seed(123)
frame1 = np.random.rand(N, N)
frame2 = frame1 + 0.02 * np.random.randn(N, N)  # small change
frame2 = np.clip(frame2, 0, 1)

enc1, v1min, v1max, b1 = encode_matrix(frame1.tolist(), bits)
enc2, v2min, v2max, b2 = encode_matrix(frame2.tolist(), bits)

total_delta_bits = 0
for p_idx in range(bits):
    l1, u1, d1 = enc1[p_idx]
    l2, u2, d2 = enc2[p_idx]
    dl = delta_encode_layers(l1, l2)
    du = delta_encode_layers(u1, u2)
    dd = [a ^ b for a, b in zip(d1, d2)]
    plane_delta = count_delta_bits(dl) + count_delta_bits(du) + sum(dd)
    total_delta_bits += plane_delta

full_bits = bits * N * N
print(f"  {N}x{N} x {bits}bit: full={full_bits} bits, delta={total_delta_bits} bits")
print(f"  Delta compression: {full_bits/max(total_delta_bits,1):.1f}x")
print(f"  (Small motion: 2% RMS noise on [0,1] values)")

# More motion levels
for noise_level in [0.01, 0.02, 0.05, 0.10, 0.20]:
    frame_b = frame1 + noise_level * np.random.randn(N, N)
    frame_b = np.clip(frame_b, 0, 1)
    enc_b, _, _, _ = encode_matrix(frame_b.tolist(), bits)
    delta_bits = 0
    for p_idx in range(bits):
        l1, u1, d1 = enc1[p_idx]
        lb, ub, db = enc_b[p_idx]
        delta_bits += count_delta_bits(delta_encode_layers(l1, lb))
        delta_bits += count_delta_bits(delta_encode_layers(u1, ub))
        delta_bits += sum(a ^ b for a, b in zip(d1, db))
    print(f"  noise={noise_level:.2f}: delta={delta_bits:5d} bits ({delta_bits/full_bits*100:.1f}%), {full_bits/max(delta_bits,1):.1f}x")

# ================================================================
print("\n  3. COMBINED PROGRESSIVE (planes + layers):")
N = 16
bits = 8
M = np.array([[np.exp(-abs(i-j)/3.0) for j in range(N)] for i in range(N)])
encoded, vmin, vmax, b = encode_matrix(M.tolist(), bits)

print(f"  {N}x{N} x {bits}bit: DOUBLY PROGRESSIVE decode")
print(f"  {'Planes':>6} {'Layers':>6} {'Bits':>8} {'Quality':>8} {'Compress':>9}")

for n_planes in [1, 2, 4, 8]:
    for n_layers in [2, N//2, N-1]:
        # Decode only first n_planes planes, each with first n_layers layers
        partial_planes = []
        total_bits = 0
        for p in range(n_planes):
            lower, upper, diag = encoded[p]
            # Truncate layers
            pl = binary_tiling_decode(lower[:n_layers], upper[:n_layers], diag, N)
            partial_planes.append(pl)
            total_bits += sum(range(1, n_layers+1)) * 2 + N  # approximate

        while len(partial_planes) < bits:
            partial_planes.append([[0]*N for _ in range(N)])

        Q_partial = combine_bitplanes(partial_planes, bits)
        M_partial = dequantize(np.array(Q_partial), vmin, vmax, bits)
        err = np.linalg.norm(M - M_partial) / np.linalg.norm(M)
        full_bits_total = bits * N * N
        compress = full_bits_total / max(total_bits, 1)

        print(f"  {n_planes:6d} {n_layers:6d} {total_bits:8d} {(1-err)*100:7.1f}% {compress:8.1f}x")

# ================================================================
print("\n  4. SPEED:")
for N in [16, 32, 64]:
    M = np.random.rand(N, N).tolist()
    reps = 100 if N <= 32 else 10
    t0 = time.time()
    for _ in range(reps):
        encode_matrix(M, 8)
    enc_time = (time.time() - t0) / reps * 1e3
    encoded, vmin, vmax, b = encode_matrix(M, 8)
    t0 = time.time()
    for _ in range(reps):
        decode_matrix(encoded, vmin, vmax, b, N)
    dec_time = (time.time() - t0) / reps * 1e3
    print(f"    {N}x{N} x 8bit: encode={enc_time:.1f}ms decode={dec_time:.1f}ms")

print("\nDONE.")
print("=" * 80)
