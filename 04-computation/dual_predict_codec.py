#!/usr/bin/env python3
"""
dual_predict_codec.py — Dual-basis prediction: diagonal + orthogonal combined
kind-pasteur-2026-03-25-S20gr

THE IDEA: Each pixel has 8 neighbors in two coordinate systems:
  Orthogonal: up, down, left, right (standard row/col)
  Diagonal: up-left, up-right, down-left, down-right (45-degree)

Standard codecs predict from orthogonal neighbors only (PNG, JPEG-LS).
Our diagonal codec predicts from diagonal neighbors only.
COMBINE: predict from ALL 8 neighbors, pick the best prediction.

MIDDLE-OUT: Start from the CENTER of the block, expand outward.
Each ring of pixels gets predicted from the already-decoded inner ring.
The center pixel is the "keyframe" (stored raw).

This captures BOTH orthogonal and diagonal correlation simultaneously.
The residual (actual - predicted) should be smaller than either alone.

COMPARE against:
  1. Left-only prediction (PNG filter 1)
  2. Up-only prediction (PNG filter 2)
  3. Average of left+up (PNG filter 3 / Paeth-like)
  4. Diagonal-only prediction (our previous method)
  5. DUAL 8-neighbor prediction (NEW)
  6. Middle-out with dual prediction (NEW)
"""

import sys
import zlib
import numpy as np
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  DUAL-BASIS PREDICTION CODEC")
print("  kind-pasteur-2026-03-25-S20gr")
print("=" * 80)


def predict_left(M):
    """PNG-style left prediction."""
    H, W = M.shape
    residual = np.zeros_like(M)
    residual[:, 0] = M[:, 0]
    for j in range(1, W):
        residual[:, j] = M[:, j] - M[:, j-1]
    return residual


def predict_up(M):
    """PNG-style up prediction."""
    H, W = M.shape
    residual = np.zeros_like(M)
    residual[0, :] = M[0, :]
    for i in range(1, H):
        residual[i, :] = M[i, :] - M[i-1, :]
    return residual


def predict_avg_lu(M):
    """Average of left and up prediction."""
    H, W = M.shape
    residual = np.zeros_like(M)
    residual[0, :] = M[0, :]
    residual[:, 0] = M[:, 0]
    for i in range(1, H):
        for j in range(1, W):
            pred = (int(M[i, j-1]) + int(M[i-1, j])) // 2
            residual[i, j] = M[i, j] - pred
    return residual


def predict_diagonal(M):
    """Diagonal (up-left) prediction."""
    H, W = M.shape
    residual = np.zeros_like(M)
    residual[0, :] = M[0, :]
    residual[:, 0] = M[:, 0]
    for i in range(1, H):
        for j in range(1, W):
            residual[i, j] = M[i, j] - M[i-1, j-1]
    return residual


def predict_dual_8(M):
    """Dual 8-neighbor prediction: average of all available neighbors."""
    H, W = M.shape
    residual = np.zeros_like(M)
    for i in range(H):
        for j in range(W):
            neighbors = []
            for di in [-1, 0, 1]:
                for dj in [-1, 0, 1]:
                    if di == 0 and dj == 0:
                        continue
                    ni, nj = i + di, j + dj
                    if 0 <= ni < H and 0 <= nj < W:
                        # Only use already-scanned neighbors (raster order)
                        if ni < i or (ni == i and nj < j):
                            neighbors.append(int(M[ni, nj]))
            if neighbors:
                pred = sum(neighbors) // len(neighbors)
            else:
                pred = 0
            residual[i, j] = M[i, j] - pred
    return residual


def predict_middle_out(M):
    """Middle-out: start from center, expand in rings, predict from inner ring."""
    H, W = M.shape
    residual = np.zeros_like(M)
    decoded = np.full_like(M, -1, dtype=int)

    ci, cj = H // 2, W // 2
    decoded[ci, cj] = M[ci, cj]
    residual[ci, cj] = M[ci, cj]  # center stored raw

    max_ring = max(H, W)
    for ring in range(1, max_ring):
        for i in range(max(0, ci - ring), min(H, ci + ring + 1)):
            for j in range(max(0, cj - ring), min(W, cj + ring + 1)):
                if decoded[i, j] >= 0:
                    continue
                # Check if at this ring distance
                if max(abs(i - ci), abs(j - cj)) != ring:
                    continue
                # Predict from decoded neighbors
                neighbors = []
                for di in [-1, 0, 1]:
                    for dj in [-1, 0, 1]:
                        if di == 0 and dj == 0:
                            continue
                        ni, nj = i + di, j + dj
                        if 0 <= ni < H and 0 <= nj < W and decoded[ni, nj] >= 0:
                            neighbors.append(decoded[ni, nj])
                if neighbors:
                    pred = sum(neighbors) // len(neighbors)
                else:
                    pred = 128
                residual[i, j] = M[i, j] - pred
                decoded[i, j] = int(M[i, j])

    return residual


def compress_residual(residual):
    """Compress residual with zlib. Map to unsigned first."""
    # Residuals are in [-255, 255]. Map to [0, 510].
    mapped = (residual.astype(int) + 255).astype(np.uint8)
    # Clip to valid range
    mapped = np.clip(residual.astype(int) + 128, 0, 255).astype(np.uint8)
    return len(zlib.compress(bytes(mapped.flatten()), 9))


# ================================================================
# BENCHMARK
# ================================================================
print(f"\n  {'Method':>20} {'32x32':>8} {'64x64':>8} {'128x128':>9} {'256x256':>9}")

methods = {
    'Raw': lambda M: M,
    'Left': predict_left,
    'Up': predict_up,
    'Avg(L+U)': predict_avg_lu,
    'Diagonal': predict_diagonal,
    'Dual-8': predict_dual_8,
    'Middle-out': predict_middle_out,
}

for name, fn in methods.items():
    results = []
    for N in [32, 64, 128, 256]:
        np.random.seed(42)
        # Natural-looking image: smooth gradient + texture + noise
        x = np.linspace(0, 4*np.pi, N)
        y = np.linspace(0, 4*np.pi, N)
        X, Y = np.meshgrid(x, y)
        M = (128 + 50*np.sin(X)*np.cos(Y) + 20*np.sin(3*X+Y) + 5*np.random.randn(N, N)).clip(0, 255).astype(np.uint8)

        if N <= 128 or name not in ['Dual-8', 'Middle-out']:
            residual = fn(M)
            compressed = compress_residual(residual)
            results.append(compressed)
        else:
            results.append(None)

    line = f"  {name:>20}"
    for r in results:
        if r is not None:
            line += f" {r:>8}B"
        else:
            line += f" {'slow':>8}"
    print(line)

# ================================================================
print(f"\n  DETAILED COMPARISON (64x64 natural image):")
N = 64
np.random.seed(42)
x = np.linspace(0, 4*np.pi, N)
y = np.linspace(0, 4*np.pi, N)
X, Y = np.meshgrid(x, y)
M = (128 + 50*np.sin(X)*np.cos(Y) + 20*np.sin(3*X+Y) + 5*np.random.randn(N, N)).clip(0, 255).astype(np.uint8)

raw_size = len(zlib.compress(bytes(M.flatten()), 9))
print(f"  Raw (no prediction) + zlib: {raw_size}B")

for name, fn in methods.items():
    if name == 'Raw':
        continue
    residual = fn(M)
    comp = compress_residual(residual)
    improvement = (1 - comp / raw_size) * 100
    mean_abs_resid = np.mean(np.abs(residual.astype(float)))
    print(f"  {name:>20}: {comp:6}B ({improvement:+5.1f}% vs raw), mean|resid|={mean_abs_resid:.1f}")

# ================================================================
print(f"\n  VIDEO DELTA (64x64, small motion):")
np.random.seed(42)
frame1 = (128 + 50*np.sin(X)*np.cos(Y)).clip(0, 255).astype(np.uint8)
frame2 = np.clip(frame1.astype(int) + np.random.randint(-3, 4, frame1.shape), 0, 255).astype(np.uint8)
diff = (frame2.astype(int) - frame1.astype(int)).astype(np.int16)

print(f"  Raw delta + zlib: {len(zlib.compress(bytes((diff + 128).clip(0,255).astype(np.uint8).flatten()), 9))}B")
for name in ['Left', 'Up', 'Avg(L+U)', 'Diagonal', 'Dual-8', 'Middle-out']:
    fn = methods[name]
    # Apply prediction to the DELTA frame
    residual = fn(np.clip(diff + 128, 0, 255).astype(np.uint8))
    comp = compress_residual(residual)
    print(f"  {name:>20}: {comp:6}B")

print(f"\n{'='*60}")
print("THE DUAL-BASIS INSIGHT")
print(f"{'='*60}")
print("""
Each pixel lives at the intersection of TWO coordinate systems:
  Orthogonal: (row, col) — captured by Left/Up prediction
  Diagonal: (row+col, row-col) — captured by Diagonal prediction

The DUAL-8 predictor uses ALL 8 neighbors (both bases).
The MIDDLE-OUT predictor additionally changes the scan ORDER
so that each pixel is predicted from its NEAREST decoded neighbor
(not just raster-order predecessors).

If Dual-8 beats both Left and Diagonal individually:
  The two bases carry INDEPENDENT information about each pixel.
  Combining them reduces the residual below either alone.

If Middle-out beats Dual-8:
  The scan order matters — center-out gives better predictions
  because the center is most constrained (most neighbors available).
""")

print("DONE.")
print("=" * 80)
