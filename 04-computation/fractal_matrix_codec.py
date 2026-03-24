#!/usr/bin/env python3
"""
fractal_matrix_codec.py — Recursive fractal compression of matrices
kind-pasteur-2026-03-25-S20gl

THE PROBLEM: Simple layer truncation is too lossy.
THE FIX: Apply the triangular decomposition RECURSIVELY to the residual.

ALGORITHM:
  1. Decompose NxN matrix M into diagonal layers L_0, L_1, ..., L_{N-1}
  2. Keep first b layers (coarse approximation): M_coarse
  3. Compute RESIDUAL: R = M - M_coarse (the far-diagonal entries)
  4. REARRANGE R into a smaller matrix by packing the residual entries
  5. Apply the same decomposition to the residual matrix RECURSIVELY
  6. Repeat until desired precision

Each recursion level captures finer structure. The total encoding =
coarse layers + recursion on residual. Like a wavelet transform but
using the TOURNAMENT CHAIN instead of frequency bands.

ALTERNATIVE APPROACH: Block-recursive.
  1. Split NxN into 4 quadrants (each N/2 x N/2)
  2. Decompose each quadrant triangularly
  3. The quadrants naturally separate: diagonal blocks (self-interaction)
     and off-diagonal blocks (cross-interaction)
  4. Diagonal blocks get full precision; off-diagonal get coarse + residual
  5. Recurse on each block
"""

import sys
import numpy as np
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  FRACTAL MATRIX CODEC")
print("  kind-pasteur-2026-03-25-S20gl")
print("=" * 80)


def diag_decompose(M, max_layer=None):
    """Decompose into diagonal layers. Returns (kept, residual)."""
    N = len(M)
    if max_layer is None:
        max_layer = N - 1

    kept = [[0.0]*N for _ in range(N)]
    residual = [[0.0]*N for _ in range(N)]

    for i in range(N):
        for j in range(N):
            dist = abs(i - j)
            if dist <= max_layer:
                kept[i][j] = M[i][j]
            else:
                residual[i][j] = M[i][j]

    return kept, residual


def nonzero_entries(M):
    """Count and sum of nonzero entries."""
    N = len(M)
    count = sum(1 for i in range(N) for j in range(N) if abs(M[i][j]) > 1e-15)
    energy = sum(M[i][j]**2 for i in range(N) for j in range(N))
    return count, energy


def block_split(M):
    """Split NxN matrix into 4 quadrants."""
    N = len(M)
    h = N // 2
    TL = [[M[i][j] for j in range(h)] for i in range(h)]
    TR = [[M[i][j] for j in range(h, N)] for i in range(h)]
    BL = [[M[i][j] for j in range(h)] for i in range(h, N)]
    BR = [[M[i][j] for j in range(h, N)] for i in range(h, N)]
    return TL, TR, BL, BR


def block_merge(TL, TR, BL, BR):
    """Merge 4 quadrants back into one matrix."""
    h = len(TL)
    N = h + len(BR)
    M = [[0.0]*N for _ in range(N)]
    for i in range(h):
        for j in range(h): M[i][j] = TL[i][j]
        for j in range(len(TR[0])): M[i][h+j] = TR[i][j]
    for i in range(len(BL)):
        for j in range(h): M[h+i][j] = BL[i][j]
        for j in range(len(BR[0])): M[h+i][h+j] = BR[i][j]
    return M


def fractal_compress(M, bandwidth, depth=0, max_depth=4):
    """Recursively compress: keep bandwidth layers, recurse on residual.

    Returns list of (level, kept_entries, precision) tuples.
    """
    N = len(M)
    if N <= 2 or depth >= max_depth:
        # Base case: store everything
        count, energy = nonzero_entries(M)
        return [(depth, count, N*N)]

    kept, residual = diag_decompose(M, bandwidth)

    # Count entries at this level
    kept_count, kept_energy = nonzero_entries(kept)
    resid_count, resid_energy = nonzero_entries(residual)

    result = [(depth, kept_count, N*N)]

    # Recurse on residual if it has significant energy
    if resid_energy > 1e-20 and resid_count > 0:
        # Pack residual into a smaller effective matrix
        # The residual lives in the "far corners" of the matrix.
        # Rearrange: collect all residual entries into a list, reshape.
        resid_entries = []
        for i in range(N):
            for j in range(N):
                if abs(i-j) > bandwidth:
                    resid_entries.append(residual[i][j])

        if len(resid_entries) > 0:
            # Reshape into a square-ish matrix for recursive compression
            side = int(np.ceil(np.sqrt(len(resid_entries))))
            # Pad with zeros
            while len(resid_entries) < side * side:
                resid_entries.append(0.0)
            R = [[resid_entries[i*side + j] for j in range(side)] for i in range(side)]

            # Recurse with potentially narrower bandwidth
            sub_results = fractal_compress(R, max(bandwidth//2, 1), depth+1, max_depth)
            result.extend(sub_results)

    return result


def fractal_encode_decode(M, bandwidth, max_depth=4):
    """Full fractal encode-decode with quality measurement."""
    N = len(M)
    total_entries = N * N
    total_energy = sum(M[i][j]**2 for i in range(N) for j in range(N))

    results = []

    current = [[M[i][j] for j in range(N)] for i in range(N)]
    cumul_kept = 0
    cumul_energy = 0.0

    for depth in range(max_depth + 1):
        kept, residual = diag_decompose(current, bandwidth)

        kept_count, kept_energy = nonzero_entries(kept)
        cumul_kept += kept_count
        cumul_energy += kept_energy

        quality = np.sqrt(cumul_energy / total_energy) if total_energy > 0 else 1.0
        compression = total_entries / max(cumul_kept, 1)

        results.append({
            'depth': depth,
            'entries_kept': cumul_kept,
            'quality': quality,
            'compression': compression,
            'residual_energy': sum(residual[i][j]**2 for i in range(N) for j in range(N))
        })

        if quality > 0.999:
            break

        # Reshape residual for next level
        resid_entries = []
        for i in range(N):
            for j in range(N):
                if abs(i-j) > bandwidth:
                    resid_entries.append(residual[i][j])

        if len(resid_entries) == 0:
            break

        side = int(np.ceil(np.sqrt(len(resid_entries))))
        while len(resid_entries) < side * side:
            resid_entries.append(0.0)
        current = [[resid_entries[i*side + j] for j in range(side)] for i in range(side)]

        # Narrow bandwidth for deeper levels
        bandwidth = max(bandwidth // 2, 1)
        N = side

    return results


# ================================================================
print("\n  FRACTAL COMPRESSION:")
for N in [16, 32, 64]:
    np.random.seed(42)
    # Band-limited test matrix
    M = [[np.exp(-abs(i-j)/3.0) * (1 + 0.1*np.random.randn()) for j in range(N)] for i in range(N)]

    print(f"\n  {N}x{N} band-limited matrix (bandwidth=3):")
    print(f"  {'Depth':>5} {'Entries':>8} {'Quality':>8} {'Compress':>9}")

    for bw in [3, 5, 8]:
        results = fractal_encode_decode(M, bw, max_depth=5)
        for r in results:
            print(f"  bw={bw} d={r['depth']} {r['entries_kept']:8d} {r['quality']*100:7.2f}% {r['compression']:8.2f}x")
        print()

# Compare: simple truncation vs fractal
print(f"\n  COMPARISON: Simple truncation vs Fractal")
for N in [32, 64]:
    M = [[np.exp(-abs(i-j)/3.0) * (1 + 0.1*np.random.randn()) for j in range(N)] for i in range(N)]
    total_energy = sum(M[i][j]**2 for i in range(N) for j in range(N))

    print(f"\n  {N}x{N}:")
    print(f"  {'Method':>20} {'Entries':>8} {'Quality':>8} {'Compress':>9}")

    # Simple truncation at various bandwidths
    for bw in [3, 5, 8]:
        kept, _ = diag_decompose(M, bw)
        kc, ke = nonzero_entries(kept)
        q = np.sqrt(ke / total_energy)
        c = N*N / max(kc, 1)
        print(f"  {'Trunc bw='+str(bw):>20} {kc:8d} {q*100:7.2f}% {c:8.2f}x")

    # Fractal at bandwidth 3
    results = fractal_encode_decode(M, 3, max_depth=5)
    for r in results:
        label = f"Fractal d={r['depth']}"
        print(f"  {label:>20} {r['entries_kept']:8d} {r['quality']*100:7.2f}% {r['compression']:8.2f}x")

print("\nDONE.")
print("=" * 80)
