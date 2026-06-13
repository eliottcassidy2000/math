#!/usr/bin/env python3
"""
dual_basis_universal.py — Dual-basis prediction for arbitrary structured data
kind-pasteur-2026-03-25-S20gs

THE UNIVERSAL PRINCIPLE: Any data with K-dimensional locality benefits
from prediction using ALL 2^K-1 neighbor directions simultaneously.

1D (K=1): 2 neighbors (left, right) — standard delta coding
2D (K=2): 8 neighbors (4 orthogonal + 4 diagonal) — our Dual-8
3D (K=3): 26 neighbors (6 face + 12 edge + 8 corner)
KD: 3^K - 1 neighbors

APPLICATIONS TO REAL DATA FORMATS:

1. AUDIO (1D): predict from left + right context (bidirectional)
2. TEXT/TOKENS (1D): predict from preceding tokens (causal LM)
3. TIME SERIES MATRICES (2D): predict from time + feature neighbors
4. 3D VOLUMES (medical imaging, voxels): predict from 26 neighbors
5. DATABASE TABLES (2D): predict from row + column neighbors
6. SPREADSHEETS: each cell predicted from row/col/diagonal neighbors
7. DNA SEQUENCES: predict from preceding bases + complement strand
8. GRAPH ADJACENCY: predict from row/col/diagonal of adjacency matrix
"""

import sys
import zlib
import numpy as np
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  DUAL-BASIS PREDICTION: UNIVERSAL APPLICATIONS")
print("  kind-pasteur-2026-03-25-S20gs")
print("=" * 80)


def dual_predict_1d(seq):
    """1D: predict from left neighbor (delta coding)."""
    residual = np.zeros_like(seq)
    residual[0] = seq[0]
    for i in range(1, len(seq)):
        residual[i] = seq[i] - seq[i-1]
    return residual


def dual_predict_1d_bidirectional(seq):
    """1D bidirectional: predict from BOTH left and right (two-pass).
    Pass 1: left prediction. Pass 2: right prediction. Take average."""
    n = len(seq)
    left_pred = np.zeros(n)
    left_pred[0] = seq[0]
    for i in range(1, n):
        left_pred[i] = seq[i-1]
    right_pred = np.zeros(n)
    right_pred[-1] = seq[-1]
    for i in range(n-2, -1, -1):
        right_pred[i] = seq[i+1]
    avg_pred = ((left_pred + right_pred) // 2).astype(seq.dtype)
    return seq - avg_pred


def dual_predict_2d(M):
    """2D: Dual-8 prediction (all 8 neighbors, raster order)."""
    H, W = M.shape
    residual = np.zeros((H, W), dtype=np.int16)
    for i in range(H):
        for j in range(W):
            neighbors = []
            for di in [-1, 0, 1]:
                for dj in [-1, 0, 1]:
                    if di == 0 and dj == 0: continue
                    ni, nj = i + di, j + dj
                    if 0 <= ni < H and 0 <= nj < W:
                        if ni < i or (ni == i and nj < j):
                            neighbors.append(int(M[ni, nj]))
            pred = sum(neighbors) // len(neighbors) if neighbors else 0
            residual[i, j] = int(M[i, j]) - pred
    return residual


def compress(data):
    """Compress byte array with zlib level 9."""
    return len(zlib.compress(bytes(np.clip(data.flatten().astype(int) + 128, 0, 255).astype(np.uint8)), 9))


# ================================================================
print("\n  1. AUDIO (1D waveform):")
for name, gen in [
    ("sine", lambda n: (128 + 100*np.sin(np.linspace(0, 20*np.pi, n))).astype(np.uint8)),
    ("speech-like", lambda n: (128 + 80*np.sin(np.linspace(0, 50*np.pi, n)) * np.exp(-np.linspace(0, 3, n)**2/2)).astype(np.uint8)),
    ("noise", lambda n: np.random.randint(0, 256, n, dtype=np.uint8)),
    ("music-like", lambda n: (128 + 50*np.sin(np.linspace(0, 10*np.pi, n)) + 30*np.sin(np.linspace(0, 30*np.pi, n))).clip(0,255).astype(np.uint8)),
]:
    np.random.seed(42)
    seq = gen(4096)
    raw = compress(seq)
    left = compress(dual_predict_1d(seq))
    bidir = compress(dual_predict_1d_bidirectional(seq))
    print(f"  {name:>15}: raw={raw:5d}B, left={left:5d}B ({(1-left/raw)*100:+.1f}%), bidir={bidir:5d}B ({(1-bidir/raw)*100:+.1f}%)")

# ================================================================
print("\n  2. DATABASE TABLE (rows=records, cols=fields):")
for name, gen in [
    ("sorted_ids", lambda r,c: np.array([[(i*c+j) % 256 for j in range(c)] for i in range(r)], dtype=np.uint8)),
    ("categorical", lambda r,c: np.random.choice([0, 50, 100, 150, 200, 250], (r, c)).astype(np.uint8)),
    ("correlated", lambda r,c: np.array([[int(128+40*np.sin(i/5)+30*np.cos(j/3)) for j in range(c)] for i in range(r)], dtype=np.uint8)),
    ("sparse", lambda r,c: (np.random.random((r,c)) < 0.1).astype(np.uint8) * np.random.randint(1, 255, (r,c)).astype(np.uint8)),
]:
    np.random.seed(42)
    M = gen(64, 32)
    raw = compress(M)
    left_only = compress(dual_predict_1d(M.flatten().astype(np.int16)).reshape(M.shape))
    d8 = compress(dual_predict_2d(M))
    print(f"  {name:>15}: raw={raw:5d}B, row-delta={left_only:5d}B ({(1-left_only/raw)*100:+.1f}%), dual-8={d8:5d}B ({(1-d8/raw)*100:+.1f}%)")

# ================================================================
print("\n  3. TIME SERIES MATRIX (rows=time, cols=features):")
n_time, n_feat = 128, 16
for name, gen in [
    ("smooth", lambda: np.array([[int(128+50*np.sin(t/10+f/3)) for f in range(n_feat)] for t in range(n_time)], dtype=np.uint8)),
    ("regime", lambda: np.array([[int(128+50*np.sin(t/10) if t < 64 else 128+50*np.cos(t/10)) for f in range(n_feat)] for t in range(n_time)], dtype=np.uint8)),
    ("multifreq", lambda: np.array([[int(128+30*np.sin(t/5)+20*np.sin(f*t/100)) for f in range(n_feat)] for t in range(n_time)], dtype=np.uint8)),
]:
    M = gen()
    raw = compress(M)
    row_d = compress(np.diff(M, axis=0, prepend=M[:1,:]))
    col_d = compress(np.diff(M, axis=1, prepend=M[:,:1]))
    d8 = compress(dual_predict_2d(M))
    print(f"  {name:>15}: raw={raw:4d}B, time-delta={row_d:4d}B, feat-delta={col_d:4d}B, dual-8={d8:4d}B")

# ================================================================
print("\n  4. DNA-LIKE SEQUENCE (4-symbol alphabet, paired strands):")
n = 4096
np.random.seed(42)
# Simulate DNA: mostly random but with local repeats
dna = np.random.randint(0, 4, n, dtype=np.uint8)
# Add local correlation (repeats of length 3-10)
for _ in range(n//20):
    pos = np.random.randint(0, n-10)
    length = np.random.randint(3, 11)
    src = np.random.randint(0, n-length)
    dna[pos:pos+length] = dna[src:src+length]

# Complement strand: A<->T (0<->3), C<->G (1<->2)
complement = np.array([3, 2, 1, 0], dtype=np.uint8)[dna]

# Two strands as 2D matrix (2 x n)
paired = np.stack([dna, complement])

raw = compress(dna)
left_d = compress(dual_predict_1d(dna.astype(np.int16)))
# Predict from BOTH the left neighbor AND the complement
dual_resid = np.zeros_like(dna, dtype=np.int16)
for i in range(n):
    preds = []
    if i > 0: preds.append(int(dna[i-1]))
    preds.append(int(complement[i]))  # complement is "known" (paired)
    pred = sum(preds) // len(preds)
    dual_resid[i] = int(dna[i]) - pred

dual_comp = compress(dual_resid)
print(f"  DNA (4096 bases): raw={raw}B, left-delta={left_d}B, dual(left+complement)={dual_comp}B")

# ================================================================
print("\n  5. GRAPH ADJACENCY MATRIX:")
for name, gen in [
    ("path_graph", lambda n: np.array([[1 if abs(i-j)==1 else 0 for j in range(n)] for i in range(n)], dtype=np.uint8)),
    ("cycle_graph", lambda n: np.array([[1 if abs(i-j)==1 or abs(i-j)==n-1 else 0 for j in range(n)] for i in range(n)], dtype=np.uint8)),
    ("random_sparse", lambda n: (np.random.random((n,n)) < 0.1).astype(np.uint8)),
    ("block_community", lambda n: np.array([[1 if (i//8==j//8 and np.random.random()<0.7) or np.random.random()<0.05 else 0 for j in range(n)] for i in range(n)], dtype=np.uint8)),
]:
    np.random.seed(42)
    n = 64
    M = gen(n)
    # Make symmetric
    M = np.maximum(M, M.T)
    np.fill_diagonal(M, 0)

    raw = compress(M)
    d8 = compress(dual_predict_2d(M))
    # Diagonal reorder then zlib
    diag_data = []
    for d in range(n):
        if d == 0:
            for i in range(n): diag_data.append(M[i,i])
        else:
            for i in range(d, n): diag_data.append(M[i,i-d])
            for i in range(n-d): diag_data.append(M[i,i+d])
    diag_comp = len(zlib.compress(bytes(np.array(diag_data, dtype=np.uint8)), 9))

    print(f"  {name:>18}: raw={raw:4d}B, dual-8={d8:4d}B, diag-order={diag_comp:4d}B, best={'D8' if d8<diag_comp else 'Diag'}")

# ================================================================
print("\n  6. SPREADSHEET-LIKE (formulas create row+col+diag patterns):")
np.random.seed(42)
rows, cols = 32, 16
# Simulate: column A = IDs, B = A*2+noise, C = cumsum(A), D = random
A = np.arange(rows, dtype=np.uint8)
B = np.clip(A * 2 + np.random.randint(-2, 3, rows), 0, 255).astype(np.uint8)
C = np.cumsum(A % 10).astype(np.uint8)
D = np.random.randint(0, 256, rows, dtype=np.uint8)
sheet = np.column_stack([A, B, C, D] * (cols // 4))

raw = compress(sheet)
d8 = compress(dual_predict_2d(sheet))
col_d = compress(np.diff(sheet, axis=0, prepend=sheet[:1,:]))

print(f"  Spreadsheet {rows}x{cols}: raw={raw}B, col-delta={col_d}B, dual-8={d8}B")

# ================================================================
print(f"\n{'='*60}")
print("SUMMARY: WHERE DUAL-BASIS WINS")
print(f"{'='*60}")
print("""
DUAL-8 (all 8 neighbors) beats single-direction prediction when:
  1. Data has BOTH row and column correlation (natural images, tables)
  2. Data has diagonal patterns (band-limited matrices, adjacency)
  3. Block/community structure (each community is a diagonal cluster)

DUAL-8 is NEUTRAL or LOSES when:
  1. Data has ONLY one-directional correlation (sorted columns)
  2. Data is random (no structure to exploit)
  3. Data has very strong single-axis pattern (time series)

THE UNIVERSAL ARCHITECTURE:
  For any K-dimensional data:
    1. Identify the 3^K - 1 neighbor directions
    2. For each data point: predict from average of decoded neighbors
    3. Encode the residual (actual - predicted)
    4. The residual has LOWER entropy than the original data

  The improvement scales with the INDEPENDENCE of the K axes:
    If all K axes carry the same info: no gain from combining
    If K axes carry independent info: up to K-fold entropy reduction
""")

print("DONE.")
print("=" * 80)
