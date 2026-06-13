#!/usr/bin/env python3
"""
dual_basis_applications.py — The 45-degree trick applied abstractly
kind-pasteur-2026-03-25-S20gq

THE ABSTRACT PRINCIPLE:
  Any NxN matrix has TWO natural coordinate systems:
    Standard: (i, j) — row and column
    Rotated:  (i+j, i-j) — sum and difference

  The rotation maps:
    Row correlation -> diagonal correlation
    Column correlation -> anti-diagonal correlation
    Diagonal correlation -> row correlation!

  This means: ANY algorithm that works well on row-structured data
  can be applied to diagonal-structured data by rotating first.

APPLICATIONS:
  1. TOURNAMENT ADJACENCY: A[i][j] in rotated coords reveals rank-gap structure
  2. ATTENTION MATRICES: rotated attention groups by "attention distance"
  3. TIME SERIES DELAY EMBEDDING: rotated = autocorrelation structure
  4. CORRELATION MATRICES: rotated = sector structure in finance
  5. DISTANCE MATRICES: rotated = metric ball structure
"""

import sys
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  THE 45-DEGREE TRICK: ABSTRACT APPLICATIONS")
print("  kind-pasteur-2026-03-25-S20gq")
print("=" * 80)


def rotate_45(M):
    """Rearrange NxN matrix by (i+j, i-j) coordinates.
    Returns a (2N-1) x (2N-1) sparse representation."""
    N = len(M)
    # In rotated coords: s = i+j ranges 0..2N-2, d = i-j ranges -(N-1)..N-1
    result = {}
    for i in range(N):
        for j in range(N):
            s, d = i+j, i-j
            result[(s, d)] = M[i][j]
    return result, N


def extract_sum_slices(rotated, N):
    """Extract slices at constant s = i+j (anti-diagonals of original)."""
    slices = {}
    for s in range(2*N - 1):
        sl = []
        for d in range(-(N-1), N):
            if (s, d) in rotated:
                sl.append(rotated[(s, d)])
        slices[s] = sl
    return slices


def extract_diff_slices(rotated, N):
    """Extract slices at constant d = i-j (diagonals of original)."""
    slices = {}
    for d in range(-(N-1), N):
        sl = []
        for s in range(2*N - 1):
            if (s, d) in rotated:
                sl.append(rotated[(s, d)])
        slices[d] = sl
    return slices


# ================================================================
print("\n  1. TOURNAMENT ADJACENCY IN ROTATED COORDINATES:")
# A tournament A[i][j] = 1 if i beats j. Antisymmetric: A[i][j] + A[j][i] = 1.
# In rotated coords:
#   s = i+j = "average rank position" of the pair
#   d = i-j = "rank gap" (positive = i has higher index)
# For antisymmetric: A(s,d) + A(s,-d) = 1 (opposite gaps are complementary)

N = 8
np.random.seed(42)
# Random tournament
A = np.zeros((N, N), dtype=int)
for i in range(N):
    for j in range(i+1, N):
        if np.random.random() < 0.5:
            A[i][j] = 1
        else:
            A[j][i] = 1

rotated, _ = rotate_45(A.tolist())
diff_slices = extract_diff_slices(rotated, N)

print(f"  Tournament on {N} vertices:")
print(f"  In standard coords: antisymmetric (A[i][j] + A[j][i] = 1)")
print(f"  In rotated coords: d=i-j slices show GAP-DEPENDENT win rates")
for d in sorted(diff_slices.keys()):
    sl = diff_slices[d]
    if len(sl) > 0:
        avg = sum(sl) / len(sl)
        print(f"    d={d:+3d} (gap {abs(d)}): avg={avg:.3f}, entries={len(sl)}")

# For a near-transitive tournament: avg increases with d
# (higher-index beats lower-index more often)
print(f"  Near-transitive should show avg increasing with d.")

# ================================================================
print(f"\n  2. ATTENTION MATRIX (transformer-style):")
# Self-attention: A[i][j] = softmax(Q_i . K_j / sqrt(d)).
# In rotated coords: s=i+j, d=i-j.
# d = "attention distance" (how far apart in sequence).
# The rotated view groups by attention distance -> reveals
# whether attention is LOCAL (high near d=0) or GLOBAL.

N = 16
# Simulate local attention: A[i][j] ~ exp(-|i-j|/tau)
tau = 3.0
A_att = np.array([[np.exp(-abs(i-j)/tau) for j in range(N)] for i in range(N)])
# Normalize rows (softmax-like)
A_att = A_att / A_att.sum(axis=1, keepdims=True)

rotated, _ = rotate_45(A_att.tolist())
diff_slices = extract_diff_slices(rotated, N)

print(f"  Local attention (tau={tau}) on {N} tokens:")
print(f"  In rotated coords: d=i-j = attention distance")
for d in range(-4, 5):
    sl = diff_slices.get(d, [])
    if sl:
        avg = sum(sl) / len(sl)
        print(f"    d={d:+3d}: avg_attention={avg:.4f}")

print(f"  -> Attention decays with |d| (local attention = band-limited in rotated)")
print(f"  -> Compressible: keep only |d| <= bandwidth, zero the rest")

# How many entries needed for 95% of total attention?
total = sum(sum(sl) for sl in diff_slices.values())
cumul = 0
for bandwidth in range(N):
    for d in range(-bandwidth, bandwidth+1):
        if d in diff_slices:
            cumul += sum(diff_slices[d])
    if cumul / total >= 0.95:
        entries = sum(len(diff_slices.get(d, [])) for d in range(-bandwidth, bandwidth+1))
        print(f"  95% attention captured by |d| <= {bandwidth}: {entries}/{N*N} entries ({entries/(N*N)*100:.0f}%)")
        break

# ================================================================
print(f"\n  3. TIME SERIES DELAY EMBEDDING:")
# Given sequence x[0], ..., x[N-1]
# Delay matrix: M[i][j] = x[i] * x[j] (outer product = rank-1 embedding)
# In rotated coords: s=i+j = "center time", d=i-j = "lag"
# Constant-d slices = autocorrelation at lag d!

N = 32
# Periodic signal + noise
x = np.array([np.sin(2*np.pi*i/8) + 0.2*np.random.randn() for i in range(N)])
M = np.outer(x, x)

rotated, _ = rotate_45(M.tolist())
diff_slices = extract_diff_slices(rotated, N)

print(f"  Periodic signal (period=8) with noise, N={N}:")
print(f"  Rotated constant-d slices = autocorrelation at lag d")
for d in [0, 1, 2, 4, 8, 12, 16]:
    if d in diff_slices and diff_slices[d]:
        avg = sum(diff_slices[d]) / len(diff_slices[d])
        print(f"    lag {d:3d}: autocorr={avg:.4f}")

print(f"  -> Period visible: autocorr peaks at lag 0, 8, 16 (multiples of period)")

# ================================================================
print(f"\n  4. FINANCIAL CORRELATION MATRIX:")
# Stocks in same sector are correlated. Order stocks by sector.
# In standard coords: block-diagonal structure.
# In rotated coords: the blocks become diagonal bands!

N = 16
# 4 sectors of 4 stocks each
sectors = [0]*4 + [1]*4 + [2]*4 + [3]*4
C = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        if sectors[i] == sectors[j]:
            C[i][j] = 0.8 + 0.2*np.random.rand()
        else:
            C[i][j] = 0.1 + 0.1*np.random.rand()

rotated, _ = rotate_45(C.tolist())
diff_slices = extract_diff_slices(rotated, N)

print(f"  4-sector correlation matrix, {N} stocks:")
print(f"  Standard: block-diagonal (sectors along diagonal)")
print(f"  Rotated: sector structure visible in d-slices")
for d in [0, 1, 2, 3, 4, 5, 8, 12]:
    if d in diff_slices and diff_slices[d]:
        avg = sum(diff_slices[d]) / len(diff_slices[d])
        print(f"    d={d:3d}: avg_corr={avg:.4f} {'(within sector)' if d < 4 else '(cross sector)' if d >= 4 else ''}")

# ================================================================
print(f"\n  5. COMPRESSION BENEFIT (rotated vs standard for each type):")
import zlib

for name, M_data in [
    ("tournament", A.tolist()),
    ("attention", (A_att * 255).astype(int).tolist()),
    ("correlation", (C * 255).astype(int).tolist()),
]:
    data_row = bytes(np.array([M_data[i][j] for i in range(len(M_data)) for j in range(len(M_data[0]))], dtype=np.uint8))
    # Diagonal scan
    N_m = len(M_data)
    data_diag = []
    for d in range(N_m):
        for i in range(d, N_m):
            data_diag.append(M_data[i][i-d])
        for i in range(N_m-d):
            data_diag.append(M_data[i][i+d])
    data_diag = bytes(np.array(data_diag[:N_m*N_m], dtype=np.uint8))

    row_z = len(zlib.compress(data_row, 9))
    diag_z = len(zlib.compress(data_diag, 9))
    winner = "DIAG" if diag_z < row_z else "ROW"
    ratio = row_z / diag_z if diag_z > 0 else 0

    print(f"  {name:>15}: row={row_z}B, diag={diag_z}B, winner={winner}, ratio={ratio:.2f}x")

print(f"\n{'='*60}")
print("THE ABSTRACT PRINCIPLE")
print(f"{'='*60}")
print("""
The 45-degree rotation is a CHANGE OF BASIS that reveals hidden structure:

  STANDARD (i,j):  Row structure visible. Good for row-correlated data.
  ROTATED (i+j, i-j): Diagonal structure visible. Good for distance-correlated data.

APPLICATION MAP:
  Standard basis good for: images (spatial), sequences (temporal), tables (records)
  Rotated basis good for: tournaments (ranking gap), attention (distance),
                           autocorrelation (lag), band-limited matrices

THE UNIVERSAL TRICK: When your data has structure along DIAGONALS,
rotate 45 degrees to make it ROW structure, then apply any row-based tool.

This is why the tournament decomposition works on band-limited matrices:
  Band structure = diagonal correlation = ROW correlation after rotation.
  zlib works on rows. Rotate -> zlib works on bands.

It's also why attention is compressible:
  Local attention = band-limited = few significant diagonal slices.
  In rotated coords: few significant ROWS = sparse representation.
""")

print("DONE.")
print("=" * 80)
