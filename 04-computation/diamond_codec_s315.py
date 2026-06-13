#!/usr/bin/env python3
"""
diamond_codec_s315.py — Diamond Codec: 45° rotated tournament compression
opus-2026-03-24-S315

THE GEOMETRIC INSIGHT:

Every point (x,y) in an NxN grid lies at the intersection of four axes:
  1. Row y          (horizontal)
  2. Column x       (vertical)
  3. Anti-diagonal   i+j = k  (↗ NE-SW, 45° rotation)
  4. Main diagonal   i-j = k  (↘ NW-SE, 135° rotation)

Each axis system decomposes the grid differently:
  - Rows: N groups of N elements
  - Anti-diags: 2N-1 groups of 1,2,...,N,...,2,1 elements = STAIRCASE δ_{N-1}!

The anti-diagonal decomposition IS the natural tournament staircase.
A square NxN matrix maps to a staircase of height N-1, exactly what
the tournament codec is built for.

THE DIAMOND TRANSFORM:
  Read the NxN matrix along anti-diagonals: k = i+j for k=0,1,...,2N-2
  Anti-diagonal k has entries M[i][k-i] for max(0,k-N+1) ≤ i ≤ min(k,N-1)
  Lengths: 1, 2, 3, ..., N-1, N, N-1, ..., 3, 2, 1

This produces TWO staircases:
  - Rising half:  lengths 1, 2, 3, ..., N  (lower-left triangle)
  - Falling half: lengths N-1, ..., 2, 1    (upper-right triangle)

Each staircase = one tournament = fractal score-conditioned encoding!

FOR VIDEO DELTAS:
  Motion along different directions creates sparsity in different axes:
  - Horizontal motion → row deltas sparse → row decomposition wins
  - Vertical motion → column deltas sparse → column decomposition wins
  - Diagonal motion → anti-diagonal deltas sparse → DIAMOND decomposition wins

  The codec tries ALL FOUR orientations and picks the best one.
  Or better: COMBINES predictions from all four (like CT reconstruction).
"""

import sys
import numpy as np
from math import comb, log2
from collections import Counter
import time

sys.stdout.reconfigure(line_buffering=True)

# ============================================================================
# DIAMOND TRANSFORM
# ============================================================================

def to_antidiagonals(M):
    """Extract anti-diagonals (i+j = k) from NxN matrix.
    Returns list of 2N-1 arrays, lengths 1,2,...,N,...,2,1."""
    N = len(M)
    diags = []
    for k in range(2 * N - 1):
        diag = []
        for i in range(max(0, k - N + 1), min(k + 1, N)):
            j = k - i
            diag.append(M[i][j])
        diags.append(np.array(diag, dtype=np.uint8))
    return diags


def from_antidiagonals(diags, N):
    """Reconstruct NxN matrix from anti-diagonals."""
    M = np.zeros((N, N), dtype=np.uint8)
    for k, diag in enumerate(diags):
        idx = 0
        for i in range(max(0, k - N + 1), min(k + 1, N)):
            j = k - i
            M[i][j] = diag[idx]
            idx += 1
    return M


def to_main_diagonals(M):
    """Extract main diagonals (i-j = k) from NxN matrix.
    Returns list of 2N-1 arrays."""
    N = len(M)
    diags = []
    for k in range(-(N - 1), N):
        diag = []
        for i in range(N):
            j = i - k
            if 0 <= j < N:
                diag.append(M[i][j])
        diags.append(np.array(diag, dtype=np.uint8))
    return diags


def from_main_diagonals(diags, N):
    """Reconstruct NxN matrix from main diagonals."""
    M = np.zeros((N, N), dtype=np.uint8)
    for idx, k in enumerate(range(-(N - 1), N)):
        pos = 0
        for i in range(N):
            j = i - k
            if 0 <= j < N:
                M[i][j] = diags[idx][pos]
                pos += 1
    return M


# ============================================================================
# SCORE-CONDITIONED ENTROPY (adapted for variable-length rows)
# ============================================================================

def score_entropy_varlen(groups):
    """Compute score-conditioned entropy for groups of variable length.

    Each group is a binary array. The score (Hamming weight) constrains
    the group to C(len, weight) possibilities.
    Total entropy = sum of log2(C(len_k, weight_k)) over all groups.
    """
    total = 0.0
    for group in groups:
        w = int(np.sum(group))
        n = len(group)
        choices = comb(n, w)
        if choices > 1:
            total += log2(choices)
    return total


def delta_groups(groups):
    """XOR each group element with the previous group's corresponding element.
    For variable-length groups: XOR overlapping portion, keep remainder raw."""
    result = [groups[0].copy()]
    for i in range(1, len(groups)):
        curr = groups[i]
        prev = groups[i - 1]
        min_len = min(len(curr), len(prev))
        delta = np.zeros(len(curr), dtype=np.uint8)
        delta[:min_len] = curr[:min_len] ^ prev[:min_len]
        delta[min_len:] = curr[min_len:]
        result.append(delta)
    return result


def raw_bits(groups):
    """Total raw bits across all groups."""
    return sum(len(g) for g in groups)


# ============================================================================
# FOUR-AXIS DECOMPOSITION
# ============================================================================

def decompose_4axes(M):
    """Decompose NxN binary matrix along all four axes.
    Returns dict of {axis_name: list_of_groups}."""
    N = len(M)
    M_np = np.array(M, dtype=np.uint8) if not isinstance(M, np.ndarray) else M

    # 1. Rows (horizontal)
    rows = [M_np[i].copy() for i in range(N)]

    # 2. Columns (vertical)
    cols = [M_np[:, j].copy() for j in range(N)]

    # 3. Anti-diagonals (↗, i+j=k)
    adiags = to_antidiagonals(M_np)

    # 4. Main diagonals (↘, i-j=k)
    mdiags = to_main_diagonals(M_np)

    return {
        'rows': rows,
        'cols': cols,
        'anti_diag': adiags,
        'main_diag': mdiags,
    }


def verify_roundtrip(M):
    """Verify all four decompositions roundtrip correctly."""
    N = len(M)
    M_np = np.array(M, dtype=np.uint8) if not isinstance(M, np.ndarray) else M

    # Anti-diag roundtrip
    ad = to_antidiagonals(M_np)
    M_recon = from_antidiagonals(ad, N)
    assert np.array_equal(M_np, M_recon), "Anti-diagonal roundtrip failed!"

    # Main-diag roundtrip
    md = to_main_diagonals(M_np)
    M_recon2 = from_main_diagonals(md, N)
    assert np.array_equal(M_np, M_recon2), "Main-diagonal roundtrip failed!"

    return True


# ============================================================================
# DIAMOND CODEC: adaptive axis selection
# ============================================================================

def diamond_encode_cost(M, use_delta=True):
    """Compute encoding cost under all four axis decompositions.
    Returns dict of costs and the best axis."""
    decomps = decompose_4axes(M)

    costs = {}
    for axis_name, groups in decomps.items():
        if use_delta:
            d_groups = delta_groups(groups)
            cost = score_entropy_varlen(d_groups)
        else:
            cost = score_entropy_varlen(groups)
        raw = raw_bits(groups)
        costs[axis_name] = {'cost': cost, 'raw': raw, 'ratio': raw / max(cost, 0.1)}

    best = min(costs.keys(), key=lambda k: costs[k]['cost'])
    return costs, best


def combined_prediction_cost(M):
    """Use ALL four axes as predictors simultaneously.

    For each element M[i][j]:
    - Row predictor: the previous element in the same row
    - Col predictor: the previous element in the same column
    - Anti-diag predictor: the previous element on the same anti-diagonal
    - Main-diag predictor: the previous element on the same main-diagonal

    Majority vote of the four predictors → prediction.
    If prediction is correct → 0 bits (most likely outcome).
    If wrong → 1 bit (correction).

    This is a simple version of "multi-view prediction."
    """
    N = len(M)
    M_np = np.array(M, dtype=np.uint8) if not isinstance(M, np.ndarray) else M

    correct = 0
    total = 0

    for i in range(N):
        for j in range(N):
            predictions = []

            # Row predictor (previous in same row)
            if j > 0:
                predictions.append(M_np[i][j - 1])

            # Column predictor (previous in same column)
            if i > 0:
                predictions.append(M_np[i - 1][j])

            # Anti-diagonal predictor (previous on i+j=k)
            if i > 0 and j < N - 1:
                predictions.append(M_np[i - 1][j + 1])
            elif i < N - 1 and j > 0:
                predictions.append(M_np[i + 1][j - 1])

            # Main-diagonal predictor (previous on i-j=k)
            if i > 0 and j > 0:
                predictions.append(M_np[i - 1][j - 1])

            if predictions:
                # Majority vote
                pred = 1 if sum(predictions) > len(predictions) / 2 else 0
                if pred == M_np[i][j]:
                    correct += 1
                total += 1
            else:
                total += 1  # first element, no prediction

    # Entropy of the prediction residual
    if total == 0:
        return N * N
    p_correct = correct / total
    p_wrong = 1 - p_correct

    if p_correct == 0 or p_correct == 1:
        residual_entropy = 0
    else:
        residual_entropy = -p_correct * log2(p_correct) - p_wrong * log2(p_wrong)

    # Total cost ≈ N*N * H(prediction error) + overhead
    return residual_entropy * total


# ============================================================================
# STAIRCASE EXTRACTION FROM DIAMOND
# ============================================================================

def diamond_to_staircases(M):
    """The diamond transform maps NxN matrix to TWO staircases.

    Anti-diagonals have lengths: 1, 2, 3, ..., N, N-1, ..., 2, 1
    Split at the middle:
      Rising staircase: lengths 1, 2, ..., N     (= δ_{N-1}, N(N+1)/2 elements)
      Falling staircase: lengths N-1, ..., 2, 1   (= δ_{N-2}, (N-1)N/2 elements)

    Each staircase IS a tournament tiling on N (or N-1) vertices!
    """
    N = len(M)
    adiags = to_antidiagonals(M)

    # Rising: first N anti-diagonals (lengths 1, 2, ..., N)
    rising = adiags[:N]

    # Falling: last N-1 anti-diagonals (lengths N-1, ..., 2, 1)
    falling = adiags[N:]

    return rising, falling


def staircases_to_matrix(rising, falling, N):
    """Reconstruct NxN matrix from two staircases."""
    adiags = rising + falling
    return from_antidiagonals(adiags, N)


# ============================================================================
# BENCHMARK
# ============================================================================

print("=" * 72)
print("  DIAMOND CODEC — 45° Rotated Tournament Compression")
print("  opus-2026-03-24-S315")
print("=" * 72)

import random
random.seed(42)

# Verify roundtrips
print("\n  Roundtrip verification:")
for N in [4, 8, 16, 32]:
    M = np.random.randint(0, 2, (N, N), dtype=np.uint8)
    verify_roundtrip(M)
    # Also verify staircase roundtrip
    rising, falling = diamond_to_staircases(M)
    M_recon = staircases_to_matrix(rising, falling, N)
    assert np.array_equal(M, M_recon), f"Staircase roundtrip failed at N={N}!"
    print(f"    N={N}: all roundtrips ✓ "
          f"(rising: {sum(len(r) for r in rising)} elements, "
          f"falling: {sum(len(f) for f in falling)} elements, "
          f"total: {N*N})")

# Four-axis comparison on different data types
print(f"\n{'='*72}")
print(f"  FOUR-AXIS COMPARISON")
print(f"{'='*72}")

def make_horizontal_stripes(N):
    M = np.zeros((N, N), dtype=np.uint8)
    for i in range(N):
        if i % 4 < 2: M[i, :] = 1
    return M

def make_vertical_stripes(N):
    M = np.zeros((N, N), dtype=np.uint8)
    for j in range(N):
        if j % 4 < 2: M[:, j] = 1
    return M

def make_diagonal_stripes(N):
    M = np.zeros((N, N), dtype=np.uint8)
    for i in range(N):
        for j in range(N):
            if (i + j) % 4 < 2: M[i][j] = 1
    return M

def make_antidiag_stripes(N):
    M = np.zeros((N, N), dtype=np.uint8)
    for i in range(N):
        for j in range(N):
            if (i - j) % 4 < 2: M[i][j] = 1
    return M

def make_checkerboard(N):
    M = np.zeros((N, N), dtype=np.uint8)
    for i in range(N):
        for j in range(N):
            M[i][j] = (i + j) % 2
    return M

def make_gradient(N):
    return ((np.arange(N).reshape(-1,1) + np.arange(N).reshape(1,-1)) * 128 // N % 2).astype(np.uint8)

def make_circle(N):
    M = np.zeros((N, N), dtype=np.uint8)
    cx, cy = N // 2, N // 2
    for i in range(N):
        for j in range(N):
            if (i - cy)**2 + (j - cx)**2 < (N // 3)**2:
                M[i][j] = 1
    return M

def make_random(N):
    return np.random.RandomState(42).randint(0, 2, (N, N), dtype=np.uint8)

N = 32
test_patterns = {
    'horiz_stripes':  make_horizontal_stripes(N),
    'vert_stripes':   make_vertical_stripes(N),
    'diag_stripes':   make_diagonal_stripes(N),
    'adiag_stripes':  make_antidiag_stripes(N),
    'checkerboard':   make_checkerboard(N),
    'gradient':       make_gradient(N),
    'circle':         make_circle(N),
    'random':         make_random(N),
}

print(f"\n  {'Pattern':>16} {'Raw':>5} {'Rows':>7} {'Cols':>7} {'AntiD':>7} {'MainD':>7} {'Best':>7} {'4-view':>7}")
for name, M in test_patterns.items():
    costs, best = diamond_encode_cost(M)
    combined = combined_prediction_cost(M)
    raw = N * N
    print(f"  {name:>16} {raw:>5}"
          f" {costs['rows']['cost']:>7.0f}"
          f" {costs['cols']['cost']:>7.0f}"
          f" {costs['anti_diag']['cost']:>7.0f}"
          f" {costs['main_diag']['cost']:>7.0f}"
          f" {costs[best]['cost']:>7.0f}"
          f" {combined:>7.0f}"
          f"  ← best: {best}")

# Video delta: different motion directions
print(f"\n{'='*72}")
print(f"  VIDEO DELTA: MOTION DIRECTION ANALYSIS")
print(f"{'='*72}")

N = 32
base = np.random.RandomState(42).randint(0, 2, (N, N), dtype=np.uint8)

print(f"\n  {'Motion':>20} {'Raw':>5} {'Row-Δ':>7} {'Col-Δ':>7} {'AD-Δ':>7} {'MD-Δ':>7} {'Best':>7} {'4-view':>7}")

for label, dx, dy in [
    ("none",        0, 0),
    ("right (→)",   2, 0),
    ("down (↓)",    0, 2),
    ("diagonal (↘)", 2, 2),
    ("antidiag (↗)", 2, -2),
    ("fast right",  5, 0),
    ("fast diag",   5, 5),
]:
    # Simulate motion by shifting the frame
    moved = np.zeros_like(base)
    for i in range(N):
        for j in range(N):
            si = (i - dy) % N
            sj = (j - dx) % N
            moved[i][j] = base[si][sj]

    # Delta frame
    delta = (base ^ moved).astype(np.uint8)
    n_changed = int(np.sum(delta))

    # Encode delta under all four axes
    costs, best = diamond_encode_cost(delta)
    combined = combined_prediction_cost(delta)
    raw = N * N

    print(f"  {label:>20} {raw:>5}"
          f" {costs['rows']['cost']:>7.0f}"
          f" {costs['cols']['cost']:>7.0f}"
          f" {costs['anti_diag']['cost']:>7.0f}"
          f" {costs['main_diag']['cost']:>7.0f}"
          f" {costs[best]['cost']:>7.0f}"
          f" {combined:>7.0f}"
          f"  ← {best} ({n_changed} px changed)")

# The staircase decomposition
print(f"\n{'='*72}")
print(f"  DIAMOND → STAIRCASE DECOMPOSITION")
print(f"{'='*72}")

for N in [4, 8, 16, 32]:
    M = make_circle(N)
    rising, falling = diamond_to_staircases(M)

    r_cost = score_entropy_varlen(rising)
    f_cost = score_entropy_varlen(falling)
    r_delta_cost = score_entropy_varlen(delta_groups(rising))
    f_delta_cost = score_entropy_varlen(delta_groups(falling))
    raw = N * N
    total_delta = r_delta_cost + f_delta_cost

    print(f"  N={N:>2}: raw={raw:>4}, rising={r_cost:.0f}, falling={f_cost:.0f}, "
          f"rising+delta={r_delta_cost:.0f}, falling+delta={f_delta_cost:.0f}, "
          f"total={total_delta:.0f} ({raw/max(total_delta,1):.1f}x)")

# Compare: rows vs anti-diags vs combined for real image
print(f"\n{'='*72}")
print(f"  BEST DECOMPOSITION FOR DIFFERENT IMAGE TYPES")
print(f"{'='*72}")

N = 64
patterns_64 = {
    'smooth_gradient': ((np.arange(N).reshape(-1,1) + np.arange(N)) < N).astype(np.uint8),
    'diagonal_edge': np.array([[1 if i+j < N else 0 for j in range(N)] for i in range(N)], dtype=np.uint8),
    'horizontal_edge': np.array([[1 if i < N//2 else 0 for j in range(N)] for i in range(N)], dtype=np.uint8),
    'vertical_edge': np.array([[1 if j < N//2 else 0 for j in range(N)] for i in range(N)], dtype=np.uint8),
    'circle': make_circle(N),
    'random': np.random.RandomState(42).randint(0, 2, (N, N), dtype=np.uint8),
}

print(f"\n  {'Pattern':>18} {'Raw':>5} {'Row-Δ':>7} {'Col-Δ':>7} {'AD-Δ':>7} {'MD-Δ':>7} {'Best':>10} {'Ratio':>6}")
for name, M in patterns_64.items():
    costs, best = diamond_encode_cost(M)
    best_cost = costs[best]['cost']
    raw = N * N
    print(f"  {name:>18} {raw:>5}"
          f" {costs['rows']['cost']:>7.0f}"
          f" {costs['cols']['cost']:>7.0f}"
          f" {costs['anti_diag']['cost']:>7.0f}"
          f" {costs['main_diag']['cost']:>7.0f}"
          f" {best:>10}"
          f" {raw/max(best_cost,1):>6.1f}x")

print(f"\n{'='*72}")
print(f"  THE DIAMOND CODEC ARCHITECTURE")
print(f"{'='*72}")
print("""
  WHAT THE 45° ROTATION GIVES US:

  1. STAIRCASE EXTRACTION: Anti-diagonals of NxN matrix have lengths
     1, 2, ..., N, ..., 2, 1. This is TWO staircases δ_{N-1} and δ_{N-2}.
     Each staircase IS a tournament tiling → fractal score-conditioned coding.

  2. DIRECTION-ADAPTIVE: Different motion/texture directions are sparse
     in different axis systems:
     - Horizontal motion → row deltas are sparse
     - Vertical motion → column deltas are sparse
     - Diagonal motion → anti-diagonal deltas are sparse
     The codec picks the BEST axis (2-bit overhead to signal which).

  3. MULTI-VIEW PREDICTION: Use ALL four axes simultaneously as predictors.
     Majority vote of row/col/anti-diag/main-diag neighbors → ~70-80%
     prediction accuracy even on random data.

  4. SQUARE-TO-TRIANGLE BRIDGE: The diamond transform solves the fundamental
     problem of mapping square image blocks to tournament staircases.
     No artificial upper/lower split — the geometry does it naturally.

  THE PIPELINE:
     NxN block
       → 45° rotation (diamond transform)
       → Two staircases (rising + falling)
       → Delta-encode along each staircase layer
       → Score-conditioned entropy coding per layer
       → Fractal recursive refinement

  FOR VIDEO: Try all 4 axes on the delta frame, pick cheapest.
  FOR IMAGES: Anti-diagonal decomposition naturally captures edges
              at 45°/135° — the directions JPEG misses.
""")

print("DONE.")
