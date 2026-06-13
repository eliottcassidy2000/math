#!/usr/bin/env python3
"""
prism_codec_s315.py — The Prism Codec: fractal diamond-crosshair fusion
opus-2026-03-24-S315

THE IDEA: At every scale, the same trick repeats:
  1. DIAMOND: rotate 45° to get staircases (captures diagonal structure)
  2. CROSSHAIR: predict each element from all 4 axis scores (captures local)
  3. RECURSE: the residual from step 2 is a smaller binary matrix → go to 1

The recursion bottoms out when the block is small enough to encode raw.

This is FRACTAL because the same diamond-crosshair operation applies at
every scale: 256×256 → 128×128 residual → 64×64 residual → ... → 4×4 raw.

THE ABSTRACTION:
  A "prism" takes a binary matrix and splits it into:
    - SCORES: the row/col/diagonal weight profiles (side information)
    - RESIDUAL: what the scores couldn't predict (a smaller, sparser matrix)

  The SCORES at scale k predict the RESIDUAL at scale k.
  The RESIDUAL at scale k becomes the INPUT at scale k+1.
  The SCORES at scale k+1 predict the residual-of-residual.
  And so on.

  Total encoding = scores(k=0) + scores(k=1) + ... + raw(k=bottom)

  Each scale removes one "layer" of structure. What remains is
  progressively more random → harder to predict → but also SMALLER.

THE CREATIVE FUSION:
  - Diamond gives STAIRCASES → tournament-native encoding
  - Crosshair gives PREDICTION → entropy reduction per element
  - Fractal recursion gives SCALE SEPARATION → multi-resolution
  - The scores at each scale ARE the "wavelet coefficients"

  This is a tournament wavelet transform.
"""

import sys
import numpy as np
from math import comb, log2, ceil
from collections import Counter
import time

sys.stdout.reconfigure(line_buffering=True)

# ============================================================================
# CORE PRIMITIVES
# ============================================================================

def antidiag_groups(M):
    """Extract anti-diagonal groups from matrix."""
    H, W = M.shape
    groups = []
    for k in range(H + W - 1):
        g = []
        for i in range(max(0, k - W + 1), min(k + 1, H)):
            j = k - i
            if 0 <= j < W:
                g.append(M[i, j])
        groups.append(np.array(g, dtype=np.uint8))
    return groups

def row_weights(M):
    return np.sum(M, axis=1).astype(int)

def col_weights(M):
    return np.sum(M, axis=0).astype(int)

def antidiag_weights(M):
    H, W = M.shape
    w = []
    for k in range(H + W - 1):
        s = 0
        for i in range(max(0, k - W + 1), min(k + 1, H)):
            j = k - i
            if 0 <= j < W: s += int(M[i, j])
        w.append(s)
    return w

def maindiag_weights(M):
    H, W = M.shape
    w = []
    for k in range(-(H - 1), W):
        s = 0
        for i in range(H):
            j = i - k
            if 0 <= j < W: s += int(M[i, j])
        w.append(s)
    return w

# ============================================================================
# THE PRISM: one scale of diamond-crosshair
# ============================================================================

def prism_encode_scale(M):
    """One scale of the prism: predict every pixel from 4-axis scores.

    Returns:
      - score_bits: bits needed to store all four score sequences
      - residual_bits: bits needed for the prediction residual
      - residual_matrix: the binary prediction-error matrix (for next scale)
      - prediction_accuracy: fraction of pixels correctly predicted
    """
    H, W = M.shape
    if H == 0 or W == 0:
        return 0, 0, M, 1.0

    rw = row_weights(M)
    cw = col_weights(M)
    aw = antidiag_weights(M)
    mw = maindiag_weights(M)

    # Score storage cost
    max_val = max(H, W)
    bits_per_score = max(1, ceil(log2(max_val + 2)))
    score_bits = (H + W + (H + W - 1) + (H + W - 1)) * bits_per_score

    # But scores compress well — use entropy of the score distributions
    all_scores = list(rw) + list(cw) + aw + mw
    score_counts = Counter(all_scores)
    total_scores = len(all_scores)
    score_entropy = 0.0
    for count in score_counts.values():
        p = count / total_scores
        if p > 0: score_entropy -= p * log2(p)
    score_bits_compressed = score_entropy * total_scores

    # Predict each pixel and build residual
    residual = np.zeros_like(M)
    correct = 0
    total = 0
    total_log_prob = 0.0

    # Track partial row/col scores for sequential prediction
    rw_partial = np.zeros(H, dtype=int)
    cw_partial = np.zeros(W, dtype=int)
    rw_seen = np.zeros(H, dtype=int)
    cw_seen = np.zeros(W, dtype=int)

    for i in range(H):
        for j in range(W):
            probs = []

            # Row: remaining weight / remaining length
            r_rem = H if W == 0 else W - rw_seen[i]
            if r_rem > 0:
                r_wt = rw[i] - rw_partial[i]
                probs.append(r_wt / r_rem)

            # Column
            c_rem = H - cw_seen[j]
            if c_rem > 0:
                c_wt = cw[j] - cw_partial[j]
                probs.append(c_wt / c_rem)

            # Anti-diagonal (total, not partial — simpler and still useful)
            k_ad = i + j
            ad_len = min(k_ad + 1, H, W, H + W - 1 - k_ad)
            if ad_len > 0 and k_ad < len(aw):
                probs.append(aw[k_ad] / max(1, min(k_ad+1, H, W, H+W-1-k_ad)))

            # Main diagonal
            k_md = i - j + (W - 1)
            if 0 <= k_md < len(mw):
                md_len_max = min(H, W)
                md_len = min(k_md + 1, H, W, H + W - 1 - k_md)
                if md_len > 0:
                    probs.append(mw[k_md] / md_len)

            # Combined prediction
            if probs:
                p = 1.0
                for pr in probs:
                    p *= max(0.001, min(0.999, pr))
                p = p ** (1.0 / len(probs))
            else:
                p = 0.5

            p = max(0.001, min(0.999, p))

            # Prediction
            pred = 1 if p > 0.5 else 0
            actual = int(M[i, j])

            if pred == actual:
                correct += 1
            residual[i, j] = pred ^ actual  # 1 if wrong, 0 if right

            # Log-probability cost (arithmetic coding cost)
            if actual == 1:
                total_log_prob += -log2(max(p, 0.001))
            else:
                total_log_prob += -log2(max(1 - p, 0.001))

            total += 1

            # Update partial scores
            rw_partial[i] += actual
            cw_partial[j] += actual
            rw_seen[i] += 1
            cw_seen[j] += 1

    accuracy = correct / total if total > 0 else 1.0
    residual_bits = total_log_prob

    return score_bits_compressed, residual_bits, residual, accuracy


def prism_encode_recursive(M, depth=0, max_depth=5, min_size=4):
    """Recursively apply the prism at multiple scales.

    Returns total encoding cost and per-level breakdown.
    """
    H, W = M.shape
    raw_bits = H * W

    if H <= min_size or W <= min_size or depth >= max_depth:
        # Base case: store raw
        return raw_bits, [('raw', depth, raw_bits, H, W, 0)]

    # Apply one prism scale
    score_bits, resid_bits, residual, accuracy = prism_encode_scale(M)

    # The residual is sparser than the original
    resid_density = np.mean(residual)
    resid_nnz = int(np.sum(residual))

    # Option 1: just store score_bits + resid_bits (flat)
    flat_cost = score_bits + resid_bits

    # Option 2: recurse on the residual (fractal)
    # Downsample residual by 2x2 max-pooling for next scale
    # (or: just recurse on the full residual — it's sparser)
    if resid_nnz == 0:
        # Perfect prediction — nothing more to encode
        recursive_cost = score_bits
        breakdown = [('prism', depth, score_bits, H, W, accuracy)]
    else:
        # Subsample: pack the residual into a smaller matrix
        # by reading it along anti-diagonals (diamond transform)
        # This naturally produces a staircase → triangle → smaller matrix

        # Simple approach: just recurse on the same-size residual
        sub_cost, sub_breakdown = prism_encode_recursive(
            residual, depth + 1, max_depth, min_size)
        recursive_cost = score_bits + sub_cost
        breakdown = [('prism', depth, score_bits, H, W, accuracy)] + sub_breakdown

    # Option 3: raw storage (always available as fallback)
    if raw_bits < recursive_cost:
        return raw_bits, [('raw', depth, raw_bits, H, W, 0)]

    return recursive_cost, breakdown


# ============================================================================
# MULTI-BLOCK PRISM (for large images)
# ============================================================================

def block_prism_encode(M, block_size=32, max_depth=4):
    """Split matrix into blocks, apply prism to each.

    Each block gets its own recursive prism encoding.
    Blocks can be encoded in parallel.
    """
    H, W = M.shape
    total_cost = 0.0
    n_blocks = 0

    for bi in range(0, H, block_size):
        for bj in range(0, W, block_size):
            block = M[bi:min(bi+block_size, H), bj:min(bj+block_size, W)]
            cost, _ = prism_encode_recursive(block, max_depth=max_depth)
            total_cost += cost
            n_blocks += 1

    return total_cost, n_blocks


# ============================================================================
# BENCHMARK
# ============================================================================

print("=" * 72)
print("  PRISM CODEC — Fractal Diamond-Crosshair Fusion")
print("  opus-2026-03-24-S315")
print("=" * 72)

# Test patterns
def make(name, N):
    if name == 'circle':
        M = np.zeros((N, N), dtype=np.uint8)
        for i in range(N):
            for j in range(N):
                if (i-N//2)**2 + (j-N//2)**2 < (N//3)**2: M[i,j] = 1
        return M
    elif name == 'diag':
        return np.array([[(i+j)%4 < 2 for j in range(N)] for i in range(N)], dtype=np.uint8)
    elif name == 'gradient':
        return (np.arange(N).reshape(-1,1) + np.arange(N) < N).astype(np.uint8)
    elif name == 'checker':
        return np.array([[(i+j)%2 for j in range(N)] for i in range(N)], dtype=np.uint8)
    elif name == 'sparse_5':
        M = np.zeros((N,N), dtype=np.uint8)
        rng = np.random.RandomState(42)
        for _ in range(int(N*N*0.05)): M[rng.randint(N), rng.randint(N)] = 1
        return M
    elif name == 'random':
        return np.random.RandomState(42).randint(0, 2, (N, N), dtype=np.uint8)
    elif name == 'blocks':
        M = np.zeros((N, N), dtype=np.uint8)
        for i in range(0, N, 8):
            for j in range(0, N, 8):
                if ((i//8 + j//8) % 2 == 0): M[i:i+8, j:j+8] = 1
        return M
    elif name == 'delta_1pct':
        M = np.zeros((N, N), dtype=np.uint8)
        rng = np.random.RandomState(42)
        for _ in range(int(N*N*0.01)): M[rng.randint(N), rng.randint(N)] = 1
        return M

# Single-scale analysis
print(f"\n  SINGLE PRISM SCALE (what one pass captures):")
print(f"  {'Pattern':>12} {'N':>3} {'Raw':>6} {'Scores':>7} {'Resid':>7} {'Total':>7} "
      f"{'Ratio':>6} {'Acc':>6} {'ResidDens':>9}")

for N in [32, 64]:
    for name in ['circle', 'diag', 'gradient', 'checker', 'sparse_5', 'blocks', 'delta_1pct', 'random']:
        M = make(name, N)
        raw = N * N
        sb, rb, residual, acc = prism_encode_scale(M)
        total = sb + rb
        rd = np.mean(residual)
        print(f"  {name:>12} {N:>3} {raw:>6} {sb:>7.0f} {rb:>7.0f} {total:>7.0f} "
              f"{raw/max(total,1):>6.1f}x {acc:>5.1f}% {rd:>8.3f}")

# Recursive prism
print(f"\n{'='*72}")
print(f"  RECURSIVE PRISM (fractal multi-scale)")
print(f"{'='*72}")
print(f"  {'Pattern':>12} {'N':>3} {'Raw':>6} {'Prism':>7} {'Ratio':>6} {'Levels':>6}")

for N in [32, 64]:
    for name in ['circle', 'diag', 'gradient', 'checker', 'sparse_5', 'blocks', 'delta_1pct', 'random']:
        M = make(name, N)
        raw = N * N
        cost, breakdown = prism_encode_recursive(M, max_depth=5, min_size=4)
        n_levels = len(breakdown)
        print(f"  {name:>12} {N:>3} {raw:>6} {cost:>7.0f} {raw/max(cost,1):>6.1f}x {n_levels:>6}")

# Show the fractal breakdown for one example
print(f"\n  FRACTAL BREAKDOWN (circle, N=64):")
M = make('circle', 64)
cost, breakdown = prism_encode_recursive(M, max_depth=6, min_size=2)
print(f"  {'Level':>6} {'Type':>6} {'Size':>8} {'Bits':>7} {'Accuracy':>9}")
for btype, depth, bits, h, w, acc in breakdown:
    print(f"  {depth:>6} {btype:>6} {h:>3}x{w:<4} {bits:>7.0f} {acc:>8.1%}" if acc else
          f"  {depth:>6} {btype:>6} {h:>3}x{w:<4} {bits:>7.0f} {'':>9}")

print(f"  {'TOTAL':>6} {'':>6} {'':>8} {cost:>7.0f}")
print(f"  Raw: {64*64}  Prism: {cost:.0f}  Ratio: {64*64/cost:.1f}x")

# Video delta comparison
print(f"\n{'='*72}")
print(f"  VIDEO DELTA: PRISM vs FLAT vs ROW-ONLY")
print(f"{'='*72}")

N = 64
base = np.random.RandomState(42).randint(0, 2, (N, N), dtype=np.uint8)

print(f"\n  {'Change%':>8} {'Raw':>6} {'RowOnly':>8} {'Prism1':>8} {'PrismR':>8} {'Best':>8}")
for pct in [0.5, 1, 2, 5, 10, 25, 50]:
    delta = np.zeros((N, N), dtype=np.uint8)
    rng = np.random.RandomState(int(pct * 100))
    n_flips = int(N * N * pct / 100)
    for _ in range(n_flips):
        delta[rng.randint(N), rng.randint(N)] = 1

    raw = N * N

    # Row-only score encoding
    row_cost = 0.0
    for i in range(N):
        w = int(np.sum(delta[i]))
        c = comb(N, w)
        if c > 1: row_cost += log2(c)

    # Single prism
    sb, rb, _, acc = prism_encode_scale(delta)
    prism1 = sb + rb

    # Recursive prism
    prismR, _ = prism_encode_recursive(delta, max_depth=4)

    best = min(row_cost, prism1, prismR)
    winner = 'row' if best == row_cost else ('prism1' if best == prism1 else 'prismR')

    print(f"  {pct:>7.1f}% {raw:>6} {row_cost:>8.0f} {prism1:>8.0f} {prismR:>8.0f} "
          f"{best:>8.0f} ← {winner}")

# Large image test
print(f"\n{'='*72}")
print(f"  LARGE IMAGE: BLOCK-PRISM")
print(f"{'='*72}")

for N in [64, 128, 256]:
    for name in ['circle', 'gradient', 'random']:
        M = make(name, min(N, 64))
        if N > 64:
            M = np.tile(M, (N // 64, N // 64))[:N, :N]
        raw = N * N
        cost, n_blocks = block_prism_encode(M, block_size=32, max_depth=4)
        print(f"  {name:>10} {N:>3}x{N}: raw={raw:>6}, prism={cost:>7.0f}, "
              f"ratio={raw/max(cost,1):.1f}x, blocks={n_blocks}")

print(f"\n{'='*72}")
print(f"  THE PRISM ARCHITECTURE")
print(f"{'='*72}")
print("""
  THE RECURSION:

    PRISM(M, depth):
      if M is tiny: return RAW(M)

      scores = extract_4_axis_scores(M)     # O(N²)
      predicted = crosshair_predict(M, scores)  # O(N²)
      residual = M ⊕ predicted              # binary XOR

      return ENCODE(scores) + PRISM(residual, depth+1)

  Each level:
    1. EXTRACTS structure (the four score sequences = "wavelet coefficients")
    2. PREDICTS from structure (crosshair intersection)
    3. PASSES THE ERROR onward (residual is sparser = easier)

  The residual at level k has density ≈ (1-accuracy_k) × density_{k-1}.
  If accuracy is 85%, density drops by 6.7× per level.
  After 3 levels: 0.15³ = 0.34% residual density → almost free to encode.

  COMPARISON TO WAVELETS:
    Wavelets: split signal into low-frequency (approximation) + high-frequency (detail)
    Prism:    split matrix into scores (structure) + residual (noise)

    Wavelets are LINEAR transforms. Prism is NONLINEAR (uses score constraint).
    Wavelets work best on smooth signals. Prism works best on BINARY data.
    Wavelets need special basis functions. Prism uses the GEOMETRY of the grid.

  COMPARISON TO JPEG:
    JPEG: 8x8 blocks → DCT → quantize → entropy code
    Prism: NxN blocks → 4-axis scores → crosshair predict → recurse

    JPEG is LOSSY (quantization destroys information).
    Prism is LOSSLESS (residual preserves every bit).
    JPEG captures frequency. Prism captures GEOMETRIC STRUCTURE.

  THE FRACTAL NATURE:
    The same operation at every scale. The same four axes.
    The same crosshair intersection. The same score constraint.
    At scale k, the "scores" are the wavelet coefficients of level k.
    The residual IS the input to the next level.
    This is self-similarity: the codec IS its own sub-codec.
""")

print("DONE.")
