#!/usr/bin/env python3
"""
tournament_wavelet.py -- Tournament Wavelet Transform for Images
kind-pasteur-2026-03-24-S20cq

THE DEEP CONNECTION: Traditional wavelets decompose by FREQUENCY (low/high).
Tournament theory decomposes by SPATIAL SKIP DISTANCE (long-range/short-range).

These are DUAL decompositions:
  - DCT: x-axis = frequency, y-axis = amplitude
  - Tournament: x-axis = spatial distance, y-axis = correlation strength

For NATURAL images, frequency decomposition (DCT/wavelet) works well because
natural images are smooth (low-frequency dominant).

For STRUCTURED images (diagrams, charts, UI, game art, maps), the tournament
wavelet excels because structure is defined by SPATIAL RELATIONSHIPS, not
frequency content.

THE TOURNAMENT WAVELET:
  Level 0: Original image (all skip distances)
  Level 1: Predict even rows/cols from odd (skip-1 interactions)
  Level 2: Predict every-4th from every-2nd (skip-2 interactions)
  Level k: Predict skip-2^k from skip-2^(k-1)

Each level captures the "detail at that spatial scale." The residual at each
level is the TOURNAMENT COEFFICIENT at that scale — analogous to a wavelet
coefficient but in tournament space.

ALSO: Lazy Wavelet (Haar-like lifting scheme) applied to pixel values.
This is lossless and reversible.

APPLICATIONS:
  1. Image compression (progressive, multi-resolution)
  2. Image denoising (threshold high-skip coefficients)
  3. Edge detection (large coefficients at skip-1 = edges)
  4. Texture analysis (coefficient distribution = texture signature)
  5. Super-resolution (predict missing skip levels from available)
"""

import sys
import os
import zlib
import bz2
import numpy as np
import time

__version__ = "1.0.0"


# ============================================================================
# TOURNAMENT WAVELET TRANSFORM
# ============================================================================

def twt_forward_1d(signal):
    """Forward Tournament Wavelet Transform (1D).

    Haar-like lifting scheme:
      1. Split into even and odd samples
      2. Detail = odd - even (prediction step)
      3. Approximation = even + detail//2 (update step)

    Returns: (approximation, detail)
    """
    n = len(signal)
    if n < 2:
        return signal, np.array([], dtype=signal.dtype)

    # Ensure even length
    if n % 2 == 1:
        signal = np.append(signal, signal[-1])
        n += 1

    even = signal[0::2].astype(np.int16)
    odd = signal[1::2].astype(np.int16)

    # Prediction: detail = odd - even (the "residual")
    detail = odd - even

    # Update: approximation = even + detail//2 (preserve mean)
    approx = even + detail // 2

    return approx.astype(np.int16), detail.astype(np.int16)


def twt_inverse_1d(approx, detail):
    """Inverse Tournament Wavelet Transform (1D)."""
    n_even = len(approx)
    n_detail = len(detail)

    even = approx.astype(np.int16) - detail[:n_even].astype(np.int16) // 2
    odd = detail[:n_even].astype(np.int16) + even

    result = np.zeros(n_even * 2, dtype=np.int16)
    result[0::2] = even
    result[1::2] = odd
    return result


def twt_forward_2d(img, levels=None):
    """Forward 2D Tournament Wavelet Transform.

    Applies 1D TWT to rows, then to columns.
    Produces multi-resolution pyramid like standard wavelet:
      LL | LH
      HL | HH

    Returns: list of (LL, LH, HL, HH) tuples for each level,
             plus final LL approximation.
    """
    if levels is None:
        levels = int(np.log2(min(img.shape))) - 1
        levels = max(1, min(levels, 6))

    current = img.astype(np.int16)
    coeffs = []

    for level in range(levels):
        H, W = current.shape

        # Apply to rows
        row_approx = np.zeros((H, W // 2), dtype=np.int16)
        row_detail = np.zeros((H, W // 2), dtype=np.int16)
        for i in range(H):
            a, d = twt_forward_1d(current[i])
            row_approx[i] = a[:W // 2]
            row_detail[i] = d[:W // 2]

        # Apply to columns of row_approx
        LL = np.zeros((H // 2, W // 2), dtype=np.int16)
        LH = np.zeros((H // 2, W // 2), dtype=np.int16)
        for j in range(W // 2):
            a, d = twt_forward_1d(row_approx[:, j])
            LL[:, j] = a[:H // 2]
            LH[:, j] = d[:H // 2]

        # Apply to columns of row_detail
        HL = np.zeros((H // 2, W // 2), dtype=np.int16)
        HH = np.zeros((H // 2, W // 2), dtype=np.int16)
        for j in range(W // 2):
            a, d = twt_forward_1d(row_detail[:, j])
            HL[:, j] = a[:H // 2]
            HH[:, j] = d[:H // 2]

        coeffs.append((LH, HL, HH))
        current = LL

    return current, coeffs


def twt_inverse_2d(LL, coeffs):
    """Inverse 2D Tournament Wavelet Transform."""
    current = LL.astype(np.int16)

    for level in range(len(coeffs) - 1, -1, -1):
        LH, HL, HH = coeffs[level]
        H2, W2 = current.shape

        # Inverse column transform on (LL, LH) -> row_approx
        row_approx = np.zeros((H2 * 2, W2), dtype=np.int16)
        for j in range(W2):
            row_approx[:, j] = twt_inverse_1d(current[:, j], LH[:, j])[:H2 * 2]

        # Inverse column transform on (HL, HH) -> row_detail
        row_detail = np.zeros((H2 * 2, W2), dtype=np.int16)
        for j in range(W2):
            row_detail[:, j] = twt_inverse_1d(HL[:, j], HH[:, j])[:H2 * 2]

        # Inverse row transform
        H, W = H2 * 2, W2 * 2
        current = np.zeros((H, W), dtype=np.int16)
        for i in range(H):
            current[i] = twt_inverse_1d(row_approx[i], row_detail[i])[:W]

    return current


# ============================================================================
# COMPRESSION USING TWT
# ============================================================================

def compress_twt(img, levels=None):
    """Compress image using Tournament Wavelet Transform.

    Strategy:
      1. Apply TWT to get multi-resolution pyramid
      2. Compress LL (low-res approximation) — small, smooth
      3. Compress each detail band — sparse, compressible
      4. Pack together with metadata
    """
    LL, coeffs = twt_forward_2d(img, levels)

    # Compress each component
    parts = []

    # LL: convert to uint8 range, delta encode
    ll_shifted = np.clip(LL + 128, 0, 255).astype(np.uint8)
    ll_bytes = ll_shifted.tobytes()
    ll_comp = _compress_best(ll_bytes)
    parts.append(ll_comp)

    # Detail bands: shift to uint8, compress
    for LH, HL, HH in coeffs:
        for band in [LH, HL, HH]:
            shifted = np.clip(band + 128, 0, 255).astype(np.uint8)
            band_comp = _compress_best(shifted.tobytes())
            parts.append(band_comp)

    # Pack: n_levels (1 byte) + H (2 bytes) + W (2 bytes) +
    #        for each part: size (4 bytes) + data
    n_levels = len(coeffs)
    header = bytes([n_levels]) + np.array([img.shape[0], img.shape[1]], dtype=np.uint16).tobytes()
    packed = header
    for part in parts:
        packed += len(part).to_bytes(4, 'big') + part

    return packed


def _compress_best(data):
    """Best of zlib/bz2."""
    results = [zlib.compress(data, 9)]
    try:
        results.append(bz2.compress(data, 9))
    except:
        pass
    return min(results, key=len)


# ============================================================================
# ANALYSIS
# ============================================================================

def analyze_coefficients(img, levels=None):
    """Analyze TWT coefficient statistics."""
    LL, coeffs = twt_forward_2d(img, levels)

    stats = {'LL_energy': np.sum(LL.astype(float)**2)}
    total_energy = stats['LL_energy']

    for i, (LH, HL, HH) in enumerate(coeffs):
        level_energy = (np.sum(LH.astype(float)**2) +
                       np.sum(HL.astype(float)**2) +
                       np.sum(HH.astype(float)**2))
        stats[f'level_{i}_energy'] = level_energy
        total_energy += level_energy

        # Sparsity: fraction of near-zero coefficients
        all_coeffs = np.concatenate([LH.flatten(), HL.flatten(), HH.flatten()])
        stats[f'level_{i}_sparsity'] = np.mean(np.abs(all_coeffs) <= 2)
        stats[f'level_{i}_max'] = np.max(np.abs(all_coeffs))

    stats['total_energy'] = total_energy
    stats['LL_fraction'] = stats['LL_energy'] / total_energy if total_energy > 0 else 0

    return stats


# ============================================================================
# BENCHMARK
# ============================================================================

def benchmark():
    """Benchmark TWT compression."""
    np.random.seed(42)

    print(f"Tournament Wavelet v{__version__}")
    print("=" * 90)

    N = 256
    tests = {}

    tests['grad_h'] = np.tile(np.arange(N, dtype=np.uint8), (N, 1))
    tests['grad_v'] = np.tile(np.arange(N, dtype=np.uint8).reshape(-1, 1), (1, N))
    tests['grad_diag'] = np.array([[(i+j) % 256 for j in range(N)] for i in range(N)], dtype=np.uint8)
    tests['checker'] = np.array([[(255 if (i//8 + j//8) % 2 == 0 else 0) for j in range(N)] for i in range(N)], dtype=np.uint8)
    tests['smooth'] = ((128 + 100 * np.sin(np.linspace(0, 4*np.pi, N)[None,:]) * np.cos(np.linspace(0, 4*np.pi, N)[:,None]))).clip(0, 255).astype(np.uint8)
    tests['circle'] = np.array([[255 if (i-N//2)**2 + (j-N//2)**2 < (N//3)**2 else 0 for j in range(N)] for i in range(N)], dtype=np.uint8)
    tests['sparse'] = np.zeros((N, N), dtype=np.uint8)
    for _ in range(50):
        x, y = np.random.randint(0, N), np.random.randint(0, N)
        tests['sparse'][max(0,x-3):min(N,x+3), max(0,y-3):min(N,y+3)] = 255
    tests['random'] = np.random.randint(0, 256, (N, N), dtype=np.uint8)

    print(f"\n  {'Name':>12} {'Raw':>8} {'TWT':>8} {'zlib9':>8} {'bz2':>8} {'TWT/best':>9} {'LL%':>6} {'L0-spar':>8} {'L2-spar':>8}")
    print("  " + "-" * 80)

    for name, img in tests.items():
        raw = img.nbytes

        # TWT compression
        twt_data = compress_twt(img)
        twt_size = len(twt_data)

        # Comparison
        zl = len(zlib.compress(img.tobytes(), 9))
        try: bz = len(bz2.compress(img.tobytes(), 9))
        except: bz = raw
        best = min(zl, bz)

        ratio = best / twt_size if twt_size > 0 else 0

        # Coefficient analysis
        stats = analyze_coefficients(img, levels=4)

        # Verify roundtrip
        LL, coeffs = twt_forward_2d(img, levels=4)
        recon = twt_inverse_2d(LL, coeffs)
        recon_clipped = np.clip(recon, 0, 255).astype(np.uint8)
        roundtrip_ok = np.array_equal(img, recon_clipped)

        tag = "WIN" if ratio > 1.02 else "TIE" if ratio > 0.98 else "LOSE"
        rt = "OK" if roundtrip_ok else "FAIL"

        print(f"  {name:>12} {raw:7d}B {twt_size:7d}B {zl:7d}B {bz:7d}B {ratio:8.2f}x "
              f"{stats['LL_fraction']:5.1%} {stats.get('level_0_sparsity', 0):7.1%} "
              f"{stats.get('level_2_sparsity', 0):7.1%} {tag} {rt}")

    print(f"\n  Note: TWT uses Haar-like lifting. Better wavelets (CDF 5/3, 9/7) would improve results.")
    print(f"  The key advantage is skip-hierarchy decomposition, not the wavelet basis.")


if __name__ == "__main__":
    benchmark()
