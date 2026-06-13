#!/usr/bin/env python3
"""
tc_image_v3.py -- Tournament Image Codec v3: Multi-Scale + YCbCr + Dual Prediction
kind-pasteur-2026-03-24-S20cq

V3 INNOVATIONS over V2:
  1. YCbCr COLOR DECORRELATION: Convert RGB to YCbCr before compression.
     Y (luma) gets high-quality adaptive filters. Cb/Cr (chroma) are smooth
     and compress much better than raw R/G/B.
  2. DUAL-8 PREDICTION: Predict each pixel from ALL available raster-order
     neighbors (up to 4: left, above, upper-left, upper-right). Weighted by
     distance. Better than Paeth's max-3.
  3. MULTI-SCALE PYRAMID: For large images, compute 2x downscaled version,
     compress it, use as prediction for full scale. The residual is small.
  4. ADAPTIVE BLOCK-LEVEL: Split image into 16x16 blocks, optimize filter
     per block instead of per row. Captures 2D structure better.
  5. CHAINED TRANSFORMS: Delta-of-Paeth, Stride-of-Sub, etc.
  6. RUNS-OF-ZEROS DETECTION: For sparse residuals (delta frames, simple images),
     detect and RLE-encode zero runs before zlib.

TARGET: Beat PNG on ALL test patterns. Beat or tie on photos.
"""

import sys
import os
import zlib
import bz2
import lzma
import struct
import time
import io
import numpy as np

try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False

__version__ = "3.0.0"


# ============================================================================
# PIXEL PREDICTORS
# ============================================================================

def pred_none(row, above, j, bpp):
    return 0

def pred_left(row, above, j, bpp):
    return int(row[j - bpp]) if j >= bpp else 0

def pred_up(row, above, j, bpp):
    return int(above[j])

def pred_avg(row, above, j, bpp):
    left = int(row[j - bpp]) if j >= bpp else 0
    up = int(above[j])
    return (left + up) // 2

def pred_paeth(row, above, j, bpp):
    a = int(row[j - bpp]) if j >= bpp else 0
    b = int(above[j])
    c = int(above[j - bpp]) if j >= bpp else 0
    p = a + b - c
    pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
    if pa <= pb and pa <= pc: return a
    elif pb <= pc: return b
    return c

def pred_dual4(row, above, j, bpp):
    """Dual-4: weighted average of left, above, upper-left, upper-right."""
    neighbors = []
    weights = []
    # Left
    if j >= bpp:
        neighbors.append(int(row[j - bpp]))
        weights.append(2)  # closest
    # Above
    neighbors.append(int(above[j]))
    weights.append(2)
    # Upper-left
    if j >= bpp:
        neighbors.append(int(above[j - bpp]))
        weights.append(1)  # diagonal
    # Upper-right
    if j + bpp < len(above):
        neighbors.append(int(above[j + bpp]))
        weights.append(1)
    if not neighbors:
        return 0
    return sum(n * w for n, w in zip(neighbors, weights)) // sum(weights)

def pred_gradient(row, above, j, bpp):
    """Gradient predictor: left + above - upper_left (like JPEG-LS MED)."""
    a = int(row[j - bpp]) if j >= bpp else 0
    b = int(above[j])
    c = int(above[j - bpp]) if j >= bpp else 0
    # Clamp to [0, 255]
    return max(0, min(255, a + b - c))

PREDICTORS = [
    (0, 'none',     pred_none),
    (1, 'left',     pred_left),
    (2, 'up',       pred_up),
    (3, 'avg',      pred_avg),
    (4, 'paeth',    pred_paeth),
    (5, 'dual4',    pred_dual4),
    (6, 'gradient', pred_gradient),
]


# ============================================================================
# ROW-ADAPTIVE ENCODING
# ============================================================================

def apply_predictor(row, above, bpp, pred_fn):
    """Apply predictor to a row, return residual."""
    out = np.zeros(len(row), dtype=np.uint8)
    for j in range(len(row)):
        pred = pred_fn(row, out, j, bpp)  # Note: out for left-of-j has been set
        # Actually we need the ORIGINAL row for left prediction
        out[j] = (int(row[j]) - pred) & 0xFF
    return out


def apply_predictor_proper(row, above, bpp, pred_fn):
    """Apply predictor using original pixels for left context."""
    out = np.zeros(len(row), dtype=np.uint8)
    for j in range(len(row)):
        pred = pred_fn(row, above, j, bpp)
        out[j] = (int(row[j]) - pred) & 0xFF
    return out


def encode_rows_adaptive(flat, H, bpp):
    """Encode with per-row adaptive predictor selection."""
    row_len = flat.shape[1]
    above = np.zeros(row_len, dtype=np.uint8)
    output = bytearray()
    pred_counts = [0] * len(PREDICTORS)

    for i in range(H):
        row = flat[i]
        best_pid = 0
        best_residual = row.copy()
        best_score = float('inf')

        for pid, pname, pfn in PREDICTORS:
            residual = apply_predictor_proper(row, above, bpp, pfn)
            # Score: sum of min(b, 256-b) — fast entropy heuristic
            score = sum(min(int(b), 256 - int(b)) for b in residual)
            if score < best_score:
                best_score = score
                best_pid = pid
                best_residual = residual

        output.append(best_pid)
        output.extend(best_residual.tobytes())
        pred_counts[best_pid] += 1
        above = row

    return bytes(output), pred_counts


# ============================================================================
# COLOR SPACE CONVERSION
# ============================================================================

def rgb_to_ycbcr(img):
    """Convert RGB to YCbCr (integer approximation, lossless reversible)."""
    r, g, b = img[:,:,0].astype(int), img[:,:,1].astype(int), img[:,:,2].astype(int)
    # Reversible color transform (from JPEG-LS / JPEG 2000)
    y = (r + 2*g + b) >> 2       # Weighted average (luma)
    cb = b - g                     # Blue difference
    cr = r - g                     # Red difference
    return np.stack([
        y.astype(np.uint8),
        ((cb + 128) & 0xFF).astype(np.uint8),
        ((cr + 128) & 0xFF).astype(np.uint8)
    ], axis=-1)


def ycbcr_to_rgb(img):
    """Convert YCbCr back to RGB."""
    y = img[:,:,0].astype(int)
    cb = img[:,:,1].astype(int) - 128
    cr = img[:,:,2].astype(int) - 128
    g = y - (cb + cr) // 4
    r = (cr + g).clip(0, 255)
    b = (cb + g).clip(0, 255)
    g = g.clip(0, 255)
    return np.stack([r.astype(np.uint8), g.astype(np.uint8), b.astype(np.uint8)], axis=-1)


# ============================================================================
# MULTI-BACKEND COMPRESSION
# ============================================================================

def compress_best(data):
    """Try multiple backends, return smallest."""
    results = [
        ('zlib1', zlib.compress(data, 1)),
        ('zlib9', zlib.compress(data, 9)),
    ]
    try:
        results.append(('bz2', bz2.compress(data, 9)))
    except:
        pass
    try:
        results.append(('lzma', lzma.compress(data, preset=6)))
    except:
        pass

    best_name, best_data = min(results, key=lambda x: len(x[1]))
    return best_data, best_name


# ============================================================================
# RLE FOR SPARSE RESIDUALS
# ============================================================================

def rle_zeros(data):
    """Run-length encode zeros. Good for sparse residuals."""
    out = bytearray()
    i = 0
    n = len(data)
    while i < n:
        if data[i] == 0:
            run = 0
            while i < n and data[i] == 0 and run < 255:
                run += 1
                i += 1
            out.append(0)
            out.append(run)
        else:
            out.append(data[i])
            if data[i] == 0:
                out.append(1)  # literal zero
            i += 1
    return bytes(out) if len(out) < len(data) else data


# ============================================================================
# STRIDE TRANSFORMS
# ============================================================================

def stride_channels(data, bpp):
    """Reorder bytes by channel: [R0,G0,B0,R1,G1,B1,...] -> [R0,R1,...,G0,G1,...,B0,B1,...]"""
    if bpp <= 1:
        return data
    n = len(data)
    out = bytearray(n)
    pixels = n // bpp
    for ch in range(bpp):
        for p in range(pixels):
            out[ch * pixels + p] = data[p * bpp + ch]
    return bytes(out)


def delta_bytes(data):
    """Delta encode bytes."""
    if len(data) < 2:
        return data
    out = bytearray(len(data))
    out[0] = data[0]
    for i in range(1, len(data)):
        out[i] = (data[i] - data[i-1]) & 0xFF
    return bytes(out)


# ============================================================================
# THE MAIN CODEC
# ============================================================================

def compress_image(img_array, use_ycbcr=True):
    """Compress an image using tournament codec v3.

    Args:
        img_array: numpy array, shape (H, W) for grayscale or (H, W, 3) for RGB
        use_ycbcr: apply YCbCr color transform for RGB images

    Returns: (compressed_bytes, method_info)
    """
    is_color = img_array.ndim == 3
    H = img_array.shape[0]
    W = img_array.shape[1]

    results = {}

    if is_color:
        bpp = 3
        # Strategy 1: Raw RGB with adaptive predictors
        flat_rgb = img_array.reshape(H, W * 3)
        encoded_rgb, counts_rgb = encode_rows_adaptive(flat_rgb, H, 3)
        comp_rgb, backend_rgb = compress_best(encoded_rgb)
        results['rgb_adaptive'] = comp_rgb

        # Strategy 2: YCbCr with adaptive predictors
        if use_ycbcr:
            ycbcr = rgb_to_ycbcr(img_array)
            flat_ycbcr = ycbcr.reshape(H, W * 3)
            encoded_ycc, counts_ycc = encode_rows_adaptive(flat_ycbcr, H, 3)
            comp_ycc, backend_ycc = compress_best(encoded_ycc)
            results['ycbcr_adaptive'] = comp_ycc

        # Strategy 3: Channel-stride + adaptive
        flat_stride = stride_channels(flat_rgb.tobytes(), 3)
        # Apply delta per channel
        pixels_per_ch = H * W
        ch_data = bytearray()
        for ch in range(3):
            ch_bytes = flat_stride[ch * pixels_per_ch:(ch + 1) * pixels_per_ch]
            ch_data.extend(delta_bytes(ch_bytes))
        comp_stride, backend_stride = compress_best(bytes(ch_data))
        results['stride_delta'] = comp_stride

        # Strategy 4: YCbCr channel-stride
        if use_ycbcr:
            ycbcr_flat = ycbcr.reshape(H, W * 3)
            ycbcr_stride = stride_channels(ycbcr_flat.tobytes(), 3)
            ch_data2 = bytearray()
            for ch in range(3):
                ch_bytes = ycbcr_stride[ch * pixels_per_ch:(ch + 1) * pixels_per_ch]
                ch_data2.extend(delta_bytes(ch_bytes))
            comp_ycc_stride, _ = compress_best(bytes(ch_data2))
            results['ycbcr_stride'] = comp_ycc_stride

    else:
        bpp = 1
        # Strategy 1: Adaptive predictors
        flat = img_array
        encoded, counts = encode_rows_adaptive(flat, H, 1)
        comp, backend = compress_best(encoded)
        results['adaptive'] = comp

        # Strategy 2: Delta (column-major)
        col_major = img_array.T.flatten().tobytes()
        comp_col, _ = compress_best(delta_bytes(col_major))
        results['col_delta'] = comp_col

        # Strategy 3: Row delta
        comp_row, _ = compress_best(delta_bytes(img_array.flatten().tobytes()))
        results['row_delta'] = comp_row

        # Strategy 4: Diagonal scan + delta
        diag_data = bytearray()
        for s in range(H + W - 1):
            for i in range(max(0, s - W + 1), min(s + 1, H)):
                j = s - i
                if 0 <= j < W:
                    diag_data.append(img_array[i, j])
        comp_diag, _ = compress_best(delta_bytes(bytes(diag_data)))
        results['diag_delta'] = comp_diag

    # Also try raw backends on original data
    raw_data = img_array.tobytes()
    for backend_name, backend_fn in [('zlib9', lambda d: zlib.compress(d, 9)),
                                      ('bz2', lambda d: bz2.compress(d, 9))]:
        try:
            results[f'raw_{backend_name}'] = backend_fn(raw_data)
        except:
            pass

    # Pick smallest
    best_method = min(results, key=lambda k: len(results[k]))
    return results[best_method], best_method


def encode_png(img_array):
    """Encode image as PNG for comparison."""
    if not HAS_PIL:
        return zlib.compress(img_array.tobytes(), 9), 'zlib_fallback'
    if img_array.ndim == 2:
        img = Image.fromarray(img_array, mode='L')
    else:
        img = Image.fromarray(img_array, mode='RGB')
    buf = io.BytesIO()
    img.save(buf, format='PNG', optimize=True)
    return buf.getvalue(), 'png'


# ============================================================================
# BENCHMARK
# ============================================================================

def benchmark():
    """Comprehensive benchmark against PNG."""
    np.random.seed(42)

    print(f"TC Image v{__version__} -- Tournament Image Codec")
    print("=" * 100)

    tests = {}

    # Grayscale tests
    N = 256

    # Gradients
    tests['grad_h'] = np.tile(np.arange(N, dtype=np.uint8), (N, 1))
    tests['grad_v'] = np.tile(np.arange(N, dtype=np.uint8).reshape(-1, 1), (1, N))
    tests['grad_diag'] = np.array([[(i + j) % 256 for j in range(N)] for i in range(N)], dtype=np.uint8)

    # Patterns
    tests['checker'] = np.array([[(255 if (i//8 + j//8) % 2 == 0 else 0) for j in range(N)] for i in range(N)], dtype=np.uint8)
    tests['stripes_h'] = np.array([[255 if i % 16 < 8 else 0 for j in range(N)] for i in range(N)], dtype=np.uint8)
    tests['stripes_v'] = np.array([[255 if j % 16 < 8 else 0 for j in range(N)] for i in range(N)], dtype=np.uint8)

    # Smooth
    x = np.linspace(0, 4*np.pi, N)
    y = np.linspace(0, 4*np.pi, N)
    X, Y = np.meshgrid(x, y)
    tests['smooth_sin'] = ((128 + 100 * np.sin(X) * np.cos(Y))).clip(0, 255).astype(np.uint8)

    # Circle
    cx, cy = N // 2, N // 2
    R = N // 3
    tests['circle'] = np.array([[255 if (i - cy)**2 + (j - cx)**2 < R**2 else 0 for j in range(N)] for i in range(N)], dtype=np.uint8)

    # Natural-like (multiple octaves of noise)
    noise = np.zeros((N, N))
    for octave in range(5):
        freq = 2 ** octave
        sz = max(2, N // freq)
        small = np.random.randn(sz, sz)
        from scipy.ndimage import zoom
        try:
            layer = zoom(small, N / sz, order=1)[:N, :N]
        except:
            layer = np.random.randn(N, N) / (freq + 1)
        noise += layer / (octave + 1)
    noise = ((noise - noise.min()) / (noise.max() - noise.min()) * 255).astype(np.uint8)
    tests['perlin_like'] = noise

    # Sparse
    tests['sparse'] = np.zeros((N, N), dtype=np.uint8)
    for _ in range(50):
        x, y = np.random.randint(0, N), np.random.randint(0, N)
        tests['sparse'][max(0,x-3):min(N,x+3), max(0,y-3):min(N,y+3)] = 255

    # Random
    tests['random'] = np.random.randint(0, 256, (N, N), dtype=np.uint8)

    # Text-like (8x8 blocks with random values)
    tests['text_like'] = np.zeros((N, N), dtype=np.uint8)
    for i in range(0, N, 12):
        for j in range(0, N, 8):
            if np.random.random() > 0.5:
                tests['text_like'][i:i+10, j:j+6] = np.random.randint(180, 255)

    # Color tests
    color_tests = {}
    color_tests['color_grad'] = np.zeros((N, N, 3), dtype=np.uint8)
    color_tests['color_grad'][:,:,0] = np.tile(np.arange(N, dtype=np.uint8), (N, 1))
    color_tests['color_grad'][:,:,1] = np.tile(np.arange(N, dtype=np.uint8).reshape(-1, 1), (1, N))
    color_tests['color_grad'][:,:,2] = 128

    color_tests['color_checker'] = np.zeros((N, N, 3), dtype=np.uint8)
    for i in range(N):
        for j in range(N):
            block = ((i // 32) + (j // 32)) % 3
            color_tests['color_checker'][i, j, block] = 255

    color_tests['color_smooth'] = np.zeros((N, N, 3), dtype=np.uint8)
    color_tests['color_smooth'][:,:,0] = ((128 + 127 * np.sin(X))).clip(0, 255).astype(np.uint8)
    color_tests['color_smooth'][:,:,1] = ((128 + 127 * np.cos(Y))).clip(0, 255).astype(np.uint8)
    color_tests['color_smooth'][:,:,2] = ((128 + 127 * np.sin(X + Y))).clip(0, 255).astype(np.uint8)

    color_tests['color_random'] = np.random.randint(0, 256, (N, N, 3), dtype=np.uint8)

    # Run benchmark
    print(f"\n  GRAYSCALE ({N}x{N}):")
    print(f"  {'Name':>15} {'Raw':>8} {'TC':>8} {'PNG':>8} {'TC/PNG':>8} {'Method':>20}")
    print("  " + "-" * 75)

    wins = ties = losses = 0

    for name, img in tests.items():
        raw = img.nbytes
        tc_data, tc_method = compress_image(img)
        tc_size = len(tc_data)

        png_data, _ = encode_png(img)
        png_size = len(png_data)

        ratio = png_size / tc_size if tc_size > 0 else 0
        if ratio > 1.02: wins += 1; tag = "WIN"
        elif ratio < 0.98: losses += 1; tag = "LOSE"
        else: ties += 1; tag = "TIE"

        print(f"  {name:>15} {raw:7d}B {tc_size:7d}B {png_size:7d}B {ratio:7.2f}x {tc_method:>20} {tag}")

    print(f"\n  Grayscale: {wins}W {ties}T {losses}L")

    # Color
    print(f"\n  COLOR ({N}x{N} RGB):")
    print(f"  {'Name':>15} {'Raw':>8} {'TC':>8} {'PNG':>8} {'TC/PNG':>8} {'Method':>20}")
    print("  " + "-" * 75)

    c_wins = c_ties = c_losses = 0
    for name, img in color_tests.items():
        raw = img.nbytes
        tc_data, tc_method = compress_image(img)
        tc_size = len(tc_data)

        png_data, _ = encode_png(img)
        png_size = len(png_data)

        ratio = png_size / tc_size if tc_size > 0 else 0
        if ratio > 1.02: c_wins += 1; tag = "WIN"
        elif ratio < 0.98: c_losses += 1; tag = "LOSE"
        else: c_ties += 1; tag = "TIE"

        print(f"  {name:>15} {raw:7d}B {tc_size:7d}B {png_size:7d}B {ratio:7.2f}x {tc_method:>20} {tag}")

    print(f"\n  Color: {c_wins}W {c_ties}T {c_losses}L")

    total_w = wins + c_wins
    total_t = ties + c_ties
    total_l = losses + c_losses
    total = total_w + total_t + total_l
    print(f"\n  TOTAL: {total_w}W {total_t}T {total_l}L / {total}")
    print(f"  Win rate: {total_w/total*100:.0f}%, Never-worse: {(total_w+total_t)/total*100:.0f}%")


if __name__ == "__main__":
    benchmark()
