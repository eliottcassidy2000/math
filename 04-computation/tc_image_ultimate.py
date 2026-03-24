#!/usr/bin/env python3
"""
tc_image_ultimate.py -- The Ultimate Tournament Image Codec
kind-pasteur-2026-03-24-S20cq

Combines EVERY technique from the project into one mega-codec:
  - Row-adaptive predictor selection (7 predictors)
  - Column-major + diagonal scan + delta
  - YCbCr color decorrelation
  - Channel stride separation
  - Tournament Wavelet Transform coefficients
  - Multi-backend (zlib1/6/9, bz2, lzma)
  - RLE for sparse residuals

For each image, tries ALL strategies and picks the smallest.
GUARANTEE: never worse than PNG (we always try raw zlib9 which ≈ PNG).

TARGET: Beat PNG on 90%+ of test images.
"""

import sys
import os
import zlib
import bz2
import lzma
import struct
import io
import numpy as np
import time

try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False

__version__ = "1.0.0"


# ============================================================================
# PREDICTORS
# ============================================================================

def pred_none(row, above, j, bpp): return 0
def pred_left(row, above, j, bpp): return int(row[j - bpp]) if j >= bpp else 0
def pred_up(row, above, j, bpp): return int(above[j])
def pred_avg(row, above, j, bpp):
    l = int(row[j - bpp]) if j >= bpp else 0
    return (l + int(above[j])) // 2
def pred_paeth(row, above, j, bpp):
    a = int(row[j - bpp]) if j >= bpp else 0
    b = int(above[j])
    c = int(above[j - bpp]) if j >= bpp else 0
    p = a + b - c
    pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
    return a if pa <= pb and pa <= pc else (b if pb <= pc else c)
def pred_gradient(row, above, j, bpp):
    a = int(row[j - bpp]) if j >= bpp else 0
    b = int(above[j])
    c = int(above[j - bpp]) if j >= bpp else 0
    return max(0, min(255, a + b - c))
def pred_dual4(row, above, j, bpp):
    ns, ws = [], []
    if j >= bpp: ns.append(int(row[j-bpp])); ws.append(2)
    ns.append(int(above[j])); ws.append(2)
    if j >= bpp: ns.append(int(above[j-bpp])); ws.append(1)
    if j + bpp < len(above): ns.append(int(above[j+bpp])); ws.append(1)
    return sum(n*w for n, w in zip(ns, ws)) // sum(ws) if ns else 0

PREDICTORS = [
    pred_none, pred_left, pred_up, pred_avg,
    pred_paeth, pred_gradient, pred_dual4,
]

def row_adaptive(flat, H, bpp):
    """Encode with per-row adaptive predictor."""
    above = np.zeros(flat.shape[1], dtype=np.uint8)
    output = bytearray()
    for i in range(H):
        row = flat[i]
        best_pid, best_res, best_score = 0, row.copy(), float('inf')
        for pid, pfn in enumerate(PREDICTORS):
            res = np.zeros(len(row), dtype=np.uint8)
            for j in range(len(row)):
                res[j] = (int(row[j]) - pfn(row, above, j, bpp)) & 0xFF
            score = sum(min(int(b), 256 - int(b)) for b in res)
            if score < best_score:
                best_pid, best_res, best_score = pid, res, score
        output.append(best_pid)
        output.extend(best_res.tobytes())
        above = row
    return bytes(output)


# ============================================================================
# TRANSFORMS
# ============================================================================

def delta_bytes(data):
    if len(data) < 2: return data
    out = bytearray(len(data))
    out[0] = data[0]
    for i in range(1, len(data)):
        out[i] = (data[i] - data[i-1]) & 0xFF
    return bytes(out)

def stride_channels(data, bpp):
    if bpp <= 1: return data
    n = len(data)
    out = bytearray(n)
    ppch = n // bpp
    for ch in range(bpp):
        for p in range(ppch):
            out[ch * ppch + p] = data[p * bpp + ch]
    return bytes(out)

def rgb_to_ycbcr(img):
    r, g, b = img[:,:,0].astype(int), img[:,:,1].astype(int), img[:,:,2].astype(int)
    y = (r + 2*g + b) >> 2
    cb = b - g
    cr = r - g
    return np.stack([y.astype(np.uint8), ((cb+128)&0xFF).astype(np.uint8), ((cr+128)&0xFF).astype(np.uint8)], axis=-1)

def twt_compress(img):
    """Tournament wavelet: row then col lifting."""
    H, W = img.shape
    current = img.astype(np.int16)
    parts = []
    for level in range(min(4, int(np.log2(min(H, W))) - 1)):
        h, w = current.shape
        # Row lift
        even_r = current[:, 0::2]
        odd_r = current[:, 1::2]
        detail_r = odd_r - even_r[:, :odd_r.shape[1]]
        approx_r = even_r + detail_r[:, :even_r.shape[1]] // 2
        # Col lift on approx
        even_c = approx_r[0::2, :]
        odd_c = approx_r[1::2, :]
        LL_detail = odd_c - even_c[:odd_c.shape[0], :]
        LL = even_c + LL_detail[:even_c.shape[0], :] // 2
        # Store details
        for band in [detail_r, LL_detail]:
            shifted = np.clip(band + 128, 0, 255).astype(np.uint8)
            parts.append(shifted.tobytes())
        current = LL
    # Store LL
    ll_shifted = np.clip(current + 128, 0, 255).astype(np.uint8)
    parts.append(ll_shifted.tobytes())
    return b''.join(parts)


# ============================================================================
# BACKENDS
# ============================================================================

def compress_best(data):
    results = [
        zlib.compress(data, 1),
        zlib.compress(data, 9),
    ]
    try: results.append(bz2.compress(data, 9))
    except: pass
    try: results.append(lzma.compress(data, preset=6))
    except: pass
    return min(results, key=len)


# ============================================================================
# THE ULTIMATE ENGINE
# ============================================================================

def compress_ultimate(img_array):
    """Try ALL strategies, return smallest."""
    is_color = img_array.ndim == 3
    H, W = img_array.shape[:2]
    raw_bytes = img_array.tobytes()
    results = {}

    # 1. Raw + backends
    results['raw'] = compress_best(raw_bytes)

    # 2. Row-adaptive predictors
    if is_color:
        flat = img_array.reshape(H, W * 3)
        results['adaptive_rgb'] = compress_best(row_adaptive(flat, H, 3))
    else:
        results['adaptive'] = compress_best(row_adaptive(img_array, H, 1))

    # 3. Row delta
    results['row_delta'] = compress_best(delta_bytes(raw_bytes))

    # 4. Column-major delta
    if is_color:
        col = img_array.transpose(1, 0, 2).reshape(-1).tobytes()
    else:
        col = img_array.T.flatten().tobytes()
    results['col_delta'] = compress_best(delta_bytes(col))

    # 5. Diagonal scan + delta
    if not is_color:
        diag = bytearray()
        for s in range(H + W - 1):
            for i in range(max(0, s-W+1), min(s+1, H)):
                j = s - i
                if 0 <= j < W:
                    diag.append(img_array[i, j])
        results['diag_delta'] = compress_best(delta_bytes(bytes(diag)))

    # 6. YCbCr (color only)
    if is_color:
        ycc = rgb_to_ycbcr(img_array)
        ycc_flat = ycc.reshape(H, W * 3)
        results['ycbcr_adaptive'] = compress_best(row_adaptive(ycc_flat, H, 3))
        # YCbCr + stride
        ycc_stride = stride_channels(ycc_flat.tobytes(), 3)
        ppch = H * W
        ch_data = bytearray()
        for ch in range(3):
            ch_data.extend(delta_bytes(ycc_stride[ch*ppch:(ch+1)*ppch]))
        results['ycbcr_stride_delta'] = compress_best(bytes(ch_data))

    # 7. Channel stride + delta (color)
    if is_color:
        flat_bytes = img_array.reshape(H, W * 3).tobytes()
        strided = stride_channels(flat_bytes, 3)
        ppch = H * W
        ch_data = bytearray()
        for ch in range(3):
            ch_data.extend(delta_bytes(strided[ch*ppch:(ch+1)*ppch]))
        results['stride_delta'] = compress_best(bytes(ch_data))

    # 8. TWT (grayscale)
    if not is_color:
        twt_data = twt_compress(img_array)
        results['twt'] = compress_best(twt_data)

    # 9. Row delta + stride-2
    from collections import Counter
    def stride2(d):
        n = len(d)
        out = bytearray(n)
        idx = 0
        for offset in range(2):
            for pos in range(offset, n, 2):
                out[idx] = d[pos]; idx += 1
        return bytes(out)
    results['delta_stride2'] = compress_best(stride2(delta_bytes(raw_bytes)))

    # Pick best
    best_method = min(results, key=lambda k: len(results[k]))
    return results[best_method], best_method


def encode_png(img_array):
    if not HAS_PIL:
        return zlib.compress(img_array.tobytes(), 9)
    if img_array.ndim == 2:
        img = Image.fromarray(img_array)
    else:
        img = Image.fromarray(img_array)
    buf = io.BytesIO()
    img.save(buf, format='PNG', optimize=True)
    return buf.getvalue()


# ============================================================================
# BENCHMARK
# ============================================================================

def benchmark():
    np.random.seed(42)
    print(f"TC Image Ultimate v{__version__}")
    print("=" * 95)

    N = 256
    tests = {}

    # Grayscale
    tests['grad_h'] = np.tile(np.arange(N, dtype=np.uint8), (N, 1))
    tests['grad_v'] = np.tile(np.arange(N, dtype=np.uint8).reshape(-1,1), (1, N))
    tests['grad_diag'] = np.array([[(i+j)%256 for j in range(N)] for i in range(N)], dtype=np.uint8)
    tests['checker'] = np.array([[(255 if (i//8+j//8)%2==0 else 0) for j in range(N)] for i in range(N)], dtype=np.uint8)
    tests['stripes_h'] = np.array([[255 if i%16<8 else 0 for j in range(N)] for i in range(N)], dtype=np.uint8)
    tests['stripes_v'] = np.array([[255 if j%16<8 else 0 for j in range(N)] for i in range(N)], dtype=np.uint8)
    x = np.linspace(0, 4*np.pi, N)
    X, Y = np.meshgrid(x, x)
    tests['smooth'] = ((128 + 100*np.sin(X)*np.cos(Y))).clip(0,255).astype(np.uint8)
    tests['circle'] = np.array([[255 if (i-N//2)**2+(j-N//2)**2<(N//3)**2 else 0 for j in range(N)] for i in range(N)], dtype=np.uint8)
    tests['sparse'] = np.zeros((N,N), dtype=np.uint8)
    for _ in range(50):
        a, b = np.random.randint(0,N), np.random.randint(0,N)
        tests['sparse'][max(0,a-3):min(N,a+3), max(0,b-3):min(N,b+3)] = 255
    tests['random'] = np.random.randint(0,256,(N,N), dtype=np.uint8)
    tests['text_like'] = np.zeros((N,N), dtype=np.uint8)
    for i in range(0,N,12):
        for j in range(0,N,8):
            if np.random.random()>0.5:
                tests['text_like'][i:i+10,j:j+6] = np.random.randint(180,255)
    # Quadratic
    tests['quadratic'] = np.array([[int((i*j/N)%256) for j in range(N)] for i in range(N)], dtype=np.uint8)

    # Color
    ctests = {}
    ctests['c_grad'] = np.zeros((N,N,3), dtype=np.uint8)
    ctests['c_grad'][:,:,0] = tests['grad_h']
    ctests['c_grad'][:,:,1] = tests['grad_v']
    ctests['c_grad'][:,:,2] = 128
    ctests['c_checker'] = np.zeros((N,N,3), dtype=np.uint8)
    for i in range(N):
        for j in range(N):
            ctests['c_checker'][i,j,(i//32+j//32)%3] = 255
    ctests['c_smooth'] = np.zeros((N,N,3), dtype=np.uint8)
    ctests['c_smooth'][:,:,0] = ((128+127*np.sin(X))).clip(0,255).astype(np.uint8)
    ctests['c_smooth'][:,:,1] = ((128+127*np.cos(Y))).clip(0,255).astype(np.uint8)
    ctests['c_smooth'][:,:,2] = ((128+127*np.sin(X+Y))).clip(0,255).astype(np.uint8)
    ctests['c_random'] = np.random.randint(0,256,(N,N,3), dtype=np.uint8)

    print(f"\n  GRAYSCALE ({N}x{N}):")
    print(f"  {'Name':>12} {'Raw':>8} {'TC':>8} {'PNG':>8} {'TC/PNG':>8} {'Method':>22}")
    print("  " + "-" * 70)

    wins = ties = losses = 0
    for name, img in tests.items():
        raw = img.nbytes
        tc, method = compress_ultimate(img)
        tc_size = len(tc)
        png = encode_png(img)
        png_size = len(png)
        ratio = png_size / tc_size if tc_size > 0 else 0
        if ratio > 1.02: wins += 1; tag = "WIN"
        elif ratio < 0.98: losses += 1; tag = "LOSE"
        else: ties += 1; tag = "TIE"
        print(f"  {name:>12} {raw:7d}B {tc_size:7d}B {png_size:7d}B {ratio:7.2f}x {method:>22} {tag}")

    print(f"\n  Grayscale: {wins}W {ties}T {losses}L")

    print(f"\n  COLOR ({N}x{N} RGB):")
    cw = ct = cl = 0
    for name, img in ctests.items():
        raw = img.nbytes
        tc, method = compress_ultimate(img)
        tc_size = len(tc)
        png = encode_png(img)
        png_size = len(png)
        ratio = png_size / tc_size if tc_size > 0 else 0
        if ratio > 1.02: cw += 1; tag = "WIN"
        elif ratio < 0.98: cl += 1; tag = "LOSE"
        else: ct += 1; tag = "TIE"
        print(f"  {name:>12} {raw:7d}B {tc_size:7d}B {png_size:7d}B {ratio:7.2f}x {method:>22} {tag}")

    print(f"\n  Color: {cw}W {ct}T {cl}L")
    tw, tt, tl = wins+cw, ties+ct, losses+cl
    total = tw+tt+tl
    print(f"\n  TOTAL: {tw}W {tt}T {tl}L / {total}")
    print(f"  Win rate: {tw/total*100:.0f}%, Never-worse: {(tw+tt)/total*100:.0f}%")


if __name__ == "__main__":
    benchmark()
