#!/usr/bin/env python3
"""
quincunx.py — Quincunx Pyramid Codec

A genuinely new multiresolution image compression approach.

QUINCUNX SUBSAMPLING: Instead of 2:1 in both dimensions (loses diagonals),
use checkerboard subsampling. At each level, keep pixels where (r+c) % 2 == 0.

Level 0: all N*N pixels (full resolution)
Level 1: N*N/2 "white" squares (checkerboard phase 0)
Level 2: N*N/4 (checkerboard of level 1)
Level K: N*N/2^K pixels

Encoding: encode from COARSEST level down to finest.
At each level, predict from the next-coarser level (all neighbors available).

Key property: at every level, prediction uses ALL 4 diagonal neighbors from
the coarser level — not just left/above like raster. This gives 2D isotropic
prediction regardless of scan order.

The quincunx lattice IS the dual of the checkerboard. And:
- The quincunx = rotation by 45 degrees + 2:1 scale
- It's used in image processing for anti-aliasing
- New here: using it as a COMPRESSION pyramid

Combined with spiral scan within each level for optimal 1D ordering.

Author: opus-2026-03-25-S342
"""

import sys, struct, zlib, lzma, io, time
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

def _quincunx_levels(H, W, max_levels=None):
    """Generate quincunx pyramid levels.
    Each level is a list of (r,c) coordinates.
    Level 0 = all pixels. Level k+1 = checkerboard of level k."""
    if max_levels is None:
        max_levels = max(1, int(np.log2(min(H, W))) - 1)

    levels = []
    current = set()
    for r in range(H):
        for c in range(W):
            current.add((r, c))
    levels.append(sorted(current))

    for k in range(max_levels):
        # Keep only (r+c) % 2 == 0 from current
        next_level = set()
        for r, c in current:
            if (r + c) % 2 == 0:
                next_level.add((r, c))
        if len(next_level) < 4:
            break
        levels.append(sorted(next_level))
        current = next_level

    return levels


def _predict_quincunx(plane, r, c, known_set, H, W):
    """Predict pixel from known neighbors. Uses all 8 neighbors weighted by distance."""
    total = 0.0
    weight = 0.0
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            if dr == 0 and dc == 0: continue
            nr, nc = r + dr, c + dc
            if 0 <= nr < H and 0 <= nc < W and (nr, nc) in known_set:
                # Cardinal neighbors weight 1, diagonal weight 1/sqrt(2)
                w = 1.0 if (abs(dr) + abs(dc) == 1) else 0.707
                total += w * int(plane[nr, nc])
                weight += w
    if weight > 0:
        return int(round(total / weight))
    return 128


def _spiral_order(coords, H, W):
    """Sort coordinates in spiral order from center."""
    cr, cc = H / 2, W / 2
    return sorted(coords, key=lambda rc: (rc[0] - cr)**2 + (rc[1] - cc)**2)


def encode_quincunx(plane):
    """Encode a single plane using quincunx pyramid."""
    H, W = plane.shape
    levels = _quincunx_levels(H, W)
    n_levels = len(levels)

    # Start from coarsest level
    coarsest = levels[-1]
    known = set()
    residuals = bytearray()

    # Encode coarsest level raw (small — spiral ordered)
    coarsest_spiral = _spiral_order(coarsest, H, W)
    for r, c in coarsest_spiral:
        residuals.append(int(plane[r, c]))
        known.add((r, c))

    # Encode each finer level's NEW pixels
    level_sizes = [len(coarsest)]
    for k in range(n_levels - 2, -1, -1):  # from second-coarsest to finest
        new_pixels = [rc for rc in levels[k] if rc not in known]
        new_spiral = _spiral_order(new_pixels, H, W)

        count = 0
        for r, c in new_spiral:
            pred = _predict_quincunx(plane, r, c, known, H, W)
            residuals.append((int(plane[r, c]) - pred) & 0xFF)
            known.add((r, c))
            count += 1
        level_sizes.append(count)

    return bytes(residuals), n_levels, level_sizes


def decode_quincunx(data, H, W, n_levels, level_sizes):
    """Decode a plane from quincunx pyramid."""
    levels = _quincunx_levels(H, W, max_levels=n_levels - 1)

    plane = np.zeros((H, W), dtype=np.uint8)
    known = set()
    offset = 0

    # Decode coarsest level (raw)
    coarsest = levels[-1]
    coarsest_spiral = _spiral_order(coarsest, H, W)
    for r, c in coarsest_spiral:
        plane[r, c] = data[offset]; offset += 1
        known.add((r, c))

    # Decode each finer level
    for k in range(n_levels - 2, -1, -1):
        new_pixels = [rc for rc in levels[k] if rc not in known]
        new_spiral = _spiral_order(new_pixels, H, W)

        for r, c in new_spiral:
            pred = _predict_quincunx(plane, r, c, known, H, W)
            plane[r, c] = (int(data[offset]) + pred) & 0xFF; offset += 1
            known.add((r, c))

    return plane


def png_size(img):
    from PIL import Image
    mode = 'L' if img.ndim == 2 else 'RGB'
    pil = Image.fromarray(img, mode)
    buf = io.BytesIO()
    pil.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()


# ============================================================
# TEST
# ============================================================

if __name__ == "__main__":
    np.random.seed(42)
    print("=" * 80)
    print("  Quincunx Pyramid Codec — Multiresolution Checkerboard Compression")
    print("  opus-2026-03-25-S342")
    print("=" * 80)

    tests = {}
    for N in [64, 128, 256]:
        tests[f"circle_{N}"] = np.array([[min(255, int(255*np.exp(-((r-N//2)**2+(c-N//2)**2)/(N**2/8)))) for c in range(N)] for r in range(N)], dtype=np.uint8)
        tests[f"gradient_{N}"] = np.tile(np.linspace(0, 255, N, dtype=np.uint8), (N, 1))
        tests[f"diag_grad_{N}"] = np.array([[int(255*(r+c)/(2*N-2)) for c in range(N)] for r in range(N)], dtype=np.uint8)
        tests[f"checker_{N}"] = np.array([[(r+c)%2*255 for c in range(N)] for r in range(N)], dtype=np.uint8)
        tests[f"noise_{N}"] = np.random.randint(0, 256, (N, N), dtype=np.uint8)
        tests[f"smooth_{N}"] = np.clip(np.array([[int(128+60*np.sin(r/5)*np.cos(c/4)) for c in range(N)] for r in range(N)]), 0, 255).astype(np.uint8)
        tests[f"radial_{N}"] = np.array([[int(128+127*np.sin(np.sqrt((r-N//2)**2+(c-N//2)**2)/5)) for c in range(N)] for r in range(N)], dtype=np.uint8)

    print(f"\n  {'Image':<25} {'Raw':>7} {'Quin':>7} {'PNG':>7} {'Q/PNG':>7} {'Levels':>7} {'OK':>4}")
    print("  " + "-" * 68)

    wins = 0; total = 0; all_ok = True
    for name, img in sorted(tests.items()):
        total += 1
        raw = img.size

        t0 = time.time()
        residuals, n_levels, level_sizes = encode_quincunx(img)
        t_enc = time.time() - t0

        # Compress residuals
        compressed = zlib.compress(residuals, 9)
        # Also try lzma
        try:
            lzma_c = lzma.compress(residuals, preset=9)
            if len(lzma_c) < len(compressed):
                compressed = lzma_c
        except: pass

        quin_size = len(compressed) + 8  # small header overhead

        # Verify roundtrip
        recovered = decode_quincunx(residuals, img.shape[0], img.shape[1], n_levels, level_sizes)
        ok = np.array_equal(img, recovered)
        if not ok: all_ok = False

        ps = png_size(img)
        ratio = ps / quin_size if quin_size > 0 else 0
        status = "WIN" if ratio > 1.001 else "LOSS"
        if ratio > 1.001: wins += 1

        print(f"  {name:<25} {raw:>7} {quin_size:>7} {ps:>7} {ratio:>6.2f}x {n_levels:>3}lev {'OK' if ok else 'FAIL':>4}")

    print(f"\n  {wins}/{total} beat PNG. All lossless: {all_ok}")

    # Real photo test
    print("\n--- Real Photos ---")
    import os
    from PIL import Image as PILImage
    for path in ['/home/e/tron.jpeg', '/home/e/t.jpg']:
        if os.path.exists(path):
            img = np.array(PILImage.open(path).convert('L'))
            H, W = img.shape
            residuals, n_levels, level_sizes = encode_quincunx(img)
            compressed = zlib.compress(residuals, 9)
            try:
                lzma_c = lzma.compress(residuals, preset=9)
                if len(lzma_c) < len(compressed): compressed = lzma_c
            except: pass
            quin_size = len(compressed) + 8
            recovered = decode_quincunx(residuals, H, W, n_levels, level_sizes)
            ok = np.array_equal(img, recovered)
            ps = png_size(img)
            print(f"  {os.path.basename(path):<25} Quin={quin_size:>7} PNG={ps:>7} {ps/quin_size:.3f}x {'OK' if ok else 'FAIL'}")
