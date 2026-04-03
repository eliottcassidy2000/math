#!/usr/bin/env python3
"""
TIC-Mosaic v4: Full-Image First, Block-Adaptive Refinement

Strategy: try multiple full-image encoding strategies first (fast),
then only do block-level adaptive encoding if the image is mixed content.

Key changes from v3:
- Full-image encoding is the primary path (much faster)
- Block-level encoding only when content is genuinely mixed
- Faster: reuse zstd compressor, skip unnecessary strategies
- Added: YCoCg-R color transform option (lossless, often better than G-RG-BG)
- Added: delta coding (row differencing) before compression

opus-2026-04-03
"""
import sys, io, os, zlib, struct, math, time
import numpy as np
from PIL import Image
from collections import Counter

sys.stdout.reconfigure(line_buffering=True)

# Global compressors (reuse)
try:
    import zstandard
    ZSTD = zstandard.ZstdCompressor(level=19)
    ZSTD_FAST = zstandard.ZstdCompressor(level=3)
    HAS_ZSTD = True
except ImportError:
    ZSTD = None
    ZSTD_FAST = None
    HAS_ZSTD = False

try:
    import brotli
    HAS_BROTLI = True
except ImportError:
    HAS_BROTLI = False

import lzma

# ================================================================
# PREDICTION FILTERS (vectorized)
# ================================================================

def med_predict(plane):
    """MED (median edge detector) prediction. Returns int16 residuals."""
    h, w = plane.shape
    p = plane.astype(np.int16)
    res = np.zeros_like(p)
    res[0, 0] = p[0, 0]
    if w > 1:
        res[0, 1:] = p[0, 1:] - p[0, :-1]
    if h > 1:
        res[1:, 0] = p[1:, 0] - p[:-1, 0]
    if h > 1 and w > 1:
        a = p[1:, :-1]
        b = p[:-1, 1:]
        c = p[:-1, :-1]
        mn = np.minimum(a, b)
        mx = np.maximum(a, b)
        pred = np.where(c >= mx, mn, np.where(c <= mn, mx, a + b - c))
        res[1:, 1:] = p[1:, 1:] - pred
    return res

def paeth_predict(plane):
    """Paeth prediction (as used in PNG)."""
    h, w = plane.shape
    p = plane.astype(np.int16)
    res = np.zeros_like(p)
    res[0, 0] = p[0, 0]
    if w > 1:
        res[0, 1:] = p[0, 1:] - p[0, :-1]
    if h > 1:
        res[1:, 0] = p[1:, 0] - p[:-1, 0]
    if h > 1 and w > 1:
        a = p[1:, :-1]
        b = p[:-1, 1:]
        c = p[:-1, :-1]
        pp = a + b - c
        pa = np.abs(pp - a)
        pb = np.abs(pp - b)
        pc = np.abs(pp - c)
        pred = np.where((pa <= pb) & (pa <= pc), a, np.where(pb <= pc, b, c))
        res[1:, 1:] = p[1:, 1:] - pred
    return res

def left_predict(plane):
    """Left (Sub) prediction — simple horizontal delta."""
    p = plane.astype(np.int16)
    res = np.zeros_like(p)
    res[:, 0] = p[:, 0]
    if p.shape[1] > 1:
        res[:, 1:] = p[:, 1:] - p[:, :-1]
    return res

def up_predict(plane):
    """Up prediction — simple vertical delta."""
    p = plane.astype(np.int16)
    res = np.zeros_like(p)
    res[0, :] = p[0, :]
    if p.shape[0] > 1:
        res[1:, :] = p[1:, :] - p[:-1, :]
    return res

def wrap_residuals(res):
    """Wrap int16 residuals to uint8."""
    return ((res + 128) % 256).astype(np.uint8)

# ================================================================
# COLOR TRANSFORMS
# ================================================================

def grd_decorrelate(arr):
    """G, R-G, B-G color decorrelation."""
    g = arr[:,:,1]
    rg = ((arr[:,:,0].astype(np.int16) - arr[:,:,1].astype(np.int16)) % 256).astype(np.uint8)
    bg = ((arr[:,:,2].astype(np.int16) - arr[:,:,1].astype(np.int16)) % 256).astype(np.uint8)
    return [g, rg, bg]

def ycocg_decorrelate(arr):
    """YCoCg-R (reversible) color transform. Often better than G-RG-BG."""
    r = arr[:,:,0].astype(np.int16)
    g = arr[:,:,1].astype(np.int16)
    b = arr[:,:,2].astype(np.int16)
    co = (r - b) & 0xFF
    tmp = (b + (((r - b) + 256) >> 1)) & 0xFF  # careful with mod arithmetic
    cg = (g - tmp) & 0xFF
    y = (tmp + (((g - tmp) + 256) >> 1)) & 0xFF
    return [y.astype(np.uint8), co.astype(np.uint8), cg.astype(np.uint8)]

# ================================================================
# COMPRESSION BACKENDS
# ================================================================

def compress_zlib(data_bytes):
    return zlib.compress(data_bytes, 9)

def compress_zstd(data_bytes):
    if HAS_ZSTD:
        return ZSTD.compress(data_bytes)
    return zlib.compress(data_bytes, 9)

def compress_brotli(data_bytes):
    if HAS_BROTLI:
        return brotli.compress(data_bytes, quality=11)
    return zlib.compress(data_bytes, 9)

def compress_lzma(data_bytes):
    return lzma.compress(data_bytes, preset=9 | lzma.PRESET_EXTREME)

# ================================================================
# FULL-IMAGE ENCODERS
# ================================================================

def encode_pipeline(planes, predictor, compressor):
    """Generic: apply predictor to each plane, concatenate, compress."""
    residuals = []
    for plane in planes:
        res = wrap_residuals(predictor(plane))
        residuals.append(res.ravel())
    combined = np.concatenate(residuals)
    return compressor(combined.tobytes())

def encode_420_pipeline(arr, predictor, compressor):
    """4:2:0 chroma subsampling + predictor + compressor."""
    h, w = arr.shape[:2]
    planes = grd_decorrelate(arr)
    ph = h if h % 2 == 0 else h - 1
    pw = w if w % 2 == 0 else w - 1

    # Full-res luma
    luma_res = wrap_residuals(predictor(planes[0]))

    # Subsampled chroma
    chroma_res = []
    for ch_plane in planes[1:]:
        sub = ch_plane[:ph, :pw].reshape(ph//2, 2, pw//2, 2).mean(axis=(1,3)).astype(np.uint8)
        chroma_res.append(wrap_residuals(predictor(sub)))

    combined = np.concatenate([luma_res.ravel()] + [c.ravel() for c in chroma_res])
    return compressor(combined.tobytes())

def encode_palette(arr, compressor):
    """Global palette + packed indices."""
    h, w = arr.shape[:2]
    n = h * w

    # Quick reject for photo-like images
    if n > 500000:
        sample_idx = np.random.default_rng(42).choice(n, min(50000, n), replace=False)
        flat_sample = arr.reshape(-1, 3)[sample_idx]
        packed_sample = flat_sample[:, 0].astype(np.uint32) * 65536 + \
                        flat_sample[:, 1].astype(np.uint32) * 256 + \
                        flat_sample[:, 2].astype(np.uint32)
        if len(np.unique(packed_sample)) / len(packed_sample) > 0.5:
            raise ValueError("Too many colors")

    flat = arr.reshape(-1, 3)
    packed = flat[:, 0].astype(np.uint32) * 65536 + flat[:, 1].astype(np.uint32) * 256 + flat[:, 2].astype(np.uint32)
    unique_packed = np.unique(packed)
    n_colors = len(unique_packed)

    if n_colors > 65536:
        raise ValueError("Too many colors")

    color_to_idx = {int(c): i for i, c in enumerate(unique_packed)}
    if n_colors <= 256:
        indices = np.array([color_to_idx[int(c)] for c in packed], dtype=np.uint8)
    else:
        indices = np.array([color_to_idx[int(c)] for c in packed], dtype=np.uint16)

    header = struct.pack('<H', n_colors)
    palette_bytes = bytearray()
    for c in unique_packed:
        palette_bytes.extend([(int(c) >> 16) & 0xFF, (int(c) >> 8) & 0xFF, int(c) & 0xFF])

    compressed_idx = compressor(indices.tobytes())
    return header + bytes(palette_bytes) + compressed_idx


# ================================================================
# PNG-STYLE ROW FILTERS
# ================================================================

def png_filter_rows(plane):
    """Apply PNG-style adaptive row filtering to a plane.
    Returns bytes with filter-type prefix per row."""
    h, w = plane.shape
    p = plane.astype(np.int16)
    rows = []
    for r in range(h):
        row = p[r]
        # Filter 0: None
        f0 = row.copy()
        # Filter 1: Sub (left delta)
        f1 = np.zeros_like(row)
        f1[0] = row[0]
        if w > 1:
            f1[1:] = row[1:] - row[:-1]
        # Filter 2: Up (above delta)
        if r > 0:
            f2 = row - p[r-1]
        else:
            f2 = row.copy()
        # Filter 3: Average
        if r > 0:
            above = p[r-1]
        else:
            above = np.zeros_like(row)
        left = np.zeros_like(row)
        if w > 1:
            left[1:] = p[r, :-1]
        f3 = row - ((above + left) >> 1)

        # Pick filter with minimum absolute sum (PNG minimum sum heuristic)
        candidates = [(0, f0), (1, f1), (2, f2), (3, f3)]
        best_filter, best_data = min(candidates, key=lambda x: int(np.sum(np.abs(x[1]))))
        rows.append(bytes([best_filter]) + ((best_data % 256).astype(np.uint8)).tobytes())
    return b''.join(rows)

def encode_pngrow(arr, compressor):
    """PNG-style row filter on decorrelated planes + compressor."""
    planes = grd_decorrelate(arr)
    all_filtered = []
    for plane in planes:
        all_filtered.append(png_filter_rows(plane))
    return compressor(b''.join(all_filtered))

def encode_pngrow_raw(arr, compressor):
    """PNG-style row filter on raw RGB interleaved scanlines + compressor."""
    h, w = arr.shape[:2]
    # Interleave RGB channels per scanline (like PNG does it)
    rows = []
    for r in range(h):
        scanline = arr[r].ravel()  # R,G,B,R,G,B,...
        prev = arr[r-1].ravel() if r > 0 else np.zeros_like(scanline)
        # Filter 0: None
        f0 = scanline.astype(np.int16)
        # Filter 1: Sub (left, stride 3 for RGB)
        f1 = scanline.astype(np.int16).copy()
        if w > 1:
            f1[3:] = scanline[3:].astype(np.int16) - scanline[:-3].astype(np.int16)
        # Filter 2: Up
        f2 = scanline.astype(np.int16) - prev.astype(np.int16)

        candidates = [(0, f0), (1, f1), (2, f2)]
        best_filter, best_data = min(candidates, key=lambda x: int(np.sum(np.abs(x[1]))))
        rows.append(bytes([best_filter]) + ((best_data % 256).astype(np.uint8)).tobytes())
    return compressor(b''.join(rows))

# ================================================================
# THE CODEC — try all full-image strategies, pick smallest
# ================================================================

def mosaic_encode(arr):
    """Try multiple full-image encoding strategies. Return smallest."""
    h, w = arr.shape[:2]
    n_pixels = h * w

    results = {}

    # Compute color planes once
    grd_planes = grd_decorrelate(arr)

    # MED + GRD + zstd (primary strategy)
    if HAS_ZSTD:
        results['med_grd_zstd'] = len(encode_pipeline(grd_planes, med_predict, compress_zstd))

    # MED + GRD + zlib
    results['med_grd_zlib'] = len(encode_pipeline(grd_planes, med_predict, compress_zlib))

    # 4:2:0 + MED + zstd
    if h >= 4 and w >= 4 and HAS_ZSTD:
        try:
            results['420_med_zstd'] = len(encode_420_pipeline(arr, med_predict, compress_zstd))
        except Exception:
            pass

    # 4:2:0 + MED + zlib
    if h >= 4 and w >= 4:
        try:
            results['420_med_zlib'] = len(encode_420_pipeline(arr, med_predict, compress_zlib))
        except Exception:
            pass

    # 4:2:0 + brotli (potentially best for medium photos)
    if h >= 4 and w >= 4 and HAS_BROTLI and n_pixels < 4000000:
        try:
            results['420_med_brotli'] = len(encode_420_pipeline(arr, med_predict, compress_brotli))
        except Exception:
            pass

    # 4:2:0 + lzma (potentially best ratio for small-medium photos)
    if h >= 4 and w >= 4 and n_pixels < 2200000:
        try:
            results['420_med_lzma'] = len(encode_420_pipeline(arr, med_predict, compress_lzma))
        except Exception:
            pass

    # Paeth + GRD + zstd
    if HAS_ZSTD:
        results['paeth_grd_zstd'] = len(encode_pipeline(grd_planes, paeth_predict, compress_zstd))

    # Left (Sub) + GRD + zstd
    if HAS_ZSTD:
        results['left_grd_zstd'] = len(encode_pipeline(grd_planes, left_predict, compress_zstd))

    # Up + GRD + zstd
    if HAS_ZSTD:
        results['up_grd_zstd'] = len(encode_pipeline(grd_planes, up_predict, compress_zstd))

    # YCoCg-R variants (if helpful)
    if HAS_ZSTD:
        try:
            ycocg_planes = ycocg_decorrelate(arr)
            results['med_ycocg_zstd'] = len(encode_pipeline(ycocg_planes, med_predict, compress_zstd))
        except Exception:
            pass

    # Palette (only for low-color images)
    if HAS_ZSTD:
        try:
            results['pal_zstd'] = len(encode_palette(arr, compress_zstd))
        except Exception:
            pass
    try:
        results['pal_zlib'] = len(encode_palette(arr, compress_zlib))
    except Exception:
        pass

    # Raw bytes + various compressors (surprisingly competitive on structured content)
    raw_bytes = arr.tobytes()
    if HAS_ZSTD and n_pixels < 4000000:
        results['raw_zstd'] = len(compress_zstd(raw_bytes))
    if HAS_BROTLI and n_pixels < 4000000:
        results['raw_brotli'] = len(compress_brotli(raw_bytes))
    if n_pixels < 2200000:
        results['raw_lzma'] = len(compress_lzma(raw_bytes))

    # Brotli alternatives (excellent on small structured data)
    if HAS_BROTLI and n_pixels < 4000000:
        results['med_grd_brotli'] = len(encode_pipeline(grd_planes, med_predict, compress_brotli))
        try:
            results['pal_brotli'] = len(encode_palette(arr, compress_brotli))
        except Exception:
            pass
        results['paeth_grd_brotli'] = len(encode_pipeline(grd_planes, paeth_predict, compress_brotli))

    # LZMA alternatives (best ratio, slowest — only for small images)
    if n_pixels < 2200000:
        results['med_grd_lzma'] = len(encode_pipeline(grd_planes, med_predict, compress_lzma))
        try:
            results['pal_lzma'] = len(encode_palette(arr, compress_lzma))
        except Exception:
            pass

    # PNG-style: per-row filter byte + filtered scanlines + zstd
    # This captures the same spatial redundancy that makes PNG/WebP good on screenshots
    # Only for images small enough that per-row Python loop is tolerable
    if HAS_ZSTD and n_pixels < 4000000:
        try:
            results['pngrow_grd_zstd'] = len(encode_pngrow(arr, compress_zstd))
        except Exception:
            pass

        # PNG-style with raw RGB (no color transform — sometimes better on palette content)
        try:
            results['pngrow_raw_zstd'] = len(encode_pngrow_raw(arr, compress_zstd))
        except Exception:
            pass

    # --- Block-adaptive encoding (only if worth it) ---
    # Do block-level if image has mixed content (check variance of block energies)
    if n_pixels >= 100000 and n_pixels < 8000000:  # Skip very small and very large
        block_result = block_adaptive_encode(arr)
        if block_result is not None:
            results['mosaic'] = block_result

    if not results:
        results['raw_zlib'] = len(compress_zlib(arr.tobytes()))

    best_name = min(results, key=results.get)
    return results[best_name], best_name, results


# ================================================================
# BLOCK-ADAPTIVE PATH (for mixed content)
# ================================================================

BLOCK_SIZE = 64

def block_adaptive_encode(arr):
    """Block-level encoding: classify blocks, merge strips, encode per-strip."""
    h, w = arr.shape[:2]
    bh = (h + BLOCK_SIZE - 1) // BLOCK_SIZE
    bw = (w + BLOCK_SIZE - 1) // BLOCK_SIZE

    # Quick content analysis: sample some blocks
    block_energies = []
    block_colors = []
    for by in range(0, bh, max(1, bh // 8)):
        for bx in range(0, bw, max(1, bw // 8)):
            y0, y1 = by * BLOCK_SIZE, min((by + 1) * BLOCK_SIZE, h)
            x0, x1 = bx * BLOCK_SIZE, min((bx + 1) * BLOCK_SIZE, w)
            block = arr[y0:y1, x0:x1]
            # Quick energy
            res = med_predict(block[:,:,1])  # just green channel
            energy = float(np.mean(np.abs(res)))
            block_energies.append(energy)
            # Quick color count
            packed = block[:,:,0].astype(np.uint32).ravel() * 65536 + \
                     block[:,:,1].astype(np.uint32).ravel() * 256 + \
                     block[:,:,2].astype(np.uint32).ravel()
            n_unique = len(np.unique(packed))
            block_colors.append(n_unique)

    # Only do block-adaptive if there's content diversity
    if len(block_energies) >= 4:
        energy_std = np.std(block_energies)
        color_std = np.std(block_colors)
        # If all blocks look similar, skip block-adaptive
        if energy_std < 2.0 and color_std < 50:
            return None

    # Actually encode block by block
    total_data = 0
    header_size = 12  # magic + dims + block config
    strip_overhead_per = 5  # type + count + size per strip

    n_strips = 0
    for by in range(bh):
        bx = 0
        while bx < bw:
            y0, y1 = by * BLOCK_SIZE, min((by + 1) * BLOCK_SIZE, h)
            x0, x1 = bx * BLOCK_SIZE, min((bx + 1) * BLOCK_SIZE, w)
            block = arr[y0:y1, x0:x1]

            # Classify
            packed = block[:,:,0].astype(np.uint32).ravel() * 65536 + \
                     block[:,:,1].astype(np.uint32).ravel() * 256 + \
                     block[:,:,2].astype(np.uint32).ravel()
            n_unique = len(np.unique(packed))
            res = med_predict(block[:,:,1])
            energy = float(np.mean(np.abs(res)))

            # Merge adjacent blocks of similar type
            strip_end = bx + 1
            is_flat = n_unique <= 4 and energy < 2.0
            is_photo = energy >= 5.0

            while strip_end < bw:
                nx0 = strip_end * BLOCK_SIZE
                nx1 = min((strip_end + 1) * BLOCK_SIZE, w)
                next_block = arr[y0:y1, nx0:nx1]
                np2 = next_block[:,:,0].astype(np.uint32).ravel() * 65536 + \
                      next_block[:,:,1].astype(np.uint32).ravel() * 256 + \
                      next_block[:,:,2].astype(np.uint32).ravel()
                nu2 = len(np.unique(np2))
                nr = med_predict(next_block[:,:,1])
                ne = float(np.mean(np.abs(nr)))
                next_flat = nu2 <= 4 and ne < 2.0
                next_photo = ne >= 5.0
                if is_flat == next_flat and is_photo == next_photo:
                    strip_end += 1
                else:
                    break

            # Encode the strip
            strip_x0 = bx * BLOCK_SIZE
            strip_x1 = min(strip_end * BLOCK_SIZE, w)
            strip = arr[y0:y1, strip_x0:strip_x1]

            # Try a few strategies
            best_sz = float('inf')
            if is_flat:
                # Check if solid
                sp = strip[:,:,0].astype(np.uint32).ravel() * 65536 + \
                     strip[:,:,1].astype(np.uint32).ravel() * 256 + \
                     strip[:,:,2].astype(np.uint32).ravel()
                su = np.unique(sp)
                if len(su) == 1:
                    best_sz = 4  # solid: just the color
                else:
                    # Palette
                    try:
                        pal_data = encode_palette(strip, compress_zstd if HAS_ZSTD else compress_zlib)
                        best_sz = min(best_sz, len(pal_data))
                    except Exception:
                        pass

            # MED + zstd
            if HAS_ZSTD:
                grd = grd_decorrelate(strip)
                data = encode_pipeline(grd, med_predict, compress_zstd)
                best_sz = min(best_sz, len(data))

            # 4:2:0 + zstd (for photo strips)
            if is_photo and strip.shape[0] >= 4 and strip.shape[1] >= 4 and HAS_ZSTD:
                try:
                    data = encode_420_pipeline(strip, med_predict, compress_zstd)
                    best_sz = min(best_sz, len(data))
                except Exception:
                    pass

            if best_sz == float('inf'):
                grd = grd_decorrelate(strip)
                data = encode_pipeline(grd, med_predict, compress_zlib)
                best_sz = len(data)

            total_data += best_sz
            n_strips += 1
            bx = strip_end

    # Total with overhead
    map_bits = bh * bw * 2  # 2 bits per block for type
    total = header_size + n_strips * strip_overhead_per + (map_bits + 7) // 8 + total_data
    return total


# ================================================================
# BENCHMARK
# ================================================================

def png_sz(arr):
    buf = io.BytesIO()
    Image.fromarray(arr).save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()

def webp_ll_sz(arr):
    buf = io.BytesIO()
    Image.fromarray(arr).save(buf, format='WEBP', lossless=True)
    return buf.tell()

def tic_sz(arr):
    import subprocess
    h, w = arr.shape[:2]
    arr.flatten().astype(np.uint8).tofile('/tmp/tic_t.raw')
    try:
        r = subprocess.run(
            ['./tic_adaptive', 'bench', '/tmp/tic_t.raw', str(w), str(h), '1'],
            capture_output=True, text=True, timeout=300,
            cwd='/Users/e/Documents/GitHub/math/04-computation/tic'
        )
        for l in r.stdout.split('\n'):
            if 'Compressed:' in l:
                return int(l.split()[1])
    except Exception:
        pass
    return 0


def run_benchmark():
    print("=" * 135)
    print("  TIC-Mosaic v4: Full-Image First, Block-Adaptive Refinement")
    print("  Multiple predictors × color transforms × compressors | Block-adaptive for mixed content only")
    print("=" * 135)

    dirs = ['/Users/e/Documents/GitHub/math/inbox/processed/2026-04-02/new/',
            '/Users/e/Documents/GitHub/math/test_screenshots/']

    imgs = {}
    for d in dirs:
        if not os.path.isdir(d):
            continue
        for f in sorted(os.listdir(d)):
            if f.startswith('.'):
                continue
            if f.lower().endswith(('.mp4', '.gif', '.cr2', '.dng', '.heic', '.webp')):
                continue
            path = os.path.join(d, f)
            try:
                arr = np.array(Image.open(path).convert('RGB'))
                if arr.shape[0] < 100 or arr.shape[1] < 100:
                    continue
                label = f[:42]
                if f.lower().endswith('.jpg'):
                    label = 'PHOTO:' + label
                imgs[label] = arr
            except Exception:
                pass

    n = len(imgs)
    print(f"\nLoaded {n} images\n")
    hdr = f"{'#':>2} {'Image':<44} {'WxH':>11} {'PNG':>10} {'TIC-A':>10} {'Mos4':>10} {'WebP-LL':>10} {'Best':>8} {'Method':>16} {'vs PNG':>7} {'vs WebP':>7}"
    print(hdr)
    print("-" * 150)

    wins = Counter()
    total = {'PNG': 0, 'TIC-A': 0, 'Mos4': 0, 'WebP-LL': 0}

    for idx, (name, arr) in enumerate(sorted(imgs.items()), 1):
        h, w = arr.shape[:2]

        t0 = time.time()
        mosaic_sz, method, all_results = mosaic_encode(arr)
        encode_sec = time.time() - t0

        ps = png_sz(arr)
        ts = tic_sz(arr)
        ws = webp_ll_sz(arr)

        total['PNG'] += ps
        total['TIC-A'] += ts
        total['Mos4'] += mosaic_sz
        total['WebP-LL'] += ws

        candidates = {'PNG': ps, 'Mos4': mosaic_sz, 'WebP-LL': ws}
        if ts > 0:
            candidates['TIC-A'] = ts
        best = min(candidates, key=candidates.get)
        wins[best] += 1

        ts_str = f"{ts:>10,}" if ts > 0 else f"{'N/A':>10}"
        vs_png = f"{mosaic_sz/ps*100:5.1f}%" if ps > 0 else "N/A"
        vs_webp = f"{mosaic_sz/ws*100:5.1f}%" if ws > 0 else "N/A"

        print(f"{idx:>2} {name:<44} {w}x{h:>4} {ps:>10,} {ts_str} {mosaic_sz:>10,} {ws:>10,} {best:>8} {method:>16} {vs_png:>7} {vs_webp:>7}  ({encode_sec:.1f}s)")

    print("-" * 150)
    tp, tt, tm, tw = total['PNG'], total['TIC-A'], total['Mos4'], total['WebP-LL']
    print(f"   {'TOTALS':<44} {'':>11} {tp:>10,} {tt:>10,} {tm:>10,} {tw:>10,}")
    print(f"\n   vs PNG:  Mos4 = {tm/tp*100:.1f}%  |  WebP-LL = {tw/tp*100:.1f}%  |  TIC-A = {tt/tp*100:.1f}%")
    print(f"   vs WebP: Mos4 = {tm/tw*100:.1f}%")
    print(f"\n   WINS: {dict(wins)} (out of {n} images)")

    # Show where we lost
    print(f"\n   Legend: method = best internal encoding strategy chosen")


if __name__ == '__main__':
    run_benchmark()
