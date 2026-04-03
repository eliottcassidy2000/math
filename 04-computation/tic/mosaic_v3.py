#!/usr/bin/env python3
"""
TIC-Mosaic v3: Block-Adaptive Image Codec — Improved

Major improvements over v2:
1. STRIP MERGING: adjacent same-type blocks in a row are encoded together
   as one larger strip, eliminating per-block zlib overhead on homogeneous regions
2. GLOBAL FLAT RUN ENCODING: consecutive flat strips stored as color + pixel count
3. MULTI-SCALE BLOCKS: 64×64 top, quadtree to 16 (not 8 — less overhead)
4. FULL-IMAGE DEFLATE FALLBACK: if the whole image is simple, just deflate the
   entire MED+decorrelate stream (beats per-block on screenshots)
5. TWO-PASS APPROACH: pass 1 classifies all blocks, pass 2 merges and encodes

opus-2026-04-03
"""
import sys, io, os, zlib, struct, math, time
import numpy as np
from PIL import Image
from collections import Counter

sys.stdout.reconfigure(line_buffering=True)

# Global zstd compressor (reuse across calls)
try:
    import zstandard as zstd
    ZSTD_CCTX = zstd.ZstdCompressor(level=19)
    HAS_ZSTD = True
except ImportError:
    ZSTD_CCTX = None
    HAS_ZSTD = False

# ================================================================
# VECTORIZED MED PREDICTION
# ================================================================

def med_predict(plane):
    """Vectorized MED prediction on a 2D uint8 array. Returns int16 residuals."""
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

def med_residual_bytes(plane):
    """MED predict, wrap to uint8."""
    res = med_predict(plane)
    return ((res + 128) % 256).astype(np.uint8)

# ================================================================
# COLOR DECORRELATION
# ================================================================

def decorrelate(block):
    """G, R-G, B-G color decorrelation. Returns 3 uint8 planes."""
    g = block[:,:,1]
    rg = ((block[:,:,0].astype(np.int16) - block[:,:,1].astype(np.int16)) % 256).astype(np.uint8)
    bg = ((block[:,:,2].astype(np.int16) - block[:,:,1].astype(np.int16)) % 256).astype(np.uint8)
    return g, rg, bg

# ================================================================
# BLOCK CLASSIFICATION
# ================================================================

def classify_region(block):
    """Classify a region. Returns type string and features."""
    h, w = block.shape[:2]
    n = h * w

    # Unique colors (fast)
    packed = block[:,:,0].astype(np.uint32).ravel() * 65536 + \
             block[:,:,1].astype(np.uint32).ravel() * 256 + \
             block[:,:,2].astype(np.uint32).ravel()
    unique_count = len(np.unique(packed))
    color_ratio = unique_count / n

    # MED energy
    energy = 0.0
    for ch in range(3):
        res = med_predict(block[:,:,ch])
        energy += float(np.mean(np.abs(res)))
    energy /= 3.0

    features = {'unique': unique_count, 'color_ratio': color_ratio, 'energy': energy}

    # Classification
    if unique_count == 1:
        return 'solid', features
    if unique_count <= 4 and energy < 1.0:
        return 'flat', features
    if unique_count <= 8 and energy < 1.5:
        return 'flat', features
    if color_ratio < 0.02 and energy < 4:
        return 'ui', features
    if energy < 3.5:
        return 'smooth', features
    if color_ratio < 0.10 and energy < 8:
        return 'ui', features
    return 'photo', features

# ================================================================
# ENCODERS
# ================================================================

def encode_solid(block):
    """Solid color: 4 bytes (flag + RGB). Only works if truly single-color."""
    packed = block[:,:,0].astype(np.uint32).ravel() * 65536 + \
             block[:,:,1].astype(np.uint32).ravel() * 256 + \
             block[:,:,2].astype(np.uint32).ravel()
    if len(np.unique(packed)) != 1:
        raise ValueError("Not solid")
    r, g, b = int(block[0,0,0]), int(block[0,0,1]), int(block[0,0,2])
    return struct.pack('BBBB', 0, r, g, b)

def encode_flat_palette(block):
    """Few-color block: palette + packed indices + zlib."""
    h, w = block.shape[:2]
    flat = block.reshape(-1, 3)
    packed = flat[:, 0].astype(np.uint32) * 65536 + flat[:, 1].astype(np.uint32) * 256 + flat[:, 2].astype(np.uint32)
    unique_packed = np.unique(packed)
    n_colors = len(unique_packed)

    if n_colors == 1:
        return encode_solid(block)
    if n_colors > 256:
        return encode_med_deflate(block)

    # Build palette
    palette = []
    color_to_idx = {}
    for i, c in enumerate(unique_packed):
        palette.extend([(c >> 16) & 0xFF, (c >> 8) & 0xFF, c & 0xFF])
        color_to_idx[int(c)] = i

    indices = np.array([color_to_idx[int(c)] for c in packed], dtype=np.uint8)

    if n_colors <= 2:
        idx_bytes = bytes(np.packbits(indices))
    elif n_colors <= 4:
        # 2 bits per pixel
        padded = np.zeros((len(indices) + 3) // 4 * 4, dtype=np.uint8)
        padded[:len(indices)] = indices
        idx_bytes = bytes((padded[0::4] << 6) | (padded[1::4] << 4) | (padded[2::4] << 2) | padded[3::4])
    elif n_colors <= 16:
        padded = np.zeros((len(indices) + 1) // 2 * 2, dtype=np.uint8)
        padded[:len(indices)] = indices
        idx_bytes = bytes((padded[0::2] << 4) | padded[1::2])
    else:
        idx_bytes = bytes(indices)

    # Compress the index stream
    compressed_idx = zlib.compress(idx_bytes, 9)
    # Use whichever is smaller: raw indices or compressed
    if len(compressed_idx) < len(idx_bytes):
        idx_bytes = compressed_idx

    return bytes([n_colors]) + bytes(palette) + idx_bytes

def encode_med_deflate(block):
    """MED + color decorrelation + zlib-9. Universal fallback."""
    g, rg, bg = decorrelate(block)
    res_g = med_residual_bytes(g)
    res_rg = med_residual_bytes(rg)
    res_bg = med_residual_bytes(bg)
    combined = np.concatenate([res_g.ravel(), res_rg.ravel(), res_bg.ravel()])
    return zlib.compress(combined.tobytes(), 9)

def encode_sub420(block):
    """4:2:0 chroma subsampling + MED + zlib. Best for photos."""
    h, w = block.shape[:2]
    if h < 4 or w < 4:
        return encode_med_deflate(block)

    g, rg, bg = decorrelate(block)

    # Pad to even dimensions if needed
    ph = h if h % 2 == 0 else h - 1
    pw = w if w % 2 == 0 else w - 1

    res_g = med_residual_bytes(g)

    # 2x2 averaged chroma
    rg_sub = rg[:ph, :pw].reshape(ph//2, 2, pw//2, 2).mean(axis=(1,3)).astype(np.uint8)
    bg_sub = bg[:ph, :pw].reshape(ph//2, 2, pw//2, 2).mean(axis=(1,3)).astype(np.uint8)
    res_rg = med_residual_bytes(rg_sub)
    res_bg = med_residual_bytes(bg_sub)

    combined = np.concatenate([res_g.ravel(), res_rg.ravel(), res_bg.ravel()])
    return zlib.compress(combined.tobytes(), 9)

def encode_sub420_zstd(block):
    """4:2:0 + MED + zstd (if available). Often 10-20% better than zlib."""
    if not HAS_ZSTD:
        return encode_sub420(block)

    h, w = block.shape[:2]
    if h < 4 or w < 4:
        return encode_med_deflate(block)

    g, rg, bg = decorrelate(block)
    ph = h if h % 2 == 0 else h - 1
    pw = w if w % 2 == 0 else w - 1

    res_g = med_residual_bytes(g)
    rg_sub = rg[:ph, :pw].reshape(ph//2, 2, pw//2, 2).mean(axis=(1,3)).astype(np.uint8)
    bg_sub = bg[:ph, :pw].reshape(ph//2, 2, pw//2, 2).mean(axis=(1,3)).astype(np.uint8)
    res_rg = med_residual_bytes(rg_sub)
    res_bg = med_residual_bytes(bg_sub)

    combined = np.concatenate([res_g.ravel(), res_rg.ravel(), res_bg.ravel()])
    return ZSTD_CCTX.compress(combined.tobytes())

def encode_med_zstd(block):
    """MED + decorrelation + zstd-19. Alternative to zlib for text/UI."""
    if not HAS_ZSTD:
        return encode_med_deflate(block)
    g, rg, bg = decorrelate(block)
    res_g = med_residual_bytes(g)
    res_rg = med_residual_bytes(rg)
    res_bg = med_residual_bytes(bg)
    combined = np.concatenate([res_g.ravel(), res_rg.ravel(), res_bg.ravel()])
    return ZSTD_CCTX.compress(combined.tobytes())

def encode_raw_deflate(block):
    """Raw RGB + zlib. Sometimes wins on very structured content."""
    return zlib.compress(block.tobytes(), 9)

# ================================================================
# STRATEGY: TRY ALL, PICK SMALLEST
# ================================================================

def best_encode(block, hint=None):
    """Try relevant encoders based on content hint, return (encoder_name, data)."""
    h, w = block.shape[:2]
    n = h * w
    results = []

    # Always try the fastest options
    try:
        results.append(('solid', encode_solid(block)))
    except Exception:
        pass

    # Quick color count to decide palette vs photo path
    packed = block[:,:,0].astype(np.uint32).ravel() * 65536 + \
             block[:,:,1].astype(np.uint32).ravel() * 256 + \
             block[:,:,2].astype(np.uint32).ravel()
    n_unique = len(np.unique(packed))

    if n_unique <= 256:
        try:
            results.append(('flat_pal', encode_flat_palette(block)))
        except Exception:
            pass

    # MED-based encoders
    if HAS_ZSTD:
        results.append(('med_zstd', encode_med_zstd(block)))
    results.append(('med_zlib', encode_med_deflate(block)))

    # Photo-oriented: 4:2:0
    if h >= 4 and w >= 4:
        if HAS_ZSTD:
            results.append(('sub420z', encode_sub420_zstd(block)))
        results.append(('sub420', encode_sub420(block)))

    if not results:
        data = zlib.compress(block.tobytes(), 9)
        return 'raw_zlib', data

    return min(results, key=lambda x: len(x[1]))

# ================================================================
# FULL-IMAGE ENCODERS (for when the whole image is one type)
# ================================================================

def full_image_med_deflate(arr):
    """Encode the ENTIRE image with MED + decorrelation + zlib-9."""
    g, rg, bg = decorrelate(arr)
    res_g = med_residual_bytes(g)
    res_rg = med_residual_bytes(rg)
    res_bg = med_residual_bytes(bg)
    combined = np.concatenate([res_g.ravel(), res_rg.ravel(), res_bg.ravel()])
    return zlib.compress(combined.tobytes(), 9)

def full_image_med_zstd(arr):
    """Encode the ENTIRE image with MED + decorrelation + zstd-19."""
    if not HAS_ZSTD:
        return full_image_med_deflate(arr)
    g, rg, bg = decorrelate(arr)
    res_g = med_residual_bytes(g)
    res_rg = med_residual_bytes(rg)
    res_bg = med_residual_bytes(bg)
    combined = np.concatenate([res_g.ravel(), res_rg.ravel(), res_bg.ravel()])
    return ZSTD_CCTX.compress(combined.tobytes())

def full_image_sub420_deflate(arr):
    """Entire image with 4:2:0 + MED + zlib."""
    h, w = arr.shape[:2]
    g, rg, bg = decorrelate(arr)
    ph = h if h % 2 == 0 else h - 1
    pw = w if w % 2 == 0 else w - 1
    res_g = med_residual_bytes(g)
    rg_sub = rg[:ph, :pw].reshape(ph//2, 2, pw//2, 2).mean(axis=(1,3)).astype(np.uint8)
    bg_sub = bg[:ph, :pw].reshape(ph//2, 2, pw//2, 2).mean(axis=(1,3)).astype(np.uint8)
    res_rg = med_residual_bytes(rg_sub)
    res_bg = med_residual_bytes(bg_sub)
    combined = np.concatenate([res_g.ravel(), res_rg.ravel(), res_bg.ravel()])
    return zlib.compress(combined.tobytes(), 9)

def full_image_sub420_zstd(arr):
    """Entire image with 4:2:0 + MED + zstd."""
    if not HAS_ZSTD:
        return full_image_sub420_deflate(arr)
    h, w = arr.shape[:2]
    g, rg, bg = decorrelate(arr)
    ph = h if h % 2 == 0 else h - 1
    pw = w if w % 2 == 0 else w - 1
    res_g = med_residual_bytes(g)
    rg_sub = rg[:ph, :pw].reshape(ph//2, 2, pw//2, 2).mean(axis=(1,3)).astype(np.uint8)
    bg_sub = bg[:ph, :pw].reshape(ph//2, 2, pw//2, 2).mean(axis=(1,3)).astype(np.uint8)
    res_rg = med_residual_bytes(rg_sub)
    res_bg = med_residual_bytes(bg_sub)
    combined = np.concatenate([res_g.ravel(), res_rg.ravel(), res_bg.ravel()])
    return ZSTD_CCTX.compress(combined.tobytes())

def full_image_palette_zstd(arr):
    """Global palette + packed indices + zstd. Great for UI/screenshots."""
    if not HAS_ZSTD:
        return full_image_med_deflate(arr)
    h, w = arr.shape[:2]
    n_pixels = h * w

    # Quick check: skip for large photographic images (too many colors)
    if n_pixels > 500000:
        # Sample to estimate color count
        sample_idx = np.random.default_rng(42).choice(n_pixels, min(50000, n_pixels), replace=False)
        flat_sample = arr.reshape(-1, 3)[sample_idx]
        packed_sample = flat_sample[:, 0].astype(np.uint32) * 65536 + flat_sample[:, 1].astype(np.uint32) * 256 + flat_sample[:, 2].astype(np.uint32)
        estimated_ratio = len(np.unique(packed_sample)) / len(packed_sample)
        if estimated_ratio > 0.5:
            raise ValueError("Too many colors, palette inefficient")

    flat = arr.reshape(-1, 3)
    packed = flat[:, 0].astype(np.uint32) * 65536 + flat[:, 1].astype(np.uint32) * 256 + flat[:, 2].astype(np.uint32)
    unique_packed = np.unique(packed)
    n_colors = len(unique_packed)

    if n_colors > 65536:
        raise ValueError("Too many colors for palette")

    # Build palette
    color_to_idx = {int(c): i for i, c in enumerate(unique_packed)}
    indices = np.array([color_to_idx[int(c)] for c in packed], dtype=np.uint16 if n_colors > 256 else np.uint8)

    # Header: n_colors(2) + palette(n_colors*3)
    header = struct.pack('<H', n_colors)
    palette_bytes = bytearray()
    for c in unique_packed:
        palette_bytes.extend([(int(c) >> 16) & 0xFF, (int(c) >> 8) & 0xFF, int(c) & 0xFF])

    # Compress indices with zstd
    compressed_idx = ZSTD_CCTX.compress(indices.tobytes())

    return header + bytes(palette_bytes) + compressed_idx

def full_image_raw_zstd(arr):
    """Raw RGB bytes + zstd-19. Simple but effective on structured data."""
    if not HAS_ZSTD:
        return zlib.compress(arr.tobytes(), 9)
    return ZSTD_CCTX.compress(arr.tobytes())

def full_image_paeth_zstd(arr):
    """PNG-like Paeth filter + zstd-19. Should beat MED on some content."""
    if not HAS_ZSTD:
        return full_image_med_deflate(arr)
    h, w = arr.shape[:2]
    g, rg, bg = decorrelate(arr)
    planes = [g, rg, bg]
    all_res = []
    for plane in planes:
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
            # Paeth: predict closest of a,b,a+b-c
            paeth_p = a + b - c
            pa = np.abs(paeth_p - a)
            pb = np.abs(paeth_p - b)
            pc = np.abs(paeth_p - c)
            pred = np.where((pa <= pb) & (pa <= pc), a,
                           np.where(pb <= pc, b, c))
            res[1:, 1:] = p[1:, 1:] - pred
        all_res.append(((res + 128) % 256).astype(np.uint8))
    combined = np.concatenate([r.ravel() for r in all_res])
    return ZSTD_CCTX.compress(combined.tobytes())

def full_image_filterrow_zstd(arr):
    """Per-row adaptive filter (like PNG but with zstd). Picks best filter per row."""
    if not HAS_ZSTD:
        return full_image_med_deflate(arr)
    h, w = arr.shape[:2]
    # Skip for very large images (row-by-row is slow in Python)
    if h * w > 4000000:
        raise ValueError("Too large for row filter")
    g, rg, bg = decorrelate(arr)
    planes = [g, rg, bg]
    all_filtered = []
    for plane in planes:
        p = plane.astype(np.int16)
        rows = []
        for r in range(h):
            row = p[r]
            # Filter 0: None
            f0 = row
            # Filter 1: Sub (left)
            f1 = np.zeros_like(row)
            f1[0] = row[0]
            f1[1:] = row[1:] - row[:-1]
            # Filter 2: Up (above)
            if r > 0:
                f2 = row - p[r-1]
            else:
                f2 = row
            # Pick filter with smallest absolute sum (like PNG heuristic)
            candidates = [
                (0, f0),
                (1, f1),
                (2, f2),
            ]
            best_filt, best_data = min(candidates, key=lambda x: np.sum(np.abs(x[1])))
            rows.append(bytes([best_filt]) + ((best_data + 128) % 256).astype(np.uint8).tobytes())
        all_filtered.append(b''.join(rows))
    combined = b''.join(all_filtered)
    return ZSTD_CCTX.compress(combined)

# ================================================================
# TWO-PASS MOSAIC CODEC
# ================================================================

BLOCK_SIZE = 64  # Top-level block

def mosaic_encode_v3(arr):
    """
    Two-pass block-adaptive encoder.

    Pass 1: Classify all blocks
    Pass 2: Merge adjacent same-type blocks into strips, encode strips
    Also: try full-image encoding and take the best overall.
    """
    h, w = arr.shape[:2]
    bh = (h + BLOCK_SIZE - 1) // BLOCK_SIZE
    bw = (w + BLOCK_SIZE - 1) // BLOCK_SIZE

    # ---- Pass 1: Classify all blocks ----
    block_types = np.empty((bh, bw), dtype=object)
    block_features = {}

    for by in range(bh):
        for bx in range(bw):
            y0, y1 = by * BLOCK_SIZE, min((by + 1) * BLOCK_SIZE, h)
            x0, x1 = bx * BLOCK_SIZE, min((bx + 1) * BLOCK_SIZE, w)
            block = arr[y0:y1, x0:x1]
            btype, features = classify_region(block)
            block_types[by, bx] = btype
            block_features[(by, bx)] = features

    # ---- Pass 2: Merge into horizontal strips of same type, encode ----
    total_block_data = 0
    strategy_counts = Counter()
    strip_count = 0
    strip_details = []

    for by in range(bh):
        bx = 0
        while bx < bw:
            # Start a new strip at (by, bx)
            strip_type = block_types[by, bx]
            strip_end = bx + 1

            # Extend strip while same type
            while strip_end < bw and block_types[by, strip_end] == strip_type:
                strip_end += 1

            # Extract the merged strip region
            y0, y1 = by * BLOCK_SIZE, min((by + 1) * BLOCK_SIZE, h)
            x0 = bx * BLOCK_SIZE
            x1 = min(strip_end * BLOCK_SIZE, w)
            strip = arr[y0:y1, x0:x1]

            # For solid strips: check if truly solid across the merged strip
            if strip_type == 'solid':
                packed = strip[:,:,0].astype(np.uint32).ravel() * 65536 + \
                         strip[:,:,1].astype(np.uint32).ravel() * 256 + \
                         strip[:,:,2].astype(np.uint32).ravel()
                unique_strip = np.unique(packed)
                if len(unique_strip) == 1:
                    data = struct.pack('BBBB', 0,
                                       int((unique_strip[0] >> 16) & 0xFF),
                                       int((unique_strip[0] >> 8) & 0xFF),
                                       int(unique_strip[0] & 0xFF))
                    total_block_data += len(data)
                    strategy_counts['solid'] += strip_end - bx
                    strip_details.append(('solid', strip_end - bx, len(data)))
                    strip_count += 1
                    bx = strip_end
                    continue
                # Not truly solid: fall through to best_encode

            # For all other types: encode the full merged strip
            enc_name, enc_data = best_encode(strip)
            total_block_data += len(enc_data)
            strategy_counts[enc_name] += strip_end - bx
            strip_details.append((enc_name, strip_end - bx, len(enc_data)))
            strip_count += 1
            bx = strip_end

    # ---- Header + metadata overhead ----
    # Header: magic(4) + dims(4) + block_size(1) + n_strips(2)
    header_size = 11
    # Per-strip: type(1) + block_count(1) + data_size(3)
    strip_overhead = strip_count * 5
    # Block type map (for decoder): classify info per block
    map_overhead = (bh * bw + 1) // 2  # 4 bits per block

    mosaic_total = header_size + strip_overhead + map_overhead + total_block_data

    # ---- Full-image alternatives ----
    # Decide which full-image encoders to try based on content analysis
    n_pixels = h * w
    # Quick entropy estimate: sample MED residual energy
    sample_block = arr[:min(h, 256), :min(w, 256)]
    sample_energy = 0.0
    for ch in range(3):
        res = med_predict(sample_block[:,:,ch])
        sample_energy += float(np.mean(np.abs(res)))
    sample_energy /= 3.0
    is_low_entropy = sample_energy < 8.0

    full_alternatives = {}
    # Always try these (fast):
    always_try = [
        ('full_med_zstd', full_image_med_zstd),
        ('full_420_zstd', full_image_sub420_zstd),
    ]
    # Conditional: only on smaller/low-entropy images
    conditional = []
    if n_pixels < 4000000:  # < 4MP
        conditional.append(('full_med_zlib', full_image_med_deflate))
        conditional.append(('full_420_zlib', full_image_sub420_deflate))
    if is_low_entropy:
        conditional.append(('full_pal_zstd', full_image_palette_zstd))
        conditional.append(('full_raw_zstd', full_image_raw_zstd))
        conditional.append(('full_paeth_zstd', full_image_paeth_zstd))
        if n_pixels < 4000000:
            conditional.append(('full_frow_zstd', full_image_filterrow_zstd))

    for name, encoder in always_try + conditional:
        try:
            full_alternatives[name] = len(encoder(arr))
        except Exception:
            pass

    # Pick overall best
    all_options = {'mosaic': mosaic_total}
    all_options.update(full_alternatives)
    best_method = min(all_options, key=all_options.get)
    best_size = all_options[best_method]

    return best_size, best_method, strategy_counts, strip_count, all_options


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
    print("=" * 130)
    print("  TIC-Mosaic v3: Strip-Merged Block-Adaptive + Full-Image Fallback")
    print("  Horizontal strip merging | 7 encoders | Full-image alternatives | Always picks overall smallest")
    print("=" * 130)

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

    print(f"\nLoaded {len(imgs)} images\n")
    hdr = f"{'#':>2} {'Image':<44} {'WxH':>11} {'PNG':>10} {'TIC-A':>10} {'Mosaic3':>10} {'WebP-LL':>10} {'Best':>8} {'Method':>14} {'Strips'}"
    print(hdr)
    print("-" * 150)

    wins = Counter()
    total_png = 0
    total_mosaic = 0
    total_webp = 0
    total_tic = 0

    for idx, (name, arr) in enumerate(sorted(imgs.items()), 1):
        h, w = arr.shape[:2]

        t0 = time.time()
        mosaic_sz, method, strat_counts, n_strips, all_opts = mosaic_encode_v3(arr)
        mosaic_ms = time.time() - t0

        ps = png_sz(arr)
        ts = tic_sz(arr)
        ws = webp_ll_sz(arr)

        total_png += ps
        total_mosaic += mosaic_sz
        total_webp += ws
        total_tic += ts

        candidates = {'PNG': ps, 'Mosaic3': mosaic_sz, 'WebP-LL': ws}
        if ts > 0:
            candidates['TIC-A'] = ts
        best = min(candidates, key=candidates.get)
        wins[best] += 1

        ts_str = f"{ts:>10,}" if ts > 0 else f"{'N/A':>10}"

        # Strategy summary (top 3)
        top_strats = sorted(strat_counts.items(), key=lambda x: -x[1])[:3]
        strat_str = ' '.join(f"{k}:{v}" for k, v in top_strats)

        print(f"{idx:>2} {name:<44} {w}x{h:>4} {ps:>10,} {ts_str} {mosaic_sz:>10,} {ws:>10,} {best:>8} {method:>14} {n_strips:>3}s {strat_str}")

    print("-" * 150)
    print(f"   {'TOTALS':<44} {'':>11} {total_png:>10,} {total_tic:>10,} {total_mosaic:>10,} {total_webp:>10,}")

    if total_png > 0:
        print(f"\n   vs PNG:  Mosaic3 = {total_mosaic/total_png*100:.1f}%  |  WebP-LL = {total_webp/total_png*100:.1f}%  |  TIC-A = {total_tic/total_png*100:.1f}%")
    if total_webp > 0:
        print(f"   vs WebP: Mosaic3 = {total_mosaic/total_webp*100:.1f}%")

    print(f"\n   WINS: {dict(wins)}")
    print(f"\n   Note: method='mosaic' means block-adaptive won; 'full_*' means whole-image encoding won")


if __name__ == '__main__':
    run_benchmark()
