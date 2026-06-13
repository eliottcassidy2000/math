#!/usr/bin/env python3
"""
TIC-Mosaic: block-adaptive codec that classifies image regions
and applies the optimal compression strategy per block.

opus-2026-04-02-S25

The key insight: real images are COMPOSITES.
- A screenshot has flat UI + anti-aliased text + embedded photos
- A photo editor has a photo center + dark UI panels + slider widgets
- A web page has text blocks + images + gradients + icons

Instead of one global codec, we:
1. Divide into 32×32 blocks
2. Classify each block by content type
3. Apply the best encoder per block
4. Store a strategy map (2 bits/block ≈ 500 bytes overhead for 1080p)

This should combine TIC's photo strength with deflate's UI strength.
"""
import sys, io, os, zlib, struct, math, time
import numpy as np
from PIL import Image
from collections import Counter

sys.stdout.reconfigure(line_buffering=True)

BLOCK = 32  # Block size

# ================================================================
# BLOCK CLASSIFICATION
# ================================================================

def classify_block(block):
    """Classify a block into one of 4 types based on its content.
    Returns: ('flat'|'ui'|'gradient'|'photo', features_dict)
    """
    h, w = block.shape[:2]
    n = h * w

    # Feature 1: unique colors
    flat = block.reshape(-1, 3)
    unique = len(set(map(tuple, flat)))

    # Feature 2: gradient energy (MED residual magnitude)
    energy = 0
    for ch in range(3):
        p = block[:,:,ch].astype(int)
        for r in range(h):
            for c in range(w):
                a = p[r, c-1] if c > 0 else 0
                b = p[r-1, c] if r > 0 else 0
                cv = p[r-1, c-1] if (r > 0 and c > 0) else 0
                mn, mx = min(a,b), max(a,b)
                pred = mn if cv >= mx else (mx if cv <= mn else a + b - cv)
                energy += abs(int(p[r,c]) - pred)
    energy_per_pixel = energy / (n * 3)

    # Feature 3: is it completely flat?
    if unique == 1:
        return 'flat', {'unique': unique, 'energy': energy_per_pixel}

    # Feature 4: unique color ratio
    color_ratio = unique / n

    # Classification rules:
    if unique <= 4 and energy_per_pixel < 2:
        return 'flat', {'unique': unique, 'energy': energy_per_pixel}
    elif color_ratio < 0.05 and energy_per_pixel < 10:
        return 'ui', {'unique': unique, 'energy': energy_per_pixel, 'ratio': color_ratio}
    elif energy_per_pixel < 5:
        return 'gradient', {'unique': unique, 'energy': energy_per_pixel}
    else:
        return 'photo', {'unique': unique, 'energy': energy_per_pixel}

# ================================================================
# PER-BLOCK ENCODERS
# ================================================================

def med_predict_block(block_ch):
    """MED predict a single-channel block, return residuals as bytes."""
    h, w = block_ch.shape
    res = np.zeros((h, w), dtype=np.int16)
    p = block_ch.astype(np.int16)
    for r in range(h):
        for c in range(w):
            a = p[r, c-1] if c > 0 else 0
            b = p[r-1, c] if r > 0 else 0
            cv = p[r-1, c-1] if (r > 0 and c > 0) else 0
            mn, mx = min(a,b), max(a,b)
            pred = mn if cv >= mx else (mx if cv <= mn else a + b - cv)
            res[r, c] = p[r, c] - pred
    return ((res + 128) % 256).astype(np.uint8)

def encode_flat(block):
    """Flat block: store just the color(s). 3-12 bytes."""
    flat = block.reshape(-1, 3)
    unique = list(set(map(tuple, flat)))
    if len(unique) == 1:
        return bytes([0, unique[0][0], unique[0][1], unique[0][2]])
    # Few colors: store palette + indices
    palette = unique[:16]
    idx_map = {c: i for i, c in enumerate(palette)}
    bits = len(palette).bit_length()
    data = bytearray([len(palette)])
    for c in palette:
        data.extend(c)
    # Pack indices
    for pixel in flat:
        t = tuple(pixel)
        data.append(idx_map.get(t, 0))
    return bytes(data)

def encode_deflate(block):
    """UI/text block: MED predict + decorrelate + deflate. Best for repeated patterns."""
    h, w = block.shape[:2]
    # G, R-G, B-G decorrelation
    g = block[:,:,1]
    rg = ((block[:,:,0].astype(int) - block[:,:,1].astype(int)) & 0xFF).astype(np.uint8)
    bg = ((block[:,:,2].astype(int) - block[:,:,1].astype(int)) & 0xFF).astype(np.uint8)
    # MED predict each
    res = np.concatenate([med_predict_block(g), med_predict_block(rg), med_predict_block(bg)])
    return zlib.compress(res.tobytes(), 9)

def encode_gr_residuals(residuals):
    """Golomb-Rice encode residual bytes. Returns compressed bytes."""
    # Simple GR: zigzag + adaptive k
    out = bytearray()
    avg = 4.0
    for val in residuals:
        # Unwrap from [0,255] to [-128,127]
        v = int(val) - 128
        # Zigzag
        zz = 2 * v if v >= 0 else 2 * (-v) - 1
        # k from average
        k = 0
        if avg >= 0.5:
            k = max(0, min(12, round(math.log2(avg / 0.693))))
        # Encode: store q and r as bytes (simplified — real GR is bitwise)
        q = zz >> k
        r = zz & ((1 << k) - 1)
        # Pack: 1 byte for q (capped), k bits for r
        if q > 254:
            out.append(255)
            out.extend(struct.pack('<H', zz))
        else:
            out.append(q)
            if k > 0:
                out.extend(r.to_bytes((k + 7) // 8, 'little'))
        avg = 0.9 * avg + 0.1 * abs(v)
    return bytes(out)

def encode_photo(block):
    """Photo block: G-RG-BG decorrelation + MED + GR encoding."""
    g = block[:,:,1]
    rg = ((block[:,:,0].astype(int) - block[:,:,1].astype(int)) & 0xFF).astype(np.uint8)
    bg = ((block[:,:,2].astype(int) - block[:,:,1].astype(int)) & 0xFF).astype(np.uint8)
    res_g = med_predict_block(g).flatten()
    res_rg = med_predict_block(rg).flatten()
    res_bg = med_predict_block(bg).flatten()
    # Try both GR and deflate, pick smaller
    gr_data = encode_gr_residuals(np.concatenate([res_g, res_rg, res_bg]))
    zl_data = zlib.compress(np.concatenate([res_g, res_rg, res_bg]).tobytes(), 9)
    return gr_data if len(gr_data) <= len(zl_data) else zl_data

def encode_block(block, strategy):
    """Encode a block with the given strategy."""
    if strategy == 'flat':
        return encode_flat(block)
    elif strategy == 'ui':
        return encode_deflate(block)
    elif strategy == 'gradient':
        return encode_deflate(block)  # deflate works well for smooth gradients too
    else:  # photo
        return encode_photo(block)

# ================================================================
# MOSAIC CODEC
# ================================================================

def mosaic_encode(arr):
    """Encode image with block-adaptive strategy selection."""
    h, w = arr.shape[:2]
    bh = (h + BLOCK - 1) // BLOCK
    bw = (w + BLOCK - 1) // BLOCK

    strategy_map = []
    block_data = []
    type_counts = Counter()

    for by in range(bh):
        for bx in range(bw):
            y0, y1 = by * BLOCK, min((by + 1) * BLOCK, h)
            x0, x1 = bx * BLOCK, min((bx + 1) * BLOCK, w)
            block = arr[y0:y1, x0:x1]

            # Classify
            btype, features = classify_block(block)
            type_counts[btype] += 1

            # Try ALL strategies, pick smallest
            results = {}
            for strat in ['flat', 'ui', 'photo']:
                try:
                    data = encode_block(block, strat)
                    results[strat] = data
                except:
                    pass

            best_strat = min(results, key=lambda k: len(results[k]))
            best_data = results[best_strat]

            strat_id = {'flat': 0, 'ui': 1, 'gradient': 2, 'photo': 3}[best_strat]
            strategy_map.append(strat_id)
            block_data.append(best_data)

    # Pack: header + strategy map + block sizes + block data
    header = struct.pack('<4sHHHH', b'TICM', w, h, bw, bh)
    # Strategy map: 2 bits per block, packed into bytes
    map_bytes = bytearray()
    for i in range(0, len(strategy_map), 4):
        byte = 0
        for j in range(4):
            if i + j < len(strategy_map):
                byte |= (strategy_map[i + j] << (j * 2))
        map_bytes.append(byte)

    # Block sizes (4 bytes each)
    size_data = b''.join(struct.pack('<I', len(d)) for d in block_data)

    # All block data concatenated
    all_blocks = b''.join(block_data)

    total = len(header) + len(map_bytes) + len(size_data) + len(all_blocks)
    return total, type_counts, strategy_map

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
    r = subprocess.run(['./tic_adaptive', 'bench', '/tmp/tic_t.raw', str(w), str(h), '1'],
                       capture_output=True, text=True, timeout=300)
    for l in r.stdout.split('\n'):
        if 'Compressed:' in l: return int(l.split()[1])
    return 0

print("="*110)
print("  TIC-MOSAIC: Block-Adaptive Compression")
print("  Classifies each 32×32 block and applies optimal strategy")
print("="*110)

# Test on diverse images
dirs = ['/Users/e/Documents/GitHub/math/inbox/processed/2026-04-02/new/',
        '/Users/e/Documents/GitHub/math/test_screenshots/']

imgs = {}
for d in dirs:
    for f in sorted(os.listdir(d)):
        if f.startswith('.') or f.endswith(('.mp4','.gif','.cr2','.dng','.HEIC','.jpg')): continue
        if not f.lower().endswith('.png'): continue
        path = os.path.join(d, f)
        try:
            arr = np.array(Image.open(path).convert('RGB'))
            if arr.shape[0] < 100 or arr.shape[1] < 100: continue
            imgs[f[:35]] = arr
        except: pass

# Add a couple photos for comparison
for f in ['IMG_8669.jpg', 'IMG_1470.jpg']:
    path = os.path.join(dirs[0], f)
    if os.path.exists(path):
        imgs['PHOTO:'+f[:25]] = np.array(Image.open(path).convert('RGB'))

print(f"\nLoaded {len(imgs)} images\n")
print(f"{'#':>2} {'Image':<37} {'WxH':>10} {'PNG':>10} {'TIC':>10} {'Mosaic':>10} {'WebP-LL':>10} {'Best':>7} {'Composition'}")
print("-"*130)

wins = Counter()

for idx, (name, arr) in enumerate(sorted(imgs.items()), 1):
    h, w = arr.shape[:2]

    t0 = time.time()
    mosaic_sz, type_counts, strat_map = mosaic_encode(arr)
    mosaic_ms = (time.time() - t0) * 1000

    ps = png_sz(arr)
    ts = tic_sz(arr)
    ws = webp_ll_sz(arr)

    candidates = {'PNG': ps, 'TIC': ts, 'Mosaic': mosaic_sz, 'WebP-LL': ws}
    best = min(candidates, key=candidates.get)
    wins[best] += 1

    # Composition summary
    total_blocks = sum(type_counts.values())
    comp = ' '.join(f"{k}:{v}" for k, v in sorted(type_counts.items(), key=lambda x: -x[1]))

    print(f"{idx:>2} {name:<37} {w}x{h:>4} {ps:>10,} {ts:>10,} {mosaic_sz:>10,} {ws:>10,} {best:>7} {comp}")

print(f"\nWINS: {dict(wins)}")
print(f"\nMosaic approach: try ALL block strategies, pick smallest per block.")
print(f"If Mosaic beats TIC, it means block-adaptive helps.")
print(f"If Mosaic beats WebP-LL, we have a genuine advance.")
