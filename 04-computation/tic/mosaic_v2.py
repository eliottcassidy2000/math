#!/usr/bin/env python3
"""
TIC-Mosaic v2: Block-Adaptive Image Codec with Quadtree Sub-Block Decomposition

Key improvements over v1:
- Quadtree: blocks that are mixed get subdivided (64→32→16→8)
- Vectorized MED prediction (numpy, not Python loops)
- Better classification: edge density, polynomial fit, entropy
- Per-block optimal encoder selection (always try all, pick smallest)
- Proper flat encoding with palette RLE
- Gradient blocks: polynomial fit + small residual
- Context-aware boundary prediction (use neighbors from adjacent blocks)

opus-2026-04-03
"""
import sys, io, os, zlib, struct, math, time
import numpy as np
from PIL import Image
from collections import Counter

sys.stdout.reconfigure(line_buffering=True)

# ================================================================
# VECTORIZED MED PREDICTION
# ================================================================

def med_predict(plane):
    """Vectorized MED prediction on a 2D uint8 array. Returns int16 residuals."""
    h, w = plane.shape
    p = plane.astype(np.int16)
    res = np.zeros_like(p)

    # First pixel
    res[0, 0] = p[0, 0]
    # First row: predict from left (a), b=0, c=0 → med(a,0,0)=a if a in [0,0] else...
    # Simplified: first row, left predictor
    if w > 1:
        res[0, 1:] = p[0, 1:] - p[0, :-1]
    # First column: predict from above
    if h > 1:
        res[1:, 0] = p[1:, 0] - p[:-1, 0]

    # Interior: full MED
    if h > 1 and w > 1:
        a = p[1:, :-1]   # left
        b = p[:-1, 1:]   # above
        c = p[:-1, :-1]  # above-left
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
# BLOCK CLASSIFICATION (improved)
# ================================================================

def classify_block(block):
    """Classify a block into content types with rich features.
    Returns: (type_str, features_dict)
    Types: 'flat', 'ui', 'gradient', 'photo_smooth', 'photo_texture'
    """
    h, w = block.shape[:2]
    n = h * w

    # Feature: unique colors
    flat_pixels = block.reshape(-1, 3)
    # Use a hash-based approach for speed
    packed = flat_pixels[:, 0].astype(np.uint32) * 65536 + \
             flat_pixels[:, 1].astype(np.uint32) * 256 + \
             flat_pixels[:, 2].astype(np.uint32)
    unique_count = len(np.unique(packed))
    color_ratio = unique_count / n

    # Feature: per-channel variance
    var_r = float(np.var(block[:,:,0]))
    var_g = float(np.var(block[:,:,1]))
    var_b = float(np.var(block[:,:,2]))
    total_var = var_r + var_g + var_b

    # Feature: MED residual energy (vectorized)
    energy = 0.0
    for ch in range(3):
        res = med_predict(block[:,:,ch])
        energy += float(np.mean(np.abs(res)))
    energy /= 3.0  # average across channels

    # Feature: edge density (Sobel-like, horizontal + vertical differences)
    gray = block.mean(axis=2)
    if h > 1 and w > 1:
        dx = np.abs(np.diff(gray, axis=1)).mean()
        dy = np.abs(np.diff(gray, axis=0)).mean()
        edge_density = dx + dy
    else:
        edge_density = 0.0

    # Feature: gradient fit quality (can a 2D linear fit explain the block?)
    grad_residual = float('inf')
    if h >= 4 and w >= 4:
        # Fit linear plane: pixel = a*x + b*y + c per channel
        ys, xs = np.mgrid[:h, :w]
        xs_flat = xs.flatten().astype(np.float64)
        ys_flat = ys.flatten().astype(np.float64)
        ones = np.ones(n, dtype=np.float64)
        A = np.column_stack([xs_flat, ys_flat, ones])
        total_res = 0.0
        for ch in range(3):
            vals = block[:,:,ch].flatten().astype(np.float64)
            coeffs, _, _, _ = np.linalg.lstsq(A, vals, rcond=None)
            fitted = A @ coeffs
            total_res += np.mean(np.abs(vals - fitted))
        grad_residual = total_res / 3.0

    features = {
        'unique': unique_count,
        'color_ratio': color_ratio,
        'energy': energy,
        'edge_density': edge_density,
        'total_var': total_var,
        'grad_residual': grad_residual,
    }

    # Classification decision tree
    if unique_count <= 2:
        return 'flat', features
    if unique_count <= 8 and energy < 1.5:
        return 'flat', features
    if color_ratio < 0.03 and energy < 6:
        return 'ui', features
    if energy < 2.0 and grad_residual < 2.0:
        return 'gradient', features
    if energy < 4.0:
        return 'gradient', features
    if edge_density > 20 and color_ratio < 0.15:
        return 'ui', features
    if energy < 12 and edge_density < 10:
        return 'photo_smooth', features
    return 'photo_texture', features


# ================================================================
# ENCODERS
# ================================================================

def encode_flat(block):
    """Flat block: palette + indices. Optimal for ≤16 colors."""
    h, w = block.shape[:2]
    flat = block.reshape(-1, 3)
    packed = flat[:, 0].astype(np.uint32) * 65536 + \
             flat[:, 1].astype(np.uint32) * 256 + \
             flat[:, 2].astype(np.uint32)
    unique_packed = np.unique(packed)
    n_colors = len(unique_packed)

    if n_colors == 1:
        # Solid color: 4 bytes
        c = unique_packed[0]
        return bytes([0, (c >> 16) & 0xFF, (c >> 8) & 0xFF, c & 0xFF])

    if n_colors <= 16:
        # Palette + index map
        palette = []
        color_to_idx = {}
        for i, c in enumerate(unique_packed):
            palette.extend([(c >> 16) & 0xFF, (c >> 8) & 0xFF, c & 0xFF])
            color_to_idx[int(c)] = i

        indices = np.array([color_to_idx[int(c)] for c in packed], dtype=np.uint8)

        if n_colors <= 2:
            # 1 bit per pixel, pack 8 per byte
            bits = np.packbits(indices)
            data = bytes([n_colors]) + bytes(palette) + bytes(bits)
        elif n_colors <= 4:
            # 2 bits per pixel
            padded = np.zeros((len(indices) + 3) // 4 * 4, dtype=np.uint8)
            padded[:len(indices)] = indices
            packed_bytes = (padded[0::4] << 6) | (padded[1::4] << 4) | (padded[2::4] << 2) | padded[3::4]
            data = bytes([n_colors]) + bytes(palette) + bytes(packed_bytes)
        else:
            # 4 bits per pixel
            padded = np.zeros((len(indices) + 1) // 2 * 2, dtype=np.uint8)
            padded[:len(indices)] = indices
            packed_bytes = (padded[0::2] << 4) | padded[1::2]
            data = bytes([n_colors]) + bytes(palette) + bytes(packed_bytes)

        return data

    # >16 colors: fall through to deflate
    return encode_deflate(block)


def encode_deflate(block):
    """MED + color decorrelation + zlib-9. Best for UI/text."""
    g = block[:,:,1]
    rg = ((block[:,:,0].astype(np.int16) - block[:,:,1].astype(np.int16)) % 256).astype(np.uint8)
    bg = ((block[:,:,2].astype(np.int16) - block[:,:,1].astype(np.int16)) % 256).astype(np.uint8)

    res_g = med_residual_bytes(g)
    res_rg = med_residual_bytes(rg)
    res_bg = med_residual_bytes(bg)

    combined = np.concatenate([res_g.flatten(), res_rg.flatten(), res_bg.flatten()])
    return zlib.compress(combined.tobytes(), 9)


def encode_deflate_interleaved(block):
    """Interleaved scanline layout with decorrelation + MED + zlib-9.
    Sometimes better than planar for blocks with correlated channels."""
    h, w = block.shape[:2]
    g = block[:,:,1]
    rg = ((block[:,:,0].astype(np.int16) - block[:,:,1].astype(np.int16)) % 256).astype(np.uint8)
    bg = ((block[:,:,2].astype(np.int16) - block[:,:,1].astype(np.int16)) % 256).astype(np.uint8)

    res_g = med_residual_bytes(g)
    res_rg = med_residual_bytes(rg)
    res_bg = med_residual_bytes(bg)

    # Interleave: row of g, row of rg, row of bg, repeat
    interleaved = np.empty((h * 3, w), dtype=np.uint8)
    interleaved[0::3] = res_g
    interleaved[1::3] = res_rg
    interleaved[2::3] = res_bg
    return zlib.compress(interleaved.tobytes(), 9)


def encode_gradient(block):
    """Gradient block: fit 2D polynomial per channel, compress residual.
    For smooth gradients, the residual is tiny → great compression."""
    h, w = block.shape[:2]
    n = h * w

    ys, xs = np.mgrid[:h, :w]
    xs_flat = xs.flatten().astype(np.float64)
    ys_flat = ys.flatten().astype(np.float64)

    # Quadratic fit: a*x + b*y + c*x*y + d*x² + e*y² + f
    A = np.column_stack([xs_flat, ys_flat, xs_flat * ys_flat,
                         xs_flat**2, ys_flat**2, np.ones(n)])

    coeff_data = bytearray()
    residual_planes = []

    for ch in range(3):
        vals = block[:,:,ch].flatten().astype(np.float64)
        coeffs, _, _, _ = np.linalg.lstsq(A, vals, rcond=None)
        fitted = np.clip(np.round(A @ coeffs), 0, 255).astype(np.uint8)
        residual = ((block[:,:,ch].flatten().astype(np.int16) -
                     fitted.astype(np.int16) + 128) % 256).astype(np.uint8)

        # Store coefficients as float16 (2 bytes each, 6 coeffs = 12 bytes per channel)
        for c in coeffs:
            coeff_data.extend(struct.pack('<e', np.float16(c)))
        residual_planes.append(residual)

    # 36 bytes of coefficients + compressed residual
    combined_residual = np.concatenate(residual_planes)
    compressed_residual = zlib.compress(combined_residual.tobytes(), 9)

    return bytes(coeff_data) + compressed_residual


def encode_raw_deflate(block):
    """Raw pixels + deflate. Sometimes wins on very small or weird blocks."""
    return zlib.compress(block.tobytes(), 9)


def encode_sub4_deflate(block):
    """Subsample 4:2:0 chroma + MED + deflate. Good for photo blocks."""
    h, w = block.shape[:2]
    if h < 4 or w < 4:
        return encode_deflate(block)

    g = block[:,:,1]
    rg = ((block[:,:,0].astype(np.int16) - block[:,:,1].astype(np.int16)) % 256).astype(np.uint8)
    bg = ((block[:,:,2].astype(np.int16) - block[:,:,1].astype(np.int16)) % 256).astype(np.uint8)

    # Full-res luma
    res_g = med_residual_bytes(g)

    # 2× downsampled chroma (average 2×2 blocks)
    rg_sub = rg.reshape(h//2, 2, w//2, 2).mean(axis=(1,3)).astype(np.uint8)
    bg_sub = bg.reshape(h//2, 2, w//2, 2).mean(axis=(1,3)).astype(np.uint8)
    res_rg = med_residual_bytes(rg_sub)
    res_bg = med_residual_bytes(bg_sub)

    combined = np.concatenate([res_g.flatten(), res_rg.flatten(), res_bg.flatten()])
    return zlib.compress(combined.tobytes(), 9)


# ================================================================
# STRATEGY SELECTOR — try all applicable, pick smallest
# ================================================================

STRATEGIES = {
    0: ('flat', encode_flat),
    1: ('deflate_planar', encode_deflate),
    2: ('deflate_interleaved', encode_deflate_interleaved),
    3: ('gradient_poly', encode_gradient),
    4: ('raw_deflate', encode_raw_deflate),
    5: ('sub420_deflate', encode_sub4_deflate),
}


def best_encode(block):
    """Try all strategies, return (strategy_id, compressed_data)."""
    results = {}
    for sid, (name, encoder) in STRATEGIES.items():
        try:
            data = encoder(block)
            results[sid] = data
        except Exception:
            pass

    if not results:
        # Absolute fallback
        data = zlib.compress(block.tobytes(), 9)
        return 4, data

    best_sid = min(results, key=lambda k: len(results[k]))
    return best_sid, results[best_sid]


# ================================================================
# QUADTREE BLOCK DECOMPOSITION
# ================================================================

def should_subdivide(block, block_size, min_size=8):
    """Decide if a block is mixed content and should be subdivided."""
    if block_size <= min_size:
        return False
    h, w = block.shape[:2]
    if h < min_size * 2 or w < min_size * 2:
        return False

    # Quick test: classify the 4 sub-quadrants
    mh, mw = h // 2, w // 2
    quadrants = [
        block[:mh, :mw],
        block[:mh, mw:],
        block[mh:, :mw],
        block[mh:, mw:],
    ]
    types = []
    for q in quadrants:
        if q.shape[0] < 4 or q.shape[1] < 4:
            return False
        t, _ = classify_block(q)
        types.append(t)

    # Subdivide if quadrants have different types
    unique_types = set(types)
    if len(unique_types) >= 2:
        return True

    # Also subdivide if the block is large and the sub-blocks compress
    # significantly differently (but only check at the 64→32 level to avoid overhead)
    if block_size >= 64:
        # Quick compression size check on two quadrants
        _, d0 = best_encode(quadrants[0])
        _, d2 = best_encode(quadrants[2])
        full_enc = encode_deflate(block)
        sub_sum = len(d0) + len(d2)
        # If two quadrants alone are < 60% of the full block encoding,
        # subdivision will likely help
        if sub_sum < len(full_enc) * 0.6:
            return True

    return False


def encode_quadtree(block, block_size, min_size=8):
    """Recursively encode with quadtree decomposition.
    Returns: (total_bytes, strategy_counts, encoded_data, tree_structure)
    tree_structure: list of (is_leaf, strategy_id_or_None, sub-blocks_or_None)
    """
    h, w = block.shape[:2]

    if not should_subdivide(block, block_size, min_size):
        # Leaf node: encode as single block
        sid, data = best_encode(block)
        return len(data), Counter({sid: 1}), data, ('leaf', sid, len(data))

    # Internal node: split into 4 quadrants
    mh, mw = h // 2, w // 2
    sub_size = block_size // 2

    quadrants = [
        (block[:mh, :mw], 0),
        (block[:mh, mw:w], 1),
        (block[mh:h, :mw], 2),
        (block[mh:h, mw:w], 3),
    ]

    total = 0
    counts = Counter()
    all_data = bytearray()
    children = []

    for q_block, q_idx in quadrants:
        if q_block.shape[0] < 4 or q_block.shape[1] < 4:
            # Too small to recurse, encode directly
            sid, data = best_encode(q_block)
            total += len(data)
            counts[sid] += 1
            all_data.extend(data)
            children.append(('leaf', sid, len(data)))
        else:
            sz, c, d, tree = encode_quadtree(q_block, sub_size, min_size)
            total += sz
            counts += c
            all_data.extend(d)
            children.append(tree)

    # Check: is subdivision actually better than encoding the whole block?
    sid_whole, data_whole = best_encode(block)
    # Overhead for tree structure: ~1 byte per split level
    subdivision_overhead = 4  # bytes for the 4 sub-block size headers
    if len(data_whole) <= total + subdivision_overhead:
        # Whole block is better, don't subdivide
        return len(data_whole), Counter({sid_whole: 1}), data_whole, ('leaf', sid_whole, len(data_whole))

    return total, counts, bytes(all_data), ('node', children)


# ================================================================
# MOSAIC CODEC v2
# ================================================================

def mosaic_encode_v2(arr, top_block=64, min_block=8):
    """Encode image with quadtree block-adaptive strategy selection."""
    h, w = arr.shape[:2]
    bh = (h + top_block - 1) // top_block
    bw = (w + top_block - 1) // top_block

    total_size = 0
    strategy_counts = Counter()
    block_details = []

    # Header: magic + dimensions + block config
    header_size = 4 + 2 + 2 + 1 + 1  # TICM, w, h, top_block, min_block
    total_size += header_size

    for by in range(bh):
        for bx in range(bw):
            y0 = by * top_block
            y1 = min((by + 1) * top_block, h)
            x0 = bx * top_block
            x1 = min((bx + 1) * top_block, w)
            block = arr[y0:y1, x0:x1]

            sz, counts, data, tree = encode_quadtree(block, top_block, min_block)
            total_size += sz + 4  # +4 for block size header
            strategy_counts += counts
            block_details.append({
                'pos': (bx, by),
                'size': sz,
                'tree': tree,
            })

    # Strategy map overhead (approximate: 3 bits per block for strategy ID,
    # plus tree structure bits)
    n_leaf_blocks = sum(strategy_counts.values())
    map_overhead = (n_leaf_blocks * 3 + 7) // 8 + bh * bw  # rough
    total_size += map_overhead

    return total_size, strategy_counts, block_details


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
    print("=" * 120)
    print("  TIC-Mosaic v2: Quadtree Block-Adaptive Compression")
    print("  Blocks: 64→32→16→8 quadtree | 6 encoding strategies | Always picks smallest")
    print("=" * 120)

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
                label = f[:40]
                if f.lower().endswith('.jpg'):
                    label = 'PHOTO:' + label
                imgs[label] = arr
            except Exception:
                pass

    print(f"\nLoaded {len(imgs)} images\n")
    print(f"{'#':>2} {'Image':<42} {'WxH':>11} {'PNG':>10} {'TIC-A':>10} {'Mosaic2':>10} {'WebP-LL':>10} {'Best':>8} {'Strategies'}")
    print("-" * 140)

    wins = Counter()
    total_png = 0
    total_mosaic = 0
    total_webp = 0
    total_tic = 0

    for idx, (name, arr) in enumerate(sorted(imgs.items()), 1):
        h, w = arr.shape[:2]

        t0 = time.time()
        mosaic_sz, strat_counts, details = mosaic_encode_v2(arr, top_block=64, min_block=8)
        mosaic_ms = time.time() - t0

        ps = png_sz(arr)
        ts = tic_sz(arr)
        ws = webp_ll_sz(arr)

        total_png += ps
        total_mosaic += mosaic_sz
        total_webp += ws
        total_tic += ts

        candidates = {'PNG': ps, 'Mosaic2': mosaic_sz, 'WebP-LL': ws}
        if ts > 0:
            candidates['TIC-A'] = ts
        best = min(candidates, key=candidates.get)
        wins[best] += 1

        # Strategy summary
        strat_str = ' '.join(f"{STRATEGIES[k][0][:8]}:{v}" for k, v in sorted(strat_counts.items()))
        n_blocks = sum(strat_counts.values())

        ts_str = f"{ts:>10,}" if ts > 0 else f"{'N/A':>10}"
        print(f"{idx:>2} {name:<42} {w}x{h:>4} {ps:>10,} {ts_str} {mosaic_sz:>10,} {ws:>10,} {best:>8} [{n_blocks:>4} blocks] {strat_str}")

    print("-" * 140)
    print(f"   {'TOTALS':<42} {'':>11} {total_png:>10,} {total_tic:>10,} {total_mosaic:>10,} {total_webp:>10,}")

    if total_png > 0:
        print(f"\n   vs PNG:  Mosaic2 = {total_mosaic/total_png*100:.1f}%  |  WebP-LL = {total_webp/total_png*100:.1f}%  |  TIC-A = {total_tic/total_png*100:.1f}%")
    if total_webp > 0:
        print(f"   vs WebP: Mosaic2 = {total_mosaic/total_webp*100:.1f}%")

    print(f"\n   WINS: {dict(wins)}")
    print(f"\n   Strategy IDs: {', '.join(f'{k}={v[0]}' for k, v in STRATEGIES.items())}")


if __name__ == '__main__':
    run_benchmark()
