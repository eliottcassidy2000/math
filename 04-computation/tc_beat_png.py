#!/usr/bin/env python3
"""
Tournament Codec: BEAT PNG ON EVERYTHING.

Strategy: PNG uses 5 scanline filters + Deflate. We use:
1. More filters (8+ including diagonal, GAP, LOCO-I style)
2. Per-row optimal selection (try all, pick smallest compressed)
3. Color transforms (YCoCg better than RGB)
4. Per-plane adaptive encoding
5. Same Deflate backend (zlib-9) — fair comparison

The ONLY way to beat PNG with same backend is BETTER PREDICTION.

kind-pasteur-2026-03-25-S2
"""
import sys, io, struct, zlib, time, hashlib
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

# ============================================================
# PREDICTORS (more than PNG's 5)
# ============================================================

def pred_none(row, prev, x):
    """PNG filter 0: no prediction."""
    return 0

def pred_sub(row, prev, x):
    """PNG filter 1: predict from left."""
    return row[x-1] if x > 0 else 0

def pred_up(row, prev, x):
    """PNG filter 2: predict from above."""
    return prev[x] if prev is not None else 0

def pred_avg(row, prev, x):
    """PNG filter 3: average of left and above."""
    a = row[x-1] if x > 0 else 0
    b = prev[x] if prev is not None else 0
    return (a + b) // 2

def pred_paeth(row, prev, x):
    """PNG filter 4: Paeth predictor."""
    a = row[x-1] if x > 0 else 0
    b = prev[x] if prev is not None else 0
    c = prev[x-1] if (prev is not None and x > 0) else 0
    p = a + b - c
    pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
    if pa <= pb and pa <= pc:
        return a
    elif pb <= pc:
        return b
    return c

def pred_med(row, prev, x):
    """MED (Median Edge Detector) — better than Paeth for edges.
    Used in JPEG-LS / LOCO-I."""
    a = row[x-1] if x > 0 else 0
    b = prev[x] if prev is not None else 0
    c = prev[x-1] if (prev is not None and x > 0) else 0
    mn, mx = min(a, b), max(a, b)
    if c >= mx:
        return mn
    elif c <= mn:
        return mx
    return a + b - c

def pred_gap(row, prev, x):
    """GAP (Gradient-Adjusted Predictor) — from CALIC.
    Uses gradient to adaptively weight horizontal vs vertical."""
    a = row[x-1] if x > 0 else 0
    b = prev[x] if prev is not None else 0
    c = prev[x-1] if (prev is not None and x > 0) else 0
    d = prev[x+1] if (prev is not None and x + 1 < len(prev)) else b
    # Gradient
    dh = abs(a - c) + abs(b - d)  # horizontal discontinuity
    dv = abs(b - c) + abs(a - (row[x-2] if x >= 2 else 0))  # vertical discontinuity
    if dh > dv + dv:      # strong horizontal edge → predict from above
        return b
    elif dv > dh + dh:    # strong vertical edge → predict from left
        return a
    else:                  # smooth → average
        return (a + b) // 2

def pred_diagonal(row, prev, x):
    """Diagonal predictor: above-left."""
    return prev[x-1] if (prev is not None and x > 0) else 0

def pred_aboveright(row, prev, x):
    """Above-right predictor: useful for diagonal edges going the other way."""
    if prev is not None and x + 1 < len(prev):
        return prev[x+1]
    return prev[x] if prev is not None else 0

def pred_2d_avg(row, prev, x):
    """2D average: average of left, above, above-left, above-right."""
    vals = []
    if x > 0:
        vals.append(row[x-1])
    if prev is not None:
        vals.append(prev[x])
        if x > 0:
            vals.append(prev[x-1])
        if x + 1 < len(prev):
            vals.append(prev[x+1])
    return sum(vals) // len(vals) if vals else 0

def pred_loco(row, prev, x):
    """LOCO-I predictor with context-based correction.
    Uses gradient context to bias MED prediction."""
    a = row[x-1] if x > 0 else 0
    b = prev[x] if prev is not None else 0
    c = prev[x-1] if (prev is not None and x > 0) else 0
    d = prev[x+1] if (prev is not None and x + 1 < len(prev)) else b

    # MED base prediction
    mn, mx = min(a, b), max(a, b)
    if c >= mx:
        base = mn
    elif c <= mn:
        base = mx
    else:
        base = a + b - c

    # Gradient context correction
    # Quantize gradient into 3 levels
    gh = a - c   # horizontal gradient
    gv = b - c   # vertical gradient
    # Slight correction toward gradient direction
    if abs(gh) > 2 * abs(gv) and abs(gh) > 8:
        # Strong horizontal gradient → lean toward vertical (above)
        correction = (b - base) // 4
    elif abs(gv) > 2 * abs(gh) and abs(gv) > 8:
        # Strong vertical gradient → lean toward horizontal (left)
        correction = (a - base) // 4
    else:
        correction = 0

    return np.clip(base + correction, 0, 255)

PREDICTORS = [
    ("none", pred_none),
    ("sub", pred_sub),
    ("up", pred_up),
    ("avg", pred_avg),
    ("paeth", pred_paeth),
    ("med", pred_med),
    ("gap", pred_gap),
    ("diag", pred_diagonal),
    ("aboveR", pred_aboveright),
    ("2davg", pred_2d_avg),
    ("loco", pred_loco),
]

# ============================================================
# ENCODING ENGINE
# ============================================================

def apply_filter(plane, pred_func):
    """Apply prediction filter to entire plane, return residuals (mod 256)."""
    h, w = plane.shape
    out = np.zeros_like(plane)
    for y in range(h):
        row = plane[y].copy()
        prev = plane[y-1] if y > 0 else None
        for x in range(w):
            p = pred_func(row, prev, x)
            out[y, x] = (int(plane[y, x]) - int(p)) & 0xFF
    return out

def apply_filter_row(plane_row, prev_row, pred_func, width):
    """Apply filter to a single row."""
    out = np.zeros(width, dtype=np.uint8)
    row = plane_row.copy()
    for x in range(width):
        p = pred_func(row, prev_row, x)
        out[x] = (int(plane_row[x]) - int(p)) & 0xFF
    return out

def compress_plane_perrow(plane, level=9):
    """Per-row optimal filter selection (like PNG but with more filters).
    Try all filters for each row, pick the one that compresses best."""
    h, w = plane.shape
    # For each row, try each filter, compress with zlib, pick smallest
    filtered_data = bytearray()

    for y in range(h):
        row = plane[y]
        prev = plane[y-1] if y > 0 else None
        best_size = float('inf')
        best_fdata = None
        best_fid = 0

        for fid, (fname, pfunc) in enumerate(PREDICTORS):
            frow = apply_filter_row(row, prev, pfunc, w)
            # Prepend filter ID byte + compressed row
            trial = bytes([fid]) + bytes(frow)
            csize = len(zlib.compress(trial, level))
            if csize < best_size:
                best_size = csize
                best_fdata = trial
                best_fid = fid

        filtered_data.extend(best_fdata)

    return zlib.compress(bytes(filtered_data), level)

def compress_plane_wholeframe(plane, level=9):
    """Try each filter for the whole frame, pick best."""
    best_size = float('inf')
    best_data = None
    best_name = None

    for fid, (fname, pfunc) in enumerate(PREDICTORS):
        filtered = apply_filter(plane, pfunc)
        # Prepend filter ID
        data = bytes([fid]) + filtered.tobytes()
        cdata = zlib.compress(data, level)
        if len(cdata) < best_size:
            best_size = len(cdata)
            best_data = cdata
            best_name = fname

    return best_data, best_name

def rgb_to_ycocg(img):
    """Convert RGB to YCoCg (reversible color transform)."""
    r, g, b = img[:,:,0].astype(int), img[:,:,1].astype(int), img[:,:,2].astype(int)
    Co = r - b
    tmp = b + (Co >> 1)
    Cg = g - tmp
    Y = tmp + (Cg >> 1)
    return np.stack([Y & 0xFF, Co & 0xFF, Cg & 0xFF], axis=-1).astype(np.uint8)

def compress_image(img_array):
    """Compress RGB image using our best codec.
    Returns compressed bytes."""
    h, w = img_array.shape[:2]
    channels = img_array.shape[2] if img_array.ndim == 3 else 1

    if channels == 3:
        # Try both RGB and YCoCg
        planes_rgb = [img_array[:,:,c] for c in range(3)]
        planes_ycocg = [rgb_to_ycocg(img_array)[:,:,c] for c in range(3)]

        # Compress each set
        rgb_sizes = []
        ycocg_sizes = []
        rgb_data = []
        ycocg_data = []

        for planes, sizes, datas in [(planes_rgb, rgb_sizes, rgb_data),
                                      (planes_ycocg, ycocg_sizes, ycocg_data)]:
            for p in planes:
                # Try per-row and whole-frame, pick best
                pr_data = compress_plane_perrow(p)
                wf_data, wf_name = compress_plane_wholeframe(p)
                if len(pr_data) < len(wf_data):
                    sizes.append(len(pr_data))
                    datas.append(pr_data)
                else:
                    sizes.append(len(wf_data))
                    datas.append(wf_data)

        rgb_total = sum(rgb_sizes)
        ycocg_total = sum(ycocg_sizes)

        # Header: 12 bytes (magic + dims + color flag)
        header = struct.pack('<4sHHB', b'TC02', w, h, 1 if ycocg_total < rgb_total else 0)

        if ycocg_total < rgb_total:
            return header + b''.join(struct.pack('<I', len(d)) + d for d in ycocg_data)
        else:
            return header + b''.join(struct.pack('<I', len(d)) + d for d in rgb_data)
    else:
        # Grayscale
        plane = img_array if img_array.ndim == 2 else img_array[:,:,0]
        pr_data = compress_plane_perrow(plane)
        wf_data, wf_name = compress_plane_wholeframe(plane)
        data = pr_data if len(pr_data) < len(wf_data) else wf_data
        header = struct.pack('<4sHHB', b'TC02', w, h, 2)
        return header + struct.pack('<I', len(data)) + data

def png_size(img_array):
    """Get PNG compressed size using optimal settings."""
    img = Image.fromarray(img_array)
    buf = io.BytesIO()
    img.save(buf, format='PNG', optimize=True)
    return buf.tell()

# ============================================================
# TEST IMAGE GENERATORS
# ============================================================

def gen_gradient_h(w, h):
    """Horizontal gradient."""
    return np.tile(np.linspace(0, 255, w, dtype=np.uint8), (h, 1))

def gen_gradient_v(w, h):
    """Vertical gradient."""
    return np.tile(np.linspace(0, 255, h, dtype=np.uint8).reshape(-1, 1), (1, w))

def gen_gradient_diag(w, h):
    """Diagonal gradient."""
    x = np.arange(w)
    y = np.arange(h).reshape(-1, 1)
    return ((x + y) * 255 // (w + h - 2)).astype(np.uint8)

def gen_solid(w, h, val=128):
    """Solid color."""
    return np.full((h, w), val, dtype=np.uint8)

def gen_checkerboard(w, h, block=8):
    """Checkerboard pattern."""
    x = np.arange(w) // block
    y = np.arange(h).reshape(-1, 1) // block
    return ((x + y) % 2 * 255).astype(np.uint8)

def gen_random(w, h):
    """Random noise."""
    return np.random.randint(0, 256, (h, w), dtype=np.uint8)

def gen_smooth_noise(w, h, scale=32):
    """Smooth noise (photo-like)."""
    small = np.random.randint(0, 256, (h // scale + 1, w // scale + 1))
    # Bilinear upscale
    from PIL import Image as PILImage
    img = PILImage.fromarray(small.astype(np.uint8), 'L')
    img = img.resize((w, h), PILImage.BILINEAR)
    return np.array(img)

def gen_edges(w, h):
    """Sharp edges: blocks of different gray levels."""
    img = np.zeros((h, w), dtype=np.uint8)
    np.random.seed(42)
    for _ in range(20):
        x0, y0 = np.random.randint(0, w-32), np.random.randint(0, h-32)
        bw, bh = np.random.randint(16, 64), np.random.randint(16, 64)
        val = np.random.randint(0, 256)
        img[y0:min(y0+bh, h), x0:min(x0+bw, w)] = val
    return img

def gen_text_like(w, h):
    """Text-like: mostly white with thin dark lines."""
    img = np.full((h, w), 255, dtype=np.uint8)
    np.random.seed(123)
    for _ in range(h // 8):
        y = np.random.randint(0, h)
        x_start = np.random.randint(0, w // 2)
        x_end = x_start + np.random.randint(10, w // 2)
        img[y, x_start:min(x_end, w)] = np.random.randint(0, 50)
    return img

def gen_natural_like(w, h):
    """Natural photo simulation: smooth regions + noise + edges."""
    # Base: smooth noise
    base = gen_smooth_noise(w, h, scale=16).astype(float)
    # Add Gaussian noise
    noise = np.random.normal(0, 8, (h, w))
    # Add some edges
    edges = np.zeros((h, w))
    for _ in range(5):
        y = np.random.randint(h // 4, 3 * h // 4)
        edges[y:y+2, :] = np.random.randint(-50, 50)
    return np.clip(base + noise + edges, 0, 255).astype(np.uint8)

def gen_dithered(w, h):
    """Floyd-Steinberg dithered gradient."""
    gradient = np.linspace(0, 255, w).reshape(1, -1).repeat(h, axis=0).astype(float)
    out = np.zeros((h, w), dtype=np.uint8)
    for y in range(h):
        for x in range(w):
            old = gradient[y, x]
            new = 255 if old > 128 else 0
            out[y, x] = new
            err = old - new
            if x + 1 < w:
                gradient[y, x+1] += err * 7/16
            if y + 1 < h:
                if x > 0:
                    gradient[y+1, x-1] += err * 3/16
                gradient[y+1, x] += err * 5/16
                if x + 1 < w:
                    gradient[y+1, x+1] += err * 1/16
    return out

def gen_circles(w, h):
    """Concentric circles."""
    cx, cy = w // 2, h // 2
    x = np.arange(w) - cx
    y = (np.arange(h) - cy).reshape(-1, 1)
    r = np.sqrt(x*x + y*y)
    return (np.sin(r / 5) * 127 + 128).astype(np.uint8)

def gen_stripes_h(w, h, period=4):
    """Horizontal stripes."""
    return (np.arange(h).reshape(-1, 1) % period < period // 2).astype(np.uint8) * 255

def gen_sparse(w, h, density=0.01):
    """Sparse image: mostly zero with random dots."""
    img = np.zeros((h, w), dtype=np.uint8)
    n_dots = int(w * h * density)
    xs = np.random.randint(0, w, n_dots)
    ys = np.random.randint(0, h, n_dots)
    vals = np.random.randint(1, 256, n_dots)
    img[ys, xs] = vals
    return img

def gen_rgb_gradient(w, h):
    """RGB color gradient."""
    img = np.zeros((h, w, 3), dtype=np.uint8)
    img[:,:,0] = np.tile(np.linspace(0, 255, w, dtype=np.uint8), (h, 1))
    img[:,:,1] = np.tile(np.linspace(0, 255, h, dtype=np.uint8).reshape(-1, 1), (1, w))
    img[:,:,2] = 128
    return img

def gen_rgb_natural(w, h):
    """RGB natural photo simulation."""
    np.random.seed(777)
    channels = []
    for c in range(3):
        channels.append(gen_natural_like(w, h))
    return np.stack(channels, axis=-1)

def gen_rgb_edges(w, h):
    """RGB with sharp color edges."""
    img = np.zeros((h, w, 3), dtype=np.uint8)
    np.random.seed(99)
    for _ in range(15):
        x0, y0 = np.random.randint(0, w-32), np.random.randint(0, h-32)
        bw, bh = np.random.randint(16, 80), np.random.randint(16, 80)
        color = np.random.randint(0, 256, 3)
        img[y0:min(y0+bh, h), x0:min(x0+bw, w)] = color
    return img

# ============================================================
# MAIN BENCHMARK
# ============================================================

print("=" * 72)
print("  TOURNAMENT CODEC vs PNG: COMPREHENSIVE BENCHMARK")
print("  kind-pasteur-2026-03-25-S2")
print("=" * 72)

SIZE = 128  # Use 128x128 for speed; scale results

# Define ALL test images
TESTS = [
    # Grayscale
    ("gray_solid", lambda: gen_solid(SIZE, SIZE)),
    ("gray_grad_h", lambda: gen_gradient_h(SIZE, SIZE)),
    ("gray_grad_v", lambda: gen_gradient_v(SIZE, SIZE)),
    ("gray_grad_diag", lambda: gen_gradient_diag(SIZE, SIZE)),
    ("gray_checker_8", lambda: gen_checkerboard(SIZE, SIZE, 8)),
    ("gray_checker_2", lambda: gen_checkerboard(SIZE, SIZE, 2)),
    ("gray_circles", lambda: gen_circles(SIZE, SIZE)),
    ("gray_stripes_4", lambda: gen_stripes_h(SIZE, SIZE, 4)),
    ("gray_stripes_16", lambda: gen_stripes_h(SIZE, SIZE, 16)),
    ("gray_edges", lambda: gen_edges(SIZE, SIZE)),
    ("gray_text", lambda: gen_text_like(SIZE, SIZE)),
    ("gray_natural", lambda: gen_natural_like(SIZE, SIZE)),
    ("gray_dithered", lambda: gen_dithered(SIZE, SIZE)),
    ("gray_smooth_noise", lambda: gen_smooth_noise(SIZE, SIZE)),
    ("gray_random", lambda: gen_random(SIZE, SIZE)),
    ("gray_sparse_1pct", lambda: gen_sparse(SIZE, SIZE, 0.01)),
    ("gray_sparse_10pct", lambda: gen_sparse(SIZE, SIZE, 0.10)),
    # RGB
    ("rgb_gradient", lambda: gen_rgb_gradient(SIZE, SIZE)),
    ("rgb_natural", lambda: gen_rgb_natural(SIZE, SIZE)),
    ("rgb_edges", lambda: gen_rgb_edges(SIZE, SIZE)),
]

results = []
wins, ties, losses = 0, 0, 0

print(f"\n{'Test Image':<22} {'Raw':>7} {'PNG':>7} {'TC':>7} {'Ratio':>7} {'Result':>8}")
print("-" * 72)

for name, gen_func in TESTS:
    np.random.seed(hash(name) % (2**31))  # Reproducible
    img = gen_func()

    raw_size = img.nbytes

    # PNG size
    psize = png_size(img)

    # Our codec
    t0 = time.time()
    tc_data = compress_image(img)
    tc_time = time.time() - t0
    tsize = len(tc_data)

    ratio = psize / tsize if tsize > 0 else 0
    if ratio > 1.005:
        result = "WIN"
        wins += 1
    elif ratio < 0.995:
        result = "LOSS"
        losses += 1
    else:
        result = "TIE"
        ties += 1

    print(f"{name:<22} {raw_size:>7} {psize:>7} {tsize:>7} {ratio:>7.3f} {result:>8}")
    results.append((name, raw_size, psize, tsize, ratio, result))

print("-" * 72)
total_png = sum(r[2] for r in results)
total_tc = sum(r[3] for r in results)
agg_ratio = total_png / total_tc if total_tc > 0 else 0
print(f"{'AGGREGATE':<22} {'':>7} {total_png:>7} {total_tc:>7} {agg_ratio:>7.3f}")
print(f"\nWins: {wins}  Ties: {ties}  Losses: {losses}  Total: {len(results)}")
print(f"Win rate: {wins/len(results)*100:.0f}%")

# Detail the losses
if losses > 0:
    print(f"\n--- LOSSES (we need to fix these) ---")
    for name, raw, psize, tsize, ratio, result in results:
        if result == "LOSS":
            print(f"  {name}: PNG={psize}, TC={tsize}, ratio={ratio:.3f} ({(1-ratio)*100:+.1f}%)")

# Detail the biggest wins
print(f"\n--- BIGGEST WINS ---")
for name, raw, psize, tsize, ratio, result in sorted(results, key=lambda r: -r[4])[:5]:
    print(f"  {name}: PNG={psize}, TC={tsize}, ratio={ratio:.3f} ({(ratio-1)*100:+.1f}%)")
