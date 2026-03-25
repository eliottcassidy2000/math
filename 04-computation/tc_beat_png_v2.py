#!/usr/bin/env python3
"""
Tournament Codec v2: BEAT PNG ON EVERYTHING.

Fixes from v1:
1. Integer overflow in predictors (use int, clamp to 0-255)
2. Cleaner byte stream: filter map separate from data
3. Try BOTH per-row and whole-frame, pick best per plane
4. Multiple color transforms for RGB
5. Smarter selection: use residual entropy to pick, not just compressed size

kind-pasteur-2026-03-25-S2
"""
import sys, io, struct, zlib, time
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

# ============================================================
# PREDICTORS — use Python int to avoid overflow
# ============================================================

def pix(arr, x, default=0):
    """Safe pixel access."""
    if arr is not None and 0 <= x < len(arr):
        return int(arr[x])
    return default

def pred_none(row, prev, x, w):
    return 0

def pred_sub(row, prev, x, w):
    return pix(row, x-1)

def pred_up(row, prev, x, w):
    return pix(prev, x)

def pred_avg(row, prev, x, w):
    a = pix(row, x-1)
    b = pix(prev, x)
    return (a + b) >> 1

def pred_paeth(row, prev, x, w):
    a = pix(row, x-1)
    b = pix(prev, x)
    c = pix(prev, x-1)
    p = a + b - c
    pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
    if pa <= pb and pa <= pc: return a
    if pb <= pc: return b
    return c

def pred_med(row, prev, x, w):
    """LOCO-I MED predictor."""
    a = pix(row, x-1)
    b = pix(prev, x)
    c = pix(prev, x-1)
    if c >= max(a, b): return min(a, b)
    if c <= min(a, b): return max(a, b)
    return a + b - c

def pred_gap(row, prev, x, w):
    """Gradient-Adjusted Predictor."""
    a = pix(row, x-1)
    b = pix(prev, x)
    c = pix(prev, x-1)
    d = pix(prev, x+1, b)
    e = pix(row, x-2)
    dh = abs(a - c) + abs(b - d)
    dv = abs(b - c) + abs(a - e)
    if dh > 2 * dv: return b        # horizontal edge → vertical predictor
    if dv > 2 * dh: return a        # vertical edge → horizontal predictor
    return (a + b) >> 1

def pred_diag(row, prev, x, w):
    return pix(prev, x-1)

def pred_avg4(row, prev, x, w):
    """Average of 4 neighbors."""
    vals = [pix(row, x-1)]
    if prev is not None:
        vals.append(pix(prev, x))
        vals.append(pix(prev, x-1))
        if x + 1 < w:
            vals.append(pix(prev, x+1))
    return sum(vals) // len(vals)

def pred_loco_ctx(row, prev, x, w):
    """LOCO-I with gradient context correction."""
    a = pix(row, x-1)
    b = pix(prev, x)
    c = pix(prev, x-1)
    # MED
    if c >= max(a, b): base = min(a, b)
    elif c <= min(a, b): base = max(a, b)
    else: base = a + b - c
    # Gradient context
    gh = a - c
    gv = b - c
    if abs(gh) > 2 * abs(gv) + 10:
        base = (base + b) >> 1
    elif abs(gv) > 2 * abs(gh) + 10:
        base = (base + a) >> 1
    return max(0, min(255, base))

def pred_weighted(row, prev, x, w):
    """Weighted average: give more weight to closer neighbors."""
    a = pix(row, x-1)
    b = pix(prev, x)
    c = pix(prev, x-1)
    # Weights: a=3, b=3, c=1 (closer gets more)
    if prev is not None and x > 0:
        return (3*a + 3*b - c) // 5
    elif prev is not None:
        return b
    return a

FILTERS = [
    pred_none, pred_sub, pred_up, pred_avg, pred_paeth,
    pred_med, pred_gap, pred_diag, pred_avg4, pred_loco_ctx,
    pred_weighted,
]
NFILTERS = len(FILTERS)

# ============================================================
# FAST FILTERING WITH NUMPY (vectorized where possible)
# ============================================================

def filter_plane_perrow(plane):
    """Per-row optimal filter selection. Returns (filter_ids, filtered_data)."""
    h, w = plane.shape
    filter_ids = np.zeros(h, dtype=np.uint8)
    filtered_rows = []

    for y in range(h):
        row = plane[y]
        prev = plane[y-1] if y > 0 else None
        best_fid = 0
        best_abs_sum = float('inf')

        # Try each filter, pick lowest absolute sum of residuals
        for fid, pfunc in enumerate(FILTERS):
            residuals = np.zeros(w, dtype=np.uint8)
            for x in range(w):
                p = pfunc(row, prev, x, w)
                residuals[x] = (int(row[x]) - p) & 0xFF
            # Heuristic: minimum absolute sum (maps large residuals to small via mod 256)
            # Convert to signed for better heuristic
            signed = residuals.astype(np.int16)
            signed[signed > 128] -= 256
            abs_sum = np.sum(np.abs(signed))
            if abs_sum < best_abs_sum:
                best_abs_sum = abs_sum
                best_fid = fid
                best_residuals = residuals

        filter_ids[y] = best_fid
        filtered_rows.append(best_residuals)

    filtered_data = np.vstack(filtered_rows)
    return filter_ids, filtered_data

def filter_plane_whole(plane, fid):
    """Apply single filter to entire plane."""
    h, w = plane.shape
    pfunc = FILTERS[fid]
    out = np.zeros_like(plane)
    for y in range(h):
        row = plane[y]
        prev = plane[y-1] if y > 0 else None
        for x in range(w):
            p = pfunc(row, prev, x, w)
            out[y, x] = (int(row[x]) - p) & 0xFF
    return out

# ============================================================
# COMPRESSION STRATEGIES
# ============================================================

def compress_strategy_perrow(plane, level=9):
    """Strategy 1: Per-row filter selection."""
    fids, fdata = filter_plane_perrow(plane)
    # Pack: filter_ids (h bytes) + filtered data (h*w bytes)
    payload = bytes(fids) + fdata.tobytes()
    return struct.pack('<BH', 1, len(fids)) + zlib.compress(payload, level)

def compress_strategy_wholeframe(plane, level=9):
    """Strategy 2: Try each filter for whole frame, pick best."""
    best_size = float('inf')
    best_data = None
    for fid in range(NFILTERS):
        fdata = filter_plane_whole(plane, fid)
        payload = struct.pack('<B', fid) + fdata.tobytes()
        cdata = zlib.compress(payload, level)
        if len(cdata) < best_size:
            best_size = len(cdata)
            best_data = struct.pack('<B', 2) + cdata
    return best_data

def compress_strategy_interleaved(plane, level=9):
    """Strategy 3: PNG-style interleaved (filter byte per row)."""
    h, w = plane.shape
    fids, fdata = filter_plane_perrow(plane)
    # Interleave: for each row, filter_id byte + filtered pixels
    buf = bytearray()
    for y in range(h):
        buf.append(fids[y])
        buf.extend(fdata[y])
    return struct.pack('<B', 3) + zlib.compress(bytes(buf), level)

def compress_strategy_delta_row(plane, level=9):
    """Strategy 4: Row-difference (predict each row from previous)."""
    h, w = plane.shape
    out = np.zeros_like(plane)
    out[0] = plane[0]
    for y in range(1, h):
        out[y] = (plane[y].astype(int) - plane[y-1].astype(int)) & 0xFF
    return struct.pack('<B', 4) + zlib.compress(out.tobytes(), level)

def compress_strategy_raw(plane, level=9):
    """Strategy 5: Raw (no filter, just zlib)."""
    return struct.pack('<B', 5) + zlib.compress(plane.tobytes(), level)

STRATEGIES = [
    compress_strategy_perrow,
    compress_strategy_wholeframe,
    compress_strategy_interleaved,
    compress_strategy_delta_row,
    compress_strategy_raw,
]

def compress_plane_best(plane, level=9):
    """Try all strategies, pick smallest."""
    best = None
    for sfunc in STRATEGIES:
        data = sfunc(plane, level)
        if best is None or len(data) < len(best):
            best = data
    return best

# ============================================================
# COLOR TRANSFORMS
# ============================================================

def rgb_to_ycocg(img):
    r, g, b = img[:,:,0].astype(int), img[:,:,1].astype(int), img[:,:,2].astype(int)
    Co = r - b
    tmp = b + (Co >> 1)
    Cg = g - tmp
    Y = tmp + (Cg >> 1)
    return np.stack([Y & 0xFF, Co & 0xFF, Cg & 0xFF], axis=-1).astype(np.uint8)

def rgb_to_rct(img):
    """Reversible Color Transform (JPEG 2000 style)."""
    r, g, b = img[:,:,0].astype(int), img[:,:,1].astype(int), img[:,:,2].astype(int)
    Y = (r + 2*g + b) >> 2
    Cb = b - g
    Cr = r - g
    return np.stack([Y & 0xFF, Cb & 0xFF, Cr & 0xFF], axis=-1).astype(np.uint8)

def rgb_to_grd(img):
    """Green-Red-Delta transform."""
    r, g, b = img[:,:,0].astype(int), img[:,:,1].astype(int), img[:,:,2].astype(int)
    return np.stack([g & 0xFF, (r - g) & 0xFF, (b - g) & 0xFF], axis=-1).astype(np.uint8)

COLOR_TRANSFORMS = [
    ("RGB", lambda img: img),
    ("YCoCg", rgb_to_ycocg),
    ("RCT", rgb_to_rct),
    ("GRD", rgb_to_grd),
]

# ============================================================
# MAIN COMPRESS FUNCTION
# ============================================================

def compress_image(img_array):
    """Compress image using best combination of color transform + per-plane strategy."""
    h, w = img_array.shape[:2]
    is_rgb = img_array.ndim == 3 and img_array.shape[2] == 3

    if is_rgb:
        best_total = float('inf')
        best_result = None
        for ct_id, (ct_name, ct_func) in enumerate(COLOR_TRANSFORMS):
            transformed = ct_func(img_array)
            planes_data = []
            total = 0
            for c in range(3):
                pdata = compress_plane_best(transformed[:,:,c])
                planes_data.append(pdata)
                total += len(pdata)
            if total < best_total:
                best_total = total
                best_result = (ct_id, planes_data)

        ct_id, planes_data = best_result
        header = struct.pack('<4sHHBB', b'TC03', w, h, 3, ct_id)
        body = b''.join(struct.pack('<I', len(d)) + d for d in planes_data)
        return header + body
    else:
        plane = img_array if img_array.ndim == 2 else img_array[:,:,0]
        pdata = compress_plane_best(plane)
        header = struct.pack('<4sHHBB', b'TC03', w, h, 1, 0)
        return header + struct.pack('<I', len(pdata)) + pdata

def png_size(img_array):
    """Get PNG compressed size."""
    if img_array.ndim == 2:
        img = Image.fromarray(img_array, 'L')
    else:
        img = Image.fromarray(img_array, 'RGB')
    buf = io.BytesIO()
    img.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()

# ============================================================
# TEST IMAGES (comprehensive)
# ============================================================

def gen_solid(w, h): return np.full((h, w), 128, dtype=np.uint8)
def gen_grad_h(w, h): return np.tile(np.linspace(0, 255, w, dtype=np.uint8), (h, 1))
def gen_grad_v(w, h): return np.tile(np.linspace(0, 255, h, dtype=np.uint8).reshape(-1,1), (1, w))
def gen_grad_d(w, h):
    x, y = np.meshgrid(np.arange(w), np.arange(h))
    return ((x + y) * 255 // (w + h - 2)).astype(np.uint8)
def gen_checker(w, h, b=8):
    x, y = np.meshgrid(np.arange(w)//b, np.arange(h)//b)
    return ((x + y) % 2 * 255).astype(np.uint8)
def gen_circles(w, h):
    cx, cy = w//2, h//2
    x, y = np.meshgrid(np.arange(w)-cx, np.arange(h)-cy)
    r = np.sqrt(x*x + y*y)
    return (np.sin(r/5)*127+128).astype(np.uint8)
def gen_random(w, h): return np.random.randint(0, 256, (h, w), dtype=np.uint8)
def gen_sparse(w, h, d=0.01):
    img = np.zeros((h,w), dtype=np.uint8)
    n = int(w*h*d)
    img[np.random.randint(0,h,n), np.random.randint(0,w,n)] = np.random.randint(1,256,n).astype(np.uint8)
    return img

def gen_smooth(w, h, sc=16):
    """Smooth noise — simulates natural photo."""
    small = np.random.randint(0, 256, (max(h//sc,2), max(w//sc,2)), dtype=np.uint8)
    img = Image.fromarray(small).resize((w, h), Image.BILINEAR)
    return np.array(img)

def gen_natural(w, h):
    """Natural photo simulation with noise."""
    base = gen_smooth(w, h, 16).astype(float)
    noise = np.random.normal(0, 10, (h, w))
    return np.clip(base + noise, 0, 255).astype(np.uint8)

def gen_edges(w, h):
    """Sharp blocks."""
    img = np.zeros((h,w), dtype=np.uint8)
    for _ in range(20):
        x0, y0 = np.random.randint(0, max(w-16,1)), np.random.randint(0, max(h-16,1))
        bw, bh = np.random.randint(8, 48), np.random.randint(8, 48)
        img[y0:min(y0+bh,h), x0:min(x0+bw,w)] = np.random.randint(0, 256)
    return img

def gen_text(w, h):
    """Text-like."""
    img = np.full((h,w), 245, dtype=np.uint8)
    for _ in range(h//6):
        y = np.random.randint(0, h)
        xs = np.random.randint(0, w//2)
        xe = xs + np.random.randint(10, w//2)
        img[y, xs:min(xe,w)] = np.random.randint(0, 30)
    return img

def gen_dithered(w, h):
    """Floyd-Steinberg dithered gradient."""
    g = np.linspace(0, 255, w).reshape(1,-1).repeat(h, axis=0).astype(float)
    out = np.zeros((h,w), dtype=np.uint8)
    for y in range(h):
        for x in range(w):
            old = g[y,x]; new = 255 if old > 128 else 0; out[y,x] = new
            err = old - new
            if x+1<w: g[y,x+1] += err*7/16
            if y+1<h:
                if x>0: g[y+1,x-1] += err*3/16
                g[y+1,x] += err*5/16
                if x+1<w: g[y+1,x+1] += err/16
    return out

def gen_stripes(w, h, p=4):
    return (np.arange(h).reshape(-1,1) % p < p//2).astype(np.uint8) * 255

def gen_rgb_grad(w, h):
    img = np.zeros((h,w,3), dtype=np.uint8)
    img[:,:,0] = gen_grad_h(w,h)
    img[:,:,1] = gen_grad_v(w,h)
    img[:,:,2] = 128
    return img

def gen_rgb_natural(w, h):
    np.random.seed(777)
    return np.stack([gen_natural(w,h) for _ in range(3)], axis=-1)

def gen_rgb_edges(w, h):
    img = np.zeros((h,w,3), dtype=np.uint8)
    np.random.seed(99)
    for _ in range(15):
        x0, y0 = np.random.randint(0, max(w-16,1)), np.random.randint(0, max(h-16,1))
        bw, bh = np.random.randint(16, 64), np.random.randint(16, 64)
        img[y0:min(y0+bh,h), x0:min(x0+bw,w)] = np.random.randint(0, 256, 3)
    return img

def gen_rgb_smooth(w, h):
    """Smooth RGB — hardest natural test."""
    np.random.seed(888)
    return np.stack([gen_smooth(w, h, 8) for _ in range(3)], axis=-1)

# ============================================================
# BENCHMARK
# ============================================================

print("=" * 76)
print("  TOURNAMENT CODEC v2 vs PNG: COMPREHENSIVE BENCHMARK")
print("  kind-pasteur-2026-03-25-S2")
print("=" * 76)

SZ = 128

TESTS = [
    ("gray_solid",       lambda: gen_solid(SZ, SZ)),
    ("gray_grad_h",      lambda: gen_grad_h(SZ, SZ)),
    ("gray_grad_v",      lambda: gen_grad_v(SZ, SZ)),
    ("gray_grad_diag",   lambda: gen_grad_d(SZ, SZ)),
    ("gray_checker8",    lambda: gen_checker(SZ, SZ, 8)),
    ("gray_checker2",    lambda: gen_checker(SZ, SZ, 2)),
    ("gray_circles",     lambda: gen_circles(SZ, SZ)),
    ("gray_stripes4",    lambda: gen_stripes(SZ, SZ, 4)),
    ("gray_edges",       lambda: gen_edges(SZ, SZ)),
    ("gray_text",        lambda: gen_text(SZ, SZ)),
    ("gray_natural",     lambda: gen_natural(SZ, SZ)),
    ("gray_dithered",    lambda: gen_dithered(SZ, SZ)),
    ("gray_smooth",      lambda: gen_smooth(SZ, SZ)),
    ("gray_random",      lambda: gen_random(SZ, SZ)),
    ("gray_sparse1",     lambda: gen_sparse(SZ, SZ, 0.01)),
    ("gray_sparse10",    lambda: gen_sparse(SZ, SZ, 0.10)),
    ("rgb_gradient",     lambda: gen_rgb_grad(SZ, SZ)),
    ("rgb_natural",      lambda: gen_rgb_natural(SZ, SZ)),
    ("rgb_edges",        lambda: gen_rgb_edges(SZ, SZ)),
    ("rgb_smooth",       lambda: gen_rgb_smooth(SZ, SZ)),
]

results = []
W, T, L = 0, 0, 0

print(f"\n{'Test':<20} {'Raw':>7} {'PNG':>7} {'TC':>7} {'Ratio':>7} {'Better':>8}")
print("-" * 76)

for name, gf in TESTS:
    np.random.seed(abs(hash(name)) % (2**31))
    img = gf()
    raw = img.nbytes
    ps = png_size(img)
    tc = len(compress_image(img))
    r = ps / tc if tc > 0 else 0
    if r > 1.005: res = "WIN"; W += 1
    elif r < 0.995: res = "LOSS"; L += 1
    else: res = "TIE"; T += 1
    print(f"{name:<20} {raw:>7} {ps:>7} {tc:>7} {r:>7.3f} {res:>8}")
    results.append((name, raw, ps, tc, r, res))

print("-" * 76)
tp = sum(r[2] for r in results)
tt = sum(r[3] for r in results)
print(f"{'AGGREGATE':<20} {'':>7} {tp:>7} {tt:>7} {tp/tt:>7.3f}")
print(f"\n  Wins: {W}  Ties: {T}  Losses: {L}  Total: {len(results)}")
print(f"  Win rate: {W}/{len(results)} = {W/len(results)*100:.0f}%")

if L > 0:
    print(f"\n--- LOSSES ---")
    for n, _, p, t, r, res in results:
        if res == "LOSS":
            print(f"  {n}: PNG={p}, TC={t}, overhead={t-p} bytes ({(1-r)*100:.1f}%)")

print(f"\n--- TOP 5 WINS ---")
for n, _, p, t, r, res in sorted(results, key=lambda x: -x[4])[:5]:
    print(f"  {n}: TC is {r:.2f}x smaller than PNG (saves {p-t} bytes)")
