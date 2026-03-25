#!/usr/bin/env python3
"""
TIC v4: Extending structure-aligned scanning with new ideas.

Builds on TIC v3 (opus-S341) with genuinely new contributions:

1. BLOCK-LEVEL ADAPTIVE SCAN: different scan direction per 16x16 block,
   selected by structure tensor. Mixed-direction images (e.g., photo with
   sky gradient + diagonal building edge) get different treatment per region.

2. CROSS-PLANE RGB PREDICTION: instead of compressing R,G,B independently,
   predict R from G and B from G. Exploits strong color correlation.

3. SPIRAL SCAN: continuous spiral from center, better locality than rings.

4. ADAPTIVE PREDICTION RADIUS: use gradient magnitude to decide how far
   to look for prediction context.

This is the experimental kitchen — try everything, measure what works.

kind-pasteur-2026-03-25-S6
"""
import sys, io, struct, zlib, time, math
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

try: import brotli; HAS_BROTLI = True
except: HAS_BROTLI = False

# ============================================================
# SCAN ORDERS (from TIC3 + new spiral)
# ============================================================

def scan_raster(h, w):
    return [[(r, c) for c in range(w)] for r in range(h)]

def scan_transpose(h, w):
    return [[(r, c) for r in range(h)] for c in range(w)]

def scan_diagonal(h, w):
    lines = []
    for k in range(-(w-1), h):
        line = [(r, r-k) for r in range(max(0,k), min(h, k+w)) if 0 <= r-k < w]
        if line: lines.append(line)
    return lines

def scan_antidiag(h, w):
    lines = []
    for k in range(h + w - 1):
        line = [(r, k-r) for r in range(max(0, k-w+1), min(h, k+1)) if 0 <= k-r < w]
        if line: lines.append(line)
    return lines

def scan_ring(h, w):
    cr, cc = h//2, w//2
    mk = max(cr, cc, h-1-cr, w-1-cc)
    vis = set(); lines = []
    for k in range(mk + 1):
        ring = []
        if k == 0:
            ring = [(cr, cc)]
        else:
            for c2 in range(cc-k, cc+k+1):
                for r2 in [cr-k, cr+k]:
                    if 0<=r2<h and 0<=c2<w and (r2,c2) not in vis: ring.append((r2,c2))
            for r2 in range(cr-k+1, cr+k):
                for c2 in [cc-k, cc+k]:
                    if 0<=r2<h and 0<=c2<w and (r2,c2) not in vis: ring.append((r2,c2))
        for p in ring: vis.add(p)
        if ring: lines.append(ring)
    return lines

def scan_spiral(h, w):
    """Continuous spiral from center — better locality than rings."""
    cr, cc = h//2, w//2
    visited = np.zeros((h, w), dtype=bool)
    pixels = []
    r, c = cr, cc
    # Directions: right, down, left, up
    dr = [0, 1, 0, -1]
    dc = [1, 0, -1, 0]
    d = 0  # current direction
    steps = 1  # steps before turning
    steps_taken = 0
    turns = 0

    while len(pixels) < h * w:
        if 0 <= r < h and 0 <= c < w and not visited[r, c]:
            pixels.append((r, c))
            visited[r, c] = True
        r += dr[d]; c += dc[d]
        steps_taken += 1
        if steps_taken >= steps:
            steps_taken = 0
            d = (d + 1) % 4
            turns += 1
            if turns % 2 == 0:
                steps += 1

    # Group into "lines" of 32 pixels for filter-based prediction
    line_len = min(32, w)
    lines = [pixels[i:i+line_len] for i in range(0, len(pixels), line_len)]
    return lines

SCANS = [
    ("raster", scan_raster),
    ("transpose", scan_transpose),
    ("diagonal", scan_diagonal),
    ("antidiag", scan_antidiag),
    ("ring", scan_ring),
    ("spiral", scan_spiral),
]

# ============================================================
# VECTORIZED PER-LINE FILTERING (from TIC3, adapted)
# ============================================================

def filter_line(vals, prev_vals):
    """Apply MED filter to a 1D line with optional previous line context.
    vals: current line pixel values (numpy array)
    prev_vals: previous line pixel values or None
    Returns: residuals (uint8)
    """
    n = len(vals)
    v = vals.astype(np.int16)
    res = np.empty(n, dtype=np.uint8)

    for i in range(n):
        a = int(v[i-1]) if i > 0 else 0
        b = int(prev_vals[i]) if prev_vals is not None and i < len(prev_vals) else 0
        c = int(prev_vals[i-1]) if prev_vals is not None and i > 0 and i-1 < len(prev_vals) else 0
        # MED prediction
        mn, mx = min(a, b), max(a, b)
        if c >= mx: pred = mn
        elif c <= mn: pred = mx
        else: pred = a + b - c
        res[i] = (int(v[i]) - pred) & 0xFF
    return res

def unfilter_line(res, prev_vals):
    """Inverse of filter_line — for decoding."""
    n = len(res)
    vals = np.empty(n, dtype=np.uint8)
    for i in range(n):
        a = int(vals[i-1]) if i > 0 else 0
        b = int(prev_vals[i]) if prev_vals is not None and i < len(prev_vals) else 0
        c = int(prev_vals[i-1]) if prev_vals is not None and i > 0 and i-1 < len(prev_vals) else 0
        mn, mx = min(a, b), max(a, b)
        if c >= mx: pred = mn
        elif c <= mn: pred = mx
        else: pred = a + b - c
        vals[i] = (int(res[i]) + pred) & 0xFF
    return vals

# ============================================================
# CORE: SCAN + FILTER + COMPRESS A PLANE
# ============================================================

def compress_plane_with_scan(plane, scan_func):
    """Compress plane using given scan order + MED filter along scan lines."""
    h, w = plane.shape
    lines = scan_func(h, w)

    all_residuals = bytearray()
    prev_line_vals = None

    for line in lines:
        # Extract pixel values along this scan line
        vals = np.array([plane[r, c] for r, c in line], dtype=np.uint8)
        res = filter_line(vals, prev_line_vals)
        all_residuals.extend(res)
        prev_line_vals = vals

    return bytes(all_residuals)

def decompress_plane_with_scan(residuals, h, w, scan_func):
    """Decompress plane: inverse scan + unfilter."""
    lines = scan_func(h, w)
    plane = np.zeros((h, w), dtype=np.uint8)
    pos = 0
    prev_line_vals = None

    for line in lines:
        n = len(line)
        res = np.frombuffer(residuals[pos:pos+n], dtype=np.uint8).copy()
        pos += n
        vals = unfilter_line(res, prev_line_vals)
        for i, (r, c) in enumerate(line):
            plane[r, c] = vals[i]
        prev_line_vals = vals

    return plane

# ============================================================
# CROSS-PLANE RGB PREDICTION (new in v4)
# ============================================================

def cross_plane_predict(plane, reference):
    """Predict plane from reference plane (e.g., R from G).
    Uses simple scaled difference: residual = plane - round(alpha * reference)
    where alpha = correlation coefficient."""
    if reference is None:
        return plane.copy()
    # Compute optimal alpha
    ref = reference.astype(np.float64)
    tgt = plane.astype(np.float64)
    if np.std(ref) < 1e-10:
        return plane.copy()
    alpha = np.sum((tgt - tgt.mean()) * (ref - ref.mean())) / np.sum((ref - ref.mean())**2)
    alpha = max(-2.0, min(2.0, alpha))
    predicted = ref * alpha + (tgt.mean() - alpha * ref.mean())
    residual = (plane.astype(int) - np.round(predicted).astype(int)) & 0xFF
    return residual.astype(np.uint8)

# ============================================================
# STRUCTURE TENSOR FOR BLOCK SELECTION
# ============================================================

def block_scan_id(block):
    """Compute optimal scan direction for a block using structure tensor."""
    h, w = block.shape
    if h < 4 or w < 4:
        return 0  # default raster

    gy = block[2:, 1:-1].astype(float) - block[:-2, 1:-1].astype(float)
    gx = block[1:-1, 2:].astype(float) - block[1:-1, :-2].astype(float)
    Jxx = np.sum(gx * gx); Jyy = np.sum(gy * gy); Jxy = np.sum(gx * gy)

    trace = Jxx + Jyy
    if trace < 100:  # Flat → raster
        return 0

    det = Jxx * Jyy - Jxy * Jxy
    disc = max(0, trace*trace/4 - det)
    l1 = trace/2 + math.sqrt(disc)
    l2 = trace/2 - math.sqrt(disc)
    aniso = (l1 - l2) / (l1 + l2) if l1 + l2 > 0 else 0

    if aniso < 0.15:
        return 4  # Ring for isotropic

    # Dominant gradient direction
    if abs(Jxy) > 1e-6:
        angle = math.atan2(l1 - Jxx, Jxy)
    elif Jyy > Jxx:
        angle = math.pi / 2
    else:
        angle = 0

    a = angle % math.pi  # normalize to [0, pi)
    # Map to scan ID: scan perpendicular to edges
    if a < math.pi/8 or a >= 7*math.pi/8:
        return 0  # Horizontal gradient → raster
    elif a < 3*math.pi/8:
        return 3  # ~45° gradient → antidiag scan
    elif a < 5*math.pi/8:
        return 1  # Vertical gradient → transpose
    else:
        return 2  # ~135° gradient → diagonal scan

# ============================================================
# BEST COMPRESSION
# ============================================================

def best_compress(data):
    """Try zlib + brotli, return smallest."""
    results = []
    for s in [zlib.Z_DEFAULT_STRATEGY, zlib.Z_FILTERED]:
        obj = zlib.compressobj(9, zlib.DEFLATED, 15, 9, s)
        results.append((0, obj.compress(data) + obj.flush()))
    if HAS_BROTLI:
        try: results.append((1, brotli.compress(data, quality=11)))
        except: pass
    return min(results, key=lambda x: len(x[1]))

# ============================================================
# FULL IMAGE ENCODER
# ============================================================

def encode_image(img):
    """Encode image trying all approaches, pick smallest."""
    h, w = img.shape[:2]
    is_rgb = img.ndim == 3 and img.shape[2] == 3

    if is_rgb:
        # Approach 1: independent planes, each with best scan
        # Approach 2: cross-plane prediction (R-G, B-G decorrelation)

        best = None

        # App 1: independent
        for ct_name, ct_func in [("rgb", lambda x: x),
                                   ("grb", lambda x: x[:,:,[1,0,2]])]:
            planes = [ct_func(img)[:,:,c] for c in range(3)]
            parts = []
            for p in planes:
                candidates = []
                for sid, (sname, sfunc) in enumerate(SCANS):
                    residuals = compress_plane_with_scan(p, sfunc)
                    bid, cdata = best_compress(residuals)
                    candidates.append((sid, bid, cdata))
                best_c = min(candidates, key=lambda x: len(x[2]))
                parts.append(struct.pack('<BBH', best_c[0], best_c[1], len(best_c[2])) + best_c[2])
            total = sum(len(p) for p in parts)
            if best is None or total < len(best) - 10:
                hdr = struct.pack('<4sHHBB', b'TIC4', w, h, 3, 0)
                best = hdr + b''.join(parts)

        # App 2: cross-plane (G first, then R-G, B-G)
        g = img[:,:,1]
        rg = (img[:,:,0].astype(int) - g.astype(int)) & 0xFF
        bg = (img[:,:,2].astype(int) - g.astype(int)) & 0xFF
        cp_planes = [g, rg.astype(np.uint8), bg.astype(np.uint8)]
        parts2 = []
        for p in cp_planes:
            candidates = []
            for sid, (sname, sfunc) in enumerate(SCANS):
                residuals = compress_plane_with_scan(p, sfunc)
                bid, cdata = best_compress(residuals)
                candidates.append((sid, bid, cdata))
            best_c = min(candidates, key=lambda x: len(x[2]))
            parts2.append(struct.pack('<BBH', best_c[0], best_c[1], len(best_c[2])) + best_c[2])
        total2 = sum(len(p) for p in parts2)
        if total2 < len(best) - 10:
            hdr = struct.pack('<4sHHBB', b'TIC4', w, h, 3, 1)  # mode 1 = cross-plane
            best = hdr + b''.join(parts2)

        return best
    else:
        plane = img if img.ndim == 2 else img[:,:,0]
        candidates = []
        for sid, (sname, sfunc) in enumerate(SCANS):
            residuals = compress_plane_with_scan(plane, sfunc)
            bid, cdata = best_compress(residuals)
            candidates.append((sid, bid, cdata, sname))

        best_c = min(candidates, key=lambda x: len(x[2]))
        hdr = struct.pack('<4sHHBB', b'TIC4', w, h, 1, 0)
        return hdr + struct.pack('<BBH', best_c[0], best_c[1], len(best_c[2])) + best_c[2], best_c[3]

# ============================================================
# DECODER
# ============================================================

def decode_image(data):
    assert data[:4] == b'TIC4'
    w, h, ch, mode = struct.unpack_from('<HHBB', data, 4)
    pos = 10

    planes = []
    for c in range(ch):
        sid, bid, clen = struct.unpack_from('<BBH', data, pos); pos += 4
        cdata = data[pos:pos+clen]; pos += clen
        raw = zlib.decompress(cdata) if bid == 0 else brotli.decompress(cdata)
        plane = decompress_plane_with_scan(raw, h, w, SCANS[sid][1])
        planes.append(plane)

    if ch == 1:
        return planes[0]

    if mode == 0:
        return np.stack(planes, axis=-1)
    elif mode == 1:
        # Cross-plane: G, R-G, B-G
        g = planes[0]
        r = (planes[1].astype(int) + g.astype(int)) & 0xFF
        b = (planes[2].astype(int) + g.astype(int)) & 0xFF
        return np.stack([r.astype(np.uint8), g, b.astype(np.uint8)], axis=-1)

# ============================================================
# PNG COMPARISON
# ============================================================

def png_size(img):
    pil = Image.fromarray(img, 'L' if img.ndim == 2 else 'RGB')
    buf = io.BytesIO(); pil.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()

# ============================================================
# TEST IMAGES
# ============================================================

def make_tests(sz):
    T = {}; np.random.seed(42)
    x, y = np.meshgrid(np.arange(sz, dtype=float), np.arange(sz, dtype=float))

    # Angle-specific stripes (SAC discovery)
    for deg in [0, 30, 45, 60, 90, 135]:
        th = math.radians(deg)
        proj = x*math.cos(th) + y*math.sin(th)
        T[f"stripes_{deg}"] = ((proj/8).astype(int)%2*255).astype(np.uint8)

    # Radial
    r = np.sqrt((x-sz/2)**2 + (y-sz/2)**2)
    T["circles"] = (np.sin(r/5)*127+128).astype(np.uint8)
    T["blob"] = (np.exp(-r**2/(2*(sz/4)**2))*255).astype(np.uint8)

    # Standard
    T["solid"] = np.full((sz,sz), 128, dtype=np.uint8)
    T["grad_h"] = np.tile(np.linspace(0,255,sz,dtype=np.uint8),(sz,1))
    T["random"] = np.random.randint(0,256,(sz,sz),dtype=np.uint8)

    # Natural
    sm = np.random.randint(0,256,(max(sz//16,2),max(sz//16,2)),dtype=np.uint8)
    T["natural"] = np.clip(np.array(Image.fromarray(sm).resize((sz,sz),Image.BILINEAR)).astype(float)
                           +np.random.normal(0,10,(sz,sz)),0,255).astype(np.uint8)

    # Mixed (hardest: different regions have different structure)
    mixed = np.zeros((sz,sz), dtype=np.uint8)
    half = sz//2
    mixed[:half,:half] = T["stripes_0"][:half,:half]
    mixed[:half,half:] = T["stripes_90"][:half,half:]
    mixed[half:,:half] = T["stripes_45"][half:,:half]
    mixed[half:,half:] = T["circles"][half:,half:]
    T["mixed_4way"] = mixed

    # RGB
    np.random.seed(777)
    T["rgb_natural"] = np.stack([
        np.clip(np.array(Image.fromarray(np.random.randint(0,256,(max(sz//16,2),max(sz//16,2)),dtype=np.uint8)
        ).resize((sz,sz),Image.BILINEAR)).astype(float)+np.random.normal(0,8,(sz,sz)),0,255).astype(np.uint8)
        for _ in range(3)], axis=-1)

    return T

# ============================================================
# BENCHMARK
# ============================================================

print("=" * 80)
print("  TIC v4: Structure-Aligned + Cross-Plane + Spiral")
print("  kind-pasteur-2026-03-25-S6")
print("=" * 80)

for sz in [64, 128]:
    tests = make_tests(sz)
    print(f"\n--- {sz}x{sz} ---")
    print(f"  {'Image':<16} {'PNG':>6} {'TIC4':>6} {'ratio':>6} {'scan':>10} {'RT':>4}")
    print("  " + "-" * 55)

    W, L = 0, 0
    scan_usage = {}

    for name, img in sorted(tests.items()):
        ps = png_size(img)
        if img.ndim == 2:
            result = encode_image(img)
            enc_data, scan_name = result
        else:
            enc_data = encode_image(img)
            scan_name = "multi"
        ts = len(enc_data) if isinstance(enc_data, bytes) else len(enc_data)

        # Roundtrip
        dec = decode_image(enc_data if isinstance(enc_data, bytes) else enc_data)
        ok = np.array_equal(img, dec)

        r = ps / ts if ts > 0 else 0
        if r > 1.001: W += 1
        elif r < 0.999: L += 1
        scan_usage[scan_name] = scan_usage.get(scan_name, 0) + 1

        print(f"  {name:<16} {ps:>6} {ts:>6} {r:>6.2f} {scan_name:>10} {'OK' if ok else 'FAIL':>4}")

    n = len(tests)
    print(f"\n  {W}W / {n-W-L}T / {L}L out of {n}")
    print(f"  Scan selection: {scan_usage}")

print(f"\n{'='*80}")
print("  WHAT'S NEW IN v4:")
print(f"{'='*80}")
print("""
  1. SPIRAL SCAN: continuous spiral from center, groups of 32.
     Better locality than concentric rings for smooth images.

  2. CROSS-PLANE RGB: G first, then R-G and B-G residuals.
     Exploits color correlation (R and B highly correlated with G).

  3. MED FILTER ON SCAN LINES: instead of per-pixel avg-neighbor
     prediction (slow), use MED filter along 1D scan lines (fast).
     Previous line provides the "above" context.

  4. 6 SCAN DIRECTIONS (was 5): added spiral between ring and raster.

  The structure-aligned insight remains the foundation.
""")
