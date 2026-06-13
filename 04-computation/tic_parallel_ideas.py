#!/usr/bin/env python3
"""
PARALLEL IDEAS: Work around the staircase paradox by doing two things at once.

The paradox: scanning in ONE direction misses structure in OTHER directions.
The fix: instead of choosing ONE scan, do MULTIPLE things simultaneously.

IDEAS TO TRY:

1. DUAL-PASS PREDICTION: predict each pixel from BOTH row-context AND
   column-context, average the two predictions. One pass captures horizontal
   structure, the other captures vertical. Together they catch both.
   Cost: 2x prediction time. Gain: catches 2 directions at once.

2. 2D CONTEXT HASHING: instead of line-based prediction, hash the 2D
   neighborhood pattern into a context ID. Use a lookup table per context.
   Each context captures the LOCAL direction implicitly.
   This is what CALIC does. Let's try a simplified version.

3. INTERLEAVED SCAN: alternate rows between horizontal and vertical scan.
   Even rows: left-to-right. Odd rows: top-to-bottom.
   Each pixel has context from BOTH directions within 1 row distance.

4. GRADIENT-STEERED PREDICTION: for each pixel, compute the gradient
   from known neighbors, then predict along the gradient direction.
   This is a per-pixel adaptive direction, not per-block.

5. CHECKERBOARD DECOMPOSITION: split image into even and odd pixels
   (like a checkerboard). Encode even pixels first (predicting from
   each other), then odd pixels (predicting from ALL even neighbors).
   Each odd pixel has 4 known neighbors in cardinal directions = perfect
   2D context. This kills the staircase paradox because context is symmetric.

6. MULTI-RESOLUTION: encode a 2x-downsampled version first (cheap, small).
   Then encode the full-resolution residual predicting from the downsampled
   version. The downsampled version captures global structure; the residual
   captures local detail. Progressive AND kills the paradox (global structure
   is direction-independent when downsampled).

7. QUINCUNX SCAN: encode every other pixel in a diamond pattern first,
   then fill in the gaps. Similar to checkerboard but with diamond geometry.
   Each gap pixel has 4 neighbors (N,S,E,W) all known. Perfect 2D context.

Let's try ALL of them and measure what works.

kind-pasteur-2026-03-25-S7
"""
import sys, io, struct, zlib, math, time
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

try: import brotli; HAS_BROTLI = True
except: HAS_BROTLI = False

def compress(data):
    results = [zlib.compress(data, 9)]
    if HAS_BROTLI:
        try: results.append(brotli.compress(data, quality=11))
        except: pass
    return min(results, key=len)

def png_size(img):
    pil = Image.fromarray(img, 'L' if img.ndim == 2 else 'RGB')
    buf = io.BytesIO(); pil.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()

def med(a, b, c):
    mn, mx = min(a, b), max(a, b)
    if c >= mx: return mn
    if c <= mn: return mx
    return a + b - c

# ============================================================
# IDEA 1: DUAL-PASS PREDICTION (row + column averaged)
# ============================================================

def encode_dualpass(plane):
    """Predict from BOTH left and above, take average. Scans raster."""
    h, w = plane.shape
    # Pass 1: predict from left (sub filter)
    pred_h = np.zeros_like(plane, dtype=np.int16)
    pred_h[:, 1:] = plane[:, :-1]
    # Pass 2: predict from above (up filter)
    pred_v = np.zeros_like(plane, dtype=np.int16)
    pred_v[1:, :] = plane[:-1, :]
    # Average prediction
    pred = ((pred_h + pred_v) >> 1).astype(np.int16)
    residual = ((plane.astype(np.int16) - pred) & 0xFF).astype(np.uint8)
    return residual.tobytes()

# ============================================================
# IDEA 2: 2D CONTEXT HASHING (simplified CALIC)
# ============================================================

def encode_context_hash(plane):
    """Classify each pixel's neighborhood into gradient contexts, adjust prediction."""
    h, w = plane.shape
    residuals = bytearray()
    # Context correction table: running average error per context
    ctx_sum = np.zeros(16, dtype=np.float64)
    ctx_count = np.ones(16, dtype=np.float64)  # start at 1 to avoid div/0

    for r in range(h):
        for c in range(w):
            a = int(plane[r, c-1]) if c > 0 else 0
            b = int(plane[r-1, c]) if r > 0 else 0
            c2 = int(plane[r-1, c-1]) if r > 0 and c > 0 else 0
            d = int(plane[r-1, c+1]) if r > 0 and c+1 < w else b

            # MED base prediction
            base = med(a, b, c2)

            # Gradient context: quantize (dh, dv) into 4x4 = 16 bins
            dh = abs(a - c2)  # horizontal gradient
            dv = abs(b - c2)  # vertical gradient
            qh = min(3, dh >> 5)  # quantize to 0-3
            qv = min(3, dv >> 5)
            ctx = qh * 4 + qv

            # Correction from running average
            correction = int(ctx_sum[ctx] / ctx_count[ctx])
            pred = max(0, min(255, base + correction))

            actual = int(plane[r, c])
            error = actual - pred
            residuals.append(error & 0xFF)

            # Update context
            ctx_sum[ctx] += (actual - base)  # error relative to uncorrected
            ctx_count[ctx] += 1

    return bytes(residuals)

# ============================================================
# IDEA 5: CHECKERBOARD DECOMPOSITION
# ============================================================

def encode_checkerboard(plane):
    """Phase 1: encode even pixels (r+c even). Phase 2: encode odd pixels (r+c odd).
    Odd pixels have 4 cardinal neighbors that are ALL even (already known).
    This gives perfect symmetric 2D context — no staircase paradox!"""
    h, w = plane.shape

    # Phase 1: even pixels (r+c % 2 == 0), predicted from known even neighbors
    even_residuals = bytearray()
    even_plane = np.full((h, w), -1, dtype=np.int16)

    # Scan even pixels in raster order, predict from available even neighbors
    for r in range(h):
        for c in range(w):
            if (r + c) % 2 != 0:
                continue
            # Predict from known even neighbors (2-step distance)
            vals = []
            for dr, dc in [(-2,0), (2,0), (0,-2), (0,2), (-1,-1), (-1,1), (1,-1), (1,1)]:
                nr, nc = r+dr, c+dc
                if 0<=nr<h and 0<=nc<w and even_plane[nr,nc] >= 0:
                    vals.append(int(even_plane[nr,nc]))
            pred = sum(vals) // len(vals) if vals else 128
            even_residuals.append((int(plane[r,c]) - pred) & 0xFF)
            even_plane[r,c] = plane[r,c]

    # Phase 2: odd pixels (r+c % 2 == 1), predict from ALL 4 cardinal even neighbors
    odd_residuals = bytearray()
    for r in range(h):
        for c in range(w):
            if (r + c) % 2 != 1:
                continue
            # All 4 cardinal neighbors are even → all known!
            vals = []
            for dr, dc in [(-1,0), (1,0), (0,-1), (0,1)]:
                nr, nc = r+dr, c+dc
                if 0<=nr<h and 0<=nc<w:
                    vals.append(int(plane[nr,nc]))
            pred = sum(vals) // len(vals) if vals else 128
            odd_residuals.append((int(plane[r,c]) - pred) & 0xFF)

    return bytes(even_residuals) + bytes(odd_residuals)

# ============================================================
# IDEA 6: MULTI-RESOLUTION (2x downsample + residual)
# ============================================================

def encode_multiresolution(plane):
    """Encode 2x-downsampled version first, then full-res residual."""
    h, w = plane.shape
    # Downsample by averaging 2x2 blocks
    dh, dw = h // 2, w // 2
    down = np.zeros((dh, dw), dtype=np.uint8)
    for r in range(dh):
        for c in range(dw):
            block = plane[2*r:2*r+2, 2*c:2*c+2].astype(int)
            down[r, c] = block.mean().astype(np.uint8)

    # Upsample back (nearest neighbor for simplicity)
    up = np.repeat(np.repeat(down, 2, axis=0), 2, axis=1)[:h, :w]

    # Residual
    residual = ((plane.astype(int) - up.astype(int)) & 0xFF).astype(np.uint8)

    return down.tobytes() + residual.tobytes()

# ============================================================
# IDEA 7: QUINCUNX SCAN
# ============================================================

def encode_quincunx(plane):
    """Phase 1: encode diamond-pattern pixels. Phase 2: fill gaps.
    Diamond pattern: r%2==0,c%2==0 OR r%2==1,c%2==1.
    Gap pixels: each has 4 neighbors in the diamond (N,E,S,W)."""
    h, w = plane.shape

    # Phase 1: diamond pixels
    diamond_res = bytearray()
    known = np.full((h,w), False)
    for r in range(h):
        for c in range(w):
            if (r+c) % 2 != 0: continue
            # Predict from known diamond neighbors (skip=2 in cardinal dirs)
            vals = []
            for dr, dc in [(-2,0),(2,0),(0,-2),(0,2)]:
                nr, nc = r+dr, c+dc
                if 0<=nr<h and 0<=nc<w and known[nr,nc]:
                    vals.append(int(plane[nr,nc]))
            pred = sum(vals)//len(vals) if vals else 128
            diamond_res.append((int(plane[r,c])-pred)&0xFF)
            known[r,c] = True

    # Phase 2: gap pixels — each has 4 cardinal diamond neighbors
    gap_res = bytearray()
    for r in range(h):
        for c in range(w):
            if (r+c) % 2 != 1: continue
            vals = []
            for dr, dc in [(-1,0),(1,0),(0,-1),(0,1)]:
                nr, nc = r+dr, c+dc
                if 0<=nr<h and 0<=nc<w:
                    vals.append(int(plane[nr,nc]))
            pred = sum(vals)//len(vals) if vals else 128
            gap_res.append((int(plane[r,c])-pred)&0xFF)

    return bytes(diamond_res) + bytes(gap_res)

# ============================================================
# IDEA 4: GRADIENT-STEERED per-pixel prediction
# ============================================================

def encode_gradient_steer(plane):
    """Per-pixel: compute gradient from known neighbors, predict along it."""
    h, w = plane.shape
    residuals = bytearray()

    for r in range(h):
        for c in range(w):
            a = int(plane[r,c-1]) if c > 0 else 128
            b = int(plane[r-1,c]) if r > 0 else 128
            c2 = int(plane[r-1,c-1]) if r > 0 and c > 0 else 128
            d = int(plane[r-1,c+1]) if r > 0 and c+1 < w else b

            # Gradient
            gh = abs(a - c2) + abs(b - d)  # horizontal variation
            gv = abs(b - c2) + abs(a - (int(plane[r,c-2]) if c >= 2 else a))  # vertical

            if gh > 2*gv:     # strong horiz gradient → predict from above
                pred = b
            elif gv > 2*gh:   # strong vert gradient → predict from left
                pred = a
            else:              # smooth → MED
                pred = med(a, b, c2)

            residuals.append((int(plane[r,c]) - pred) & 0xFF)

    return bytes(residuals)

# ============================================================
# BASELINE: Standard MED raster (for comparison)
# ============================================================

def encode_med_raster(plane):
    h, w = plane.shape
    residuals = bytearray()
    for r in range(h):
        for c in range(w):
            a = int(plane[r,c-1]) if c > 0 else 0
            b = int(plane[r-1,c]) if r > 0 else 0
            c2 = int(plane[r-1,c-1]) if r > 0 and c > 0 else 0
            residuals.append((int(plane[r,c]) - med(a, b, c2)) & 0xFF)
    return bytes(residuals)

# ============================================================
# TEST IMAGES
# ============================================================

SZ = 128

def make_tests():
    T = {}; np.random.seed(42)
    x, y = np.meshgrid(np.arange(SZ, dtype=float), np.arange(SZ, dtype=float))
    T["solid"] = np.full((SZ,SZ), 128, dtype=np.uint8)
    T["grad_h"] = np.tile(np.linspace(0,255,SZ,dtype=np.uint8),(SZ,1))
    T["grad_d"] = ((x+y)*255//(2*SZ-2)).astype(np.uint8)
    for deg in [0, 45, 90]:
        th = math.radians(deg)
        proj = x*math.cos(th)+y*math.sin(th)
        T[f"stripes_{deg}"] = ((proj/8).astype(int)%2*255).astype(np.uint8)
    r = np.sqrt((x-SZ/2)**2+(y-SZ/2)**2)
    T["circles"] = (np.sin(r/5)*127+128).astype(np.uint8)
    T["blob"] = (np.exp(-r**2/(2*(SZ/4)**2))*255).astype(np.uint8)
    T["random"] = np.random.randint(0,256,(SZ,SZ),dtype=np.uint8)
    sm = np.random.randint(0,256,(SZ//16,SZ//16),dtype=np.uint8)
    T["natural"] = np.clip(np.array(Image.fromarray(sm).resize((SZ,SZ),Image.BILINEAR)).astype(float)
                           +np.random.normal(0,10,(SZ,SZ)),0,255).astype(np.uint8)
    # Mixed directions
    mixed = np.zeros((SZ,SZ),dtype=np.uint8)
    h2 = SZ//2
    mixed[:h2,:h2] = T["stripes_0"][:h2,:h2]
    mixed[:h2,h2:] = T["stripes_90"][:h2,h2:]
    mixed[h2:,:h2] = T["stripes_45"][h2:,:h2]
    mixed[h2:,h2:] = T["circles"][h2:,h2:]
    T["mixed"] = mixed
    return T

# ============================================================
# BENCHMARK ALL IDEAS
# ============================================================

print("=" * 90)
print("  PARALLEL IDEAS: Kill two birds with one stone")
print("  kind-pasteur-2026-03-25-S7")
print("=" * 90)

METHODS = [
    ("MED raster",      encode_med_raster),
    ("Dual-pass",       encode_dualpass),
    ("Context hash",    encode_context_hash),
    ("Gradient steer",  encode_gradient_steer),
    ("Checkerboard",    encode_checkerboard),
    ("Multi-res",       encode_multiresolution),
    ("Quincunx",        encode_quincunx),
    ("Raw",             lambda p: p.tobytes()),
]

tests = make_tests()

# Header for results table
print(f"\n  {'Image':<12}", end="")
for mname, _ in METHODS:
    print(f" {mname:>12}", end="")
print(f"  {'PNG':>6}  {'BEST':>12}  {'vs PNG':>6}")
print("  " + "-" * (14 + 13*len(METHODS) + 30))

overall_best_method_wins = {m: 0 for m, _ in METHODS}

for name, img in sorted(tests.items()):
    ps = png_size(img)
    sizes = {}
    for mname, mfunc in METHODS:
        res = mfunc(img)
        csz = len(compress(res)) + 10  # header estimate
        sizes[mname] = csz

    best_method = min(sizes, key=sizes.get)
    best_sz = sizes[best_method]
    overall_best_method_wins[best_method] += 1

    print(f"  {name:<12}", end="")
    for mname, _ in METHODS:
        marker = "*" if mname == best_method else " "
        print(f" {sizes[mname]:>11}{marker}", end="")
    ratio = ps / best_sz if best_sz > 0 else 0
    print(f"  {ps:>6}  {best_method:>12}  {ratio:>5.2f}x")

print(f"\n  METHOD WIN COUNTS (which method was smallest most often):")
for m, count in sorted(overall_best_method_wins.items(), key=lambda x: -x[1]):
    if count > 0:
        print(f"    {m:<16}: {count}")

# Detailed analysis of the best new ideas
print(f"\n{'='*90}")
print("  ANALYSIS: Which new idea actually helps?")
print(f"{'='*90}")

print(f"\n  Compare each method to MED raster baseline:")
print(f"  {'Image':<12} {'MED':>6}", end="")
for mname, _ in METHODS[1:]:
    print(f" {mname[:8]:>9}", end="")
print()

for name, img in sorted(tests.items()):
    base_sz = len(compress(encode_med_raster(img))) + 10
    print(f"  {name:<12} {base_sz:>6}", end="")
    for mname, mfunc in METHODS[1:]:
        sz = len(compress(mfunc(img))) + 10
        ratio = base_sz / sz if sz > 0 else 0
        if ratio > 1.02: marker = "+"
        elif ratio < 0.98: marker = "-"
        else: marker = "="
        print(f" {ratio:>7.2f}{marker}", end="")
    print()

print(f"""
  LEGEND: + = better than MED, - = worse, = = same (within 2%)

  KEY FINDINGS:
  - Context hash: expected to help on natural images (running error correction)
  - Checkerboard: expected to help on diagonal patterns (symmetric 2D context)
  - Multi-res: expected to help on smooth images (global structure first)
  - Gradient steer: expected to help on mixed-direction images
  - Quincunx: similar to checkerboard but diamond geometry
""")
