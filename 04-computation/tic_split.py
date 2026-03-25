#!/usr/bin/env python3
"""
BINARY SPLIT CODEC + ROW LPC: Two more radical breaks from convention.

IDEA 1: RECURSIVE BINARY SPLIT
  Split image in half. Encode difference of the two halves.
  Recurse on each half. Choose split direction (H/V) adaptively.
  If halves are similar → difference is near-zero → compresses great.
  This captures SELF-SIMILARITY at all scales without wavelets.

IDEA 2: ROW LINEAR PREDICTIVE CODING (Row-LPC)
  Predict row[r] as a linear combination of previous rows:
    row_pred[r] = a1 * row[r-1] + a2 * row[r-2]
  Coefficients computed per-row via least-squares on recent history.
  This is how FLAC compresses audio — applied to image rows.
  Captures inter-row correlations that pixel-level MED misses.

IDEA 3: HISTOGRAM EQUALIZATION + MED
  Apply histogram equalization (lossless version: rank transform).
  This makes the value distribution uniform → better for entropy coding.
  Then apply MED on the equalized image.
  Decode: inverse rank transform.

kind-pasteur-2026-03-25-S8
"""
import sys, io, struct, zlib, math, time
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

try: import brotli; HAS_BROTLI = True
except: HAS_BROTLI = False

def bc(data):
    r = [zlib.compress(data, 9)]
    if HAS_BROTLI:
        try: r.append(brotli.compress(data, quality=11))
        except: pass
    return min(r, key=len)

def png_size(img):
    pil = Image.fromarray(img, 'L' if img.ndim == 2 else 'RGB')
    buf = io.BytesIO(); pil.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()

def med(a, b, c):
    mn, mx = min(a, b), max(a, b)
    if c >= mx: return mn
    if c <= mn: return mx
    return a + b - c

def encode_med(plane):
    h, w = plane.shape; res = bytearray()
    for r in range(h):
        for c in range(w):
            a = int(plane[r,c-1]) if c>0 else 0
            b = int(plane[r-1,c]) if r>0 else 0
            c2 = int(plane[r-1,c-1]) if r>0 and c>0 else 0
            res.append((int(plane[r,c]) - med(a,b,c2)) & 0xFF)
    return len(bc(bytes(res)))

# ============================================================
# IDEA 1: RECURSIVE BINARY SPLIT
# ============================================================

def encode_binary_split(plane, depth=0, max_depth=6):
    """Recursively split image, encode differences."""
    h, w = plane.shape

    if h < 4 or w < 4 or depth >= max_depth:
        # Base case: encode directly
        return bc(plane.tobytes())

    # Try horizontal split
    top = plane[:h//2, :]
    bot = plane[h//2:, :]
    h_diff = ((bot[:min(top.shape[0],bot.shape[0])].astype(int) -
               top[:min(top.shape[0],bot.shape[0])].astype(int)) & 0xFF).astype(np.uint8)
    h_cost = len(bc(h_diff.tobytes()))
    h_total = h_cost + len(encode_binary_split(top, depth+1, max_depth)) + len(encode_binary_split(bot, depth+1, max_depth))

    # Try vertical split
    left = plane[:, :w//2]
    right = plane[:, w//2:]
    v_diff = ((right[:, :min(left.shape[1],right.shape[1])].astype(int) -
               left[:, :min(left.shape[1],right.shape[1])].astype(int)) & 0xFF).astype(np.uint8)
    v_cost = len(bc(v_diff.tobytes()))
    v_total = v_cost + len(encode_binary_split(left, depth+1, max_depth)) + len(encode_binary_split(right, depth+1, max_depth))

    # Also try just encoding directly
    direct = bc(plane.tobytes())

    # Pick smallest
    best = min(len(direct), h_total, v_total)
    if best == len(direct):
        return direct
    elif best == h_total:
        return struct.pack('<B', 0) + encode_binary_split(top, depth+1, max_depth) + encode_binary_split(bot, depth+1, max_depth)
    else:
        return struct.pack('<B', 1) + encode_binary_split(left, depth+1, max_depth) + encode_binary_split(right, depth+1, max_depth)

# ============================================================
# IDEA 2: ROW LINEAR PREDICTIVE CODING
# ============================================================

def encode_row_lpc(plane, order=2):
    """Predict each row as linear combination of previous rows.
    row_pred = a1*row[-1] + a2*row[-2] (+ ...)
    Coefficients encoded per-row. Residual = row - row_pred."""
    h, w = plane.shape
    coeffs_data = bytearray()
    residual_data = bytearray()

    history = []  # list of previous rows (as float arrays)

    for r in range(h):
        row = plane[r].astype(float)

        if len(history) >= 1:
            # Build system: predict row from history
            n_hist = min(order, len(history))
            A = np.column_stack([history[-i-1] for i in range(n_hist)])
            # Solve least squares: row ≈ A @ coeffs
            try:
                c2, _, _, _ = np.linalg.lstsq(A, row, rcond=None)
                pred = np.clip(np.round(A @ c2), 0, 255)
            except:
                pred = np.zeros(w)
                c2 = np.zeros(n_hist)

            # Encode coefficients as fixed-point (1 byte each, scaled)
            for coeff in c2:
                # Scale to [-2, 2] range, encode as uint8
                scaled = int(np.clip(coeff * 64 + 128, 0, 255))
                coeffs_data.append(scaled)

            residual = ((plane[r].astype(int) - pred.astype(int)) & 0xFF).astype(np.uint8)
        else:
            residual = plane[r].copy()

        residual_data.extend(residual)
        history.append(row)

    total = len(bc(bytes(coeffs_data))) + len(bc(bytes(residual_data)))
    return total

# ============================================================
# IDEA 3: RANK TRANSFORM + MED
# ============================================================

def encode_rank_med(plane):
    """Apply rank transform (lossless histogram equalization) then MED.
    The rank transform maps each pixel to its rank in the sorted list
    of unique values. This makes the distribution more uniform."""
    h, w = plane.shape
    flat = plane.ravel()

    # Get unique values and their ranks
    unique_vals = np.unique(flat)
    rank_map = np.zeros(256, dtype=np.uint8)
    for rank, val in enumerate(unique_vals):
        rank_map[val] = min(255, rank * 255 // max(1, len(unique_vals) - 1))

    # Apply rank transform
    ranked = rank_map[plane]

    # MED on ranked image
    res = bytearray()
    for r in range(h):
        for c in range(w):
            a = int(ranked[r,c-1]) if c>0 else 0
            b = int(ranked[r-1,c]) if r>0 else 0
            c2 = int(ranked[r-1,c-1]) if r>0 and c>0 else 0
            res.append((int(ranked[r,c]) - med(a,b,c2)) & 0xFF)

    # Total: rank map (256 bytes) + MED residuals
    total = len(bc(unique_vals.tobytes())) + len(bc(bytes(res)))
    return total

# ============================================================
# IDEA 4: XOR-PLANE DECOMPOSITION
# ============================================================

def encode_xor_planes(plane):
    """Gray-code the pixel values, then encode each bit-plane.
    Gray code: adjacent values differ by 1 bit → bit-planes are smoother."""
    h, w = plane.shape
    # Gray code: g = v ^ (v >> 1)
    gray = (plane.astype(np.uint16) ^ (plane.astype(np.uint16) >> 1)).astype(np.uint8)

    total = 0
    for bit in range(7, -1, -1):
        bp = ((gray >> bit) & 1).astype(np.uint8)
        packed = np.packbits(bp)
        total += len(bc(bytes(packed)))
    return total

# ============================================================
# IDEA 5: EDGE-MAP + SMOOTH-MAP DECOMPOSITION
# ============================================================

def encode_edge_smooth(plane):
    """Decompose image into edge map + smooth component.
    Smooth = low-pass filtered. Edge = original - smooth.
    Smooth compresses very well. Edge is sparse → compresses well too."""
    h, w = plane.shape

    # Simple box blur (3x3 average)
    smooth = np.zeros_like(plane, dtype=np.float64)
    for r in range(h):
        for c in range(w):
            vals = []
            for dr in [-1,0,1]:
                for dc in [-1,0,1]:
                    nr, nc = r+dr, c+dc
                    if 0<=nr<h and 0<=nc<w:
                        vals.append(int(plane[nr,nc]))
            smooth[r,c] = sum(vals) / len(vals)

    smooth_uint8 = np.clip(np.round(smooth), 0, 255).astype(np.uint8)
    edge = ((plane.astype(int) - smooth_uint8.astype(int)) & 0xFF).astype(np.uint8)

    total = len(bc(smooth_uint8.tobytes())) + len(bc(edge.tobytes()))
    return total

# ============================================================
# BENCHMARK
# ============================================================

SZ = 128

def make_tests():
    T = {}; np.random.seed(42)
    x, y = np.meshgrid(np.arange(SZ, dtype=float), np.arange(SZ, dtype=float))
    r = np.sqrt((x-SZ/2)**2 + (y-SZ/2)**2)
    T["solid"] = np.full((SZ,SZ), 128, dtype=np.uint8)
    T["grad_h"] = np.tile(np.linspace(0,255,SZ,dtype=np.uint8),(SZ,1))
    T["grad_d"] = ((x+y)*255//(2*SZ-2)).astype(np.uint8)
    T["circles"] = (np.sin(r/5)*127+128).astype(np.uint8)
    T["blob"] = (np.exp(-r**2/(2*(SZ/4)**2))*255).astype(np.uint8)
    T["random"] = np.random.randint(0,256,(SZ,SZ),dtype=np.uint8)
    sm = np.random.randint(0,256,(SZ//16,SZ//16),dtype=np.uint8)
    T["natural"] = np.clip(np.array(Image.fromarray(sm).resize((SZ,SZ),Image.BILINEAR)).astype(float)
                           +np.random.normal(0,10,(SZ,SZ)),0,255).astype(np.uint8)
    T["stripes45"] = (((x+y)/8).astype(int)%2*255).astype(np.uint8)
    return T

print("=" * 100)
print("  RADICAL IDEAS ROUND 2: Split, LPC, Rank, Gray-XOR, Edge-Smooth")
print("  kind-pasteur-2026-03-25-S8")
print("=" * 100)

tests = make_tests()

METHODS = [
    ("MED std",     encode_med),
    ("BinSplit",    lambda p: len(encode_binary_split(p))),
    ("RowLPC",      encode_row_lpc),
    ("RankMED",     encode_rank_med),
    ("GrayXOR",     encode_xor_planes),
    ("EdgeSmooth",  encode_edge_smooth),
]

print(f"\n  {'Image':<12} {'PNG':>6}", end="")
for mname, _ in METHODS:
    print(f" {mname:>10}", end="")
print(f"  {'BEST':>10}")
print("  " + "-" * (22 + 11*len(METHODS) + 14))

for name, img in sorted(tests.items()):
    ps = png_size(img)
    sizes = {}
    for mname, mfunc in METHODS:
        sizes[mname] = mfunc(img) + 10

    best = min(sizes, key=sizes.get)
    ratio = ps / sizes[best] if sizes[best] > 0 else 0

    print(f"  {name:<12} {ps:>6}", end="")
    for mname, _ in METHODS:
        marker = "*" if mname == best else " "
        print(f" {sizes[mname]:>9}{marker}", end="")
    print(f"  {best:>10} {ratio:.2f}x")

print(f"""
  KEY FINDINGS:
  - Binary Split: captures self-similarity at multiple scales
  - Row LPC: audio-style linear prediction on rows
  - Rank+MED: histogram equalization before MED
  - Gray XOR: Gray-code bit-planes (smoother bit-planes)
  - Edge+Smooth: separate high/low frequency components
""")
