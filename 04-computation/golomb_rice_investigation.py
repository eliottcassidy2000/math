#!/usr/bin/env python3
"""
GOLOMB-RICE INVESTIGATION: Can we compete with JPEG-LS?

JPEG-LS uses: MED prediction + Golomb-Rice coding + context modeling.
We use:       MED prediction + zlib (Deflate).

The gap is in the ENTROPY CODER, not the predictor.

WHY GOLOMB-RICE BEATS ZLIB ON MED RESIDUALS:

1. Golomb-Rice is OPTIMAL for geometric/Laplacian distributions.
   MED residuals ARE approximately Laplacian (peaked at 0, exponential tails).
   Golomb-Rice encodes each residual in ~log₂(1 + |r|/k) + k bits.
   For k=4 and residual 0: 1 bit. For residual 3: 4 bits.

2. Zlib (Deflate) uses LZ77 + Huffman. It finds REPEATED PATTERNS in
   the byte stream. For i.i.d. Laplacian residuals, there are no repeated
   patterns — each residual is independent. Zlib falls back to pure Huffman,
   which is suboptimal for Laplacian distributions (byte-level, not bit-level).

3. The gap: Golomb-Rice operates at BIT level (variable-length per residual).
   Zlib operates at BYTE level (fixed 8-bit symbols). For small residuals
   (0, ±1, ±2), Golomb-Rice uses 1-3 bits. Zlib uses 8 bits minimum.

THIS IS THE FUNDAMENTAL BOTTLENECK: zlib can't encode a residual of 0 in
less than ~1 byte. Golomb-Rice can encode it in 1 bit.

CAN WE BUILD GOLOMB-RICE IN PYTHON? Yes, it's simple:
  1. Map signed residual to unsigned via zigzag: 0→0, -1→1, 1→2, -2→3, ...
  2. Divide by 2^k: quotient q, remainder r.
  3. Encode: q ones + 1 zero (unary) + r in k bits (binary).

The parameter k is chosen per-row or per-block to minimize total bits.

kind-pasteur-2026-03-25-S31
"""
import sys, io, zlib, math, time, os, glob
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

try: import brotli; HAS_BROTLI = True
except: HAS_BROTLI = False

def png_size(img):
    if isinstance(img, np.ndarray):
        pil = Image.fromarray(img, 'L' if img.ndim == 2 else 'RGB')
    else:
        pil = img
    buf = io.BytesIO(); pil.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()

def med(a, b, c):
    mn, mx = min(a, b), max(a, b)
    if c >= mx: return mn
    if c <= mn: return mx
    return a + b - c

# ============================================================
# GOLOMB-RICE ENCODER/DECODER
# ============================================================

def zigzag_encode(signed):
    """Map signed to unsigned: 0→0, -1→1, 1→2, -2→3, 2→4, ..."""
    if signed >= 0: return 2 * signed
    return 2 * (-signed) - 1

def zigzag_decode(unsigned):
    """Inverse of zigzag_encode."""
    if unsigned % 2 == 0: return unsigned // 2
    return -(unsigned + 1) // 2

def golomb_rice_encode_stream(values, k):
    """Encode a stream of unsigned integers with Golomb-Rice parameter k.
    Returns a bytearray of packed bits."""
    bits = []
    for v in values:
        q = v >> k  # quotient
        r = v & ((1 << k) - 1)  # remainder
        # Unary: q ones + 1 zero
        bits.extend([1] * q)
        bits.append(0)
        # Binary: k bits for remainder
        for i in range(k-1, -1, -1):
            bits.append((r >> i) & 1)

    # Pack to bytes
    n_bytes = (len(bits) + 7) // 8
    out = bytearray(n_bytes)
    for i, b in enumerate(bits):
        if b: out[i >> 3] |= (1 << (7 - (i & 7)))
    return bytes(out), len(bits)

def golomb_rice_decode_stream(data, n_bits, n_values, k):
    """Decode n_values unsigned integers from Golomb-Rice bitstream."""
    values = []
    bit_pos = 0

    def read_bit():
        nonlocal bit_pos
        byte_idx = bit_pos >> 3
        bit_idx = 7 - (bit_pos & 7)
        bit_pos += 1
        return (data[byte_idx] >> bit_idx) & 1

    for _ in range(n_values):
        # Read unary: count ones until zero
        q = 0
        while read_bit() == 1:
            q += 1
        # Read binary: k bits
        r = 0
        for _ in range(k):
            r = (r << 1) | read_bit()
        values.append((q << k) | r)

    return values

def optimal_k(values):
    """Find the optimal Golomb-Rice parameter k for a set of values.
    Optimal k ≈ max(0, round(log2(mean/ln(2))))."""
    if not values: return 0
    mean_val = sum(values) / len(values)
    if mean_val < 1: return 0
    k = max(0, round(math.log2(mean_val / math.log(2))))
    return min(k, 15)  # cap at 15

# ============================================================
# FULL CODEC: G-RG-BG + MED + Golomb-Rice
# ============================================================

def encode_plane_gr(plane):
    """Encode one plane with MED prediction + Golomb-Rice coding."""
    h, w = plane.shape

    # MED prediction → signed residuals
    residuals = []
    for r in range(h):
        for c in range(w):
            a = int(plane[r,c-1]) if c > 0 else 0
            b = int(plane[r-1,c]) if r > 0 else 0
            d = int(plane[r-1,c-1]) if r > 0 and c > 0 else 0
            pred = med(a, b, d)
            res = int(plane[r, c]) - pred
            # Map to signed in [-128, 127]
            if res > 128: res -= 256
            if res < -128: res += 256
            residuals.append(res)

    # Zigzag encode → unsigned
    unsigned = [zigzag_encode(r) for r in residuals]

    # Find optimal k
    k = optimal_k(unsigned)

    # Golomb-Rice encode
    gr_data, n_bits = golomb_rice_encode_stream(unsigned, k)

    return gr_data, n_bits, k, residuals

def encode_plane_zlib(plane):
    """Encode one plane with MED prediction + zlib."""
    h, w = plane.shape
    res = bytearray()
    for r in range(h):
        for c in range(w):
            a = int(plane[r,c-1]) if c > 0 else 0
            b = int(plane[r-1,c]) if r > 0 else 0
            d = int(plane[r-1,c-1]) if r > 0 and c > 0 else 0
            res.append((int(plane[r,c]) - med(a, b, d)) & 0xFF)
    return zlib.compress(bytes(res), 9)

# ============================================================
# BENCHMARK: Golomb-Rice vs zlib on real images
# ============================================================

print("=" * 80)
print("  GOLOMB-RICE vs ZLIB: Why JPEG-LS beats our codec")
print("  kind-pasteur-2026-03-25-S31")
print("=" * 80)

# Test on Kodak + real photo
test_files = sorted(glob.glob("test_images/kodim*.png"))
test_files.append("inbox/vlcsnap-2026-03-10-17h12m34s408.png")

print(f"\n  {'Image':<18} {'PNG':>8} {'MED+zlib':>10} {'MED+GR':>10} {'GR bits':>10} "
      f"{'Shannon':>10} {'k':>3} {'GR/zlib':>8}")
print("  " + "-" * 85)

total_png = total_zlib = total_gr = total_shannon = 0

for fpath in test_files:
    name = os.path.basename(fpath)[:16]
    pil = Image.open(fpath).convert('RGB')
    w0, h0 = pil.size
    crop = pil.crop((w0//2-128, h0//2-128, w0//2+128, h0//2+128))

    # Work on green channel (representative)
    green = np.array(crop)[:,:,1]

    ps = png_size(Image.fromarray(green, 'L'))

    # MED + zlib
    zlib_data = encode_plane_zlib(green)
    zlib_sz = len(zlib_data)

    # MED + Golomb-Rice
    gr_data, gr_bits, k, residuals = encode_plane_gr(green)
    gr_bytes = (gr_bits + 7) // 8 + 1  # +1 for k parameter

    # Shannon entropy (theoretical minimum)
    unsigned = [zigzag_encode(r) for r in residuals]
    from collections import Counter
    counts = Counter(unsigned)
    total_n = len(unsigned)
    shannon_bits = -sum(c * math.log2(c / total_n) for c in counts.values())
    shannon_bytes = int(math.ceil(shannon_bits / 8))

    gr_ratio = zlib_sz / gr_bytes if gr_bytes > 0 else 0

    total_png += ps; total_zlib += zlib_sz; total_gr += gr_bytes; total_shannon += shannon_bytes

    print(f"  {name:<18} {ps:>8} {zlib_sz:>10} {gr_bytes:>10} {gr_bits:>10} "
          f"{shannon_bytes:>10} {k:>3} {gr_ratio:>7.3f}x")

# RGB totals (3 planes)
print(f"\n  {'AGGREGATE (G only)':<18} {total_png:>8} {total_zlib:>10} {total_gr:>10} "
      f"{'':>10} {total_shannon:>10}")

agg_gr_vs_zlib = total_zlib / total_gr if total_gr > 0 else 0
agg_gr_vs_png = total_png / total_gr if total_gr > 0 else 0
agg_shannon_vs_gr = total_gr / total_shannon if total_shannon > 0 else 0

print(f"""
  ANALYSIS:
    GR vs zlib:    {agg_gr_vs_zlib:.3f}x ({'GR better' if agg_gr_vs_zlib > 1 else 'zlib better'})
    GR vs PNG:     {agg_gr_vs_png:.3f}x
    GR vs Shannon: {agg_shannon_vs_gr:.3f}x from optimal

  WHY THE GAP EXISTS:
    zlib encodes at BYTE level (min 1 byte per symbol).
    Golomb-Rice encodes at BIT level (can use <8 bits per symbol).
    For MED residuals peaked at 0: most residuals are small (0, ±1, ±2).
    GR encodes 0 in {k+1} bits. At k=4: 5 bits. At k=2: 3 bits.
    zlib encodes 0 in ~8 bits minimum (one byte).

  THE FIX: Replace zlib with Golomb-Rice in our C codec.
    This would make TIC competitive with JPEG-LS.
    Golomb-Rice IS lossless (it's just a different way to encode integers).
    JPEG-LS adds context modeling on top (per-context k selection),
    which gives another ~5% over plain Golomb-Rice.
""")

# Verify Golomb-Rice roundtrip
print("  ROUNDTRIP VERIFICATION:")
for k_test in [2, 4, 6]:
    test_vals = list(range(20)) + [100, 200, 255]
    encoded, n_bits = golomb_rice_encode_stream(test_vals, k_test)
    decoded = golomb_rice_decode_stream(encoded, n_bits, len(test_vals), k_test)
    ok = (test_vals == decoded)
    print(f"    k={k_test}: encode {len(test_vals)} values → {n_bits} bits → decode {'PASS' if ok else 'FAIL'}")
