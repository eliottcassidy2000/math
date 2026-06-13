#!/usr/bin/env python3
"""
RADICAL IDEAS: Break every assumption about image compression.

ASSUMPTION 1: "I must scan pixels in some order"
  BREAK IT: Don't scan. Describe the image as a SET of constraints.
  Idea: express the image as a system of equations that a solver can
  reconstruct. E.g., "row 5 has sum 1247" + "column 3 has sum 983" +
  enough constraints to uniquely determine the image.
  This is like TOMOGRAPHIC compression — project along multiple angles.

ASSUMPTION 2: "I must predict pixel from neighbors"
  BREAK IT: Predict BLOCKS from OTHER BLOCKS, or predict the entire
  image from a TEMPLATE. If the image is "gradient + noise", encode
  the gradient parameters (2 numbers) + the noise separately.

ASSUMPTION 3: "Residuals go to a general-purpose compressor"
  BREAK IT: If residuals are approximately Laplacian distributed
  (which they are for MED), use a CUSTOM entropy coder optimized
  for Laplacian distributions. Golomb-Rice codes are optimal for
  geometric distributions. Let's USE that.

ASSUMPTION 4: "The image is a grid of pixels"
  BREAK IT: Treat the image as a FUNCTION f(x,y) and find its most
  compact representation. Fourier? No — too many coefficients.
  POLYNOMIAL? Maybe for smooth regions. PIECEWISE LINEAR? Yes —
  triangulate the image, encode vertices + values, interpolate.

ASSUMPTION 5: "Compress the image as given"
  BREAK IT: TRANSFORM the image into a domain where it's simpler.
  Not DCT (lossy) — but what about the WAVELET domain? Lossless
  wavelet transform (like JPEG 2000's 5/3 wavelet) makes the image
  into a sparse representation that compresses better.

ASSUMPTION 6: "Each pixel is independent of distant pixels"
  BREAK IT: Find GLOBAL patterns. If row 10 is similar to row 50,
  use a COPY reference (like LZ77 but in 2D). This is what zlib does
  in 1D — extend it to 2D matching.

Let me try the most radical ideas and measure.

kind-pasteur-2026-03-25-S8
"""
import sys, io, struct, zlib, math, time
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

try: import brotli; HAS_BROTLI = True
except: HAS_BROTLI = False

def bc(data):
    """Best compress."""
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
# RADICAL IDEA 1: GOLOMB-RICE CODING of MED residuals
# ============================================================

def golomb_rice_encode(residuals, k=4):
    """Encode residuals using Golomb-Rice code with parameter k.
    Optimal for geometric/Laplacian distributed data.
    Maps signed residuals to unsigned via zigzag encoding,
    then encodes with unary(q) + binary(r) where q = val >> k, r = val & ((1<<k)-1)."""
    # Zigzag encode: 0 -> 0, -1 -> 1, 1 -> 2, -2 -> 3, 2 -> 4, ...
    bits = []
    for res in residuals:
        # Convert from unsigned mod-256 to signed
        signed = res if res <= 128 else res - 256
        # Zigzag
        zz = (signed << 1) ^ (signed >> 31) if signed >= 0 else ((-signed) << 1) - 1
        zz = abs(signed) * 2 - (1 if signed > 0 else 0) if signed != 0 else 0
        # Simpler zigzag
        if signed >= 0:
            zz = 2 * signed
        else:
            zz = 2 * (-signed) - 1
        q = zz >> k
        r = zz & ((1 << k) - 1)
        # Unary: q ones followed by a zero
        bits.extend([1] * q)
        bits.append(0)
        # Binary: k bits for remainder
        for i in range(k-1, -1, -1):
            bits.append((r >> i) & 1)

    # Pack bits into bytes
    n_bytes = (len(bits) + 7) // 8
    out = bytearray(n_bytes)
    for i, b in enumerate(bits):
        if b:
            out[i >> 3] |= (1 << (7 - (i & 7)))
    return bytes(out), len(bits)

def golomb_rice_size(plane, k=None):
    """Compute Golomb-Rice coded size of MED residuals."""
    h, w = plane.shape
    residuals = []
    for r in range(h):
        for c in range(w):
            a = int(plane[r,c-1]) if c>0 else 0
            b = int(plane[r-1,c]) if r>0 else 0
            c2 = int(plane[r-1,c-1]) if r>0 and c>0 else 0
            pred = med(a, b, c2)
            residuals.append((int(plane[r,c]) - pred) & 0xFF)

    if k is None:
        # Try multiple k values, pick best
        best_bits = float('inf')
        best_k = 0
        for kk in range(0, 8):
            _, nbits = golomb_rice_encode(residuals, kk)
            if nbits < best_bits:
                best_bits = nbits
                best_k = kk
        return (best_bits + 7) // 8 + 1, best_k  # +1 for k parameter
    else:
        _, nbits = golomb_rice_encode(residuals, k)
        return (nbits + 7) // 8 + 1, k

# ============================================================
# RADICAL IDEA 2: PIECEWISE-FLAT REGIONS (run-length on regions)
# ============================================================

def encode_flat_regions(plane, threshold=3):
    """Segment image into nearly-flat regions. Encode each region by
    its bounding box + mean value. Residual is only the per-pixel
    deviation from the region mean.
    This is NOT like any standard codec — it's more like vector graphics."""
    h, w = plane.shape
    visited = np.zeros((h, w), dtype=bool)
    regions = []

    for r in range(h):
        for c in range(w):
            if visited[r, c]:
                continue
            # Flood-fill from (r,c) to find flat region
            val = int(plane[r, c])
            stack = [(r, c)]
            region = []
            while stack:
                rr, cc = stack.pop()
                if rr < 0 or rr >= h or cc < 0 or cc >= w:
                    continue
                if visited[rr, cc]:
                    continue
                if abs(int(plane[rr, cc]) - val) <= threshold:
                    visited[rr, cc] = True
                    region.append((rr, cc))
                    stack.extend([(rr-1,cc),(rr+1,cc),(rr,cc-1),(rr,cc+1)])

            if region:
                mean_val = sum(int(plane[rr2,cc2]) for rr2,cc2 in region) // len(region)
                regions.append((mean_val, region))

    # Encode: number of regions + per-region (mean + bitmap of which pixels)
    # This is compact when regions are large (smooth images) and expensive
    # when regions are small (noisy images).
    # For measurement: encode region means + residuals
    region_data = bytearray()
    residual_data = bytearray()
    region_map = np.zeros((h,w), dtype=np.int16)

    for rid, (mean_val, region) in enumerate(regions):
        region_data.append(mean_val & 0xFF)
        for rr, cc in region:
            region_map[rr, cc] = rid
            residual_data.append((int(plane[rr, cc]) - mean_val) & 0xFF)

    # Total size: region means + residuals (compressed)
    total = len(bc(bytes(region_data))) + len(bc(bytes(residual_data)))
    return total, len(regions)

# ============================================================
# RADICAL IDEA 3: POLYNOMIAL FIT + RESIDUAL
# ============================================================

def encode_polyfit(plane, degree=2):
    """Fit a 2D polynomial to the image, encode coefficients + residual.
    For smooth images, the polynomial captures global structure in ~10 numbers.
    Residual is just local noise/detail."""
    h, w = plane.shape
    x = np.arange(w, dtype=float) / w
    y = np.arange(h, dtype=float) / h
    X, Y = np.meshgrid(x, y)

    # Build feature matrix for polynomial of given degree
    features = []
    for dy in range(degree + 1):
        for dx in range(degree + 1 - dy):
            features.append((X**dx * Y**dy).ravel())
    A = np.column_stack(features)

    # Least squares fit
    vals = plane.ravel().astype(float)
    coeffs, _, _, _ = np.linalg.lstsq(A, vals, rcond=None)

    # Reconstruct polynomial surface
    surface = (A @ coeffs).reshape(h, w)
    surface_uint8 = np.clip(np.round(surface), 0, 255).astype(np.uint8)

    # Residual
    residual = ((plane.astype(int) - surface_uint8.astype(int)) & 0xFF).astype(np.uint8)

    # Size: coefficients (as float16 = 2 bytes each) + compressed residual
    n_coeffs = len(coeffs)
    coeff_bytes = np.array(coeffs, dtype=np.float32).tobytes()  # 4 bytes each
    residual_compressed = bc(residual.tobytes())

    total = len(coeff_bytes) + len(residual_compressed)
    return total, n_coeffs, np.mean(np.abs(residual.astype(np.int8)))

# ============================================================
# RADICAL IDEA 4: 2D LZ MATCHING (block-copy references)
# ============================================================

def encode_2d_lz(plane, block_size=8):
    """Find repeated blocks in the image. Encode novel blocks directly,
    repeated blocks as references (x, y offset to match).
    This is 2D pattern matching — NOT done by any standard image codec."""
    h, w = plane.shape
    bh = (h + block_size - 1) // block_size
    bw = (w + block_size - 1) // block_size

    blocks = {}
    references = []
    novel_data = bytearray()
    n_novel = 0
    n_ref = 0

    for by in range(bh):
        for bx in range(bw):
            r0, c0 = by * block_size, bx * block_size
            r1, c1 = min(r0 + block_size, h), min(c0 + block_size, w)
            block = plane[r0:r1, c0:c1]

            # Normalize block size
            if block.shape != (block_size, block_size):
                padded = np.zeros((block_size, block_size), dtype=np.uint8)
                padded[:block.shape[0], :block.shape[1]] = block
                block = padded

            block_key = block.tobytes()

            # Check for exact match in previous blocks
            found = False
            if block_key in blocks:
                # Reference to first occurrence
                ref_by, ref_bx = blocks[block_key]
                references.append((1, ref_by, ref_bx))  # 1 = reference
                n_ref += 1
                found = True

            if not found:
                # Check for close match (within threshold)
                best_dist = float('inf')
                best_ref = None
                for key, (pby, pbx) in blocks.items():
                    prev_block = np.frombuffer(key, dtype=np.uint8).reshape(block_size, block_size)
                    dist = np.sum(np.abs(block.astype(int) - prev_block.astype(int)))
                    if dist < best_dist:
                        best_dist = dist
                        best_ref = (pby, pbx)

                if best_dist < block_size * block_size * 2 and best_ref is not None:
                    # Close match: encode as reference + small residual
                    ref_block_key = None
                    for key, (pby, pbx) in blocks.items():
                        if (pby, pbx) == best_ref:
                            ref_block_key = key
                            break
                    if ref_block_key:
                        prev_block = np.frombuffer(ref_block_key, dtype=np.uint8).reshape(block_size, block_size)
                        residual = ((block.astype(int) - prev_block.astype(int)) & 0xFF).astype(np.uint8)
                        references.append((2, best_ref[0], best_ref[1]))  # 2 = ref+residual
                        novel_data.extend(residual.tobytes())
                        n_ref += 1
                    else:
                        references.append((0,))  # 0 = novel
                        novel_data.extend(block.tobytes())
                        n_novel += 1
                else:
                    references.append((0,))  # 0 = novel
                    novel_data.extend(block.tobytes())
                    n_novel += 1

            blocks[block_key] = (by, bx)

    # Encode reference map + novel data
    ref_bytes = bytearray()
    for ref in references:
        ref_bytes.append(ref[0])
        if ref[0] > 0:
            ref_bytes.extend(struct.pack('<BB', ref[1], ref[2]))

    total = len(bc(bytes(ref_bytes))) + len(bc(bytes(novel_data)))
    return total, n_novel, n_ref

# ============================================================
# RADICAL IDEA 5: BIT-PLANE DECOMPOSITION
# ============================================================

def encode_bitplanes(plane):
    """Decompose image into 8 bit-planes (MSB to LSB).
    Higher bit-planes are smoother → compress better.
    Lower bit-planes are noisier → compress less.
    This separates signal from noise without wavelets."""
    h, w = plane.shape
    total = 0
    for bit in range(7, -1, -1):
        bitplane = ((plane >> bit) & 1).astype(np.uint8)
        # Pack bits into bytes
        packed = np.packbits(bitplane)
        total += len(bc(bytes(packed)))
    return total

# ============================================================
# RADICAL IDEA 6: SORT-TRANSFORM (BWT-inspired)
# ============================================================

def encode_sort_transform(plane):
    """Sort rows by similarity, then compress. Similar rows become adjacent,
    making row-differences smaller. Like BWT for images.
    This breaks the spatial order but creates compressible sequences."""
    h, w = plane.shape
    # Sort rows by their hash (groups similar rows together)
    rows = [plane[r] for r in range(h)]
    indices = list(range(h))
    indices.sort(key=lambda i: rows[i].tobytes())

    sorted_plane = np.stack([rows[i] for i in indices])

    # Encode: sorted indices (permutation) + delta-compressed sorted plane
    perm_data = np.array(indices, dtype=np.uint16 if h > 256 else np.uint8).tobytes()
    delta = np.zeros_like(sorted_plane)
    delta[0] = sorted_plane[0]
    delta[1:] = ((sorted_plane[1:].astype(int) - sorted_plane[:-1].astype(int)) & 0xFF).astype(np.uint8)

    total = len(bc(perm_data)) + len(bc(delta.tobytes()))
    return total

# ============================================================
# RADICAL IDEA 7: FREQUENCY-SORTED PIXELS
# ============================================================

def encode_freq_sort(plane):
    """Sort pixels by value frequency. Most common values get shortest codes.
    Encode as: value histogram + sorted pixel indices.
    This is essentially Huffman without building a tree."""
    h, w = plane.shape
    flat = plane.ravel()

    # Count frequencies
    hist = np.bincount(flat, minlength=256)

    # Sort values by frequency (most common first)
    value_order = np.argsort(-hist)
    # Map each pixel to its rank in frequency order
    rank_map = np.zeros(256, dtype=np.uint8)
    for rank, val in enumerate(value_order):
        rank_map[val] = rank

    ranked = rank_map[flat]

    # Encode: histogram (256 bytes) + ranked values (compressed)
    # The ranked values cluster low ranks → compress well
    total = len(bc(hist.astype(np.uint8).tobytes())) + len(bc(ranked.tobytes()))
    return total

# ============================================================
# TEST
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

print("=" * 90)
print("  RADICAL IDEAS: Break every assumption")
print("  kind-pasteur-2026-03-25-S8")
print("=" * 90)

tests = make_tests()

print(f"\n  {'Image':<12} {'PNG':>6} {'MED':>6} {'GolRice':>8} {'FlatRgn':>8} {'PolyFit':>8} {'2D-LZ':>7} {'BitPln':>7} {'RowSort':>8} {'FreqSrt':>8}")
print("  " + "-" * 95)

for name, img in sorted(tests.items()):
    ps = png_size(img)
    med_sz = encode_med(img) + 10

    gr_sz, gr_k = golomb_rice_size(img)
    gr_sz += 10

    flat_sz, n_regions = encode_flat_regions(img)
    flat_sz += 10

    poly_sz, n_coeffs, mean_res = encode_polyfit(img, degree=3)
    poly_sz += 10

    lz_sz, n_novel, n_ref = encode_2d_lz(img)
    lz_sz += 10

    bp_sz = encode_bitplanes(img) + 10

    sort_sz = encode_sort_transform(img) + 10

    freq_sz = encode_freq_sort(img) + 10

    # Find best
    all_methods = [("PNG", ps), ("MED", med_sz), ("GolRice", gr_sz),
                   ("FlatRgn", flat_sz), ("PolyFit", poly_sz),
                   ("2D-LZ", lz_sz), ("BitPln", bp_sz),
                   ("RowSort", sort_sz), ("FreqSrt", freq_sz)]
    best_name, best_sz = min(all_methods, key=lambda x: x[1])

    print(f"  {name:<12} {ps:>6} {med_sz:>6} {gr_sz:>8} {flat_sz:>8} {poly_sz:>8} {lz_sz:>7} {bp_sz:>7} {sort_sz:>8} {freq_sz:>8}  BEST: {best_name}")

print(f"""
  ANALYSIS OF RADICAL IDEAS:

  Golomb-Rice: Custom entropy coding for Laplacian residuals.
    Should beat zlib on MED residuals because it's distribution-optimal.

  Flat Regions: Flood-fill segmentation → mean + tiny residuals.
    Great for smooth images, terrible for noisy ones.

  Polynomial Fit: 2D polynomial surface + residual.
    Captures gradients in ~20 coefficients. Good for smooth, bad for texture.

  2D-LZ: Block-copy matching in 2D.
    Finds repeated 8x8 blocks. Good for screenshots/tiled content.

  Bit-Planes: MSB planes are smooth, LSB planes are noise.
    Separates signal from noise without wavelets.

  Row-Sort: Reorder rows by similarity, then delta-compress.
    Groups similar rows → smaller deltas. Good for horizontal patterns.

  Freq-Sort: Map pixel values to frequency ranks.
    Makes common values cluster at low ranks → better compression.
""")
