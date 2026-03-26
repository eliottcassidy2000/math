#!/usr/bin/env python3
"""
REVERSE COMPRESSION: Start with the answer, find the question.

Standard compression: image → predict → residual → encode
Reverse compression:  image → SUBTRACT STRUCTURE → noise → encode structure + noise

The insight: an image = structure + noise.
Structure is COMPACTLY describable (few parameters).
Noise is INCOMPRESSIBLE (but has known statistics).

If we can perfectly separate structure from noise, we get:
  compressed = describe(structure) + raw(noise)
  size = |structure description| + |noise|

For natural images:
  structure ≈ piecewise smooth surface (described by edges + smooth values)
  noise ≈ i.i.d. Gaussian with variance σ² per pixel

The OPTIMAL separation minimizes |structure description| + |noise|.
This is the Minimum Description Length (MDL) principle.

APPROACH: PROGRESSIVE STRUCTURE REMOVAL
  1. Remove the BIGGEST structure first (global mean/gradient) — costs 2-6 bytes
  2. Remove the NEXT biggest (block means at 16×16) — costs ~64 bytes
  3. Remove the NEXT (block means at 4×4) — costs ~1024 bytes
  4. What remains is LOCAL noise — incompressible but minimal

Each removal step is a DELETION that inserts information into the compressed stream.
The noise at each scale has decreasing energy → decreasing entropy.

THIS IS THE WAVELET IDEA but without wavelets.
We use PLAIN SUBTRACTION at each scale.
The innovation: choose the scale decomposition ADAPTIVELY based on image content.

kind-pasteur-2026-03-25-S12
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

# ============================================================
# THE REVERSE APPROACH: Subtract structure, encode what's left
# ============================================================

def subtract_global_mean(plane):
    """Remove global mean. Structure cost: 1 byte. Residual: centered data."""
    mean_val = int(round(plane.mean()))
    residual = ((plane.astype(int) - mean_val) & 0xFF).astype(np.uint8)
    return bytes([mean_val]), residual

def subtract_row_means(plane):
    """Remove per-row means. Structure cost: h bytes."""
    h, w = plane.shape
    means = np.array([int(round(plane[r].mean())) for r in range(h)], dtype=np.uint8)
    residual = np.zeros_like(plane)
    for r in range(h):
        residual[r] = ((plane[r].astype(int) - int(means[r])) & 0xFF).astype(np.uint8)
    return bytes(means), residual

def subtract_block_means(plane, bs):
    """Remove per-block means at given block size."""
    h, w = plane.shape
    bh, bw = (h + bs - 1) // bs, (w + bs - 1) // bs
    means = np.zeros((bh, bw), dtype=np.uint8)
    residual = plane.copy()

    for by in range(bh):
        for bx in range(bw):
            r0, c0 = by*bs, bx*bs
            r1, c1 = min(r0+bs, h), min(c0+bs, w)
            block = plane[r0:r1, c0:c1]
            m = int(round(block.mean()))
            means[by, bx] = m
            residual[r0:r1, c0:c1] = ((block.astype(int) - m) & 0xFF).astype(np.uint8)

    return means.tobytes(), residual

def subtract_bilinear(plane, grid_size=16):
    """Fit a bilinear surface on a coarse grid, subtract it.
    Structure: grid values. Residual: local deviation.
    This is like encoding the DC component at multiple scales."""
    h, w = plane.shape
    gh = h // grid_size + 1
    gw = w // grid_size + 1

    # Sample grid points
    grid = np.zeros((gh, gw), dtype=np.uint8)
    for gy in range(gh):
        for gx in range(gw):
            r = min(gy * grid_size, h - 1)
            c = min(gx * grid_size, w - 1)
            # Average in neighborhood
            r0, r1 = max(0, r-2), min(h, r+3)
            c0, c1 = max(0, c-2), min(w, c+3)
            grid[gy, gx] = int(round(plane[r0:r1, c0:c1].mean()))

    # Upsample grid with bilinear interpolation
    surface = np.array(Image.fromarray(grid).resize((w, h), Image.BILINEAR))
    residual = ((plane.astype(int) - surface.astype(int)) & 0xFF).astype(np.uint8)

    return grid.tobytes(), residual

# ============================================================
# PROGRESSIVE STRUCTURE REMOVAL (the main idea)
# ============================================================

def encode_progressive_removal(plane):
    """Remove structure at decreasing scales:
    1. Global mean (1 byte)
    2. 32×32 block means (~16 bytes)
    3. 8×8 block means (~256 bytes)
    4. Remaining noise (per-pixel residual)
    Each step removes structure, reducing the residual's entropy."""
    h, w = plane.shape
    total_structure = bytearray()

    # Level 0: global mean
    current = plane.copy()
    mean_val = int(round(current.mean()))
    total_structure.append(mean_val)
    current = ((current.astype(int) - mean_val) & 0xFF).astype(np.uint8)

    # Level 1: 32×32 block means of residual
    bs = 32
    bh, bw = (h+bs-1)//bs, (w+bs-1)//bs
    means1 = np.zeros((bh, bw), dtype=np.uint8)
    for by in range(bh):
        for bx in range(bw):
            r0, c0 = by*bs, bx*bs
            r1, c1 = min(r0+bs, h), min(c0+bs, w)
            block = current[r0:r1, c0:c1]
            # Convert to signed for proper mean
            signed = block.astype(np.int16)
            signed[signed > 128] -= 256
            m = int(round(signed.mean()))
            means1[by, bx] = m & 0xFF
            current[r0:r1, c0:c1] = ((block.astype(int) - m) & 0xFF).astype(np.uint8)
    total_structure.extend(means1.tobytes())

    # Level 2: 8×8 block means of residual
    bs = 8
    bh2, bw2 = (h+bs-1)//bs, (w+bs-1)//bs
    means2 = np.zeros((bh2, bw2), dtype=np.uint8)
    for by in range(bh2):
        for bx in range(bw2):
            r0, c0 = by*bs, bx*bs
            r1, c1 = min(r0+bs, h), min(c0+bs, w)
            block = current[r0:r1, c0:c1]
            signed = block.astype(np.int16)
            signed[signed > 128] -= 256
            m = int(round(signed.mean()))
            means2[by, bx] = m & 0xFF
            current[r0:r1, c0:c1] = ((block.astype(int) - m) & 0xFF).astype(np.uint8)
    total_structure.extend(means2.tobytes())

    # Level 3: remaining per-pixel noise
    structure_compressed = bc(bytes(total_structure))
    noise_compressed = bc(current.tobytes())

    return len(structure_compressed) + len(noise_compressed)

# ============================================================
# THE DUAL VIEW: encode the DECODER (generative model)
# ============================================================

def encode_as_generator(plane):
    """Instead of encoding the image, encode a PROGRAM that generates it.
    The program: "start with bilinear surface, add corrections."

    Level 0: bilinear surface on 16×16 grid (cheap, captures global shape)
    Level 1: corrections at 4×4 grid (medium, captures edges)
    Level 2: per-pixel corrections (expensive, captures noise)

    Each level refines the previous. Encode each level independently."""
    h, w = plane.shape

    # Level 0: coarse grid (16×16 spacing)
    g0_h, g0_w = h//16+1, w//16+1
    g0 = np.zeros((g0_h, g0_w), dtype=np.uint8)
    for gy in range(g0_h):
        for gx in range(g0_w):
            r, c = min(gy*16, h-1), min(gx*16, w-1)
            g0[gy, gx] = plane[r, c]

    # Upsample to full resolution
    surface0 = np.array(Image.fromarray(g0).resize((w, h), Image.BILINEAR))

    # Level 1: corrections at 4×4 grid
    g1_h, g1_w = h//4+1, w//4+1
    resid0 = ((plane.astype(int) - surface0.astype(int)) & 0xFF).astype(np.uint8)
    g1 = np.zeros((g1_h, g1_w), dtype=np.uint8)
    for gy in range(g1_h):
        for gx in range(g1_w):
            r, c = min(gy*4, h-1), min(gx*4, w-1)
            g1[gy, gx] = resid0[r, c]

    surface1 = np.array(Image.fromarray(g1).resize((w, h), Image.BILINEAR))
    resid1 = ((resid0.astype(int) - surface1.astype(int)) & 0xFF).astype(np.uint8)

    # Level 2: per-pixel corrections
    total = (len(bc(g0.tobytes())) +
             len(bc(g1.tobytes())) +
             len(bc(resid1.tobytes())))

    return total

# ============================================================
# THE INVERSION: what if the IMAGE is the compressed form?
# ============================================================

def encode_self_describing(plane):
    """The image IS its own most compact description in a different basis.

    KEY INSIGHT: the autocorrelation matrix R of the image defines
    an eigenbasis (the KLT/PCA basis). In this basis, the image is
    a set of INDEPENDENT coefficients with DECREASING variance.

    We don't compute the full KLT (too expensive). Instead, we use
    the ROWS of the image as vectors and compute their PCA.
    Store: mean row + principal components + coefficients.

    For h rows of width w:
      - Mean row: w bytes
      - k principal components: k×w bytes (but they're orthogonal → compress well)
      - Coefficients: h×k bytes (independent → minimal entropy)
      - Residual: h×w - k components worth of noise

    This is literally solving the compression problem in reverse:
    find the basis that makes the data simplest."""
    h, w = plane.shape
    rows = plane.astype(float)

    # Compute mean row
    mean_row = rows.mean(axis=0)
    centered = rows - mean_row

    # PCA: compute covariance, find top-k eigenvectors
    # For speed, use SVD on centered rows
    # centered = U @ S @ Vt, where Vt rows are principal components
    try:
        U, S, Vt = np.linalg.svd(centered, full_matrices=False)
    except:
        return len(bc(plane.tobytes()))  # fallback

    # How many components to keep? Use MDL criterion:
    # Keep k components if storing k components + residual is smaller
    total_var = np.sum(S**2)
    best_k = 0
    best_size = len(bc(plane.tobytes()))  # baseline: raw

    for k in [1, 2, 4, 8, 16, 32, min(64, min(h, w))]:
        if k > len(S): break

        # Reconstruct from k components
        recon = U[:, :k] @ np.diag(S[:k]) @ Vt[:k, :]
        recon_full = recon + mean_row
        recon_uint8 = np.clip(np.round(recon_full), 0, 255).astype(np.uint8)

        residual = ((plane.astype(int) - recon_uint8.astype(int)) & 0xFF).astype(np.uint8)

        # Size: mean_row + k components (Vt[:k]) + k coefficients per row (U[:,:k]*S[:k]) + residual
        mean_bytes = bc(np.round(mean_row).astype(np.uint8).tobytes())
        # Components as int8 (scaled)
        vt_scaled = np.clip(np.round(Vt[:k] * 127 / (np.abs(Vt[:k]).max() + 1e-10)), -128, 127).astype(np.int8)
        comp_bytes = bc(vt_scaled.tobytes())
        # Coefficients: U[:,:k] * S[:k], quantized to int16
        coeffs = U[:, :k] * S[:k]
        coeffs_scaled = np.clip(np.round(coeffs), -32768, 32767).astype(np.int16)
        coeff_bytes = bc(coeffs_scaled.tobytes())
        # Scale factors
        scale_bytes = bc(np.array(S[:k], dtype=np.float32).tobytes())

        resid_bytes = bc(residual.tobytes())

        total = len(mean_bytes) + len(comp_bytes) + len(coeff_bytes) + len(resid_bytes) + len(scale_bytes)

        if total < best_size:
            best_size = total
            best_k = k

    return best_size

# ============================================================
# THE DELETION APPROACH: remove pixels until remainder is simple
# ============================================================

def encode_by_deletion(plane):
    """Instead of predicting missing pixels, DELETE redundant pixels.

    Start with the full image. Repeatedly find and remove the pixel
    whose value is MOST PREDICTABLE from its neighbors.
    Store only the UNPREDICTABLE pixels + their positions.

    The predictable pixels can be reconstructed by the decoder
    using the same prediction rule.

    This is like building a MINIMAL SPANNING SET: the smallest
    subset of pixels from which the full image can be reconstructed."""
    h, w = plane.shape

    # Simpler version: flag each pixel as "predictable" or "novel"
    # A pixel is predictable if |actual - MED(neighbors)| <= threshold
    # Store: threshold + flag bitmap + novel pixel values
    img = plane.astype(int)
    best_size = len(bc(plane.tobytes()))  # baseline
    best_threshold = -1

    for threshold in [0]:  # ONLY threshold=0 is lossless
        flags = np.zeros(h * w, dtype=np.uint8)
        novels = bytearray()

        for r in range(h):
            for c in range(w):
                a = img[r, c-1] if c > 0 else 128
                b = img[r-1, c] if r > 0 else 128
                c2 = img[r-1, c-1] if r > 0 and c > 0 else 128

                mn, mx = min(a, b), max(a, b)
                if c2 >= mx: pred = mn
                elif c2 <= mn: pred = mx
                else: pred = a + b - c2

                if abs(img[r, c] - pred) <= threshold:
                    flags[r * w + c] = 0  # predictable
                else:
                    flags[r * w + c] = 1  # novel
                    novels.append(img[r, c] & 0xFF)

        packed_flags = np.packbits(flags)
        total = 1 + len(bc(bytes(packed_flags))) + len(bc(bytes(novels)))

        if total < best_size:
            best_size = total
            best_threshold = threshold

    return best_size

# ============================================================
# BASELINES
# ============================================================

def m_med(plane):
    h,w=plane.shape; res=bytearray()
    for r in range(h):
        for c in range(w):
            a=int(plane[r,c-1]) if c>0 else 0; b=int(plane[r-1,c]) if r>0 else 0
            c3=int(plane[r-1,c-1]) if r>0 and c>0 else 0
            res.append((int(plane[r,c])-med(a,b,c3))&0xFF)
    return len(bc(bytes(res)))

def lpc_coeffs(sig, order):
    n=len(sig)
    if n<order+1: return np.zeros(order)
    r2=np.correlate(sig,sig,mode='full'); r2=r2[n-1:n+order]
    if r2[0]<1e-10: return np.zeros(order)
    a=np.zeros(order); e=r2[0]
    for i in range(order):
        if e<1e-10: break
        lam=r2[i+1]
        for j in range(i): lam-=a[j]*r2[i-j]
        k=lam/e; a_new=np.zeros(order); a_new[i]=k
        for j in range(i): a_new[j]=a[j]-k*a[i-1-j]
        a=a_new; e*=(1-k*k)
    return a

def m_lpc8(plane):
    sig=plane.ravel().astype(float); n=len(sig)
    coeffs=lpc_coeffs(sig-sig.mean(),8)
    res=np.zeros(n,dtype=np.uint8)
    for i in range(n):
        pred=sum(coeffs[j]*sig[i-1-j] for j in range(min(i,8)))
        res[i]=(int(sig[i])-int(np.clip(round(pred),0,255)))&0xFF
    return len(np.array(coeffs,dtype=np.float32).tobytes())+len(bc(bytes(res)))

# ============================================================
# TEST
# ============================================================

SZ = 128

def make_tests():
    T={}; np.random.seed(42)
    x,y=np.meshgrid(np.arange(SZ,dtype=float),np.arange(SZ,dtype=float))
    r=np.sqrt((x-SZ/2)**2+(y-SZ/2)**2)
    T["solid"]=np.full((SZ,SZ),128,dtype=np.uint8)
    T["grad_h"]=np.tile(np.linspace(0,255,SZ,dtype=np.uint8),(SZ,1))
    T["grad_d"]=((x+y)*255//(2*SZ-2)).astype(np.uint8)
    T["circles"]=(np.sin(r/5)*127+128).astype(np.uint8)
    T["blob"]=(np.exp(-r**2/(2*(SZ/4)**2))*255).astype(np.uint8)
    T["random"]=np.random.randint(0,256,(SZ,SZ),dtype=np.uint8)
    sm=np.random.randint(0,256,(SZ//16,SZ//16),dtype=np.uint8)
    T["natural"]=np.clip(np.array(Image.fromarray(sm).resize((SZ,SZ),Image.BILINEAR)).astype(float)
                         +np.random.normal(0,10,(SZ,SZ)),0,255).astype(np.uint8)
    T["stripes45"]=(((x+y)/8).astype(int)%2*255).astype(np.uint8)
    return T

print("="*80)
print("  REVERSE COMPRESSION: Start with the answer, find the question")
print("  kind-pasteur-2026-03-25-S12")
print("="*80)

tests=make_tests()

METHODS=[
    ("MED",           m_med),
    ("LPC8",          m_lpc8),
    ("ProgRemoval",   encode_progressive_removal),
    ("Generator",     encode_as_generator),
    ("PCA-rows",      encode_self_describing),
    ("Deletion",      encode_by_deletion),
    ("Raw",           lambda p: len(bc(p.tobytes()))),
]

print(f"\n  {'Image':<12} {'PNG':>6}", end="")
for mname, _ in METHODS:
    print(f" {mname:>10}", end="")
print(f"  {'BEST':>12}")
print("  "+"-"*(20+11*len(METHODS)+16))

for name,img in sorted(tests.items()):
    ps=png_size(img); sizes={}
    for mname,mfunc in METHODS:
        try: sizes[mname]=mfunc(img)+10
        except: sizes[mname]=999999
    best=min(sizes,key=sizes.get)
    ratio=ps/sizes[best] if sizes[best]>0 and sizes[best]<999999 else 0
    print(f"  {name:<12} {ps:>6}", end="")
    for mname,_ in METHODS:
        v=sizes[mname]; marker="*" if mname==best else " "
        print(f" {v if v<999999 else 'ERR':>9}{marker}", end="")
    print(f"  {best:>12} {ratio:.2f}x")

print(f"""
  THE REVERSE INSIGHTS:

  1. PROGRESSIVE REMOVAL: subtract structure at decreasing scales.
     Like wavelets but with plain block means.
     Should excel at multi-scale images (natural photos).

  2. GENERATOR MODEL: encode the image as a multi-scale generator.
     Coarse grid → upsample → add corrections → repeat.
     The decoder IS the generative model.

  3. PCA ON ROWS: find the eigenbasis of the row-space.
     Represent each row as coefficients in the eigenbasis.
     For images with repeated row patterns: massive compression.

  4. DELETION: remove predictable pixels, encode only surprises.
     The decoder reconstructs predictable pixels from neighbors.
     For smooth images: most pixels are predictable → few surprises.
""")
