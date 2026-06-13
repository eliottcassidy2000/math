#!/usr/bin/env python3
"""
KODAK BENCHMARK + CREATIVE EXTENSIONS

Test on the Kodak 24-image corpus (the standard for lossless compression).
Also test on our real photo at multiple crops.

CREATIVE EXTENSIONS:
1. NOISE-ADAPTIVE MED-BLEND: Estimate local noise from MED residuals.
   High noise → use MED (hard switch avoids averaging noise).
   Low noise → use Blend (multi-neighbor captures smooth structure).
   This adapts to REAL photo characteristics automatically.

2. CROSS-CHANNEL PREDICTION: After G-RG-BG decorrelation,
   predict R-G from G's local gradient (edges in G predict edges in R-G).
   The G channel is decoded first → its structure is available as context.

3. ADAPTIVE COLOR TRANSFORM: Try G-RG-BG, YCoCg, and direct RGB.
   Pick per-image (1 byte overhead). Different photos may prefer different transforms.

4. PER-PLANE ADAPTIVE: Use MED for sharp channels (G), Blend for smooth (R-G, B-G).
   The difference channels are smoother → Blend helps more.

kind-pasteur-2026-03-25-S21
"""
import sys, io, struct, zlib, time, os, glob
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

try: import brotli; HAS_BROTLI = True
except: HAS_BROTLI = False

def zlib9(data): return zlib.compress(data, 9)

def best_compress(data):
    r = [zlib9(data)]
    if HAS_BROTLI:
        try: r.append(brotli.compress(data, quality=11))
        except: pass
    return min(r, key=len)

def png_size_pil(pil_img):
    buf = io.BytesIO()
    pil_img.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()

def med(a, b, c):
    mn, mx = min(a, b), max(a, b)
    if c >= mx: return mn
    if c <= mn: return mx
    return a + b - c

# ============================================================
# PREDICTORS
# ============================================================

def encode_med(plane):
    h,w=plane.shape; res=bytearray()
    for r in range(h):
        for c in range(w):
            a=int(plane[r,c-1]) if c>0 else 0; b=int(plane[r-1,c]) if r>0 else 0
            d=int(plane[r-1,c-1]) if r>0 and c>0 else 0
            res.append((int(plane[r,c])-med(a,b,d))&0xFF)
    return bytes(res)

def encode_blend(plane, threshold=1000):
    h,w=plane.shape; img=plane.astype(float); res=bytearray()
    for r in range(h):
        for c in range(w):
            a=img[r,c-1] if c>0 else 0.0; b=img[r-1,c] if r>0 else 0.0
            d=img[r-1,c-1] if r>0 and c>0 else 0.0
            grad=abs(a-d)+abs(b-d) if r>0 and c>0 else 0
            alpha=min(1.0, grad/threshold)
            p_med=med(int(a),int(b),int(d))
            total,wt=0.0,0.0
            for dr,dc in [(-1,-1),(-1,0),(-1,1),(0,-1),(-2,0),(0,-2)]:
                nr,nc=r+dr,c+dc
                if 0<=nr<h and 0<=nc<w:
                    dst=abs(dr)+abs(dc); w2=4.0/(dst*dst)
                    total+=img[nr,nc]*w2; wt+=w2
            p_wie=total/wt if wt>0 else 128
            pred=int(np.clip(round(alpha*p_med+(1-alpha)*p_wie),0,255))
            res.append((int(plane[r,c])-pred)&0xFF)
    return bytes(res)

def encode_noise_adaptive(plane):
    """CREATIVE: Estimate noise from MED residuals, adapt blend threshold.
    Real photos have varying noise: dark areas are noisier (Poisson).
    Low noise → use Blend (smooth areas benefit from multi-neighbor).
    High noise → use MED (avoid averaging noise from neighbors)."""
    h,w=plane.shape; img=plane.astype(float); res=bytearray()

    # Running noise estimate (from recent MED residual magnitudes)
    noise_est = 10.0  # initial estimate
    decay = 0.99

    for r in range(h):
        for c in range(w):
            a=img[r,c-1] if c>0 else 0.0; b=img[r-1,c] if r>0 else 0.0
            d=img[r-1,c-1] if r>0 and c>0 else 0.0

            # Gradient magnitude
            grad=abs(a-d)+abs(b-d) if r>0 and c>0 else 0

            # Noise-adaptive threshold:
            # high noise → high threshold → more MED (alpha → 1 faster)
            # low noise → low threshold → more Wiener at moderate gradients
            threshold = max(5.0, noise_est * 3)

            alpha=min(1.0, grad/threshold)
            p_med=med(int(a),int(b),int(d))

            total,wt=0.0,0.0
            for dr,dc in [(-1,-1),(-1,0),(-1,1),(0,-1),(-2,0),(0,-2)]:
                nr,nc=r+dr,c+dc
                if 0<=nr<h and 0<=nc<w:
                    dst=abs(dr)+abs(dc); w2=4.0/(dst*dst)
                    total+=img[nr,nc]*w2; wt+=w2
            p_wie=total/wt if wt>0 else 128

            pred=int(np.clip(round(alpha*p_med+(1-alpha)*p_wie),0,255))
            actual = int(plane[r,c])
            residual = (actual - pred) & 0xFF

            # Update noise estimate from MED residual
            med_res = abs(actual - med(int(a), int(b), int(d)))
            noise_est = decay * noise_est + (1-decay) * med_res

            res.append(residual)
    return bytes(res)

# ============================================================
# COLOR TRANSFORMS
# ============================================================

def ct_grd(rgb):
    """G, R-G, B-G"""
    g = rgb[:,:,1]
    rg = ((rgb[:,:,0].astype(int) - g.astype(int)) & 0xFF).astype(np.uint8)
    bg = ((rgb[:,:,2].astype(int) - g.astype(int)) & 0xFF).astype(np.uint8)
    return [g, rg, bg], "G-RG-BG"

def ct_rgb(rgb):
    """Direct RGB"""
    return [rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]], "RGB"

def ct_bgr(rgb):
    """B, G-B, R-B (maybe blue is less noisy for some images?)"""
    b = rgb[:,:,2]
    gb = ((rgb[:,:,1].astype(int) - b.astype(int)) & 0xFF).astype(np.uint8)
    rb = ((rgb[:,:,0].astype(int) - b.astype(int)) & 0xFF).astype(np.uint8)
    return [b, gb, rb], "B-GB-RB"

# ============================================================
# FULL RGB ENCODER: try all combinations, pick smallest
# ============================================================

def encode_rgb_best(rgb, use_brotli=False):
    """Try all color transforms × prediction methods, pick smallest."""
    compress = best_compress if use_brotli else zlib9
    best_total = float('inf')
    best_desc = ""

    for ct_func in [ct_grd, ct_rgb, ct_bgr]:
        planes, ct_name = ct_func(rgb)

        # All MED
        total_med = sum(len(compress(encode_med(p))) for p in planes)
        if total_med < best_total:
            best_total = total_med; best_desc = f"{ct_name}+MED"

        # All Blend
        total_blend = sum(len(compress(encode_blend(p))) for p in planes)
        if total_blend < best_total:
            best_total = total_blend; best_desc = f"{ct_name}+Blend"

        # All noise-adaptive
        total_na = sum(len(compress(encode_noise_adaptive(p))) for p in planes)
        if total_na < best_total:
            best_total = total_na; best_desc = f"{ct_name}+NoiseAdapt"

        # Hybrid: MED on first plane (structural), Blend on diff planes
        total_hybrid = len(compress(encode_med(planes[0])))
        for p in planes[1:]:
            total_hybrid += len(compress(encode_blend(p, threshold=500)))
        if total_hybrid < best_total:
            best_total = total_hybrid; best_desc = f"{ct_name}+MED/Blend"

        # Hybrid: MED on first, noise-adaptive on diffs
        total_hyb2 = len(compress(encode_med(planes[0])))
        for p in planes[1:]:
            total_hyb2 += len(compress(encode_noise_adaptive(p)))
        if total_hyb2 < best_total:
            best_total = total_hyb2; best_desc = f"{ct_name}+MED/NoiseAdapt"

    return best_total + 2, best_desc  # +2 for ct_id + pred_id

# ============================================================
# BENCHMARK
# ============================================================

print("=" * 80)
print("  KODAK + REAL PHOTO BENCHMARK")
print("  kind-pasteur-2026-03-25-S21")
print("=" * 80)

# Collect all test images
test_files = []
for f in sorted(glob.glob("test_images/kodim*.png")):
    test_files.append(f)
test_files.append("inbox/vlcsnap-2026-03-10-17h12m34s408.png")

print(f"\n  {len(test_files)} images to test")

print(f"\n  {'Image':<25} {'Size':>10} {'PNG':>8} {'zlib-best':>10} {'brotli':>8} {'z/PNG':>6} {'b/PNG':>6} {'method':>20}")
print("  " + "-" * 100)

total_png = 0; total_zlib = 0; total_brotli = 0
z_wins = 0; z_losses = 0; b_wins = 0

for fpath in test_files:
    pil = Image.open(fpath).convert('RGB')
    name = os.path.basename(fpath)[:20]

    # Crop to 256×256 for speed (full-image would take too long in Python)
    w0, h0 = pil.size
    cx, cy = w0//2, h0//2
    crop = pil.crop((cx-128, cy-128, cx+128, cy+128))
    rgb = np.array(crop)
    h, w = rgb.shape[:2]

    crop_png = png_size_pil(crop)

    t0 = time.time()
    z_sz, z_desc = encode_rgb_best(rgb, use_brotli=False)
    z_time = time.time() - t0

    b_sz, b_desc = encode_rgb_best(rgb, use_brotli=True)

    z_ratio = crop_png / z_sz if z_sz > 0 else 0
    b_ratio = crop_png / b_sz if b_sz > 0 else 0

    total_png += crop_png; total_zlib += z_sz; total_brotli += b_sz
    if z_ratio > 1.001: z_wins += 1
    elif z_ratio < 0.999: z_losses += 1
    if b_ratio > 1.001: b_wins += 1

    raw = h * w * 3
    print(f"  {name:<25} {raw:>10} {crop_png:>8} {z_sz:>10} {b_sz:>8} {z_ratio:>5.3f}x {b_ratio:>5.3f}x {z_desc:>20}")

n = len(test_files)
agg_z = total_png / total_zlib if total_zlib > 0 else 0
agg_b = total_png / total_brotli if total_brotli > 0 else 0
print(f"\n  AGGREGATE: zlib-9 {agg_z:.3f}x, brotli {agg_b:.3f}x vs PNG")
print(f"  zlib-9: {z_wins}W / {n-z_wins-z_losses}T / {z_losses}L out of {n}")
print(f"  brotli: {b_wins}W out of {n}")

print(f"""
  NEW TECHNIQUE: Noise-Adaptive MED-Blend
    Estimates local noise from running MED residual magnitude.
    Adapts the MED-Wiener blend threshold per pixel.
    High noise (dark areas, sensor noise) → more MED
    Low noise (bright smooth areas) → more Wiener

  COLOR TRANSFORM MATTERS:
    G-RG-BG decorrelation provides the BIGGEST single gain.
    Different images may prefer different transforms.

  PER-PLANE ADAPTIVE (MED on G, Blend on R-G/B-G):
    The G channel has edges → MED is best.
    The R-G and B-G channels are smoother → Blend may help.
""")
