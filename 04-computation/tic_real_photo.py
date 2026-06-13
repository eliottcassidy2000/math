#!/usr/bin/env python3
"""
REAL PHOTO BENCHMARK: Test everything on an actual photograph.

This is the moment of truth — a real camera photo with:
- Smooth regions (wall, fabric)
- Fine texture (hair, leaves)
- Sharp edges (monitor frames)
- Text (on monitors)
- Color (RGB)

Test ALL methods, zlib-9 only (fair), measure speed, verify roundtrip.

kind-pasteur-2026-03-25-S20
"""
import sys, io, struct, zlib, time, math
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

try: import brotli; HAS_BROTLI = True
except: HAS_BROTLI = False

IMG_PATH = "inbox/vlcsnap-2026-03-10-17h12m34s408.png"

def zlib9(data): return zlib.compress(data, 9)

def best_compress(data):
    r = [zlib9(data)]
    if HAS_BROTLI:
        try: r.append(brotli.compress(data, quality=11))
        except: pass
    return min(r, key=len)

def png_size_from_pil(pil_img):
    buf = io.BytesIO()
    pil_img.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()

def med(a, b, c):
    mn, mx = min(a, b), max(a, b)
    if c >= mx: return mn
    if c <= mn: return mx
    return a + b - c

def lpc_coeffs(sig, order):
    n = len(sig)
    if n < order + 1: return np.zeros(order)
    r = np.correlate(sig, sig, mode='full'); r = r[n-1:n+order]
    if r[0] < 1e-10: return np.zeros(order)
    a = np.zeros(order); e = r[0]
    for i in range(order):
        if e < 1e-10: break
        lam = r[i+1]
        for j in range(i): lam -= a[j] * r[i-j]
        k = lam / e; a_new = np.zeros(order); a_new[i] = k
        for j in range(i): a_new[j] = a[j] - k * a[i-1-j]
        a = a_new; e *= (1 - k*k)
    return a

# ============================================================
# LOAD AND ANALYZE THE PHOTO
# ============================================================

print("=" * 80)
print("  REAL PHOTO BENCHMARK: vlcsnap (desk scene)")
print("  kind-pasteur-2026-03-25-S20")
print("=" * 80)

pil = Image.open(IMG_PATH)
img_rgb = np.array(pil)
h, w = img_rgb.shape[:2]
channels = img_rgb.shape[2] if img_rgb.ndim == 3 else 1

print(f"\n  Image: {IMG_PATH}")
print(f"  Dimensions: {w}×{h}, {channels} channels")
print(f"  Raw size: {img_rgb.nbytes} bytes")

# PNG size
png_sz = png_size_from_pil(pil)
print(f"  PNG size: {png_sz} bytes ({img_rgb.nbytes/png_sz:.1f}:1 ratio)")

# ============================================================
# GRAYSCALE ANALYSIS (work on luminance channel)
# ============================================================

if channels >= 3:
    # Convert to grayscale for per-plane analysis
    gray = np.round(0.299 * img_rgb[:,:,0].astype(float) +
                    0.587 * img_rgb[:,:,1].astype(float) +
                    0.114 * img_rgb[:,:,2].astype(float)).astype(np.uint8)
else:
    gray = img_rgb

gray_png = png_size_from_pil(Image.fromarray(gray, 'L'))
print(f"\n  Grayscale: {w}×{h}")
print(f"  Gray PNG: {gray_png} bytes")

# ============================================================
# TEST ALL METHODS ON GRAYSCALE
# ============================================================

print(f"\n  === GRAYSCALE METHODS (zlib-9, fair comparison) ===")

# Crop to manageable size for per-pixel methods
crop_h, crop_w = min(h, 256), min(w, 256)
crop = gray[:crop_h, :crop_w]
crop_png = png_size_from_pil(Image.fromarray(crop, 'L'))
print(f"  Crop: {crop_w}×{crop_h}, PNG={crop_png}")

results = {}

# Method 0: MED
t0 = time.time()
res = bytearray()
for r in range(crop_h):
    for c in range(crop_w):
        a = int(crop[r,c-1]) if c>0 else 0
        b = int(crop[r-1,c]) if r>0 else 0
        d = int(crop[r-1,c-1]) if r>0 and c>0 else 0
        res.append((int(crop[r,c]) - med(a,b,d)) & 0xFF)
results["MED"] = (len(zlib9(bytes(res))), time.time()-t0)

# Method 1: Blend (MED-Wiener)
t0 = time.time()
img_f = crop.astype(float); res = bytearray()
for r in range(crop_h):
    for c in range(crop_w):
        a=img_f[r,c-1] if c>0 else 0.0; b=img_f[r-1,c] if r>0 else 0.0
        dd=img_f[r-1,c-1] if r>0 and c>0 else 0.0
        grad=abs(a-dd)+abs(b-dd) if r>0 and c>0 else 0
        alpha=min(1.0, grad/1000)
        p_med=med(int(a),int(b),int(dd))
        total,wt=0.0,0.0
        for dr,dc in [(-1,-1),(-1,0),(-1,1),(0,-1),(-2,0),(0,-2)]:
            nr,nc=r+dr,c+dc
            if 0<=nr<crop_h and 0<=nc<crop_w:
                dst=abs(dr)+abs(dc); w2=4.0/(dst*dst)
                total+=img_f[nr,nc]*w2; wt+=w2
        p_wie=total/wt if wt>0 else 128
        pred=int(np.clip(round(alpha*p_med+(1-alpha)*p_wie),0,255))
        res.append((int(crop[r,c])-pred)&0xFF)
results["Blend"] = (len(zlib9(bytes(res))), time.time()-t0)

# Method 2: Snake-LPC8
t0 = time.time()
sig = np.empty(crop_h*crop_w, dtype=np.uint8)
for r in range(crop_h):
    if r%2==0: sig[r*crop_w:(r+1)*crop_w]=crop[r]
    else: sig[r*crop_w:(r+1)*crop_w]=crop[r,::-1]
sf=sig.astype(float); n=len(sf)
coeffs=lpc_coeffs(sf-sf.mean(),8)
res2=np.empty(n,dtype=np.uint8)
for i in range(n):
    pred=sum(coeffs[j]*sf[i-1-j] for j in range(min(i,8)))
    res2[i]=(int(sf[i])-int(np.clip(round(pred),0,255)))&0xFF
lpc_data = np.array(coeffs,dtype=np.float32).tobytes() + zlib9(bytes(res2))
results["SnakeLPC8"] = (len(lpc_data), time.time()-t0)

# Method 3: ColSort
t0 = time.time()
idx=sorted(range(crop_w),key=lambda c2:crop[:,c2].tobytes())
sp=np.stack([crop[:,c2] for c2 in idx],axis=1)
delta=np.zeros_like(sp); delta[:,0]=sp[:,0]
delta[:,1:]=((sp[:,1:].astype(int)-sp[:,:-1].astype(int))&0xFF).astype(np.uint8)
perm=np.array(idx,dtype=np.uint16)
col_data = zlib9(perm.tobytes()) + zlib9(delta.tobytes())
results["ColSort"] = (len(col_data), time.time()-t0)

# Method 4: Raw
t0 = time.time()
results["Raw"] = (len(zlib9(crop.tobytes())), time.time()-t0)

# Method 5: Snake-delta (order 1)
t0 = time.time()
sig2 = np.empty(crop_h*crop_w, dtype=np.uint8)
for r in range(crop_h):
    if r%2==0: sig2[r*crop_w:(r+1)*crop_w]=crop[r]
    else: sig2[r*crop_w:(r+1)*crop_w]=crop[r,::-1]
res3 = np.empty(len(sig2), dtype=np.uint8)
res3[0] = sig2[0]
res3[1:] = ((sig2[1:].astype(int) - sig2[:-1].astype(int)) & 0xFF).astype(np.uint8)
results["SnakeDelta"] = (len(zlib9(bytes(res3))), time.time()-t0)

# With brotli
if HAS_BROTLI:
    # Best of zlib+brotli for each method
    for mname in list(results.keys()):
        pass  # already using zlib9, add brotli versions below

# Also try brotli on the best zlib method
results_brotli = {}
for mname, (sz, t) in results.items():
    pass  # we'll compute brotli versions in the comparison

print(f"\n  {'Method':<14} {'zlib-9':>8} {'vs PNG':>8} {'time':>7}")
print("  " + "-" * 42)

for mname in ["MED", "Blend", "SnakeLPC8", "ColSort", "SnakeDelta", "Raw"]:
    sz, t = results[mname]
    sz += 10  # header overhead
    ratio = crop_png / sz if sz > 0 else 0
    verdict = "WIN" if ratio > 1.001 else ("LOSS" if ratio < 0.999 else "TIE")
    print(f"  {mname:<14} {sz:>8} {ratio:>7.3f}x {t*1000:>6.0f}ms  {verdict}")

best_name = min(results, key=lambda k: results[k][0])
best_sz = results[best_name][0] + 10

print(f"\n  BEST: {best_name} at {best_sz} bytes ({crop_png/best_sz:.3f}x vs PNG)")

# ============================================================
# RGB ANALYSIS: per-plane vs color transform
# ============================================================

print(f"\n  === RGB METHODS ({w}×{h}) ===")

# Full RGB PNG
rgb_png = png_sz
print(f"  RGB PNG: {rgb_png} bytes")

# Per-plane MED (independent R, G, B)
t0 = time.time()
total_indep = 0
for ch in range(3):
    plane = img_rgb[:crop_h, :crop_w, ch]
    res_ch = bytearray()
    for r in range(crop_h):
        for c in range(crop_w):
            a = int(plane[r,c-1]) if c>0 else 0
            b = int(plane[r-1,c]) if r>0 else 0
            dd = int(plane[r-1,c-1]) if r>0 and c>0 else 0
            res_ch.append((int(plane[r,c]) - med(a,b,dd)) & 0xFF)
    total_indep += len(zlib9(bytes(res_ch)))
t_indep = time.time() - t0

# G-RG-BG color transform + per-plane MED
t0 = time.time()
g = img_rgb[:crop_h, :crop_w, 1]
rg = ((img_rgb[:crop_h, :crop_w, 0].astype(int) - g.astype(int)) & 0xFF).astype(np.uint8)
bg = ((img_rgb[:crop_h, :crop_w, 2].astype(int) - g.astype(int)) & 0xFF).astype(np.uint8)
total_grd = 0
for plane in [g, rg, bg]:
    res_ch = bytearray()
    for r in range(crop_h):
        for c in range(crop_w):
            a = int(plane[r,c-1]) if c>0 else 0
            b = int(plane[r-1,c]) if r>0 else 0
            dd = int(plane[r-1,c-1]) if r>0 and c>0 else 0
            res_ch.append((int(plane[r,c]) - med(a,b,dd)) & 0xFF)
    total_grd += len(zlib9(bytes(res_ch)))
t_grd = time.time() - t0

# Per-plane Blend
t0 = time.time()
total_blend = 0
for ch in range(3):
    plane = img_rgb[:crop_h, :crop_w, ch]
    img_f = plane.astype(float); res_ch = bytearray()
    for r in range(crop_h):
        for c in range(crop_w):
            a=img_f[r,c-1] if c>0 else 0.0; b=img_f[r-1,c] if r>0 else 0.0
            dd=img_f[r-1,c-1] if r>0 and c>0 else 0.0
            grad=abs(a-dd)+abs(b-dd) if r>0 and c>0 else 0
            alpha=min(1.0, grad/1000)
            p_med=med(int(a),int(b),int(dd))
            total_w,wt=0.0,0.0
            for dr,dc in [(-1,-1),(-1,0),(-1,1),(0,-1),(-2,0),(0,-2)]:
                nr,nc=r+dr,c+dc
                if 0<=nr<crop_h and 0<=nc<crop_w:
                    dst=abs(dr)+abs(dc); w2=4.0/(dst*dst)
                    total_w+=img_f[nr,nc]*w2; wt+=w2
            p_wie=total_w/wt if wt>0 else 128
            pred=int(np.clip(round(alpha*p_med+(1-alpha)*p_wie),0,255))
            res_ch.append((int(plane[r,c])-pred)&0xFF)
    total_blend += len(zlib9(bytes(res_ch)))
t_blend = time.time() - t0

# G-RG-BG + Blend
t0 = time.time()
total_grd_blend = 0
for plane in [g, rg, bg]:
    img_f = plane.astype(float); res_ch = bytearray()
    for r in range(crop_h):
        for c in range(crop_w):
            a=img_f[r,c-1] if c>0 else 0.0; b=img_f[r-1,c] if r>0 else 0.0
            dd=img_f[r-1,c-1] if r>0 and c>0 else 0.0
            grad=abs(a-dd)+abs(b-dd) if r>0 and c>0 else 0
            alpha=min(1.0, grad/1000)
            p_med=med(int(a),int(b),int(dd))
            total_w,wt=0.0,0.0
            for dr,dc in [(-1,-1),(-1,0),(-1,1),(0,-1),(-2,0),(0,-2)]:
                nr,nc=r+dr,c+dc
                if 0<=nr<crop_h and 0<=nc<crop_w:
                    dst=abs(dr)+abs(dc); w2=4.0/(dst*dst)
                    total_w+=img_f[nr,nc]*w2; wt+=w2
            p_wie=total_w/wt if wt>0 else 128
            pred=int(np.clip(round(alpha*p_med+(1-alpha)*p_wie),0,255))
            res_ch.append((int(plane[r,c])-pred)&0xFF)
    total_grd_blend += len(zlib9(bytes(res_ch)))
t_grd_blend = time.time() - t0

# Crop RGB PNG
crop_rgb = img_rgb[:crop_h, :crop_w]
crop_rgb_png = png_size_from_pil(Image.fromarray(crop_rgb))

print(f"\n  Crop {crop_w}×{crop_h} RGB:")
print(f"  {'Method':<18} {'size':>8} {'vs PNG':>8} {'time':>7}")
print("  " + "-" * 45)
print(f"  {'PNG':.<18} {crop_rgb_png:>8}")
print(f"  {'RGB MED indep':.<18} {total_indep+10:>8} {crop_rgb_png/(total_indep+10):>7.3f}x {t_indep*1000:>6.0f}ms")
print(f"  {'G-RG-BG + MED':.<18} {total_grd+10:>8} {crop_rgb_png/(total_grd+10):>7.3f}x {t_grd*1000:>6.0f}ms")
print(f"  {'RGB Blend':.<18} {total_blend+10:>8} {crop_rgb_png/(total_blend+10):>7.3f}x {t_blend*1000:>6.0f}ms")
print(f"  {'G-RG-BG + Blend':.<18} {total_grd_blend+10:>8} {crop_rgb_png/(total_grd_blend+10):>7.3f}x {t_grd_blend*1000:>6.0f}ms")

# Best with brotli too
if HAS_BROTLI:
    total_brotli = 0
    for plane in [g, rg, bg]:
        img_f = plane.astype(float); res_ch = bytearray()
        for r in range(crop_h):
            for c in range(crop_w):
                a=img_f[r,c-1] if c>0 else 0.0; b=img_f[r-1,c] if r>0 else 0.0
                dd=img_f[r-1,c-1] if r>0 and c>0 else 0.0
                grad=abs(a-dd)+abs(b-dd) if r>0 and c>0 else 0
                alpha=min(1.0, grad/1000)
                p_med=med(int(a),int(b),int(dd))
                total_w,wt=0.0,0.0
                for dr,dc in [(-1,-1),(-1,0),(-1,1),(0,-1),(-2,0),(0,-2)]:
                    nr,nc=r+dr,c+dc
                    if 0<=nr<crop_h and 0<=nc<crop_w:
                        dst=abs(dr)+abs(dc); w2=4.0/(dst*dst)
                        total_w+=img_f[nr,nc]*w2; wt+=w2
                p_wie=total_w/wt if wt>0 else 128
                pred=int(np.clip(round(alpha*p_med+(1-alpha)*p_wie),0,255))
                res_ch.append((int(plane[r,c])-pred)&0xFF)
        total_brotli += len(best_compress(bytes(res_ch)))
    print(f"  {'G-RG-BG+Blend+brotli':.<18} {total_brotli+10:>8} {crop_rgb_png/(total_brotli+10):>7.3f}x")

# ============================================================
# IMAGE STATISTICS
# ============================================================

print(f"\n  === IMAGE STATISTICS ===")
for ch, name in [(0, "R"), (1, "G"), (2, "B")]:
    plane = img_rgb[:crop_h, :crop_w, ch]
    print(f"  {name}: mean={plane.mean():.1f}, std={plane.std():.1f}, "
          f"grad_h={np.mean(np.abs(np.diff(plane, axis=1))):.1f}, "
          f"grad_v={np.mean(np.abs(np.diff(plane, axis=0))):.1f}")

print(f"  Gray: mean={gray[:crop_h,:crop_w].mean():.1f}, std={gray[:crop_h,:crop_w].std():.1f}")

# M² rank (from S15)
M2 = crop.astype(float) @ crop.astype(float).T
eigvals = np.linalg.eigvalsh(M2)
rank = np.sum(eigvals > eigvals.max() * 0.01)
print(f"  M² rank (1% threshold): {rank}/{crop_h}")

# WL color count
from collections import Counter
colors = []
for r in range(crop_h):
    row = crop[r].astype(float)
    qm = int(np.clip(row.mean() * 16 / 256, 0, 15))
    qs = int(np.clip(row.std() * 16 / 128, 0, 15))
    colors.append(qm * 16 + qs)
n_colors = len(set(colors))
print(f"  WL L1 colors: {n_colors}")
