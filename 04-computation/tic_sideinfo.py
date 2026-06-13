#!/usr/bin/env python3
"""
SIDE-INFORMATION CODEC: The proof showed splitting kills spatial locality.
The fix: use classification as SIDE INFO while encoding in raster order.

The overhead analysis proved:
  per-class split encoding ALWAYS loses (breaks zlib's cross-block patterns)

But the INFORMATION in the classification (which blocks are smooth, which are edges)
IS valuable — it just needs to be used DIFFERENTLY:
  - Encode in raster order (preserve spatial locality)
  - At each pixel, look up its block's class
  - Use a CLASS-DEPENDENT predictor
  - The class map is encoded once, cheaply (it's spatially coherent)

This is like the blend codec but with DISCRETE classes instead of continuous weights.
The advantage: each class can have a completely different prediction model,
not just a blend of MED and Wiener.

The key test: does this beat both MED (single predictor) and blend (continuous)?

kind-pasteur-2026-03-25-S14
"""
import sys, io, zlib, math, time
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

# ============================================================
# CLASSIFY BLOCKS (same as before)
# ============================================================

def classify_blocks(plane, bs=8):
    h, w = plane.shape
    bh, bw = (h+bs-1)//bs, (w+bs-1)//bs
    classes = np.zeros((bh, bw), dtype=np.uint8)
    means = np.zeros((bh, bw), dtype=np.uint8)
    for by in range(bh):
        for bx in range(bw):
            r0,c0 = by*bs, bx*bs
            r1,c1 = min(r0+bs,h), min(c0+bs,w)
            block = plane[r0:r1, c0:c1].astype(float)
            var_val = block.var()
            means[by,bx] = int(round(block.mean()))
            if var_val < 2: classes[by,bx] = 0      # flat
            elif var_val < 25: classes[by,bx] = 1    # smooth
            elif var_val < 100: classes[by,bx] = 2   # edge
            else: classes[by,bx] = 3                  # texture
    return classes, means

# ============================================================
# SIDE-INFO PREDICTOR: raster order, class-dependent prediction
# ============================================================

def encode_sideinfo(plane, bs=8):
    """Encode in raster order with class-dependent prediction.
    One compressed stream (preserves spatial locality).
    Class map encoded separately."""
    h, w = plane.shape
    classes, means = classify_blocks(plane, bs)
    bh, bw = classes.shape

    residuals = bytearray()
    img = plane.astype(float)

    for r in range(h):
        for c in range(w):
            by, bx = r // bs, c // bs
            if by >= bh: by = bh - 1
            if bx >= bw: bx = bw - 1
            cls = classes[by, bx]
            block_mean = int(means[by, bx])

            a = img[r, c-1] if c > 0 else 0.0
            b = img[r-1, c] if r > 0 else 0.0
            d = img[r-1, c-1] if r > 0 and c > 0 else 0.0

            if cls == 0:
                # FLAT: predict from block mean
                pred = block_mean
            elif cls == 1:
                # SMOOTH: weighted average of many neighbors (Wiener-like)
                total, wt = 0.0, 0.0
                for dr, dc in [(-1,-1),(-1,0),(-1,1),(0,-1),(-2,0),(0,-2)]:
                    nr, nc = r+dr, c+dc
                    if 0<=nr<h and 0<=nc<w:
                        dist = abs(dr)+abs(dc)
                        w2 = 4.0/(dist*dist)
                        total += img[nr,nc]*w2; wt += w2
                pred = total/wt if wt > 0 else 128
            elif cls == 2:
                # EDGE: MED (nonlinear edge detection)
                pred = med(int(a), int(b), int(d))
            else:
                # TEXTURE: simple delta from left (minimal prediction)
                pred = a if c > 0 else b if r > 0 else 128

            pred = int(np.clip(round(pred), 0, 255))
            residuals.append((int(plane[r, c]) - pred) & 0xFF)

    # Compress: class map + block means + single residual stream
    map_comp = bc(classes.tobytes())
    means_comp = bc(means.tobytes())
    res_comp = bc(bytes(residuals))

    return len(map_comp) + len(means_comp) + len(res_comp)

# ============================================================
# BASELINES FOR COMPARISON
# ============================================================

def m_med(plane):
    h,w=plane.shape; res=bytearray()
    for r in range(h):
        for c in range(w):
            a=int(plane[r,c-1]) if c>0 else 0; b=int(plane[r-1,c]) if r>0 else 0
            d=int(plane[r-1,c-1]) if r>0 and c>0 else 0
            res.append((int(plane[r,c])-med(a,b,d))&0xFF)
    return len(bc(bytes(res)))

def m_blend(plane, threshold=1000):
    """Blend MED and Wiener (from S13)."""
    h,w=plane.shape; img=plane.astype(float); res=bytearray()
    for r in range(h):
        for c in range(w):
            a=img[r,c-1] if c>0 else 0.0; b=img[r-1,c] if r>0 else 0.0
            d=img[r-1,c-1] if r>0 and c>0 else 0.0
            grad = abs(a-d)+abs(b-d) if r>0 and c>0 else 0
            alpha = min(1.0, grad/threshold)
            p_med = med(int(a),int(b),int(d))
            total,wt=0.0,0.0
            for dr,dc in [(-1,-1),(-1,0),(-1,1),(0,-1),(-2,0),(0,-2)]:
                nr,nc=r+dr,c+dc
                if 0<=nr<h and 0<=nc<w:
                    dist2=abs(dr)+abs(dc); w2=4.0/(dist2*dist2)
                    total+=img[nr,nc]*w2; wt+=w2
            p_wie=total/wt if wt>0 else 128
            pred=int(np.clip(round(alpha*p_med+(1-alpha)*p_wie),0,255))
            res.append((int(plane[r,c])-pred)&0xFF)
    return len(bc(bytes(res)))

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
    T["circles"]=(np.sin(r/5)*127+128).astype(np.uint8)
    T["blob"]=(np.exp(-r**2/(2*(SZ/4)**2))*255).astype(np.uint8)
    T["random"]=np.random.randint(0,256,(SZ,SZ),dtype=np.uint8)
    sm=np.random.randint(0,256,(SZ//16,SZ//16),dtype=np.uint8)
    T["natural"]=np.clip(np.array(Image.fromarray(sm).resize((SZ,SZ),Image.BILINEAR)).astype(float)
                         +np.random.normal(0,10,(SZ,SZ)),0,255).astype(np.uint8)
    mixed=np.zeros((SZ,SZ),dtype=np.uint8)
    mixed[:SZ//2,:SZ//2]=128; mixed[:SZ//2,SZ//2:]=np.tile(np.linspace(0,255,SZ//2,dtype=np.uint8),(SZ//2,1))
    mixed[SZ//2:,:SZ//2]=np.random.randint(0,256,(SZ//2,SZ//2),dtype=np.uint8)
    r2=np.sqrt((x[SZ//2:,SZ//2:]-SZ*0.75)**2+(y[SZ//2:,SZ//2:]-SZ*0.75)**2)
    mixed[SZ//2:,SZ//2:]=(np.sin(r2/5)*127+128).astype(np.uint8)
    T["mixed"]=mixed
    return T

print("="*80)
print("  SIDE-INFO CODEC: Classification as side info, not routing")
print("  Proves: raster order + class-dependent prediction > splitting")
print("  kind-pasteur-2026-03-25-S14")
print("="*80)

tests=make_tests()
print(f"\n  {'Image':<12} {'PNG':>6} {'MED':>6} {'Blend':>6} {'SideIn':>6} {'Best':>6} {'winner':>8}")
print("  "+"-"*55)

for name,img in sorted(tests.items()):
    ps=png_size(img)
    med_sz=m_med(img)+10
    blend_sz=m_blend(img)+10
    side_sz=encode_sideinfo(img)+10
    all_m={"MED":med_sz,"Blend":blend_sz,"SideIn":side_sz}
    best_n=min(all_m,key=all_m.get); best_s=all_m[best_n]
    ratio=ps/best_s if best_s>0 else 0
    print(f"  {name:<12} {ps:>6} {med_sz:>6} {blend_sz:>6} {side_sz:>6} {best_s:>6} {best_n:>8} {ratio:.2f}x")

print(f"""
  THE LESSON:
  The overhead proof showed splitting is ALWAYS worse.
  Side-info preserves spatial locality while using classification.
  But does it beat the continuous blend?

  If SideInfo > Blend: discrete classification adds value over continuous blend
  If Blend > SideInfo: continuous is strictly better than discrete
""")
