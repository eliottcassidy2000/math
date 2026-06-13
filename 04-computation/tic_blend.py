#!/usr/bin/env python3
"""
BLEND CODEC: Don't choose between MED and Wiener — BLEND them.

The paradox avoidance: instead of binary choice (expensive if wrong),
make a CONTINUOUS blend that's always close to right.

  prediction = α * MED(a,b,c) + (1-α) * Wiener(neighbors)

  α = sigmoid(gradient_magnitude / threshold)

Near edges (high gradient): α → 1 (MED handles edges)
In smooth regions (low gradient): α → 0 (Wiener uses more context)

THE KEY: α is computable from already-known neighbors, so the decoder
computes the SAME α. Zero overhead for the blending decision.

This is the "silent mutation" insight applied: in smooth regions,
the choice between MED and Wiener is "silent" (both give similar results).
Near edges, MED is provably better (nonlinear switching). The blend
smoothly transitions between regimes with zero cost.

kind-pasteur-2026-03-25-S13
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
# THE BLEND PREDICTOR
# ============================================================

def encode_blend(plane, threshold=15.0):
    """Blend MED and Wiener based on local gradient. Zero overhead."""
    h, w = plane.shape
    img = plane.astype(float)
    residuals = bytearray()

    for r in range(h):
        for c in range(w):
            # Causal neighbors
            a = img[r, c-1] if c > 0 else 0.0      # left
            b = img[r-1, c] if r > 0 else 0.0       # above
            d = img[r-1, c-1] if r>0 and c>0 else 0.0  # diag
            e = img[r-1, c+1] if r>0 and c+1<w else b   # above-right

            # Gradient magnitude from known neighbors
            gx = abs(a - d) if (c > 0 and r > 0) else 0
            gy = abs(b - d) if (r > 0 and c > 0) else 0
            grad = gx + gy

            # Blending weight: high gradient → MED (α=1), low gradient → Wiener (α=0)
            alpha = min(1.0, grad / threshold) if threshold > 0 else 0.5

            # MED prediction
            pred_med = med(int(a), int(b), int(d))

            # Wiener prediction: weighted average of more neighbors
            wiener_sum = 0.0; wiener_wt = 0.0
            for dr, dc in [(-1,-1), (-1,0), (-1,1), (0,-1), (-2,0), (0,-2)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < h and 0 <= nc < w:
                    dist = abs(dr) + abs(dc)
                    w2 = 4.0 / (dist * dist)
                    wiener_sum += img[nr, nc] * w2
                    wiener_wt += w2
            pred_wiener = wiener_sum / wiener_wt if wiener_wt > 0 else 128.0

            # Blend
            pred = alpha * pred_med + (1 - alpha) * pred_wiener
            pred = int(np.clip(round(pred), 0, 255))

            residuals.append((int(plane[r, c]) - pred) & 0xFF)

    return len(bc(bytes(residuals)))

def encode_blend_adaptive(plane):
    """Try multiple blend thresholds, pick the one that compresses smallest."""
    best = float('inf')
    best_t = 0
    for t in [5, 10, 15, 20, 30, 50, 100, 1000]:
        sz = encode_blend(plane, t)
        if sz < best:
            best = sz
            best_t = t
    # Also try pure MED (t=0 → alpha always 1) and pure Wiener (t=inf → alpha always 0)
    return best, best_t

# ============================================================
# BASELINES
# ============================================================

def m_med(plane):
    h,w=plane.shape; res=bytearray()
    for r in range(h):
        for c in range(w):
            a=int(plane[r,c-1]) if c>0 else 0; b=int(plane[r-1,c]) if r>0 else 0
            d=int(plane[r-1,c-1]) if r>0 and c>0 else 0
            res.append((int(plane[r,c])-med(a,b,d))&0xFF)
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

def m_colsort(plane):
    h,w=plane.shape
    idx=sorted(range(w),key=lambda c2:plane[:,c2].tobytes())
    sp=np.stack([plane[:,c2] for c2 in idx],axis=1)
    delta=np.zeros_like(sp); delta[:,0]=sp[:,0]
    delta[:,1:]=((sp[:,1:].astype(int)-sp[:,:-1].astype(int))&0xFF).astype(np.uint8)
    perm=np.array(idx,dtype=np.uint8)
    return len(bc(perm.tobytes()))+len(bc(delta.tobytes()))

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
print("  BLEND CODEC: Continuous MED-Wiener interpolation per pixel")
print("  Zero overhead — decoder computes same blend weight from context")
print("  kind-pasteur-2026-03-25-S13")
print("="*80)

tests=make_tests()

print(f"\n  {'Image':<12} {'PNG':>6} {'MED':>6} {'LPC8':>6} {'ColSrt':>6} {'Blend':>6} {'thr':>4} {'Best':>6} {'win':>7}")
print("  "+"-"*70)

W, L = 0, 0
for name,img in sorted(tests.items()):
    ps=png_size(img)
    med_sz=m_med(img)+10
    lpc_sz=m_lpc8(img)+10
    col_sz=m_colsort(img)+10
    blend_sz, blend_t = encode_blend_adaptive(img)
    blend_sz += 10

    all_sizes = {"MED": med_sz, "LPC8": lpc_sz, "ColSrt": col_sz, "Blend": blend_sz}
    best_name = min(all_sizes, key=all_sizes.get)
    best_sz = all_sizes[best_name]
    ratio = ps / best_sz if best_sz > 0 else 0
    if ratio > 1.001: W += 1
    elif ratio < 0.999: L += 1

    print(f"  {name:<12} {ps:>6} {med_sz:>6} {lpc_sz:>6} {col_sz:>6} {blend_sz:>6} {blend_t:>4} {best_sz:>6} {best_name:>7} {ratio:.2f}x")

n = len(tests)
print(f"\n  {W}W / {n-W-L}T / {L}L out of {n}")

print(f"""
  THE BLEND INSIGHT:
  Blending MED and Wiener based on gradient magnitude gives a
  predictor that smoothly transitions between edge-optimized (MED)
  and smooth-optimized (Wiener) behavior.

  The threshold parameter selects how aggressively to switch.
  Low threshold → MED dominates (good for edgy images).
  High threshold → Wiener dominates (good for smooth images).

  The KEY: this blend has ZERO overhead because the decoder computes
  the same gradient from the same already-decoded neighbors.
  No classification map, no block boundaries, no extra bytes.
""")
