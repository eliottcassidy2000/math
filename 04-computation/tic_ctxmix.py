#!/usr/bin/env python3
"""
CONTEXT MIXING: The technique that makes PAQ the world's best compressor.

The Blend codec uses 1 context feature (gradient) and 2 experts (MED, Wiener).
Context mixing uses MULTIPLE context features and learns per-context weights.

CONTEXTS (quantized features of the causal neighborhood):
  - Gradient magnitude (4 levels: flat, low, medium, high)
  - Gradient direction (4 quadrants: H, V, diag+, diag-)
  = 16 total context classes

EXPERTS (predictors):
  - Left pixel
  - Above pixel
  - MED(left, above, diag)
  - Wiener (weighted average of 6 neighbors)

PER-CONTEXT WEIGHTS: 16 × 4 = 64 weights, updated online.
The mixer learns which expert is best for EACH context.

ZERO OVERHEAD: decoder computes same contexts from same decoded pixels.
Weights are updated identically by encoder and decoder.

This is the Blend codec generalized:
  Blend: 1 context (gradient) × 2 experts (MED, Wiener) = 2 weights
  CtxMix: 16 contexts × 4 experts = 64 weights (all learned online)

kind-pasteur-2026-03-25-S19
"""
import sys, io, zlib, math, time
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

def zlib9(data): return zlib.compress(data, 9)

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
# CONTEXT CLASSIFICATION
# ============================================================

def get_context(img, r, c, h, w):
    """Classify the local context into one of 16 classes.
    Based on gradient magnitude (4 levels) × direction (4 quadrants)."""
    a = img[r, c-1] if c > 0 else 128.0
    b = img[r-1, c] if r > 0 else 128.0
    d = img[r-1, c-1] if r > 0 and c > 0 else 128.0

    gh = a - d  # horizontal gradient
    gv = b - d  # vertical gradient
    gmag = abs(gh) + abs(gv)

    # Magnitude: 4 levels
    if gmag < 5: mag = 0      # flat
    elif gmag < 15: mag = 1   # low
    elif gmag < 40: mag = 2   # medium
    else: mag = 3              # high

    # Direction: 4 quadrants
    if abs(gh) < 3 and abs(gv) < 3: direction = 0  # isotropic
    elif abs(gh) > abs(gv) * 2: direction = 1       # horizontal
    elif abs(gv) > abs(gh) * 2: direction = 2       # vertical
    else: direction = 3                               # diagonal

    return mag * 4 + direction

# ============================================================
# CONTEXT MIXER
# ============================================================

N_CONTEXTS = 16
N_EXPERTS = 4

def encode_ctxmix(plane, lr=0.1):
    """Context-mixed prediction with online weight learning."""
    h, w = plane.shape
    img = plane.astype(float)
    residuals = bytearray()

    # Per-context weights (initialized equal)
    ctx_weights = np.ones((N_CONTEXTS, N_EXPERTS)) / N_EXPERTS

    for r in range(h):
        for c in range(w):
            # Neighbor values
            a = img[r, c-1] if c > 0 else 0.0
            b = img[r-1, c] if r > 0 else 0.0
            d = img[r-1, c-1] if r > 0 and c > 0 else 0.0

            # Expert predictions
            preds = np.array([
                a,                                    # left
                b,                                    # above
                float(med(int(a), int(b), int(d))),   # MED
                0.0,                                  # Wiener (computed below)
            ])

            # Wiener: weighted average of all causal neighbors
            total, wt = 0.0, 0.0
            for dr, dc in [(-1,-1),(-1,0),(-1,1),(0,-1),(-2,0),(0,-2)]:
                nr, nc = r+dr, c+dc
                if 0 <= nr < h and 0 <= nc < w:
                    dst = abs(dr)+abs(dc); w2 = 4.0/(dst*dst)
                    total += img[nr,nc]*w2; wt += w2
            preds[3] = total/wt if wt > 0 else 128.0

            # Get context
            ctx = get_context(img, r, c, h, w)

            # Weighted prediction
            weights = ctx_weights[ctx]
            mixture = np.dot(weights, preds)
            mixture = int(np.clip(round(mixture), 0, 255))

            actual = int(plane[r, c])
            residuals.append((actual - mixture) & 0xFF)

            # Update weights (multiplicative weights / exponential gradient)
            errors = (actual - preds) ** 2
            # Reduce weight of bad experts, increase good ones
            for k in range(N_EXPERTS):
                ctx_weights[ctx, k] *= np.exp(-lr * errors[k] / (255*255))
            # Normalize
            ws = ctx_weights[ctx].sum()
            if ws > 1e-10:
                ctx_weights[ctx] /= ws

    return zlib9(bytes(residuals))

def decode_ctxmix(data, h, w, lr=0.1):
    """Decode context-mixed prediction — mirrors encoder exactly."""
    raw = zlib.decompress(data)
    plane = np.zeros((h, w), dtype=np.uint8)
    img = np.zeros((h, w), dtype=float)
    ctx_weights = np.ones((N_CONTEXTS, N_EXPERTS)) / N_EXPERTS
    pos = 0

    for r in range(h):
        for c in range(w):
            a = img[r, c-1] if c > 0 else 0.0
            b = img[r-1, c] if r > 0 else 0.0
            d = img[r-1, c-1] if r > 0 and c > 0 else 0.0

            preds = np.array([a, b, float(med(int(a),int(b),int(d))), 0.0])
            total, wt = 0.0, 0.0
            for dr, dc in [(-1,-1),(-1,0),(-1,1),(0,-1),(-2,0),(0,-2)]:
                nr, nc = r+dr, c+dc
                if 0 <= nr < h and 0 <= nc < w:
                    dst = abs(dr)+abs(dc); w2 = 4.0/(dst*dst)
                    total += img[nr,nc]*w2; wt += w2
            preds[3] = total/wt if wt > 0 else 128.0

            ctx = get_context(img, r, c, h, w)
            weights = ctx_weights[ctx]
            mixture = int(np.clip(round(np.dot(weights, preds)), 0, 255))

            actual = (int(raw[pos]) + mixture) & 0xFF
            plane[r, c] = actual
            img[r, c] = float(actual)
            pos += 1

            errors = (actual - preds) ** 2
            for k in range(N_EXPERTS):
                ctx_weights[ctx, k] *= np.exp(-lr * errors[k] / (255*255))
            ws = ctx_weights[ctx].sum()
            if ws > 1e-10: ctx_weights[ctx] /= ws

    return plane

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
    for sigma in [5, 10, 20]:
        sm=np.random.randint(0,256,(SZ//16,SZ//16),dtype=np.uint8)
        base=np.array(Image.fromarray(sm).resize((SZ,SZ),Image.BILINEAR)).astype(float)
        T[f"nat_s{sigma}"]=np.clip(base+np.random.normal(0,sigma,(SZ,SZ)),0,255).astype(np.uint8)
    T["stripes45"]=(((x+y)/8).astype(int)%2*255).astype(np.uint8)
    return T

# Baselines
def m_med(p):
    h,w=p.shape; res=bytearray()
    for r in range(h):
        for c in range(w):
            a=int(p[r,c-1]) if c>0 else 0; b=int(p[r-1,c]) if r>0 else 0
            d=int(p[r-1,c-1]) if r>0 and c>0 else 0
            res.append((int(p[r,c])-med(a,b,d))&0xFF)
    return len(zlib9(bytes(res)))

def m_blend(p):
    h,w=p.shape; img=p.astype(float); res=bytearray()
    for r in range(h):
        for c in range(w):
            a=img[r,c-1] if c>0 else 0.0; b=img[r-1,c] if r>0 else 0.0
            d=img[r-1,c-1] if r>0 and c>0 else 0.0
            grad=abs(a-d)+abs(b-d) if r>0 and c>0 else 0
            alpha=min(1.0, grad/1000)
            p_med=med(int(a),int(b),int(d))
            total,wt=0.0,0.0
            for dr,dc in [(-1,-1),(-1,0),(-1,1),(0,-1),(-2,0),(0,-2)]:
                nr,nc=r+dr,c+dc
                if 0<=nr<h and 0<=nc<w:
                    dst=abs(dr)+abs(dc); w2=4.0/(dst*dst)
                    total+=img[nr,nc]*w2; wt+=w2
            p_wie=total/wt if wt>0 else 128
            pred=int(np.clip(round(alpha*p_med+(1-alpha)*p_wie),0,255))
            res.append((int(p[r,c])-pred)&0xFF)
    return len(zlib9(bytes(res)))

print("="*80)
print("  CONTEXT MIXING: Per-context expert weights (PAQ-inspired)")
print("  16 contexts × 4 experts, online learning, zero overhead")
print("  kind-pasteur-2026-03-25-S19")
print("="*80)

tests=make_tests()

print(f"\n  {'Image':<12} {'PNG':>6} {'MED':>6} {'Blend':>6}", end="")
for lr in [0.01, 0.05, 0.1, 0.5]:
    print(f" {'CM'+str(lr):>8}", end="")
print(f"  {'Best':>6} {'RT':>4}")
print("  "+"-"*75)

for name,img in sorted(tests.items()):
    ps=png_size(img)
    med_sz=m_med(img)+10
    blend_sz=m_blend(img)+10

    cm_results={}
    for lr in [0.01, 0.05, 0.1, 0.5]:
        cm_data=encode_ctxmix(img, lr)
        cm_sz=len(cm_data)+10
        cm_results[lr]=cm_sz

    best_lr=min(cm_results, key=cm_results.get)
    best_cm=cm_results[best_lr]

    # Roundtrip check at best lr
    cm_data=encode_ctxmix(img, best_lr)
    dec=decode_ctxmix(cm_data, SZ, SZ, best_lr)
    ok=np.array_equal(img, dec)

    all_sizes={"MED":med_sz,"Blend":blend_sz,"CtxMix":best_cm}
    winner=min(all_sizes, key=all_sizes.get)
    ratio=ps/all_sizes[winner]

    print(f"  {name:<12} {ps:>6} {med_sz:>6} {blend_sz:>6}", end="")
    for lr in [0.01, 0.05, 0.1, 0.5]:
        print(f" {cm_results[lr]:>8}", end="")
    print(f"  {winner:>6} {'OK' if ok else 'FAIL':>4} {ratio:.2f}x")

print(f"""
  THE CONTEXT MIXING PRINCIPLE:
  16 contexts (gradient magnitude × direction) × 4 experts (left, above, MED, Wiener)
  Weights updated per-pixel with multiplicative weights algorithm.
  Decoder mirrors encoder exactly → zero overhead.

  Compared to Blend:
    Blend: 1 continuous context, 2 experts, hardcoded switching
    CtxMix: 16 discrete contexts, 4 experts, learned switching

  If CtxMix beats Blend: finer context granularity helps.
  If Blend beats CtxMix: continuous context is better than discrete.
""")
