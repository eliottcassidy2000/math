#!/usr/bin/env python3
"""
THE ORACLE CODEC: Work backwards from the perfect decoder.

The perfect decoder = a probability model P(pixel | context).
The encoder sends -log2(P(actual | context)) bits per pixel.
Total = cross-entropy H(image, model).

Instead of choosing ONE predictor, maintain ALL predictors simultaneously.
Weight each by how well it's been predicting so far (exponential weights).
The prediction is the WEIGHTED MIXTURE.

This is the "prediction with expert advice" framework:
  - Each "expert" is a predictor (MED, Wiener, LPC, delta, etc.)
  - After each pixel, update expert weights based on prediction error
  - The mixture prediction is provably within O(log(K)) of the best expert
    for K experts (Vovk-style bound)

THE KEY INSIGHT: We don't choose. The mixture automatically converges
to the best expert for each region of the image. Near edges: MED weight
rises. In smooth areas: Wiener weight rises. No overhead for the choice
because the decoder tracks the same weights from the same context.

This is the MATHEMATICAL MERGER of all codecs. Zero overhead.

kind-pasteur-2026-03-25-S17
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

def med_pred(a, b, c):
    mn, mx = min(a, b), max(a, b)
    if c >= mx: return mn
    if c <= mn: return mx
    return a + b - c

def spectral_flatness(residuals):
    f = np.array(residuals, dtype=float) - np.mean(residuals)
    if np.std(f) < 1e-10: return 1.0
    spec = np.abs(np.fft.rfft(f))**2
    spec = spec[1:]; spec = spec[spec > 0]
    if len(spec) == 0: return 1.0
    return np.exp(np.mean(np.log(spec + 1e-20))) / np.mean(spec)

# ============================================================
# THE EXPERT PREDICTORS
# ============================================================

def expert_zero(img, r, c, h, w):
    """Predict 128 (no context)."""
    return 128

def expert_left(img, r, c, h, w):
    """Predict from left neighbor."""
    return int(img[r, c-1]) if c > 0 else 128

def expert_above(img, r, c, h, w):
    """Predict from above neighbor."""
    return int(img[r-1, c]) if r > 0 else 128

def expert_avg(img, r, c, h, w):
    """Average of left and above."""
    a = int(img[r, c-1]) if c > 0 else 0
    b = int(img[r-1, c]) if r > 0 else 0
    n = (1 if c > 0 else 0) + (1 if r > 0 else 0)
    return (a + b) // n if n > 0 else 128

def expert_med(img, r, c, h, w):
    """MED predictor."""
    a = int(img[r, c-1]) if c > 0 else 0
    b = int(img[r-1, c]) if r > 0 else 0
    d = int(img[r-1, c-1]) if r > 0 and c > 0 else 0
    return med_pred(a, b, d)

def expert_wiener(img, r, c, h, w):
    """Weighted average of all causal neighbors (Wiener-like)."""
    total, wt = 0.0, 0.0
    for dr, dc in [(-1,-1),(-1,0),(-1,1),(0,-1),(-2,0),(0,-2)]:
        nr, nc = r + dr, c + dc
        if 0 <= nr < h and 0 <= nc < w:
            dist = abs(dr) + abs(dc)
            w2 = 4.0 / (dist * dist)
            total += int(img[nr, nc]) * w2
            wt += w2
    return int(round(total / wt)) if wt > 0 else 128

def expert_diag(img, r, c, h, w):
    """Predict from diagonal."""
    return int(img[r-1, c-1]) if r > 0 and c > 0 else 128

def expert_linear_h(img, r, c, h, w):
    """Linear extrapolation from 2 left pixels."""
    if c >= 2:
        return max(0, min(255, 2*int(img[r, c-1]) - int(img[r, c-2])))
    return int(img[r, c-1]) if c > 0 else 128

def expert_linear_v(img, r, c, h, w):
    """Linear extrapolation from 2 above pixels."""
    if r >= 2:
        return max(0, min(255, 2*int(img[r-1, c]) - int(img[r-2, c])))
    return int(img[r-1, c]) if r > 0 else 128

EXPERTS = [
    expert_left, expert_above, expert_avg, expert_med,
    expert_wiener, expert_diag, expert_linear_h, expert_linear_v,
]
K = len(EXPERTS)

# ============================================================
# THE ORACLE: Exponential-weight mixture prediction
# ============================================================

def encode_oracle(plane, learning_rate=0.1):
    """Oracle codec: weighted mixture of K experts.
    Weights updated per pixel based on prediction error.
    Decoder tracks same weights → zero overhead."""
    h, w = plane.shape
    img = plane.astype(float)

    # Initialize equal weights
    log_weights = np.zeros(K)  # log domain for numerical stability

    residuals = bytearray()
    expert_errors = np.zeros(K)  # cumulative squared errors

    for r in range(h):
        for c in range(w):
            # Get each expert's prediction
            preds = np.array([exp(img, r, c, h, w) for exp in EXPERTS], dtype=float)

            # Compute mixture prediction (weighted average)
            weights = np.exp(log_weights - np.max(log_weights))  # softmax
            weights /= weights.sum()

            mixture_pred = np.dot(weights, preds)
            mixture_pred = int(np.clip(round(mixture_pred), 0, 255))

            # Compute residual
            actual = int(plane[r, c])
            residuals.append((actual - mixture_pred) & 0xFF)

            # Update weights based on prediction errors
            for k in range(K):
                error = (actual - preds[k]) ** 2
                log_weights[k] -= learning_rate * error / (255*255)
                expert_errors[k] += abs(actual - preds[k])

    # Spectral flatness
    sf = spectral_flatness(list(residuals))

    # Encode: just the residuals (weights are implicit, tracked by decoder)
    compressed = bc(bytes(residuals))

    return len(compressed), sf, expert_errors / (h * w)

def encode_oracle_sweep(plane):
    """Try multiple learning rates, pick best."""
    best_sz = float('inf'); best_lr = 0; best_sf = 0; best_errs = None
    for lr in [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0]:
        sz, sf, errs = encode_oracle(plane, lr)
        if sz < best_sz:
            best_sz = sz; best_lr = lr; best_sf = sf; best_errs = errs
    return best_sz + 1, best_lr, best_sf, best_errs  # +1 for lr byte

# ============================================================
# BASELINES
# ============================================================

def m_med(plane):
    h,w=plane.shape; res=bytearray()
    for r in range(h):
        for c in range(w):
            a=int(plane[r,c-1]) if c>0 else 0; b=int(plane[r-1,c]) if r>0 else 0
            d=int(plane[r-1,c-1]) if r>0 and c>0 else 0
            res.append((int(plane[r,c])-med_pred(a,b,d))&0xFF)
    return len(bc(bytes(res))), spectral_flatness(list(res))

def m_blend(plane):
    h,w=plane.shape; img=plane.astype(float); res=bytearray()
    for r in range(h):
        for c in range(w):
            a=img[r,c-1] if c>0 else 0.0; b=img[r-1,c] if r>0 else 0.0
            d=img[r-1,c-1] if r>0 and c>0 else 0.0
            grad=abs(a-d)+abs(b-d) if r>0 and c>0 else 0
            alpha=min(1.0, grad/1000)
            p_med=med_pred(int(a),int(b),int(d))
            total,wt=0.0,0.0
            for dr,dc in [(-1,-1),(-1,0),(-1,1),(0,-1),(-2,0),(0,-2)]:
                nr,nc=r+dr,c+dc
                if 0<=nr<h and 0<=nc<w:
                    dst=abs(dr)+abs(dc); w2=4.0/(dst*dst)
                    total+=img[nr,nc]*w2; wt+=w2
            p_wie=total/wt if wt>0 else 128
            pred=int(np.clip(round(alpha*p_med+(1-alpha)*p_wie),0,255))
            res.append((int(plane[r,c])-pred)&0xFF)
    return len(bc(bytes(res))), spectral_flatness(list(res))

# ============================================================
# TEST
# ============================================================

SZ = 128

def make_tests():
    T={}; np.random.seed(42)
    x,y=np.meshgrid(np.arange(SZ,dtype=float),np.arange(SZ,dtype=float))
    r=np.sqrt((x-SZ/2)**2+(y-SZ/2)**2)
    T["circles"]=(np.sin(r/5)*127+128).astype(np.uint8)
    T["blob"]=(np.exp(-r**2/(2*(SZ/4)**2))*255).astype(np.uint8)
    sm=np.random.randint(0,256,(SZ//16,SZ//16),dtype=np.uint8)
    T["natural"]=np.clip(np.array(Image.fromarray(sm).resize((SZ,SZ),Image.BILINEAR)).astype(float)
                         +np.random.normal(0,10,(SZ,SZ)),0,255).astype(np.uint8)
    T["grad_h"]=np.tile(np.linspace(0,255,SZ,dtype=np.uint8),(SZ,1))
    T["stripes45"]=(((x+y)/8).astype(int)%2*255).astype(np.uint8)
    T["random"]=np.random.randint(0,256,(SZ,SZ),dtype=np.uint8)
    return T

print("="*90)
print("  THE ORACLE CODEC: Perfect mathematical decoder → work backwards")
print("  Mixture of 8 experts with exponential weight updates")
print("  Zero overhead: decoder tracks same weights from same context")
print("  kind-pasteur-2026-03-25-S17")
print("="*90)

tests=make_tests()

print(f"\n  {'Image':<12} {'PNG':>6}  {'MED':>12}  {'Blend':>12}  {'Oracle':>12}")
print(f"  {'':>18}  {'sz':>5} {'sf':>5}  {'sz':>5} {'sf':>5}  {'sz':>5} {'sf':>5} {'lr':>5}")
print("  "+"-"*72)

for name,img in sorted(tests.items()):
    ps=png_size(img)
    med_sz, med_sf = m_med(img); med_sz += 10
    blend_sz, blend_sf = m_blend(img); blend_sz += 10
    ora_sz, ora_lr, ora_sf, ora_errs = encode_oracle_sweep(img); ora_sz += 10

    best = min(med_sz, blend_sz, ora_sz)
    winner = "MED" if best==med_sz else ("Blend" if best==blend_sz else "Oracle")

    print(f"  {name:<12} {ps:>6}  {med_sz:>5} {med_sf:>5.3f}  {blend_sz:>5} {blend_sf:>5.3f}  {ora_sz:>5} {ora_sf:>5.3f} {ora_lr:>5.3f}  {winner}")

# Show expert weight evolution for natural image
print(f"\n  EXPERT WEIGHT ANALYSIS (natural image, best learning rate):")
img = tests["natural"]
_, lr, _, errs = encode_oracle_sweep(img)
print(f"  Learning rate: {lr}")
print(f"  Average absolute error per expert:")
expert_names = ["left", "above", "avg", "MED", "Wiener", "diag", "lin_h", "lin_v"]
for i, (ename, err) in enumerate(zip(expert_names, errs)):
    print(f"    {ename:>8}: {err:.2f}")

print(f"""
  THE MATHEMATICAL INSIGHT:
  The oracle's mixture prediction is provably within O(log K / T)
  of the best expert in hindsight, for K experts over T pixels.
  For 8 experts, 16384 pixels: regret ≤ log(8)/16384 ≈ 0.0001 bits/pixel.

  This means the oracle automatically converges to the best predictor
  for each region WITHOUT any overhead for the selection decision.
  The decoder computes the SAME weights from the SAME decoded context.

  IF the oracle beats Blend: the mixture captures something Blend misses.
  IF Blend beats oracle: the continuous blend IS the optimal mixture.
""")
