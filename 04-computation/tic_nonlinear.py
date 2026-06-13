#!/usr/bin/env python3
"""
NONLINEAR PREDICTION: All nonlinear functions = compositions of linear functions.

The linear ceiling is 0.54 spectral flatness. To break through, we need
the RIGHT composition of linear maps with simple nonlinearities.

MED is already piecewise-linear (3 pieces). The question: what composition
of linear functions captures the remaining 46% of structure?

APPROACH: Build prediction as a 2-layer "network" with features:
  Layer 1 (linear): extract gradient features from causal neighborhood
    f1 = a - c  (horizontal gradient)
    f2 = b - c  (vertical gradient)
    f3 = a + b - 2c  (Laplacian / curvature)
    f4 = d - c  (diagonal gradient, d = above-right)
    f5 = a - e  (second-order horizontal, e = 2-left)

  Layer 2 (nonlinear): piecewise-linear activations
    g_i = max(0, f_i)   (positive part = ReLU)
    h_i = max(0, -f_i)  (negative part)

  Layer 3 (linear): weighted combination
    pred = w0 + Σ w_i * g_i + Σ v_i * h_i

  Weights learned from the image (like LPC but with nonlinear features).
  Decoder recomputes same features from decoded pixels → zero overhead.

This is EXACTLY a 1-hidden-layer ReLU network applied to the causal neighborhood.
The universal approximation theorem guarantees: with enough features, this
captures ANY continuous function of the neighbors.

kind-pasteur-2026-03-25-S16
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

def spectral_flatness(signal):
    f = signal.astype(float) - signal.mean()
    if np.std(f) < 1e-10: return 1.0
    spec = np.abs(np.fft.rfft(f))**2
    spec = spec[1:]; spec = spec[spec > 0]
    if len(spec) == 0: return 1.0
    return np.exp(np.mean(np.log(spec + 1e-20))) / np.mean(spec)

# ============================================================
# FEATURE EXTRACTION: Linear features from causal neighborhood
# ============================================================

def extract_features(img, r, c, h, w):
    """Extract causal neighbor values and gradient features."""
    a = img[r, c-1] if c > 0 else 0.0        # left
    b = img[r-1, c] if r > 0 else 0.0         # above
    cv = img[r-1, c-1] if r>0 and c>0 else 0.0  # above-left
    d = img[r-1, c+1] if r>0 and c+1<w else b   # above-right
    e = img[r, c-2] if c > 1 else a              # 2-left
    f = img[r-2, c] if r > 1 else b              # 2-above

    # Raw neighbors (linear features)
    raw = [a, b, cv, d, e, f]

    # Gradient features
    gh = a - cv    # horizontal gradient
    gv = b - cv    # vertical gradient
    lap = a + b - 2*cv  # Laplacian
    gd = d - cv    # diagonal gradient
    gh2 = a - e    # second-order horizontal

    # ReLU activations (nonlinear features)
    nonlinear = [
        max(0, gh), max(0, -gh),    # positive/negative horizontal
        max(0, gv), max(0, -gv),    # positive/negative vertical
        max(0, lap), max(0, -lap),  # positive/negative curvature
        max(0, gd), max(0, -gd),    # positive/negative diagonal
        max(0, gh2), max(0, -gh2),  # second-order horizontal
    ]

    return np.array(raw + nonlinear)

# ============================================================
# NONLINEAR PREDICTOR: Learn weights from the image
# ============================================================

def encode_nonlinear(plane):
    """Nonlinear predictor: learn feature weights from the image.
    Features = linear neighbors + ReLU(gradients).
    Weights = least-squares fit on the causal neighborhood."""
    h, w = plane.shape
    img = plane.astype(float)

    # Collect training data: features → pixel value
    X = []; Y = []
    for r in range(2, h):
        for c in range(2, w-1):
            feats = extract_features(img, r, c, h, w)
            X.append(feats)
            Y.append(img[r, c])

    X = np.array(X); Y = np.array(Y)

    # Add bias term
    X_bias = np.column_stack([X, np.ones(len(X))])

    # Solve least squares: Y ≈ X_bias @ weights
    try:
        weights, _, _, _ = np.linalg.lstsq(X_bias, Y, rcond=None)
    except:
        weights = np.zeros(X_bias.shape[1])
        weights[0] = 1.0  # fallback: predict from left

    # Apply prediction with learned weights
    residuals = bytearray()
    for r in range(h):
        for c in range(w):
            feats = extract_features(img, r, c, h, w)
            feats_bias = np.append(feats, 1.0)
            pred = np.dot(weights, feats_bias)
            pred = int(np.clip(round(pred), 0, 255))
            residuals.append((int(plane[r, c]) - pred) & 0xFF)

    # Encode: weights + compressed residuals
    weight_bytes = np.array(weights, dtype=np.float32).tobytes()
    res_comp = bc(bytes(residuals))

    sf = spectral_flatness(np.array(list(residuals)))

    return len(weight_bytes) + len(res_comp), sf

# ============================================================
# ADAPTIVE NONLINEAR: Learn per-block (like adaptive LPC but nonlinear)
# ============================================================

def encode_nonlinear_adaptive(plane, block_rows=16):
    """Adaptive nonlinear: relearn weights every block_rows rows."""
    h, w = plane.shape
    img = plane.astype(float)
    residuals = bytearray()
    all_weights = bytearray()

    for block_start in range(0, h, block_rows):
        block_end = min(block_start + block_rows, h)

        # Collect training from this block (using pixels already processed)
        X = []; Y = []
        for r in range(max(2, block_start), block_end):
            for c in range(2, w-1):
                feats = extract_features(img, r, c, h, w)
                X.append(feats)
                Y.append(img[r, c])

        if len(X) < 20:
            # Not enough data, use previous weights or default
            weights = np.zeros(17); weights[0] = 1.0
        else:
            X = np.array(X); Y = np.array(Y)
            X_bias = np.column_stack([X, np.ones(len(X))])
            try:
                weights, _, _, _ = np.linalg.lstsq(X_bias, Y, rcond=None)
            except:
                weights = np.zeros(X_bias.shape[1]); weights[0] = 1.0

        # Store weights (quantized to float16 for compactness)
        w16 = np.array(weights, dtype=np.float16)
        all_weights.extend(w16.tobytes())

        # Apply to this block
        for r in range(block_start, block_end):
            for c in range(w):
                feats = extract_features(img, r, c, h, w)
                feats_bias = np.append(feats, 1.0)
                pred = np.dot(weights, feats_bias)
                pred = int(np.clip(round(pred), 0, 255))
                residuals.append((int(plane[r, c]) - pred) & 0xFF)

    weight_comp = bc(bytes(all_weights))
    res_comp = bc(bytes(residuals))

    sf = spectral_flatness(np.array(list(residuals)))

    return len(weight_comp) + len(res_comp), sf

# ============================================================
# DEEP COMPOSITION: 2-layer nonlinear (features of features)
# ============================================================

def encode_deep(plane):
    """2-layer nonlinear: features → ReLU → more features → ReLU → prediction.
    Layer 1: gradient features (same as above)
    Layer 2: products and interactions of Layer 1 features
    This captures QUADRATIC interactions (edges meeting at corners, etc.)"""
    h, w = plane.shape
    img = plane.astype(float)

    X = []; Y = []
    for r in range(2, h):
        for c in range(2, w-1):
            f = extract_features(img, r, c, h, w)  # 16 features
            # Layer 2: add pairwise products of the most important features
            # (a*b, a*c, b*c, gh*gv, etc.)
            a, b, cv = f[0], f[1], f[2]
            gh, gv = f[6], f[8]  # positive gradient parts
            layer2 = [
                a * b / 255,      # horizontal-vertical interaction
                gh * gv / 255,    # gradient cross-product
                abs(a - b),       # |horizontal - vertical| (edge detector)
                min(a, b),        # min neighbor (MED-like)
                max(a, b),        # max neighbor (MED-like)
            ]
            feats = np.concatenate([f, layer2])
            X.append(feats)
            Y.append(img[r, c])

    X = np.array(X); Y = np.array(Y)
    X_bias = np.column_stack([X, np.ones(len(X))])

    try:
        weights, _, _, _ = np.linalg.lstsq(X_bias, Y, rcond=None)
    except:
        weights = np.zeros(X_bias.shape[1]); weights[0] = 1.0

    residuals = bytearray()
    for r in range(h):
        for c in range(w):
            f = extract_features(img, r, c, h, w)
            a, b, cv = f[0], f[1], f[2]
            gh, gv = f[6], f[8]
            layer2 = [a*b/255, gh*gv/255, abs(a-b), min(a,b), max(a,b)]
            feats = np.concatenate([f, layer2])
            feats_bias = np.append(feats, 1.0)
            pred = np.dot(weights, feats_bias)
            pred = int(np.clip(round(pred), 0, 255))
            residuals.append((int(plane[r,c]) - pred) & 0xFF)

    weight_bytes = np.array(weights, dtype=np.float32).tobytes()
    res_comp = bc(bytes(residuals))
    sf = spectral_flatness(np.array(list(residuals)))

    return len(weight_bytes) + len(res_comp), sf

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
    sf = spectral_flatness(np.array(list(res)))
    return len(bc(bytes(res))), sf

def m_blend(plane):
    h,w=plane.shape; img=plane.astype(float); res=bytearray()
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
            res.append((int(plane[r,c])-pred)&0xFF)
    sf = spectral_flatness(np.array(list(res)))
    return len(bc(bytes(res))), sf

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
    return T

print("="*90)
print("  NONLINEAR PREDICTION: Compositions of linear functions break the ceiling")
print("  All nonlinear functions = compositions of linear functions (ReLU network)")
print("  kind-pasteur-2026-03-25-S16")
print("="*90)

tests=make_tests()

METHODS=[
    ("MED",             lambda p: m_med(p)),
    ("Blend",           lambda p: m_blend(p)),
    ("Nonlinear",       lambda p: encode_nonlinear(p)),
    ("NL-Adaptive",     lambda p: encode_nonlinear_adaptive(p)),
    ("Deep (2-layer)",  lambda p: encode_deep(p)),
]

print(f"\n  {'Image':<12} {'PNG':>6}", end="")
for mname,_ in METHODS: print(f"  {mname:>14}", end="")
print()
print(f"  {'':>18}", end="")
for mname,_ in METHODS: print(f"  {'sz':>6} {'sf':>6}", end="")
print()
print("  "+"-"*(20+16*len(METHODS)))

for name,img in sorted(tests.items()):
    ps=png_size(img)
    print(f"  {name:<12} {ps:>6}", end="")
    for mname,mfunc in METHODS:
        try:
            sz, sf = mfunc(img)
            sz += 10
            print(f"  {sz:>6} {sf:>6.3f}", end="")
        except Exception as ex:
            print(f"  {'ERR':>6} {'':>6}", end="")
    print()

print(f"""
  KEY: sz = compressed size, sf = spectral flatness (1.0 = optimal whitening)

  THE COMPOSITION THEORY:
    MED = 1-layer piecewise-linear (3 pieces, hardcoded thresholds)
    Nonlinear = 1-layer ReLU network (16 features, learned weights)
    Deep = 2-layer (16 base features + 5 interaction features, learned weights)

  If spectral flatness increases with depth → nonlinear structure IS being captured.
  If compressed size decreases → the weight overhead is worth the prediction gain.
  If sf > 0.54 → WE'VE BROKEN THE LINEAR CEILING.
""")
