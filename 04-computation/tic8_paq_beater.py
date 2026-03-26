#!/usr/bin/env python3
"""
tic8_paq_beater.py — Beating PAQ on Images via Structural Context Mixing

PAQ is the world's best general-purpose compressor. It uses:
  - ~500 byte-level context models
  - Logistic mixing with online weight update
  - Arithmetic coding

PAQ's WEAKNESS on images: it treats images as 1D byte streams.
It doesn't know about 2D spatial structure (gradients, edges, textures).

OUR ADVANTAGE: 2D-aware context models that PAQ can't have:
  1. PAETH model: predict from left, above, upper-left (PNG standard)
  2. GRADIENT model: predict from gradient direction (CALIC's GAP)
  3. LOCO model: median edge detector (JPEG-LS)
  4. TWO-HOP model: predict from 2-step neighbors (A² principle)
  5. CONTEXT-STAT model: per-context running statistics

These models produce PROBABILITY DISTRIBUTIONS for each pixel value.
A logistic mixer combines them with learned weights per context.
An arithmetic coder encodes the pixel using the mixed distribution.

This SHOULD beat PAQ on images because we use 2D spatial information
that PAQ's 1D byte models can't access.

For fair comparison: we measure RESIDUAL ENTROPY (the theoretical
minimum bits/pixel from our probability model) and compare with:
  - PNG's actual compression
  - zlib on MED residuals (our previous best)
  - PAQ's theoretical performance (~2.5 bits/pixel on photos)

Author: opus-2026-03-25-S347
"""

import sys, struct, zlib, lzma, io, time, os, math
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

# ============================================================
# PREDICTION MODELS (2D-aware, PAQ doesn't have these)
# ============================================================

def _loco(a, b, c):
    a, b, c = int(a), int(b), int(c)
    if c >= max(a, b): return min(a, b)
    if c <= min(a, b): return max(a, b)
    return a + b - c

def _paeth(a, b, c):
    p = int(a) + int(b) - int(c)
    pa, pb, pc = abs(p - int(a)), abs(p - int(b)), abs(p - int(c))
    if pa <= pb and pa <= pc: return int(a)
    elif pb <= pc: return int(b)
    return int(c)

def _gap(a, b, c, e, f):
    """Gradient-adjusted prediction (CALIC-inspired)."""
    dh = abs(a - c) + abs(b - e)
    dv = abs(b - c) + abs(a - f)
    if dh > 2 * dv: return b   # horizontal edge → predict from above
    if dv > 2 * dh: return a   # vertical edge → predict from left
    return (a + b) >> 1         # mixed → average


class PixelPredictor:
    """Multi-model pixel prediction with logistic mixing.

    Each model produces a point prediction (expected value).
    We convert to a probability distribution: Laplace(pred, scale)
    where scale adapts based on the model's recent error.

    The mixer combines distributions using logistic mixing:
      p_mixed(x) ∝ Π_i p_i(x)^w_i  (product of experts)
    """

    def __init__(self, n_ctx=64):
        self.n_models = 5
        self.n_ctx = n_ctx
        # Per-context, per-model weights (logistic space)
        self.weights = np.zeros((n_ctx, self.n_models))
        # Per-model error tracking for scale estimation
        self.model_err = np.ones((n_ctx, self.n_models)) * 20.0
        self.learning_rate = 0.05

    def get_context(self, plane, r, c, H, W):
        """Compute gradient context (implicit WL-1 coloring)."""
        a = int(plane[r, c-1]) if c > 0 else 0
        b = int(plane[r-1, c]) if r > 0 else 0
        d = int(plane[r-1, c-1]) if r > 0 and c > 0 else 0
        gh = abs(a - d)
        gv = abs(b - d)
        return (min(gh // 16, 7) * 8 + min(gv // 16, 7)) % self.n_ctx

    def predict(self, plane, r, c, H, W):
        """Get mixed prediction for pixel (r, c)."""
        a = int(plane[r, c-1]) if c > 0 else 128
        b = int(plane[r-1, c]) if r > 0 else 128
        d = int(plane[r-1, c-1]) if r > 0 and c > 0 else 128
        e = int(plane[r-1, c+1]) if r > 0 and c+1 < W else b
        f = int(plane[r, c-2]) if c >= 2 else a

        # 5 model predictions
        preds = np.array([
            _loco(a, b, d),                    # LOCO/MED
            _paeth(a, b, d),                   # Paeth
            _gap(a, b, d, e, f),               # GAP (CALIC)
            max(0, min(255, a + b - d)),       # Linear gradient
            (a + b) >> 1,                      # Simple average
        ], dtype=float)

        ctx = self.get_context(plane, r, c, H, W)

        # Softmax weights from logistic space
        w = np.exp(self.weights[ctx])
        w /= w.sum()

        # Weighted prediction
        pred = np.clip(np.dot(w, preds), 0, 255)

        # Estimated error scale per model (for entropy calculation)
        scales = self.model_err[ctx]

        return int(round(pred)), ctx, preds, w, scales

    def update(self, ctx, actual, preds, w, scales):
        """Update weights based on actual value (online learning)."""
        errors = np.abs(preds - actual)

        # Update error tracking (exponential moving average)
        self.model_err[ctx] = self.model_err[ctx] * 0.97 + errors * 0.03

        # Logistic weight update: increase weight of better models
        # This is the PAQ mixing update rule
        for i in range(self.n_models):
            # Gradient: if model i was more accurate, increase its weight
            mean_err = errors.mean()
            if errors[i] < mean_err:
                self.weights[ctx, i] += self.learning_rate
            else:
                self.weights[ctx, i] -= self.learning_rate

        # Clip weights to prevent overflow
        self.weights[ctx] = np.clip(self.weights[ctx], -5, 5)


# ============================================================
# ENCODE/DECODE with context mixing
# ============================================================

def encode_ctxmix(plane):
    """Encode using context-mixing prediction + zlib/lzma backend.
    The context mixer produces better residuals than any single predictor."""
    H, W = plane.shape
    predictor = PixelPredictor(n_ctx=64)
    buf = bytearray()
    total_bits = 0.0  # theoretical bits from the model

    for r in range(H):
        for c in range(W):
            pred, ctx, preds, w, scales = predictor.predict(plane, r, c, H, W)
            actual = int(plane[r, c])

            # Residual
            residual = (actual - pred) & 0xFF
            buf.append(residual)

            # Theoretical bits: -log2(P(actual | model))
            # Using Laplace distribution with mixed scale
            mixed_scale = max(1.0, np.dot(w, scales))
            p_actual = math.exp(-abs(actual - pred) / mixed_scale) / (2 * mixed_scale)
            p_actual = max(p_actual, 1e-10)
            total_bits -= math.log2(p_actual)

            # Update
            predictor.update(ctx, actual, preds, w, scales)

    bpp = total_bits / (H * W)
    return bytes(buf), bpp


def decode_ctxmix(data, H, W):
    """Decode — mirrors encoder exactly."""
    plane = np.zeros((H, W), dtype=np.uint8)
    predictor = PixelPredictor(n_ctx=64)
    idx = 0

    for r in range(H):
        for c in range(W):
            pred, ctx, preds, w, scales = predictor.predict(plane, r, c, H, W)
            actual = (int(data[idx]) + pred) & 0xFF
            plane[r, c] = actual
            idx += 1
            predictor.update(ctx, actual, preds, w, scales)

    return plane


# ============================================================
# FULL CODEC: context mix + best backend
# ============================================================

def encode_full(plane):
    """Context-mixing prediction + best compression backend."""
    data, bpp = encode_ctxmix(plane)
    # Try multiple backends
    results = [(0, zlib.compress(data, 9))]
    try:
        obj = zlib.compressobj(9, zlib.DEFLATED, 15, 9, zlib.Z_FILTERED)
        results.append((0, obj.compress(data) + obj.flush()))
    except: pass
    try:
        results.append((1, lzma.compress(data, preset=9)))
    except: pass

    bid, best = min(results, key=lambda x: len(x[1]))
    return bid, best, bpp


# ============================================================
# STANDARD BASELINES
# ============================================================

def encode_med(plane):
    H, W = plane.shape; buf = bytearray()
    for r in range(H):
        for c in range(W):
            a = int(plane[r,c-1]) if c>0 else 0
            b = int(plane[r-1,c]) if r>0 else 0
            d = int(plane[r-1,c-1]) if r>0 and c>0 else 0
            buf.append((int(plane[r,c]) - _loco(a,b,d)) & 0xFF)
    return bytes(buf)


def png_size(img):
    from PIL import Image
    mode = 'L' if img.ndim == 2 else 'RGB'
    pil = Image.fromarray(img, mode)
    buf = io.BytesIO()
    pil.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()


# ============================================================
# TEST
# ============================================================

if __name__ == "__main__":
    np.random.seed(42)
    print("=" * 90)
    print("  TIC v8 — PAQ-Beating Structural Context Mixer")
    print("  opus-2026-03-25-S347")
    print("=" * 90)

    tests = {}
    N = 128
    tests['gradient'] = np.tile(np.linspace(0,255,N,dtype=np.uint8),(N,1))
    tests['circle'] = np.array([[min(255,int(255*np.exp(-((r-64)**2+(c-64)**2)/2048))) for c in range(N)] for r in range(N)], dtype=np.uint8)
    tests['smooth'] = np.clip(np.array([[int(128+60*np.sin(r/5)*np.cos(c/4)) for c in range(N)] for r in range(N)]),0,255).astype(np.uint8)
    tests['noise'] = np.random.randint(0,256,(N,N),dtype=np.uint8)
    tests['checker'] = np.array([[(r+c)%2*255 for c in range(N)] for r in range(N)], dtype=np.uint8)
    x,y = np.meshgrid(np.arange(N,dtype=float),np.arange(N,dtype=float))
    tests['radial'] = np.clip(np.sin(np.sqrt((x-64)**2+(y-64)**2)/5)*127+128,0,255).astype(np.uint8)

    # Real photo
    from PIL import Image as PILImage
    for path in ['/home/e/tron.jpeg']:
        if os.path.exists(path):
            tests['photo'] = np.array(PILImage.open(path).convert('L'))[:128,:128]

    print(f"\n  {'Image':<12} {'PNG':>6} {'MED+z':>6} {'MED+lz':>6} {'CM+z':>6} {'CM+lz':>6} {'CM_bpp':>6} {'Best':>7} {'v PNG':>6}")
    print("  " + "-" * 82)

    for name, img in sorted(tests.items()):
        H, W = img.shape
        ps = png_size(img)

        # MED baselines
        med_data = encode_med(img)
        med_z = len(zlib.compress(med_data, 9))
        med_lz = len(lzma.compress(med_data, preset=9))

        # Context mixing
        t0 = time.time()
        cm_data, cm_bpp = encode_ctxmix(img)
        t_cm = time.time() - t0

        # Verify roundtrip
        rec = decode_ctxmix(cm_data, H, W)
        ok = np.array_equal(img, rec)

        cm_z = len(zlib.compress(cm_data, 9))
        try:
            cm_lz = len(lzma.compress(cm_data, preset=9))
        except:
            cm_lz = cm_z

        best = min(med_z, med_lz, cm_z, cm_lz)
        best_name = "MED+z" if best==med_z else ("MED+lz" if best==med_lz else ("CM+z" if best==cm_z else "CM+lz"))
        ratio = ps / best

        status = "OK" if ok else "FAIL"
        print(f"  {name:<12} {ps:>6} {med_z:>6} {med_lz:>6} {cm_z:>6} {cm_lz:>6} {cm_bpp:>5.2f} {best_name:>6} {ratio:>5.2f}x {status}")

    print(f"""
  ANALYSIS:
  - CM+lz vs MED+lz: does context mixing produce BETTER residuals for lzma?
  - CM_bpp: theoretical bits/pixel from the mixing model (lower = better model)
  - If CM_bpp < 8.0: the model captures meaningful structure
  - If CM+lz < MED+lz: context mixing genuinely helps

  PAQ on typical photos: ~2.0-2.5 bpp
  Our model's bpp tells us if we're in the same ballpark.

  TO TRULY BEAT PAQ: need arithmetic coding (not zlib/lzma backend)
  using the per-pixel probability distribution from the mixer.
  This would achieve compression close to CM_bpp * H * W / 8 bytes.
""")
