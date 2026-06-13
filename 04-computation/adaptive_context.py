#!/usr/bin/env python3
"""
adaptive_context.py — Context-Adaptive Prediction (WL without overhead)

THE KEY INSIGHT: Fractal WL compression showed that CLASSIFYING regions
and using class-specific prediction beats uniform prediction.
But WL coloring has overhead (storing the color map).

SOLUTION: Use the LOCAL GRADIENT CONTEXT as an IMPLICIT classifier.
The decoder can compute the same gradient context from decoded pixels.
No overhead at all — the classification is reproducible.

THIS IS CALIC DONE RIGHT:
  1. Compute gradient context from known neighbors (quantized to N contexts)
  2. For each context, maintain a running ERROR CORRECTION
  3. Apply MED + error correction = adaptive MED
  4. The error correction adapts to each context class

The contexts form a WL-1 coloring of the pixel graph.
The error correction adapts the prediction per WL class.
Combined: WL-class-adaptive prediction at zero overhead.

Also: BLENDING MULTIPLE PREDICTIONS with context-adaptive weights.
  pred = w_med(ctx) * MED + w_avg(ctx) * AVG + w_grad(ctx) * GRADIENT
  where weights adapt based on which prediction worked best in this context.

Author: opus-2026-03-25-S346
"""

import sys, struct, zlib, lzma, io, time, os
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

def _loco(a, b, c):
    a, b, c = int(a), int(b), int(c)
    if c >= max(a, b): return min(a, b)
    if c <= min(a, b): return max(a, b)
    return a + b - c

# ============================================================
# CONTEXT-ADAPTIVE MED (CALIC-inspired, zero overhead)
# ============================================================

def encode_adaptive_med(plane):
    """MED with context-adaptive error correction.
    Zero overhead — decoder reproduces same correction from decoded pixels."""
    H, W = plane.shape
    buf = bytearray()
    # Error correction table: for each quantized gradient context,
    # track the running average prediction error
    N_CTX = 128
    err_sum = np.zeros(N_CTX)
    err_count = np.zeros(N_CTX)

    for r in range(H):
        for c in range(W):
            a = int(plane[r,c-1]) if c > 0 else 0
            b = int(plane[r-1,c]) if r > 0 else 0
            d = int(plane[r-1,c-1]) if r > 0 and c > 0 else 0
            e = int(plane[r-1,c+1]) if r > 0 and c+1 < W else b

            # Base prediction
            pred_base = _loco(a, b, d)

            # Gradient context (quantized)
            gh = abs(a - d) + abs(b - e)
            gv = abs(b - d) + abs(a - (int(plane[r,c-2]) if c >= 2 else a))
            energy = (a + b + d) // 3 if (a + b + d) > 0 else 0
            ctx = (min(gh // 8, 7) * 16 + min(gv // 8, 7)) % N_CTX

            # Error correction
            if err_count[ctx] > 2:
                correction = err_sum[ctx] / err_count[ctx]
                pred = max(0, min(255, int(round(pred_base + correction))))
            else:
                pred = pred_base

            actual = int(plane[r, c])
            buf.append((actual - pred) & 0xFF)

            # Update error statistics
            error = actual - pred_base  # error relative to BASE prediction
            err_sum[ctx] = err_sum[ctx] * 0.95 + error * 0.05
            err_count[ctx] = min(err_count[ctx] + 1, 200)

    return bytes(buf)


def decode_adaptive_med(data, H, W):
    """Decode adaptive MED — reproduces exact same correction."""
    plane = np.zeros((H, W), dtype=np.uint8)
    N_CTX = 128
    err_sum = np.zeros(N_CTX)
    err_count = np.zeros(N_CTX)
    idx = 0

    for r in range(H):
        for c in range(W):
            a = int(plane[r,c-1]) if c > 0 else 0
            b = int(plane[r-1,c]) if r > 0 else 0
            d = int(plane[r-1,c-1]) if r > 0 and c > 0 else 0
            e = int(plane[r-1,c+1]) if r > 0 and c+1 < W else b

            pred_base = _loco(a, b, d)

            gh = abs(a - d) + abs(b - e)
            gv = abs(b - d) + abs(a - (int(plane[r,c-2]) if c >= 2 else a))
            ctx = (min(gh // 8, 7) * 16 + min(gv // 8, 7)) % N_CTX

            if err_count[ctx] > 2:
                correction = err_sum[ctx] / err_count[ctx]
                pred = max(0, min(255, int(round(pred_base + correction))))
            else:
                pred = pred_base

            actual = (int(data[idx]) + pred) & 0xFF
            plane[r, c] = actual
            idx += 1

            error = actual - pred_base
            err_sum[ctx] = err_sum[ctx] * 0.95 + error * 0.05
            err_count[ctx] = min(err_count[ctx] + 1, 200)

    return plane


# ============================================================
# MULTI-PREDICTOR BLEND WITH ADAPTIVE WEIGHTS
# ============================================================

def encode_adaptive_blend(plane):
    """Blend MED, AVG, and GRADIENT predictions with context-adaptive weights.
    Each context learns which predictor works best."""
    H, W = plane.shape
    buf = bytearray()
    N_CTX = 64
    # Track squared error for each predictor in each context
    pred_err = np.ones((N_CTX, 3)) * 100  # [ctx, predictor] = cumulative error²

    for r in range(H):
        for c in range(W):
            a = int(plane[r,c-1]) if c > 0 else 0
            b = int(plane[r-1,c]) if r > 0 else 0
            d = int(plane[r-1,c-1]) if r > 0 and c > 0 else 0

            # Three predictors
            p_med = _loco(a, b, d)
            p_avg = (a + b) >> 1
            p_grad = max(0, min(255, a + b - d))

            # Context
            gh = abs(a - d); gv = abs(b - d)
            ctx = (min(gh // 16, 7) * 8 + min(gv // 16, 7)) % N_CTX

            # Adaptive blend: weight inversely proportional to cumulative error²
            errs = pred_err[ctx]
            w = 1.0 / (errs + 1.0)
            w /= w.sum()

            pred = int(round(w[0] * p_med + w[1] * p_avg + w[2] * p_grad))
            pred = max(0, min(255, pred))

            actual = int(plane[r, c])
            buf.append((actual - pred) & 0xFF)

            # Update error tracking
            for i, p in enumerate([p_med, p_avg, p_grad]):
                e = (actual - p) ** 2
                pred_err[ctx, i] = pred_err[ctx, i] * 0.98 + e * 0.02

    return bytes(buf)


def decode_adaptive_blend(data, H, W):
    plane = np.zeros((H, W), dtype=np.uint8)
    N_CTX = 64
    pred_err = np.ones((N_CTX, 3)) * 100
    idx = 0

    for r in range(H):
        for c in range(W):
            a = int(plane[r,c-1]) if c > 0 else 0
            b = int(plane[r-1,c]) if r > 0 else 0
            d = int(plane[r-1,c-1]) if r > 0 and c > 0 else 0

            p_med = _loco(a, b, d)
            p_avg = (a + b) >> 1
            p_grad = max(0, min(255, a + b - d))

            gh = abs(a - d); gv = abs(b - d)
            ctx = (min(gh // 16, 7) * 8 + min(gv // 16, 7)) % N_CTX

            errs = pred_err[ctx]
            w = 1.0 / (errs + 1.0)
            w /= w.sum()

            pred = int(round(w[0] * p_med + w[1] * p_avg + w[2] * p_grad))
            pred = max(0, min(255, pred))

            actual = (int(data[idx]) + pred) & 0xFF
            plane[r, c] = actual
            idx += 1

            for i, p in enumerate([p_med, p_avg, p_grad]):
                e = (actual - p) ** 2
                pred_err[ctx, i] = pred_err[ctx, i] * 0.98 + e * 0.02

    return plane


# ============================================================
# STANDARD MED (baseline)
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


# ============================================================
# TEST
# ============================================================

def png_size(img):
    from PIL import Image
    pil = Image.fromarray(img, 'L')
    buf = io.BytesIO(); pil.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()


if __name__ == "__main__":
    np.random.seed(42)
    print("=" * 85)
    print("  Context-Adaptive Prediction (WL without overhead)")
    print("  opus-2026-03-25-S346")
    print("=" * 85)

    tests = {}
    N = 128
    tests['gradient'] = np.tile(np.linspace(0,255,N,dtype=np.uint8),(N,1))
    tests['circle'] = np.array([[min(255,int(255*np.exp(-((r-64)**2+(c-64)**2)/2048))) for c in range(N)] for r in range(N)], dtype=np.uint8)
    tests['smooth'] = np.clip(np.array([[int(128+60*np.sin(r/5)*np.cos(c/4)) for c in range(N)] for r in range(N)]),0,255).astype(np.uint8)
    tests['noise'] = np.random.randint(0,256,(N,N),dtype=np.uint8)
    tests['blob'] = np.array([[min(255,int(200*np.exp(-((r-40)**2+(c-80)**2)/800)+100*np.exp(-((r-80)**2+(c-40)**2)/1000))) for c in range(N)] for r in range(N)], dtype=np.uint8)
    x,y = np.meshgrid(np.arange(N,dtype=float),np.arange(N,dtype=float))
    tests['radial'] = np.clip(np.sin(np.sqrt((x-64)**2+(y-64)**2)/5)*127+128,0,255).astype(np.uint8)

    # Real photo
    from PIL import Image as PILImage
    for path in ['/home/e/tron.jpeg']:
        if os.path.exists(path):
            tests['photo'] = np.array(PILImage.open(path).convert('L'))[:128,:128]

    print(f"\n  {'Image':<15} {'PNG':>6} {'MED':>6} {'AdpMED':>6} {'Blend':>6} {'Best':>7} {'vs PNG':>6} {'OK':>4}")
    print("  " + "-" * 72)

    for name, img in sorted(tests.items()):
        H, W = img.shape
        ps = png_size(img)

        # MED
        med = encode_med(img)
        c_med = len(zlib.compress(med, 9))

        # Adaptive MED
        amed = encode_adaptive_med(img)
        rec_amed = decode_adaptive_med(amed, H, W)
        ok_amed = np.array_equal(img, rec_amed)
        c_amed = len(zlib.compress(amed, 9))

        # Adaptive blend
        ablend = encode_adaptive_blend(img)
        rec_ablend = decode_adaptive_blend(ablend, H, W)
        ok_ablend = np.array_equal(img, rec_ablend)
        c_ablend = len(zlib.compress(ablend, 9))

        best = min(c_med, c_amed, c_ablend)
        best_name = "MED" if best == c_med else ("AdpMED" if best == c_amed else "Blend")
        ok = "OK" if (ok_amed and ok_ablend) else "FAIL"

        print(f"  {name:<15} {ps:>6} {c_med:>6} {c_amed:>6} {c_ablend:>6} {best_name:>5} {ps/best:>5.2f}x {ok}")

    print(f"\n  KEY FINDING: Adaptive MED/Blend should beat plain MED by 1-5%")
    print(f"  on images where the prediction error has context-dependent bias.")
    print(f"  This is WL classification at ZERO overhead cost.")
