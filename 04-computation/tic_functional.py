#!/usr/bin/env python3
"""
FUNCTIONAL COMPRESSION: Compression as pure function composition.

In Haskell:
  compress   = encode . residualize . predict . reorder
  decompress = unreorder . unpredict . unresidual . decode

Each step is a BIJECTION on the data. The composition is a bijection.
The art: choose the bijections so the final representation has minimum entropy.

THE KEY SYMMETRY: compress and decompress are inverses.
This means every operation must be INVERTIBLE.

THE DEEPER SYMMETRY: The image f(x,y) has spatial autocorrelation.
The OPTIMAL transform diagonalizes this autocorrelation — making pixels
independent. For Gaussian sources, this is the KLT (Karhunen-Loève Transform).
For images, it's approximately the DCT. But for LOSSLESS compression,
we need an INTEGER version.

THE FUNCTIONAL PIPELINE:
  1. REORDER: permute pixels to maximize serial correlation
     (scan order choice — raster, diagonal, ring, etc.)
  2. PREDICT: subtract prediction from each pixel
     (MED, LPC, Wiener — a linear/nonlinear filter)
  3. RESIDUALIZE: map residuals to a form with minimum entropy
     (zigzag, score conditioning, context modeling)
  4. ENCODE: entropy code the residuals
     (Huffman, arithmetic, ANS, or just zlib)

Each step is a function. The COMPOSITION is the codec.
The PRODUCT of entropies at each step gives the total compressed size.

But here's the insight: steps 1-3 don't change the TOTAL entropy.
They only change the CONDITIONAL entropy that step 4 can exploit.
The best codec makes each residual INDEPENDENT and IDENTICALLY DISTRIBUTED.

WHAT DOES "INDEPENDENT" MEAN FOR IMAGES?
After perfect prediction, residuals should be white noise.
The spectrum of residuals should be flat.
Any remaining structure = wasted compression opportunity.

THIS SCRIPT MEASURES THE SPECTRAL FLATNESS OF RESIDUALS
for each predictor. Higher flatness = closer to optimal.
The predictor with flattest residual spectrum is provably best.

kind-pasteur-2026-03-25-S11
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

def med(a, b, c):
    mn, mx = min(a, b), max(a, b)
    if c >= mx: return mn
    if c <= mn: return mx
    return a + b - c

# ============================================================
# THE FUNCTIONAL DECOMPOSITION
# ============================================================

# STEP 1: REORDER — choose the permutation of pixels
def reorder_raster(plane):
    """Identity permutation (row-major)."""
    return plane.ravel(), plane.shape

def reorder_snake(plane):
    """Snake scan: alternate row direction for better 1D locality."""
    h, w = plane.shape
    result = np.empty(h * w, dtype=np.uint8)
    for r in range(h):
        if r % 2 == 0:
            result[r*w:(r+1)*w] = plane[r]
        else:
            result[r*w:(r+1)*w] = plane[r, ::-1]
    return result, plane.shape

def reorder_hilbert(plane):
    """Hilbert curve for maximum 2D locality preservation."""
    h, w = plane.shape
    # Simple Hilbert for power-of-2 sizes
    def hilbert_d2xy(n, d):
        x = y = 0
        s = 1
        while s < n:
            rx = 1 if (d & 2) else 0
            ry = 1 if ((d & 1) ^ rx) else 0
            if ry == 0:
                if rx == 1: x, y = s-1-x, s-1-y
                x, y = y, x
            x += s * rx; y += s * ry
            d >>= 2; s <<= 1
        return x, y

    n = max(h, w)
    # Round up to power of 2
    n2 = 1
    while n2 < n: n2 <<= 1

    result = []
    for d in range(n2 * n2):
        x, y = hilbert_d2xy(n2, d)
        if y < h and x < w:
            result.append(plane[y, x])
    return np.array(result, dtype=np.uint8), plane.shape

# STEP 2: PREDICT — subtract prediction (bijection via mod 256)
def predict_none(signal, shape):
    """No prediction."""
    return signal.copy()

def predict_delta(signal, shape):
    """First-order difference (like audio delta coding)."""
    result = np.empty_like(signal)
    result[0] = signal[0]
    result[1:] = (signal[1:].astype(int) - signal[:-1].astype(int)) & 0xFF
    return result

def predict_lpc(signal, shape, order=4):
    """Linear predictive coding with Levinson-Durbin."""
    n = len(signal)
    sig_f = signal.astype(float)

    # Compute autocorrelation
    r = np.correlate(sig_f - sig_f.mean(), sig_f - sig_f.mean(), mode='full')
    r = r[n-1:n+order]
    if r[0] < 1e-10:
        return predict_delta(signal, shape)

    # Levinson-Durbin
    a = np.zeros(order); e = r[0]
    for i in range(order):
        if e < 1e-10: break
        lam = r[i+1]
        for j in range(i): lam -= a[j] * r[i-j]
        k = lam / e
        a_new = np.zeros(order); a_new[i] = k
        for j in range(i): a_new[j] = a[j] - k * a[i-1-j]
        a = a_new; e *= (1 - k*k)

    # Apply prediction
    result = np.empty(n, dtype=np.uint8)
    for i in range(n):
        pred = sum(a[j] * sig_f[i-1-j] for j in range(min(i, order)))
        result[i] = (int(signal[i]) - int(np.clip(round(pred), 0, 255))) & 0xFF

    return result

def predict_adaptive_lpc(signal, shape, order=4, block=1024):
    """Adaptive LPC: recompute coefficients per block."""
    n = len(signal)
    sig_f = signal.astype(float)
    result = np.empty(n, dtype=np.uint8)

    for start in range(0, n, block):
        end = min(start + block, n)
        blk = sig_f[start:end]
        blen = len(blk)

        if blen < order + 1:
            # Too short for LPC, use delta
            for i in range(blen):
                idx = start + i
                pred = sig_f[idx-1] if idx > 0 else 0
                result[idx] = (int(signal[idx]) - int(np.clip(round(pred), 0, 255))) & 0xFF
            continue

        # Local autocorrelation
        r = np.correlate(blk - blk.mean(), blk - blk.mean(), mode='full')
        r = r[blen-1:blen+order]
        if r[0] < 1e-10:
            a = np.zeros(order)
        else:
            a = np.zeros(order); e = r[0]
            for ii in range(order):
                if e < 1e-10: break
                lam = r[ii+1]
                for j in range(ii): lam -= a[j] * r[ii-j]
                k = lam / e
                a_new = np.zeros(order); a_new[ii] = k
                for j in range(ii): a_new[j] = a[j] - k * a[ii-1-j]
                a = a_new; e *= (1 - k*k)

        for i in range(blen):
            idx = start + i
            pred = sum(a[j] * sig_f[idx-1-j] for j in range(min(idx, order)))
            result[idx] = (int(signal[idx]) - int(np.clip(round(pred), 0, 255))) & 0xFF

    return result

# STEP 3: RESIDUALIZE — transform residuals for minimum entropy
def residualize_none(residuals):
    """Identity."""
    return residuals

def residualize_zigzag(residuals):
    """Map signed residuals to unsigned via zigzag (0,-1,1,-2,2,...)."""
    result = np.empty_like(residuals)
    for i in range(len(residuals)):
        v = int(residuals[i])
        signed = v if v <= 128 else v - 256
        if signed >= 0:
            result[i] = min(255, 2 * signed)
        else:
            result[i] = min(255, 2 * (-signed) - 1)
    return result

def residualize_xor_prev(residuals):
    """XOR with previous residual — decorrelates if residuals are correlated."""
    result = np.empty_like(residuals)
    result[0] = residuals[0]
    for i in range(1, len(residuals)):
        result[i] = residuals[i] ^ residuals[i-1]
    return result

# STEP 4: ENCODE — entropy coding (zlib/brotli as black box)
def encode(data):
    return bc(bytes(data))

# ============================================================
# THE COMPOSITION ENGINE
# ============================================================

REORDERS = [
    ("raster", reorder_raster),
    ("snake",  reorder_snake),
]

PREDICTORS = [
    ("none",        lambda s, sh: predict_none(s, sh)),
    ("delta",       lambda s, sh: predict_delta(s, sh)),
    ("lpc4",        lambda s, sh: predict_lpc(s, sh, 4)),
    ("lpc8",        lambda s, sh: predict_lpc(s, sh, 8)),
    ("alpc4",       lambda s, sh: predict_adaptive_lpc(s, sh, 4, 512)),
    ("alpc4_1k",    lambda s, sh: predict_adaptive_lpc(s, sh, 4, 1024)),
]

RESIDUALIZERS = [
    ("none",    residualize_none),
    ("zigzag",  residualize_zigzag),
    ("xor",     residualize_xor_prev),
]

def try_all_compositions(plane):
    """Try all combinations of reorder × predict × residualize × encode.
    Return the smallest and its composition."""
    best_size = float('inf')
    best_comp = None

    for rname, rfunc in REORDERS:
        signal, shape = rfunc(plane)
        for pname, pfunc in PREDICTORS:
            residuals = pfunc(signal, shape)
            for zname, zfunc in RESIDUALIZERS:
                transformed = zfunc(residuals)
                compressed = encode(transformed)
                size = len(compressed)
                if size < best_size:
                    best_size = size
                    best_comp = f"{rname}|{pname}|{zname}"

    return best_size, best_comp

# ============================================================
# SPECTRAL FLATNESS: measure how close to optimal each predictor is
# ============================================================

def spectral_flatness(signal):
    """Spectral flatness = geometric mean / arithmetic mean of power spectrum.
    1.0 = perfectly white noise (optimal residual).
    <1.0 = colored noise (prediction can still improve)."""
    if len(signal) < 16:
        return 0.0
    f_signal = signal.astype(float)
    f_signal -= f_signal.mean()
    if np.std(f_signal) < 1e-10:
        return 1.0  # constant signal

    spectrum = np.abs(np.fft.rfft(f_signal))**2
    spectrum = spectrum[1:]  # remove DC
    spectrum = spectrum[spectrum > 0]  # remove zeros

    if len(spectrum) == 0:
        return 1.0

    log_mean = np.mean(np.log(spectrum + 1e-20))
    arith_mean = np.mean(spectrum)

    if arith_mean < 1e-20:
        return 1.0

    return np.exp(log_mean) / arith_mean

# ============================================================
# BENCHMARK
# ============================================================

SZ = 128

def make_tests():
    T = {}; np.random.seed(42)
    x, y = np.meshgrid(np.arange(SZ, dtype=float), np.arange(SZ, dtype=float))
    r = np.sqrt((x-SZ/2)**2 + (y-SZ/2)**2)
    T["solid"] = np.full((SZ,SZ), 128, dtype=np.uint8)
    T["grad_h"] = np.tile(np.linspace(0,255,SZ,dtype=np.uint8),(SZ,1))
    T["grad_d"] = ((x+y)*255//(2*SZ-2)).astype(np.uint8)
    T["circles"] = (np.sin(r/5)*127+128).astype(np.uint8)
    T["blob"] = (np.exp(-r**2/(2*(SZ/4)**2))*255).astype(np.uint8)
    T["random"] = np.random.randint(0,256,(SZ,SZ),dtype=np.uint8)
    sm = np.random.randint(0,256,(SZ//16,SZ//16),dtype=np.uint8)
    T["natural"] = np.clip(np.array(Image.fromarray(sm).resize((SZ,SZ),Image.BILINEAR)).astype(float)
                           +np.random.normal(0,10,(SZ,SZ)),0,255).astype(np.uint8)
    T["stripes45"] = (((x+y)/8).astype(int)%2*255).astype(np.uint8)
    return T

print("=" * 80)
print("  FUNCTIONAL COMPRESSION: Pure composition of bijections")
print("  compress = encode . residualize . predict . reorder")
print("  kind-pasteur-2026-03-25-S11")
print("=" * 80)

tests = make_tests()

# Part 1: Try all compositions
print(f"\n  PART 1: Best composition per image (2 reorders × 6 predictors × 3 residualizers = 36 combos)")
print(f"  {'Image':<12} {'PNG':>6} {'BEST':>6} {'ratio':>6} {'composition':>30}")
print("  " + "-" * 65)

W = 0
for name, img in sorted(tests.items()):
    ps = png_size(img)
    csz, comp = try_all_compositions(img)
    total = csz + 10
    r = ps / total if total > 0 else 0
    if r > 1.001: W += 1
    print(f"  {name:<12} {ps:>6} {total:>6} {r:>6.2f} {comp:>30}")

print(f"\n  {W}/{len(tests)} beat PNG")

# Part 2: Spectral flatness of residuals
print(f"\n  PART 2: SPECTRAL FLATNESS of residuals (1.0 = optimal white noise)")
print(f"  {'Image':<12}", end="")
for pname, _ in PREDICTORS:
    print(f" {pname:>8}", end="")
print(f"  {'raw':>8}")
print("  " + "-" * (14 + 9 * (len(PREDICTORS) + 1)))

for name, img in sorted(tests.items()):
    signal, _ = reorder_raster(img)
    print(f"  {name:<12}", end="")
    for pname, pfunc in PREDICTORS:
        res = pfunc(signal, img.shape)
        sf = spectral_flatness(res)
        print(f" {sf:>8.4f}", end="")
    sf_raw = spectral_flatness(signal)
    print(f"  {sf_raw:>8.4f}")

print(f"""
  INTERPRETATION:
  - Spectral flatness close to 1.0 = residuals are white noise = OPTIMAL
  - Below 1.0 = residuals have structure = room for better prediction
  - The BEST predictor for each image makes flattest residuals

  The functional insight: compression quality = how well the composition
  (reorder . predict . residualize) whitens the signal.
  Perfect whitening → entropy coding is optimal → minimum size.
""")

# Part 3: Which step contributes most?
print(f"  PART 3: Contribution of each step (bits per pixel)")
print(f"  {'Image':<12} {'raw':>6} {'reord':>6} {'pred':>6} {'resid':>6} {'total':>6}")
print("  " + "-" * 50)

for name, img in sorted(tests.items()):
    h, w = img.shape
    n = h * w
    raw_bpp = len(bc(img.tobytes())) * 8 / n

    # Best reorder only
    best_reord = float('inf')
    for rname, rfunc in REORDERS:
        sig, _ = rfunc(img)
        sz = len(bc(bytes(sig)))
        if sz < best_reord: best_reord = sz
    reord_bpp = best_reord * 8 / n

    # Best reorder + predict
    best_pred = float('inf')
    for rname, rfunc in REORDERS:
        sig, sh = rfunc(img)
        for pname, pfunc in PREDICTORS:
            res = pfunc(sig, sh)
            sz = len(bc(bytes(res)))
            if sz < best_pred: best_pred = sz
    pred_bpp = best_pred * 8 / n

    # Best full composition
    best_full, _ = try_all_compositions(img)
    full_bpp = best_full * 8 / n

    print(f"  {name:<12} {raw_bpp:>6.2f} {reord_bpp:>6.2f} {pred_bpp:>6.2f} {full_bpp:>6.2f} {full_bpp:>6.2f}")
