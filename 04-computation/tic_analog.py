#!/usr/bin/env python3
"""
ANALOG THINKING: The digital-analog divide is only as real as you lack math.

An image is a sampled 2D continuous function f(x,y).
Audio compression treats sound as a sampled 1D function g(t).
The math doesn't care about dimension. Apply EVERY audio insight to images.

IDEAS:

1. FULL LPC ON RASTER SCAN: Scan the entire image as a 1D signal
   (raster order = unrolled to a vector). Apply LPC order-8 prediction.
   This captures inter-ROW and inter-PIXEL correlations in one model.
   Audio compressors (FLAC) use LPC order 1-12. Why not for images?

2. ADAPTIVE LPC PER BLOCK: Like FLAC's per-frame adaptation.
   Divide image into blocks (like audio frames of 4096 samples).
   Compute optimal LPC coefficients per block.
   Encode coefficients + residual.

3. SPECTRAL PREDICTION: Compute DFT of each row, predict DFT coefficients
   from previous row's DFT. Smooth images have smooth spectra → small
   spectral differences.

4. SINUSOIDAL MODEL: Decompose each row into a sum of sinusoids
   (like sinusoidal audio coding). For smooth images, very few sinusoids
   suffice. Encode frequencies + amplitudes + residual.

5. CONTINUOUS INTERPOLATION: Treat pixels as samples of a continuous
   function. Use SINC INTERPOLATION to predict pixel values from
   non-adjacent known samples. This uses the FULL bandwidth of the signal,
   not just nearest neighbors.

6. AUTOCORRELATION-OPTIMAL PREDICTOR: Compute the autocorrelation
   function of the image (or a local window). Design the optimal
   Wiener predictor from the autocorrelation. This is provably the
   best linear predictor for stationary signals.

7. Z-TRANSFORM VIEW: The MED predictor has a specific transfer function
   in the z-domain. Can we design a BETTER transfer function by
   analyzing the image's spectral properties?

8. DELTA-SIGMA MODULATION: Instead of predicting each pixel independently,
   use a feedback loop (like delta-sigma ADC). The integrator accumulates
   error, making the noise spectrum shaped (more error at high frequencies,
   less at low frequencies where it matters).

kind-pasteur-2026-03-25-S9
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

def m_med(plane):
    h, w = plane.shape; res = bytearray()
    for r in range(h):
        for c in range(w):
            a = int(plane[r,c-1]) if c>0 else 0
            b = int(plane[r-1,c]) if r>0 else 0
            c2 = int(plane[r-1,c-1]) if r>0 and c>0 else 0
            res.append((int(plane[r,c]) - med(a,b,c2)) & 0xFF)
    return bc(bytes(res))

# ============================================================
# IDEA 1: FULL LPC ON 1D RASTER SCAN (like FLAC for the whole image)
# ============================================================

def lpc_coeffs(signal, order):
    """Compute LPC coefficients using autocorrelation method (Levinson-Durbin)."""
    n = len(signal)
    if n < order + 1:
        return np.zeros(order)
    # Autocorrelation
    r = np.correlate(signal, signal, mode='full')
    r = r[n-1:n+order]  # r[0], r[1], ..., r[order]
    if r[0] < 1e-10:
        return np.zeros(order)
    # Levinson-Durbin
    a = np.zeros(order)
    e = r[0]
    for i in range(order):
        if e < 1e-10:
            break
        lam = r[i+1]
        for j in range(i):
            lam -= a[j] * r[i-j]
        k = lam / e
        a_new = np.zeros(order)
        a_new[i] = k
        for j in range(i):
            a_new[j] = a[j] - k * a[i-1-j]
        a = a_new
        e *= (1 - k*k)
    return a

def encode_lpc_1d(plane, order=4):
    """Flatten image to 1D (raster), apply LPC prediction."""
    h, w = plane.shape
    signal = plane.ravel().astype(float)
    n = len(signal)

    # Compute LPC coefficients on the whole signal
    coeffs = lpc_coeffs(signal - signal.mean(), order)

    # Predict
    residuals = np.zeros(n, dtype=np.uint8)
    for i in range(n):
        pred = 0.0
        for j in range(min(i, order)):
            pred += coeffs[j] * signal[i-1-j]
        pred = np.clip(round(pred), 0, 255)
        residuals[i] = (int(signal[i]) - int(pred)) & 0xFF

    # Encode: coefficients (as float32) + compressed residuals
    coeff_bytes = np.array(coeffs, dtype=np.float32).tobytes()
    total = len(coeff_bytes) + len(bc(bytes(residuals)))
    return total

# ============================================================
# IDEA 2: ADAPTIVE LPC PER BLOCK (FLAC-style frames)
# ============================================================

def encode_lpc_adaptive(plane, block_size=1024, order=4):
    """Split raster signal into blocks, compute LPC per block."""
    h, w = plane.shape
    signal = plane.ravel().astype(float)
    n = len(signal)

    all_coeffs = []
    all_residuals = bytearray()

    for start in range(0, n, block_size):
        end = min(start + block_size, n)
        block = signal[start:end]

        coeffs = lpc_coeffs(block - block.mean(), order)
        all_coeffs.extend(coeffs)

        for i in range(len(block)):
            pred = 0.0
            for j in range(min(i, order)):
                pred += coeffs[j] * block[i-1-j]
            idx = start + i
            # Also use pixels from previous block as context
            if i == 0 and start > 0:
                pred = signal[start - 1]
            pred = np.clip(round(pred), 0, 255)
            all_residuals.append((int(block[i]) - int(pred)) & 0xFF)

    coeff_bytes = np.array(all_coeffs, dtype=np.float32).tobytes()
    total = len(coeff_bytes) + len(bc(bytes(all_residuals)))
    return total

# ============================================================
# IDEA 3: 2D LPC (use BOTH horizontal and vertical neighbors)
# ============================================================

def encode_lpc_2d(plane, order=3):
    """2D LPC: predict pixel from 2D causal neighborhood.
    Neighborhood: (r,c-1), (r,c-2), (r-1,c), (r-1,c-1), (r-1,c+1), (r-2,c).
    This is the OPTIMAL linear predictor for the given autocorrelation."""
    h, w = plane.shape
    img = plane.astype(float)

    # Collect training data: each pixel and its 2D neighborhood
    # Neighborhood: a1=(r,c-1), a2=(r,c-2), b1=(r-1,c), b2=(r-1,c-1), b3=(r-1,c+1), b4=(r-2,c)
    n_feats = 6
    X = []
    Y = []
    for r in range(2, h):
        for c in range(2, w-1):
            feats = [img[r,c-1], img[r,c-2], img[r-1,c], img[r-1,c-1], img[r-1,c+1], img[r-2,c]]
            X.append(feats)
            Y.append(img[r,c])

    X = np.array(X)
    Y = np.array(Y)

    # Solve least squares for optimal 2D LPC coefficients
    try:
        coeffs, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
    except:
        coeffs = np.array([1, 0, 0, 0, 0, 0], dtype=float)  # fallback: just use left

    # Apply prediction
    residuals = bytearray()
    for r in range(h):
        for c in range(w):
            feats = [
                img[r,c-1] if c>0 else 0,
                img[r,c-2] if c>1 else 0,
                img[r-1,c] if r>0 else 0,
                img[r-1,c-1] if r>0 and c>0 else 0,
                img[r-1,c+1] if r>0 and c+1<w else 0,
                img[r-2,c] if r>1 else 0,
            ]
            pred = np.clip(round(sum(f*k for f,k in zip(feats, coeffs))), 0, 255)
            residuals.append((int(plane[r,c]) - pred) & 0xFF)

    coeff_bytes = np.array(coeffs, dtype=np.float32).tobytes()
    total = len(coeff_bytes) + len(bc(bytes(residuals)))
    return total

# ============================================================
# IDEA 4: AUTOCORRELATION-OPTIMAL PREDICTOR (Wiener filter)
# ============================================================

def encode_wiener(plane, radius=3):
    """Optimal linear predictor from local autocorrelation.
    For each pixel, predict from ALL pixels within radius r.
    Weights computed from the image's autocorrelation function."""
    h, w = plane.shape
    img = plane.astype(float)

    # Compute 2D autocorrelation from the image itself
    from scipy.signal import correlate2d
    acf = correlate2d(img - img.mean(), img - img.mean(), mode='full')
    # Normalize
    acf /= (acf.max() + 1e-10)

    # The causal neighborhood within radius
    cr, cc = acf.shape[0]//2, acf.shape[1]//2
    causal_offsets = []
    for dr in range(-radius, 1):
        for dc in range(-radius, radius+1):
            if dr == 0 and dc >= 0:
                continue  # not yet known
            if dr*dr + dc*dc <= radius*radius:
                causal_offsets.append((dr, dc))

    n_feats = len(causal_offsets)

    # Build autocorrelation matrix R and cross-correlation vector p
    R = np.zeros((n_feats, n_feats))
    p = np.zeros(n_feats)
    for i, (dr1, dc1) in enumerate(causal_offsets):
        p[i] = acf[cr + 0 - dr1, cc + 0 - dc1]  # r(target - neighbor_i)
        for j, (dr2, dc2) in enumerate(causal_offsets):
            R[i, j] = acf[cr + dr1 - dr2, cc + dc1 - dc2]

    # Solve Wiener-Hopf equation: R @ w = p
    try:
        weights = np.linalg.solve(R + 1e-8 * np.eye(n_feats), p)
    except:
        weights = np.zeros(n_feats)
        if causal_offsets:
            weights[0] = 1.0  # fallback: just use nearest

    # Apply predictor
    residuals = bytearray()
    for r in range(h):
        for c in range(w):
            pred = 0.0
            for k, (dr, dc) in enumerate(causal_offsets):
                nr, nc = r + dr, c + dc
                if 0 <= nr < h and 0 <= nc < w:
                    pred += weights[k] * img[nr, nc]
            pred = np.clip(round(pred), 0, 255)
            residuals.append((int(plane[r, c]) - pred) & 0xFF)

    weight_bytes = np.array(weights, dtype=np.float32).tobytes()
    total = len(weight_bytes) + len(bc(bytes(residuals)))
    return total

# ============================================================
# IDEA 5: IMAGE-AS-AUDIO (flatten to 1D, apply FLAC-style encoding)
# ============================================================

def encode_image_as_audio(plane):
    """Scan image as 1D signal with SNAKE scan (alternating row direction).
    Apply adaptive LPC (order varies per block to minimize residual).
    This is literally how FLAC would encode the image."""
    h, w = plane.shape

    # Snake scan: even rows L→R, odd rows R→L (better 1D locality)
    signal = []
    for r in range(h):
        if r % 2 == 0:
            signal.extend(plane[r])
        else:
            signal.extend(plane[r][::-1])
    signal = np.array(signal, dtype=float)

    # Adaptive LPC: try orders 1-8, pick best per block
    block_size = w * 2  # 2 rows per block (natural for snake scan)
    n = len(signal)
    all_residuals = bytearray()
    orders_used = bytearray()

    for start in range(0, n, block_size):
        end = min(start + block_size, n)
        block = signal[start:end]

        # Try multiple orders
        best_size = float('inf')
        best_res = None
        best_order = 0

        for order in [0, 1, 2, 4, 8]:
            if order == 0:
                # Raw
                res = block.astype(int) & 0xFF
            else:
                coeffs = lpc_coeffs(block - block.mean(), order)
                res = np.zeros(len(block), dtype=int)
                for i in range(len(block)):
                    pred = 0.0
                    for j in range(min(i, order)):
                        pred += coeffs[j] * block[i-1-j]
                    if i == 0 and start > 0:
                        pred = signal[start-1]
                    res[i] = (int(block[i]) - int(np.clip(round(pred), 0, 255))) & 0xFF

            trial = len(zlib.compress(bytes(res.astype(np.uint8)), 9))
            if trial < best_size:
                best_size = trial
                best_res = res.astype(np.uint8)
                best_order = order

        orders_used.append(best_order)
        all_residuals.extend(best_res)

    total = len(bc(bytes(orders_used))) + len(bc(bytes(all_residuals)))
    return total

# ============================================================
# IDEA 6: INTER-ROW DIFFERENCE SPECTRUM (predict in frequency domain)
# ============================================================

def encode_spectral_delta(plane):
    """Compute DCT of each row. Delta-encode DCT coefficients between rows.
    For smooth images, DCT coefficients change slowly → small deltas."""
    h, w = plane.shape
    from scipy.fft import dct, idct

    prev_dct = None
    all_deltas = bytearray()

    for r in range(h):
        row_dct = dct(plane[r].astype(float), type=2, norm='ortho')
        # Quantize to integer (lossless: scale and round)
        # Actually for LOSSLESS, we can't quantize. Store exact differences.
        # Compromise: store DCT coefficients as int16
        row_int = np.clip(np.round(row_dct), -32768, 32767).astype(np.int16)

        if prev_dct is not None:
            delta = ((row_int - prev_dct + 32768) & 0xFFFF).astype(np.uint16)
            all_deltas.extend(delta.tobytes())
        else:
            all_deltas.extend(row_int.astype(np.uint16).tobytes())

        prev_dct = row_int

    total = len(bc(bytes(all_deltas)))
    return total

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

print("=" * 100)
print("  ANALOG THINKING: Audio compression meets image compression")
print("  The digital-analog divide is only as real as you lack math skills")
print("  kind-pasteur-2026-03-25-S9")
print("=" * 100)

tests = make_tests()

METHODS = [
    ("MED std",       lambda p: len(m_med(p))),
    ("LPC-1D ord4",   lambda p: encode_lpc_1d(p, 4)),
    ("LPC-1D ord8",   lambda p: encode_lpc_1d(p, 8)),
    ("LPC adaptive",  lambda p: encode_lpc_adaptive(p, 1024, 4)),
    ("LPC-2D",        encode_lpc_2d),
    ("Wiener r=3",    lambda p: encode_wiener(p, 3)),
    ("FLAC-style",    encode_image_as_audio),
    ("Spectral dlt",  encode_spectral_delta),
]

print(f"\n  {'Image':<12} {'PNG':>6}", end="")
for mname, _ in METHODS:
    print(f" {mname:>12}", end="")
print(f"  {'BEST':>12}")
print("  " + "-" * (22 + 13*len(METHODS) + 16))

for name, img in sorted(tests.items()):
    ps = png_size(img)
    sizes = {}
    for mname, mfunc in METHODS:
        try:
            sizes[mname] = mfunc(img) + 10
        except Exception as e:
            sizes[mname] = 999999

    best = min(sizes, key=sizes.get)
    ratio = ps / sizes[best] if sizes[best] > 0 else 0

    print(f"  {name:<12} {ps:>6}", end="")
    for mname, _ in METHODS:
        marker = "*" if mname == best else " "
        print(f" {sizes[mname]:>11}{marker}", end="")
    print(f"  {best:>12} {ratio:.2f}x")

print(f"""
  KEY QUESTION: Does the optimal LINEAR predictor (Wiener/2D-LPC) beat MED?

  MED is a NONLINEAR predictor (it selects between left, above, and linear).
  The Wiener filter is the optimal LINEAR predictor for stationary signals.
  2D-LPC finds optimal linear coefficients from the image autocorrelation.

  If Wiener/2D-LPC beats MED: there's room for better linear predictors.
  If MED beats them: MED's nonlinearity is essential, linear methods plateau.
""")
