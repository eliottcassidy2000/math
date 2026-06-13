#!/usr/bin/env python3
"""
HIGH-ORDER LPC: How much of the 48% untapped spectral potential can we capture?

From S11: spectral flatness of LPC8 residuals on natural = 0.52 (halfway to optimal).
Theory: higher-order LPC captures more of the autocorrelation structure.
LPC order 8 on raster ≈ 2D predictor using 8 nearest causal neighbors.
LPC order 16 ≈ using 2 full rows of context.
LPC order 128 ≈ using entire previous row as explicit context.

Question: does flatness approach 1.0 as order increases?
If yes: we just need higher-order prediction.
If no: there's a NONLINEAR component that linear methods can't capture.

Also test: SNAKE SCAN + HIGH ORDER (continuity at row boundaries → better LPC).

kind-pasteur-2026-03-25-S15
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

def lpc_coeffs(sig, order):
    n = len(sig)
    if n < order + 1: return np.zeros(order)
    r = np.correlate(sig, sig, mode='full')
    r = r[n-1:n+order]
    if r[0] < 1e-10: return np.zeros(order)
    a = np.zeros(order); e = r[0]
    for i in range(order):
        if e < 1e-10: break
        lam = r[i+1]
        for j in range(i): lam -= a[j] * r[i-j]
        k = lam / e
        a_new = np.zeros(order); a_new[i] = k
        for j in range(i): a_new[j] = a[j] - k * a[i-1-j]
        a = a_new; e *= (1 - k*k)
    return a

def spectral_flatness(signal):
    f = signal.astype(float) - signal.mean()
    if np.std(f) < 1e-10: return 1.0
    spec = np.abs(np.fft.rfft(f))**2
    spec = spec[1:]  # remove DC
    spec = spec[spec > 0]
    if len(spec) == 0: return 1.0
    return np.exp(np.mean(np.log(spec + 1e-20))) / np.mean(spec)

def encode_lpc(plane, order, snake=False):
    """LPC with given order on raster (or snake) scan."""
    h, w = plane.shape
    if snake:
        signal = []
        for r in range(h):
            if r % 2 == 0: signal.extend(plane[r])
            else: signal.extend(plane[r][::-1])
        signal = np.array(signal, dtype=float)
    else:
        signal = plane.ravel().astype(float)

    n = len(signal)
    coeffs = lpc_coeffs(signal - signal.mean(), order)

    res = np.empty(n, dtype=np.uint8)
    for i in range(n):
        pred = sum(coeffs[j] * signal[i-1-j] for j in range(min(i, order)))
        res[i] = (int(signal[i]) - int(np.clip(round(pred), 0, 255))) & 0xFF

    # Spectral flatness of residuals
    sf = spectral_flatness(res)

    coeff_bytes = len(np.array(coeffs, dtype=np.float32).tobytes())
    compressed = len(bc(bytes(res)))

    return compressed + coeff_bytes, sf

# ============================================================
# TEST: Sweep LPC order from 1 to 256
# ============================================================

SZ = 128

def make_tests():
    T = {}; np.random.seed(42)
    x, y = np.meshgrid(np.arange(SZ, dtype=float), np.arange(SZ, dtype=float))
    r = np.sqrt((x-SZ/2)**2 + (y-SZ/2)**2)
    T["circles"] = (np.sin(r/5)*127+128).astype(np.uint8)
    T["blob"] = (np.exp(-r**2/(2*(SZ/4)**2))*255).astype(np.uint8)
    sm = np.random.randint(0, 256, (SZ//16, SZ//16), dtype=np.uint8)
    T["natural"] = np.clip(np.array(Image.fromarray(sm).resize((SZ, SZ), Image.BILINEAR)).astype(float)
                           + np.random.normal(0, 10, (SZ, SZ)), 0, 255).astype(np.uint8)
    return T

print("=" * 90)
print("  HIGH-ORDER LPC: Chasing the 48% untapped spectral potential")
print("  kind-pasteur-2026-03-25-S15")
print("=" * 90)

tests = make_tests()

for name, img in sorted(tests.items()):
    ps = png_size(img)
    print(f"\n  {name} (PNG={ps}):")
    print(f"  {'order':>6} {'raster_sz':>10} {'raster_sf':>10} {'snake_sz':>10} {'snake_sf':>10} {'vs_PNG':>7}")
    print("  " + "-" * 60)

    for order in [1, 2, 4, 8, 16, 32, 64, 128]:
        r_sz, r_sf = encode_lpc(img, order, snake=False)
        s_sz, s_sf = encode_lpc(img, order, snake=True)
        best = min(r_sz, s_sz) + 10
        ratio = ps / best if best > 0 else 0
        print(f"  {order:>6} {r_sz:>10} {r_sf:>10.4f} {s_sz:>10} {s_sf:>10.4f} {ratio:>7.2f}x")

print(f"""
  KEY QUESTIONS ANSWERED:
  1. Does spectral flatness increase with LPC order?
     → If yes, higher order = closer to optimal.
  2. Does compressed size decrease monotonically with order?
     → If no, there's an optimal order (coefficient overhead vs prediction gain).
  3. Does snake scan help LPC?
     → Snake eliminates the discontinuity at row boundaries.
  4. Can we reach spectral flatness > 0.8?
     → If yes, we're capturing 80%+ of available structure.
""")
