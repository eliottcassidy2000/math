#!/usr/bin/env python3
"""
Investigate:
1. Why GR4 loses to GR2 on real photos
2. Whether recursive compression helps
3. What the per-plane breakdown looks like

opus-2026-04-02-S2
"""
import sys, io, os, time, subprocess, tempfile
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)
DIR = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(DIR, "..", "..")

IMG = os.path.join(BASE, "inbox/processed/2026-04-02/new/IMG_2811.jpg")
pil = Image.open(IMG).convert('RGB')

# Use a 1024x1024 crop for fast iteration
cx, cy = pil.size[0]//2, pil.size[1]//2
crop = pil.crop((cx-512, cy-512, cx+512, cy+512))
arr = np.array(crop)
h, w = arr.shape[:2]
print(f"Test image: {w}x{h} crop of IMG_2811.jpg")

# ============================================================
# EXPERIMENT 1: Per-plane analysis — which plane is GR4 losing on?
# ============================================================
print("\n" + "="*80)
print("  EXPERIMENT 1: Per-plane entropy analysis")
print("="*80)

# MED prediction for a plane
def med_predict(plane):
    h, w = plane.shape
    residuals = np.zeros_like(plane, dtype=np.int16)
    for r in range(h):
        for c in range(w):
            a = int(plane[r, c-1]) if c > 0 else 0
            b = int(plane[r-1, c]) if r > 0 else 0
            cv = int(plane[r-1, c-1]) if (r > 0 and c > 0) else 0
            mn, mx = min(a, b), max(a, b)
            if cv >= mx: pred = mn
            elif cv <= mn: pred = mx
            else: pred = a + b - cv
            residuals[r, c] = int(plane[r, c]) - pred
    return residuals

R, G, B = arr[:,:,0], arr[:,:,1], arr[:,:,2]

# G-RG-BG transform
grd_g = G.copy()
grd_rg = ((R.astype(int) - G.astype(int)) & 0xFF).astype(np.uint8)
grd_bg = ((B.astype(int) - G.astype(int)) & 0xFF).astype(np.uint8)

# YCoCg-R transform (mod 256)
co = ((R.astype(int) - B.astype(int)) & 0xFF).astype(np.uint8)
t = ((B.astype(int) + (co.astype(int) >> 1)) & 0xFF).astype(np.uint8)
cg = ((G.astype(int) - t.astype(int)) & 0xFF).astype(np.uint8)
y = ((t.astype(int) + (cg.astype(int) >> 1)) & 0xFF).astype(np.uint8)

print(f"\n  Color transform comparison:")
print(f"  {'Plane':<10} {'GRD range':>12} {'YCoCg range':>12} {'GRD mean|diff|':>16} {'YCoCg mean|diff|':>18}")

for name, grd_p, ycc_p in [("Luma", grd_g, y), ("Chroma1", grd_rg, co), ("Chroma2", grd_bg, cg)]:
    grd_res = med_predict(grd_p)
    ycc_res = med_predict(ycc_p)

    grd_wrapped = grd_res.copy()
    grd_wrapped[grd_wrapped > 128] -= 256
    grd_wrapped[grd_wrapped < -128] += 256
    ycc_wrapped = ycc_res.copy()
    ycc_wrapped[ycc_wrapped > 128] -= 256
    ycc_wrapped[ycc_wrapped < -128] += 256

    grd_mad = np.mean(np.abs(grd_wrapped))
    ycc_mad = np.mean(np.abs(ycc_wrapped))

    print(f"  {name:<10} [{grd_wrapped.min():>4},{grd_wrapped.max():>4}] [{ycc_wrapped.min():>4},{ycc_wrapped.max():>4}] {grd_mad:>16.2f} {ycc_mad:>18.2f}")

# ============================================================
# EXPERIMENT 2: What is the theoretical entropy?
# ============================================================
print("\n" + "="*80)
print("  EXPERIMENT 2: Residual entropy (bits per pixel)")
print("="*80)

def empirical_entropy_bpp(residuals):
    """Compute empirical entropy of wrapped residuals in bits per pixel."""
    wrapped = residuals.copy().flatten()
    wrapped[wrapped > 128] -= 256
    wrapped[wrapped < -128] += 256
    # Zigzag
    zz = np.where(wrapped >= 0, 2*wrapped, 2*(-wrapped)-1)
    # Histogram
    vals, counts = np.unique(zz, return_counts=True)
    probs = counts / counts.sum()
    H = -np.sum(probs * np.log2(probs))
    return H

for name, grd_p, ycc_p in [("Luma", grd_g, y), ("Chroma1", grd_rg, co), ("Chroma2", grd_bg, cg)]:
    grd_H = empirical_entropy_bpp(med_predict(grd_p))
    ycc_H = empirical_entropy_bpp(med_predict(ycc_p))
    print(f"  {name:<10}  GRD: {grd_H:.3f} bpp   YCoCg: {ycc_H:.3f} bpp   diff: {ycc_H-grd_H:+.3f}")

grd_total = sum(empirical_entropy_bpp(med_predict(p)) for p in [grd_g, grd_rg, grd_bg])
ycc_total = sum(empirical_entropy_bpp(med_predict(p)) for p in [y, co, cg])
print(f"  {'TOTAL':<10}  GRD: {grd_total:.3f} bpp   YCoCg: {ycc_total:.3f} bpp   diff: {ycc_total-grd_total:+.3f}")

# ============================================================
# EXPERIMENT 3: Recursive compression
# ============================================================
print("\n" + "="*80)
print("  EXPERIMENT 3: Recursive compression")
print("="*80)

# Idea: after MED prediction, the residual plane itself has structure.
# Can we predict the residuals from THEIR neighbors and compress the
# "residuals of residuals"?

plane = grd_g  # Use luma
res1 = med_predict(plane)
res1_wrapped = res1.copy()
res1_wrapped[res1_wrapped > 128] -= 256
res1_wrapped[res1_wrapped < -128] += 256
res1_uint8 = (res1_wrapped & 0xFF).astype(np.uint8)

# Second round: predict residuals from their neighbors
res2 = med_predict(res1_uint8)
res2_wrapped = res2.copy()
res2_wrapped[res2_wrapped > 128] -= 256
res2_wrapped[res2_wrapped < -128] += 256

H1 = empirical_entropy_bpp(res1)
H2 = empirical_entropy_bpp(res2)

print(f"\n  Luma plane ({w}x{h}):")
print(f"    Raw entropy:                            8.000 bpp")
print(f"    After MED (1st round residuals):        {H1:.3f} bpp")
print(f"    After 2nd MED (residuals of residuals): {H2:.3f} bpp")
print(f"    Gain from 2nd round:                    {H1-H2:+.3f} bpp")

# Third round
res2_uint8 = (res2_wrapped & 0xFF).astype(np.uint8)
res3 = med_predict(res2_uint8)
H3 = empirical_entropy_bpp(res3)
print(f"    After 3rd MED (res of res of res):      {H3:.3f} bpp")
print(f"    Gain from 3rd round:                    {H2-H3:+.3f} bpp")

# What about a different predictor on the residuals?
# Residuals tend to be ~0 with Laplacian distribution.
# A simple "predict 0" (no predictor) might be best for residuals:
zero_pred_H = empirical_entropy_bpp(res1)
print(f"\n    'Predict 0' on raw residuals:           {zero_pred_H:.3f} bpp  (same as 1st round, obviously)")

# Multi-resolution: encode at half-res, use it to predict full-res
print(f"\n  Multi-resolution approach:")
half = plane[::2, ::2]  # Downsample 2x
half_res = med_predict(half)
half_H = empirical_entropy_bpp(half_res)
print(f"    Half-res ({half.shape[1]}x{half.shape[0]}) entropy: {half_H:.3f} bpp")

# Predict full-res from upsampled half-res
upsampled = np.repeat(np.repeat(half, 2, axis=0), 2, axis=1)[:h, :w]
detail = (plane.astype(int) - upsampled.astype(int))
detail_wrapped = detail.copy()
detail_wrapped[detail_wrapped > 128] -= 256
detail_wrapped[detail_wrapped < -128] += 256
detail_H = empirical_entropy_bpp(detail.astype(np.int16))
print(f"    Detail residual entropy:                {detail_H:.3f} bpp")
# Total: half_H * (1/4 of pixels) + detail_H * (3/4 of pixels)
multi_total = half_H * 0.25 + detail_H * 0.75
print(f"    Multi-res total: {half_H:.3f}/4 + {detail_H:.3f}*3/4 = {multi_total:.3f} bpp")
print(f"    vs single-pass MED:                     {H1:.3f} bpp")
print(f"    Gain from multi-res:                    {H1-multi_total:+.3f} bpp")

# ============================================================
# EXPERIMENT 4: Sign correction value
# ============================================================
print("\n" + "="*80)
print("  EXPERIMENT 4: Sign correction impact")
print("="*80)

# Sign correction: when dominant gradient is negative, flip residual sign.
# How much does this help on real data?
res_flat = res1_wrapped.flatten()
pos_count = np.sum(res_flat > 0)
neg_count = np.sum(res_flat < 0)
zero_count = np.sum(res_flat == 0)
print(f"\n  Luma residual sign distribution:")
print(f"    Positive: {pos_count:>8,} ({100*pos_count/len(res_flat):.1f}%)")
print(f"    Zero:     {zero_count:>8,} ({100*zero_count/len(res_flat):.1f}%)")
print(f"    Negative: {neg_count:>8,} ({100*neg_count/len(res_flat):.1f}%)")
print(f"    Skewness: {np.mean(res_flat**3) / (np.std(res_flat)**3 + 1e-10):.3f}")

# After sign correction, what's the distribution?
print(f"\n  Zigzag encoding efficiency:")
zz = np.where(res_flat >= 0, 2*res_flat, 2*(-res_flat)-1)
print(f"    Mean zigzag value:    {np.mean(zz):.2f}")
print(f"    Median zigzag value:  {np.median(zz):.1f}")
mean_bits_k0 = np.mean(zz) + 1  # Unary coding with k=0
mean_bits_k2 = np.mean(zz >> 2) + 1 + 2  # With k=2
mean_bits_k4 = np.mean(zz >> 4) + 1 + 4  # With k=4
print(f"    Mean bits at k=0:     {mean_bits_k0:.2f}")
print(f"    Mean bits at k=2:     {mean_bits_k2:.2f}")
print(f"    Mean bits at k=4:     {mean_bits_k4:.2f}")
print(f"    Optimal (entropy):    {H1:.2f}")

# ============================================================
# EXPERIMENT 5: What about using zlib on residuals?
# ============================================================
print("\n" + "="*80)
print("  EXPERIMENT 5: Backend comparison on MED residuals")
print("="*80)

import zlib

for name, p in [("Luma(G)", grd_g), ("Chroma(R-G)", grd_rg), ("Chroma(B-G)", grd_bg)]:
    res = med_predict(p)
    wrapped = res.copy()
    wrapped[wrapped > 128] -= 256
    wrapped[wrapped < -128] += 256
    res_bytes = wrapped.astype(np.int8).tobytes()

    raw_sz = len(res_bytes)
    zlib_sz = len(zlib.compress(res_bytes, 9))

    H = empirical_entropy_bpp(res)
    entropy_bound = int(H * w * h / 8)

    print(f"  {name:<14} raw:{raw_sz:>8,}  zlib-9:{zlib_sz:>8,}  entropy:{entropy_bound:>8,}  "
          f"zlib/entropy:{zlib_sz/entropy_bound:.3f}x  zlib bpp:{8*zlib_sz/(w*h):.3f}")

print("\n  (zlib/entropy > 1.0 means zlib wastes bits beyond the theoretical minimum)")
print("  (GR typically matches or slightly exceeds zlib for Laplacian residuals)")
