#!/usr/bin/env python3
"""
tic_real_benchmark.py — Test TIC vs PNG on REAL photographs + screenshots

Addresses the devil's advocate critique (S338) that our benchmarks are cherry-picked.
Tests on actual camera photos, screenshots, and crops at various sizes.

Author: opus-2026-03-25-S340
"""

import sys, io, time, os
import numpy as np
from PIL import Image

sys.path.insert(0, os.path.dirname(__file__))
from tic2 import encode, decode

sys.stdout.reconfigure(line_buffering=True)

def png_size(img_array, mode='L'):
    """Get optimized PNG size for an image array."""
    if img_array.ndim == 3:
        mode = 'RGB'
    pil = Image.fromarray(img_array, mode)
    buf = io.BytesIO()
    pil.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()

def test_image(name, img, verify=True):
    """Test TIC vs PNG on a single image. Returns (name, png, tic, ratio, ok)."""
    t0 = time.time()
    tic_data = encode(img)
    t_enc = time.time() - t0

    # Verify roundtrip
    ok = True
    if verify:
        recovered = decode(tic_data)
        ok = np.array_equal(img, recovered)

    png = png_size(img)
    tic = len(tic_data)
    ratio = png / tic if tic > 0 else 0
    status = "WIN" if ratio > 1.001 else ("TIE" if ratio > 0.999 else "LOSS")
    lossless = "OK" if ok else "FAIL"

    print(f"  {name:<40} {img.shape[0]:>4}x{img.shape[1]:<4} PNG={png:>8} TIC={tic:>8} {ratio:>6.3f}x {status:>4} {lossless:>4} {t_enc:>6.2f}s")
    return (name, png, tic, ratio, status, ok)

print("=" * 100)
print("  TIC vs PNG: REAL IMAGE BENCHMARK")
print("  opus-2026-03-25-S340")
print("=" * 100)

results = []

# ============================================================
# REAL PHOTOGRAPHS
# ============================================================
print("\n--- REAL PHOTOGRAPHS ---")
print(f"  {'Name':<40} {'Size':>10} {'PNG':>8} {'TIC':>8} {'Ratio':>7} {'':>4} {'':>4} {'Time':>7}")

real_images = {
    'tron.jpeg': '/home/e/tron.jpeg',
    'baby.jpeg': '/home/e/baby.jpeg',
    't.jpg': '/home/e/t.jpg',
}

for name, path in real_images.items():
    if os.path.exists(path):
        try:
            img = np.array(Image.open(path).convert('RGB'))
            results.append(test_image(name, img))
        except Exception as e:
            print(f"  {name}: {e}")

# ============================================================
# SCREENSHOTS (mixed text/graphics)
# ============================================================
print("\n--- SCREENSHOTS ---")

screenshots = []
ss_dir = '/home/e/Pictures/Screenshots/'
if os.path.isdir(ss_dir):
    for f in sorted(os.listdir(ss_dir))[:4]:
        if f.endswith('.png'):
            screenshots.append(os.path.join(ss_dir, f))

for path in screenshots:
    try:
        img = np.array(Image.open(path).convert('RGB'))
        name = os.path.basename(path)[:35]
        H, W = img.shape[:2]
        if H >= 512 and W >= 512:
            cy, cx = H // 2, W // 2
            crop = img[cy-256:cy+256, cx-256:cx+256]
            results.append(test_image(f"{name}_512crop", crop))
        if H >= 256 and W >= 256:
            cy, cx = H // 2, W // 2
            crop2 = img[cy-128:cy+128, cx-128:cx+128]
            results.append(test_image(f"{name}_256crop", crop2))
        elif H >= 64 and W >= 64:
            results.append(test_image(name, img))
    except Exception as e:
        print(f"  {os.path.basename(path)}: {e}")

# ============================================================
# CROPS AT DIFFERENT SIZES (from first photo)
# ============================================================
print("\n--- SIZE SCALING (tron.jpeg crops) ---")

try:
    base = np.array(Image.open('/home/e/tron.jpeg').convert('RGB'))
    for sz in [32, 64, 128, 225]:
        if sz <= min(base.shape[:2]):
            crop = base[:sz, :sz]
            results.append(test_image(f"tron_{sz}x{sz}", crop))
except Exception as e:
    print(f"  tron crops: {e}")

# ============================================================
# GRAYSCALE CONVERSIONS OF PHOTOS
# ============================================================
print("\n--- GRAYSCALE PHOTOS ---")

for name, path in real_images.items():
    if os.path.exists(path):
        try:
            img = np.array(Image.open(path).convert('L'))
            results.append(test_image(f"{name}_gray", img))
        except Exception as e:
            print(f"  {name}: {e}")

# ============================================================
# HARD SYNTHETIC (from devil's advocate)
# ============================================================
print("\n--- HARD SYNTHETIC (from devil's advocate S338) ---")

np.random.seed(42)
N = 256

# Pure noise
results.append(test_image("pure_noise_256",
    np.random.randint(0, 256, (N, N), dtype=np.uint8)))

# High-noise simulated photo
results.append(test_image("high_noise_photo_256",
    np.clip(np.array([[int(128+60*np.sin(r/5)*np.cos(c/4)+np.random.normal(0,40))
        for c in range(N)] for r in range(N)]), 0, 255).astype(np.uint8)))

# JPEG artifact simulation
results.append(test_image("jpeg_artifacts_256",
    np.array([[(r//8*8 + c//8*8) % 256 + np.random.randint(-3, 4)
        for c in range(N)] for r in range(N)], dtype=np.uint8).clip(0, 255)))

# Random walk
results.append(test_image("random_walk_256",
    np.clip(np.cumsum(np.random.normal(0, 1, N*N)).reshape(N,N) * 3 + 128,
            0, 255).astype(np.uint8)))

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*100}")

wins = sum(1 for r in results if r[4] == "WIN")
ties = sum(1 for r in results if r[4] == "TIE")
losses = sum(1 for r in results if r[4] == "LOSS")
fails = sum(1 for r in results if not r[5])
total = len(results)

total_png = sum(r[1] for r in results)
total_tic = sum(r[2] for r in results)

print(f"  {wins} WINS / {ties} TIES / {losses} LOSSES out of {total}")
print(f"  Win rate: {wins/total*100:.0f}%")
print(f"  Lossless: {total - fails}/{total}")
print(f"  Aggregate: TIC is {total_png/total_tic:.3f}x vs PNG")

if losses > 0:
    print(f"\n  LOSSES:")
    for n, p, t, r, s, ok in sorted(results, key=lambda x: x[3]):
        if s == "LOSS":
            print(f"    {n:<40} PNG={p:>8} TIC={t:>8} ({r:.4f}x)")

if wins > 0:
    print(f"\n  TOP 5 WINS:")
    for n, p, t, r, s, ok in sorted(results, key=lambda x: -x[3])[:5]:
        if s == "WIN":
            print(f"    {n:<40} {r:.2f}x better ({p-t} bytes saved)")

print(f"\nDone. opus-2026-03-25-S340")
