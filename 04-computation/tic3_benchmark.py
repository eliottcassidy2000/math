#!/usr/bin/env python3
"""
tic3_benchmark.py — TIC v3 vs TIC v2 vs PNG on real and synthetic images

Tests the structure-aligned scanning innovation.
Focuses on images where diagonal scan should help.

Author: opus-2026-03-25-S341
"""

import sys, io, time, os
import numpy as np
from PIL import Image

sys.path.insert(0, os.path.dirname(__file__))
from tic3 import encode as enc3, decode as dec3, MAGIC as M3
from tic2 import encode as enc2, decode as dec2
import struct

sys.stdout.reconfigure(line_buffering=True)

def png_size(img_array):
    mode = 'L' if img_array.ndim == 2 else 'RGB'
    pil = Image.fromarray(img_array, mode)
    buf = io.BytesIO()
    pil.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()

def get_scan(data):
    if data[:4] == b'TIC3':
        _, _, _, _, _, scan_id, _, _ = struct.unpack_from('<4sBHHBBBB', data, 0)
        if scan_id < 5:
            return ["raster","transp","diag","adiag","ring"][scan_id]
        if 32 <= scan_id < 36: return f"D/{['n','s','a','d'][scan_id-32]}"
        if 48 <= scan_id < 52: return f"A/{['n','s','a','d'][scan_id-48]}"
        return f"s{scan_id}"
    return "?"

def test(name, img):
    t0 = time.time()
    d3 = enc3(img)
    t3 = time.time() - t0
    d2 = enc2(img)
    png = png_size(img)
    s3, s2 = len(d3), len(d2)
    r3 = png / s3 if s3 > 0 else 0
    r2 = png / s2 if s2 > 0 else 0

    # Verify v3 roundtrip
    ok = np.array_equal(img, dec3(d3))
    scan = get_scan(d3)

    winner = "v3" if s3 < s2 else ("v2" if s2 < s3 else "tie")
    dim = f"{img.shape[0]}x{img.shape[1]}" + ("x3" if img.ndim == 3 else "")

    print(f"  {name:<30} {dim:>10} PNG={png:>7} v2={s2:>7} v3={s3:>7} v3/PNG={r3:.3f}x {winner:>3} {'OK' if ok else 'FAIL'} [{scan}]")
    return (name, png, s2, s3, r3, winner, ok, scan)

print("=" * 115)
print("  TIC v3 vs v2 vs PNG — Structure-Aligned Scanning Benchmark")
print("  opus-2026-03-25-S341")
print("=" * 115)

results = []

# ============================================================
# DIAGONAL-DOMINATED PATTERNS (where v3 should excel)
# ============================================================
print("\n--- DIAGONAL PATTERNS (v3 target) ---")
np.random.seed(42)
for N in [128, 256]:
    # 45-degree stripes
    results.append(test(f"diag_stripes_{N}",
        np.array([[(r+c)//8%2*255 for c in range(N)] for r in range(N)], dtype=np.uint8)))
    # 135-degree stripes
    results.append(test(f"antidiag_stripes_{N}",
        np.array([[(r-c)//8%2*255 for c in range(N)] for r in range(N)], dtype=np.uint8)))
    # Diagonal gradient
    results.append(test(f"diag_gradient_{N}",
        np.array([[int(255*(r+c)/(2*N-2)) for c in range(N)] for r in range(N)], dtype=np.uint8)))
    # Diagonal sine wave
    results.append(test(f"diag_sine_{N}",
        np.array([[int(128+127*np.sin((r+c)/10)) for c in range(N)] for r in range(N)], dtype=np.uint8)))
    # Diagonal edge
    results.append(test(f"diag_edge_{N}",
        np.array([[255 if r > c else 0 for c in range(N)] for r in range(N)], dtype=np.uint8)))
    # Checkerboard (has both diagonal directions)
    results.append(test(f"checker_{N}",
        np.array([[(r+c)%2*255 for c in range(N)] for r in range(N)], dtype=np.uint8)))

# ============================================================
# STANDARD PATTERNS (where raster should win)
# ============================================================
print("\n--- STANDARD PATTERNS ---")
for N in [128, 256]:
    results.append(test(f"h_gradient_{N}",
        np.tile(np.linspace(0, 255, N, dtype=np.uint8), (N, 1))))
    results.append(test(f"v_gradient_{N}",
        np.tile(np.linspace(0, 255, N, dtype=np.uint8).reshape(-1, 1), (1, N))))
    results.append(test(f"circle_{N}",
        np.array([[min(255, int(255*np.exp(-((r-N//2)**2+(c-N//2)**2)/(N**2/8)))) for c in range(N)] for r in range(N)], dtype=np.uint8)))
    results.append(test(f"solid_{N}",
        np.full((N, N), 128, dtype=np.uint8)))
    results.append(test(f"noise_{N}",
        np.random.randint(0, 256, (N, N), dtype=np.uint8)))

# ============================================================
# REAL PHOTOS
# ============================================================
print("\n--- REAL PHOTOS ---")
photo_files = {
    'tron.jpeg': '/home/e/tron.jpeg',
    'baby.jpeg': '/home/e/baby.jpeg',
    't.jpg': '/home/e/t.jpg',
}
for name, path in photo_files.items():
    if os.path.exists(path):
        img = np.array(Image.open(path).convert('RGB'))
        results.append(test(name, img))
        # Also grayscale
        gray = np.array(Image.open(path).convert('L'))
        results.append(test(f"{name}_gray", gray))

# Screenshots
ss_dir = '/home/e/Pictures/Screenshots/'
if os.path.isdir(ss_dir):
    for f in sorted(os.listdir(ss_dir))[:2]:
        if f.endswith('.png'):
            img = np.array(Image.open(os.path.join(ss_dir, f)).convert('RGB'))
            H, W = img.shape[:2]
            if H >= 256 and W >= 256:
                cy, cx = H // 2, W // 2
                crop = img[cy-128:cy+128, cx-128:cx+128]
                results.append(test(f"ss_{f[:20]}_256", crop))

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*115}")

v3_wins_vs_png = sum(1 for r in results if r[4] > 1.001)
v3_wins_vs_v2 = sum(1 for r in results if r[5] == "v3")
v2_wins_vs_v3 = sum(1 for r in results if r[5] == "v2")
ties = sum(1 for r in results if r[5] == "tie")
fails = sum(1 for r in results if not r[6])
total = len(results)

total_png = sum(r[1] for r in results)
total_v2 = sum(r[2] for r in results)
total_v3 = sum(r[3] for r in results)

print(f"  v3 vs PNG: {v3_wins_vs_png}/{total} wins ({v3_wins_vs_png/total*100:.0f}%)")
print(f"  v3 vs v2:  {v3_wins_vs_v2} wins / {v2_wins_vs_v3} losses / {ties} ties")
print(f"  Lossless:  {total - fails}/{total}")
print(f"  Aggregate: v3 = {total_png/total_v3:.3f}x vs PNG, v2 = {total_png/total_v2:.3f}x vs PNG")
print(f"  v3 saves {total_v2 - total_v3} bytes vs v2 total")

# Where v3 wins most over v2
v3_gains = [(r[0], r[2]-r[3], r[2]/r[3] if r[3]>0 else 0) for r in results if r[5] == "v3"]
if v3_gains:
    v3_gains.sort(key=lambda x: -x[2])
    print(f"\n  Top v3 gains over v2:")
    for name, saved, ratio in v3_gains[:5]:
        print(f"    {name:<30} saved {saved:>6}B ({ratio:.3f}x)")

# Where v2 wins over v3
v2_gains = [(r[0], r[3]-r[2], r[3]/r[2] if r[2]>0 else 0) for r in results if r[5] == "v2"]
if v2_gains:
    v2_gains.sort(key=lambda x: -x[2])
    print(f"\n  v2 still better on:")
    for name, extra, ratio in v2_gains[:5]:
        print(f"    {name:<30} v3 adds {extra:>6}B ({ratio:.3f}x)")

# Scan mode distribution
scan_dist = {}
for r in results:
    scan_dist[r[7]] = scan_dist.get(r[7], 0) + 1
print(f"\n  Scan mode distribution:")
for scan, count in sorted(scan_dist.items(), key=lambda x: -x[1]):
    print(f"    {scan:<15} {count:>3} images")

print(f"\nDone. opus-2026-03-25-S341")
