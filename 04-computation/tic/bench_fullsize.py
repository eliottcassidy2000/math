#!/usr/bin/env python3
"""
Full-size natural image benchmark: TIC-GR4 vs GR2 vs PNG.
Tests at multiple sizes (crop + full) to see how compression scales.

opus-2026-04-02-S2
"""
import sys, io, os, time, subprocess, tempfile
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

DIR = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(DIR, "..", "..")

CODECS = {}
for name, exe in [("GR2", "tic_gr2_mac"), ("GR4", "tic_gr4")]:
    path = os.path.join(DIR, exe)
    if os.path.exists(path):
        CODECS[name] = path

IMAGES = [
    os.path.join(BASE, "inbox/processed/2026-04-02/new/vlcsnap-2026-03-10-17h12m34s408.png"),
    os.path.join(BASE, "inbox/processed/2026-04-02/new/IMG_2811.jpg"),
]

CROP_SIZES = [256, 512, 1024, 2048, 0]  # 0 = full size

def png_size(arr):
    pil = Image.fromarray(arr)
    buf = io.BytesIO()
    pil.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()

def tic_bench(exe, rgb_flat, w, h):
    raw_path = os.path.join(tempfile.gettempdir(), "tic_fullbench.raw")
    rgb_flat.tofile(raw_path)
    try:
        r = subprocess.run([exe, "bench", raw_path, str(w), str(h), "1"],
                           capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        return 0, 0, 0, False
    sz = 0; enc_ms = 0; dec_ms = 0; ok = False
    for line in r.stdout.split('\n'):
        if 'Compressed:' in line:
            sz = int(line.split()[1])
        if 'Encode:' in line:
            enc_ms = float(line.split()[1])
        if 'Decode:' in line:
            dec_ms = float(line.split()[1])
        if 'PASS' in line:
            ok = True
    if not ok and r.stderr:
        print(f"    STDERR: {r.stderr.strip()[:200]}")
    try: os.remove(raw_path)
    except: pass
    return sz, enc_ms, dec_ms, ok

print("=" * 110)
print("  FULL-SIZE NATURAL IMAGE BENCHMARK")
print("  opus-2026-04-02-S2")
print("=" * 110)
print(f"\n  Codecs: {', '.join(CODECS.keys())}")

for img_path in IMAGES:
    if not os.path.exists(img_path):
        print(f"\n  SKIP: {img_path} not found")
        continue

    fname = os.path.basename(img_path)[:20]
    pil = Image.open(img_path).convert('RGB')
    W0, H0 = pil.size
    print(f"\n{'='*110}")
    print(f"  IMAGE: {fname} ({W0}×{H0})")
    print(f"{'='*110}")

    print(f"\n  {'Size':<12} {'Pixels':>10} {'Raw':>10} {'PNG':>10}", end="")
    for c in CODECS: print(f" {c+' sz':>10} {c+' enc':>8} {c+' dec':>8}", end="")
    print(f" {'best/PNG':>8} {'RT':>4}")
    print("  " + "-" * (60 + 26 * len(CODECS)))

    for crop in CROP_SIZES:
        if crop == 0:
            # Full size
            arr = np.array(pil)
            label = f"{W0}x{H0}"
        else:
            if crop > min(W0, H0):
                continue
            cx, cy = W0 // 2, H0 // 2
            half = crop // 2
            arr = np.array(pil.crop((cx - half, cy - half, cx + half, cy + half)))
            label = f"{crop}x{crop}"

        h, w = arr.shape[:2]
        raw_sz = h * w * 3
        rgb_flat = arr.flatten().astype(np.uint8)

        ps = png_size(arr)

        results = {}
        for cname, exe in CODECS.items():
            sz, enc_ms, dec_ms, ok = tic_bench(exe, rgb_flat, w, h)
            results[cname] = (sz, enc_ms, dec_ms, ok)

        best_sz = min(r[0] for r in results.values() if r[0] > 0) if any(r[0] > 0 for r in results.values()) else ps
        ratio = ps / best_sz if best_sz > 0 else 0
        all_ok = all(r[3] for r in results.values())

        print(f"  {label:<12} {w*h:>10,} {raw_sz:>10,} {ps:>10,}", end="")
        for cname in CODECS:
            sz, enc, dec, ok = results[cname]
            print(f" {sz:>10,} {enc:>7.0f}ms {dec:>7.0f}ms", end="")
        print(f" {ratio:>7.3f}x {'OK' if all_ok else 'FAIL':>4}")

    # Summary for this image
    print()
    if 0 in CROP_SIZES or CROP_SIZES[-1] == 0:
        arr_full = np.array(pil)
        h, w = arr_full.shape[:2]
        ps_full = png_size(arr_full)
        rgb_full = arr_full.flatten().astype(np.uint8)
        print(f"  FULL IMAGE DETAIL ({w}×{h}):")
        print(f"    PNG:  {ps_full:>12,} bytes  ({ps_full*8/(w*h):.2f} bpp)")
        for cname, exe in CODECS.items():
            sz, enc, dec, ok = tic_bench(exe, rgb_full, w, h)
            bpp = sz * 8 / (w * h) if sz > 0 else 0
            ratio = ps_full / sz if sz > 0 else 0
            print(f"    {cname}: {sz:>12,} bytes  ({bpp:.2f} bpp)  {ratio:.3f}x vs PNG  RT:{'PASS' if ok else 'FAIL'}  enc:{enc:.0f}ms  dec:{dec:.0f}ms")

print(f"\n{'='*110}")
print("  DONE")
print(f"{'='*110}")
