#!/usr/bin/env python3
"""
Comprehensive lossless format comparison on ALL test images.
TIC vs PNG vs QOI (via qoi module or manual) vs WebP-lossless vs BMP vs entropy bound.

opus-2026-04-02-S8
"""
import sys, io, os, subprocess, time, struct, zlib, math
import numpy as np
from PIL import Image
from collections import Counter

sys.stdout.reconfigure(line_buffering=True)

DIR = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(DIR, '..', '..')
IMG_DIR = os.path.join(BASE, 'inbox/processed/2026-04-02/new/')
TIC = os.path.join(DIR, 'tic_adaptive')

# ================================================================
# FORMAT ENCODERS
# ================================================================

def png_encode(arr):
    buf = io.BytesIO()
    Image.fromarray(arr).save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.getvalue()

def webp_lossless(arr):
    buf = io.BytesIO()
    Image.fromarray(arr).save(buf, format='WEBP', lossless=True)
    return buf.getvalue()

def bmp_size(arr):
    h, w = arr.shape[:2]
    row_sz = ((w * 3 + 3) // 4) * 4
    return 54 + row_sz * h

def tic_encode(arr):
    """Encode with our adaptive TIC codec (C binary)."""
    h, w = arr.shape[:2]
    rgb = arr.flatten().astype(np.uint8)
    rgb.tofile('/tmp/tic_fmt.raw')
    r = subprocess.run([TIC, 'bench', '/tmp/tic_fmt.raw', str(w), str(h), '1'],
                       capture_output=True, text=True, timeout=600)
    sz = 0; ok = False
    for line in r.stdout.split('\n'):
        if 'Compressed:' in line: sz = int(line.split()[1])
        if 'PASS' in line: ok = True
    return sz, ok

def qoi_encode(arr):
    """QOI encoder in Python (matches the JS implementation)."""
    h, w = arr.shape[:2]
    out = bytearray()
    # Header
    out += b'qoif'
    out += struct.pack('>II', w, h)
    out += bytes([3, 0])  # channels=3, srgb
    idx = [0] * 64
    pr, pg, pb = 0, 0, 0
    run = 0
    for i in range(w * h):
        r, g, b = int(arr[i // w, i % w, 0]), int(arr[i // w, i % w, 1]), int(arr[i // w, i % w, 2])
        if r == pr and g == pg and b == pb:
            run += 1
            if run == 62 or i == w * h - 1:
                out.append(0xC0 | (run - 1)); run = 0
        else:
            if run > 0:
                out.append(0xC0 | (run - 1)); run = 0
            h_idx = (r * 3 + g * 5 + b * 7 + 255 * 11) & 63
            sr, sg, sb = (idx[h_idx] >> 16) & 0xFF, (idx[h_idx] >> 8) & 0xFF, idx[h_idx] & 0xFF
            if sr == r and sg == g and sb == b:
                out.append(h_idx)  # OP_INDEX
            else:
                idx[h_idx] = (r << 16) | (g << 8) | b
                dr, dg, db = (r - pr) & 0xFF, (g - pg) & 0xFF, (b - pb) & 0xFF
                if dr > 127: dr -= 256
                if dg > 127: dg -= 256
                if db > 127: db -= 256
                if -2 <= dr <= 1 and -2 <= dg <= 1 and -2 <= db <= 1:
                    out.append(0x40 | ((dr + 2) << 4) | ((dg + 2) << 2) | (db + 2))
                else:
                    drdg, dbdg = dr - dg, db - dg
                    if -32 <= dg <= 31 and -8 <= drdg <= 7 and -8 <= dbdg <= 7:
                        out.append(0x80 | (dg + 32))
                        out.append(((drdg + 8) << 4) | (dbdg + 8))
                    else:
                        out.append(0xFE); out.append(r); out.append(g); out.append(b)
            pr, pg, pb = r, g, b
    out += b'\x00' * 7 + b'\x01'
    return bytes(out)

def entropy_bound(arr):
    """Theoretical minimum: Shannon entropy on MED residuals per channel."""
    h, w = arr.shape[:2]
    n = w * h
    total_bits = 0
    def med(a, b, c):
        mn, mx = min(a, b), max(a, b)
        return mn if c >= mx else (mx if c <= mn else a + b - c)
    for ch in range(3):
        counts = Counter()
        plane = arr[:, :, ch].astype(int)
        for r in range(h):
            for c in range(w):
                a = plane[r, c-1] if c > 0 else 0
                b = plane[r-1, c] if r > 0 else 0
                cv = plane[r-1, c-1] if (r > 0 and c > 0) else 0
                res = plane[r, c] - med(a, b, cv)
                if res > 128: res -= 256
                elif res < -128: res += 256
                counts[res] += 1
        H = -sum((cnt/n) * math.log2(cnt/n) for cnt in counts.values())
        total_bits += H * n
    return int(math.ceil(total_bits / 8))

# ================================================================
# LOAD ALL IMAGES
# ================================================================

print("=" * 120)
print("  COMPREHENSIVE LOSSLESS FORMAT COMPARISON — ALL IMAGES")
print("  opus-2026-04-02-S8")
print("=" * 120)

imgs = {}
for f in sorted(os.listdir(IMG_DIR)):
    if f.startswith('.'): continue
    path = os.path.join(IMG_DIR, f)
    fsize = os.path.getsize(path)
    try:
        if f.endswith(('.cr2', '.dng')):
            import rawpy
            raw = rawpy.imread(path)
            arr = raw.postprocess(output_bps=8, no_auto_bright=True, use_camera_wb=True)
            imgs[f[:32]] = (arr, 'RAW', fsize)
        elif f.endswith('.gif'):
            arr = np.array(Image.open(path).convert('RGB'))
            imgs[f[:32]] = (arr, 'GIF', fsize)
        elif f.endswith('.mp4'):
            # Skip large videos, process small ones
            if fsize > 10_000_000:
                print(f"  SKIP {f}: {fsize/1e6:.0f}MB video (too large)")
                continue
            # Extract frame with ffmpeg or skip
            print(f"  SKIP {f}: video frame extraction not implemented in Python")
            continue
        else:
            arr = np.array(Image.open(path).convert('RGB'))
            imgs[f[:32]] = (arr, f.split('.')[-1].upper(), fsize)
    except Exception as e:
        print(f"  SKIP {f}: {e}")

print(f"\n  Loaded {len(imgs)} images\n")

# ================================================================
# BENCHMARK
# ================================================================

header = f"{'#':>2} {'Image':<34} {'Type':<4} {'WxH':>10} {'Raw':>10} {'Entropy':>10} {'TIC':>10} {'PNG':>10} {'QOI':>10} {'WebP-LL':>10} {'BMP':>10} {'Best':>6}"
print(header)
print("-" * len(header))

totals = {'raw': 0, 'entropy': 0, 'tic': 0, 'png': 0, 'qoi': 0, 'webpll': 0, 'bmp': 0}
tic_wins = 0
total_count = 0

for idx, (name, (arr, kind, fsize)) in enumerate(sorted(imgs.items()), 1):
    h, w = arr.shape[:2]
    npix = w * h
    raw_sz = npix * 3

    # Skip entropy computation for large images (slow in Python)
    if npix > 2_000_000:
        ent_sz = 0  # skip
    else:
        ent_sz = entropy_bound(arr)

    tic_sz, tic_ok = tic_encode(arr)
    png_data = png_encode(arr)
    png_sz = len(png_data)

    # QOI — skip for very large images
    if npix > 4_000_000:
        qoi_sz = 0
    else:
        qoi_data = qoi_encode(arr)
        qoi_sz = len(qoi_data)

    webp_data = webp_lossless(arr)
    webp_sz = len(webp_data)

    bmp_sz = bmp_size(arr)

    # Find best lossless (excluding entropy bound and BMP)
    candidates = {'TIC': tic_sz, 'PNG': png_sz}
    if qoi_sz > 0: candidates['QOI'] = qoi_sz
    candidates['WebP-LL'] = webp_sz
    best_name = min(candidates, key=candidates.get)
    if best_name == 'TIC': tic_wins += 1
    total_count += 1

    totals['raw'] += raw_sz
    totals['entropy'] += ent_sz
    totals['tic'] += tic_sz
    totals['png'] += png_sz
    totals['qoi'] += qoi_sz
    totals['webpll'] += webp_sz
    totals['bmp'] += bmp_sz

    ent_str = f"{ent_sz:>10,}" if ent_sz > 0 else f"{'(skip)':>10}"
    qoi_str = f"{qoi_sz:>10,}" if qoi_sz > 0 else f"{'(skip)':>10}"
    rt = "✓" if tic_ok else "✗"

    print(f"{idx:>2} {name:<34} {kind:<4} {w}x{h:>4} {raw_sz:>10,} {ent_str} {tic_sz:>10,}{rt} {png_sz:>10,} {qoi_str} {webp_sz:>10,} {bmp_sz:>10,} {best_name:>6}")

# ================================================================
# SUMMARY
# ================================================================

print(f"\n{'':>2} {'TOTALS':<34} {'':4} {'':>10} {totals['raw']:>10,} {totals['entropy']:>10,} {totals['tic']:>10,}  {totals['png']:>10,} {totals['qoi']:>10,} {totals['webpll']:>10,} {totals['bmp']:>10,}")

print(f"\n  === LOSSLESS RANKING (aggregate) ===")
agg = [('TIC', totals['tic']), ('PNG', totals['png']), ('WebP-LL', totals['webpll'])]
if totals['qoi'] > 0: agg.append(('QOI', totals['qoi']))
agg.sort(key=lambda x: x[1])
for i, (name, sz) in enumerate(agg, 1):
    vs_best = f"+{(sz/agg[0][1]-1)*100:.1f}%" if i > 1 else "BEST"
    print(f"    {i}. {name:<10} {sz:>12,} bytes   {vs_best}")

print(f"\n  TIC wins: {tic_wins}/{total_count} images ({100*tic_wins/total_count:.0f}%)")
print(f"\n  NOTE: WebP lossless via Pillow is TRUE lossless (unlike canvas API).")
print(f"  NOTE: Entropy bound computed on small images only (MED residuals, Shannon entropy).")
print(f"  NOTE: JPEG-LS / JPEG-XL not available in Python without extra native libraries.")
