#!/usr/bin/env python3
"""
Comprehensive benchmark on genuine raw camera data.
opus-2026-04-02-S3

Tests:
1. TIC vs PNG on raw-derived 8-bit RGB at multiple sizes
2. Raw Bayer mosaic compression analysis
3. Information-theoretic breakdown
4. Cross-channel correlation in raw vs JPEG-decoded data
5. Creative compression ideas for raw sensor data
"""
import sys, io, os, time, subprocess, zlib
import numpy as np
from PIL import Image
from collections import Counter

sys.stdout.reconfigure(line_buffering=True)

DIR = os.path.dirname(os.path.abspath(__file__))

def png_size(arr):
    buf = io.BytesIO()
    Image.fromarray(arr).save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()

def tic_size(arr):
    h, w = arr.shape[:2]
    arr.flatten().tofile('/tmp/tic_tmp.raw')
    r = subprocess.run([os.path.join(DIR, 'tic_prod'), 'bench', '/tmp/tic_tmp.raw', str(w), str(h), '1'],
                       capture_output=True, text=True, timeout=600)
    for line in r.stdout.split('\n'):
        if 'Compressed:' in line: return int(line.split()[1])
    return 0

def entropy_bpp(data):
    flat = data.flatten()
    counts = Counter(flat)
    total = len(flat)
    return -sum((c/total) * np.log2(c/total) for c in counts.values())

def med_predict(plane):
    """MED prediction, returns int32 residuals."""
    h, w = plane.shape
    res = np.zeros((h, w), dtype=np.int32)
    p = plane.astype(np.int32)
    for r in range(h):
        for c in range(w):
            a = p[r, c-1] if c > 0 else 0
            b = p[r-1, c] if r > 0 else 0
            cv = p[r-1, c-1] if (r > 0 and c > 0) else 0
            mn, mx = min(a, b), max(a, b)
            pred = mn if cv >= mx else (mx if cv <= mn else a + b - cv)
            res[r, c] = p[r, c] - pred
    return res

print("="*100)
print("  GENUINE RAW CAMERA DATA BENCHMARK")
print("  opus-2026-04-02-S3")
print("="*100)

# ================================================================
# Load raw files
# ================================================================
import rawpy

raws = {}
for name, path in [("CR2_EOS60D", "/Users/e/Documents/GitHub/math/inbox/sample1.cr2"),
                    ("DNG_EOS350D", "/Users/e/Documents/GitHub/math/inbox/processed/2026-04-02/new/sample1.dng")]:
    raw = rawpy.imread(path)
    rgb8 = raw.postprocess(output_bps=8, no_auto_bright=True, use_camera_wb=True)
    rgb16 = raw.postprocess(output_bps=16, no_auto_bright=True, use_camera_wb=True)
    bayer = raw.raw_image_visible.copy()
    raws[name] = {
        'rgb8': rgb8, 'rgb16': rgb16, 'bayer': bayer,
        'black': raw.black_level_per_channel, 'white': raw.white_level,
        'pattern': raw.raw_pattern,
        'h': rgb8.shape[0], 'w': rgb8.shape[1],
    }
    print(f"\n  {name}: {rgb8.shape[1]}×{rgb8.shape[0]}, Bayer {bayer.shape[1]}×{bayer.shape[0]}, "
          f"{bayer.dtype} [{bayer.min()}-{bayer.max()}]")

# Also load JPEG-decoded images for comparison
jpegs = {}
for name, path in [("JPEG_Canon", "/Users/e/Documents/GitHub/math/inbox/processed/2026-04-02/new/IMG_2811.jpg"),
                    ("JPEG_iPhone", "/Users/e/Documents/GitHub/math/inbox/processed/2026-04-02/new/IMG_8669.jpg")]:
    if os.path.exists(path):
        arr = np.array(Image.open(path).convert('RGB'))
        jpegs[name] = arr
        print(f"  {name}: {arr.shape[1]}×{arr.shape[0]}")

# ================================================================
# EXPERIMENT 1: TIC vs PNG at all sizes, RAW vs JPEG-derived
# ================================================================
print(f"\n{'='*100}")
print(f"  EXPERIMENT 1: TIC vs PNG — raw-derived vs JPEG-derived")
print(f"{'='*100}")

print(f"\n  {'Image':<20} {'Size':<12} {'PNG':>10} {'TIC':>10} {'TIC/PNG':>8} {'bpp':>6} {'Δ':>6}")
print(f"  {'-'*75}")

for name, info in raws.items():
    rgb = info['rgb8']
    h, w = rgb.shape[:2]
    for sz in [512, 1024, 0]:
        if sz == 0:
            crop = rgb
            label = f"{w}x{h}"
        else:
            if sz > min(w, h): continue
            cy, cx = h//2, w//2
            crop = rgb[cy-sz//2:cy+sz//2, cx-sz//2:cx+sz//2]
            label = f"{sz}x{sz}"
        ch, cw = crop.shape[:2]
        npix = cw * ch
        ps = png_size(crop)
        ts = tic_size(crop)
        ratio = ts / ps
        bpp = 8.0 * ts / npix
        win = "WIN" if ts < ps else "loss"
        print(f"  {name:<20} {label:<12} {ps:>10,} {ts:>10,} {ratio:>7.3f}x {bpp:>5.2f} {win:>6}")

for name, rgb in jpegs.items():
    h, w = rgb.shape[:2]
    for sz in [512, 1024, 0]:
        if sz == 0:
            crop = rgb
            label = f"{w}x{h}"
        else:
            if sz > min(w, h): continue
            cy, cx = h//2, w//2
            crop = rgb[cy-sz//2:cy+sz//2, cx-sz//2:cx+sz//2]
            label = f"{sz}x{sz}"
        ch, cw = crop.shape[:2]
        npix = cw * ch
        ps = png_size(crop)
        ts = tic_size(crop)
        ratio = ts / ps
        bpp = 8.0 * ts / npix
        win = "WIN" if ts < ps else "loss"
        print(f"  {name:<20} {label:<12} {ps:>10,} {ts:>10,} {ratio:>7.3f}x {bpp:>5.2f} {win:>6}")

# ================================================================
# EXPERIMENT 2: Cross-channel correlation — RAW vs JPEG
# ================================================================
print(f"\n{'='*100}")
print(f"  EXPERIMENT 2: Cross-channel correlation (RAW vs JPEG-decoded)")
print(f"{'='*100}")

for name, rgb in list(raws.items()):
    crop = rgb['rgb8']
    h, w = crop.shape[:2]
    if min(w,h) > 1024:
        cy,cx = h//2, w//2
        crop = crop[cy-512:cy+512, cx-512:cx+512]
    R, G, B = crop[:,:,0].flatten().astype(float), crop[:,:,1].flatten().astype(float), crop[:,:,2].flatten().astype(float)
    rg = np.corrcoef(R, G)[0,1]
    rb = np.corrcoef(R, B)[0,1]
    gb = np.corrcoef(G, B)[0,1]
    # Decorrelation benefit
    rg_plane = ((crop[:,:,0].astype(int) - crop[:,:,1].astype(int)) & 0xFF).astype(np.uint8)
    bg_plane = ((crop[:,:,2].astype(int) - crop[:,:,1].astype(int)) & 0xFF).astype(np.uint8)
    H_R = entropy_bpp(crop[:,:,0])
    H_G = entropy_bpp(crop[:,:,1])
    H_B = entropy_bpp(crop[:,:,2])
    H_RG = entropy_bpp(rg_plane)
    H_BG = entropy_bpp(bg_plane)
    print(f"\n  {name} (RAW-derived 8-bit):")
    print(f"    corr(R,G)={rg:.3f}  corr(R,B)={rb:.3f}  corr(G,B)={gb:.3f}")
    print(f"    H(R)={H_R:.2f}  H(G)={H_G:.2f}  H(B)={H_B:.2f}  → total {H_R+H_G+H_B:.2f} bpp")
    print(f"    H(R-G)={H_RG:.2f}  H(B-G)={H_BG:.2f}  → decorrelated {H_G+H_RG+H_BG:.2f} bpp")
    print(f"    Decorrelation saves: {(H_R+H_G+H_B)-(H_G+H_RG+H_BG):.2f} bpp")

for name, crop in jpegs.items():
    h, w = crop.shape[:2]
    if min(w,h) > 1024:
        cy,cx = h//2, w//2
        crop = crop[cy-512:cy+512, cx-512:cx+512]
    R, G, B = crop[:,:,0].flatten().astype(float), crop[:,:,1].flatten().astype(float), crop[:,:,2].flatten().astype(float)
    rg = np.corrcoef(R, G)[0,1]
    rb = np.corrcoef(R, B)[0,1]
    gb = np.corrcoef(G, B)[0,1]
    rg_plane = ((crop[:,:,0].astype(int) - crop[:,:,1].astype(int)) & 0xFF).astype(np.uint8)
    bg_plane = ((crop[:,:,2].astype(int) - crop[:,:,1].astype(int)) & 0xFF).astype(np.uint8)
    H_R = entropy_bpp(crop[:,:,0])
    H_G = entropy_bpp(crop[:,:,1])
    H_B = entropy_bpp(crop[:,:,2])
    H_RG = entropy_bpp(rg_plane)
    H_BG = entropy_bpp(bg_plane)
    print(f"\n  {name} (JPEG-decoded 8-bit):")
    print(f"    corr(R,G)={rg:.3f}  corr(R,B)={rb:.3f}  corr(G,B)={gb:.3f}")
    print(f"    H(R)={H_R:.2f}  H(G)={H_G:.2f}  H(B)={H_B:.2f}  → total {H_R+H_G+H_B:.2f} bpp")
    print(f"    H(R-G)={H_RG:.2f}  H(B-G)={H_BG:.2f}  → decorrelated {H_G+H_RG+H_BG:.2f} bpp")
    print(f"    Decorrelation saves: {(H_R+H_G+H_B)-(H_G+H_RG+H_BG):.2f} bpp")

# ================================================================
# EXPERIMENT 3: Bayer-level compression analysis
# ================================================================
print(f"\n{'='*100}")
print(f"  EXPERIMENT 3: Raw Bayer mosaic — compression frontier")
print(f"{'='*100}")

for name, info in raws.items():
    bayer = info['bayer']
    black = info['black']
    white = info['white']
    h, w = bayer.shape
    bits = int(np.ceil(np.log2(white + 1)))

    print(f"\n  {name}: {w}×{h} Bayer, {bits}-bit, black={black}, white={white}")

    # Crop for speed
    crop = bayer[h//2-256:h//2+256, w//2-256:w//2+256].astype(np.int32)
    crop = np.clip(crop - min(black), 0, white - min(black))
    ch, cw = crop.shape

    # Raw entropy
    H_raw = entropy_bpp(crop)
    print(f"    Pixel entropy (512x512 crop): {H_raw:.2f} bpp (vs {bits} raw)")

    # MED prediction on interleaved Bayer
    print(f"    Computing MED residuals...")
    res_full = med_predict(crop.astype(np.uint16) if bits <= 16 else crop)
    H_med = entropy_bpp(res_full)
    print(f"    MED residual entropy (interleaved): {H_med:.2f} bpp")

    # Per-channel MED
    total_H_perchan = 0
    for ri in range(2):
        for ci in range(2):
            ch_data = crop[ri::2, ci::2].astype(np.uint16)
            ch_res = med_predict(ch_data)
            H_ch = entropy_bpp(ch_res)
            ch_name = ['R','G1','G2','B'][ri*2+ci]
            frac = ch_data.size / crop.size
            total_H_perchan += H_ch * frac
            print(f"    MED on {ch_name} channel: {H_ch:.2f} bpp")
    print(f"    Per-channel MED total: {total_H_perchan:.2f} bpp")

    # Decorrelated Bayer (G-avg, R-G, B-G, G-diff)
    R_ch = crop[0::2, 0::2]; G1_ch = crop[0::2, 1::2]
    G2_ch = crop[1::2, 0::2]; B_ch = crop[1::2, 1::2]
    G_avg = ((G1_ch.astype(np.int64) + G2_ch.astype(np.int64)) // 2).astype(np.int32)
    G_diff = (G1_ch.astype(np.int32) - G2_ch.astype(np.int32))
    R_dec = R_ch.astype(np.int32) - G_avg
    B_dec = B_ch.astype(np.int32) - G_avg

    H_gavg = entropy_bpp(G_avg)
    H_gdiff = entropy_bpp(G_diff)
    H_rdec = entropy_bpp(R_dec)
    H_bdec = entropy_bpp(B_dec)
    # Each covers 1/4 of pixels
    total_dec = (H_gavg + H_gdiff + H_rdec + H_bdec) / 4 * 4  # same total pixels
    print(f"\n    Decorrelated Bayer:")
    print(f"      G_avg:  {H_gavg:.2f} bpp")
    print(f"      G_diff: {H_gdiff:.2f} bpp (green channel noise)")
    print(f"      R-G:    {H_rdec:.2f} bpp")
    print(f"      B-G:    {H_bdec:.2f} bpp")
    print(f"      Total (weighted): {(H_gavg + H_gdiff + H_rdec + H_bdec)/4:.2f} bpp per Bayer pixel")

    # With MED on each decorrelated channel
    print(f"    Decorrelated + MED:")
    for ch_name, ch_data in [("G_avg", G_avg), ("G_diff", G_diff), ("R-G", R_dec), ("B-G", B_dec)]:
        # Need uint-compatible data for MED
        shifted = (ch_data - ch_data.min()).astype(np.uint16)
        ch_res = med_predict(shifted)
        H = entropy_bpp(ch_res)
        print(f"      {ch_name}: {H:.2f} bpp after MED")

print(f"\n{'='*100}")
print(f"  DONE")
print(f"{'='*100}")
