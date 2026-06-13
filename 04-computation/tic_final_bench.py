#!/usr/bin/env python3
"""
TIC FINAL BENCHMARK: 6-filter + adaptive color transform.

Tests on Kodak 8 + real photo. Compares against PNG.
Uses compiled C binary (tic_final.exe).

kind-pasteur-2026-03-25-S23
"""
import sys, io, os, time, struct, subprocess, tempfile, glob, math
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

TIC_EXE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tic_final.exe")

def png_benchmark(pil_img, repeats=10):
    buf = io.BytesIO()
    pil_img.save(buf, format='PNG', optimize=True, compress_level=9)
    sz = buf.tell()
    t0 = time.perf_counter()
    for _ in range(repeats):
        buf = io.BytesIO()
        pil_img.save(buf, format='PNG', optimize=True, compress_level=9)
    t = (time.perf_counter() - t0) / repeats
    return sz, t

def tic_benchmark(rgb_array, w, h, iters=50):
    raw_path = os.path.join(tempfile.gettempdir(), "tic_input.raw")
    rgb_array.tofile(raw_path)

    result = subprocess.run(
        [TIC_EXE, "bench", raw_path, str(w), str(h), str(iters)],
        capture_output=True, text=True, timeout=120
    )

    enc_ms = dec_ms = 0; comp_sz = 0; rt_ok = False; ct = "?"
    for line in result.stdout.strip().split('\n'):
        if 'Compressed:' in line: comp_sz = int(line.split()[1])
        if 'Encode:' in line: enc_ms = float(line.split()[1])
        if 'Decode:' in line: dec_ms = float(line.split()[1])
        if 'Roundtrip:' in line: rt_ok = 'PASS' in line
        if 'ct=' in line and 'TIC FINAL' in line:
            ct = line.split('ct=')[1].strip() if 'ct=' in line else "?"

    try: os.remove(raw_path)
    except: pass

    return comp_sz, enc_ms / 1000, dec_ms / 1000, rt_ok, ct

print("=" * 100)
print("  TIC FINAL: 6-filter adaptive + YCoCg-R/G-RG-BG selection — Kodak + Real Photo")
print("  kind-pasteur-2026-03-25-S23")
print("=" * 100)

test_files = sorted(glob.glob("test_images/kodim*.png"))
test_files.append("inbox/vlcsnap-2026-03-10-17h12m34s408.png")

print(f"\n  {'Image':<18} {'Raw':>8} {'PNG':>8} {'TIC':>8} {'PNG_r':>6} {'TIC_r':>6} "
      f"{'PNG_ms':>7} {'TIC_ms':>7} {'gain':>6} {'speed':>6} {'CT':>8} {'RT':>4}")
print("  " + "-" * 110)

totals = {'raw': 0, 'png': 0, 'tic': 0, 'png_t': 0, 'tic_t': 0}
all_pass = True; wins = 0

for fpath in test_files:
    name = os.path.basename(fpath)[:16]
    pil = Image.open(fpath).convert('RGB')

    # 256×256 center crop
    w0, h0 = pil.size
    cx, cy = w0//2, h0//2
    crop = pil.crop((cx-128, cy-128, cx+128, cy+128))
    rgb = np.array(crop)
    h, w = rgb.shape[:2]
    raw_sz = h * w * 3

    png_sz, png_t = png_benchmark(crop, repeats=10)
    tic_sz, tic_enc_t, tic_dec_t, rt_ok, ct = tic_benchmark(rgb, w, h, iters=50)
    if not rt_ok: all_pass = False

    png_r = raw_sz / png_sz
    tic_r = raw_sz / tic_sz if tic_sz > 0 else 0
    gain = png_sz / tic_sz if tic_sz > 0 else 0
    speed = png_t / tic_enc_t if tic_enc_t > 0 else 0

    if gain > 1.001: wins += 1

    totals['raw'] += raw_sz; totals['png'] += png_sz; totals['tic'] += tic_sz
    totals['png_t'] += png_t; totals['tic_t'] += tic_enc_t

    print(f"  {name:<18} {raw_sz:>8} {png_sz:>8} {tic_sz:>8} {png_r:>5.2f}x {tic_r:>5.2f}x "
          f"{png_t*1000:>6.1f}ms {tic_enc_t*1000:>5.1f}ms {gain:>5.3f}x {speed:>5.2f}x {ct:>8} {'OK' if rt_ok else 'FAIL':>4}")

n = len(test_files)
agg_gain = totals['png'] / totals['tic'] if totals['tic'] > 0 else 0
agg_speed = totals['png_t'] / totals['tic_t'] if totals['tic_t'] > 0 else 0
agg_png_r = totals['raw'] / totals['png']
agg_tic_r = totals['raw'] / totals['tic']

print(f"\n  {'AGGREGATE':<18} {totals['raw']:>8} {totals['png']:>8} {totals['tic']:>8} "
      f"{agg_png_r:>5.2f}x {agg_tic_r:>5.2f}x "
      f"{totals['png_t']*1000:>6.1f}ms {totals['tic_t']*1000:>5.1f}ms "
      f"{agg_gain:>5.3f}x {agg_speed:>5.2f}x {'':>8} {'ALL' if all_pass else '!!!'}")

# Weissman
W = agg_gain * (math.log(totals['png_t']) / math.log(totals['tic_t'])) if totals['tic_t'] > 0 and totals['png_t'] > 0 else agg_gain

print(f"""
  ╔══════════════════════════════════════════════════════════════╗
  ║  TIC FINAL SCORES                                            ║
  ║                                                              ║
  ║  Images:       {n} ({n-1} Kodak + 1 real photo, 256×256 crops)    ║
  ║  Win rate:     {wins}/{n} ({wins/n*100:.0f}%)                                         ║
  ║  Compression:  {agg_gain:.3f}x better than PNG ({(agg_gain-1)*100:+.1f}%)              ║
  ║  Speed:        {agg_speed:.2f}x vs PNG                                   ║
  ║  Weissman:     {W:.2f}                                            ║
  ║  Roundtrip:    {'ALL PASS' if all_pass else 'FAILURES'}                                      ║
  ║                                                              ║
  ║  Algorithm: color decorrelation + 6-filter adaptive + zlib-9 ║
  ║  Key: MED as 6th filter + adaptive YCoCg-R / G-RG-BG        ║
  ╚══════════════════════════════════════════════════════════════╝
""")
