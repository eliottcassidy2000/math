#!/usr/bin/env python3
"""
FULL C BENCHMARK: True Weissman scores on Kodak + real photo.

Uses the compiled tic_bench.exe for pure C encode/decode speed.
Compares against PIL's PNG (which uses libpng, also C).

kind-pasteur-2026-03-25-S22
"""
import sys, io, os, time, struct, subprocess, tempfile
import numpy as np
from PIL import Image
import glob

sys.stdout.reconfigure(line_buffering=True)

TIC_EXE = os.path.join(os.path.dirname(__file__), "tic_bench.exe")

def png_benchmark(pil_img, repeats=10):
    """Benchmark PNG encode: size and time."""
    buf = io.BytesIO()
    pil_img.save(buf, format='PNG', optimize=True, compress_level=9)
    sz = buf.tell()
    t0 = time.perf_counter()
    for _ in range(repeats):
        buf = io.BytesIO()
        pil_img.save(buf, format='PNG', optimize=True, compress_level=9)
    t = (time.perf_counter() - t0) / repeats
    return sz, t

def tic_benchmark(rgb_array, w, h, repeats=50):
    """Benchmark TIC encode using the C executable."""
    # Write raw RGB to temp file
    raw_path = os.path.join(tempfile.gettempdir(), "tic_bench_input.raw")
    tic_path = os.path.join(tempfile.gettempdir(), "tic_bench_output.tic")
    dec_path = os.path.join(tempfile.gettempdir(), "tic_bench_decoded.raw")

    rgb_array.tofile(raw_path)

    # Run benchmark
    result = subprocess.run(
        [TIC_EXE, "bench", raw_path, str(w), str(h), str(repeats)],
        capture_output=True, text=True, timeout=60
    )

    # Parse output
    lines = result.stdout.strip().split('\n')
    enc_ms = dec_ms = 0
    comp_sz = 0
    rt_ok = False

    for line in lines:
        if 'Compressed:' in line:
            comp_sz = int(line.split()[1])
        if 'Encode:' in line:
            enc_ms = float(line.split()[1])
        if 'Decode:' in line:
            dec_ms = float(line.split()[1])
        if 'Roundtrip:' in line:
            rt_ok = 'PASS' in line

    # Also do a single encode to get the actual compressed file
    subprocess.run([TIC_EXE, "encode", raw_path, str(w), str(h), tic_path],
                   capture_output=True, timeout=30)
    if os.path.exists(tic_path):
        comp_sz = os.path.getsize(tic_path)

    # Verify roundtrip
    subprocess.run([TIC_EXE, "decode", tic_path, str(w), str(h), dec_path],
                   capture_output=True, timeout=30)
    if os.path.exists(dec_path):
        decoded = np.fromfile(dec_path, dtype=np.uint8)
        rt_ok = np.array_equal(rgb_array.ravel(), decoded)

    # Clean up
    for p in [raw_path, tic_path, dec_path]:
        try: os.remove(p)
        except: pass

    return comp_sz, enc_ms / 1000, dec_ms / 1000, rt_ok

# ============================================================
# BENCHMARK ALL IMAGES
# ============================================================

print("=" * 95)
print("  TIC FULL C BENCHMARK: True Weissman Scores on Kodak + Real Photo")
print("  kind-pasteur-2026-03-25-S22")
print("=" * 95)

test_files = sorted(glob.glob("test_images/kodim*.png"))
test_files.append("inbox/vlcsnap-2026-03-10-17h12m34s408.png")

print(f"\n  {'Image':<18} {'Raw':>8} {'PNG':>8} {'TIC':>8} {'PNG_r':>6} {'TIC_r':>6} "
      f"{'PNG_ms':>7} {'TIC_ms':>7} {'size+':>6} {'speed':>6} {'W':>5} {'RT':>4}")
print("  " + "-" * 105)

total_raw = total_png = total_tic = 0
total_png_t = total_tic_t = 0
all_pass = True
weissman_sum = 0

import math

for fpath in test_files:
    name = os.path.basename(fpath)[:16]
    pil = Image.open(fpath).convert('RGB')

    # Full-size crop (256×256 center)
    w0, h0 = pil.size
    cx, cy = w0//2, h0//2
    crop = pil.crop((cx-128, cy-128, cx+128, cy+128))
    rgb = np.array(crop)
    h, w = rgb.shape[:2]
    raw_sz = h * w * 3

    # PNG benchmark
    png_sz, png_t = png_benchmark(crop, repeats=10)

    # TIC C benchmark
    tic_sz, tic_enc_t, tic_dec_t, rt_ok = tic_benchmark(rgb, w, h, repeats=100)
    if not rt_ok: all_pass = False

    png_r = raw_sz / png_sz
    tic_r = raw_sz / tic_sz if tic_sz > 0 else 0
    size_plus = png_sz / tic_sz if tic_sz > 0 else 0  # >1 = we're smaller
    speed = png_t / tic_enc_t if tic_enc_t > 0 else 0  # >1 = we're faster

    # Weissman score: W = r * log(T_ref) / log(T) (simplified)
    # Using compression ratio relative to raw, times speed relative
    if tic_enc_t > 0 and png_t > 0:
        W = size_plus * (math.log(png_t) / math.log(tic_enc_t)) if tic_enc_t != 1 else size_plus
    else:
        W = size_plus

    total_raw += raw_sz; total_png += png_sz; total_tic += tic_sz
    total_png_t += png_t; total_tic_t += tic_enc_t
    weissman_sum += W

    print(f"  {name:<18} {raw_sz:>8} {png_sz:>8} {tic_sz:>8} {png_r:>5.2f}x {tic_r:>5.2f}x "
          f"{png_t*1000:>6.1f}ms {tic_enc_t*1000:>5.1f}ms {size_plus:>5.3f}x {speed:>5.2f}x {W:>5.2f} {'OK' if rt_ok else 'FAIL':>4}")

n = len(test_files)
agg_size = total_png / total_tic if total_tic > 0 else 0
agg_speed = total_png_t / total_tic_t if total_tic_t > 0 else 0
agg_W = weissman_sum / n

print(f"\n  {'AGGREGATE':>18} {total_raw:>8} {total_png:>8} {total_tic:>8} "
      f"{total_raw/total_png:>5.2f}x {total_raw/total_tic:>5.2f}x "
      f"{total_png_t*1000:>6.1f}ms {total_tic_t*1000:>5.1f}ms "
      f"{agg_size:>5.3f}x {agg_speed:>5.2f}x {agg_W:>5.2f} {'ALL' if all_pass else '!!!'}")

print(f"""
  ╔════════════════════════════════════════════════════════════╗
  ║  FINAL SCORES (Full C Implementation)                      ║
  ║                                                            ║
  ║  Compression:  {agg_size:.1%} of PNG  ({(agg_size-1)*100:+.1f}% better)         ║
  ║  Speed:        {agg_speed:.2f}x vs PNG  ({'faster' if agg_speed > 1 else 'slower':>6})                  ║
  ║  Weissman:     {agg_W:.2f}  ({'BEATS PNG' if agg_W > 1 else 'loses to PNG'})                          ║
  ║  Roundtrip:    {'ALL PASS' if all_pass else 'FAILURES'}                                ║
  ║  Images:       {n} ({n-1} Kodak + 1 real photo)                     ║
  ╚════════════════════════════════════════════════════════════╝
""")
