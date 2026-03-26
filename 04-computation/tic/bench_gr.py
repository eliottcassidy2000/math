#!/usr/bin/env python3
"""
Benchmark TIC-GR (Golomb-Rice) vs TIC-zlib vs PNG on all 9 images.
kind-pasteur-2026-03-26-S32
"""
import sys, io, os, time, subprocess, tempfile, glob
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

TIC_ZLIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tic_cli.exe")
TIC_GR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tic_gr2.exe")

def png_bench(crop, repeats=10):
    buf = io.BytesIO()
    crop.save(buf, format='PNG', optimize=True, compress_level=9)
    sz = buf.tell()
    t0 = time.perf_counter()
    for _ in range(repeats):
        buf = io.BytesIO(); crop.save(buf, format='PNG', optimize=True, compress_level=9)
    return sz, (time.perf_counter() - t0) / repeats

def tic_bench(exe, rgb, w, h, iters=50):
    raw_path = os.path.join(tempfile.gettempdir(), "tic_gr_in.raw")
    rgb.tofile(raw_path)
    r = subprocess.run([exe, "bench", raw_path, str(w), str(h), str(iters)],
                       capture_output=True, text=True, timeout=120)
    sz = enc_ms = dec_ms = 0; ok = False
    for line in r.stdout.split('\n'):
        if 'Compressed:' in line: sz = int(line.split()[1])
        if 'Encode:' in line: enc_ms = float(line.split()[1])
        if 'Decode:' in line: dec_ms = float(line.split()[1])
        if 'PASS' in line: ok = True
    try: os.remove(raw_path)
    except: pass
    return sz, enc_ms/1000, dec_ms/1000, ok

print("=" * 95)
print("  TIC-GR vs TIC-zlib vs PNG: Kodak + Real Photo")
print("  kind-pasteur-2026-03-26-S32")
print("=" * 95)

files = sorted(glob.glob("../../test_images/kodim*.png"))
files.append("../../inbox/vlcsnap-2026-03-10-17h12m34s408.png")

print(f"\n  {'Image':<16} {'PNG':>7} {'TIC-z':>7} {'TIC-GR':>7} {'z/PNG':>6} {'GR/PNG':>6} "
      f"{'GR/z':>6} {'GR_enc':>7} {'GR_dec':>7} {'RT':>4}")
print("  " + "-" * 95)

tp = tz = tg = 0

for fpath in files:
    name = os.path.basename(fpath)[:14]
    pil = Image.open(fpath).convert('RGB')
    w0, h0 = pil.size
    crop = pil.crop((w0//2-128, h0//2-128, w0//2+128, h0//2+128))
    rgb = np.array(crop); h, w = rgb.shape[:2]

    ps, _ = png_bench(crop)
    zs, _, _, _ = tic_bench(TIC_ZLIB, rgb, w, h, iters=50)
    gs, ge, gd, gok = tic_bench(TIC_GR, rgb, w, h, iters=50)

    tp += ps; tz += zs; tg += gs
    rz = ps/zs if zs > 0 else 0
    rg = ps/gs if gs > 0 else 0
    rgz = zs/gs if gs > 0 else 0

    print(f"  {name:<16} {ps:>7} {zs:>7} {gs:>7} {rz:>5.3f}x {rg:>5.3f}x "
          f"{rgz:>5.3f}x {ge*1000:>6.1f}ms {gd*1000:>6.1f}ms {'OK' if gok else 'FAIL':>4}")

print(f"\n  {'AGGREGATE':<16} {tp:>7} {tz:>7} {tg:>7} {tp/tz:>5.3f}x {tp/tg:>5.3f}x {tz/tg:>5.3f}x")

print(f"""
  SUMMARY:
    TIC-zlib vs PNG: {tp/tz:.3f}x ({(tp/tz-1)*100:.1f}% better)
    TIC-GR vs PNG:   {tp/tg:.3f}x ({(tp/tg-1)*100:.1f}% better)
    TIC-GR vs zlib:  {tz/tg:.3f}x ({(tz/tg-1)*100:.1f}% better entropy coding)

  GR has NO ZLIB DEPENDENCY — pure C, no external libraries.
  GR encoding is typically faster than zlib-9 (no LZ77 search).
""")
