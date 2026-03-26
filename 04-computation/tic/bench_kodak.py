#!/usr/bin/env python3
"""
Final Kodak benchmark for TIC v1.0.
Single canonical C implementation. Honest numbers. All corrections applied.

kind-pasteur-2026-03-25-S24
"""
import sys, io, os, time, subprocess, tempfile, glob
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

TIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tic_cli.exe")

def png_bench(crop, repeats=10):
    buf = io.BytesIO()
    crop.save(buf, format='PNG', optimize=True, compress_level=9)
    sz = buf.tell()
    t0 = time.perf_counter()
    for _ in range(repeats):
        buf = io.BytesIO()
        crop.save(buf, format='PNG', optimize=True, compress_level=9)
    return sz, (time.perf_counter() - t0) / repeats

def tic_bench(rgb, w, h, iters=50):
    raw_path = os.path.join(tempfile.gettempdir(), "tic_in.raw")
    rgb.tofile(raw_path)
    r = subprocess.run([TIC, "bench", raw_path, str(w), str(h), str(iters)],
                       capture_output=True, text=True, timeout=120)
    comp_sz = enc_ms = dec_ms = 0; ok = False
    for line in r.stdout.split('\n'):
        if 'Compressed:' in line: comp_sz = int(line.split()[1])
        if 'Encode:' in line: enc_ms = float(line.split()[1])
        if 'Decode:' in line: dec_ms = float(line.split()[1])
        if 'PASS' in line: ok = True
    try: os.remove(raw_path)
    except: pass
    return comp_sz, enc_ms/1000, dec_ms/1000, ok

print("=" * 88)
print("  TIC v1.0 — Final Kodak Benchmark (single C implementation, honest numbers)")
print("=" * 88)

files = sorted(glob.glob("../../test_images/kodim*.png"))
files.append("../../inbox/vlcsnap-2026-03-10-17h12m34s408.png")

print(f"\n  {'Image':<18} {'Raw':>8} {'PNG':>8} {'TIC':>8} {'gain':>6} "
      f"{'PNG_ms':>7} {'TIC_enc':>8} {'TIC_dec':>8} {'RT':>4}")
print("  " + "-" * 90)

tp = tt = tr = 0; tpt = ttt = 0; wins = 0; allok = True

for fpath in files:
    name = os.path.basename(fpath)[:16]
    pil = Image.open(fpath).convert('RGB')
    w0, h0 = pil.size
    crop = pil.crop((w0//2-128, h0//2-128, w0//2+128, h0//2+128))
    rgb = np.array(crop); h, w = rgb.shape[:2]; raw = h*w*3

    ps, pt = png_bench(crop)
    ts, te, td, ok = tic_bench(rgb, w, h)
    if not ok: allok = False
    g = ps/ts if ts > 0 else 0
    if g > 1.001: wins += 1

    tp += ps; tt += ts; tr += raw; tpt += pt; ttt += te

    print(f"  {name:<18} {raw:>8} {ps:>8} {ts:>8} {g:>5.3f}x "
          f"{pt*1000:>6.1f}ms {te*1000:>7.1f}ms {td*1000:>7.1f}ms {'OK' if ok else 'FAIL':>4}")

n = len(files)
ag = tp/tt if tt > 0 else 0

print(f"\n  {'AGGREGATE':<18} {tr:>8} {tp:>8} {tt:>8} {ag:>5.3f}x "
      f"{tpt*1000:>6.1f}ms {ttt*1000:>7.1f}ms")

print(f"""
  ╔═══════════════════════════════════════════════════════╗
  ║  TIC v1.0 FINAL RESULTS                               ║
  ║                                                        ║
  ║  Images:      {n} ({n-1} Kodak + 1 real photo)              ║
  ║  Wins:        {wins}/{n} ({wins*100//n}%)                                   ║
  ║  Compression: {ag:.3f}x better than PNG                    ║
  ║  Roundtrip:   {'ALL PASS' if allok else 'FAILURES'}                               ║
  ║                                                        ║
  ║  Algorithm:   G-RG-BG + MED + zlib-9 (single pass)    ║
  ║  Code:        tic.c + tic.h (one file, ~250 lines)    ║
  ║  Build:       gcc -O2 -o tic tic.c -lz                ║
  ╚═══════════════════════════════════════════════════════╝
""")
