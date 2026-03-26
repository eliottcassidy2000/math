#!/usr/bin/env python3
"""
Benchmark ALL TIC variants vs PNG on Kodak + real photo.
kind-pasteur-2026-03-26-S33
"""
import sys, io, os, time, subprocess, tempfile, glob
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

DIR = os.path.dirname(os.path.abspath(__file__))
CODECS = {
    "TIC-zlib": os.path.join(DIR, "tic_cli.exe"),
    "TIC-GR2":  os.path.join(DIR, "tic_gr2.exe"),
    "TIC-GR3":  os.path.join(DIR, "tic_gr3.exe"),
}

def png_bench(crop, repeats=10):
    buf = io.BytesIO()
    crop.save(buf, format='PNG', optimize=True, compress_level=9)
    sz = buf.tell()
    t0 = time.perf_counter()
    for _ in range(repeats):
        buf = io.BytesIO(); crop.save(buf, format='PNG', optimize=True, compress_level=9)
    return sz, (time.perf_counter() - t0) / repeats

def tic_bench(exe, rgb, w, h, iters=50):
    raw_path = os.path.join(tempfile.gettempdir(), "tic_bench_in.raw")
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

print("=" * 100)
print("  TIC CODEC COMPARISON: zlib vs GR2 vs GR3 vs PNG")
print("  kind-pasteur-2026-03-26-S33")
print("=" * 100)

files = sorted(glob.glob("../../test_images/kodim*.png"))
files.append("../../inbox/vlcsnap-2026-03-10-17h12m34s408.png")

print(f"\n  {'Image':<14} {'PNG':>7}", end="")
for name in CODECS: print(f" {name:>9}", end="")
print(f" {'best/PNG':>8} {'RT':>6}")
print("  " + "-" * 75)

totals = {'PNG': 0}
for name in CODECS: totals[name] = 0
all_ok = True

for fpath in files:
    fname = os.path.basename(fpath)[:12]
    pil = Image.open(fpath).convert('RGB')
    w0, h0 = pil.size
    crop = pil.crop((w0//2-128, h0//2-128, w0//2+128, h0//2+128))
    rgb = np.array(crop); h, w = rgb.shape[:2]

    ps, _ = png_bench(crop)
    totals['PNG'] += ps

    sizes = {}; oks = {}
    for cname, exe in CODECS.items():
        if os.path.exists(exe):
            sz, _, _, ok = tic_bench(exe, rgb, w, h)
            sizes[cname] = sz; oks[cname] = ok
            totals[cname] += sz
            if not ok: all_ok = False
        else:
            sizes[cname] = 0; oks[cname] = False

    best = min(sizes.values())
    ratio = ps / best if best > 0 else 0
    rt = "OK" if all(oks.values()) else "FAIL"

    print(f"  {fname:<14} {ps:>7}", end="")
    for cname in CODECS: print(f" {sizes[cname]:>9}", end="")
    print(f" {ratio:>7.3f}x {rt:>6}")

# Aggregates
print(f"\n  {'AGGREGATE':<14} {totals['PNG']:>7}", end="")
for name in CODECS: print(f" {totals[name]:>9}", end="")
best_total = min(totals[name] for name in CODECS if totals[name] > 0)
print(f" {totals['PNG']/best_total:>7.3f}x {'ALL OK' if all_ok else 'FAIL':>6}")

print(f"\n  vs PNG:")
for name in CODECS:
    if totals[name] > 0:
        r = totals['PNG'] / totals[name]
        print(f"    {name}: {r:.3f}x ({(r-1)*100:+.1f}%)")

print(f"\n  Roundtrip: {'ALL PASS' if all_ok else 'FAILURES DETECTED'}")
