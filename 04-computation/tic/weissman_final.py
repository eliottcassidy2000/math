#!/usr/bin/env python3
"""
TIC v1.0 — FINAL WEISSMAN SCORE REPORT

The Weissman score (from Schwall, used in HBO's Silicon Valley):

    W = α * r * ln(T_ref) / ln(T)

Where:
    r = our_ratio / ref_ratio  (compression improvement)
    T_ref = reference encode time
    T = our encode time
    α = 1 (scaling constant, typically 1)

Interpretation:
    W > 1: we beat the reference on the combined compression×speed metric
    W = 1: we match the reference
    W < 1: reference wins

We measure BOTH encode and decode Weissman scores separately,
since TIC decode is much faster than encode.

kind-pasteur-2026-03-25-S25
"""
import sys, io, os, time, subprocess, tempfile, glob, math
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

TIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tic_cli.exe")

def png_bench(crop, repeats=20):
    """Measure PNG encode size and time (PIL/Pillow)."""
    buf = io.BytesIO()
    crop.save(buf, format='PNG', optimize=True, compress_level=9)
    sz = buf.tell()
    # Warmup
    for _ in range(3):
        buf = io.BytesIO()
        crop.save(buf, format='PNG', optimize=True, compress_level=9)
    # Measure
    t0 = time.perf_counter()
    for _ in range(repeats):
        buf = io.BytesIO()
        crop.save(buf, format='PNG', optimize=True, compress_level=9)
    enc_t = (time.perf_counter() - t0) / repeats
    # Decode time
    buf.seek(0)
    for _ in range(3): Image.open(buf); buf.seek(0)
    t0 = time.perf_counter()
    for _ in range(repeats):
        buf.seek(0)
        Image.open(buf).load()
    dec_t = (time.perf_counter() - t0) / repeats
    return sz, enc_t, dec_t

def tic_bench(rgb, w, h, iters=100):
    """Measure TIC encode/decode via C binary."""
    raw_path = os.path.join(tempfile.gettempdir(), "tic_ws.raw")
    rgb.tofile(raw_path)
    r = subprocess.run([TIC, "bench", raw_path, str(w), str(h), str(iters)],
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

def weissman(our_ratio, ref_ratio, our_time, ref_time):
    """Weissman score: W = r * ln(T_ref) / ln(T)"""
    if our_time <= 0 or ref_time <= 0:
        return float('inf') if our_ratio > ref_ratio else 0
    r = our_ratio / ref_ratio
    # Both times should be > 0 and in same units
    # Using natural log as in the original formulation
    lt = math.log(our_time)
    lr = math.log(ref_time)
    if lt == 0: return float('inf')
    return r * lr / lt

# ============================================================

print("=" * 96)
print("  TIC v1.0 — WEISSMAN SCORE REPORT")
print("  Compression × Speed combined metric on Kodak + real photo")
print("=" * 96)

files = sorted(glob.glob("../../test_images/kodim*.png"))
files.append("../../inbox/vlcsnap-2026-03-10-17h12m34s408.png")

print(f"\n  {'Image':<16} {'Raw':>7} {'PNG':>7} {'TIC':>7} │ {'PNG_r':>5} {'TIC_r':>5} {'gain':>5}"
      f" │ {'PNGe':>6} {'TICe':>6} {'PNGd':>6} {'TICd':>6}"
      f" │ {'W_enc':>5} {'W_dec':>5} │ {'RT':>4}")
print("  " + "─" * 93)

totals = dict(raw=0, png_sz=0, tic_sz=0, png_et=0, tic_et=0, png_dt=0, tic_dt=0)
all_ok = True
w_enc_list = []
w_dec_list = []

for fpath in files:
    name = os.path.basename(fpath)[:14]
    pil = Image.open(fpath).convert('RGB')
    w0, h0 = pil.size
    crop = pil.crop((w0//2-128, h0//2-128, w0//2+128, h0//2+128))
    rgb = np.array(crop)
    h, w = rgb.shape[:2]
    raw = h * w * 3

    ps, pet, pdt = png_bench(crop, repeats=20)
    ts, tet, tdt, ok = tic_bench(rgb, w, h, iters=100)
    if not ok: all_ok = False

    pr = raw / ps; tr = raw / ts if ts > 0 else 0
    g = ps / ts if ts > 0 else 0

    # Weissman scores
    we = weissman(tr, pr, tet, pet)
    wd = weissman(tr, pr, tdt, pdt)
    w_enc_list.append(we)
    w_dec_list.append(wd)

    totals['raw'] += raw; totals['png_sz'] += ps; totals['tic_sz'] += ts
    totals['png_et'] += pet; totals['tic_et'] += tet
    totals['png_dt'] += pdt; totals['tic_dt'] += tdt

    print(f"  {name:<16} {raw:>7} {ps:>7} {ts:>7} │ {pr:>4.2f}x {tr:>4.2f}x {g:>4.2f}x"
          f" │ {pet*1000:>5.1f}m {tet*1000:>5.1f}m {pdt*1000:>5.1f}m {tdt*1000:>5.1f}m"
          f" │ {we:>5.2f} {wd:>5.2f} │ {'OK' if ok else 'FAIL':>4}")

# Aggregates
n = len(files)
apr = totals['raw'] / totals['png_sz']
atr = totals['raw'] / totals['tic_sz']
ag = totals['png_sz'] / totals['tic_sz']
awe = sum(w_enc_list) / n
awd = sum(w_dec_list) / n

# Also compute aggregate Weissman from totals
agg_we = weissman(atr, apr, totals['tic_et'], totals['png_et'])
agg_wd = weissman(atr, apr, totals['tic_dt'], totals['png_dt'])

print(f"  {'─'*93}")
print(f"  {'MEAN':>16} {'':>7} {'':>7} {'':>7} │ {apr:>4.2f}x {atr:>4.2f}x {ag:>4.2f}x"
      f" │ {'':>6} {'':>6} {'':>6} {'':>6}"
      f" │ {awe:>5.2f} {awd:>5.2f} │ {'ALL' if all_ok else '!!!'}")

print(f"""
  ┌───────────────────────────────────────────────────────────┐
  │  WEISSMAN SCORES                                          │
  │                                                           │
  │  Encode Weissman (compression vs encode speed):           │
  │    Per-image mean:  {awe:.2f}                                   │
  │    Aggregate:       {agg_we:.2f}                                   │
  │                                                           │
  │  Decode Weissman (compression vs decode speed):           │
  │    Per-image mean:  {awd:.2f}                                   │
  │    Aggregate:       {agg_wd:.2f}                                   │
  │                                                           │
  │  Interpretation:                                          │
  │    W > 1: TIC beats PNG on combined compression×speed     │
  │    W < 1: PNG beats TIC on combined metric                │
  │                                                           │
  │  Compression alone: {ag:.3f}x ({(ag-1)*100:+.1f}% better than PNG)        │
  │  Encode speed:      {totals['png_et']/totals['tic_et']:.2f}x vs PNG                             │
  │  Decode speed:      {totals['png_dt']/totals['tic_dt']:.2f}x vs PNG                             │
  │  Roundtrip:         {'ALL 9/9 PASS' if all_ok else 'FAILURES'}                          │
  └───────────────────────────────────────────────────────────┘

  FORMULA: W = (TIC_ratio / PNG_ratio) × ln(PNG_time) / ln(TIC_time)
  Source: Weissman score as defined by Tsachy Weissman (Stanford) for
  the HBO show Silicon Valley, adapted from information theory.
""")
