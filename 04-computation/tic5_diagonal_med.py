#!/usr/bin/env python3
"""
TIC v5: LINE-BASED MED prediction along diagonal scan lines.

THE KEY INSIGHT FROM EXPERIMENTS:
  - Per-pixel avg-neighbor: scan direction doesn't matter (adapts to any angle)
  - Raw zlib: beats prediction on highly structured patterns
  - LINE-BASED prediction: THIS is where structure alignment matters

  MED prediction uses "left" (a) and "above" (b) context. In raster scan,
  left = (r, c-1), above = (r-1, c). For DIAGONAL scan, we redefine:
  - left = previous pixel on same diagonal
  - above = corresponding pixel on adjacent diagonal

  This gives PERFECT context for 45° patterns: the "above" pixel is on
  the other side of the same edge, and "left" is along the edge.

IMPLEMENTATION:
  For each scan direction, define left/above in that direction's coordinates.
  Then apply MED prediction using those coordinates.

  Raster:    a = (r, c-1),     b = (r-1, c),     c = (r-1, c-1)
  Transpose: a = (r-1, c),     b = (r, c-1),     c = (r-1, c-1)
  Anti-diag: a = prev on same adiag, b = same position on prev adiag
  Diagonal:  a = prev on same diag, b = same position on prev diag

kind-pasteur-2026-03-25-S6
"""
import sys, io, struct, zlib, time, math
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

try: import brotli; HAS_BROTLI = True
except: HAS_BROTLI = False

# ============================================================
# SCAN GENERATORS (return list of lines, each line is list of (r,c))
# ============================================================

def lines_raster(h, w):
    return [[(r,c) for c in range(w)] for r in range(h)]

def lines_transpose(h, w):
    return [[(r,c) for r in range(h)] for c in range(w)]

def lines_antidiag(h, w):
    lines = []
    for d in range(h + w - 1):
        line = [(r, d-r) for r in range(max(0,d-w+1), min(h,d+1))]
        if line: lines.append(line)
    return lines

def lines_diagonal(h, w):
    lines = []
    for d in range(-(w-1), h):
        line = [(r, r-d) for r in range(max(0,d), min(h,d+w))]
        if line: lines.append(line)
    return lines

SCANS = [
    ("raster",    lines_raster),
    ("transpose", lines_transpose),
    ("antidiag",  lines_antidiag),
    ("diagonal",  lines_diagonal),
]

# ============================================================
# MED ON SCAN LINES: "left" = previous on line, "above" = same idx on prev line
# ============================================================

def med(a, b, c):
    mn, mx = min(a, b), max(a, b)
    if c >= mx: return mn
    if c <= mn: return mx
    return a + b - c

def paeth(a, b, c):
    p = a + b - c
    pa, pb, pc = abs(p-a), abs(p-b), abs(p-c)
    if pa <= pb and pa <= pc: return a
    if pb <= pc: return b
    return c

def encode_lines_med(plane, lines):
    """Encode plane along scan lines using MED prediction.
    a = previous pixel on current line.
    b = same-index pixel on previous line.
    c = previous-index pixel on previous line."""
    residuals = bytearray()
    prev_line_vals = None

    for line in lines:
        n = len(line)
        vals = np.array([int(plane[r,c]) for r,c in line])
        res = np.empty(n, dtype=np.uint8)

        for i in range(n):
            a = int(vals[i-1]) if i > 0 else 0
            b = int(prev_line_vals[i]) if prev_line_vals is not None and i < len(prev_line_vals) else 0
            c = int(prev_line_vals[i-1]) if prev_line_vals is not None and i > 0 and i-1 < len(prev_line_vals) else 0
            pred = med(a, b, c)
            res[i] = (int(vals[i]) - pred) & 0xFF

        residuals.extend(res)
        prev_line_vals = vals

    return bytes(residuals)

def decode_lines_med(residuals, h, w, lines):
    """Decode: inverse of encode_lines_med."""
    plane = np.zeros((h, w), dtype=np.uint8)
    pos = 0
    prev_line_vals = None

    for line in lines:
        n = len(line)
        res = np.frombuffer(residuals[pos:pos+n], dtype=np.uint8)
        pos += n
        vals = np.empty(n, dtype=np.int32)

        for i in range(n):
            a = int(vals[i-1]) if i > 0 else 0
            b = int(prev_line_vals[i]) if prev_line_vals is not None and i < len(prev_line_vals) else 0
            c = int(prev_line_vals[i-1]) if prev_line_vals is not None and i > 0 and i-1 < len(prev_line_vals) else 0
            pred = med(a, b, c)
            vals[i] = (int(res[i]) + pred) & 0xFF

        for i, (r, c2) in enumerate(line):
            plane[r, c2] = vals[i]
        prev_line_vals = vals

    return plane

# ============================================================
# TRY ALL 4 FILTERS (none, sub, up, MED) per line, pick best
# ============================================================

def encode_lines_adaptive(plane, lines):
    """Per-line adaptive filter selection: try none/sub/up/MED, pick best."""
    residuals = bytearray()
    prev_line_vals = None

    for line in lines:
        n = len(line)
        vals = np.array([int(plane[r,c]) for r,c in line])

        best_abs = float('inf')
        best_res = None
        best_fid = 0

        for fid in range(4):
            res = np.empty(n, dtype=np.uint8)
            for i in range(n):
                a = int(vals[i-1]) if i > 0 else 0
                b = int(prev_line_vals[i]) if prev_line_vals is not None and i < len(prev_line_vals) else 0
                c = int(prev_line_vals[i-1]) if prev_line_vals is not None and i > 0 and i-1 < len(prev_line_vals) else 0

                if fid == 0: pred = 0        # none
                elif fid == 1: pred = a       # sub (left in scan direction)
                elif fid == 2: pred = b       # up (same pos on prev line)
                else: pred = med(a, b, c)     # MED

                res[i] = (int(vals[i]) - pred) & 0xFF

            # Score: sum of signed absolute residuals
            signed = res.astype(np.int16); signed[signed > 128] -= 256
            score = int(np.sum(np.abs(signed)))
            if score < best_abs:
                best_abs = score; best_fid = fid; best_res = res.copy()

        residuals.append(best_fid)
        residuals.extend(best_res)
        prev_line_vals = vals

    return bytes(residuals)

def decode_lines_adaptive(residuals, h, w, lines):
    """Decode adaptive-filtered lines."""
    plane = np.zeros((h, w), dtype=np.uint8)
    pos = 0
    prev_line_vals = None

    for line in lines:
        n = len(line)
        fid = residuals[pos]; pos += 1
        res = np.frombuffer(residuals[pos:pos+n], dtype=np.uint8)
        pos += n
        vals = np.empty(n, dtype=np.int32)

        for i in range(n):
            a = int(vals[i-1]) if i > 0 else 0
            b = int(prev_line_vals[i]) if prev_line_vals is not None and i < len(prev_line_vals) else 0
            c = int(prev_line_vals[i-1]) if prev_line_vals is not None and i > 0 and i-1 < len(prev_line_vals) else 0

            if fid == 0: pred = 0
            elif fid == 1: pred = a
            elif fid == 2: pred = b
            else: pred = med(a, b, c)
            vals[i] = (int(res[i]) + pred) & 0xFF

        for i, (r, c2) in enumerate(line):
            plane[r, c2] = vals[i]
        prev_line_vals = vals

    return plane

# ============================================================
# COMPRESS + BEST BACKEND
# ============================================================

def compress(data):
    results = [zlib.compress(data, 9)]
    if HAS_BROTLI:
        try: results.append(brotli.compress(data, quality=11))
        except: pass
    return min(results, key=len)

def png_size(img):
    pil = Image.fromarray(img, 'L' if img.ndim == 2 else 'RGB')
    buf = io.BytesIO(); pil.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()

# ============================================================
# ENCODE IMAGE: try all scan directions × filter modes
# ============================================================

def encode_best(plane):
    """Try every scan × filter combo, also raw, pick smallest."""
    h, w = plane.shape
    candidates = []

    for sid, (sname, sfunc) in enumerate(SCANS):
        lines = sfunc(h, w)

        # MED only
        res_med = encode_lines_med(plane, lines)
        c_med = compress(res_med)
        candidates.append((sid, 'med', c_med, res_med))

        # Adaptive (per-line filter selection)
        res_adap = encode_lines_adaptive(plane, lines)
        c_adap = compress(res_adap)
        candidates.append((sid, 'adaptive', c_adap, res_adap))

    # Raw
    c_raw = compress(plane.tobytes())
    candidates.append((-1, 'raw', c_raw, plane.tobytes()))

    best = min(candidates, key=lambda x: len(x[2]))
    return best

# ============================================================
# TEST IMAGES
# ============================================================

def make_tests(sz):
    T = {}; np.random.seed(42)
    x, y = np.meshgrid(np.arange(sz, dtype=float), np.arange(sz, dtype=float))

    T["solid"] = np.full((sz,sz), 128, dtype=np.uint8)
    T["grad_h"] = np.tile(np.linspace(0,255,sz,dtype=np.uint8),(sz,1))
    T["grad_v"] = np.tile(np.linspace(0,255,sz,dtype=np.uint8).reshape(-1,1),(1,sz))
    T["grad_d"] = ((x+y)*255//(2*sz-2)).astype(np.uint8)

    for deg in [0, 30, 45, 60, 90, 135]:
        th = math.radians(deg)
        proj = x*math.cos(th) + y*math.sin(th)
        T[f"stripes_{deg}"] = ((proj/8).astype(int)%2*255).astype(np.uint8)

    r = np.sqrt((x-sz/2)**2+(y-sz/2)**2)
    T["circles"] = (np.sin(r/5)*127+128).astype(np.uint8)
    T["blob"] = (np.exp(-r**2/(2*(sz/4)**2))*255).astype(np.uint8)
    T["random"] = np.random.randint(0,256,(sz,sz),dtype=np.uint8)

    sm = np.random.randint(0,256,(max(sz//16,2),max(sz//16,2)),dtype=np.uint8)
    T["natural"] = np.clip(np.array(Image.fromarray(sm).resize((sz,sz),Image.BILINEAR)).astype(float)
                           +np.random.normal(0,10,(sz,sz)),0,255).astype(np.uint8)

    # Edges at specific angles
    for deg in [0, 45, 90]:
        th = math.radians(deg)
        proj = x*math.cos(th)+y*math.sin(th)
        T[f"edge_{deg}"] = (proj > sz*0.6).astype(np.uint8)*200+28

    return T

# ============================================================
# MAIN BENCHMARK
# ============================================================

print("=" * 80)
print("  TIC v5: LINE-BASED MED ON ALL SCAN DIRECTIONS")
print("  The right combination: fast line-based prediction + structure-aligned scan")
print("  kind-pasteur-2026-03-25-S6")
print("=" * 80)

for sz in [64, 128]:
    tests = make_tests(sz)
    print(f"\n--- {sz}x{sz} ---")
    print(f"  {'Image':<14} {'PNG':>6} {'TIC5':>6} {'ratio':>6} {'scan':>10} {'mode':>8} {'RT':>4}")
    print("  " + "-" * 62)

    W, L, n = 0, 0, 0
    for name, img in sorted(tests.items()):
        n += 1
        ps = png_size(img)
        sid, mode, cdata, raw_res = encode_best(img)
        ts = len(cdata) + 10  # header estimate

        sname = SCANS[sid][0] if sid >= 0 else "raw"

        # Roundtrip check
        if sid >= 0:
            if mode == 'med':
                dec = decode_lines_med(raw_res, sz, sz, SCANS[sid][1](sz, sz))
            else:
                dec = decode_lines_adaptive(raw_res, sz, sz, SCANS[sid][1](sz, sz))
        else:
            dec = np.frombuffer(raw_res, dtype=np.uint8).reshape(sz, sz)
        ok = np.array_equal(img, dec)

        r = ps / ts if ts > 0 else 0
        if r > 1.001: W += 1
        elif r < 0.999: L += 1

        print(f"  {name:<14} {ps:>6} {ts:>6} {r:>6.2f} {sname:>10} {mode:>8} {'OK' if ok else 'FAIL':>4}")

    print(f"\n  {W}W / {n-W-L}T / {L}L out of {n}")
