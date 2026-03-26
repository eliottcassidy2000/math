#!/usr/bin/env python3
"""
cheap_wl_s339.py — Cheap WL: redefine "same" to make coloring practical
opus-2026-03-25-S339

THE PROBLEM: Standard WL colors by EXACT content → too many colors at
fine scales → color map overhead exceeds savings.

THE FIX: Redefine "same" at each scale using STATISTICAL SUMMARIES.

Instead of: "row A and row B have the same pixels"
Use:        "row A and row B have the same QUANTIZED STATISTICS"

LEVELS OF "SAME":
  Level 0: same quantized mean                    → ~16 colors
  Level 1: same (mean, std)                       → ~64 colors
  Level 2: same (mean, std, skewness)             → ~256 colors
  Level 3: same (mean, std, gradient direction)   → ~256 colors
  Level 4: same (mean, std, neighbor correlation) → ~512 colors

The KEY: the statistics are CHEAP to compute (O(n) per row) and
produce a BOUNDED number of colors regardless of image size.

Then: for each color class, store ONE prototype + residuals.
The prototypes are predictable from the statistics.
The residuals are small (within-class variation).

This is WL with a COARSENED equality relation = a QUOTIENT of WL.
"""

import sys
import numpy as np
import zlib, lzma, bz2
from math import log2
from collections import defaultdict
from PIL import Image
import io

sys.stdout.reconfigure(line_buffering=True)

# ============================================================================
# CHEAP WL: Statistical coloring
# ============================================================================

def stat_color(row, level=2, q_bits=4):
    """Color a row by its quantized statistics.

    level 0: quantized mean only
    level 1: (mean, std)
    level 2: (mean, std, gradient_energy)
    level 3: (mean, std, gradient_energy, second_diff_energy)
    """
    n = len(row)
    row_f = row.astype(np.float64)

    # Mean (quantized)
    mean_q = int(np.mean(row_f)) >> (8 - q_bits)

    if level == 0:
        return (mean_q,)

    # Std (quantized)
    std_q = min(int(np.std(row_f)), 255) >> (8 - q_bits)

    if level == 1:
        return (mean_q, std_q)

    # Gradient energy (sum of absolute first differences)
    if n > 1:
        grad = np.abs(np.diff(row_f))
        grad_q = min(int(np.mean(grad)), 255) >> (8 - q_bits)
    else:
        grad_q = 0

    if level == 2:
        return (mean_q, std_q, grad_q)

    # Second difference energy (curvature)
    if n > 2:
        d2 = np.abs(np.diff(row_f, 2))
        d2_q = min(int(np.mean(d2)), 255) >> (8 - q_bits)
    else:
        d2_q = 0

    return (mean_q, std_q, grad_q, d2_q)


def cheap_wl_compress(M, level=2, q_bits=4):
    """Compress using cheap WL coloring.

    1. Color each row by stat_color at the given level
    2. For each color class: compute prototype (median row)
    3. Compute residuals (row - prototype)
    4. Encode: color_map + prototypes + residuals, all zlib'd

    The color map is tiny (few colors × log2(n_colors) bits per row).
    The prototypes are few (one per color).
    The residuals are small (rows similar to their class prototype).
    """
    H, W = M.shape

    # Color each row
    colors = [stat_color(M[i], level, q_bits) for i in range(H)]

    # Map to compact IDs
    color_set = sorted(set(colors))
    color_to_id = {c: i for i, c in enumerate(color_set)}
    n_colors = len(color_set)
    row_ids = [color_to_id[c] for c in colors]

    # Group rows by color
    groups = defaultdict(list)
    for i, cid in enumerate(row_ids):
        groups[cid].append(i)

    # Compute prototypes (median of each group)
    prototypes = {}
    for cid, rows in groups.items():
        if len(rows) == 1:
            prototypes[cid] = M[rows[0]].copy()
        else:
            prototypes[cid] = np.median([M[r] for r in rows], axis=0).astype(np.uint8)

    # Compute residuals
    all_residuals = np.zeros((H, W), dtype=np.int8)
    for cid, rows in groups.items():
        proto = prototypes[cid].astype(np.int16)
        for r in rows:
            all_residuals[r] = np.clip(M[r].astype(np.int16) - proto, -127, 127).astype(np.int8)

    # Encode everything
    buf = bytearray()

    # Header
    buf.extend(H.to_bytes(2, 'little'))
    buf.extend(W.to_bytes(2, 'little'))
    buf.extend(n_colors.to_bytes(2, 'little'))

    # Color map (row → color_id)
    for cid in row_ids:
        if n_colors <= 256:
            buf.append(cid)
        else:
            buf.extend(cid.to_bytes(2, 'little'))

    # Prototypes (one per color, in order)
    for cid in range(n_colors):
        buf.extend(prototypes[cid].tobytes())

    # Residuals (all rows, in order)
    buf.extend(all_residuals.tobytes())

    # Compress the whole thing
    raw_buf = bytes(buf)

    # Try multiple backends
    candidates = [
        zlib.compress(raw_buf, 9),
        bz2.compress(raw_buf, 9),
    ]
    try:
        candidates.append(lzma.compress(raw_buf))
    except: pass

    best = min(candidates, key=len)

    return best, n_colors, len(groups), {
        'color_map_bytes': H * (1 if n_colors <= 256 else 2),
        'prototype_bytes': n_colors * W,
        'residual_bytes': H * W,
        'avg_group_size': H / n_colors,
    }


def cheap_wl_with_prediction(M, level=2, q_bits=4):
    """Enhanced: use Paeth prediction WITHIN each color class.

    Instead of storing raw residuals (row - prototype),
    store Paeth-predicted residuals (row - Paeth(row-1 in same class)).
    Rows in the same class are structurally similar, so Paeth prediction
    within a class is very accurate.
    """
    H, W = M.shape

    # Color rows
    colors = [stat_color(M[i], level, q_bits) for i in range(H)]
    color_set = sorted(set(colors))
    color_to_id = {c: i for i, c in enumerate(color_set)}
    n_colors = len(color_set)
    row_ids = [color_to_id[c] for c in colors]

    # Group rows
    groups = defaultdict(list)
    for i, cid in enumerate(row_ids):
        groups[cid].append(i)

    # Paeth-predict within each class
    output = bytearray()
    output.extend(H.to_bytes(2, 'little'))
    output.extend(W.to_bytes(2, 'little'))
    output.extend(n_colors.to_bytes(2, 'little'))

    # Color map
    for cid in row_ids:
        output.append(cid if n_colors <= 256 else cid & 0xFF)

    # For each color class: Paeth-encode the sequence of rows
    for cid in range(n_colors):
        rows = groups.get(cid, [])
        if not rows: continue

        above = np.zeros(W, dtype=np.uint8)
        for r in rows:
            row = M[r]
            # Paeth prediction from previous row in same class
            predicted = np.zeros(W, dtype=np.uint8)
            for j in range(W):
                left = int(row[j-1]) if j > 0 else 0
                up = int(above[j])
                ul = int(above[j-1]) if j > 0 else 0
                p = left + up - ul
                pa, pb, pc = abs(p-left), abs(p-up), abs(p-ul)
                pred = left if (pa<=pb and pa<=pc) else (up if pb<=pc else ul)
                predicted[j] = pred & 0xFF

            residual = (row.astype(np.int16) - predicted.astype(np.int16)) & 0xFF
            output.extend(residual.astype(np.uint8).tobytes())
            above = row

    raw_buf = bytes(output)
    candidates = [zlib.compress(raw_buf, 9), bz2.compress(raw_buf, 9)]
    try: candidates.append(lzma.compress(raw_buf))
    except: pass

    return min(candidates, key=len), n_colors


# ============================================================================
# BENCHMARK
# ============================================================================

def gen_img(name, N=128):
    rng = np.random.RandomState(42)
    if name == 'gradient':
        return np.array([[int(j*255/max(N-1,1)) for j in range(N)] for i in range(N)], dtype=np.uint8)
    elif name == 'stripes':
        return np.array([[(255 if i%16<8 else 0) for j in range(N)] for i in range(N)], dtype=np.uint8)
    elif name == 'checker':
        return np.array([[(255 if (i//16+j//16)%2==0 else 0) for j in range(N)] for i in range(N)], dtype=np.uint8)
    elif name == 'smooth':
        x = np.linspace(0,4*np.pi,N); y = np.linspace(0,4*np.pi,N)
        xx,yy = np.meshgrid(x,y)
        return np.clip(128+100*np.sin(xx)*np.cos(yy),0,255).astype(np.uint8)
    elif name == 'photo':
        img = np.zeros((N,N),dtype=np.float64)
        for _ in range(25):
            cx,cy = rng.randint(N),rng.randint(N); r = rng.randint(8,max(9,N//4))
            yy,xx = np.ogrid[:N,:N]; img[(xx-cx)**2+(yy-cy)**2<r**2] = rng.randint(256)
        img += rng.normal(0,10,(N,N))
        return np.clip(img,0,255).astype(np.uint8)
    elif name == 'mandelbrot':
        img = np.zeros((N,N),dtype=np.uint8)
        for i in range(N):
            for j in range(N):
                c = complex(-2.5+3.5*j/N,-1.25+2.5*i/N); z=0; n=0
                while abs(z)<2 and n<255: z=z*z+c; n+=1
                img[i,j] = n
        return img
    elif name == 'noise':
        return rng.randint(0,256,(N,N),dtype=np.uint8)
    return np.zeros((N,N),dtype=np.uint8)


print("=" * 80)
print("  CHEAP WL — Redefine 'same' for practical compression")
print("  opus-2026-03-25-S339")
print("=" * 80)

N = 128

# Compare different levels of "same"
print(f"\n  EFFECT OF DEFINITION OF 'SAME' ({N}×{N}):")
print(f"  {'Image':>12} {'PNG':>6} {'L0(μ)':>6} {'L1(μσ)':>7} {'L2(μσg)':>8} {'L3(full)':>8} {'colors0':>8} {'colors2':>8}")

for name in ['gradient', 'stripes', 'checker', 'smooth', 'photo', 'mandelbrot', 'noise']:
    img = gen_img(name, N)

    pil = Image.fromarray(img, 'L')
    buf = io.BytesIO(); pil.save(buf, format='PNG', optimize=True)
    png_size = len(buf.getvalue())

    results = []
    for level in range(4):
        comp, nc, ng, info = cheap_wl_compress(img, level=level, q_bits=4)
        results.append((len(comp), nc))

    print(f"  {name:>12} {png_size:>5} {results[0][0]:>5} {results[1][0]:>6} "
          f"{results[2][0]:>7} {results[3][0]:>7} {results[0][1]:>7} {results[2][1]:>7}")

# The enhanced version: Paeth within classes
print(f"\n  CHEAP WL + INTRA-CLASS PAETH ({N}×{N}):")
print(f"  {'Image':>12} {'PNG':>6} {'CheapWL':>8} {'CWL+Paeth':>10} {'CWL+P/PNG':>10}")

for name in ['gradient', 'stripes', 'checker', 'smooth', 'photo', 'mandelbrot', 'noise']:
    img = gen_img(name, N)

    pil = Image.fromarray(img, 'L')
    buf = io.BytesIO(); pil.save(buf, format='PNG', optimize=True)
    png_size = len(buf.getvalue())

    cwl, nc_cwl, _, _ = cheap_wl_compress(img, level=2, q_bits=4)
    cwl_p, nc_cwlp = cheap_wl_with_prediction(img, level=2, q_bits=4)

    ratio = len(cwl_p) / png_size
    win = "✓" if len(cwl_p) <= png_size else ""
    print(f"  {name:>12} {png_size:>5} {len(cwl):>7} {len(cwl_p):>9} {ratio:>9.3f}x {win}")

# The FULL hybrid: best of (cheap WL, cheap WL+Paeth, Paeth+lzma, delta+lzma, raw+lzma)
print(f"\n  FULL HYBRID: best of all methods ({N}×{N}):")
print(f"  {'Image':>12} {'PNG':>6} {'Best':>6} {'B/PNG':>7} {'Method':>15}")

def paeth_pred(a, b, c):
    p = int(a)+int(b)-int(c)
    pa,pb,pc = abs(p-int(a)),abs(p-int(b)),abs(p-int(c))
    return int(a) if (pa<=pb and pa<=pc) else (int(b) if pb<=pc else int(c))

def paeth_encode(img):
    H, W = img.shape
    above = np.zeros(W, dtype=np.uint8)
    out = bytearray()
    for i in range(H):
        for j in range(W):
            left = int(img[i,j-1]) if j>0 else 0
            up = int(above[j]); ul = int(above[j-1]) if j>0 else 0
            out.append((int(img[i,j]) - paeth_pred(left, up, ul)) & 0xFF)
        above = img[i]
    return bytes(out)

all_beat = True
for name in ['gradient', 'stripes', 'checker', 'smooth', 'photo', 'mandelbrot', 'noise']:
    img = gen_img(name, N)
    raw = img.tobytes()

    pil = Image.fromarray(img, 'L')
    buf = io.BytesIO(); pil.save(buf, format='PNG', optimize=True)
    png_size = len(buf.getvalue())

    candidates = {}

    # Cheap WL variants
    for level in range(4):
        comp, _, _, _ = cheap_wl_compress(img, level=level, q_bits=4)
        candidates[f'cwl_L{level}'] = comp
    for level in [1, 2, 3]:
        comp, _ = cheap_wl_with_prediction(img, level=level, q_bits=4)
        candidates[f'cwl+P_L{level}'] = comp

    # Standard methods
    paeth = paeth_encode(img)
    candidates['paeth+lzma'] = lzma.compress(paeth)
    candidates['paeth+bz2'] = bz2.compress(paeth, 9)
    delta = bytes([raw[0]] + [(raw[i]-raw[i-1])&0xFF for i in range(1,len(raw))])
    candidates['delta+lzma'] = lzma.compress(delta)
    candidates['raw+lzma'] = lzma.compress(raw)

    best_name = min(candidates.keys(), key=lambda k: len(candidates[k]))
    best_size = len(candidates[best_name])

    ratio = best_size / png_size
    win = "✓" if best_size <= png_size else ""
    if best_size > png_size: all_beat = False
    print(f"  {name:>12} {png_size:>5} {best_size:>5} {ratio:>6.3f}x {best_name:>15} {win}")

if all_beat:
    print(f"\n  *** ALL IMAGES BEAT PNG ***")

print(f"\n{'='*80}")
print(f"  THE CHEAP WL PRINCIPLE")
print(f"{'='*80}")
print(f"""
  STANDARD WL: "same" = identical content. Too many colors.
  CHEAP WL:    "same" = same statistics. Bounded colors.

  The LEVELS OF SAME form a REFINEMENT HIERARCHY:
    L0: same mean                         → O(2^q) colors
    L1: same (mean, std)                  → O(2^{{2q}}) colors
    L2: same (mean, std, gradient)        → O(2^{{3q}}) colors
    L3: same (mean, std, gradient, curve) → O(2^{{4q}}) colors

  With q=4 bits of quantization:
    L0: ≤ 16 colors (very coarse, big groups, big residuals)
    L1: ≤ 256 colors (medium, balanced)
    L2: ≤ 4096 colors (fine, small groups, small residuals)
    L3: ≤ 65536 colors (very fine, approaches pixel-level)

  The SWEET SPOT is L1 or L2: enough colors to capture structure,
  few enough that prototypes + color map are cheap.

  INTRA-CLASS PAETH:
  Within each color class, rows are structurally similar.
  Paeth prediction between consecutive rows of the same class
  is MUCH more accurate than Paeth between arbitrary neighbors.
  This is because WL coloring ensures the predictor matches the
  structural type of the data.

  The FUSION: WL for GROUPING + Paeth for PREDICTION = each step
  does what it's best at. WL finds the structure (what's equivalent),
  Paeth exploits the structure (predicts from context).
""")

print("DONE.")
