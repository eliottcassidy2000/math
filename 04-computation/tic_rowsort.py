#!/usr/bin/env python3
"""
ROWSORT CODEC: Sort rows by similarity, then delta-compress.

THIS IS GENUINELY NEW: No standard image codec sorts rows.
It's like BWT (Burrows-Wheeler) applied to image rows.

Why it works: sorting groups similar rows adjacent.
Delta-coding adjacent similar rows gives near-zero residuals.
For images with rotational/concentric symmetry, this is devastating.

Extensions to try:
1. Sort by content (full row) vs sort by feature (sum, variance, etc.)
2. Sort rows then MED on sorted image
3. Sort rows AND columns (2D sort)
4. Cluster rows (k-means) then encode per-cluster + residuals
5. Sort then apply Paeth/MED on sorted image

kind-pasteur-2026-03-25-S8
"""
import sys, io, struct, zlib, math, time
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

try: import brotli; HAS_BROTLI = True
except: HAS_BROTLI = False

def bc(data):
    r = [zlib.compress(data, 9)]
    if HAS_BROTLI:
        try: r.append(brotli.compress(data, quality=11))
        except: pass
    return min(r, key=len)

def png_size(img):
    pil = Image.fromarray(img, 'L' if img.ndim == 2 else 'RGB')
    buf = io.BytesIO(); pil.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()

def med(a, b, c):
    mn, mx = min(a, b), max(a, b)
    if c >= mx: return mn
    if c <= mn: return mx
    return a + b - c

# ============================================================
# ROWSORT VARIANTS
# ============================================================

def rowsort_delta(plane):
    """V1: Sort rows by content, delta-compress."""
    h, w = plane.shape
    indices = sorted(range(h), key=lambda i: plane[i].tobytes())
    sorted_plane = np.stack([plane[i] for i in indices])
    delta = np.zeros_like(sorted_plane)
    delta[0] = sorted_plane[0]
    delta[1:] = ((sorted_plane[1:].astype(int) - sorted_plane[:-1].astype(int)) & 0xFF).astype(np.uint8)
    perm = np.array(indices, dtype=np.uint8 if h <= 256 else np.uint16)
    return len(bc(perm.tobytes())) + len(bc(delta.tobytes()))

def rowsort_med(plane):
    """V2: Sort rows by content, then MED on sorted image."""
    h, w = plane.shape
    indices = sorted(range(h), key=lambda i: plane[i].tobytes())
    sp = np.stack([plane[i] for i in indices])
    res = bytearray()
    for r in range(h):
        for c in range(w):
            a = int(sp[r,c-1]) if c>0 else 0
            b = int(sp[r-1,c]) if r>0 else 0
            c2 = int(sp[r-1,c-1]) if r>0 and c>0 else 0
            res.append((int(sp[r,c]) - med(a,b,c2)) & 0xFF)
    perm = np.array(indices, dtype=np.uint8 if h <= 256 else np.uint16)
    return len(bc(perm.tobytes())) + len(bc(bytes(res)))

def rowsort_bysum(plane):
    """V3: Sort rows by SUM (brightness). Cheapest sort key."""
    h, w = plane.shape
    sums = [int(plane[i].sum()) for i in range(h)]
    indices = sorted(range(h), key=lambda i: sums[i])
    sp = np.stack([plane[i] for i in indices])
    delta = np.zeros_like(sp)
    delta[0] = sp[0]
    delta[1:] = ((sp[1:].astype(int) - sp[:-1].astype(int)) & 0xFF).astype(np.uint8)
    perm = np.array(indices, dtype=np.uint8 if h <= 256 else np.uint16)
    return len(bc(perm.tobytes())) + len(bc(delta.tobytes()))

def colsort_delta(plane):
    """V4: Sort COLUMNS by content, then delta-compress."""
    h, w = plane.shape
    indices = sorted(range(w), key=lambda c: plane[:,c].tobytes())
    sp = np.stack([plane[:,c] for c in indices], axis=1)
    delta = np.zeros_like(sp)
    delta[:,0] = sp[:,0]
    delta[:,1:] = ((sp[:,1:].astype(int) - sp[:,:-1].astype(int)) & 0xFF).astype(np.uint8)
    perm = np.array(indices, dtype=np.uint8 if w <= 256 else np.uint16)
    return len(bc(perm.tobytes())) + len(bc(delta.tobytes()))

def rowsort_then_colsort(plane):
    """V5: 2D SORT — sort rows first, then sort columns."""
    h, w = plane.shape
    # Sort rows
    ridx = sorted(range(h), key=lambda i: plane[i].tobytes())
    rp = np.stack([plane[i] for i in ridx])
    # Sort columns of sorted image
    cidx = sorted(range(w), key=lambda c: rp[:,c].tobytes())
    cp = np.stack([rp[:,c] for c in cidx], axis=1)
    # Delta in both directions
    delta = np.zeros_like(cp)
    delta[0,:] = cp[0,:]
    delta[1:,:] = ((cp[1:,:].astype(int) - cp[:-1,:].astype(int)) & 0xFF).astype(np.uint8)
    rperm = np.array(ridx, dtype=np.uint8 if h <= 256 else np.uint16)
    cperm = np.array(cidx, dtype=np.uint8 if w <= 256 else np.uint16)
    return len(bc(rperm.tobytes())) + len(bc(cperm.tobytes())) + len(bc(delta.tobytes()))

def rowsort_hilbert_key(plane):
    """V6: Sort rows by a Hilbert-curve-inspired key that captures
    both brightness AND spatial frequency."""
    h, w = plane.shape
    keys = []
    for i in range(h):
        row = plane[i]
        mean_val = int(row.mean())
        variance = int(np.var(row))
        # Hilbert-like interleaving of mean and variance bits
        key = (mean_val << 8) | min(255, variance)
        keys.append((key, i))
    keys.sort()
    indices = [k[1] for k in keys]
    sp = np.stack([plane[i] for i in indices])
    delta = np.zeros_like(sp)
    delta[0] = sp[0]
    delta[1:] = ((sp[1:].astype(int) - sp[:-1].astype(int)) & 0xFF).astype(np.uint8)
    perm = np.array(indices, dtype=np.uint8 if h <= 256 else np.uint16)
    return len(bc(perm.tobytes())) + len(bc(delta.tobytes()))

def rowsort_nearest_neighbor(plane):
    """V7: Greedy nearest-neighbor sort — each row followed by its most
    similar unused row. Like TSP on rows."""
    h, w = plane.shape
    used = [False] * h
    order = [0]; used[0] = True

    for _ in range(h - 1):
        cur = order[-1]
        best_dist = float('inf')
        best_next = -1
        for j in range(h):
            if used[j]: continue
            dist = int(np.sum(np.abs(plane[cur].astype(int) - plane[j].astype(int))))
            if dist < best_dist:
                best_dist = dist
                best_next = j
        order.append(best_next)
        used[best_next] = True

    sp = np.stack([plane[i] for i in order])
    delta = np.zeros_like(sp)
    delta[0] = sp[0]
    delta[1:] = ((sp[1:].astype(int) - sp[:-1].astype(int)) & 0xFF).astype(np.uint8)
    perm = np.array(order, dtype=np.uint8 if h <= 256 else np.uint16)
    return len(bc(perm.tobytes())) + len(bc(delta.tobytes()))

def encode_med_standard(plane):
    h, w = plane.shape; res = bytearray()
    for r in range(h):
        for c in range(w):
            a = int(plane[r,c-1]) if c>0 else 0
            b = int(plane[r-1,c]) if r>0 else 0
            c2 = int(plane[r-1,c-1]) if r>0 and c>0 else 0
            res.append((int(plane[r,c]) - med(a,b,c2)) & 0xFF)
    return len(bc(bytes(res)))

# ============================================================
# TEST
# ============================================================

SZ = 128

def make_tests():
    T = {}; np.random.seed(42)
    x, y = np.meshgrid(np.arange(SZ, dtype=float), np.arange(SZ, dtype=float))
    r = np.sqrt((x-SZ/2)**2 + (y-SZ/2)**2)

    T["solid"] = np.full((SZ,SZ), 128, dtype=np.uint8)
    T["grad_h"] = np.tile(np.linspace(0,255,SZ,dtype=np.uint8),(SZ,1))
    T["grad_v"] = np.tile(np.linspace(0,255,SZ,dtype=np.uint8).reshape(-1,1),(1,SZ))
    T["grad_d"] = ((x+y)*255//(2*SZ-2)).astype(np.uint8)
    T["circles"] = (np.sin(r/5)*127+128).astype(np.uint8)
    T["blob"] = (np.exp(-r**2/(2*(SZ/4)**2))*255).astype(np.uint8)
    T["radial"] = np.clip(r*255/(SZ/2),0,255).astype(np.uint8)
    T["random"] = np.random.randint(0,256,(SZ,SZ),dtype=np.uint8)
    sm = np.random.randint(0,256,(SZ//16,SZ//16),dtype=np.uint8)
    T["natural"] = np.clip(np.array(Image.fromarray(sm).resize((SZ,SZ),Image.BILINEAR)).astype(float)
                           +np.random.normal(0,10,(SZ,SZ)),0,255).astype(np.uint8)
    T["stripes0"] = (y%8<4).astype(np.uint8)*255
    T["stripes45"] = (((x+y)/8).astype(int)%2*255).astype(np.uint8)
    T["stripes90"] = (x%8<4).astype(np.uint8)*255
    return T

print("=" * 100)
print("  ROWSORT CODEC: 7 variants of row-sorting for compression")
print("  kind-pasteur-2026-03-25-S8")
print("=" * 100)

tests = make_tests()

METHODS = [
    ("MED std",     encode_med_standard),
    ("RS delta",    rowsort_delta),
    ("RS MED",      rowsort_med),
    ("RS bysum",    rowsort_bysum),
    ("CS delta",    colsort_delta),
    ("RS+CS 2D",    rowsort_then_colsort),
    ("RS Hilbert",  rowsort_hilbert_key),
    ("RS NN-TSP",   rowsort_nearest_neighbor),
]

print(f"\n  {'Image':<12} {'PNG':>6}", end="")
for mname, _ in METHODS:
    print(f" {mname:>9}", end="")
print(f"  {'BEST':>10}")
print("  " + "-" * (20 + 10*len(METHODS) + 15))

for name, img in sorted(tests.items()):
    ps = png_size(img)
    sizes = {}
    for mname, mfunc in METHODS:
        sizes[mname] = mfunc(img) + 10

    best = min(sizes, key=sizes.get)
    best_vs = f"{ps/sizes[best]:.2f}x" if sizes[best] > 0 else "inf"

    print(f"  {name:<12} {ps:>6}", end="")
    for mname, _ in METHODS:
        marker = "*" if mname == best else " "
        print(f" {sizes[mname]:>8}{marker}", end="")
    print(f"  {best:>10} {best_vs}")

# Show what each sort looks like for the blob
print(f"\n  === ROW SIMILARITY AFTER SORTING (blob image) ===")
img = tests["blob"]
indices = sorted(range(SZ), key=lambda i: img[i].tobytes())
sp = np.stack([img[i] for i in indices])
# Row-to-row difference
diffs = np.abs(sp[1:].astype(int) - sp[:-1].astype(int))
mean_diffs = diffs.mean(axis=1)
print(f"  Original image: mean row-to-row diff = {np.abs(img[1:].astype(int) - img[:-1].astype(int)).mean():.2f}")
print(f"  Sorted image: mean row-to-row diff = {mean_diffs.mean():.2f}")
print(f"  Reduction: {np.abs(img[1:].astype(int) - img[:-1].astype(int)).mean() / mean_diffs.mean():.2f}x")

print(f"\n  === ROW SIMILARITY AFTER NN-TSP SORT (blob) ===")
h = SZ
used = [False] * h
order = [0]; used[0] = True
for _ in range(h - 1):
    cur = order[-1]; bd = float('inf'); bn = -1
    for j in range(h):
        if used[j]: continue
        d = int(np.sum(np.abs(img[cur].astype(int) - img[j].astype(int))))
        if d < bd: bd = d; bn = j
    order.append(bn); used[bn] = True
sp_nn = np.stack([img[i] for i in order])
d_nn = np.abs(sp_nn[1:].astype(int) - sp_nn[:-1].astype(int)).mean(axis=1)
print(f"  NN-TSP sorted: mean row-to-row diff = {d_nn.mean():.2f}")
print(f"  Reduction vs original: {np.abs(img[1:].astype(int) - img[:-1].astype(int)).mean() / d_nn.mean():.2f}x")
