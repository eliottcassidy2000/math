#!/usr/bin/env python3
"""
MST-WALK CODEC: The most radical idea yet.

Break the assumption: "visit pixels in a predetermined order."
Instead: build the MINIMUM SPANNING TREE of the pixel adjacency graph,
then walk the tree. Each pixel is predicted from its parent.

WHY THIS IS PROFOUND:
  - The MST naturally follows smooth gradients (minimum-weight path)
  - It avoids edges (high-weight edges are excluded)
  - For radial patterns, the MST spirals outward automatically
  - For gradients, the MST follows the gradient direction
  - The tree adapts to the CONTENT, not a fixed scan pattern

THE STAIRCASE PARADOX IS ELIMINATED: the MST doesn't scan in any
fixed direction — it flows along the image's natural structure.

COST: need to encode the tree structure (parent pointers).
GAIN: prediction residuals are minimized (each pixel predicted from
its closest similar neighbor).

For a 128×128 image: 16384 pixels, each needs a parent pointer
(2 bits for direction if only 4-connected: N/S/E/W).
Tree structure cost: 16384 * 2 / 8 = 4096 bytes.
But residuals should be MUCH smaller since every prediction follows
the minimum-difference path.

kind-pasteur-2026-03-25-S8
"""
import sys, io, struct, zlib, math, time
import numpy as np
from PIL import Image
import heapq

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

def encode_med(plane):
    h, w = plane.shape; res = bytearray()
    for r in range(h):
        for c in range(w):
            a = int(plane[r,c-1]) if c>0 else 0
            b = int(plane[r-1,c]) if r>0 else 0
            c2 = int(plane[r-1,c-1]) if r>0 and c>0 else 0
            res.append((int(plane[r,c]) - med(a,b,c2)) & 0xFF)
    return len(bc(bytes(res)))

# ============================================================
# MST-WALK CODEC
# ============================================================

def build_mst_4connected(plane):
    """Build MST of 4-connected pixel graph using Prim's algorithm.
    Edge weight = absolute difference between adjacent pixels.
    Returns: parent array where parent[r,c] = (pr, pc) or (-1,-1) for root."""
    h, w = plane.shape
    parent = np.full((h, w, 2), -1, dtype=np.int16)
    visited = np.zeros((h, w), dtype=bool)

    # Start from center (or top-left)
    sr, sc = h // 2, w // 2
    heap = [(0, sr, sc, -1, -1)]  # (weight, r, c, parent_r, parent_c)

    while heap:
        weight, r, c, pr, pc = heapq.heappop(heap)
        if visited[r, c]:
            continue
        visited[r, c] = True
        parent[r, c] = [pr, pc]

        # Add neighbors
        for dr, dc in [(-1,0),(1,0),(0,-1),(0,1)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < h and 0 <= nc < w and not visited[nr, nc]:
                w2 = abs(int(plane[nr, nc]) - int(plane[r, c]))
                heapq.heappush(heap, (w2, nr, nc, r, c))

    return parent

def mst_walk_encode(plane):
    """Encode image by walking the MST. Each pixel predicted from its parent."""
    h, w = plane.shape
    parent = build_mst_4connected(plane)

    # BFS walk from root
    sr, sc = h // 2, w // 2
    visited = np.zeros((h, w), dtype=bool)
    queue = [(sr, sc)]
    visited[sr, sc] = True

    # Walk order
    walk_order = [(sr, sc)]
    while queue:
        r, c = queue.pop(0)
        for dr, dc in [(-1,0),(1,0),(0,-1),(0,1)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < h and 0 <= nc < w and not visited[nr, nc]:
                # Check if (nr, nc) is a child of (r, c) in the MST
                if parent[nr, nc, 0] == r and parent[nr, nc, 1] == c:
                    visited[nr, nc] = True
                    walk_order.append((nr, nc))
                    queue.append((nr, nc))

    # Some pixels might not be reachable through BFS children only
    # Add remaining in raster order
    for r in range(h):
        for c in range(w):
            if not visited[r, c]:
                walk_order.append((r, c))
                visited[r, c] = True

    # Encode: parent directions (2 bits each) + residuals
    # Direction encoding: 0=N(-1,0), 1=S(+1,0), 2=W(0,-1), 3=E(0,+1), 4=root
    dir_map = {(-1,0): 0, (1,0): 1, (0,-1): 2, (0,1): 3}
    directions = bytearray()
    residuals = bytearray()

    for r, c in walk_order:
        pr, pc = parent[r, c]
        if pr == -1 and pc == -1:
            # Root
            directions.append(4)
            residuals.append(plane[r, c])
        else:
            dr, dc = pr - r, pc - c
            d = dir_map.get((dr, dc), 4)
            directions.append(d)
            pred = int(plane[pr, pc])
            residuals.append((int(plane[r, c]) - pred) & 0xFF)

    dir_compressed = bc(bytes(directions))
    res_compressed = bc(bytes(residuals))

    return len(dir_compressed) + len(res_compressed), directions, residuals

def mst_walk_stats(plane):
    """Compute MST edge weight statistics."""
    h, w = plane.shape
    parent = build_mst_4connected(plane)

    total_weight = 0
    n_edges = 0
    for r in range(h):
        for c in range(w):
            pr, pc = parent[r, c]
            if pr >= 0:
                total_weight += abs(int(plane[r,c]) - int(plane[pr,pc]))
                n_edges += 1

    mean_weight = total_weight / n_edges if n_edges > 0 else 0
    return mean_weight

# ============================================================
# BENCHMARK
# ============================================================

SZ = 64  # Start small (MST is O(n log n))

def make_tests():
    T = {}; np.random.seed(42)
    x, y = np.meshgrid(np.arange(SZ, dtype=float), np.arange(SZ, dtype=float))
    r = np.sqrt((x-SZ/2)**2 + (y-SZ/2)**2)

    T["solid"] = np.full((SZ,SZ), 128, dtype=np.uint8)
    T["grad_h"] = np.tile(np.linspace(0,255,SZ,dtype=np.uint8),(SZ,1))
    T["grad_d"] = ((x+y)*255//(2*SZ-2)).astype(np.uint8)
    T["circles"] = (np.sin(r/5)*127+128).astype(np.uint8)
    T["blob"] = (np.exp(-r**2/(2*(SZ/4)**2))*255).astype(np.uint8)
    T["radial"] = np.clip(r*255/(SZ/2),0,255).astype(np.uint8)
    T["random"] = np.random.randint(0,256,(SZ,SZ),dtype=np.uint8)
    sm = np.random.randint(0,256,(SZ//16,SZ//16),dtype=np.uint8)
    T["natural"] = np.clip(np.array(Image.fromarray(sm).resize((SZ,SZ),Image.BILINEAR)).astype(float)
                           +np.random.normal(0,10,(SZ,SZ)),0,255).astype(np.uint8)
    return T

print("=" * 80)
print("  MST-WALK CODEC: Adaptive scan order from minimum spanning tree")
print("  kind-pasteur-2026-03-25-S8")
print("=" * 80)

tests = make_tests()

print(f"\n  {'Image':<12} {'PNG':>6} {'MED':>6} {'MST':>6} {'MST stats':>12} {'BEST':>10}")
print("  " + "-" * 60)

for name, img in sorted(tests.items()):
    ps = png_size(img)
    med_sz = encode_med(img) + 10
    mst_sz, dirs, resids = mst_walk_encode(img)
    mst_sz += 10

    mst_mean = mst_walk_stats(img)

    best_name = "MST" if mst_sz < med_sz else "MED"
    best_sz = min(mst_sz, med_sz)
    ratio = ps / best_sz if best_sz > 0 else 0

    # Analyze MST residual distribution
    r_signed = np.array(list(resids), dtype=np.int16)
    r_signed[r_signed > 128] -= 256
    mean_res = np.mean(np.abs(r_signed))

    print(f"  {name:<12} {ps:>6} {med_sz:>6} {mst_sz:>6} mw={mst_mean:>5.1f} mr={mean_res:>5.1f} {best_name:>6} {ratio:.2f}x")

print(f"""
  ANALYSIS:
  mw = mean MST edge weight (lower = better MST prediction)
  mr = mean absolute residual in MST walk

  MST ADVANTAGE:
  - The MST naturally finds the path of minimum prediction error
  - Every pixel is predicted from its closest (most similar) neighbor
  - The tree adapts to any image structure (radial, gradient, edges)

  MST DISADVANTAGE:
  - Tree structure costs ~n * log2(5) / 8 bytes (3 bits per pixel: 4 dirs + root)
  - For 64×64 = 4096 pixels: ~1536 bytes of tree structure
  - This overhead must be less than the residual savings

  The MST codec wins when:
  - Residual savings > tree overhead
  - Happens for images with smooth regions separated by sharp edges
  - The MST follows smooth paths, avoiding edges
""")

# Compare residual distributions
print("  RESIDUAL DISTRIBUTION (blob):")
img = tests["blob"]
# MED residuals
med_res = []
for r in range(SZ):
    for c in range(SZ):
        a = int(img[r,c-1]) if c>0 else 0
        b = int(img[r-1,c]) if r>0 else 0
        c2 = int(img[r-1,c-1]) if r>0 and c>0 else 0
        med_res.append((int(img[r,c]) - med(a,b,c2)) & 0xFF)
med_signed = np.array(med_res, dtype=np.int16); med_signed[med_signed>128] -= 256

_, _, mst_res = mst_walk_encode(img)
mst_s = np.array(list(mst_res), dtype=np.int16); mst_s[mst_s>128] -= 256

print(f"  MED: mean|res| = {np.mean(np.abs(med_signed)):.2f}, max|res| = {np.max(np.abs(med_signed))}")
print(f"  MST: mean|res| = {np.mean(np.abs(mst_s)):.2f}, max|res| = {np.max(np.abs(mst_s))}")
print(f"  MST residuals are {np.mean(np.abs(med_signed))/np.mean(np.abs(mst_s)):.2f}x smaller")
print(f"  But tree structure costs {len(bc(bytes(list(_[1] for _ in [(mst_walk_encode(img))])[0])))} bytes")
