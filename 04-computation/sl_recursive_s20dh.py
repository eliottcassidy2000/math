#!/usr/bin/env python3
"""
sl_recursive_s20dh.py — SL_orbits via recursive tiling decomposition
kind-pasteur-2026-03-23-S20dh

THE RECURSIVE STRUCTURE (Mode B, fully unwound):

  Level 0: The full staircase delta_{n-2}
    = overlap(n-2) + bottom_wiring(n-3 bits) + top_wiring(n-3 bits) + apex(1 bit)

  Level 1: The overlap delta_{n-4}
    = overlap(n-4) + bottom_wiring(n-5 bits) + top_wiring(n-5 bits) + apex(1 bit)

  Level k: overlap delta_{n-2k-2}
    = ... until n-2k <= 2 (base case: 0 or 1 tiles)

Each self-loop comes from flipping ONE tile. That tile lives at a specific
LEVEL and REGION of the recursion. By tracking which level and region
generates each self-loop, we can find the recurrence for SL_orbits.

We compute: for each (tiling, tile) self-loop pair, classify the tile as:
  - overlap_k (in the k-th level overlap, i.e., the n-2k sub-staircase)
  - bottom_k (in the k-th level bottom wiring)
  - top_k (in the k-th level top wiring)
  - apex_k (the k-th level apex)
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  SL_ORBITS RECURSIVE DECOMPOSITION")
print("  kind-pasteur-2026-03-23-S20dh")
print("=" * 80)

def analyze(n):
    t0 = time.time()

    # Explorer-order tiles
    VERTS = list(range(n, 0, -1))
    TILES = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            TILES.append((x, y))
    m = len(TILES)

    # Build the RECURSIVE decomposition of tiles into levels
    # Level 0: full staircase. Vertices = {1,...,n} (1-indexed in explorer).
    # Apex_0 = tile (n, 1) = arc between vertex n and vertex 1
    # Bottom_0 = tiles (x, 1) for x from n-1 down to 3 = arcs from vertex 1 to internals
    #            Actually: bottom = arcs involving vertex 1 (the "low" endpoint)
    #            In the staircase: tiles with y=1 and x >= 3 (non-apex, non-base)
    # Wait, I need to be more careful with the recursive structure.

    # The staircase for n vertices uses tiles (x,y) with y < x, x >= y+2.
    # The APEX at level 0 = tile (n, 1) = the arc between the two extreme vertices.
    # The OVERLAP at level 0 = tiles involving only internal vertices {2,...,n-1}
    #   = tiles (x,y) with 2 <= y < x <= n-1, x >= y+2
    #   = the staircase for n-2 vertices {2,...,n-1} (shifted)
    # Bottom_0 = tiles (x, 1) for x >= 3 and x <= n-1 = arcs from vertex 1 to internals
    #   but NOT the apex (n, 1). So: (3,1), (4,1), ..., (n-1, 1) = n-3 tiles.
    # Top_0 = tiles (n, y) for y >= 2 and y <= n-2 = arcs from vertex n to internals
    #   = (n, 2), (n, 3), ..., (n, n-2) = n-3 tiles.

    # Classify each tile into its recursive level and region
    tile_level_region = {}

    def classify_tiles(vertices, level):
        """Classify tiles among the given vertex set into overlap/bottom/top/apex at this level."""
        vmin = min(vertices)
        vmax = max(vertices)
        internals = sorted(set(vertices) - {vmin, vmax})

        for i, (x, y) in enumerate(TILES):
            if i in tile_level_region:
                continue
            if x not in vertices or y not in vertices:
                continue
            # This tile involves vertices in our set
            if x == vmax and y == vmin:
                tile_level_region[i] = (level, 'apex')
            elif y == vmin and x in internals:
                tile_level_region[i] = (level, 'bottom')
            elif x == vmax and y in internals:
                tile_level_region[i] = (level, 'top')
            # else: it's in the overlap (both endpoints are internal)
            # Will be classified at the next level

        # Recurse on internals if there are enough
        if len(internals) >= 3:
            classify_tiles(internals, level + 1)
        elif len(internals) >= 2:
            # Base case: 2 internal vertices, 1 tile between them (if distance >= 2)
            for i, (x, y) in enumerate(TILES):
                if i not in tile_level_region and x in internals and y in internals:
                    tile_level_region[i] = (level + 1, 'base')
        # Any remaining unclassified tiles at this level must be overlap
        # that couldn't recurse further
        for i, (x, y) in enumerate(TILES):
            if i not in tile_level_region and x in vertices and y in vertices:
                if x in internals and y in internals:
                    tile_level_region[i] = (level + 1, 'base')

    all_vertices = list(range(1, n+1))
    classify_tiles(all_vertices, 0)

    # Check all tiles classified
    unclassified = [i for i in range(m) if i not in tile_level_region]
    if unclassified:
        print(f"  WARNING: {len(unclassified)} unclassified tiles!")
        for i in unclassified:
            print(f"    tile {i} = {TILES[i]}")

    # Print the classification
    print(f"\n  n={n}: {m} tiles, recursive decomposition:")
    level_region_tiles = defaultdict(list)
    for i, (level, region) in sorted(tile_level_region.items()):
        level_region_tiles[(level, region)].append(TILES[i])

    for (level, region) in sorted(level_region_tiles.keys()):
        tiles_here = level_region_tiles[(level, region)]
        print(f"    Level {level} {region:>6}: {len(tiles_here)} tiles = {tiles_here}")

    # Now compute self-loops per tile, decomposed by level/region
    def bitsToAdj(bits):
        A = [[0]*n for _ in range(n)]
        for k in range(n-1): A[k][k+1] = 1
        for i in range(m):
            xL, yL = TILES[i]
            xi = VERTS.index(xL); yi = VERTS.index(yL)
            if bits[i] == 0: A[xi][yi] = 1
            else: A[yi][xi] = 1
        return A

    perms = list(permutations(range(n)))
    def canon(A):
        best = None
        for p in perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(n) for j in range(n))
            if best is None or s < best: best = s
        return best

    # Count self-loops per level/region
    sl_by_lr = Counter()  # (level, region) -> count of (tiling, tile) self-loop pairs
    total_sl = 0

    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = bitsToAdj(bits)
        cn = canon(A)

        for k in range(m):
            # Flip tile k
            bits2 = bits.copy()
            bits2[k] = 1 - bits2[k]
            A2 = bitsToAdj(bits2)
            cn2 = canon(A2)

            if cn == cn2:  # self-loop!
                lr = tile_level_region.get(k, (-1, '?'))
                sl_by_lr[lr] += 1
                total_sl += 1

        if (mask + 1) % 5000 == 0:
            print(f"    {mask+1}/{1<<m} ({time.time()-t0:.0f}s)")

    # SL_orbits = total_sl / n!
    sl_orbits = total_sl // factorial(n)

    print(f"\n  SELF-LOOP DECOMPOSITION (total SL_raw = {total_sl}, SL_orbits = {sl_orbits}):")
    print(f"  {'Level':>6} {'Region':>8} {'#tiles':>6} {'SL_raw':>8} {'SL/n!':>8} {'SL/tile':>8}")

    for (level, region) in sorted(sl_by_lr.keys()):
        count = sl_by_lr[(level, region)]
        ntiles = len(level_region_tiles.get((level, region), []))
        print(f"  {level:6d} {region:>8} {ntiles:6d} {count:8d} {count/factorial(n):8.1f} {count/ntiles if ntiles > 0 else 0:8.1f}")

    print(f"  {'':>6} {'TOTAL':>8} {m:6d} {total_sl:8d} {sl_orbits:8d}")

    return {
        'n': n, 'sl_orbits': sl_orbits, 'sl_by_lr': dict(sl_by_lr),
        'total_sl': total_sl
    }

# ============================================================================
# MAIN
# ============================================================================

results = {}
for n in [4, 5, 6, 7]:
    print(f"\n{'#'*60}")
    results[n] = analyze(n)

# ============================================================================
# CROSS-n RECURSIVE PATTERN
# ============================================================================

print(f"\n\n{'='*60}")
print(f"  CROSS-n: SL BY RECURSIVE LEVEL")
print(f"{'='*60}")

# For each level, track SL contribution across n
max_level = 3
for level in range(max_level + 1):
    for region in ['apex', 'bottom', 'top', 'base']:
        vals = []
        for n in [4, 5, 6, 7]:
            if n in results:
                v = results[n]['sl_by_lr'].get((level, region), 0)
                vals.append(v // factorial(n))
            else:
                vals.append('?')
        if any(v != 0 and v != '?' for v in vals):
            print(f"  Level {level} {region:>6}: SL_orbits = {vals}")

print(f"\n  Total SL_orbits: {[results[n]['sl_orbits'] for n in [4,5,6,7]]}")
print(f"  Known sequence: 5, 20, 86, 490")

print(f"\n  DONE.")
print("=" * 80)
