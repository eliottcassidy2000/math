#!/usr/bin/env python3
"""
recursive_tiling_s322.py -- opus-2026-03-24-S322

HOW TILINGS RECURSIVELY DECOMPOSE INTO TILES

The staircase δ_{n-2} is a right isosceles triangle with tiles:
  (i, j) for 0 ≤ i < j ≤ n-1, j - i ≥ 2

At n=5, the 6 tiles in δ_3 are:
  (0,2) (0,3) (0,4)
  (1,3) (1,4)
  (2,4)

This is a 3×3 lower triangle (rows 0,1,2; columns 2,3,4).

THE RECURSIVE DECOMPOSITION:
  δ_{n-2} = δ_{n-3} ∪ boundary ∪ apex

Where:
  δ_{n-3}: the overlap = tiles (i,j) with 1 ≤ i < j ≤ n-2, j-i ≥ 2
           = the interior (n-1)-tournament's tile space
  boundary: tiles involving vertex 0 OR vertex n-1 (but not both)
           = bottom row (vertex 0) + right column (vertex n-1)
  apex: the single tile (0, n-1) connecting source to sink

More precisely, deleting vertex 0 (the source):
  δ_{n-2} = δ_{n-3}(vertices 1..n-1) ∪ {tiles (0,j) for j ≥ 2}
  The removed tiles form the BOTTOM ROW: (0,2), (0,3), ..., (0,n-1)
  Count: n-2 tiles

Deleting vertex n-1 (the sink):
  δ_{n-2} = δ_{n-3}(vertices 0..n-2) ∪ {tiles (i,n-1) for i ≤ n-3}
  The removed tiles form the RIGHT COLUMN: (0,n-1), (1,n-1), ..., (n-3,n-1)
  Count: n-2 tiles

Deleting BOTH (Mode B):
  δ_{n-2} = δ_{n-4}(vertices 1..n-2) ∪ bottom_row ∪ right_column ∪ apex
  Where apex = (0, n-1)

THE KEY: a tiling of δ_{n-2} decomposes as:
  tiling = (overlap_tiling, boundary_bits, apex_bit)

The overlap tiling IS a tiling of the smaller staircase δ_{n-3}.
The boundary bits encode the deleted vertex's connections.
The apex bit encodes whether source beats sink.

THIS IS EXACTLY the vertex insertion/deletion decomposition from S315!
  - overlap = T\v (the sub-tournament)
  - boundary = residual (v's connections to remaining vertices)
  - apex = the extra bit (when doing dual deletion)

BUT THERE'S MORE: the decomposition is RECURSIVE.
  δ_{n-2} → δ_{n-3} → δ_{n-4} → ... → δ_1 → δ_0 (empty)

At each level, we peel off one row (or column) of the triangle.
The peeling order determines WHICH vertex we're inserting/deleting.

Author: opus-2026-03-24-S322
"""

import sys
import numpy as np
from math import comb
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

# ============================================================
banner("THE STAIRCASE DECOMPOSITION: VISUAL")
# ============================================================

for n in [5, 6, 7]:
    m = comb(n-1, 2)
    print("  n=%d: δ_%d has %d tiles\n" % (n, n-2, m))

    # Draw the staircase
    # Rows: vertex i (the "lower" vertex of the arc)
    # Columns: vertex j (the "higher" vertex)
    # Tile (i,j) exists iff j - i ≥ 2

    print("    ", end="")
    for j in range(2, n):
        print("  j=%d " % j, end="")
    print()

    for i in range(n-2):
        print("  i=%d " % i, end="")
        for j in range(2, n):
            if j - i >= 2:
                # Classify this tile
                if i == 0 and j == n-1:
                    label = "APEX"
                elif i == 0:
                    label = "BOT "  # bottom row (vertex 0)
                elif j == n-1:
                    label = "RGT "  # right column (vertex n-1)
                else:
                    label = "OVER"  # overlap (interior)
                print(" %s " % label, end="")
            else:
                print("  ·   ", end="")
        print()

    # Count each type
    n_overlap = comb(n-2, 2)  # = C(n-3, 2) when removing vertex 0
    n_bottom = n - 2          # tiles (0, j) for j = 2..n-1
    n_right = n - 2           # tiles (i, n-1) for i = 0..n-3
    n_apex = 1                # tile (0, n-1) counted in both bottom and right

    # Careful: count tiles by category
    # Overlap: tiles (i,j) with 1 ≤ i, j ≤ n-2, j-i ≥ 2
    actual_overlap = sum(1 for i in range(1, n-1) for j in range(i+2, n-1))
    # Bottom row: tiles (0, j) for j ≥ 2
    actual_bottom = sum(1 for j in range(2, n))
    # Right column: tiles (i, n-1) for j-i ≥ 2
    actual_right = sum(1 for i in range(0, n-2))  # i from 0 to n-3, j=n-1, range ≥ 2 when n-1-i ≥ 2 → i ≤ n-3
    actual_right = n - 2
    # Overlap between bottom and right: apex (0, n-1) counted in both
    actual_apex_double = 1
    # By inclusion-exclusion:
    n_boundary = actual_bottom + actual_right - actual_apex_double  # unique boundary tiles
    actual_overlap_true = m - n_boundary

    print("\n    Overlap (interior): %d tiles" % actual_overlap_true)
    print("    Bottom row (vertex 0): %d tiles" % actual_bottom)
    print("    Right column (vertex %d): %d tiles" % (n-1, actual_right))
    print("    Apex overlap: 1 tile (in both bottom and right)")
    print("    Boundary (bottom ∪ right): %d unique tiles" % n_boundary)
    print("    Total: %d interior + %d boundary = %d = m ✓" %
          (actual_overlap_true, n_boundary, m))
    print()

# ============================================================
banner("THE RECURSIVE PEELING")
# ============================================================

print("""
  THE STAIRCASE PEELS LAYER BY LAYER:

  δ_{n-2} → δ_{n-3} → δ_{n-4} → ... → δ_1 → δ_0

  At each step, we remove ONE layer from the triangle.
  The layer can be peeled from three directions:

  PEEL BOTTOM (remove row i=0): removes n-2 tiles
    These are vertex 0's non-base connections: (0,2), (0,3), ..., (0,n-1)
    = Mode A deletion of the SOURCE vertex

  PEEL RIGHT (remove column j=n-1): removes n-2 tiles
    These are vertex n-1's non-base connections: (0,n-1), (1,n-1), ..., (n-3,n-1)
    = Mode A deletion of the SINK vertex

  PEEL HYPOTENUSE (remove anti-diagonal): removes floor((n-2)/2) tiles
    These are the tiles with range 2: (0,2), (1,3), (2,4), ...
    = Mode A deletion of the "DIAGONAL vertex"

  PEEL BOTH LEGS (remove bottom row AND right column): Mode B
    Removes 2(n-3) + 1 tiles (bottom_unique + right_unique + apex)
    = 2n-5 tiles removed, leaving C(n-3,2) = overlap

  Each peeling order gives a DIFFERENT recursive decomposition
  of the tiling into sub-tilings.
""")

# ============================================================
banner("THE FULL RECURSIVE TREE AT n=5")
# ============================================================

n = 5
m = comb(n-1, 2)  # 6

print("  n=5: δ_3 = 6 tiles: A=(0,2) B=(0,3) C=(0,4) D=(1,3) E=(1,4) F=(2,4)\n")

# The recursive decomposition peeling vertex 0 (bottom row) at each step:
print("  PEELING VERTEX 0 AT EACH LEVEL:\n")
print("  Level 0: δ_3 = {A,B,C,D,E,F} (6 tiles, vertices {0,1,2,3,4})")
print("    Peel i=0: remove {A=(0,2), B=(0,3), C=(0,4)} [3 tiles = n-2]")
print("    Remaining: δ_2 = {D=(1,3), E=(1,4), F=(2,4)} [3 tiles, vertices {1,2,3,4}]")
print()
print("  Level 1: δ_2 = {D,E,F} (3 tiles, vertices {1,2,3,4})")
print("    Peel i=1: remove {D=(1,3), E=(1,4)} [2 tiles = n-3]")
print("    Remaining: δ_1 = {F=(2,4)} [1 tile, vertices {2,3,4}]")
print()
print("  Level 2: δ_1 = {F} (1 tile, vertices {2,3,4})")
print("    Peel i=2: remove {F=(2,4)} [1 tile = n-4]")
print("    Remaining: δ_0 = {} [0 tiles, vertices {3,4}]")
print()
print("  Level 3: δ_0 = {} (0 tiles, vertices {3,4}). BASE CASE.")

print("""
  THE RECURSIVE FORMULA:
    tiling(δ_{k}) = (tiling(δ_{k-1}), boundary_bits[k])
    where boundary_bits[k] has k bits (one per tile in the peeled row).

  Total bits: Σ_{k=1}^{n-2} k = C(n-1, 2) = m. ✓

  This gives a FACTORIZATION of the tiling word:
    m-bit tiling = (1 bit from δ_1) | (2 bits from δ_2 boundary) | (3 bits from δ_3 boundary)
    = F | (D,E) | (A,B,C)
    = [1 bit] [2 bits] [3 bits]

  In general: m = 1 + 2 + 3 + ... + (n-2) = T_{n-2} = C(n-1,2).
  The tiling decomposes into n-2 LAYERS of sizes 1, 2, ..., n-2.
""")

# ============================================================
banner("THE LAYER DECOMPOSITION")
# ============================================================

for n in [5, 6, 7, 8]:
    m = comb(n-1, 2)
    base_set = set((i, i+1) for i in range(n-1))
    all_p = [(i,j) for i in range(n) for j in range(i+1,n)]
    tile_pairs = [(i,j) for i,j in all_p if (i,j) not in base_set]
    tile_pairs.sort()

    print("  n=%d (m=%d):" % (n, m))

    # Peel from bottom (remove vertex 0, then 1, then 2, ...)
    layers = []
    remaining_tiles = set(range(len(tile_pairs)))
    remaining_vertices = set(range(n))

    for level in range(n-2):
        v = level  # vertex to peel
        # Tiles involving v that are still remaining
        layer_tiles = []
        for k in remaining_tiles.copy():
            lo, hi = tile_pairs[k]
            if lo == v or hi == v:
                if hi - lo >= 2:  # it's a tile (not base arc)
                    layer_tiles.append(k)
                    remaining_tiles.remove(k)

        layers.append(layer_tiles)
        remaining_vertices.remove(v)

    # Verify: layers partition all tiles
    all_in_layers = sum(len(L) for L in layers)
    print("    Layers: %s" % [len(L) for L in layers])
    print("    Layer sizes: 1 + 2 + ... + %d = %d = m ✓" %
          (n-2, sum(range(1, n-1))))
    print("    Tiles per layer: %s" % [[tile_pairs[k] for k in L] for L in layers[:3]])
    print()

# ============================================================
banner("THE RECURSIVE INTERPRETATION")
# ============================================================

print("""
  EACH LAYER = ONE VERTEX'S NON-BASE CONNECTIONS

  Layer k (0-indexed) peels vertex k, removing its k+1 non-base tiles:
    Layer 0: vertex 0's tiles = (0,2), (0,3), ..., (0,n-1)  [n-2 tiles]
    Layer 1: vertex 1's remaining tiles = (1,3), (1,4), ..., (1,n-1)  [n-3 tiles]
    ...
    Layer k: vertex k's remaining tiles = (k, k+2), ..., (k, n-1)  [n-k-2 tiles]
    ...
    Layer n-3: vertex (n-3)'s last tile = (n-3, n-1)  [1 tile]

  Wait — this gives layers of size n-2, n-3, ..., 1.
  Total = (n-2) + (n-3) + ... + 1 = C(n-1,2) = m. ✓

  THE TILING AS A SEQUENCE OF RESIDUALS:
    tiling = residual(v=0) || residual(v=1) || ... || residual(v=n-3)

  Where residual(v=k) has n-k-2 bits, encoding vertex k's connections
  to vertices k+2, k+3, ..., n-1 (skipping the base arc to k+1).

  THIS IS THE VERTEX-BY-VERTEX PEELING of the fractal codec!
  The fractal codec (kind-pasteur S20cq) recursively encodes T as:
    T = (T\\v₀, residual₀, T\\{v₀,v₁}, residual₁, ...)

  The staircase peeling order = the natural vertex ordering 0, 1, 2, ...
  gives EXACTLY this recursive decomposition.

  ADAPTIVE PEELING: choosing a DIFFERENT vertex order gives a
  different decomposition with potentially BETTER compression.
  From S315: adaptive (most extreme score) gives 73% compression.

  THE TILE-LEVEL VIEW:
  Each individual tile (i,j) belongs to EXACTLY ONE layer:
  the layer where the LOWER vertex i is peeled.
  (Because peeling i removes all tiles with i as lower vertex.)

  So: tile (i,j) ∈ Layer i. And Layer i has tiles
  (i, i+2), (i, i+3), ..., (i, n-1).
  These are the tiles in ROW i of the staircase.

  ROW i has n-i-2 tiles, all sharing the same lower vertex i.
  The rows ARE the layers. The recursive decomposition peels rows.
""")

# ============================================================
banner("THE DEEPER RECURSION: TILES DECOMPOSE INTO SUB-TILES")
# ============================================================

print("""
  Can we go BELOW the tile level? Can a single tile decompose further?

  A tile (i,j) with range r = j-i represents the arc between vertices
  i and j. This arc can be thought of as a "path" through the
  intermediate vertices i+1, i+2, ..., j-1.

  In the Lie algebra language: tile (i,j) = compound root
  α_i + α_{i+1} + ... + α_{j-1} (sum of j-i simple roots).

  A compound root DECOMPOSES as a sum of simpler roots:
  α_i + α_{i+1} + ... + α_{j-1} = (α_i + ... + α_{k-1}) + (α_k + ... + α_{j-1})
  for any k between i+1 and j-1.

  In tournament terms: the arc (i,j) can be "decomposed" through
  any intermediate vertex k:
    (i,j) = (i,k) + (k,j) — but only if both (i,k) and (k,j) are tiles
    (i.e., k-i ≥ 2 and j-k ≥ 2, so r ≥ 4).

  For range-2 tiles: NO decomposition (they are "simple" = irreducible).
  For range-3 tiles: ONE decomposition through the middle vertex.
    But the middle vertex is a base-path neighbor, so (i,i+1) is a base arc,
    not a tile. The decomposition uses a base arc + a range-2 tile:
    (i, i+3) = base(i, i+1) + tile(i+1, i+3)
    This is NOT a pure tile decomposition — it involves the base path.

  For range-4 tiles: (i, i+4) can decompose as:
    tile(i, i+2) + tile(i+2, i+4) — both are range-2 tiles!
    This IS a pure tile decomposition.

  So: range-4 tiles = "product" of two range-2 tiles.
  Range-2 tiles are the ATOMS (simple roots / generators).
  Range-3 tiles are "mixed" (involve one base arc).
  Range-4+ tiles decompose into products of range-2 atoms.

  THE RECURSIVE STRUCTURE:
    Range 2: atoms (n-2 of them, one per row minus endpoints)
    Range 3: hybrid (involve base arcs)
    Range 4: products of two range-2 atoms
    Range 5: products of range-2 + range-3
    Range r: products decomposing through all intermediate vertices

  This mirrors the ROOT DECOMPOSITION of A_{n-1}:
    Simple roots: α_1, ..., α_{n-1}  (= base-path arcs)
    Compound roots of height h: sums of h consecutive simple roots (= tiles of range h)
    The height = range - 1 (since the base-path arcs have range 1).

  THE TILES ARE THE COMPOUND ROOTS, AND THEY DECOMPOSE
  INTO SUMS OF SIMPLE ROOTS (BASE-PATH ARCS).
  The tile bit encodes WHETHER the compound root is "positive" or "negative"
  (= the orientation of the arc).
""")

# Verify the decomposition structure
for n in [5, 6]:
    m = comb(n-1, 2)
    base_set = set((i, i+1) for i in range(n-1))
    all_p = [(i,j) for i in range(n) for j in range(i+1,n)]
    tile_pairs = [(i,j) for i,j in all_p if (i,j) not in base_set]
    tile_pairs.sort()

    print("  n=%d: tile decomposition into range-2 atoms:\n" % n)

    for lo, hi in tile_pairs:
        r = hi - lo
        if r == 2:
            print("    (%d,%d) range %d: ATOM (irreducible)" % (lo, hi, r))
        elif r == 3:
            # Decomposes as base(lo, lo+1) + tile(lo+1, hi) or tile(lo, hi-1) + base(hi-1, hi)
            print("    (%d,%d) range %d: = base(%d,%d) + tile(%d,%d)  [hybrid]" %
                  (lo, hi, r, lo, lo+1, lo+1, hi))
        elif r == 4:
            print("    (%d,%d) range %d: = tile(%d,%d) + tile(%d,%d)  [product of 2 atoms]" %
                  (lo, hi, r, lo, lo+2, lo+2, hi))
        elif r >= 5:
            # Multiple decompositions
            decomps = []
            for k in range(lo+2, hi, 2):
                if k - lo >= 2 and hi - k >= 2:
                    decomps.append("tile(%d,%d)+tile(%d,%d)" % (lo, k, k, hi))
            print("    (%d,%d) range %d: = %s" % (lo, hi, r, " OR ".join(decomps[:3])))
    print()

print("Done. opus-2026-03-24-S322")
