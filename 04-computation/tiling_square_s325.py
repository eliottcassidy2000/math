#!/usr/bin/env python3
"""
tiling_square_s325.py -- opus-2026-03-24-S325

THE TILING SQUARE: δ_{n-2} + δ_{n-3} = (n-2)² SQUARE

A tiling of size n (tournament on n vertices) fills the lower triangle δ_{n-2}
with C(n-1,2) = (n-1)(n-2)/2 bits.

A tiling of size n-1 (tournament on n-1 vertices) fills δ_{n-3}
with C(n-2,2) = (n-2)(n-3)/2 bits.

Together: C(n-1,2) + C(n-2,2) = (n-2)².

THE TWO TRIANGLES TILE A PERFECT SQUARE OF SIDE (n-2).

CONSTRUCTION:
  The (n-2) × (n-2) square has rows r = 0..n-3 and columns c = 0..n-3.
  Lower-left triangle: cells (r,c) with r > c. Count: C(n-2, 2).
  Upper-right triangle: cells (r,c) with r < c. Count: C(n-2, 2).
  Diagonal: cells (r,r). Count: n-2.

  Wait: C(n-2,2) + C(n-2,2) + (n-2) = (n-2)(n-3)/2 + (n-2)(n-3)/2 + (n-2)
  = (n-2)(n-3) + (n-2) = (n-2)(n-2) = (n-2)². ✓

  But the n-tiling has C(n-1,2) = (n-1)(n-2)/2 bits.
  And the (n-1)-tiling has C(n-2,2) = (n-2)(n-3)/2 bits.
  Sum: (n-1)(n-2)/2 + (n-2)(n-3)/2 = (n-2)(n-1+n-3)/2 = (n-2)(2n-4)/2 = (n-2)². ✓

  So: n-tiling = lower triangle + one extra row/column
  = C(n-2,2) cells in the lower triangle + (n-2) cells on the diagonal
  + the difference...

  Actually, let me think about this more carefully.

  The STAIRCASE δ_{n-2} has tiles (i,j) with 0 ≤ i < j ≤ n-1, j-i ≥ 2.
  There are C(n-1,2) = (n-1)(n-2)/2 such tiles.

  The staircase δ_{n-3} has tiles (i,j) with 0 ≤ i < j ≤ n-2, j-i ≥ 2.
  There are C(n-2,2) = (n-2)(n-3)/2 such tiles.

  The DIFFERENCE: δ_{n-2} \ δ_{n-3} = tiles involving vertex n-1
  = {(i, n-1) : 0 ≤ i ≤ n-3} = the rightmost column = n-2 tiles.

  So: δ_{n-2} = δ_{n-3} ∪ (rightmost column of n-2 tiles).

  The square interpretation:
  Take the (n-2) × (n-2) grid.
  Lower triangle (below diagonal): this IS δ_{n-3} (the overlap tiling).
  The rightmost column (n-2 cells): this IS the new layer from S322.
  The upper triangle (above diagonal): flip the (n-1)-tiling here.

  Lower triangle: C(n-2,2) = (n-2)(n-3)/2 cells
  Diagonal: n-2 cells
  Upper triangle: C(n-2,2) = (n-2)(n-3)/2 cells
  Total: (n-2)(n-3)/2 + (n-2) + (n-2)(n-3)/2 = (n-2)². ✓

  BUT: C(n-1,2) = (n-1)(n-2)/2 = (n-2)(n-3)/2 + (n-2)
  = lower triangle + diagonal = C(n-2,2) + (n-2).

  And C(n-2,2) = (n-2)(n-3)/2 = upper triangle.

  So: n-tiling fills lower triangle + diagonal = C(n-1,2) cells.
  (n-1)-tiling fills upper triangle = C(n-2,2) cells.
  Together: (n-2)² = the full square!

  THE SQUARE MATRIX is:
    M[r][c] = { n-tiling bit if r ≥ c (lower + diagonal)
              { (n-1)-tiling bit if r < c (upper, flipped)

Author: opus-2026-03-24-S325
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
banner("THE ARITHMETIC: δ_{n-2} + δ_{n-3} = (n-2)²")
# ============================================================

print("  %-4s %-8s %-8s %-10s %-8s %-8s" %
      ("n", "δ_{n-2}", "δ_{n-3}", "sum", "(n-2)²", "match"))
print("  " + "-"*50)

for n in range(3, 15):
    dn2 = comb(n-1, 2)  # tiles in n-tiling
    dn3 = comb(n-2, 2)  # tiles in (n-1)-tiling
    sq = (n-2) ** 2
    print("  %-4d %-8d %-8d %-10d %-8d %-8s" %
          (n, dn2, dn3, dn2 + dn3, sq, "✓" if dn2 + dn3 == sq else "✗"))

# ============================================================
banner("THE SQUARE AT n=5")
# ============================================================

n = 5
side = n - 2  # 3
print("  n=%d: side=%d, square=%d cells\n" % (n, side, side*side))

# The n=5 tiling has 6 tiles: (0,2),(0,3),(0,4),(1,3),(1,4),(2,4)
# The n=4 tiling has 3 tiles: (0,2),(0,3),(1,3)
# Together: 6 + 3 = 9 = 3² ✓

# Map tiles to square positions
# The (n-2)×(n-2) grid has rows 0..2, columns 0..2
# The n-tiling fills the lower triangle + diagonal

# Tile (i,j) in the n-tiling maps to grid position...
# We need a bijection between tiles and grid cells.

# THE NATURAL MAP:
# n-tiling tile (i, j) with j-i ≥ 2:
#   row = j - 2 (so j=2→row 0, j=3→row 1, j=4→row 2)
#   col = i     (so i=0→col 0, i=1→col 1, i=2→col 2)
#   This maps to the LOWER TRIANGLE + DIAGONAL of the 3×3 grid.

print("  n-TILING (δ_%d) → lower triangle + diagonal:\n" % (n-2))
grid_n = [['.' for _ in range(side)] for _ in range(side)]

base_set = set((i, i+1) for i in range(n-1))
all_p = [(i,j) for i in range(n) for j in range(i+1,n)]
tiles_n = [(i,j) for i,j in all_p if (i,j) not in base_set]
tiles_n.sort()

for idx, (lo, hi) in enumerate(tiles_n):
    r = hi - 2  # row
    c = lo       # column
    if r < side and c < side:
        grid_n[r][c] = chr(65 + idx)  # A, B, C, ...

for r in range(side):
    print("    " + "  ".join(grid_n[r]))

print("\n  Tiles: %s" % [(chr(65+i), t) for i, t in enumerate(tiles_n)])

# (n-1)-tiling tile (i,j) with j-i ≥ 2 and j ≤ n-2:
# row = i        (0-indexed)
# col = j - 2    (j=2→col 0, j=3→col 1)
# But we need to map to the UPPER TRIANGLE.

print("\n  (n-1)-TILING (δ_%d) → upper triangle:\n" % (n-3))

# The n-1 = 4 tiling has tiles (i,j) with 0 ≤ i < j ≤ 3, j-i ≥ 2:
# (0,2), (0,3), (1,3) = 3 tiles
base_set_4 = set((i, i+1) for i in range(n-2))
tiles_n1 = [(i,j) for i in range(n-1) for j in range(i+1,n-1) if (i,j) not in base_set_4 and j-i >= 2]
tiles_n1.sort()

grid_n1 = [['.' for _ in range(side)] for _ in range(side)]

# Map (n-1)-tiling to upper triangle of the 3×3 grid
# We flip: tile (i,j) → grid position (row=i, col=j-2) but in the UPPER part
# Upper triangle: r < c
for idx, (lo, hi) in enumerate(tiles_n1):
    # Map to upper triangle: row = lo, col = hi - 2
    r = lo
    c = hi - 2
    if r < side and c < side and r < c:  # upper triangle
        grid_n1[r][c] = chr(97 + idx)  # a, b, c, ...

for r in range(side):
    print("    " + "  ".join(grid_n1[r]))

print("\n  Tiles: %s" % [(chr(97+i), t) for i, t in enumerate(tiles_n1)])

# THE COMBINED SQUARE
print("\n  COMBINED SQUARE (n-tiling + (n-1)-tiling):\n")
grid_combined = [['.' for _ in range(side)] for _ in range(side)]

for r in range(side):
    for c in range(side):
        if grid_n[r][c] != '.':
            grid_combined[r][c] = grid_n[r][c]
        if grid_n1[r][c] != '.':
            grid_combined[r][c] = grid_n1[r][c]

for r in range(side):
    print("    " + "  ".join(grid_combined[r]))

# Count filled cells
n_filled = sum(1 for r in range(side) for c in range(side) if grid_combined[r][c] != '.')
print("\n  Filled: %d / %d cells" % (n_filled, side*side))

# ============================================================
banner("THE CHAIN: n=3 → n=4 → n=5 → n=6 → ...")
# ============================================================

print("""
  THE DIRECTED CHAIN OF TILINGS:

  Each consecutive pair (n, n-1) fills a (n-2)² square.
  The chain looks like:

  n=3 (δ_1, 1 tile) + n=2 (δ_0, 0 tiles) → 1² = 1 cell square
  n=4 (δ_2, 3 tiles) + n=3 (δ_1, 1 tile) → 2² = 4 cell square
  n=5 (δ_3, 6 tiles) + n=4 (δ_2, 3 tiles) → 3² = 9 cell square
  n=6 (δ_4, 10 tiles) + n=5 (δ_3, 6 tiles) → 4² = 16 cell square
  n=7 (δ_5, 15 tiles) + n=6 (δ_4, 10 tiles) → 5² = 25 cell square

  OVERLAPS: consecutive squares SHARE the smaller tiling.
  The n=4 tiling appears in BOTH the (4,3) square and the (5,4) square.
  This creates a CHAIN where each tiling links two squares.

  The chain structure:
    [1²] ← δ_1 → [2²] ← δ_2 → [3²] ← δ_3 → [4²] ← δ_4 → [5²] ...

  Each LINK (δ_k) is a tournament on k+2 vertices.
  Each SQUARE ((k+1)²) contains two consecutive tournaments.

  THE SQUARE MATRIX:
  The (n-2)×(n-2) matrix M has:
    M[r][c] = n-tiling bit for r ≥ c (lower + diagonal)
    M[r][c] = (n-1)-tiling bit for r < c (upper)

  This matrix encodes BOTH tournaments simultaneously.
  It has (n-2)² entries = C(n-1,2) + C(n-2,2) bits.

  WHAT DOES THIS MATRIX REPRESENT?

  The lower triangle encodes the n-tournament's non-base arcs.
  The upper triangle encodes the (n-1)-tournament's non-base arcs.
  The diagonal encodes the n-tournament's "rightmost column"
  = the new vertex's connections when going from n-1 to n.

  So: THE DIAGONAL IS THE INSERTION LAYER.
  M[r][r] = the bit that was added when extending from n-1 to n.
  There are n-2 diagonal bits = the insertion layer from S322.
""")

# Verify the diagonal = insertion layer
print("  DIAGONAL = INSERTION LAYER:\n")

for n in [5, 6, 7]:
    side = n - 2
    tiles_n_list = [(i,j) for i in range(n) for j in range(i+1,n) if j-i >= 2]
    tiles_n_list.sort()

    # The rightmost column tiles: (i, n-1) for i = 0..n-3
    insertion_layer = [(i, n-1) for i in range(n-2)]

    # In the grid, these map to: row = (n-1)-2 = n-3, col = i
    # Wait: (i, n-1) → row = n-1-2 = n-3, col = i.
    # But n-3 = side - 1. So these are ALL in the LAST ROW.
    # That's not the diagonal!

    # Let me reconsider the mapping.
    # Tile (lo, hi) → grid (row, col) with row = hi-2, col = lo.
    # The diagonal has row = col, so hi-2 = lo, so hi = lo+2.
    # These are the RANGE-2 TILES: (0,2), (1,3), (2,4), ...

    diagonal_tiles = [(lo, lo+2) for lo in range(n-2) if lo+2 <= n-1]

    print("  n=%d: diagonal tiles = %s (range-2 tiles)" % (n, diagonal_tiles))
    print("    insertion layer = %s (rightmost column)" % insertion_layer)
    print("    These are DIFFERENT: diagonal = range-2, insertion = rightmost column.")
    print()

print("""
  CORRECTION: The diagonal is NOT the insertion layer.
  The diagonal = RANGE-2 TILES (the atoms from S322).
  The insertion layer = RIGHTMOST COLUMN (all ranges involving vertex n-1).

  THE STRUCTURE OF THE SQUARE:
    Upper triangle: (n-1)-tiling bits
    Diagonal: range-2 tiles of the n-tiling (the ATOMS)
    Lower triangle (excl diagonal): higher-range tiles of n-tiling

  This decomposition separates the n-tiling into:
    ATOMS (diagonal): n-2 range-2 tiles = the most local interactions
    COMPOUNDS (lower triangle): C(n-2,2) higher-range tiles = longer-range interactions

  And the (n-1)-tiling fills the UPPER triangle:
    C(n-2,2) tiles = the previous level's compound tiles.

  THE ATOMS of level n ARE the COMPOUNDS of level n+1's upper triangle!
  Wait — that's not right either. Let me think about what the upper
  triangle actually represents.

  In the n-tiling: tiles (i,j) for j-i ≥ 2, 0 ≤ i < j ≤ n-1.
  Map: row = j-2, col = i.
  Lower triangle (row > col): j-2 > i → j > i+2 → range ≥ 3.
  Diagonal (row = col): j-2 = i → j = i+2 → range = 2.
  Upper part: row < col → j-2 < i → j < i+2 → range < 2 → impossible for tiles!

  So the n-tiling fills the lower triangle + diagonal EXACTLY.
  The upper triangle IS the (n-1)-tiling.

  THE KEY: the (n-1)-tiling's tiles (i,j) with j ≤ n-2 map to
  upper-triangle cells via some coordinate transformation.
  With the FLIP: (n-1)-tile (i,j) → grid (col, row) = (j-2, i)
  which gives row = i, col = j-2. For this to be upper: i < j-2,
  i.e., j > i+2, i.e., range ≥ 3. And for range = 2: i = j-2 → diagonal.

  So the (n-1)-tiling ALSO fills lower+diagonal of a grid, but
  when TRANSPOSED it fills the upper triangle.

  THE SQUARE IS:
    M = lower(n-tiling) + upper((n-1)-tiling transposed)
    = δ_{n-2} + δ_{n-3}^T
""")

# ============================================================
banner("THE SQUARE MATRIX AS A BINARY MATRIX")
# ============================================================

# Build actual square matrices at n=5

n = 5; side = n - 2  # 3×3 square
tiles_n5 = [(i,j) for i in range(n) for j in range(i+1,n) if j-i >= 2]
tiles_n5.sort()
tiles_n4 = [(i,j) for i in range(n-1) for j in range(i+1,n-1) if j-i >= 2]
tiles_n4.sort()

print("  n=5: 3×3 square matrix for ALL 64 × 8 = 512 tiling pairs:\n")
print("  n=5 tiling has %d bits, n=4 tiling has %d bits" % (len(tiles_n5), len(tiles_n4)))
print("  Square: %d × %d = %d cells ✓\n" % (side, side, side*side))

# For a specific (n5_bits, n4_bits) pair: build the square matrix
def build_square(n5_bits, n4_bits, n):
    side = n - 2
    M = np.zeros((side, side), dtype=int)

    # n-tiling → lower triangle + diagonal
    tiles_n_list = [(i,j) for i in range(n) for j in range(i+1,n) if j-i >= 2]
    tiles_n_list.sort()
    for k, (lo, hi) in enumerate(tiles_n_list):
        r = hi - 2; c = lo
        if r < side and c < side:
            M[r][c] = (n5_bits >> k) & 1

    # (n-1)-tiling → upper triangle (transposed)
    tiles_n1_list = [(i,j) for i in range(n-1) for j in range(i+1,n-1) if j-i >= 2]
    tiles_n1_list.sort()
    for k, (lo, hi) in enumerate(tiles_n1_list):
        # Transpose: (lo, hi) → grid (lo, hi-2) in the upper part
        r = lo; c = hi - 2
        if r < side and c < side and r < c:
            M[r][c] = (n4_bits >> k) & 1

    return M

# Show a few examples
for n5_bits in [0, 7, 42, 63]:
    for n4_bits in [0, 5, 7]:
        M = build_square(n5_bits, n4_bits, 5)
        det = int(round(np.linalg.det(M)))
        trace = int(np.trace(M))
        rank = np.linalg.matrix_rank(M)
        print("  n5=%s n4=%s → M = %s det=%d tr=%d rank=%d" %
              (format(n5_bits, '06b'), format(n4_bits, '03b'),
               M.flatten().tolist(), det, trace, rank))

# ============================================================
banner("PROPERTIES OF THE SQUARE MATRIX")
# ============================================================

# For all 64 × 8 = 512 pairs: compute determinant distribution
det_dist = Counter()
rank_dist = Counter()
trace_dist = Counter()

for n5_bits in range(64):
    for n4_bits in range(8):
        M = build_square(n5_bits, n4_bits, 5)
        det = int(round(np.linalg.det(M)))
        rank_val = np.linalg.matrix_rank(M)
        trace_val = int(np.trace(M))
        det_dist[det] += 1
        rank_dist[rank_val] += 1
        trace_dist[trace_val] += 1

print("  n=5: 512 square matrices (64 × 8 tiling pairs)\n")
print("  DETERMINANT distribution:")
for d in sorted(det_dist.keys()):
    pct = 100 * det_dist[d] / 512
    print("    det=%2d: %3d matrices (%.1f%%)" % (d, det_dist[d], pct))

print("\n  RANK distribution:")
for r in sorted(rank_dist.keys()):
    print("    rank=%d: %d matrices (%.1f%%)" % (r, rank_dist[r], 100*rank_dist[r]/512))

print("\n  TRACE distribution (= # range-2 tiles set to 1 in n-tiling):")
for t in sorted(trace_dist.keys()):
    print("    trace=%d: %d matrices (%.1f%%)" % (t, trace_dist[t], 100*trace_dist[t]/512))

# ============================================================
banner("THE CHAIN OF SQUARES")
# ============================================================

print("""
  THE CHAIN:

  Level:    3    4    5    6    7    8
  δ size:   1    3    6   10   15   21
  Square: 1×1  2×2  3×3  4×4  5×5  6×6

  Each square links consecutive levels:
    S(n) = δ_{n-2} (lower+diag) ⊕ δ_{n-3}^T (upper)

  The CHAIN:
    ... → δ_1 → [S(4)=2²] → δ_2 → [S(5)=3²] → δ_3 → [S(6)=4²] → δ_4 → ...

  Each δ_k appears as:
    - Lower+diagonal of S(k+2)
    - Upper (transposed) of S(k+3)

  So each tiling is shared between TWO consecutive squares.
  The chain has one degree of freedom per link (= per level).
  Total DOF for chain from n=3 to n=N:
    Σ_{k=1}^{N-2} k = C(N-1, 2)
  which IS the tiling of δ_{N-2}. The chain encodes the SAME bits
  as the largest tiling, just arranged differently.

  BUT: the square matrix M = δ_{n-2} ⊕ δ_{n-3}^T has ALGEBRAIC
  properties (determinant, rank, eigenvalues) that depend on
  BOTH tilings. This creates CORRELATIONS between consecutive levels
  that the individual tilings don't see.

  THE DETERMINANT: det(M) links the n and n-1 tilings algebraically.
  A chain of determinants det(M(3)), det(M(4)), det(M(5)), ...
  gives a SEQUENCE that encodes the inter-level structure.
""")

print("Done. opus-2026-03-24-S325")
