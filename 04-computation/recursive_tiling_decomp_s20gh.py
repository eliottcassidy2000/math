#!/usr/bin/env python3
"""
recursive_tiling_decomp_s20gh.py — How tilings recursively decompose into tiles
kind-pasteur-2026-03-25-S20gh

The staircase delta_k has T(k) = k(k+1)/2 tiles.
It decomposes recursively:

  delta_k = delta_{k-1} + (k new tiles forming the k-th row/column)

The new tiles at level k: the bottom row y=k and left column x=k+2.
Wait, let me be precise with the explorer's convention.

In the explorer (n vertices, staircase delta_{n-2}):
  Tile (x,y) where 1 <= y <= n-2, y+2 <= x <= n.
  Row y has (n-1-y) tiles.

Recursive decomposition: delta_{n-2} -> delta_{n-3} + new tiles.
  delta_{n-3}: tiles (x,y) with 2 <= y <= n-2, y+2 <= x <= n-1.
    (shift: inner tournament on vertices {2,...,n-1})
  New tiles: those involving vertex 1 or vertex n.
    Row y=1: (x,1) for x=3,...,n — these are (n-2) tiles [bottom wiring]
    Column x=n: (n,y) for y=1,...,n-2 — these are (n-2) tiles [top wiring]
    But (n,1) is counted in BOTH! So new tiles = 2(n-2) - 1 = 2n-3.

Overlap of bottom and top = {(n,1)} = the APEX.
New = bottom + top - apex = (n-2) + (n-2) - 1 = 2n-3.
Inner = delta_{n-3} with C(n-2,2) - (something)... let me just count.

Actually: tiles of delta_{n-2} = C(n-1,2).
         tiles of delta_{n-3} = C(n-2,2).
         new tiles = C(n-1,2) - C(n-2,2) = n-2.

Wait: C(n-1,2) - C(n-2,2) = (n-1)(n-2)/2 - (n-2)(n-3)/2 = (n-2)[(n-1)-(n-3)]/2 = (n-2)*2/2 = n-2.

So: delta_{n-2} = delta_{n-3} + (n-2) new tiles. The Mode A recursion.

But in the THREE-PART decomposition (from the tiling model):
  delta_{n-2} = overlap(delta_{n-4}) + bottom(n-2) + top(n-2) + apex(1).
  overlap has C(n-2,2) - 2(n-3) - 1 = ? tiles.

Let me just enumerate all the decompositions precisely.
"""

import sys

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  RECURSIVE TILING DECOMPOSITION")
print("  kind-pasteur-2026-03-25-S20gh")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in range(4, 9):
    tiles = get_tiles(n)
    m = len(tiles)

    # Classify each tile
    bottom = [(x,y) for (x,y) in tiles if y == 1]  # row 1: arcs from vertex 1
    top = [(x,y) for (x,y) in tiles if x == n]     # column n: arcs from vertex n
    apex = [(x,y) for (x,y) in tiles if x == n and y == 1]  # (n,1)
    overlap = [(x,y) for (x,y) in tiles if y >= 2 and x <= n-1]  # inner staircase

    # Mode A: remove one "layer" (anti-diagonal or row+column)
    # Mode A removes tiles with x = n (top column): these are the n-2 tiles involving vertex n
    mode_a_new = [(x,y) for (x,y) in tiles if x == n]  # = top
    mode_a_inner = [(x,y) for (x,y) in tiles if x < n]  # = everything else

    # Mode B: remove BOTH vertex 1 AND vertex n
    mode_b_removed = [(x,y) for (x,y) in tiles if x == n or y == 1]
    mode_b_inner = [(x,y) for (x,y) in tiles if x < n and y > 1]  # = overlap

    # Three-part: bottom + top + apex + overlap (with bottom ∩ top = apex)
    bottom_only = [t for t in bottom if t not in top]
    top_only = [t for t in top if t not in bottom]

    print(f"\n  n = {n}, delta_{n-2}: {m} tiles = C({n-1},2)")
    print(f"    Bottom (y=1):  {len(bottom)} tiles = {bottom}")
    print(f"    Top (x={n}):    {len(top)} tiles = {top}")
    print(f"    Apex ({n},1):   {len(apex)} tile")
    print(f"    Overlap:       {len(overlap)} tiles = C({n-2},2) = {(n-2)*(n-3)//2}")
    print(f"    Check: {len(bottom)} + {len(top)} - {len(apex)} + {len(overlap)} = {len(bottom)+len(top)-len(apex)+len(overlap)} = {m}")

    print(f"\n    MODE A (remove vertex {n}):")
    print(f"      New:   {len(mode_a_new)} tiles (= {n}-2 = {n-2})")
    print(f"      Inner: {len(mode_a_inner)} tiles (= C({n-2},2) = {(n-2)*(n-3)//2})")
    print(f"      Inner IS delta_{n-3}")

    print(f"\n    MODE B (remove vertices 1 and {n}):")
    print(f"      Removed: {len(mode_b_removed)} tiles (= 2*{n-2}-1 = {2*(n-2)-1})")
    print(f"      Inner:   {len(mode_b_inner)} tiles (= C({n-3},2) = {(n-3)*(n-4)//2})")
    print(f"      Inner IS delta_{n-4}")

    # FRACTAL VIEW: the full recursion tree
    # delta_k decomposes as delta_{k-1} + k new tiles (Mode A)
    # or as delta_{k-2} + (2k-1) border tiles (Mode B)
    print(f"\n    FRACTAL DECOMPOSITION (Mode A chain):")
    k = n - 2
    total = 0
    while k >= 0:
        new_tiles = k  # Mode A adds k tiles at level k
        inner = k * (k-1) // 2  # delta_{k-1} has this many
        print(f"      Level {k}: delta_{k} = delta_{k-1} + {new_tiles} new tiles")
        total += new_tiles
        k -= 1
    print(f"      Total tiles: {total} = {m}")

    # Each tile is delta_0 (a single binary cell)
    # The recursion: delta_k = delta_{k-1} + row_k
    # row_k = the k tiles in the new row/column
    # These k tiles are INDEPENDENT of delta_{k-1} (free binary choices)
    print(f"\n    INFORMATION CONTENT:")
    print(f"      delta_{n-2} = {m} free bits = {' + '.join(str(k) for k in range(n-2, 0, -1))} + 0")
    print(f"      = sum 1..{n-2} = C({n-1},2) = {m}")

    # The TILE at each level is a CHOICE: forward or backward arc.
    # At level k: the k new tiles each connect a new vertex to k existing vertices.
    # The TOURNAMENT on k+2 vertices is determined by:
    #   tournament on k+1 vertices (delta_{k-1}) + k new arc choices.

print(f"\n{'='*60}")
print("THE RECURSIVE VIEW")
print(f"{'='*60}")
print("""
A TILING IS A NESTED SEQUENCE OF CHOICES:

  Level 0: delta_0 = 0 tiles (tournament on 2 vertices = trivial)
  Level 1: delta_1 = 1 tile (tournament on 3 vertices = 1 choice)
  Level 2: delta_2 = 1 + 2 = 3 tiles (4 vertices = 2 new choices)
  Level 3: delta_3 = 3 + 3 = 6 tiles (5 vertices = 3 new choices)
  Level k: delta_k = delta_{k-1} + k tiles (k+2 vertices = k new choices)

Each level adds one vertex and k arcs connecting it to all existing vertices.
The k arcs are the k new tiles. They are FREE binary choices,
independent of the previous levels.

THIS IS THE FRACTAL CODEC: the tiling is literally a sequence of
independent layers, each adding one vertex with k free bits.

The "tile" at position (x,y) in the staircase belongs to level x-y-1:
  (x,y) is in the (x-y-1)-th layer of new tiles.
  Its "skip" = x-y-1 = the level at which it was added.

So: skip=0 tiles (base path) were "added" at level 0 (trivially fixed).
    skip=1 tiles were added at level 1 (first non-trivial choices).
    skip=k tiles were added at level k.

The SKIP HIERARCHY corresponds exactly to the RECURSION DEPTH:
  High-skip tiles were added EARLY (deep in the recursion).
  Low-skip tiles were added LATE (near the surface).
  skip=1 tiles are the MOST RECENT additions (the "hypotenuse").

This explains WHY skip=1 has low SL%: the most recent tiles are
the least constrained by the existing structure, so flipping them
is most likely to change the isomorphism class.
""")

print("DONE.")
print("=" * 80)
