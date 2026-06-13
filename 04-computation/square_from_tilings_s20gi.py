#!/usr/bin/env python3
"""
square_from_tilings_s20gi.py — Two tilings tile a square
kind-pasteur-2026-03-25-S20gi

A tiling of "size n" = delta_{n-2}, has C(n-1,2) cells.
A tiling of "size n-1" = delta_{n-3}, has C(n-2,2) cells.

Together: C(n-1,2) + C(n-2,2) = (n-1)(n-2)/2 + (n-2)(n-3)/2
= (n-2)[(n-1)+(n-3)]/2 = (n-2)(2n-4)/2 = (n-2)^2.

An (n-2) x (n-2) SQUARE has (n-2)^2 cells!

So: delta_{n-2} (lower-left) + delta_{n-3} flipped (upper-right) = square_{n-2}.

The square encodes a PAIR of tournaments: one on n vertices, one on n-1 vertices.
The n-tournament and the (n-1)-tournament share the DIAGONAL of the square.

This creates a CHAIN: ...delta_{k-1} + delta_k tile square_k...
  square_1 = delta_0 + delta_1 = 0 + 1 = 1 cell
  square_2 = delta_1 + delta_2 = 1 + 3 = 4 cells
  square_3 = delta_2 + delta_3 = 3 + 6 = 9 cells
  square_k = delta_{k-1} + delta_k = C(k,2) + C(k+1,2) = k^2 cells

The DIRECTED CHAIN of tilings: each level k uses delta_k (lower-left)
paired with delta_{k-1} (upper-right, flipped). The shared diagonal
is the BRIDGE between consecutive levels.

Verify this geometrically and find the cell correspondence.
"""

import sys

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  TWO TILINGS TILE A SQUARE")
print("  kind-pasteur-2026-03-25-S20gi")
print("=" * 80)

for k in range(1, 8):
    n = k + 2  # tournament on n vertices uses delta_k
    nm1 = k + 1  # tournament on n-1 vertices uses delta_{k-1}

    lower = k * (k + 1) // 2  # C(k+1, 2) = delta_k cells
    upper = k * (k - 1) // 2  # C(k, 2) = delta_{k-1} cells
    total = lower + upper
    square = k * k

    print(f"\n  k={k} (n={n}):")
    print(f"    delta_{k} (lower-left): {lower} cells = C({k+1},2)")
    print(f"    delta_{k-1} (upper-right): {upper} cells = C({k},2)")
    print(f"    Total: {total} = {square} = {k}^2 = square_{k}")
    print(f"    Match: {total == square}")

    # Build the actual cell layout
    # The square is k x k, rows and cols indexed 0..k-1.
    # Lower-left triangle: cells (row, col) with row >= col (including diagonal)
    # Upper-right triangle: cells (row, col) with row < col
    # Lower triangle has k(k+1)/2 cells = C(k+1,2) = delta_k. ✓
    # Upper triangle has k(k-1)/2 cells = C(k,2) = delta_{k-1}. ✓

    if k <= 5:
        print(f"    Square layout ({k}x{k}):")
        for row in range(k-1, -1, -1):
            line = "    "
            for col in range(k):
                if row >= col:
                    # Lower-left: this is a cell of delta_k
                    # In the staircase, this corresponds to tile (?, ?)
                    line += f" L{row},{col}"
                elif row < col:
                    # Upper-right: this is a cell of delta_{k-1} (flipped)
                    line += f" U{row},{col}"
            print(line)

    # The DIAGONAL: cells (i, i) for i=0..k-1.
    # These k cells are part of the LOWER triangle (row >= col includes row = col).
    # The diagonal represents the base path of delta_k!
    # In the tournament on n=k+2 vertices, the base path has k+1 arcs... wait, k arcs.
    # Hmm, delta_k has C(k+1,2) cells and the lower triangle has k(k+1)/2 cells.
    # The diagonal has k cells. The strictly lower has k(k-1)/2 cells.
    # k + k(k-1)/2 = k(1 + (k-1)/2) = k(k+1)/2 = C(k+1,2). ✓

    # So: the lower triangle = diagonal (k cells) + strictly lower (C(k,2) cells).
    # The strictly lower part = delta_{k-1} (the OVERLAP from Mode B).
    # The diagonal = the new tiles added in the last Mode A step.

    # And the upper triangle = delta_{k-1} flipped = the companion tiling.

    # INTERPRETATION: The diagonal is SHARED information.
    # Below diagonal: the n-tournament's "old" tiles (delta_{k-1}).
    # Above diagonal: the (n-1)-tournament's tiles (also delta_{k-1}).
    # The diagonal itself: the n-tournament's "new" tiles (k bits at level k).

    print(f"    Diagonal ({k} cells): NEW tiles at level k (vertex {n-1}'s arcs)")
    print(f"    Below diagonal ({k*(k-1)//2} cells): delta_{k-1} = tournament on {n-1} verts")
    print(f"    Above diagonal ({k*(k-1)//2} cells): delta_{k-1} flipped = companion")

print(f"\n{'='*60}")
print("THE CHAIN STRUCTURE")
print(f"{'='*60}")
print("""
THE SQUARE DECOMPOSITION:

  square_k = delta_k (lower) + delta_{k-1} (upper, flipped)

  Lower triangle = diagonal + below-diagonal
                  = (k new bits) + delta_{k-1}
                  = level-k tiles + sub-tournament

  Upper triangle = delta_{k-1} flipped
                  = the COMPANION tiling at level k-1

THE DIRECTED CHAIN:
  ...-> delta_{k-1} -> square_k -> delta_k -> square_{k+1} -> delta_{k+1} -> ...

Each step: delta_k = delta_{k-1} + k new bits (the diagonal of square_k).
The companion tiling in the upper triangle is delta_{k-1} viewed from
the OTHER direction (flipped staircase).

THIS IS THE TOURNAMENT MATRIX:
  The (n-2) x (n-2) adjacency matrix of the inner tournament on {1,...,n-2}
  IS the square. The lower triangle is A[i][j] for i > j.
  The upper triangle is A[i][j] for i < j = 1 - A[j][i].
  The diagonal is A[i][i] = 0 (no self-arcs).

So: the square IS the adjacency matrix of the (n-2)-tournament.
    The lower triangle IS the tiling (delta_{k} for the n-tournament).
    The upper triangle IS the complement tiling (= (n-1)-tournament data).

The two tilings (size n and size n-1) together give the COMPLETE
adjacency matrix of the overlap sub-tournament. The diagonal is
the shared information — it's where the n-tournament adds its new vertex.
""")

print("DONE.")
print("=" * 80)
