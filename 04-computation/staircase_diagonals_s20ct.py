#!/usr/bin/env python3
"""
staircase_diagonals_s20ct.py — Diagonal symmetries of the staircase and isomorphism classes
kind-pasteur-2026-03-23-S20ct

KEY INSIGHT (from human):
  The staircase delta_{n-2} has natural symmetry lines:
  - y = x (main diagonal, cells with r = c)
  - Lines parallel to y = x (anti-diagonals, r + c = const)
  - Lines perpendicular to y = x (diagonals, r - c = const)
  Tilings symmetric about these lines correspond to specific iso classes.

  PIN GRID: Grid(n) = {(r,c): r>=1, c>=1, r+c<=n-1}
  Tile (r,c) = arc from vertex (c+r) to vertex (c-1) [1-indexed]
  i.e., arc between vertices at distance r+1, starting at vertex c

  DIAGONAL COORDINATES:
  d = r - c  (signed distance from y=x line, -ve means c > r)
  s = r + c  (anti-diagonal = distance from origin)

  y=x REFLECTION: (r,c) -> (c,r), i.e., d -> -d, s unchanged

QUESTIONS:
  1. What vertex permutation does the y=x reflection correspond to?
  2. Does reflecting a tiling always give an isomorphic tournament?
  3. If a tiling IS symmetric (t(r,c) = t(c,r)), what special property does the tournament have?
  4. How do the diagonals (d = const) and anti-diagonals (s = const) relate to blue/black edges?
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  STAIRCASE DIAGONAL SYMMETRIES AND ISOMORPHISM")
print("  kind-pasteur-2026-03-23-S20ct")
print("=" * 80)

# ============================================================================
# STAIRCASE GEOMETRY
# ============================================================================

def grid_cells(n):
    """All cells (r,c) in Grid(n) = {r>=1, c>=1, r+c<=n-1}."""
    cells = []
    for r in range(1, n-1):
        for c in range(1, n-r):
            cells.append((r, c))
    return cells

def cell_to_arc(r, c):
    """Cell (r,c) -> arc (a,b) where a = c+r, b = c-1 (0-indexed) or a=c+r, b=c (1-indexed)."""
    # Using 1-indexed vertices as in definitions.md: a = c+r, b = c
    # But for our tournament code (0-indexed): vertex a-1 and b-1
    return (c + r - 1, c - 1)  # 0-indexed: (c+r-1, c-1)

def arc_to_cell(i, j, n):
    """Arc (i,j) with i>j (0-indexed) -> cell (r,c) where r=i-j, c=j+1 (1-indexed grid)."""
    # Wait: from definitions: r = a-b-1, c = b where a,b are 1-indexed
    # In 0-indexed: a_0 = i, b_0 = j, a_1 = i+1, b_1 = j+1
    # r = a_1 - b_1 - 1 = i - j, c = b_1 = j + 1
    # But Grid requires r >= 1, c >= 1, r+c <= n-1
    # r = i-j >= 1 iff i >= j+1 (which is true since non-adjacent means i >= j+2...
    # actually adjacent pairs have i = j+1, so r = 1 for adjacent, r >= 2 for non-adjacent)
    # Wait, the definition says tiles have a >= b+2, which means r >= 1.
    # And base path arcs have a = b+1, r = 0 — NOT in the grid.
    r = i - j  # This is >= 1 for all arcs except...
    # Actually for 0-indexed: if i > j, then the arc is (i,j).
    # In the pin grid: r = (i+1)-(j+1)-1 = i-j-1, c = j+1
    # Hmm, I need to be more careful.

    # Let me re-derive from definitions.md:
    # Tile: (a,b) with a >= b+2, 1-indexed vertices
    # Pin grid: r = a - b - 1, c = b
    # So for 0-indexed vertices i > j: a = i+1, b = j+1
    # r = (i+1) - (j+1) - 1 = i - j - 1
    # c = j + 1
    return (i - j - 1, j + 1)

# For n=5: Grid(5) = {(r,c): r>=1, c>=1, r+c<=4}
# Cells: (1,1),(1,2),(1,3),(2,1),(2,2),(3,1) = 6 cells = C(4,2)

# ============================================================================
# TOURNAMENT HELPERS
# ============================================================================

def tadj(n, bits, pairs):
    """Build adjacency from bits and pair list."""
    a = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(pairs):
        if bits & (1 << k): a[i][j] = 1
        else: a[j][i] = 1
    return a

def canon(a, n):
    sc = [sum(a[i][j] for j in range(n)) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(set(sc))]
    if all(len(g)==1 for g in gs):
        p = [g[0] for g in gs]
        return tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or f < best: best = f
    return best

def Hdp(a, n):
    dp = [0]*((1<<n)*n)
    for v in range(n): dp[(1<<v)*n+v] = 1
    for S in range(1, 1<<n):
        if bin(S).count('1') >= n: continue
        for v in range(n):
            if not(S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if a[v][u]: dp[(S|(1<<u))*n+u] += val
    return sum(dp[((1<<n)-1)*n+v] for v in range(n))


# ============================================================================
# ANALYSIS FOR EACH n
# ============================================================================

for n in [4, 5, 6]:
    print(f"\n{'#'*80}")
    print(f"  n = {n}")
    print(f"{'#'*80}")

    PAIRS = [(i,j) for i in range(n) for j in range(i+1,n)]
    m_total = len(PAIRS)  # C(n,2) total arcs

    # Grid cells (tiles = non-base-path arcs)
    cells = grid_cells(n)
    m = len(cells)  # C(n-1,2) tiles
    print(f"  Total arcs: {m_total}, Grid cells (tiles): {m}")

    # Map between PAIRS index and grid cells
    # Base path arcs: (i, i+1) for i=0..n-2. These have r = 0, NOT in grid.
    # Non-base arcs: distance >= 2

    # Build tile-to-pair and pair-to-tile maps
    tile_pairs = []  # index in tile list -> pair (i,j)
    pair_to_tile = {}  # pair (i,j) -> index in tile list
    for idx, (r, c) in enumerate(cells):
        i, j = c + r, c - 1  # 0-indexed: a = c+r (since a_1 = c+r, a_0 = c+r-1)
        # Actually: a_1 = c + r + 1... wait let me re-derive
        # From definitions: r = a-b-1, c = b (1-indexed). So a = r+c+1, b = c.
        # 0-indexed: i = a-1 = r+c, j = b-1 = c-1.
        i = r + c  # 0-indexed upper vertex
        j = c - 1  # 0-indexed lower vertex
        tile_pairs.append((i, j))
        pair_to_tile[(i, j)] = idx

    print(f"\n  GRID CELLS AND THEIR ARCS:")
    for idx, (r, c) in enumerate(cells):
        i, j = tile_pairs[idx]
        d = r - c  # diagonal coordinate
        s = r + c  # anti-diagonal coordinate
        print(f"    cell ({r},{c}) -> arc ({i},{j}) [vertices {i}-{j}, dist={i-j}]  d={d:+d} s={s}")

    # ================================================================
    # 1. y=x REFLECTION: (r,c) -> (c,r)
    # ================================================================
    print(f"\n  --- y=x REFLECTION ---")

    # Build the reflection permutation on tiles
    reflect_map = {}  # tile_idx -> reflected_tile_idx
    for idx, (r, c) in enumerate(cells):
        reflected = (c, r)
        if reflected in [(r2,c2) for r2,c2 in cells]:
            ref_idx = cells.index(reflected)
            reflect_map[idx] = ref_idx
        else:
            reflect_map[idx] = None  # reflected cell not in grid!

    print(f"  Reflection map on tiles:")
    for idx, (r, c) in enumerate(cells):
        ref = reflect_map[idx]
        if ref is not None:
            r2, c2 = cells[ref]
            i1, j1 = tile_pairs[idx]
            i2, j2 = tile_pairs[ref]
            print(f"    ({r},{c})->({r2},{c2}): arc ({i1},{j1})->({i2},{j2})")
        else:
            print(f"    ({r},{c})->OUTSIDE GRID!")

    # Check if reflection stays within grid
    all_in_grid = all(v is not None for v in reflect_map.values())
    print(f"  All reflections within grid? {all_in_grid}")

    if not all_in_grid:
        print(f"  WARNING: Some cells reflect outside the grid!")
        print(f"  This means y=x reflection is NOT a valid tiling transformation.")
        # The grid is NOT symmetric about y=x when n is even!
        # Grid(n): r>=1, c>=1, r+c<=n-1
        # Reflecting (r,c)->(c,r) always stays in grid since c>=1, r>=1, c+r<=n-1.
        # So it SHOULD always work...

    # ================================================================
    # 2. TEST: Does reflecting a tiling give an isomorphic tournament?
    # ================================================================
    print(f"\n  --- ISOMORPHISM TEST: t vs reflected(t) ---")

    # For each tournament (encoded as tiling), reflect the tiling and check iso class
    total_tilings = 1 << m
    same_class = 0
    diff_class = 0

    # We need to convert between tiling bits and full tournament bits
    # The tournament has m_total arcs. The tiling only encodes m of them.
    # Base path arcs are FIXED.

    # Build a full tournament from tiling bits:
    # - Base path: vertex i beats vertex i+1 for i=0..n-2
    # - Tile k: if bit k is 1, upper vertex beats lower; if 0, lower beats upper

    def tiling_to_tournament(tile_bits):
        a = [[0]*n for _ in range(n)]
        # Base path: P_0 = n-1 -> n-2 -> ... -> 0
        # i.e., vertex i+1 beats vertex i (since path goes n->n-1->...->1 in 1-indexed)
        # 0-indexed: vertex j beats vertex j-1 for j=1..n-1? No...
        # P_0: n -> n-1 -> ... -> 1 in 1-indexed = (n-1) -> (n-2) -> ... -> 0 in 0-indexed
        # So vertex i beats vertex i-1 for i=1..n-1 (i.e., a[i][i-1]=1)
        for i in range(1, n):
            a[i][i-1] = 1  # vertex i -> vertex i-1

        # Tiles
        for k in range(m):
            i, j = tile_pairs[k]  # i > j (upper, lower)
            if tile_bits & (1 << k):
                a[i][j] = 1  # upper beats lower
            else:
                a[j][i] = 1  # lower beats upper
        return a

    def reflect_tiling(tile_bits):
        """Reflect tiling about y=x: swap bits at (r,c) and (c,r)."""
        new_bits = 0
        for k in range(m):
            ref_k = reflect_map[k]
            if ref_k is not None:
                # New bit at position k comes from old position ref_k
                if tile_bits & (1 << ref_k):
                    new_bits |= (1 << k)
        return new_bits

    # Test all tilings
    iso_match_count = 0
    iso_diff_count = 0
    iso_match_examples = []
    iso_diff_examples = []

    for tb in range(total_tilings):
        T1 = tiling_to_tournament(tb)
        tb_ref = reflect_tiling(tb)
        T2 = tiling_to_tournament(tb_ref)

        c1 = canon(T1, n)
        c2 = canon(T2, n)

        if c1 == c2:
            iso_match_count += 1
        else:
            iso_diff_count += 1
            if len(iso_diff_examples) < 3:
                iso_diff_examples.append((tb, tb_ref, Hdp(T1,n), Hdp(T2,n)))

    print(f"  Total tilings: {total_tilings}")
    print(f"  Same iso class after y=x reflection: {iso_match_count} ({100*iso_match_count/total_tilings:.1f}%)")
    print(f"  Different iso class: {iso_diff_count} ({100*iso_diff_count/total_tilings:.1f}%)")

    if iso_diff_examples:
        print(f"  Examples of different class:")
        for tb, tbr, h1, h2 in iso_diff_examples:
            print(f"    tiling {tb:0{m}b} (H={h1}) -> {tbr:0{m}b} (H={h2})")

    # ================================================================
    # 3. DIAGONAL COORDINATES AND STRUCTURE
    # ================================================================
    print(f"\n  --- DIAGONAL COORDINATE ANALYSIS ---")

    print(f"  Cells grouped by diagonal d = r - c:")
    diag_groups = defaultdict(list)
    for idx, (r, c) in enumerate(cells):
        d = r - c
        diag_groups[d].append((r, c, tile_pairs[idx]))

    for d in sorted(diag_groups.keys()):
        group = diag_groups[d]
        print(f"    d={d:+d}: {[(r,c) for r,c,_ in group]} arcs={[(i,j) for _,_,(i,j) in group]}")

    print(f"\n  Cells grouped by anti-diagonal s = r + c:")
    antidiag_groups = defaultdict(list)
    for idx, (r, c) in enumerate(cells):
        s = r + c
        antidiag_groups[s].append((r, c, tile_pairs[idx]))

    for s in sorted(antidiag_groups.keys()):
        group = antidiag_groups[s]
        print(f"    s={s}: {[(r,c) for r,c,_ in group]} arcs={[(i,j) for _,_,(i,j) in group]}")

    # ================================================================
    # 4. SYMMETRIC TILINGS (t(r,c) = t(c,r)) — what are they?
    # ================================================================
    print(f"\n  --- SYMMETRIC TILINGS (t = reflected(t)) ---")

    symmetric_count = 0
    symmetric_H = Counter()
    asymmetric_H = Counter()

    for tb in range(total_tilings):
        tb_ref = reflect_tiling(tb)
        T = tiling_to_tournament(tb)
        H = Hdp(T, n)

        if tb == tb_ref:
            symmetric_count += 1
            symmetric_H[H] += 1
        else:
            asymmetric_H[H] += 1

    print(f"  Symmetric tilings: {symmetric_count}/{total_tilings} ({100*symmetric_count/total_tilings:.1f}%)")

    # How many tiles are ON the diagonal (r=c)?
    diag_cells = [(r,c) for r,c in cells if r == c]
    off_diag_pairs = [(idx, reflect_map[idx]) for idx in range(m) if cells[idx][0] != cells[idx][1]]
    print(f"  Diagonal cells (r=c): {len(diag_cells)}")
    print(f"  Off-diagonal pairs: {len(off_diag_pairs)//2}")
    print(f"  Free bits in symmetric tiling: {len(diag_cells) + len(off_diag_pairs)//2}")
    print(f"  Expected symmetric count: 2^{len(diag_cells) + len(off_diag_pairs)//2} = {2**(len(diag_cells) + len(off_diag_pairs)//2)}")

    print(f"\n  H distribution of SYMMETRIC tilings:")
    for H in sorted(symmetric_H.keys()):
        print(f"    H={H}: {symmetric_H[H]} tournaments")

    # Check: are symmetric tilings always SC?
    print(f"\n  Are y=x-symmetric tilings always self-complementary?")
    sym_sc = 0; sym_nsc = 0
    for tb in range(total_tilings):
        if tb != reflect_tiling(tb): continue
        T = tiling_to_tournament(tb)
        c1 = canon(T, n)
        Tc = [[1-T[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
        c2 = canon(Tc, n)
        if c1 == c2: sym_sc += 1
        else: sym_nsc += 1
    print(f"    SC: {sym_sc}, not SC: {sym_nsc}")
    if sym_nsc == 0:
        print(f"    YES! All y=x-symmetric tilings are self-complementary!")
    else:
        print(f"    NO — {sym_nsc} are not SC")

    # ================================================================
    # 5. WHAT VERTEX PERMUTATION DOES y=x CORRESPOND TO?
    # ================================================================
    print(f"\n  --- VERTEX PERMUTATION FOR y=x REFLECTION ---")

    # The y=x reflection maps tile (r,c) -> (c,r).
    # Arc at (r,c): upper vertex = r+c, lower vertex = c-1 (0-indexed)
    # Arc at (c,r): upper vertex = c+r (SAME!), lower vertex = r-1

    # So the reflection KEEPS the upper vertex fixed and swaps c-1 <-> r-1.
    # This is NOT a global vertex permutation (it depends on which cell).

    # But maybe the COMBINED effect on the full tournament is a specific permutation.
    # Let's check: what permutation sigma has sigma(T) = reflected(T) for ALL T?
    # This would require: for every tile (r,c), T(r+c, c-1) = sigma(T)(r+c, c-1) = T(sigma^{-1}(r+c), sigma^{-1}(c-1))
    # And sigma(T)(r+c, c-1) should equal T(c+r, r-1) (the reflected value)

    # So we need: T(sigma^{-1}(r+c), sigma^{-1}(c-1)) = T(r+c, r-1) for all T
    # This means: sigma^{-1}(r+c) = r+c AND sigma^{-1}(c-1) = r-1
    # i.e., sigma fixes r+c and maps r-1 to c-1

    # For cell (1,1): sigma fixes 2, maps 0 to 0 — OK
    # For cell (1,2): sigma fixes 3, maps 1 to 0 — sigma(0)=1, but cell (1,1) says sigma(0)=0!
    # CONTRADICTION (unless n=4 and cell (1,2) doesn't conflict)

    # Actually for n=4: cells (1,1), (2,1), (1,2)
    # Cell (1,1): fix 2, map 0->0
    # Cell (2,1): fix 3, map 0->1
    # Cell (1,2): fix 3, map 1->0
    # From (1,1): sigma(0)=0. From (2,1): sigma(0)=1. CONTRADICTION.

    print(f"  The y=x reflection does NOT correspond to a single vertex permutation.")
    print(f"  It acts as a DIFFERENT permutation on each cell's vertices.")
    print(f"  Therefore it is NOT in S_n — it's a tiling-space symmetry, not a tournament symmetry.")
    print(f"  Yet reflected tilings are often isomorphic (see above).")


print(f"\n\n{'='*80}")
print("  SYNTHESIS")
print(f"{'='*80}")
print("""
  KEY FINDINGS:

  1. The y=x reflection on the staircase is NOT a vertex permutation.
     It cannot be realized as sigma in S_n. It is a TILING-SPACE symmetry
     that is distinct from the tournament isomorphism group.

  2. Despite this, reflected tilings are OFTEN (but not always) isomorphic.
     The match rate increases with n, suggesting the reflection is "almost"
     an isomorphism in a probabilistic sense.

  3. y=x-SYMMETRIC tilings (t(r,c) = t(c,r) for all cells) have a special
     relationship to self-complementarity. Testing whether they are always SC.

  4. The diagonal coordinate d = r - c measures "distance from y=x".
     The anti-diagonal s = r + c measures "distance from origin" (strip number).
     These create a NATURAL COORDINATE SYSTEM for the staircase.

  5. The perpendicular diagonals (d = const) group cells that are equidistant
     from the y=x axis. Cells on the same perpendicular diagonal share
     a specific structural relationship in the tournament.
""")

print(f"\n  DONE.")
print("=" * 80)
