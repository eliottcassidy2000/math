#!/usr/bin/env python3
"""
wiggly_pairing_n7_s20ff.py — Grid-symmetry pairing of wiggly classes at n=4..7
kind-pasteur-2026-03-24-S20ff

At n=6: 10 wiggly classes form 4 identical pairs + 2 fixed points under
grid symmetry (x,y) -> (n+1-y, n+1-x). Jaccard = 1.000 for each pair.
The quotient = 6 = C(n-2,2) = tiles at n-1.

Does this pattern extend to n=7? The 15 tiles should form pairs + fixed
under grid symmetry, and the quotient should give C(n-2,2) = C(5,2) = 10
effective classes = tiles at n-1 = 6.

Wait: C(5,2) = 10. But the n=6 quotient has 6 = C(4,2). So the quotient
at n=7 should have C(5,2) = 10? Let me check.

Grid symmetry at n=7: (x,y) -> (8-y, 8-x). Fixed: x+y = 8.
Tiles with x+y=8: (7,1), (6,2), (5,3) — 3 fixed tiles.
Paired: (7-3=4 non-diagonal)/2... let me compute directly.
"""

import sys
from math import factorial
from itertools import permutations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  GRID-SYMMETRY PAIRING AT n=4..7")
print("  kind-pasteur-2026-03-24-S20ff")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in range(4, 8):
    TILES = get_tiles(n)
    m = len(TILES)
    labels = [chr(65+i) if i < 26 else chr(65+i-26)+str(i//26) for i in range(m)]

    print(f"\n{'#'*60}")
    print(f"  n = {n}, m = {m} tiles")
    print(f"{'#'*60}")

    # Grid symmetry map: (x,y) -> (n+1-y, n+1-x)
    tile_idx = {(t[0],t[1]): i for i, t in enumerate(TILES)}
    grid_map = []
    for i, (x, y) in enumerate(TILES):
        gx, gy = n+1-y, n+1-x
        gi = tile_idx.get((gx, gy), -1)
        grid_map.append(gi)

    # Fixed points (on anti-diagonal x+y = n+1)
    fixed = [i for i in range(m) if grid_map[i] == i]
    # Pairs
    pairs = []
    seen = set()
    for i in range(m):
        if i in seen: continue
        gi = grid_map[i]
        if gi == i:
            seen.add(i)
        elif gi >= 0:
            pairs.append((i, gi))
            seen.add(i)
            seen.add(gi)

    n_eff = len(fixed) + len(pairs)  # effective wiggly classes after quotient

    print(f"  Grid symmetry: (x,y) -> ({n+1}-y, {n+1}-x)")
    print(f"  Fixed tiles (anti-diagonal x+y={n+1}): {len(fixed)}")
    for i in fixed:
        print(f"    {labels[i]}: ({TILES[i][0]},{TILES[i][1]})")
    print(f"  Paired tiles: {len(pairs)} pairs")
    for i, j in pairs:
        print(f"    {{{labels[i]}: ({TILES[i][0]},{TILES[i][1]}), {labels[j]}: ({TILES[j][0]},{TILES[j][1]})}}")
    print(f"  Effective classes (quotient): {n_eff}")
    print(f"  C(n-2, 2) = C({n-2}, 2) = {(n-2)*(n-3)//2} = tiles at n-1")
    print(f"  Match? {n_eff == (n-2)*(n-3)//2}")

    # For n <= 7: compute Jaccard similarity to verify pairing
    if n <= 7:
        t0 = time.time()
        N = n
        VERTS = list(range(n, 0, -1))
        all_perms = list(permutations(range(N)))

        def bits_to_adj(bits):
            A = [[0]*N for _ in range(N)]
            for k in range(N-1): A[k][k+1] = 1
            for i in range(m):
                x, y = TILES[i]
                xi, yi = VERTS.index(x), VERTS.index(y)
                if bits[i] == 0: A[xi][yi] = 1
                else: A[yi][xi] = 1
            return A

        def adj_complement(A):
            return [[1-A[i][j] if i!=j else 0 for j in range(N)] for i in range(N)]

        def canonicalize(A):
            best = None
            for p in all_perms:
                s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
                if best is None or s < best: best = s
            return best

        # Build canonical map
        canon_map = {}
        for mask in range(1 << m):
            bits = [(mask >> k) & 1 for k in range(m)]
            A = bits_to_adj(bits)
            canon_map[mask] = canonicalize(A)

        # True complement merging
        comp_cn = {}
        for cn in set(canon_map.values()):
            for mask, c in canon_map.items():
                if c == cn:
                    A = bits_to_adj([(mask >> k) & 1 for k in range(m)])
                    comp_cn[cn] = canonicalize(adj_complement(A))
                    break
        merged_map = {cn: min(cn, comp_cn.get(cn, cn)) for cn in set(canon_map.values())}

        # Per wiggly class: which merged edges
        class_edges = [set() for _ in range(m)]
        class_sl = [0 for _ in range(m)]
        for mask in range(1 << m):
            mcn1 = merged_map[canon_map[mask]]
            for wi in range(m):
                mcn2 = merged_map[canon_map[mask ^ (1 << wi)]]
                if mcn1 != mcn2:
                    class_edges[wi].add((min(mcn1,mcn2), max(mcn1,mcn2)))
                else:
                    class_sl[wi] += 1

        # Verify Jaccard = 1.0 for grid-symmetry pairs
        print(f"\n  JACCARD VERIFICATION (paired classes):")
        all_perfect = True
        for i, j in pairs:
            common = len(class_edges[i] & class_edges[j])
            union = len(class_edges[i] | class_edges[j])
            jaccard = common / union if union > 0 else 0
            perfect = abs(jaccard - 1.0) < 0.0001
            if not perfect: all_perfect = False
            print(f"    {labels[i]}-{labels[j]}: J={jaccard:.4f} {'PERFECT' if perfect else 'NOT PERFECT'}, |edges|={len(class_edges[i])},{len(class_edges[j])}")

        print(f"\n  ALL pairs have Jaccard = 1.0? {all_perfect}")

        # SL rates for paired classes
        print(f"\n  SL RATES (paired classes should match):")
        for i, j in pairs:
            sl_i = class_sl[i] / (2**m) * 100
            sl_j = class_sl[j] / (2**m) * 100
            print(f"    {labels[i]}: {sl_i:.2f}%, {labels[j]}: {sl_j:.2f}%, match: {abs(sl_i-sl_j)<0.01}")

        for i in fixed:
            sl_i = class_sl[i] / (2**m) * 100
            print(f"    {labels[i]} (fixed): {sl_i:.2f}%")

        # Effective edge counts
        print(f"\n  EFFECTIVE EDGE COUNTS (after quotient):")
        eff_edges = set()
        for i in fixed:
            eff_edges.update(class_edges[i])
        for i, j in pairs:
            eff_edges.update(class_edges[i])  # = class_edges[j]

        all_merged = set()
        for edges in class_edges:
            all_merged.update(edges)

        print(f"    Total metagraph edges: {len(all_merged)}")
        print(f"    Edges reachable by quotient classes: {len(eff_edges)} (should be same)")

        elapsed = time.time() - t0
        print(f"\n  Time: {elapsed:.1f}s")

print("\n" + "=" * 60)
print("SYNTHESIS: GRID-SYMMETRY QUOTIENT PATTERN")
print("=" * 60)
print("""
At every n=4..7:
  m = C(n-1,2) tiles
  Grid symmetry pairs: (m - anti-diag) / 2 pairs + anti-diag fixed
  Quotient = (m + anti-diag) / 2 effective classes

  Anti-diagonal count = floor((n-1)/2)
  Quotient size = (C(n-1,2) + floor((n-1)/2)) / 2 = C(n-2,2) + floor((n-1)/2) ... ?

  n=4: m=3, fixed=1, pairs=1, quotient=2. C(2,2)=1. NOT matching C(n-2,2).
  n=5: m=6, fixed=2, pairs=2, quotient=4. C(3,2)=3.
  n=6: m=10, fixed=2, pairs=4, quotient=6. C(4,2)=6. MATCH!
  n=7: m=15, fixed=3, pairs=6, quotient=9. C(5,2)=10.

  Actually: quotient = (m + fixed) / 2 = (C(n-1,2) + floor((n-1)/2)) / 2
  This is EXACTLY the blue exponent! = dimension of grid-sym tiling space.
  Grid-sym tilings = 2^quotient.

  The quotient of wiggly classes IS the free parameter space of blue tilings.
""")

print("DONE.")
print("=" * 80)
