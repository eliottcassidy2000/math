#!/usr/bin/env python3
"""
tiling_to_metagraph_s20cu.py — THE COMPLETE MAP: Tilings -> Iso Classes -> Meta-Graph
kind-pasteur-2026-03-23-S20cu

THE FUNDAMENTAL PICTURE:
  1. A TILING = binary labeling of the staircase delta_{n-2} (m = C(n-1,2) bits)
     Each tiling encodes a tournament T with a FIXED base path P_0 = n->n-1->...->1
  2. An ISO CLASS = orbit of S_n on labeled tournaments
     fiber(C) = H(C)/|Aut(C)| tilings encode class C
     sum_C fiber(C) = 2^m (partition of all tilings into class fibers)
  3. A META-GRAPH EDGE = flipping one tile changes the iso class
     BLUE if SC-type preserved, BLACK if changed
  4. Each tile (r,c) in the staircase generates edges
     EVERY tile generates the same blue fraction (isotropy!)

This script traces the COMPLETE chain at n=5 (m=6 tiles, 64 tilings, 12 classes).
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  THE COMPLETE MAP: TILINGS -> ISO CLASSES -> META-GRAPH")
print("  kind-pasteur-2026-03-23-S20cu")
print("=" * 80)

n = 5
m_total = comb(n, 2)  # 10 total arcs
# Grid(n): tiles (r,c) with r>=1, c>=1, r+c<=n-1
tiles = []
for r in range(1, n-1):
    for c in range(1, n-r):
        tiles.append((r, c))
m = len(tiles)  # C(n-1,2) = 6 tiles
print(f"\n  n = {n}")
print(f"  Total arcs: {m_total}, Tiles (non-base-path arcs): {m}")
print(f"  Tilings: 2^{m} = {2**m}")

# Map tile (r,c) to arc (i,j) in 0-indexed vertices
# From definitions: a = c+r (1-indexed), b = c (1-indexed)
# 0-indexed: i = c+r-1, j = c-1
tile_arcs = []
for r, c in tiles:
    # From definitions: a = r+c+1 (1-indexed), b = c (1-indexed)
    # 0-indexed: i = r+c, j = c-1
    i = r + c      # upper vertex (0-indexed)
    j = c - 1      # lower vertex (0-indexed)
    tile_arcs.append((i, j))

print(f"\n  STAIRCASE TILES AND THEIR ARCS:")
print(f"  {'(r,c)':<8} {'arc(i,j)':<10} {'d=r-c':<8} {'s=r+c':<8} {'dist':<6}")
for idx, ((r,c), (i,j)) in enumerate(zip(tiles, tile_arcs)):
    print(f"  ({r},{c}){'':<4} ({i},{j}){'':<6} {r-c:+d}{'':<5} {r+c:<8} {i-j}")

# ============================================================================
# BUILD ALL TILINGS AND MAP TO ISO CLASSES
# ============================================================================

print(f"\n{'='*70}")
print(f"  STEP 1: TILING -> TOURNAMENT -> ISO CLASS")
print(f"{'='*70}")

# Base path: P_0 = 4->3->2->1->0 (0-indexed)
# i.e., vertex j beats vertex j-1 for j=1..n-1
def tiling_to_tournament(tile_bits):
    a = [[0]*n for _ in range(n)]
    # Base path: vertex j -> vertex j-1 for j=1..n-1
    for j in range(1, n):
        a[j][j-1] = 1
    # Tiles
    for k in range(m):
        i, j = tile_arcs[k]
        if tile_bits & (1 << k):
            a[i][j] = 1  # upper beats lower (backward arc)
        else:
            a[j][i] = 1  # lower beats upper (forward arc)
    return a

def canon(a):
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

def Hdp(a):
    dp = {}
    for v in range(n): dp[(1<<v, v)] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S & (1<<v)): continue
            val = dp.get((S,v), 0)
            if val == 0: continue
            for u in range(n):
                if S & (1<<u): continue
                if a[v][u]:
                    key = (S|(1<<u), u)
                    dp[key] = dp.get(key, 0) + val
    return sum(dp.get(((1<<n)-1, v), 0) for v in range(n))

def comp(a):
    return [[1-a[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def ssq(a):
    return tuple(sorted(sum(a[i][j] for j in range(n)) for i in range(n)))

# Build all iso classes via tiling enumeration
tiling_to_class = {}  # tiling_bits -> class_canon
class_info = {}  # canon -> {H, score, aut, sc, tilings}

for tb in range(2**m):
    T = tiling_to_tournament(tb)
    cn = canon(T)
    tiling_to_class[tb] = cn
    if cn not in class_info:
        H = Hdp(T)
        sc = ssq(T)
        cn_comp = canon(comp(T))
        is_sc = (cn == cn_comp)
        aut = sum(1 for p in permutations(range(n))
                  if all(T[p[i]][p[j]] == T[i][j] for i in range(n) for j in range(n)))
        class_info[cn] = {'H': H, 'score': sc, 'aut': aut, 'sc': is_sc,
                          'tilings': [], 'comp_canon': cn_comp}
    class_info[cn]['tilings'].append(tb)

# Assign class IDs
classes = sorted(class_info.keys(), key=lambda cn: (class_info[cn]['H'], cn))
class_id = {cn: i for i, cn in enumerate(classes)}

print(f"\n  ISO CLASSES ({len(classes)} total):")
print(f"  {'ID':>3} {'H':>4} {'|Aut|':>5} {'SC':>3} {'#tilings':>8} {'H/|Aut|':>8} {'score'}")
for cn in classes:
    ci = class_info[cn]
    cid = class_id[cn]
    fib = ci['H'] / ci['aut']
    print(f"  {cid:3d} {ci['H']:4d} {ci['aut']:5d} {'Y' if ci['sc'] else 'N':>3} "
          f"{len(ci['tilings']):8d} {fib:8.1f} {ci['score']}")

# Verify fundamental identity: sum fiber(C) = 2^m
total_fibers = sum(len(ci['tilings']) for ci in class_info.values())
print(f"\n  sum(#tilings) = {total_fibers} = 2^{m} = {2**m} {'CHECK' if total_fibers == 2**m else 'FAIL!'}")

# Verify: #tilings = H/|Aut|?
print(f"\n  ORBIT-STABILIZER VERIFICATION: #tilings = H/|Aut| ?")
all_match = True
for cn in classes:
    ci = class_info[cn]
    predicted = ci['H'] / ci['aut']
    actual = len(ci['tilings'])
    match = abs(predicted - actual) < 0.01
    if not match: all_match = False
    if not match or len(classes) <= 15:
        print(f"    Class {class_id[cn]}: predicted={predicted:.0f}, actual={actual} "
              f"{'OK' if match else 'MISMATCH!'}")
print(f"  All match: {all_match}")

# ============================================================================
# STEP 2: TILING FLIPS -> META-GRAPH EDGES
# ============================================================================

print(f"\n{'='*70}")
print(f"  STEP 2: TILE FLIPS -> META-GRAPH EDGES")
print(f"{'='*70}")

# Build merged classes
merge_map = {}  # canon -> merged_id
merged_id = 0
for cn in classes:
    if cn in merge_map: continue
    ci = class_info[cn]
    merge_map[cn] = merged_id
    if not ci['sc']:
        merge_map[ci['comp_canon']] = merged_id
    merged_id += 1
V_merged = merged_id

merged_info = {}
for cn in classes:
    mid = merge_map[cn]
    if mid not in merged_info:
        ci = class_info[cn]
        merged_info[mid] = {'H': ci['H'], 'sc': ci['sc'], 'score': ci['score'],
                            'aut': ci['aut']}

print(f"  V_merged = {V_merged}")

# For each tiling, flip each tile, find the resulting merged class
# Track: which tile generated the edge, and its color
edge_by_tile = defaultdict(lambda: defaultdict(int))  # (tile_idx) -> {(src_mid, tgt_mid): count}
edge_colors = {}
self_loops_by_tile = defaultdict(int)

for tb in range(2**m):
    cn = tiling_to_class[tb]
    src_mid = merge_map[cn]

    for k in range(m):
        tb_flip = tb ^ (1 << k)
        cn_flip = tiling_to_class[tb_flip]
        tgt_mid = merge_map[cn_flip]

        if tgt_mid == src_mid:
            self_loops_by_tile[k] += 1
        else:
            e = (min(src_mid, tgt_mid), max(src_mid, tgt_mid))
            edge_by_tile[k][e] += 1
            if e not in edge_colors:
                edge_colors[e] = 'blue' if merged_info[src_mid]['sc'] == merged_info[tgt_mid]['sc'] else 'black'

# Total edges
all_edges = set()
for k in range(m):
    for e in edge_by_tile[k]:
        all_edges.add(e)

print(f"  Total merged edges: {len(all_edges)}")
print(f"  Blue: {sum(1 for e in all_edges if edge_colors.get(e) == 'blue')}")
print(f"  Black: {sum(1 for e in all_edges if edge_colors.get(e) == 'black')}")

# ============================================================================
# STEP 3: TILE-BY-TILE DECOMPOSITION
# ============================================================================

print(f"\n{'='*70}")
print(f"  STEP 3: TILE-BY-TILE DECOMPOSITION")
print(f"{'='*70}")

print(f"\n  For each tile (r,c) in the staircase:")
print(f"  {'tile':>8} {'(r,c)':>6} {'d':>4} {'s':>4} {'SL':>5} {'edges':>6} "
      f"{'blue':>5} {'black':>6} {'blue%':>7}")

for k in range(m):
    r, c = tiles[k]
    d = r - c
    s = r + c
    sl = self_loops_by_tile[k]
    edges_from_tile = edge_by_tile[k]
    n_edges = len(edges_from_tile)
    n_blue = sum(1 for e in edges_from_tile if edge_colors.get(e) == 'blue')
    n_black = n_edges - n_blue
    total_trans = sl + sum(edges_from_tile.values())
    blue_frac = n_blue / n_edges if n_edges > 0 else 0

    # Weighted blue fraction (by transition count)
    blue_trans = sum(c for e, c in edges_from_tile.items() if edge_colors.get(e) == 'blue')
    total_cross = sum(edges_from_tile.values())
    wt_blue_frac = blue_trans / total_cross if total_cross > 0 else 0

    print(f"  tile {k:2d}  ({r},{c}) {d:+3d} {s:4d} {sl:5d} {n_edges:6d} "
          f"{n_blue:5d} {n_black:6d} {wt_blue_frac:7.3f}")

# ============================================================================
# STEP 4: WHICH EDGES DOES EACH TILE GENERATE?
# ============================================================================

print(f"\n{'='*70}")
print(f"  STEP 4: TILE -> EDGE MAPPING (complete)")
print(f"{'='*70}")

for k in range(m):
    r, c = tiles[k]
    i, j = tile_arcs[k]
    edges = edge_by_tile[k]
    if not edges: continue

    print(f"\n  Tile ({r},{c}) = arc ({i},{j}), d={r-c:+d}, s={r+c}:")
    for e in sorted(edges.keys()):
        a, b = e
        color = edge_colors.get(e, '?')
        weight = edges[e]
        ha = merged_info[a]['H']
        hb = merged_info[b]['H']
        print(f"    -> edge ({a},{b}) [{color:5s}] weight={weight:3d}  "
              f"H: {ha}->{hb} (Delta={abs(ha-hb)})")

# ============================================================================
# STEP 5: THE DIAGONAL STRUCTURE
# ============================================================================

print(f"\n{'='*70}")
print(f"  STEP 5: DIAGONAL AND ANTI-DIAGONAL STRUCTURE")
print(f"{'='*70}")

# Group edges by the diagonal (d) of the generating tile
print(f"\n  Edges grouped by DIAGONAL d = r-c of generating tile:")
for d_val in sorted(set(r-c for r,c in tiles)):
    tiles_at_d = [k for k in range(m) if tiles[k][0] - tiles[k][1] == d_val]
    all_edges_at_d = set()
    for k in tiles_at_d:
        for e in edge_by_tile[k]:
            all_edges_at_d.add(e)
    blue_at_d = sum(1 for e in all_edges_at_d if edge_colors.get(e) == 'blue')
    black_at_d = len(all_edges_at_d) - blue_at_d
    print(f"    d={d_val:+d}: {len(tiles_at_d)} tiles, {len(all_edges_at_d)} distinct edges "
          f"({blue_at_d} blue, {black_at_d} black)")

# Group edges by anti-diagonal s = r+c
print(f"\n  Edges grouped by ANTI-DIAGONAL s = r+c of generating tile:")
for s_val in sorted(set(r+c for r,c in tiles)):
    tiles_at_s = [k for k in range(m) if tiles[k][0] + tiles[k][1] == s_val]
    all_edges_at_s = set()
    for k in tiles_at_s:
        for e in edge_by_tile[k]:
            all_edges_at_s.add(e)
    blue_at_s = sum(1 for e in all_edges_at_s if edge_colors.get(e) == 'blue')
    black_at_s = len(all_edges_at_s) - blue_at_s
    print(f"    s={s_val}: {len(tiles_at_s)} tiles, {len(all_edges_at_s)} distinct edges "
          f"({blue_at_s} blue, {black_at_s} black)")

# ============================================================================
# STEP 6: THE COMPLETE CHAIN VERIFICATION
# ============================================================================

print(f"\n{'='*70}")
print(f"  STEP 6: MASTER IDENTITIES")
print(f"{'='*70}")

total_sl = sum(self_loops_by_tile.values())
total_cross = sum(sum(edges.values()) for edges in edge_by_tile.values())
total_trans = total_sl + total_cross

print(f"\n  TILING-LEVEL IDENTITIES:")
print(f"    Total tilings: {2**m}")
print(f"    Total (tiling, tile) pairs: {2**m * m}")
print(f"    Self-loops: {total_sl}")
print(f"    Cross-transitions: {total_cross}")
print(f"    Check: {total_sl} + {total_cross} = {total_trans} = {2**m * m} {'OK' if total_trans == 2**m * m else 'FAIL'}")

# Edge thickness
thicknesses = Counter()
for e in all_edges:
    total_weight = sum(edge_by_tile[k].get(e, 0) for k in range(m))
    thicknesses[total_weight] += 1

print(f"\n  EDGE THICKNESS (tiling model):")
for t in sorted(thicknesses.keys()):
    ct = thicknesses[t]
    color_dist = Counter()
    for e in all_edges:
        w = sum(edge_by_tile[k].get(e, 0) for k in range(m))
        if w == t:
            color_dist[edge_colors.get(e, '?')] += 1
    print(f"    thickness={t}: {ct} edges ({dict(color_dist)})")

avg_thick = total_cross / len(all_edges) if all_edges else 0
print(f"\n  Avg thickness: {avg_thick:.2f}")
print(f"  E = (total_trans - SL) / avg_thick = ({total_trans} - {total_sl}) / {avg_thick:.2f} "
      f"= {(total_trans - total_sl) / avg_thick:.1f} (actual: {len(all_edges)})")

# The fiber identity
print(f"\n  FIBER IDENTITY: sum H(C)/|Aut(C)| = 2^m")
fib_sum = sum(ci['H'] / ci['aut'] for ci in class_info.values())
print(f"    sum = {fib_sum:.0f}, 2^m = {2**m} {'OK' if abs(fib_sum - 2**m) < 0.01 else 'FAIL'}")

# Per-class neutral arcs
print(f"\n  PER-CLASS NEUTRAL ARCS (in tiling model):")
class_sl = defaultdict(int)
for tb in range(2**m):
    cn = tiling_to_class[tb]
    mid = merge_map[cn]
    for k in range(m):
        tb_flip = tb ^ (1 << k)
        cn_flip = tiling_to_class[tb_flip]
        tgt_mid = merge_map[cn_flip]
        if tgt_mid == mid:
            class_sl[mid] += 1

print(f"  {'mid':>4} {'H':>4} {'SC':>3} {'fiber':>6} {'SL':>6} {'SL/fiber':>10} {'neutral/m':>10}")
for mid in sorted(merged_info.keys()):
    mi = merged_info[mid]
    fib = sum(len(class_info[cn]['tilings']) for cn in class_info
              if merge_map[cn] == mid)
    sl = class_sl.get(mid, 0)
    sl_per_tiling = sl / fib if fib > 0 else 0
    neutral_frac = sl_per_tiling / m
    print(f"  {mid:4d} {mi['H']:4d} {'Y' if mi['sc'] else 'N':>3} {fib:6d} {sl:6d} "
          f"{sl_per_tiling:10.2f} {neutral_frac:10.4f}")


print(f"\n{'='*70}")
print(f"  SUMMARY: THE COMPLETE MAP")
print(f"{'='*70}")
print(f"""
  TILING (binary string on staircase)
       |
       | each tiling encodes a tournament with base path P_0
       v
  LABELED TOURNAMENT
       |
       | S_n orbits: fiber(C) = H(C)/|Aut(C)| tilings per class
       v
  ISO CLASS (one of {len(classes)} classes)
       |
       | complement pairing: SC classes = self, NS = paired
       v
  MERGED CLASS (one of {V_merged} merged classes)
       |
       | tile flip: changes one bit of tiling -> possibly new class
       v
  META-GRAPH EDGE (blue if same SC-type, black if different)
       |
       | {len(all_edges)} edges total ({sum(1 for e in all_edges if edge_colors.get(e)=='blue')} blue, {sum(1 for e in all_edges if edge_colors.get(e)=='black')} black)
       v
  G_{n}/Z_2 = SPINE + RIBS + BULK
""")

print(f"  DONE.")
print("=" * 80)
