#!/usr/bin/env python3
"""
unified_weights_s20et.py — Unified view: blue/black counts + wiggly weights
kind-pasteur-2026-03-24-S20et

THE DEEP SIMILARITY:
  Both blue/black (node decoration) and wiggly (edge decoration) have
  PERFECT FACTORIZATION:
    blue(node) = tilings(node) × 2^{-(n-2)}   [verified]
    wiggly(edge, class X) = base(edge) × factor(X)  [opus S277]

  This script computes BOTH decorations simultaneously and checks:
  1. Does the wiggly factorization hold in the TILING model (tiles only, not arcs)?
  2. What is base(edge) in terms of node properties?
  3. Do blue/black and wiggly interact? (wiggly lines from blue tilings vs black)
"""

import sys
from math import factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  UNIFIED BLUE/BLACK + WIGGLY WEIGHT ANALYSIS")
print("  kind-pasteur-2026-03-24-S20et")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in range(4, 7):
    t0 = time.time()
    TILES = get_tiles(n)
    num_tiles = len(TILES)
    N = n
    VERTS = list(range(n, 0, -1))
    all_perms = list(permutations(range(N)))

    # Tile skip (= range - 1 in opus's terms, or x - y - 1)
    tile_skip = [x - y - 1 for x, y in TILES]
    tile_range = [x - y for x, y in TILES]  # opus's "range" = skip + 1

    def bits_to_adj(bits):
        A = [[0]*N for _ in range(N)]
        for k in range(N-1): A[k][k+1] = 1
        for i in range(num_tiles):
            x, y = TILES[i]
            xi, yi = VERTS.index(x), VERTS.index(y)
            if bits[i] == 0: A[xi][yi] = 1
            else: A[yi][xi] = 1
        return A

    def canonicalize(A):
        best = None
        for p in all_perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
            if best is None or s < best: best = s
        return best

    def is_grid_sym(bits):
        tile_idx = {(t[0],t[1]): i for i, t in enumerate(TILES)}
        for i, (x, y) in enumerate(TILES):
            tx, ty = N - y + 1, N - x + 1
            ti = tile_idx.get((tx, ty), -1)
            if ti >= 0 and ti != i and bits[i] != bits[ti]:
                return False
        return True

    # Build all tilings
    canon_map = {}
    grid_sym_map = {}
    for mask in range(1 << num_tiles):
        bits = [(mask >> k) & 1 for k in range(num_tiles)]
        A = bits_to_adj(bits)
        canon_map[mask] = canonicalize(A)
        grid_sym_map[mask] = is_grid_sym(bits)

    # Merged classes
    comp_map_cn = {}
    for cn in set(canon_map.values()):
        for mask, c in canon_map.items():
            if c == cn:
                comp_mask = mask ^ ((1 << num_tiles) - 1)
                comp_map_cn[cn] = canon_map[comp_mask]
                break

    merged_map = {cn: min(cn, comp_map_cn.get(cn, cn)) for cn in set(canon_map.values())}
    merged_classes = sorted(set(merged_map.values()))
    V_merged = len(merged_classes)

    # ================================================================
    # NODE DECORATION: blue/black counts per merged class
    # ================================================================
    node_tilings = defaultdict(int)
    node_blue = defaultdict(int)
    node_black = defaultdict(int)

    for mask in range(1 << num_tiles):
        mcn = merged_map[canon_map[mask]]
        node_tilings[mcn] += 1
        if grid_sym_map[mask]:
            node_blue[mcn] += 1
        else:
            node_black[mcn] += 1

    # ================================================================
    # EDGE DECORATION: wiggly weights per wiggly class
    # ================================================================
    # edge_weights[(mcn1, mcn2)][wi] = count of wiggly lines of class wi
    edge_weights = defaultdict(lambda: defaultdict(int))
    # self_weights[mcn][wi] = count of wiggly self-loop lines of class wi
    self_weights = defaultdict(lambda: defaultdict(int))

    # Also decompose by blue/black SOURCE
    edge_from_blue = defaultdict(lambda: defaultdict(int))
    edge_from_black = defaultdict(lambda: defaultdict(int))
    self_from_blue = defaultdict(lambda: defaultdict(int))
    self_from_black = defaultdict(lambda: defaultdict(int))

    for mask in range(1 << num_tiles):
        mcn1 = merged_map[canon_map[mask]]
        is_blue = grid_sym_map[mask]

        for wi in range(num_tiles):
            flip_mask = mask ^ (1 << wi)
            mcn2 = merged_map[canon_map[flip_mask]]

            if mcn1 == mcn2:
                self_weights[mcn1][wi] += 1
                if is_blue: self_from_blue[mcn1][wi] += 1
                else: self_from_black[mcn1][wi] += 1
            else:
                edge = (min(mcn1, mcn2), max(mcn1, mcn2))
                edge_weights[edge][wi] += 1
                if is_blue: edge_from_blue[edge][wi] += 1
                else: edge_from_black[edge][wi] += 1

    all_edges = sorted(edge_weights.keys())
    E = len(all_edges)

    print(f"\n{'#'*70}")
    print(f"  n = {n}, V_merged = {V_merged}, E = {E}, tiles = {num_tiles}")
    print(f"{'#'*70}")

    # ================================================================
    # CHECK 1: Blue fraction per node = 2^{-(n-2)} ?
    # ================================================================
    expected_blue_frac = 2**(-(n-2))
    print(f"\n  NODE BLUE FRACTION (expected = 2^{{-{n-2}}} = {expected_blue_frac:.6f}):")
    blue_fracs = []
    for mcn in merged_classes:
        total = node_tilings[mcn]
        blue = node_blue[mcn]
        frac = blue / total if total > 0 else 0
        blue_fracs.append(frac)
    print(f"    Actual: min={min(blue_fracs):.6f}, max={max(blue_fracs):.6f}")
    print(f"    All equal to expected? {all(abs(f - expected_blue_frac) < 1e-10 for f in blue_fracs)}")

    # ================================================================
    # CHECK 2: Wiggly factorization by SKIP
    # For tiling model: weight(edge, class X) = base(edge) × f(skip(X))
    # Group classes by skip and check
    # ================================================================
    skip_values = sorted(set(tile_skip))

    print(f"\n  WIGGLY FACTORIZATION BY SKIP:")
    # For each edge, compute weight per skip group
    factorization_holds = True
    for edge in all_edges[:5]:  # show first 5
        weights_by_skip = {}
        for s in skip_values:
            w = sum(edge_weights[edge][wi] for wi in range(num_tiles) if tile_skip[wi] == s)
            count = sum(1 for wi in range(num_tiles) if tile_skip[wi] == s)
            weights_by_skip[s] = (w, count)

        # Check: w/count is proportional to the same factor across all edges?
        per_tile = {s: w/c if c > 0 else 0 for s, (w,c) in weights_by_skip.items()}
        print(f"    Edge {edge}: per-tile weight by skip: {dict(per_tile)}")

    # For ALL edges: compute per-tile weight by skip, check if ratio is constant
    skip_ratios = defaultdict(list)  # skip -> list of (weight_per_tile / baseline)
    for edge in all_edges:
        weights = {}
        for s in skip_values:
            tiles_in_skip = [wi for wi in range(num_tiles) if tile_skip[wi] == s]
            if tiles_in_skip:
                w = sum(edge_weights[edge][wi] for wi in tiles_in_skip)
                weights[s] = w / len(tiles_in_skip)  # per-tile weight

        # Normalize by skip=min (baseline)
        baseline_skip = min(skip_values)
        baseline = weights.get(baseline_skip, 1)
        if baseline > 0:
            for s in skip_values:
                if s in weights:
                    skip_ratios[s].append(weights[s] / baseline)

    print(f"\n  SKIP RATIO (normalized to skip={min(skip_values)}):")
    for s in skip_values:
        ratios = skip_ratios[s]
        if ratios:
            print(f"    skip={s}: min={min(ratios):.4f}, max={max(ratios):.4f}, constant={max(ratios)-min(ratios)<0.001}")

    # ================================================================
    # CHECK 3: Wiggly lines from BLUE vs BLACK tilings
    # ================================================================
    print(f"\n  WIGGLY FROM BLUE vs BLACK TILINGS:")
    for edge in all_edges[:5]:
        total = sum(edge_weights[edge].values())
        from_blue = sum(edge_from_blue[edge].values())
        from_black = sum(edge_from_black[edge].values())
        blue_frac = from_blue / total if total > 0 else 0
        print(f"    Edge {edge}: total={total}, from_blue={from_blue} ({blue_frac:.4f}), from_black={from_black}")

    # Is the blue fraction of wiggly lines = node blue fraction?
    all_blue_fracs = []
    for edge in all_edges:
        total = sum(edge_weights[edge].values())
        from_blue = sum(edge_from_blue[edge].values())
        if total > 0:
            all_blue_fracs.append(from_blue / total)

    print(f"\n    Blue fraction of wiggly lines across all edges:")
    print(f"    min={min(all_blue_fracs):.4f}, max={max(all_blue_fracs):.4f}, expected={expected_blue_frac:.4f}")
    print(f"    All equal to node blue fraction? {all(abs(f - expected_blue_frac) < 0.001 for f in all_blue_fracs)}")

    # ================================================================
    # CHECK 4: Base weight formula
    # ================================================================
    print(f"\n  BASE WEIGHT ANALYSIS:")
    # Base weight = total wiggly weight / sum(tiles per skip × factor)
    # If factorization holds: base = total_weight / sum_skip(count(skip) × factor(skip))
    for edge in all_edges[:8]:
        total = sum(edge_weights[edge].values())
        mcn1, mcn2 = edge
        t1 = node_tilings[mcn1]
        t2 = node_tilings[mcn2]
        ratio = total / (t1 * t2) if t1*t2 > 0 else 0
        print(f"    {edge}: total_weight={total}, tilings={t1}*{t2}={t1*t2}, weight/(t1*t2)={ratio:.6f}")

    # Check: weight / (tilings product) constant across edges?
    weight_per_tp = []
    for edge in all_edges:
        total = sum(edge_weights[edge].values())
        t1, t2 = node_tilings[edge[0]], node_tilings[edge[1]]
        if t1*t2 > 0:
            weight_per_tp.append(total / (t1*t2))

    if weight_per_tp:
        print(f"\n    weight/(t1*t2): min={min(weight_per_tp):.6f}, max={max(weight_per_tp):.6f}")
        print(f"    Constant? {max(weight_per_tp) - min(weight_per_tp) < 0.0001}")

    # ================================================================
    # SUMMARY TABLE: unified node + edge decoration
    # ================================================================
    print(f"\n  UNIFIED VIEW:")
    print(f"    NODES: tilings × blue_fraction = blue_count")
    print(f"    EDGES: base_weight × skip_factor = class_weight")
    print(f"    Blue fraction of wiggly lines = blue fraction of nodes (same 2^{{-{n-2}}})")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
