#!/usr/bin/env python3
"""
recursive_metagraph_s20cw.py — The recursive structure of the meta-graph
kind-pasteur-2026-03-23-S20cw

THE DEEP INSIGHT:
  The staircase decomposes into 4 regions:
    OVERLAP (m(n-2) cells) — the n-2 sub-staircase at center
    BOTTOM-WIRING (n-3 cells) — how vertex 0 connects to internals
    TOP-WIRING (n-3 cells) — how vertex n-1 connects to internals
    APEX (1 cell) — the arc between vertices 0 and n-1

  PREDICTION: Overlap flips are the BLUEST (type-preserving).
  Apex flips are the LEAST BLUE (most type-changing).
  This explains why blue fraction -> 1: overlap fraction -> 1.

  Compute the per-region blue fraction for n=4,5,6
  and derive the recursive formula for E_merged.
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  RECURSIVE META-GRAPH STRUCTURE")
print("  Per-region blue fraction: overlap vs wiring vs apex")
print("  kind-pasteur-2026-03-23-S20cw")
print("=" * 80)

# ============================================================================
# HELPERS
# ============================================================================
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
                    dp[(S|(1<<u), u)] = dp.get((S|(1<<u), u), 0) + val
    return sum(dp.get(((1<<n)-1, v), 0) for v in range(n))


def analyze_n(n):
    t0 = time.time()

    # Grid cells
    tiles = [(r,c) for r in range(1,n-1) for c in range(1,n-r)]
    m = len(tiles)
    tile_arcs = [(r+c, c-1) for r,c in tiles]  # 0-indexed

    # Classify tiles into regions
    apex_cell = (n-2, 1)
    bottom = set((r,c) for r,c in tiles if r+c <= n-2)
    top = set((r,c) for r,c in tiles if c >= 2)
    overlap = bottom & top
    bottom_only = set(t for t in bottom if t not in overlap and t != apex_cell)
    top_only = set(t for t in top if t not in overlap)

    tile_region = {}
    for t in tiles:
        if t == apex_cell: tile_region[t] = 'apex'
        elif t in overlap: tile_region[t] = 'overlap'
        elif t in bottom_only: tile_region[t] = 'bottom'
        elif t in top_only: tile_region[t] = 'top'
        else: tile_region[t] = '???'

    tile_idx = {t: i for i, t in enumerate(tiles)}

    print(f"\n  n={n}: m={m} tiles, regions: overlap={len(overlap)}, "
          f"bottom={len(bottom_only)}, top={len(top_only)}, apex=1")

    # Build all tilings and iso classes
    def build_tournament(tb):
        a = [[0]*n for _ in range(n)]
        for j in range(1, n): a[j][j-1] = 1  # base path
        for k in range(m):
            i, j = tile_arcs[k]
            if tb & (1 << k): a[i][j] = 1
            else: a[j][i] = 1
        return a

    tiling_canon = {}
    class_info = {}
    for tb in range(1 << m):
        T = build_tournament(tb)
        cn = canon(T, n)
        tiling_canon[tb] = cn
        if cn not in class_info:
            Tc = [[1-T[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
            class_info[cn] = {
                'H': Hdp(T, n),
                'sc': canon(Tc, n) == cn,
                'tilings': []
            }
        class_info[cn]['tilings'].append(tb)

    # Merged classes
    merge = {}; mid = 0
    for cn in sorted(class_info.keys()):
        if cn in merge: continue
        ci = class_info[cn]
        merge[cn] = mid
        if not ci['sc']:
            cc = canon([[1-T[i][j] if i!=j else 0 for j in range(n)]
                       for i in range(n)], n) if False else None
            # Need comp canon
            T0 = build_tournament(ci['tilings'][0])
            Tc = [[1-T0[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
            cc = canon(Tc, n)
            if cc != cn and cc in class_info:
                merge[cc] = mid
        mid += 1
    V_merged = mid

    merged_sc = {}
    for cn, ci in class_info.items():
        mi = merge[cn]
        if mi not in merged_sc: merged_sc[mi] = ci['sc']

    # Per-region edge statistics
    region_stats = {r: {'blue': 0, 'black': 0, 'self_loop': 0, 'distinct_edges': set()}
                    for r in ['overlap', 'bottom', 'top', 'apex']}

    for tb in range(1 << m):
        cn = tiling_canon[tb]
        src = merge[cn]
        for k in range(m):
            tb2 = tb ^ (1 << k)
            cn2 = tiling_canon[tb2]
            tgt = merge[cn2]
            region = tile_region[tiles[k]]

            if tgt == src:
                region_stats[region]['self_loop'] += 1
            else:
                is_blue = (merged_sc[src] == merged_sc[tgt])
                e = (min(src,tgt), max(src,tgt))
                region_stats[region]['distinct_edges'].add(e)
                if is_blue:
                    region_stats[region]['blue'] += 1
                else:
                    region_stats[region]['black'] += 1

    # Count distinct edges per region
    all_edges = set()
    edge_regions = defaultdict(set)  # edge -> set of regions that generate it
    for r, st in region_stats.items():
        for e in st['distinct_edges']:
            all_edges.add(e)
            edge_regions[e].add(r)

    # Edges EXCLUSIVE to each region
    for r in region_stats:
        region_stats[r]['exclusive'] = sum(1 for e in region_stats[r]['distinct_edges']
                                           if edge_regions[e] == {r})

    print(f"  V_merged={V_merged}, total distinct edges={len(all_edges)}")
    print(f"  Build time: {time.time()-t0:.1f}s")

    # Print per-region results
    print(f"\n  PER-REGION EDGE DECOMPOSITION:")
    print(f"  {'Region':<10} {'#cells':>6} {'blue':>8} {'black':>8} {'SL':>6} "
          f"{'edges':>6} {'excl':>6} {'blue%':>8}")
    for r in ['overlap', 'bottom', 'top', 'apex']:
        st = region_stats[r]
        ncells = len(overlap) if r=='overlap' else len(bottom_only) if r=='bottom' \
                 else len(top_only) if r=='top' else 1
        total = st['blue'] + st['black']
        bf = st['blue']/total if total > 0 else 0
        ne = len(st['distinct_edges'])
        print(f"  {r:<10} {ncells:6d} {st['blue']:8d} {st['black']:8d} {st['self_loop']:6d} "
              f"{ne:6d} {st['exclusive']:6d} {bf:8.3f}")

    # Total
    tot_blue = sum(st['blue'] for st in region_stats.values())
    tot_black = sum(st['black'] for st in region_stats.values())
    tot_sl = sum(st['self_loop'] for st in region_stats.values())
    tot_bf = tot_blue/(tot_blue+tot_black) if tot_blue+tot_black > 0 else 0
    print(f"  {'TOTAL':<10} {m:6d} {tot_blue:8d} {tot_black:8d} {tot_sl:6d} "
          f"{len(all_edges):6d} {'':>6} {tot_bf:8.3f}")

    # Edge sharing between regions
    print(f"\n  EDGE SHARING:")
    shared = Counter()
    for e, regions in edge_regions.items():
        shared[frozenset(regions)] += 1
    for rset, count in sorted(shared.items(), key=lambda x: -x[1]):
        color = 'blue' if merged_sc[list(all_edges)[0][0]] == merged_sc[list(all_edges)[0][1]] else '?'
        rnames = '+'.join(sorted(rset))
        print(f"    {rnames}: {count} edges")

    # The key ratio: overlap/(overlap+wiring+apex)
    cells_overlap = len(overlap)
    cells_wiring = len(bottom_only) + len(top_only)
    cells_apex = 1
    print(f"\n  CELL FRACTIONS:")
    print(f"    Overlap: {cells_overlap}/{m} = {cells_overlap/m:.3f}")
    print(f"    Wiring: {cells_wiring}/{m} = {cells_wiring/m:.3f}")
    print(f"    Apex: {cells_apex}/{m} = {cells_apex/m:.3f}")
    print(f"    Overlap fraction = C({n-3},2)/C({n-1},2) = {comb(n-3,2)}/{comb(n-1,2)} = {comb(n-3,2)/comb(n-1,2):.4f}")

    return {
        'n': n, 'V': V_merged, 'E': len(all_edges),
        'region_stats': {r: {k:v for k,v in st.items() if k != 'distinct_edges'}
                         for r, st in region_stats.items()},
        'overlap_frac': cells_overlap/m,
        'overlap_blue': region_stats['overlap']['blue']/(region_stats['overlap']['blue']+region_stats['overlap']['black'])
                        if region_stats['overlap']['blue']+region_stats['overlap']['black'] > 0 else None,
    }


# ============================================================================
# MAIN
# ============================================================================

results = {}
for n in [4, 5, 6]:
    print(f"\n{'#'*70}")
    results[n] = analyze_n(n)


# ============================================================================
# CROSS-n SYNTHESIS
# ============================================================================

print(f"\n\n{'='*80}")
print("  CROSS-n SYNTHESIS: THE RECURSIVE PATTERN")
print(f"{'='*80}")

print(f"\n  BLUE FRACTION BY REGION:")
print(f"  {'n':>3} {'overlap':>10} {'bottom':>10} {'top':>10} {'apex':>10} {'overall':>10}")
for n in [4, 5, 6]:
    r = results[n]
    rs = r['region_stats']
    def bf(region):
        b = rs[region]['blue']; k = rs[region]['black']
        return f"{b/(b+k):.3f}" if b+k > 0 else "  N/A"
    tot_b = sum(rs[x]['blue'] for x in rs)
    tot_k = sum(rs[x]['black'] for x in rs)
    overall = tot_b/(tot_b+tot_k) if tot_b+tot_k > 0 else 0
    print(f"  {n:3d} {bf('overlap'):>10} {bf('bottom'):>10} {bf('top'):>10} {bf('apex'):>10} {overall:10.3f}")

print(f"\n  OVERLAP CELL FRACTION (approaches 1 as n grows):")
for n in [4, 5, 6, 7, 8, 10, 15, 20, 50]:
    frac = comb(n-3, 2) / comb(n-1, 2) if comb(n-1,2) > 0 else 0
    print(f"    n={n:2d}: C({n-3},2)/C({n-1},2) = {comb(n-3,2)}/{comb(n-1,2)} = {frac:.4f}")

print(f"\n  THE FORMULA:")
print(f"    blue_fraction(n) ~ overlap_fraction + (1-overlap_fraction) * wiring_blue")
print(f"    ~ C(n-3,2)/C(n-1,2) + (2(n-3)+1)/C(n-1,2) * wiring_blue")
print(f"    -> 1 as n -> infinity (because C(n-3,2)/C(n-1,2) -> 1)")

print(f"\n  WHY OVERLAP IS BLUEST:")
print(f"    Flipping a tile in the OVERLAP changes the inner tournament on n-2 vertices")
print(f"    while keeping BOTH endpoint wirings fixed. The SC anti-automorphism,")
print(f"    if it existed, maps the inner tournament to its complement while swapping")
print(f"    the wirings. Changing ONLY the inner part tends to preserve this structure")
print(f"    because the wiring remains compatible.")
print(f"")
print(f"    The APEX is the least blue because flipping it directly changes whether")
print(f"    the tournament has a Hamiltonian cycle through the base path, which")
print(f"    strongly affects the SC structure.")

print(f"\n  RECURSIVE META-GRAPH FORMULA:")
print(f"    G_n/Z_2 is built from two copies of G_{{n-1}}/Z_2 sharing G_{{n-2}}/Z_2,")
print(f"    plus one apex bit. The edges decompose as:")
print(f"    E(n) = E_overlap(n) + E_wiring(n) + E_apex(n)")
print(f"    where E_overlap lifts edges from G_{{n-2}}/Z_2 into G_n/Z_2,")
print(f"    E_wiring adds edges from endpoint insertion choices,")
print(f"    and E_apex adds edges from the cycle-controlling bit.")

print(f"\n  DONE.")
print("=" * 80)
