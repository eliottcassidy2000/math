#!/usr/bin/env python3
"""
edge_labelings_s20cz.py — All meaningful edge labelings of G_n/Z_2
kind-pasteur-2026-03-23-S20cz

Compute EVERY possible edge classification and find which ones
reveal structure, correlate with each other, or yield formulas.

LABELINGS:
  1. BLUE/BLACK (explorer): grid-symmetry of endpoint tilings
  2. SC-SC/SC-NS/NS-NS (spine/ribs/sea): class-level SC-type
  3. H-direction: ASCENDING / LEVEL (DAG, no descending)
  4. Delta-H magnitude: STEP (|dH|=2) / LEAP (|dH|>2) / LEVEL
  5. Score change: AMBER (same score seq) / VIOLET (different)
  6. c3 change: dc3 value
  7. Staircase region: OVERLAP / WIRING / APEX
  8. Fiber (Mode A): parent class at n-1
  9. |Aut| change

For each pair of labelings, compute the contingency table.
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  ALL EDGE LABELINGS OF G_n/Z_2")
print("  kind-pasteur-2026-03-23-S20cz")
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

# ============================================================================
# BUILD G_n/Z_2 WITH ALL LABELS
# ============================================================================

for n in [5, 6]:
    print(f"\n{'#'*80}")
    print(f"  n = {n}")
    print(f"{'#'*80}")
    t0 = time.time()

    # Tiles
    tiles = [(r,c) for r in range(1,n-1) for c in range(1,n-r)]
    m = len(tiles)
    tile_arcs = [(r+c, c-1) for r,c in tiles]

    # TRANS_MAP (explorer's grid reflection): (x,y) -> (n-y+1, n-x+1) in explorer coords
    # In our (r,c): tile (r,c) = arc from vertex r+c to vertex c-1 (0-indexed)
    # Explorer: tile (x,y) where x=r+c+1, y=c (1-indexed)
    # TRANS_MAP: (x,y) -> (n-y+1, n-x+1) = (n-c+1, n-r-c) in 1-indexed
    # In (r,c): reflected r' = n - (r+c+1) - (n-c) + 1... this is getting complicated.
    # Let me just build it the explorer's way.
    explorer_tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            explorer_tiles.append((x, y))
    tile_key = {(x,y): i for i, (x,y) in enumerate(explorer_tiles)}
    trans_map = [tile_key.get((n-y+1, n-x+1), -1) for x,y in explorer_tiles]

    def is_grid_sym(bits):
        for i in range(len(explorer_tiles)):
            ti = trans_map[i]
            if ti != i and ti >= 0 and ((bits >> i) & 1) != ((bits >> ti) & 1):
                return False
        return True

    # Staircase regions
    apex_cell_rc = (n-2, 1)
    overlap_set = set((r,c) for r,c in tiles if c >= 2 and r+c <= n-2)
    bottom_only_set = set((r,c) for r,c in tiles if c < 2 and (r,c) != apex_cell_rc and r+c <= n-2)
    top_only_set = set((r,c) for r,c in tiles if c >= 2 and r+c > n-2)
    # Handle edge cases for small n
    bottom_extra = set((r,c) for r,c in tiles if c == 1 and r+c <= n-2 and (r,c) != apex_cell_rc)

    def tile_region(r, c):
        if (r,c) == apex_cell_rc: return 'apex'
        if (r,c) in overlap_set: return 'overlap'
        if c == 1 and r+c <= n-2: return 'bottom'
        if c >= 2 and r+c > n-2: return 'top'
        if c == 1 and r+c > n-2 and (r,c) != apex_cell_rc: return 'top'  # shouldn't happen
        return 'bottom'

    # Build tournaments from tilings
    def build_tournament(tb):
        a = [[0]*n for _ in range(n)]
        for j in range(1, n): a[j][j-1] = 1
        for k in range(m):
            i, j = tile_arcs[k]
            if tb & (1 << k): a[i][j] = 1
            else: a[j][i] = 1
        return a

    # Enumerate all tilings
    tiling_data = {}
    class_info = {}
    for tb in range(1 << m):
        T = build_tournament(tb)
        cn = canon(T, n)
        gs = is_grid_sym(tb)
        sc_seq = tuple(sorted(sum(T[i][j] for j in range(n)) for i in range(n)))

        if cn not in class_info:
            Tc = [[1-T[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
            c3 = 0
            for i in range(n):
                for j in range(i+1,n):
                    for k in range(j+1,n):
                        if T[i][j] and T[j][k] and T[k][i]: c3 += 1
                        if T[i][k] and T[k][j] and T[j][i]: c3 += 1
            # Parent class (delete vertex n-1)
            T_del = [[T[i][j] for j in range(n-1)] for i in range(n-1)]
            parent_cn = canon(T_del, n-1)
            aut = sum(1 for p in permutations(range(n))
                      if all(T[p[i]][p[j]] == T[i][j] for i in range(n) for j in range(n)))
            class_info[cn] = {
                'H': Hdp(T, n), 'sc': canon(Tc, n) == cn,
                'score': sc_seq, 'c3': c3, 'aut': aut,
                'parent': parent_cn, 'tilings': []
            }
        class_info[cn]['tilings'].append(tb)
        tiling_data[tb] = {'canon': cn, 'grid_sym': gs}

    # Merge map
    merge = {}; mid = 0
    for cn in sorted(class_info.keys()):
        if cn in merge: continue
        ci = class_info[cn]
        merge[cn] = mid
        if not ci['sc']:
            T0 = build_tournament(ci['tilings'][0])
            Tc = [[1-T0[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
            cc = canon(Tc, n)
            if cc != cn and cc in class_info:
                merge[cc] = mid
        mid += 1
    V_merged = mid

    minfo = {}
    for cn, ci in class_info.items():
        mi = merge[cn]
        if mi not in minfo:
            minfo[mi] = {'H': ci['H'], 'sc': ci['sc'], 'score': ci['score'],
                         'c3': ci['c3'], 'aut': ci['aut'], 'parent': ci['parent']}

    print(f"  V_merged = {V_merged}, {len(class_info)} iso classes, built in {time.time()-t0:.1f}s")

    # ================================================================
    # COMPUTE ALL EDGE LABELS
    # ================================================================

    edges = {}  # (a,b) -> dict of labels

    for tb in range(1 << m):
        td = tiling_data[tb]
        src_cn = td['canon']
        src_mid = merge[src_cn]

        for k in range(m):
            tb2 = tb ^ (1 << k)
            td2 = tiling_data[tb2]
            tgt_cn = td2['canon']
            tgt_mid = merge[tgt_cn]

            if tgt_mid == src_mid: continue

            e = (min(src_mid, tgt_mid), max(src_mid, tgt_mid))
            if e not in edges:
                a, b = e
                ma, mb = minfo[a], minfo[b]
                dH = abs(ma['H'] - mb['H'])

                # SC-type classification
                if ma['sc'] and mb['sc']: sc_type = 'SC-SC'
                elif not ma['sc'] and not mb['sc']: sc_type = 'NS-NS'
                else: sc_type = 'SC-NS'

                edges[e] = {
                    'sc_type': sc_type,
                    'dH': dH,
                    'H_dir': 'level' if dH == 0 else 'ascending',
                    'dH_class': 'level' if dH == 0 else ('step' if dH == 2 else 'leap'),
                    'score_change': 'same' if ma['score'] == mb['score'] else 'diff',
                    'dc3': abs(ma['c3'] - mb['c3']),
                    'daut': abs(ma['aut'] - mb['aut']),
                    'same_parent': ma['parent'] == mb['parent'],
                    'blue_count': 0, 'black_count': 0,
                    'regions': Counter(),
                    'total_weight': 0,
                }

            # Blue/Black (explorer definition)
            both_sym = td['grid_sym'] and td2['grid_sym']
            if both_sym:
                edges[e]['blue_count'] += 1
            else:
                edges[e]['black_count'] += 1

            # Region
            r, c = tiles[k]
            edges[e]['regions'][tile_region(r, c)] += 1
            edges[e]['total_weight'] += 1

    # Classify blue/black at edge level
    for e, info in edges.items():
        if info['blue_count'] > 0 and info['black_count'] == 0:
            info['color'] = 'BLUE'
        elif info['blue_count'] == 0:
            info['color'] = 'BLACK'
        else:
            info['color'] = 'MIXED'

    E = len(edges)
    print(f"  Edges: {E}")

    # ================================================================
    # CROSS-TABULATIONS
    # ================================================================

    print(f"\n  --- CROSS-TABULATION: Explorer color vs SC-type ---")
    ct = Counter()
    for e, info in edges.items():
        ct[(info['color'], info['sc_type'])] += 1
    print(f"  {'':>8} {'SC-SC':>8} {'SC-NS':>8} {'NS-NS':>8} {'total':>8}")
    for color in ['BLUE', 'BLACK', 'MIXED']:
        row = [ct.get((color, t), 0) for t in ['SC-SC', 'SC-NS', 'NS-NS']]
        if sum(row) > 0:
            print(f"  {color:>8} {row[0]:8d} {row[1]:8d} {row[2]:8d} {sum(row):8d}")
    totals = [sum(ct.get((c, t), 0) for c in ['BLUE','BLACK','MIXED']) for t in ['SC-SC','SC-NS','NS-NS']]
    print(f"  {'total':>8} {totals[0]:8d} {totals[1]:8d} {totals[2]:8d} {sum(totals):8d}")

    print(f"\n  --- CROSS-TABULATION: Explorer color vs dH class ---")
    ct2 = Counter()
    for e, info in edges.items():
        ct2[(info['color'], info['dH_class'])] += 1
    print(f"  {'':>8} {'step':>8} {'leap':>8} {'level':>8}")
    for color in ['BLUE', 'BLACK', 'MIXED']:
        row = [ct2.get((color, t), 0) for t in ['step', 'leap', 'level']]
        if sum(row) > 0:
            print(f"  {color:>8} {row[0]:8d} {row[1]:8d} {row[2]:8d}")

    print(f"\n  --- CROSS-TABULATION: Explorer color vs score change ---")
    ct3 = Counter()
    for e, info in edges.items():
        ct3[(info['color'], info['score_change'])] += 1
    print(f"  {'':>8} {'same':>8} {'diff':>8}")
    for color in ['BLUE', 'BLACK', 'MIXED']:
        row = [ct3.get((color, t), 0) for t in ['same', 'diff']]
        if sum(row) > 0:
            print(f"  {color:>8} {row[0]:8d} {row[1]:8d}")

    print(f"\n  --- CROSS-TABULATION: SC-type vs dH class ---")
    ct4 = Counter()
    for e, info in edges.items():
        ct4[(info['sc_type'], info['dH_class'])] += 1
    print(f"  {'':>8} {'step':>8} {'leap':>8} {'level':>8}")
    for sc in ['SC-SC', 'SC-NS', 'NS-NS']:
        row = [ct4.get((sc, t), 0) for t in ['step', 'leap', 'level']]
        if sum(row) > 0:
            print(f"  {sc:>8} {row[0]:8d} {row[1]:8d} {row[2]:8d}")

    print(f"\n  --- CROSS-TABULATION: SC-type vs score change ---")
    ct5 = Counter()
    for e, info in edges.items():
        ct5[(info['sc_type'], info['score_change'])] += 1
    print(f"  {'':>8} {'same':>8} {'diff':>8}")
    for sc in ['SC-SC', 'SC-NS', 'NS-NS']:
        row = [ct5.get((sc, t), 0) for t in ['same', 'diff']]
        print(f"  {sc:>8} {row[0]:8d} {row[1]:8d}")

    print(f"\n  --- CROSS-TABULATION: SC-type vs fiber (same parent at n-1) ---")
    ct6 = Counter()
    for e, info in edges.items():
        ct6[(info['sc_type'], info['same_parent'])] += 1
    print(f"  {'':>8} {'same_p':>8} {'diff_p':>8}")
    for sc in ['SC-SC', 'SC-NS', 'NS-NS']:
        row = [ct6.get((sc, True), 0), ct6.get((sc, False), 0)]
        print(f"  {sc:>8} {row[0]:8d} {row[1]:8d}")

    print(f"\n  --- CROSS-TABULATION: Explorer color vs region (dominant) ---")
    ct7 = Counter()
    for e, info in edges.items():
        dom_region = info['regions'].most_common(1)[0][0] if info['regions'] else '?'
        ct7[(info['color'], dom_region)] += 1
    regions_seen = sorted(set(r for _, r in ct7.keys()))
    header = ''.join(f'{r:>8}' for r in regions_seen)
    print(f"  {'':>8}{header}")
    for color in ['BLUE', 'BLACK', 'MIXED']:
        row = [ct7.get((color, r), 0) for r in regions_seen]
        if sum(row) > 0:
            print(f"  {color:>8}" + ''.join(f'{v:8d}' for v in row))

    # ================================================================
    # dH DISTRIBUTION
    # ================================================================
    print(f"\n  --- Delta-H DISTRIBUTION ---")
    dh_dist = Counter(info['dH'] for info in edges.values())
    print(f"  {'dH':>4} {'count':>6} {'%':>7}")
    for dh in sorted(dh_dist.keys()):
        print(f"  {dh:4d} {dh_dist[dh]:6d} {100*dh_dist[dh]/E:7.1f}%")

    # ================================================================
    # SCORE-PRESERVING EDGES
    # ================================================================
    same_score = sum(1 for info in edges.values() if info['score_change'] == 'same')
    diff_score = E - same_score
    print(f"\n  SCORE-PRESERVING (AMBER) edges: {same_score}/{E} = {100*same_score/E:.1f}%")
    print(f"  SCORE-CHANGING (VIOLET) edges: {diff_score}/{E} = {100*diff_score/E:.1f}%")

    # Which score pairs are connected?
    if same_score > 0 and n <= 5:
        print(f"  Same-score edge examples:")
        for e, info in sorted(edges.items()):
            if info['score_change'] == 'same':
                a, b = e
                print(f"    ({a},{b}): H={minfo[a]['H']},{minfo[b]['H']} "
                      f"score={minfo[a]['score']} c3={minfo[a]['c3']},{minfo[b]['c3']}")

    # ================================================================
    # FIBER STRUCTURE (same parent at n-1)
    # ================================================================
    same_fiber = sum(1 for info in edges.values() if info['same_parent'])
    diff_fiber = E - same_fiber
    print(f"\n  INTRA-FIBER (TEAL, same parent) edges: {same_fiber}/{E} = {100*same_fiber/E:.1f}%")
    print(f"  INTER-FIBER (MAGENTA, diff parent) edges: {diff_fiber}/{E} = {100*diff_fiber/E:.1f}%")

    # ================================================================
    # SUMMARY TABLE
    # ================================================================
    print(f"\n  EDGE LABELING SUMMARY:")
    print(f"  {'Labeling':<25} {'Categories':>30} {'Counts':>20}")
    print(f"  {'-'*75}")

    color_counts = Counter(info['color'] for info in edges.values())
    print(f"  {'Explorer (blue/black)':<25} {'BLUE/BLACK/MIXED':>30} "
          f"{'/'.join(str(color_counts.get(c,0)) for c in ['BLUE','BLACK','MIXED']):>20}")

    sc_counts = Counter(info['sc_type'] for info in edges.values())
    print(f"  {'SC-type (spine/ribs/sea)':<25} {'SC-SC/SC-NS/NS-NS':>30} "
          f"{'/'.join(str(sc_counts.get(c,0)) for c in ['SC-SC','SC-NS','NS-NS']):>20}")

    dh_counts = Counter(info['dH_class'] for info in edges.values())
    print(f"  {'dH (step/leap/level)':<25} {'step/leap/level':>30} "
          f"{'/'.join(str(dh_counts.get(c,0)) for c in ['step','leap','level']):>20}")

    score_counts = Counter(info['score_change'] for info in edges.values())
    print(f"  {'Score (amber/violet)':<25} {'same/diff':>30} "
          f"{'/'.join(str(score_counts.get(c,0)) for c in ['same','diff']):>20}")

    fiber_counts = Counter(info['same_parent'] for info in edges.values())
    print(f"  {'Fiber (teal/magenta)':<25} {'intra/inter':>30} "
          f"{'/'.join(str(fiber_counts.get(c,0)) for c in [True,False]):>20}")


print(f"\n  DONE.")
print("=" * 80)
