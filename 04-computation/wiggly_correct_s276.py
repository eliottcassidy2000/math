#!/usr/bin/env python3
"""
wiggly_correct_s276.py — Wiggly classes with CORRECT definition + connections
opus-2026-03-24-S276

WIGGLY LINE = clicking one TILE in the explorer = flipping one non-base-path arc.
WIGGLY CLASS = which tile was clicked. There are C(n-1,2) classes.

This script computes per-wiggly-class:
1. Self-loop rate (blue-self + black-self) vs opaque edges
2. Connection to hook lengths h(i,j) = 2(k-i-j)-1
3. Which metagraph edges each class generates (spine/ribs/sea?)
4. Per-class H change distribution (|ΔH| per flip)
5. Per-class contribution to the Burnside transition formula
"""

import sys, time, subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter

sys.stdout.reconfigure(line_buffering=True)

def build_data(n):
    """Build tournament class data with tiling model."""
    m_total = comb(n, 2)  # total arcs
    m_tiles = comb(n-1, 2)  # tiles (non-base-path arcs)

    # Base path: n → n-1 → ... → 1
    # Tiles: (x,y) with x-y ≥ 2, x,y ∈ {1,...,n}
    # The base-path arcs are (k, k-1) for k=2,...,n

    # Enumerate tiles in explorer order: for y=1..n-2: for x=n down to y+2
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))

    # All arcs: (i,j) with i > j, using 1-indexed vertices
    all_arcs = [(i, j) for i in range(1, n+1) for j in range(1, i)]
    base_path_arcs = [(k, k-1) for k in range(2, n+1)]

    # Pair list for gentourng format (0-indexed, i < j)
    P = [(i, j) for i in range(n) for j in range(i+1, n)]

    # Map from tile (x,y) 1-indexed to arc position in P (0-indexed)
    # Arc (x,y) with x > y maps to pair (y-1, x-1) with y-1 < x-1
    tile_to_arc_idx = {}
    for t_idx, (x, y) in enumerate(tiles):
        # 1-indexed (x,y) → 0-indexed pair (y-1, x-1) where y-1 < x-1
        pair = (y-1, x-1)
        arc_idx = P.index(pair)
        tile_to_arc_idx[t_idx] = arc_idx

    # Get class reps from nauty
    r = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
    m = len(P)
    lines = [l.strip() for l in (r.stdout or '').split('\n')
             if len(l.strip()) == m and all(c in '01' for c in l.strip())]

    def b2a(bits):
        a = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(P):
            if bits & (1 << (m-1-k)): a[i][j] = 1
            else: a[j][i] = 1
        return a
    def ss(a): return tuple(sorted(sum(a[i][j] for j in range(n)) for i in range(n)))
    def c3(a):
        c = 0
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if a[i][j] and a[j][k] and a[k][i]: c += 1
                    if a[i][k] and a[k][j] and a[j][i]: c += 1
        return c
    def H(a):
        dp = {}
        for v in range(n): dp[(1<<v,v)] = 1
        for S in range(1, 1<<n):
            for v in range(n):
                if not (S&(1<<v)): continue
                val = dp.get((S,v), 0)
                if val == 0: continue
                for u in range(n):
                    if S&(1<<u): continue
                    if a[v][u]: dp[(S|(1<<u),u)] = dp.get((S|(1<<u),u),0) + val
        return sum(dp.get(((1<<n)-1,v),0) for v in range(n))
    def cp(a):
        sc = [sum(a[i][j] for j in range(n)) for i in range(n)]
        sg = defaultdict(list)
        for v in range(n): sg[sc[v]].append(v)
        gs = [sg[s] for s in sorted(sg.keys())]
        best = None; cnt = 0
        def gp(gs2):
            if not gs2: yield []; return
            for p in permutations(gs2[0]):
                for rest in gp(gs2[1:]): yield list(p)+rest
        for p in gp(gs):
            f = tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or f < best: best = f
            cnt += 1
            if cnt > 500000: break
        return best

    cls = []
    hm = defaultdict(list); cm = {}
    for i, line in enumerate(lines):
        bits = int(line, 2); adj = b2a(bits)
        s = ss(adj); cc = c3(adj); h = H(adj); canon = cp(adj)
        comp = [[1-adj[a][b] if a!=b else 0 for b in range(n)] for a in range(n)]
        cc2 = cp(comp); sc = (canon == cc2)
        cls.append({'cid':i, 'adj':adj, 'score':s, 'c3':cc, 'H':h,
                    'sc':sc, 'comp_canon':cc2, 'canon':canon})
        hm[(s,cc,h)].append(i); cm[canon] = i
    for cl in cls: cl['comp_cid'] = cm.get(cl['comp_canon'], -1)

    def lk(adj):
        s = ss(adj); cc = c3(adj); h = H(adj)
        cids = hm.get((s,cc,h))
        if not cids: return -1
        if len(cids) == 1: return cids[0]
        return cm.get(cp(adj), -1)

    # Merge complement pairs
    mn = {}; mid = 0
    for cl in cls:
        ci = cl['cid']
        if ci in mn: continue
        if cl['sc']: mn[ci] = mid; mid += 1
        else:
            mn[ci] = mid; comp = cl['comp_cid']
            if comp >= 0 and comp != ci: mn[comp] = mid
            mid += 1
    mi = {}
    for mm in range(mid):
        cids = [cl['cid'] for cl in cls if mn[cl['cid']] == mm]
        cl0 = cls[cids[0]]
        mi[mm] = {'H': cl0['H'], 'sc': cl0['sc']}

    return {
        'n': n, 'm': m, 'm_tiles': m_tiles, 'tiles': tiles,
        'tile_to_arc': tile_to_arc_idx,
        'P': P, 'cls': cls, 'lk': lk, 'mn': mn, 'V_merged': mid, 'mi': mi,
        'b2a': b2a,
    }

def analyze_wiggly_classes(data):
    """Per-wiggly-class analysis."""
    n = data['n']; m = data['m']; tiles = data['tiles']
    lk = data['lk']; mn = data['mn']; mi = data['mi']
    b2a = data['b2a']; V = data['V_merged']
    tile_to_arc = data['tile_to_arc']
    m_tiles = data['m_tiles']

    print(f"\n{'#'*72}")
    print(f"  n = {n}: WIGGLY CLASS ANALYSIS ({m_tiles} classes, {2**m_tiles} tilings)")
    print(f"{'#'*72}")

    # Label classes A, B, C, ...
    labels = [chr(65 + i) for i in range(m_tiles)]  # A, B, C, ...

    print(f"\n  TILE → WIGGLY CLASS MAPPING:")
    print(f"  {'Class':>6} {'Tile(x,y)':>10} {'Skip':>5} {'Hook':>5}")
    for t_idx, (x, y) in enumerate(tiles):
        skip = x - y - 1
        # Hook in staircase δ_{n-2}: h(r,c) = 2(k-r-c)-1 where k=n-2
        # With r = x-y-1, c = y-1 (mapping from tile to staircase cell)
        r = x - y - 1  # row in staircase (0-indexed)
        c = y - 1       # column in staircase (0-indexed)
        k = n - 2
        hook = 2*(k - r - c) - 1 if k - r - c > 0 else 1
        print(f"  {labels[t_idx]:>6} ({x},{y}):>10 {skip:>5} {hook:>5}")

    # Iterate over all 2^{m_tiles} tilings
    # A tiling = a binary vector of length m_tiles (the tile values)
    # Combined with the fixed base path, this gives a full tournament

    per_class_sl = [0] * m_tiles      # self-loops per class
    per_class_opaque = [0] * m_tiles   # class-changing flips per class
    per_class_dH = [[] for _ in range(m_tiles)]  # |ΔH| values
    per_class_edges = [set() for _ in range(m_tiles)]  # metagraph edges found
    per_class_spine = [0] * m_tiles    # SC-SC edges
    per_class_ribs = [0] * m_tiles     # SC-NS edges
    per_class_sea = [0] * m_tiles      # NS-NS edges

    total_tilings = 2 ** m_tiles

    for tiling_bits in range(total_tilings):
        # Build full tournament from tiling
        # Base path arcs are fixed: (k, k-1) → k beats k-1 → bit=1 at position (k-1-1, k-1) in P
        # Actually, we need to construct the full m-bit tournament vector
        full_bits = 0
        for k, (i, j) in enumerate(data['P']):
            # Is this a base-path arc? i < j (0-indexed)
            # Base path: j → i when j = i+1 (j beats i, i.e., i+1 → i)
            # In P format: pair (i, i+1) with i < i+1
            # The arc direction: if j = i+1, then j→i means bit = 0
            # (because bit=1 means i→j, bit=0 means j→i)
            if j == i + 1:
                # Base path: j→i, so bit = 0
                pass  # bit stays 0
            else:
                # This is a tile. Find its tile index
                # Tile (x,y) with x=j+1, y=i+1 (converting 0-indexed to 1-indexed)
                # Wait: P has (i,j) with i<j (0-indexed)
                # The tile is (j+1, i+1) in 1-indexed with (j+1)-(i+1) = j-i ≥ 2
                tile_1idx = (j+1, i+1)  # 1-indexed
                # Find this tile in the tiles list
                try:
                    t_idx = tiles.index(tile_1idx)
                except ValueError:
                    # Try reversed
                    tile_1idx = (i+1, j+1)
                    try:
                        t_idx = tiles.index(tile_1idx)
                    except ValueError:
                        continue

                # Get bit from tiling
                bit = (tiling_bits >> (m_tiles - 1 - t_idx)) & 1
                if bit:
                    full_bits |= (1 << (m - 1 - k))

        adj = b2a(full_bits)
        cid = lk(adj)
        if cid < 0: continue
        ma = mn[cid]
        H_a = mi[ma]['H']

        # For each wiggly class (tile), flip it and check
        for t_idx in range(m_tiles):
            flipped_bits = tiling_bits ^ (1 << (m_tiles - 1 - t_idx))

            # Rebuild tournament with flipped tile
            full_bits_f = 0
            for k, (i, j) in enumerate(data['P']):
                if j == i + 1:
                    pass
                else:
                    tile_1idx = (j+1, i+1)
                    try:
                        ti = tiles.index(tile_1idx)
                    except ValueError:
                        tile_1idx = (i+1, j+1)
                        try:
                            ti = tiles.index(tile_1idx)
                        except ValueError:
                            continue
                    bit = (flipped_bits >> (m_tiles - 1 - ti)) & 1
                    if bit:
                        full_bits_f |= (1 << (m - 1 - k))

            adj_f = b2a(full_bits_f)
            cid_f = lk(adj_f)
            if cid_f < 0: continue
            mb = mn[cid_f]

            if ma == mb:
                per_class_sl[t_idx] += 1
            else:
                per_class_opaque[t_idx] += 1
                per_class_edges[t_idx].add((min(ma, mb), max(ma, mb)))
                H_b = mi[mb]['H']
                per_class_dH[t_idx].append(abs(H_a - H_b))

                # Classify edge type
                sc_a = mi[ma]['sc']
                sc_b = mi[mb]['sc']
                if sc_a and sc_b:
                    per_class_spine[t_idx] += 1
                elif sc_a != sc_b:
                    per_class_ribs[t_idx] += 1
                else:
                    per_class_sea[t_idx] += 1

    # Results
    print(f"\n  PER-CLASS STATISTICS:")
    print(f"  {'Cls':>4} {'Tile':>6} {'Skip':>5} {'Hook':>5} {'SL%':>7} "
          f"{'#Edges':>7} {'<|ΔH|>':>7} {'Spine%':>7} {'Ribs%':>7} {'Sea%':>7}")

    for t_idx in range(m_tiles):
        x, y = tiles[t_idx]
        skip = x - y - 1
        r = x - y - 1; c = y - 1; k = n - 2
        hook = 2*(k - r - c) - 1 if k - r - c > 0 else 1
        total = per_class_sl[t_idx] + per_class_opaque[t_idx]
        sl_pct = 100 * per_class_sl[t_idx] / total if total > 0 else 0
        n_edges = len(per_class_edges[t_idx])
        avg_dH = sum(per_class_dH[t_idx]) / len(per_class_dH[t_idx]) if per_class_dH[t_idx] else 0
        opq = per_class_opaque[t_idx]
        sp_pct = 100 * per_class_spine[t_idx] / opq if opq > 0 else 0
        ri_pct = 100 * per_class_ribs[t_idx] / opq if opq > 0 else 0
        se_pct = 100 * per_class_sea[t_idx] / opq if opq > 0 else 0

        print(f"  {labels[t_idx]:>4} ({x},{y}):>6 {skip:>5} {hook:>5} {sl_pct:>6.1f}% "
              f"{n_edges:>7} {avg_dH:>7.2f} {sp_pct:>6.1f}% {ri_pct:>6.1f}% {se_pct:>6.1f}%")

    # Correlation: hook vs SL rate
    import numpy as np
    hooks = []
    sl_rates = []
    for t_idx in range(m_tiles):
        x, y = tiles[t_idx]
        r = x - y - 1; c = y - 1; k = n - 2
        hook = 2*(k - r - c) - 1 if k - r - c > 0 else 1
        hooks.append(hook)
        total = per_class_sl[t_idx] + per_class_opaque[t_idx]
        sl_rates.append(per_class_sl[t_idx] / total if total > 0 else 0)

    if len(hooks) >= 3:
        r_hook = np.corrcoef(hooks, sl_rates)[0, 1]
        print(f"\n  Correlation(hook, SL_rate) = {r_hook:.4f}")

    skips = [tiles[i][0] - tiles[i][1] - 1 for i in range(m_tiles)]
    r_skip = np.corrcoef(skips, sl_rates)[0, 1] if len(skips) >= 3 else 0
    print(f"  Correlation(skip, SL_rate) = {r_skip:.4f}")

    # Which edges does EVERY class see?
    common_edges = per_class_edges[0]
    for t_idx in range(1, m_tiles):
        common_edges = common_edges & per_class_edges[t_idx]
    all_edges = set()
    for t_idx in range(m_tiles):
        all_edges |= per_class_edges[t_idx]

    print(f"\n  EDGE COVERAGE:")
    print(f"  Total metagraph edges: {len(all_edges)}")
    print(f"  Edges in ALL classes: {len(common_edges)} ({100*len(common_edges)/len(all_edges):.1f}%)")
    for t_idx in range(m_tiles):
        print(f"    Class {labels[t_idx]}: {len(per_class_edges[t_idx])} edges")

if __name__ == '__main__':
    print("=" * 72)
    print("  WIGGLY CLASSES — CORRECT DEFINITION")
    print("  opus-2026-03-24-S276")
    print("=" * 72)

    for n in [4, 5]:
        t0 = time.time()
        data = build_data(n)
        analyze_wiggly_classes(data)
        print(f"\n  Time: {time.time()-t0:.1f}s")

    print("\nDONE.")
