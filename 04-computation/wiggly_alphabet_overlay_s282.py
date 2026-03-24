#!/usr/bin/env python3
"""
wiggly_alphabet_overlay_s282.py — Which wiggly classes generate each metagraph edge?
opus-2026-03-24-S282

For each metagraph edge (i,j), compute the SET of wiggly classes (tiles)
that generate it. This gives an "alphabet signature" per edge.

Questions:
1. How many tiles generate each edge? (all m? some subset?)
2. Which tiles are shared across many edges?
3. Do spine/ribs/sea edges have different alphabet signatures?
4. Is there a pattern by tile position (skip, row, hook)?
5. What is the incidence matrix I[tile][edge] structure?
"""

import sys, time, subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

def build_data(n):
    m_total = comb(n, 2)
    P = [(i,j) for i in range(n) for j in range(i+1, n)]
    r = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
    lines_raw = [l.strip() for l in (r.stdout or '').split('\n')
                 if len(l.strip()) == m_total and all(c in '01' for c in l.strip())]

    def b2a(bits):
        a = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(P):
            if bits & (1 << (m_total-1-k)): a[i][j] = 1
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
    for i, line in enumerate(lines_raw):
        bits = int(line, 2); adj = b2a(bits)
        s = ss(adj); cc = c3(adj); h = H(adj); canon = cp(adj)
        comp = [[1-adj[a][b] if a!=b else 0 for b in range(n)] for a in range(n)]
        cc2 = cp(comp); sc = (canon == cc2)
        cls.append({'cid':i, 'score':s, 'c3':cc, 'H':h, 'sc':sc,
                    'comp_canon':cc2, 'canon':canon})
        hm[(s,cc,h)].append(i); cm[canon] = i
    for cl in cls: cl['comp_cid'] = cm.get(cl['comp_canon'], -1)

    def lk(adj):
        s = ss(adj); cc = c3(adj); h = H(adj)
        cids = hm.get((s,cc,h))
        if not cids: return -1
        if len(cids) == 1: return cids[0]
        return cm.get(cp(adj), -1)

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

    return P, cls, lk, mn, mid, mi, b2a, m_total

def analyze_alphabet(n):
    print(f"\n{'#'*72}")
    print(f"  n = {n}: WIGGLY ALPHABET OVERLAY")
    print(f"{'#'*72}")

    P, cls, lk, mn, V, mi, b2a, m_total = build_data(n)
    m_tiles = comb(n-1, 2)

    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    labels = [chr(65 + i) if i < 26 else chr(65 + i - 26) + "'" for i in range(m_tiles)]

    tile_to_P = {}; base_P = set()
    for k, (i, j) in enumerate(P):
        if j == i + 1: base_P.add(k)
        else:
            for t_idx, (tx, ty) in enumerate(tiles):
                if (ty-1, tx-1) == (i, j):
                    tile_to_P[t_idx] = k; break

    def tiling_to_bits(tb):
        full_bits = 0
        for k, (i, j) in enumerate(P):
            if k not in base_P:
                for ti in range(m_tiles):
                    if tile_to_P.get(ti) == k:
                        bit = (tb >> (m_tiles - 1 - ti)) & 1
                        if bit: full_bits |= (1 << (m_total - 1 - k))
                        break
        return full_bits

    total = 2 ** m_tiles

    # Precompute tiling → merged class
    tiling_class = [0] * total
    for tb in range(total):
        fb = tiling_to_bits(tb)
        adj = b2a(fb)
        cid = lk(adj)
        tiling_class[tb] = mn[cid] if cid >= 0 else -1

    # For each (edge, tile): count how many tilings generate that (edge, tile) pair
    # edge_tile[edge][tile] = count of (tiling, tile) pairs
    edge_tile_count = defaultdict(lambda: defaultdict(int))
    edge_tile_exists = defaultdict(set)  # edge → set of tiles that generate it

    for tb in range(total):
        ma = tiling_class[tb]
        if ma < 0: continue
        for ti in range(m_tiles):
            flipped = tb ^ (1 << (m_tiles - 1 - ti))
            mb = tiling_class[flipped]
            if mb >= 0 and ma != mb:
                edge = (min(ma, mb), max(ma, mb))
                edge_tile_count[edge][ti] += 1
                edge_tile_exists[edge].add(ti)

    all_edges = sorted(edge_tile_exists.keys())
    E = len(all_edges)

    print(f"  Tiles: {m_tiles}, Edges: {E}")
    print(f"\n  TILE MAP:")
    for ti, (x, y) in enumerate(tiles):
        skip = x - y - 1
        print(f"    {labels[ti]}: tile ({x},{y}), skip={skip}")

    # How many tiles generate each edge?
    tiles_per_edge = [len(edge_tile_exists[e]) for e in all_edges]
    print(f"\n  TILES PER EDGE:")
    print(f"    Min: {min(tiles_per_edge)}, Max: {max(tiles_per_edge)}, "
          f"Mean: {np.mean(tiles_per_edge):.2f}")
    dist = Counter(tiles_per_edge)
    for k in sorted(dist.keys()):
        print(f"    {k} tiles: {dist[k]} edges ({100*dist[k]/E:.1f}%)")

    # How many edges does each tile generate?
    edges_per_tile = [0] * m_tiles
    for e, tile_set in edge_tile_exists.items():
        for ti in tile_set:
            edges_per_tile[ti] += 1

    print(f"\n  EDGES PER TILE:")
    for ti in range(m_tiles):
        x, y = tiles[ti]
        skip = x - y - 1
        print(f"    {labels[ti]} ({x},{y}) skip={skip}: {edges_per_tile[ti]}/{E} edges "
              f"({100*edges_per_tile[ti]/E:.1f}%)")

    # Incidence matrix I[tile][edge]
    I = np.zeros((m_tiles, E), dtype=int)
    for ei, e in enumerate(all_edges):
        for ti in edge_tile_exists[e]:
            I[ti][ei] = 1

    # I patterns: which tiles have identical edge sets?
    tile_signatures = {}
    for ti in range(m_tiles):
        sig = tuple(I[ti])
        if sig not in tile_signatures:
            tile_signatures[sig] = []
        tile_signatures[sig].append(ti)

    print(f"\n  TILE EQUIVALENCE CLASSES (tiles with identical edge sets):")
    for sig, tile_group in sorted(tile_signatures.items(), key=lambda x: -len(x[1])):
        tile_labels = [labels[ti] for ti in tile_group]
        tile_descs = [f"({tiles[ti][0]},{tiles[ti][1]})" for ti in tile_group]
        skips = [tiles[ti][0] - tiles[ti][1] - 1 for ti in tile_group]
        n_edges = sum(sig)
        print(f"    {tile_labels}: {n_edges}/{E} edges, skips={skips}")

    # Edge signatures: which edges have the same tile set?
    edge_signatures = {}
    for ei, e in enumerate(all_edges):
        sig = frozenset(edge_tile_exists[e])
        if sig not in edge_signatures:
            edge_signatures[sig] = []
        edge_signatures[sig].append(e)

    print(f"\n  EDGE EQUIVALENCE CLASSES (edges with same tile set):")
    print(f"    {len(edge_signatures)} distinct tile-set patterns for {E} edges")
    for sig, edge_group in sorted(edge_signatures.items(), key=lambda x: -len(x[1]))[:15]:
        tile_labels = sorted([labels[ti] for ti in sig])
        h_pairs = [(mi[a]['H'], mi[b]['H']) for (a,b) in edge_group[:3]]
        print(f"    Tiles {tile_labels} ({len(sig)} tiles): {len(edge_group)} edges, "
              f"e.g. H={h_pairs}")

    # Spine/ribs/sea × alphabet
    spine_tiles = defaultdict(int)
    ribs_tiles = defaultdict(int)
    sea_tiles = defaultdict(int)

    for e in all_edges:
        a, b = e
        sc_a, sc_b = mi[a]['sc'], mi[b]['sc']
        for ti in edge_tile_exists[e]:
            if sc_a and sc_b: spine_tiles[ti] += 1
            elif sc_a != sc_b: ribs_tiles[ti] += 1
            else: sea_tiles[ti] += 1

    n_spine = sum(1 for e in all_edges if mi[e[0]]['sc'] and mi[e[1]]['sc'])
    n_ribs = sum(1 for e in all_edges if mi[e[0]]['sc'] != mi[e[1]]['sc'])
    n_sea = E - n_spine - n_ribs

    print(f"\n  SPINE/RIBS/SEA × WIGGLY CLASS:")
    print(f"  Total: spine={n_spine}, ribs={n_ribs}, sea={n_sea}")
    print(f"  {'Tile':>5} {'skip':>5} {'spine':>6} {'ribs':>6} {'sea':>6} {'total':>6} "
          f"{'spine%':>7} {'ribs%':>7} {'sea%':>7}")
    for ti in range(m_tiles):
        sp = spine_tiles[ti]; ri = ribs_tiles[ti]; se = sea_tiles[ti]
        tot = sp + ri + se
        if tot > 0:
            print(f"  {labels[ti]:>5} {tiles[ti][0]-tiles[ti][1]-1:>5} "
                  f"{sp:>6} {ri:>6} {se:>6} {tot:>6} "
                  f"{100*sp/tot:>6.1f}% {100*ri/tot:>6.1f}% {100*se/tot:>6.1f}%")

    return {
        'n': n, 'E': E, 'm_tiles': m_tiles,
        'tiles_per_edge': tiles_per_edge,
        'edges_per_tile': edges_per_tile,
        'n_tile_equiv': len(tile_signatures),
        'n_edge_equiv': len(edge_signatures),
    }

if __name__ == '__main__':
    print("=" * 72)
    print("  WIGGLY ALPHABET OVERLAY ON THE METAGRAPH")
    print("  opus-2026-03-24-S282")
    print("=" * 72)

    for n in [4, 5, 6]:
        t0 = time.time()
        r = analyze_alphabet(n)
        print(f"\n  Time: {time.time()-t0:.1f}s")

    print("\nDONE.")
