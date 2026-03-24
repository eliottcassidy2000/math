#!/usr/bin/env python3
"""
tile_diagonal_classify_s282.py — Classify tiles by diagonal position
opus-2026-03-24-S282

The staircase δ_{n-2} has pin grid coordinates (r,c) where r = skip = a-b-1, c = b.
The grid-symmetry reflection maps (r,c) ↔ (r, n-r-c) — it fixes r and reflects c.

TWO DIAGONAL SYSTEMS:
  MAIN DIAGONAL: c = (n-r)/2. Tiles on this diagonal are singletons (self-paired).
  ANTI-DIAGONAL (strips): r+c = k. These are the "strips" Str(k).

Each tile has:
  r = skip (distance between vertices in base path)
  c = lower vertex (y = b)
  strip = r+c (which anti-diagonal)
  d = c - (n-r)/2 (signed distance from main diagonal)
  partner = (r, n-r-c) (grid-symmetric partner)
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

    cls = []; hm = defaultdict(list); cm = {}
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
        cl0 = cls[cids[0]]; mi[mm] = {'H': cl0['H'], 'sc': cl0['sc']}

    return P, cls, lk, mn, mid, mi, b2a, m_total

def classify_tiles(n):
    """Classify all tiles by diagonal coordinates."""
    m_tiles = comb(n-1, 2)
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))

    labels = [chr(65 + i) if i < 26 else chr(65 + i - 26) + "'" for i in range(m_tiles)]

    tile_info = []
    for ti, (x, y) in enumerate(tiles):
        r = x - y - 1   # skip
        c = y            # lower vertex
        strip = r + c    # anti-diagonal
        d = c - (n - r) / 2.0  # distance from main diagonal
        partner_c = n - r - c
        is_diagonal = (abs(d) < 0.01)  # on the main diagonal
        # Find partner tile index
        partner_tile = None
        for tj, (px, py) in enumerate(tiles):
            pr = px - py - 1; pc = py
            if pr == r and pc == partner_c and tj != ti:
                partner_tile = tj
                break

        tile_info.append({
            'idx': ti, 'label': labels[ti], 'tile': (x, y),
            'r': r, 'c': c, 'strip': strip, 'd': d,
            'is_diagonal': is_diagonal,
            'partner': partner_tile,
        })

    return tiles, labels, tile_info

def analyze(n):
    print(f"\n{'#'*72}")
    print(f"  n = {n}: TILE DIAGONAL CLASSIFICATION")
    print(f"{'#'*72}")

    P, cls, lk, mn, V, mi, b2a, m_total = build_data(n)
    tiles, labels, tile_info = classify_tiles(n)
    m_tiles = len(tiles)

    # Display classification
    print(f"\n  {'Label':>6} {'Tile':>8} {'r(skip)':>8} {'c(y)':>6} {'strip':>6} "
          f"{'d':>6} {'diag?':>6} {'partner':>8}")
    for t in tile_info:
        partner_label = labels[t['partner']] if t['partner'] is not None else '—'
        print(f"  {t['label']:>6} ({t['tile'][0]},{t['tile'][1]}):>8 "
              f"{t['r']:>8} {t['c']:>6} {t['strip']:>6} "
              f"{t['d']:>+6.1f} {'YES' if t['is_diagonal'] else 'no':>6} {partner_label:>8}")

    # Group by strip (anti-diagonal)
    by_strip = defaultdict(list)
    for t in tile_info:
        by_strip[t['strip']].append(t)

    print(f"\n  BY STRIP (anti-diagonal r+c = k):")
    for k in sorted(by_strip.keys()):
        tiles_in_strip = by_strip[k]
        tile_labels = [t['label'] for t in tiles_in_strip]
        d_vals = [t['d'] for t in tiles_in_strip]
        print(f"    Strip {k}: {tile_labels}, d = {d_vals}")

    # Group by distance from diagonal
    by_d = defaultdict(list)
    for t in tile_info:
        by_d[abs(t['d'])].append(t)

    print(f"\n  BY DISTANCE FROM MAIN DIAGONAL:")
    for d_abs in sorted(by_d.keys()):
        tiles_at_d = by_d[d_abs]
        tile_labels = [t['label'] for t in tiles_at_d]
        print(f"    |d| = {d_abs:.1f}: {tile_labels}")

    # Now compute edge sets per tile and analyze by diagonal class
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
    tiling_class = [0] * total
    for tb in range(total):
        fb = tiling_to_bits(tb)
        adj = b2a(fb); cid = lk(adj)
        tiling_class[tb] = mn[cid] if cid >= 0 else -1

    # Edge sets per tile
    edge_sets = [set() for _ in range(m_tiles)]
    sl_counts = [0] * m_tiles
    for tb in range(total):
        ma = tiling_class[tb]
        if ma < 0: continue
        for ti in range(m_tiles):
            flipped = tb ^ (1 << (m_tiles - 1 - ti))
            mb = tiling_class[flipped]
            if mb >= 0:
                if ma == mb:
                    sl_counts[ti] += 1
                else:
                    edge_sets[ti].add((min(ma, mb), max(ma, mb)))

    # Analysis by diagonal class
    print(f"\n  EDGES AND SL% BY DIAGONAL POSITION:")
    print(f"  {'Label':>6} {'(x,y)':>8} {'r':>3} {'c':>3} {'strip':>6} "
          f"{'d':>6} {'#edges':>7} {'SL%':>7}")
    for t in tile_info:
        ti = t['idx']
        n_edges = len(edge_sets[ti])
        tot = sl_counts[ti] + sum(1 for tb in range(total)
                                   if tiling_class[tb] >= 0
                                   for _ in [1]  # placeholder
                                   ) // m_tiles  # approx
        sl_pct = 100 * sl_counts[ti] / total if total > 0 else 0
        print(f"  {t['label']:>6} ({t['tile'][0]},{t['tile'][1]}):>8 "
              f"{t['r']:>3} {t['c']:>3} {t['strip']:>6} {t['d']:>+6.1f} "
              f"{n_edges:>7} {sl_pct:>6.1f}%")

    # Check: do tiles at the same |d| have the same edge set?
    print(f"\n  EDGE SET EQUALITY BY |d|:")
    for d_abs in sorted(by_d.keys()):
        tiles_at_d = by_d[d_abs]
        if len(tiles_at_d) > 1:
            all_same = all(edge_sets[tiles_at_d[0]['idx']] == edge_sets[t['idx']]
                          for t in tiles_at_d)
            print(f"    |d|={d_abs:.1f}: {[t['label'] for t in tiles_at_d]} "
                  f"→ same edge set: {all_same}")

    # Check: do tiles on the same strip have the same edge set?
    print(f"\n  EDGE SET EQUALITY BY STRIP:")
    for k in sorted(by_strip.keys()):
        tiles_in_strip = by_strip[k]
        if len(tiles_in_strip) > 1:
            all_same = all(edge_sets[tiles_in_strip[0]['idx']] == edge_sets[t['idx']]
                          for t in tiles_in_strip)
            print(f"    Strip {k}: {[t['label'] for t in tiles_in_strip]} "
                  f"→ same edge set: {all_same}")

    # Check: do tiles at the same r (skip) have the same edge set?
    by_r = defaultdict(list)
    for t in tile_info: by_r[t['r']].append(t)

    print(f"\n  EDGE SET EQUALITY BY SKIP (r):")
    for r_val in sorted(by_r.keys()):
        tiles_at_r = by_r[r_val]
        if len(tiles_at_r) > 1:
            all_same = all(edge_sets[tiles_at_r[0]['idx']] == edge_sets[t['idx']]
                          for t in tiles_at_r)
            print(f"    r={r_val}: {[t['label'] for t in tiles_at_r]} "
                  f"→ same edge set: {all_same}")

if __name__ == '__main__':
    print("=" * 72)
    print("  TILE DIAGONAL CLASSIFICATION")
    print("  opus-2026-03-24-S282")
    print("=" * 72)

    for n in [4, 5, 6, 7]:
        t0 = time.time()
        analyze(n)
        print(f"\n  Time: {time.time()-t0:.1f}s")

    print("\nDONE.")
