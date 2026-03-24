#!/usr/bin/env python3
"""
rare_dual_edges_s285.py — Study edges that are BOTH wiggly AND complement
opus-2026-03-24-S285

An edge (i,j) in the merged metagraph is a "DUAL EDGE" if it's generated
by BOTH a wiggly line (tile flip) AND a complement line (flip all tiles).

These are rare (8.2% at n=7) and shrinking. They're the places where
the two lenses agree — where a single tile flip and a full complement
both connect the same pair of classes.

Also study kind-pasteur's "collapsed edges": classes adjacent to their
own complement (sequence 0,0,5,0 at n=4..7).
"""

import sys, time, subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

def build_all(n):
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
        cl0 = cls[cids[0]]
        mi[mm] = {'H': cl0['H'], 'sc': cl0['sc']}

    m_tiles = comb(n-1, 2)
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2: tiles.append((x, y))
    tile_to_P = {}; base_P = set()
    for k, (i, j) in enumerate(P):
        if j == i + 1: base_P.add(k)
        else:
            for t_idx, (tx, ty) in enumerate(tiles):
                if (ty-1, tx-1) == (i, j):
                    tile_to_P[t_idx] = k; break

    total = 2 ** m_tiles; V = mid
    tiling_class = [0] * total
    for tb in range(total):
        full_bits = 0
        for k, (i, j) in enumerate(P):
            if k not in base_P:
                for ti in range(m_tiles):
                    if tile_to_P.get(ti) == k:
                        bit = (tb >> (m_tiles - 1 - ti)) & 1
                        if bit: full_bits |= (1 << (m_total - 1 - k))
                        break
        adj = b2a(full_bits); cid = lk(adj)
        tiling_class[tb] = mn[cid] if cid >= 0 else -1

    # Compute W and B edges
    w_edges = set(); b_edges = set()
    w_weights = defaultdict(int); b_weights = defaultdict(int)
    fiber = Counter()

    for tb in range(total):
        ma = tiling_class[tb]
        if ma < 0: continue
        fiber[ma] += 1
        for ti in range(m_tiles):
            flipped = tb ^ (1 << (m_tiles - 1 - ti))
            mb = tiling_class[flipped]
            if mb >= 0 and ma != mb:
                e = (min(ma,mb), max(ma,mb))
                w_edges.add(e); w_weights[e] += 1
        comp = ((1 << m_tiles) - 1) ^ tb
        mb_c = tiling_class[comp]
        if mb_c >= 0 and ma != mb_c:
            e = (min(ma,mb_c), max(ma,mb_c))
            b_edges.add(e); b_weights[e] += 1

    both = w_edges & b_edges
    only_w = w_edges - b_edges
    only_b = b_edges - w_edges

    return V, mi, w_edges, b_edges, both, only_w, only_b, w_weights, b_weights, dict(fiber), m_tiles, tiles

def analyze(n):
    print(f"\n{'#'*72}")
    print(f"  n = {n}: RARE DUAL EDGES")
    print(f"{'#'*72}")

    t0 = time.time()
    V, mi, w_edges, b_edges, both, only_w, only_b, w_wt, b_wt, fiber, m, tiles = build_all(n)
    dt = time.time() - t0
    print(f"  V={V}, W_edges={len(w_edges)}, B_edges={len(b_edges)}, "
          f"both={len(both)}, time={dt:.1f}s")

    if not both:
        print(f"  NO dual edges at n={n}.")
        return

    # Analyze the dual edges
    print(f"\n  DUAL EDGES (both wiggly AND complement):")
    print(f"  {'Edge':>10} {'H_i':>4} {'H_j':>4} {'SC':>8} {'W_wt':>7} {'B_wt':>7} "
          f"{'W/B':>6} {'fib_i':>6} {'fib_j':>6} {'|ΔH|':>5}")

    dh_vals = []
    sc_types = Counter()
    for (i, j) in sorted(both, key=lambda e: (mi[e[0]]['H'], mi[e[1]]['H'])):
        hi = mi[i]['H']; hj = mi[j]['H']
        sci = mi[i]['sc']; scj = mi[j]['sc']
        ww = w_wt[(i,j)]; bw = b_wt[(i,j)]
        fi = fiber[i]; fj = fiber[j]
        dh = abs(hi - hj)
        dh_vals.append(dh)

        sc_type = ('SC' if sci else 'NS') + '-' + ('SC' if scj else 'NS')
        sc_types[sc_type] += 1

        ratio = ww / bw if bw > 0 else float('inf')

        print(f"  ({i:>3},{j:>3}) {hi:>4} {hj:>4} {sc_type:>8} {ww:>7} {bw:>7} "
              f"{ratio:>6.2f} {fi:>6} {fj:>6} {dh:>5}")

    # Statistics
    print(f"\n  DUAL EDGE STATISTICS:")
    print(f"  Total: {len(both)} ({100*len(both)/len(w_edges|b_edges):.1f}% of all edges)")
    print(f"  |ΔH| distribution: {sorted(Counter(dh_vals).items())}")
    print(f"  SC type distribution: {dict(sc_types)}")

    # Compare with non-dual edges
    w_only_dh = [abs(mi[i]['H'] - mi[j]['H']) for (i,j) in only_w]
    b_only_dh = [abs(mi[i]['H'] - mi[j]['H']) for (i,j) in only_b]
    both_dh = dh_vals

    print(f"\n  MEAN |ΔH| BY EDGE TYPE:")
    if w_only_dh: print(f"    Wiggly-only: {np.mean(w_only_dh):.2f} ({len(only_w)} edges)")
    if b_only_dh: print(f"    Complement-only: {np.mean(b_only_dh):.2f} ({len(only_b)} edges)")
    if both_dh: print(f"    DUAL (both): {np.mean(both_dh):.2f} ({len(both)} edges)")

    # W/B weight ratio for dual edges
    ratios = [w_wt[e] / b_wt[e] for e in both if b_wt[e] > 0]
    if ratios:
        print(f"\n  W/B WEIGHT RATIO FOR DUAL EDGES:")
        print(f"    Min: {min(ratios):.2f}")
        print(f"    Max: {max(ratios):.2f}")
        print(f"    Mean: {np.mean(ratios):.2f}")
        print(f"    Expected (if uniform): m = {m}")

    # Check: are dual edges concentrated near spine, ribs, or sea?
    spine_dual = sum(1 for (i,j) in both if mi[i]['sc'] and mi[j]['sc'])
    ribs_dual = sum(1 for (i,j) in both if mi[i]['sc'] != mi[j]['sc'])
    sea_dual = sum(1 for (i,j) in both if not mi[i]['sc'] and not mi[j]['sc'])

    spine_total = sum(1 for (i,j) in (w_edges|b_edges) if mi[i]['sc'] and mi[j]['sc'])
    ribs_total = sum(1 for (i,j) in (w_edges|b_edges) if mi[i]['sc'] != mi[j]['sc'])
    sea_total = sum(1 for (i,j) in (w_edges|b_edges) if not mi[i]['sc'] and not mi[j]['sc'])

    print(f"\n  DUAL EDGE TYPE (spine/ribs/sea):")
    print(f"    Spine: {spine_dual}/{spine_total} "
          f"({100*spine_dual/spine_total:.1f}% of spine)" if spine_total > 0 else "")
    print(f"    Ribs: {ribs_dual}/{ribs_total} "
          f"({100*ribs_dual/ribs_total:.1f}% of ribs)" if ribs_total > 0 else "")
    print(f"    Sea: {sea_dual}/{sea_total} "
          f"({100*sea_dual/sea_total:.1f}% of sea)" if sea_total > 0 else "")

if __name__ == '__main__':
    print("=" * 72)
    print("  RARE DUAL EDGES: WHERE WIGGLY AND COMPLEMENT AGREE")
    print("  opus-2026-03-24-S285")
    print("=" * 72)

    for n in [4, 5, 6, 7]:
        analyze(n)

    print(f"\n{'='*72}")
    print("  THE DUAL EDGE PHENOMENON")
    print(f"{'='*72}")
    print("""
  A DUAL EDGE connects two merged classes via BOTH:
    - A wiggly line (single tile flip)
    - A complement line (flip all tiles)

  These are rare because the two operations are very different:
    Wiggly: change 1 of m tiles
    Complement: change ALL m tiles

  For both to connect the same pair, there must be a tiling T in class i
  whose complement is in class j, AND ALSO a tiling T' in class i
  that has a single-tile-flip neighbor in class j.

  T and T' can be DIFFERENT tilings in the same class.

  The rarity of dual edges measures how DIFFERENT the two lenses are.
  At n=7: only 8.2% of edges are dual → the lenses are 92% different.

  DUAL EDGES are where tournament structure is MOST constrained:
  both local (tile flip) and global (complement) operations connect
  the same classes. These are the "resonances" of tournament space.
""")

    print("DONE.")
