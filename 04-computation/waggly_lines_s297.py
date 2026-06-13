#!/usr/bin/env python3
"""
waggly_lines_s297.py — Waggly lines: connections between NON-wiggly tilings
opus-2026-03-24-S297

WAGGLY LINES connect tilings at Hamming distance ≥ 2.
They fill the GAPS that wiggly lines leave.

ORDER-k WAGGLY: flip exactly k tiles simultaneously.
  k=1: wiggly (already studied)
  k=2: waggly-2 (flip 2 tiles) → C(m,2) classes per tiling
  k=3: waggly-3 → C(m,3) classes
  ...
  k=m: complement (flip all tiles) → 1 class = blue/black line

Blue/black lines ARE waggly-m lines (the maximum order).

STUDY:
1. Which metagraph gaps are filled by waggly-2?
2. How do waggly-2 classes (pairs of tiles) behave?
3. Is there a "waggly hierarchy" — each order fills more gaps?
"""

import sys, time, subprocess
from math import comb, factorial
from itertools import permutations, combinations
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
        cls.append({'cid':i, 'score':s, 'c3':cc, 'H':h, 'sc':sc, 'comp_canon':cc2, 'canon':canon})
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
    tc = [0] * total
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
        tc[tb] = mn[cid] if cid >= 0 else -1
    return V, mi, tc, m_tiles, total, tiles

def analyze_waggly(n):
    print(f"\n{'#'*72}")
    print(f"  n = {n}: WAGGLY LINES (multi-tile flips)")
    print(f"{'#'*72}")

    t0 = time.time()
    V, mi, tc, m, total, tiles = build_data(n)
    labels = [chr(65+i) if i < 26 else chr(65+i-26)+"'" for i in range(m)]
    print(f"  Built in {time.time()-t0:.1f}s. V={V}, m={m}")

    # Wiggly edges (order 1)
    wiggly_edges = set()
    for tb in range(total):
        ma = tc[tb]
        if ma < 0: continue
        for ti in range(m):
            mb = tc[tb ^ (1 << (m-1-ti))]
            if mb >= 0 and ma != mb:
                wiggly_edges.add((min(ma,mb), max(ma,mb)))

    # Waggly-2 edges (order 2 = flip exactly 2 tiles)
    waggly2_edges = set()
    waggly2_self = 0
    waggly2_new = set()  # edges NOT in wiggly

    for tb in range(total):
        ma = tc[tb]
        if ma < 0: continue
        for ti in range(m):
            for tj in range(ti+1, m):
                flipped = tb ^ (1 << (m-1-ti)) ^ (1 << (m-1-tj))
                mb = tc[flipped]
                if mb >= 0:
                    if ma == mb:
                        waggly2_self += 1
                    else:
                        edge = (min(ma,mb), max(ma,mb))
                        waggly2_edges.add(edge)
                        if edge not in wiggly_edges:
                            waggly2_new.add(edge)

    # Complement edges (order m = flip all tiles)
    comp_edges = set()
    for tb in range(total):
        ma = tc[tb]
        if ma < 0: continue
        comp = ((1 << m) - 1) ^ tb
        mb = tc[comp]
        if mb >= 0 and ma != mb:
            comp_edges.add((min(ma,mb), max(ma,mb)))

    # All possible pairs
    total_pairs = V * (V-1) // 2
    wiggly_gaps = total_pairs - len(wiggly_edges)

    print(f"\n  WAGGLY HIERARCHY:")
    print(f"  Order 1 (wiggly):    {len(wiggly_edges):>6} edges ({100*len(wiggly_edges)/total_pairs:.1f}%)")
    print(f"  Order 2 (waggly-2):  {len(waggly2_edges):>6} edges ({100*len(waggly2_edges)/total_pairs:.1f}%)")
    print(f"  Order m (complement):{len(comp_edges):>6} edges ({100*len(comp_edges)/total_pairs:.1f}%)")
    print(f"  Total pairs:         {total_pairs:>6}")

    # GAP FILLING
    gaps_filled_by_2 = waggly2_new  # edges in waggly-2 but NOT in wiggly
    gaps_filled_by_comp = comp_edges - wiggly_edges  # complement-only edges
    remaining_gaps = total_pairs - len(wiggly_edges | waggly2_edges | comp_edges)

    print(f"\n  GAP FILLING:")
    print(f"  Wiggly edges:              {len(wiggly_edges):>6}")
    print(f"  + Waggly-2 new edges:      {len(gaps_filled_by_2):>6}")
    print(f"  + Complement new edges:    {len(gaps_filled_by_comp - waggly2_new):>6}")
    print(f"  = Total connected:         {len(wiggly_edges | waggly2_edges | comp_edges):>6}")
    print(f"  Remaining gaps:            {remaining_gaps:>6} ({100*remaining_gaps/total_pairs:.1f}%)")

    # Is complement ⊆ waggly-2 ∪ wiggly?
    comp_in_waggly2 = comp_edges & waggly2_edges
    comp_only = comp_edges - waggly2_edges - wiggly_edges
    print(f"\n  COMPLEMENT vs WAGGLY-2:")
    print(f"  Complement edges in wiggly:  {len(comp_edges & wiggly_edges)}")
    print(f"  Complement edges in waggly-2:{len(comp_in_waggly2 - wiggly_edges)}")
    print(f"  Complement-only (not in 1 or 2): {len(comp_only)}")
    print(f"  Complement IS subset of waggly-2 ∪ wiggly: "
          f"{len(comp_only) == 0}")

    # Waggly-2 self-loop rate
    total_waggly2_lines = total * comb(m, 2)
    waggly2_sl_rate = waggly2_self / total_waggly2_lines if total_waggly2_lines > 0 else 0
    print(f"\n  WAGGLY-2 SELF-LOOP RATE: {100*waggly2_sl_rate:.2f}%")
    print(f"  (Compare wiggly self-loop rate from S296)")

    # Waggly-2 edge analysis by tile pair
    if n <= 5:
        print(f"\n  WAGGLY-2 CLASSES (tile pairs):")
        waggly2_class_edges = defaultdict(set)
        for tb in range(total):
            ma = tc[tb]
            if ma < 0: continue
            for ti in range(m):
                for tj in range(ti+1, m):
                    flipped = tb ^ (1 << (m-1-ti)) ^ (1 << (m-1-tj))
                    mb = tc[flipped]
                    if mb >= 0 and ma != mb:
                        edge = (min(ma,mb), max(ma,mb))
                        pair_label = f"{labels[ti]}{labels[tj]}"
                        waggly2_class_edges[pair_label].add(edge)

        print(f"  {'Pair':>6} {'Tiles':>12} {'#edges':>7} {'new':>5}")
        for pair_label in sorted(waggly2_class_edges.keys()):
            all_e = waggly2_class_edges[pair_label]
            new_e = all_e - wiggly_edges
            ti = ord(pair_label[0]) - 65
            tj = ord(pair_label[1]) - 65
            print(f"  {pair_label:>6} ({tiles[ti]},{tiles[tj]}):>12 "
                  f"{len(all_e):>7} {len(new_e):>5}")

    # The COMBINED graph: wiggly ∪ waggly-2
    combined_12 = wiggly_edges | waggly2_edges
    print(f"\n  COMBINED GRAPH (wiggly ∪ waggly-2):")
    print(f"  Edges: {len(combined_12)} / {total_pairs} ({100*len(combined_12)/total_pairs:.1f}%)")

    # Is combined_12 = all edges (complete graph)?
    if len(combined_12) == total_pairs:
        print(f"  *** COMPLETE! Orders 1+2 reach ALL class pairs! ***")
    else:
        still_missing = total_pairs - len(combined_12)
        print(f"  Still missing: {still_missing} pairs ({100*still_missing/total_pairs:.1f}%)")

    # ΔH for waggly-2-only edges
    if gaps_filled_by_2:
        dH_vals = [abs(mi[a]['H'] - mi[b]['H']) for (a,b) in gaps_filled_by_2]
        print(f"\n  WAGGLY-2-ONLY EDGES (not in wiggly):")
        print(f"  Count: {len(gaps_filled_by_2)}")
        print(f"  Mean |ΔH|: {np.mean(dH_vals):.2f}")
        print(f"  Range: [{min(dH_vals)}, {max(dH_vals)}]")

    return len(wiggly_edges), len(waggly2_edges), len(comp_edges), total_pairs

if __name__ == '__main__':
    print("=" * 72)
    print("  WAGGLY LINES: MULTI-TILE FLIP CONNECTIONS")
    print("  opus-2026-03-24-S297")
    print("=" * 72)

    for n in [4, 5, 6]:
        analyze_waggly(n)

    print(f"\n{'='*72}")
    print("  THE WAGGLY HIERARCHY")
    print(f"{'='*72}")
    print("""
  Order 1 (WIGGLY): flip 1 tile → m neighbors per tiling
  Order 2 (WAGGLY-2): flip 2 tiles → C(m,2) neighbors
  Order 3 (WAGGLY-3): flip 3 tiles → C(m,3) neighbors
  ...
  Order m (COMPLEMENT): flip all m tiles → 1 neighbor = blue/black line

  Each order fills more gaps. The question: at what order k does
  the combined graph (orders 1..k) become COMPLETE?

  Blue/black (complement) IS order-m. It's the most extreme waggly.
  It connects classes that are Hamming-distance-m apart.
  This is why blue/black lines often fill gaps that wiggly can't.
""")

    print("DONE.")
