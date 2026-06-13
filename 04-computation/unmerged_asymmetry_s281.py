#!/usr/bin/env python3
"""
unmerged_asymmetry_s281.py — The asymmetry hidden inside W's symmetry
opus-2026-03-24-S281

W[i][j] = W[j][i] for the MERGED metagraph. But within each merged edge,
the constituent UNMERGED classes connect asymmetrically.

If merged i = {A, A^c} and merged j = {B, B^c}, there are four sub-flows:
  A → B,  A → B^c,  A^c → B,  A^c → B^c

These four weights need not be equal. The merged W sums them:
  W[i][j] = w(A→B) + w(A→B^c) + w(A^c→B) + w(A^c→B^c)

This script computes the sub-flow structure and finds the asymmetry.
"""

import sys, time, subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

def build_data(n):
    """Build BOTH merged and unmerged class data."""
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

    # Build merge map AND keep track of which unmerged class each is
    mn = {}; mid = 0
    merge_contents = {}  # merged_id → list of unmerged cids
    for cl in cls:
        ci = cl['cid']
        if ci in mn: continue
        if cl['sc']:
            mn[ci] = mid
            merge_contents[mid] = [ci]
            mid += 1
        else:
            mn[ci] = mid
            comp = cl['comp_cid']
            merge_contents[mid] = [ci]
            if comp >= 0 and comp != ci:
                mn[comp] = mid
                merge_contents[mid].append(comp)
            mid += 1

    mi = {}
    for mm in range(mid):
        cids = merge_contents[mm]
        cl0 = cls[cids[0]]
        mi[mm] = {'H': cl0['H'], 'sc': cl0['sc'],
                  'unmerged': cids, 'n_unmerged': len(cids)}

    return P, cls, lk, mn, mid, mi, b2a, m_total, merge_contents

def analyze_asymmetry(n):
    """Compute sub-flow structure within merged edges."""
    print(f"\n{'#'*72}")
    print(f"  n = {n}: UNMERGED ASYMMETRY")
    print(f"{'#'*72}")

    P, cls, lk, mn, V, mi, b2a, m_total, merge_contents = build_data(n)

    m_tiles = comb(n-1, 2)
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))

    tile_to_P = {}
    base_P = set()
    for k, (i, j) in enumerate(P):
        if j == i + 1:
            base_P.add(k)
        else:
            for t_idx, (tx, ty) in enumerate(tiles):
                if (ty-1, tx-1) == (i, j):
                    tile_to_P[t_idx] = k
                    break

    def tiling_to_bits(tb):
        full_bits = 0
        for k, (i, j) in enumerate(P):
            if k not in base_P:
                for ti in range(m_tiles):
                    if tile_to_P.get(ti) == k:
                        bit = (tb >> (m_tiles - 1 - ti)) & 1
                        if bit:
                            full_bits |= (1 << (m_total - 1 - k))
                        break
        return full_bits

    # Precompute tiling → UNMERGED class
    total = 2 ** m_tiles
    tiling_unmerged = [0] * total
    tiling_merged = [0] * total

    for tb in range(total):
        fb = tiling_to_bits(tb)
        adj = b2a(fb)
        cid = lk(adj)
        tiling_unmerged[tb] = cid
        tiling_merged[tb] = mn[cid] if cid >= 0 else -1

    # Compute UNMERGED wiggly weight matrix
    V_unmerged = len(cls)
    W_u = np.zeros((V_unmerged, V_unmerged), dtype=int)

    for tb in range(total):
        cid_a = tiling_unmerged[tb]
        if cid_a < 0: continue
        for ti in range(m_tiles):
            flipped = tb ^ (1 << (m_tiles - 1 - ti))
            cid_b = tiling_unmerged[flipped]
            if cid_b >= 0:
                W_u[cid_a][cid_b] += 1

    # Check: is W_u symmetric?
    w_u_sym = np.allclose(W_u, W_u.T)
    print(f"\n  Unmerged W symmetric: {w_u_sym}")

    # For each MERGED edge, decompose into sub-flows
    print(f"\n  SUB-FLOW DECOMPOSITION:")
    print(f"  (For NS-NS edges: 4 sub-flows. For SC-*: 1-2 sub-flows.)")

    # Collect all merged edges with sub-flows
    merged_edges = set()
    for i in range(V):
        for j in range(i+1, V):
            # Check if any sub-flow exists
            has_flow = False
            for a in merge_contents[i]:
                for b in merge_contents[j]:
                    if W_u[a][b] > 0 or W_u[b][a] > 0:
                        has_flow = True
            if has_flow:
                merged_edges.add((i, j))

    # Classify by type
    sc_sc = [(i,j) for (i,j) in merged_edges if mi[i]['sc'] and mi[j]['sc']]
    sc_ns = [(i,j) for (i,j) in merged_edges if mi[i]['sc'] != mi[j]['sc']]
    ns_ns = [(i,j) for (i,j) in merged_edges if not mi[i]['sc'] and not mi[j]['sc']]

    print(f"\n  SC-SC edges: {len(sc_sc)}")
    print(f"  SC-NS edges: {len(sc_ns)}")
    print(f"  NS-NS edges: {len(ns_ns)}")

    # Analyze NS-NS edges in detail (4 sub-flows each)
    if ns_ns:
        print(f"\n  NS-NS EDGES (4 sub-flows each):")
        print(f"  {'Edge':>10} {'H_i':>4} {'H_j':>4} {'A→B':>6} {'A→B^c':>7} "
              f"{'A^c→B':>7} {'A^c→B^c':>8} {'total':>6} {'asymm':>6}")

        asymmetries = []
        for (i, j) in sorted(ns_ns, key=lambda e: (mi[e[0]]['H'], mi[e[1]]['H']))[:30]:
            cs_i = merge_contents[i]  # [A, A^c]
            cs_j = merge_contents[j]  # [B, B^c]
            if len(cs_i) < 2 or len(cs_j) < 2: continue

            A, Ac = cs_i[0], cs_i[1]
            B, Bc = cs_j[0], cs_j[1]

            ab = W_u[A][B]
            abc = W_u[A][Bc]
            acb = W_u[Ac][B]
            acbc = W_u[Ac][Bc]
            total_w = ab + abc + acb + acbc

            # Asymmetry: how uneven are the four sub-flows?
            vals = [ab, abc, acb, acbc]
            mean_v = np.mean(vals)
            asymm = np.std(vals) / mean_v if mean_v > 0 else 0

            asymmetries.append(asymm)

            print(f"  ({i:>3},{j:>3}) {mi[i]['H']:>4} {mi[j]['H']:>4} "
                  f"{ab:>6} {abc:>7} {acb:>7} {acbc:>8} {total_w:>6} {asymm:>6.3f}")

        if asymmetries:
            print(f"\n  Asymmetry stats: mean={np.mean(asymmetries):.3f}, "
                  f"max={max(asymmetries):.3f}, "
                  f"zero={sum(1 for a in asymmetries if a < 0.01)}/{len(asymmetries)}")

    # SC-NS edges (2 sub-flows: SC→C and SC→C^c)
    if sc_ns:
        print(f"\n  SC-NS EDGES (2 sub-flows each):")
        print(f"  {'Edge':>10} {'H_i':>4} {'H_j':>4} {'SC→C':>6} {'SC→C^c':>8} "
              f"{'total':>6} {'ratio':>7}")

        ratios = []
        for (i, j) in sorted(sc_ns, key=lambda e: (mi[e[0]]['H'], mi[e[1]]['H']))[:20]:
            # Determine which is SC
            if mi[i]['sc']:
                sc_node, ns_node = i, j
            else:
                sc_node, ns_node = j, i

            SC_cid = merge_contents[sc_node][0]
            C, Cc = merge_contents[ns_node][0], merge_contents[ns_node][1]

            w1 = W_u[SC_cid][C]
            w2 = W_u[SC_cid][Cc]
            total_w = w1 + w2
            ratio = w1 / w2 if w2 > 0 else float('inf')
            if ratio < 1: ratio = 1/ratio

            ratios.append(ratio)

            print(f"  ({i:>3},{j:>3}) {mi[i]['H']:>4} {mi[j]['H']:>4} "
                  f"{w1:>6} {w2:>8} {total_w:>6} {ratio:>7.2f}")

        if ratios:
            finite_ratios = [r for r in ratios if r < float('inf')]
            if finite_ratios:
                print(f"\n  Ratio stats: mean={np.mean(finite_ratios):.2f}, "
                      f"min={min(finite_ratios):.2f}, max={max(finite_ratios):.2f}")
                print(f"  Infinite ratios (one sub-flow = 0): "
                      f"{sum(1 for r in ratios if r == float('inf'))}/{len(ratios)}")

    # THE DEEP PATTERN: complement pairing within sub-flows
    # For NS-NS: does A→B always equal A^c→B^c? (complement preserves?)
    # And A→B^c always equal A^c→B? (complement swaps?)
    if ns_ns:
        print(f"\n  COMPLEMENT PATTERN IN SUB-FLOWS:")
        preserve_count = 0
        swap_count = 0
        neither_count = 0

        for (i, j) in ns_ns:
            cs_i = merge_contents[i]
            cs_j = merge_contents[j]
            if len(cs_i) < 2 or len(cs_j) < 2: continue

            A, Ac = cs_i[0], cs_i[1]
            B, Bc = cs_j[0], cs_j[1]

            # Check: A→B = A^c→B^c? (complement preserves the sub-edge)
            if W_u[A][B] == W_u[Ac][Bc] and W_u[A][Bc] == W_u[Ac][B]:
                preserve_count += 1
            else:
                neither_count += 1

        total_ns = preserve_count + neither_count
        print(f"  Complement-preserving: {preserve_count}/{total_ns} "
              f"({100*preserve_count/total_ns:.1f}%)")
        print(f"  Non-preserving: {neither_count}/{total_ns}")

        if preserve_count == total_ns:
            print(f"\n  *** COMPLEMENT ALWAYS PRESERVES SUB-FLOW STRUCTURE! ***")
            print(f"  For every NS-NS edge: w(A→B) = w(A^c→B^c) and w(A→B^c) = w(A^c→B)")
            print(f"  This means: the complement involution acts as a SYMMETRY on sub-flows.")
            print(f"  The merged W is symmetric BECAUSE of this complement symmetry.")

    return V

if __name__ == '__main__':
    print("=" * 72)
    print("  UNMERGED ASYMMETRY WITHIN SYMMETRIC W")
    print("  opus-2026-03-24-S281")
    print("=" * 72)

    for n in [4, 5, 6]:
        t0 = time.time()
        analyze_asymmetry(n)
        print(f"\n  Time: {time.time()-t0:.1f}s")

    print("\nDONE.")
