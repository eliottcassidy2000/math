#!/usr/bin/env python3
"""
deep_taxonomy_s236b.py -- opus-2026-03-23-S236b

DEEPER TAXONOMY: neighborhood-based node categories and more edge types.

New node categories:
  - HALO (nodes reachable from SC via BLACK edges)
  - ISOLATED_NS (NS nodes with no BLACK edges)
  - BRIDGE (nodes connecting SC and NS components)
  - AMBER_HUB (high AMBER degree)
  - TEAL_BOTTLENECK (all edges are TEAL = forced by fiber structure)
  - H_GRADIENT_SIGN (more uphill or downhill neighbors)
  - DEGREE_QUARTILE (Q1..Q4)

New edge categories:
  - THICKNESS_BIN (thin/medium/thick/very_thick)
  - H_DIRECTION (uphill/downhill/level)
  - BRIDGE_EDGE (connects different connectivity components)
  - AMBER+BLUE ("golden spine" = same score + SC backbone)
  - TEAL+LEVEL ("fiber degeneracy" = same parent + same H)

Author: opus-2026-03-23-S236b
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def canon(A, n):
    sc = [sum(A[i]) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(sg.keys())]
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f < best: best = f
    return best

def Hcount(A, n):
    dp = {}
    for v in range(n): dp[(1<<v, v)] = 1
    full = (1<<n)-1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S & (1<<v)): continue
            c = dp.get((S, v), 0)
            if c == 0: continue
            for w in range(n):
                if S & (1<<w): continue
                if A[v][w]:
                    dp[(S|(1<<w), w)] = dp.get((S|(1<<w), w), 0) + c
    return sum(dp.get((full, v), 0) for v in range(n))

def c3count(A, n):
    c = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]: c += 1
                if A[i][k] and A[k][j] and A[j][i]: c += 1
    return c

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def is_sc(A, n):
    return canon(A, n) == canon(complement(A, n), n)

for n in [5, 6]:
    banner("DEEP TAXONOMY AT n=%d" % n)
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

    cmap = {}; sizes = {}; reps = {}; tc = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canon(A, n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid; sizes[cid] = 0
            reps[cid] = [row[:] for row in A]
        sizes[cmap[c]] += 1
        tc[bits] = cmap[c]
    nc = len(cmap)

    class_info = {}
    for cid in range(nc):
        A = reps[cid]
        comp_can = canon(complement(A, n), n)
        class_info[cid] = {
            'H': Hcount(A, n), 'c3': c3count(A, n), 'score': score_seq(A, n),
            'SC': is_sc(A, n), 'comp': cmap[comp_can], 'aut': factorial(n) // sizes[cid]
        }

    # Merged graph
    merged = {}; merged_id = {}; mid = 0
    for cid in range(nc):
        if cid in merged_id: continue
        comp = class_info[cid]['comp']
        if comp == cid: merged[mid] = (cid,)
        else: merged[mid] = (cid, comp)
        merged_id[cid] = mid
        if comp != cid: merged_id[comp] = mid
        mid += 1
    n_merged = mid

    # F matrix
    F = np.zeros((nc, nc), dtype=int)
    for bits in range(2**m):
        ci = tc[bits]
        for k in range(m):
            F[ci][tc[bits^(1<<k)]] += 1

    # Merged edges with thickness
    merged_adj = {}
    for ci in range(nc):
        mi = merged_id[ci]
        for cj in range(nc):
            if F[ci][cj] > 0:
                mj = merged_id[cj]
                if mi != mj:
                    key = (min(mi, mj), max(mi, mj))
                    merged_adj[key] = merged_adj.get(key, 0) + F[ci][cj]

    merged_edges = set(merged_adj.keys())
    n_merged_edges = len(merged_edges)

    # Parent computation
    if n >= 4:
        n1 = n - 1; m1 = comb(n1, 2)
        pairs1 = [(i,j) for i in range(n1) for j in range(i+1,n1)]
        cmap1 = {}
        for bits in range(2**m1):
            A1 = [[0]*n1 for _ in range(n1)]
            for k,(ii,jj) in enumerate(pairs1):
                if (bits>>k)&1: A1[ii][jj]=1
                else: A1[jj][ii]=1
            c1 = canon(A1, n1)
            if c1 not in cmap1: cmap1[c1] = len(cmap1)
        parent_sets = {}
        for cid in range(nc):
            A = reps[cid]
            parents = set()
            for v in range(n):
                rem = [u for u in range(n) if u != v]
                A1 = [[A[rem[i]][rem[j]] for j in range(n-1)] for i in range(n-1)]
                parents.add(cmap1[canon(A1, n1)])
            parent_sets[cid] = parents

    # Basic node info for merged nodes
    node_info = {}
    for mi in range(n_merged):
        cid = merged[mi][0]
        info = class_info[cid]
        node_info[mi] = {
            'SC': info['SC'], 'H': info['H'], 'c3': info['c3'],
            'score': info['score'], 'aut': info['aut'],
            'cid': cid, 'cids': merged[mi]
        }

    # Adjacency as neighbor lists
    neighbors = defaultdict(set)
    for (a,b) in merged_edges:
        neighbors[a].add(b); neighbors[b].add(a)

    # ============================================================
    # DEEPER NODE CATEGORIES
    # ============================================================

    for mi in range(n_merged):
        ni = node_info[mi]
        nbrs = neighbors[mi]

        # Degree
        ni['deg'] = len(nbrs)

        # Self-loop
        ci = ni['cid']
        ni['self_loop'] = F[ci][ci] > 0

        # Neighbor types
        sc_nbrs = [mj for mj in nbrs if node_info[mj]['SC']]
        ns_nbrs = [mj for mj in nbrs if not node_info[mj]['SC']]
        ni['blue_deg'] = len(sc_nbrs) if ni['SC'] else 0
        ni['black_deg'] = len(ns_nbrs) if ni['SC'] else len(sc_nbrs)
        ni['green_deg'] = len(ns_nbrs) if not ni['SC'] else 0

        # H-gradient: how many neighbors have higher vs lower H
        H = ni['H']
        uphill = sum(1 for mj in nbrs if node_info[mj]['H'] > H)
        downhill = sum(1 for mj in nbrs if node_info[mj]['H'] < H)
        level = sum(1 for mj in nbrs if node_info[mj]['H'] == H)
        ni['uphill'] = uphill; ni['downhill'] = downhill; ni['level_nbrs'] = level

        # H-gradient sign: net direction
        if uphill > downhill: ni['gradient'] = 'SOURCE'  # more uphill = node is a source
        elif downhill > uphill: ni['gradient'] = 'SINK'   # more downhill = node is a sink
        else: ni['gradient'] = 'BALANCED'

        # AMBER degree (score-preserving neighbors)
        amber_nbrs = [mj for mj in nbrs if node_info[mj]['score'] == ni['score']]
        ni['amber_deg'] = len(amber_nbrs)

        # TEAL degree (same-parent neighbors)
        cid_i = ni['cid']
        teal_nbrs = [mj for mj in nbrs if parent_sets.get(cid_i, set()) & parent_sets.get(node_info[mj]['cid'], set())]
        ni['teal_deg'] = len(teal_nbrs)

        # HALO: NS node reachable from SC via ≤1 BLACK edge
        if not ni['SC'] and len(sc_nbrs) > 0:
            ni['halo'] = 'INNER_HALO'  # directly adjacent to SC
        elif not ni['SC'] and len(sc_nbrs) == 0:
            ni['halo'] = 'OUTER_SEA'   # not adjacent to SC at all
        elif ni['SC']:
            ni['halo'] = 'SPINE'

        # BRIDGE: SC node with both BLUE and BLACK edges
        if ni['SC'] and ni['blue_deg'] > 0 and ni['black_deg'] > 0:
            ni['bridge'] = True
        else:
            ni['bridge'] = False

        # Pure types
        if ni['SC'] and ni['black_deg'] == 0:
            ni['pure'] = 'ISOLATED_SC'  # SC with no NS neighbors
        elif ni['SC']:
            ni['pure'] = 'CONNECTED_SC'
        elif not ni['SC'] and ni['black_deg'] == 0:
            ni['pure'] = 'ISOLATED_NS'  # NS with no SC neighbors
        else:
            ni['pure'] = 'CONNECTED_NS'

    # Degree quartiles
    all_degs = sorted(node_info[mi]['deg'] for mi in range(n_merged))
    q1 = all_degs[n_merged // 4]
    q2 = all_degs[n_merged // 2]
    q3 = all_degs[3 * n_merged // 4]
    for mi in range(n_merged):
        d = node_info[mi]['deg']
        if d <= q1: node_info[mi]['deg_q'] = 'Q1_LOW'
        elif d <= q2: node_info[mi]['deg_q'] = 'Q2_MED'
        elif d <= q3: node_info[mi]['deg_q'] = 'Q3_HIGH'
        else: node_info[mi]['deg_q'] = 'Q4_TOP'

    # ============================================================
    # PRINT DEEP NODE CATEGORIES
    # ============================================================

    print("  DEEP NODE CATEGORIES (n=%d):\n" % n)

    deep_cats = [
        ('SPINE', lambda v: node_info[v]['halo'] == 'SPINE'),
        ('INNER_HALO', lambda v: node_info[v]['halo'] == 'INNER_HALO'),
        ('OUTER_SEA', lambda v: node_info[v]['halo'] == 'OUTER_SEA'),
        ('ISOLATED_SC', lambda v: node_info[v]['pure'] == 'ISOLATED_SC'),
        ('CONNECTED_SC', lambda v: node_info[v]['pure'] == 'CONNECTED_SC'),
        ('ISOLATED_NS', lambda v: node_info[v]['pure'] == 'ISOLATED_NS'),
        ('CONNECTED_NS', lambda v: node_info[v]['pure'] == 'CONNECTED_NS'),
        ('BRIDGE (SC w/ blue+black)', lambda v: node_info[v]['bridge']),
        ('GRADIENT=SOURCE (more uphill)', lambda v: node_info[v]['gradient'] == 'SOURCE'),
        ('GRADIENT=SINK (more downhill)', lambda v: node_info[v]['gradient'] == 'SINK'),
        ('GRADIENT=BALANCED', lambda v: node_info[v]['gradient'] == 'BALANCED'),
        ('HAS_LEVEL_NEIGHBOR', lambda v: node_info[v]['level_nbrs'] > 0),
        ('AMBER_HUB (amber_deg ≥ 3)', lambda v: node_info[v]['amber_deg'] >= 3),
        ('AMBER_ISOLATED (amber_deg = 0)', lambda v: node_info[v]['amber_deg'] == 0),
        ('SELF_LOOP', lambda v: node_info[v]['self_loop']),
        ('Q1_LOW degree', lambda v: node_info[v]['deg_q'] == 'Q1_LOW'),
        ('Q2_MED degree', lambda v: node_info[v]['deg_q'] == 'Q2_MED'),
        ('Q3_HIGH degree', lambda v: node_info[v]['deg_q'] == 'Q3_HIGH'),
        ('Q4_TOP degree', lambda v: node_info[v]['deg_q'] == 'Q4_TOP'),
        ('SC + SELF_LOOP', lambda v: node_info[v]['SC'] and node_info[v]['self_loop']),
        ('NS + OUTER_SEA', lambda v: node_info[v]['halo'] == 'OUTER_SEA'),
        ('SC + SOURCE gradient', lambda v: node_info[v]['SC'] and node_info[v]['gradient'] == 'SOURCE'),
        ('SC + SINK gradient', lambda v: node_info[v]['SC'] and node_info[v]['gradient'] == 'SINK'),
        ('NS + SOURCE gradient', lambda v: not node_info[v]['SC'] and node_info[v]['gradient'] == 'SOURCE'),
        ('NS + SINK gradient', lambda v: not node_info[v]['SC'] and node_info[v]['gradient'] == 'SINK'),
        ('AMBER_HUB + SC', lambda v: node_info[v]['amber_deg'] >= 3 and node_info[v]['SC']),
        ('AMBER_HUB + NS', lambda v: node_info[v]['amber_deg'] >= 3 and not node_info[v]['SC']),
    ]

    for name, pred in deep_cats:
        count = sum(1 for v in range(n_merged) if pred(v))
        if count > 0:
            pct = 100 * count / n_merged
            print("  %-40s: %2d / %d (%5.1f%%)" % (name, count, n_merged, pct))

    # ============================================================
    # DEEPER EDGE CATEGORIES
    # ============================================================

    banner("DEEPER EDGE CATEGORIES (n=%d)" % n)

    # Compute edge info
    edge_info = {}
    for (mi, mj) in merged_edges:
        ni = node_info[mi]; nj = node_info[mj]

        # SC-type
        if ni['SC'] and nj['SC']: sc_type = 'BLUE'
        elif ni['SC'] != nj['SC']: sc_type = 'BLACK'
        else: sc_type = 'GREEN'

        dH = abs(ni['H'] - nj['H'])
        amber = ni['score'] == nj['score']
        teal = bool(parent_sets.get(ni['cid'], set()) & parent_sets.get(nj['cid'], set()))
        c3_cross = ni['c3'] % 2 != nj['c3'] % 2

        # H direction
        if ni['H'] < nj['H']: h_dir = 'UPHILL'
        elif ni['H'] > nj['H']: h_dir = 'DOWNHILL'
        else: h_dir = 'LEVEL'

        # Thickness
        thickness = merged_adj[(mi, mj)]
        if thickness <= factorial(n) // 2: t_bin = 'WHISPER'
        elif thickness <= factorial(n): t_bin = 'THIN'
        elif thickness <= 4 * factorial(n): t_bin = 'MEDIUM'
        else: t_bin = 'HEAVY'

        # Halo crossing
        halo_i = ni['halo']; halo_j = nj['halo']
        if halo_i == halo_j: halo_edge = 'SAME_ZONE'
        else: halo_edge = 'CROSS_ZONE'

        # Gradient alignment
        if h_dir == 'UPHILL':
            # mi has lower H. Is mi a SOURCE? nj a SINK?
            aligned = ni['gradient'] == 'SOURCE' and nj['gradient'] == 'SINK'
        elif h_dir == 'DOWNHILL':
            aligned = ni['gradient'] == 'SINK' and nj['gradient'] == 'SOURCE'
        else:
            aligned = False

        edge_info[(mi, mj)] = {
            'sc_type': sc_type, 'dH': dH, 'amber': amber, 'teal': teal,
            'c3_cross': c3_cross, 'h_dir': h_dir, 't_bin': t_bin,
            'halo_edge': halo_edge, 'aligned': aligned, 'thickness': thickness
        }

    deep_edge_cats = [
        ('WHISPER (very thin)', lambda e: edge_info[e]['t_bin'] == 'WHISPER'),
        ('THIN', lambda e: edge_info[e]['t_bin'] == 'THIN'),
        ('MEDIUM', lambda e: edge_info[e]['t_bin'] == 'MEDIUM'),
        ('HEAVY', lambda e: edge_info[e]['t_bin'] == 'HEAVY'),
        ('UPHILL', lambda e: edge_info[e]['h_dir'] == 'UPHILL'),
        ('DOWNHILL', lambda e: edge_info[e]['h_dir'] == 'DOWNHILL'),
        ('LEVEL', lambda e: edge_info[e]['h_dir'] == 'LEVEL'),
        ('SAME_ZONE', lambda e: edge_info[e]['halo_edge'] == 'SAME_ZONE'),
        ('CROSS_ZONE', lambda e: edge_info[e]['halo_edge'] == 'CROSS_ZONE'),
        ('GRADIENT_ALIGNED', lambda e: edge_info[e]['aligned']),
        ('GOLDEN SPINE (BLUE+AMBER)', lambda e: edge_info[e]['sc_type'] == 'BLUE' and edge_info[e]['amber']),
        ('FIBER DEGENERATE (TEAL+LEVEL)', lambda e: edge_info[e]['teal'] and edge_info[e]['dH'] == 0),
        ('PURE LEAP (dH ≥ 8)', lambda e: edge_info[e]['dH'] >= 8),
        ('MICRO STEP (dH = 2)', lambda e: edge_info[e]['dH'] == 2),
        ('GREEN + HEAVY (NS sea bulk)', lambda e: edge_info[e]['sc_type'] == 'GREEN' and edge_info[e]['t_bin'] == 'HEAVY'),
        ('BLACK + THIN (delicate rib)', lambda e: edge_info[e]['sc_type'] == 'BLACK' and edge_info[e]['t_bin'] in ('THIN', 'WHISPER')),
        ('AMBER + c3_CROSS', lambda e: edge_info[e]['amber'] and edge_info[e]['c3_cross']),
        ('AMBER + c3_PRESERVE', lambda e: edge_info[e]['amber'] and not edge_info[e]['c3_cross']),
        ('BLUE + WHISPER', lambda e: edge_info[e]['sc_type'] == 'BLUE' and edge_info[e]['t_bin'] == 'WHISPER'),
    ]

    for name, pred in deep_edge_cats:
        count = sum(1 for e in merged_edges if pred(e))
        if count > 0:
            pct = 100 * count / n_merged_edges
            print("  %-40s: %2d / %d (%5.1f%%)" % (name, count, n_merged_edges, pct))

    # ============================================================
    # THE GRAND TAXONOMY TABLE
    # ============================================================
    banner("GRAND TAXONOMY SUMMARY (n=%d)" % n)

    # Count total distinct categories defined
    n_node_cats = len(deep_cats) + 27  # from S236 + new
    n_edge_cats = len(deep_edge_cats) + 16  # from S236 + new
    print("  TOTAL NODE CATEGORIES DEFINED: ~%d" % n_node_cats)
    print("  TOTAL EDGE CATEGORIES DEFINED: ~%d" % n_edge_cats)
    print("  V_merged = %d, E_merged = %d" % (n_merged, n_merged_edges))

    # Most informative single category (highest entropy split)
    print("\n  MOST BALANCED SPLITS (closest to 50/50):")
    all_cats = deep_cats + [
        ('SC', lambda v: node_info[v]['SC']),
        ('RIGID', lambda v: node_info[v]['aut'] == 1),
        ('SELF_LOOP', lambda v: node_info[v]['self_loop']),
        ('c3_EVEN', lambda v: node_info[v]['c3'] % 2 == 0),
        ('HAS_AMBER', lambda v: node_info[v]['amber_deg'] > 0),
        ('HAS_BLUE', lambda v: node_info[v]['blue_deg'] > 0),
    ]

    splits = []
    for name, pred in all_cats:
        count = sum(1 for v in range(n_merged) if pred(v))
        ratio = count / n_merged if n_merged > 0 else 0
        balance = 1 - 2 * abs(ratio - 0.5)  # 1 = perfect 50/50, 0 = all-or-nothing
        splits.append((balance, ratio, count, name))

    for balance, ratio, count, name in sorted(splits, reverse=True)[:10]:
        print("    %-40s: %d/%d = %.1f%% (balance=%.2f)" %
              (name, count, n_merged, 100*ratio, balance))

print("\nDone. opus-2026-03-23-S236b")
