#!/usr/bin/env python3
"""
taxonomy_merged_s236.py -- opus-2026-03-23-S236

COMPREHENSIVE TAXONOMY OF G_n/Z_2 (MERGED META-GRAPH)

Create EVERY meaningful category for nodes and edges, then cross-tabulate.
Use the 5-labeling system from S20cz:

EDGE LABELINGS:
  1. SC-type: BLUE (SC↔SC) / BLACK (SC↔NS) / GREEN (NS↔NS)
  2. Score: AMBER (same score seq) / VIOLET (different score seq)
  3. Fiber: TEAL (same parent at n-1) / MAGENTA (different parent)
  4. dH: STEP (|dH|=2) / LEAP (|dH|>2) / LEVEL (dH=0)
  5. Self-loop: LOOP (C=D) / CROSS (C≠D)

NODE LABELINGS:
  1. SC-type: SC / NS
  2. H-zone: FLOOR (H min) / BULK / CEILING (H max)
  3. Aut: RIGID (|Aut|=1) / SYMMETRIC (|Aut|>1)
  4. Self-adjacent: has self-loop or not
  5. Degree: LOW / MEDIUM / HIGH
  6. Blue-degree: has any BLUE edge or not
  7. Black-degree: has any BLACK edge or not
  8. Green-degree: has any GREEN edge or not
  9. Amber-degree: has any AMBER edge or not
  10. Teal-degree: has any TEAL edge or not
  11. Parent count: how many distinct parents at n-1
  12. c3 parity: EVEN / ODD
  13. Level membership: how many nodes share the same H
  14. Score type: unique score or shared score

Author: opus-2026-03-23-S236
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
    """Check if A is isomorphic to its complement."""
    c_A = canon(A, n)
    c_C = canon(complement(A, n), n)
    return c_A == c_C

# ============================================================
# BUILD FOR n=5 AND n=6
# ============================================================

for n in [5, 6]:
    banner("COMPREHENSIVE TAXONOMY AT n=%d" % n)
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

    # Build iso classes
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

    # Compute per-class invariants
    class_info = {}
    for cid in range(nc):
        A = reps[cid]
        H = Hcount(A, n)
        c3 = c3count(A, n)
        sc = score_seq(A, n)
        is_self_comp = is_sc(A, n)

        # Find complement class
        comp_can = canon(complement(A, n), n)
        comp_cid = cmap[comp_can]

        # Automorphism count
        aut = factorial(n) // sizes[cid]

        class_info[cid] = {
            'H': H, 'c3': c3, 'score': sc, 'SC': is_self_comp,
            'comp': comp_cid, 'aut': aut, 'size': sizes[cid]
        }

    # Build merged graph: merge complement pairs
    merged = {}  # merged_id -> (cid1, cid2) or (cid,) for SC
    merged_id = {}  # cid -> merged_node_id
    mid = 0
    for cid in range(nc):
        if cid in merged_id: continue
        comp = class_info[cid]['comp']
        if comp == cid:  # SC
            merged[mid] = (cid,)
            merged_id[cid] = mid
        else:
            merged[mid] = (cid, comp)
            merged_id[cid] = mid
            merged_id[comp] = mid
        mid += 1
    n_merged = len(merged)

    # Build F matrix
    F = np.zeros((nc, nc), dtype=int)
    for bits in range(2**m):
        ci = tc[bits]
        for k in range(m):
            F[ci][tc[bits^(1<<k)]] += 1

    # Build merged adjacency
    merged_adj = defaultdict(int)  # (mi, mj) -> weight
    for ci in range(nc):
        mi = merged_id[ci]
        for cj in range(nc):
            if F[ci][cj] > 0:
                mj = merged_id[cj]
                if mi != mj:
                    key = (min(mi, mj), max(mi, mj))
                    merged_adj[key] += F[ci][cj]

    merged_edges = set()
    for (mi, mj) in merged_adj:
        merged_edges.add((mi, mj))
    n_merged_edges = len(merged_edges)

    # Compute parent at n-1 for each class
    if n >= 4:
        n1 = n - 1
        m1 = comb(n1, 2)
        pairs1 = [(i,j) for i in range(n1) for j in range(i+1,n1)]
        cmap1 = {}
        for bits in range(2**m1):
            A1 = [[0]*n1 for _ in range(n1)]
            for k,(ii,jj) in enumerate(pairs1):
                if (bits>>k)&1: A1[ii][jj]=1
                else: A1[jj][ii]=1
            c1 = canon(A1, n1)
            if c1 not in cmap1: cmap1[c1] = len(cmap1)

        # For each n-class, delete each vertex and find parent class
        parent_sets = {}
        for cid in range(nc):
            A = reps[cid]
            parents = set()
            for v in range(n):
                rem = [u for u in range(n) if u != v]
                A1 = [[A[rem[i]][rem[j]] for j in range(n-1)] for i in range(n-1)]
                parents.add(cmap1[canon(A1, n1)])
            parent_sets[cid] = parents

    # ============================================================
    # NODE CATEGORIES
    # ============================================================
    banner("NODE CATEGORIES (n=%d, V_merged=%d)" % (n, n_merged))

    # For each merged node, compute all categories
    node_cats = {}
    H_all = []
    for mid_val in range(n_merged):
        cids = merged[mid_val]
        cid = cids[0]
        info = class_info[cid]

        # SC type
        is_sc_node = info['SC']

        # H value
        H = info['H']
        H_all.append(H)

        # c3 parity
        c3_par = 'EVEN' if info['c3'] % 2 == 0 else 'ODD'

        # Aut
        aut = info['aut']
        rigid = aut == 1

        # Self-loop (does this class have F[C][C] > 0?)
        has_self_loop = F[cid][cid] > 0

        # Degree in merged graph
        deg = sum(1 for (a,b) in merged_edges if a == mid_val or b == mid_val)

        # Blue degree (connections to other SC nodes)
        blue_deg = 0; black_deg = 0; green_deg = 0
        for (a,b) in merged_edges:
            if a != mid_val and b != mid_val: continue
            other = b if a == mid_val else a
            other_sc = class_info[merged[other][0]]['SC']
            if is_sc_node and other_sc:
                blue_deg += 1
            elif is_sc_node != other_sc:
                black_deg += 1
            else:
                green_deg += 1

        # AMBER degree (connections preserving score)
        amber_deg = 0; teal_deg = 0
        for (a,b) in merged_edges:
            if a != mid_val and b != mid_val: continue
            other = b if a == mid_val else a
            other_cid = merged[other][0]
            if info['score'] == class_info[other_cid]['score']:
                amber_deg += 1
            if n >= 4 and parent_sets.get(cid) and parent_sets.get(other_cid):
                if parent_sets[cid] & parent_sets[other_cid]:
                    teal_deg += 1

        # Level membership (how many merged nodes share this H)
        H_counts = Counter(class_info[merged[m][0]]['H'] for m in range(n_merged))
        level_size = H_counts[H]

        # Parent count
        n_parents = len(parent_sets.get(cid, set())) if n >= 4 else 0

        node_cats[mid_val] = {
            'SC': is_sc_node,
            'H': H,
            'c3': info['c3'],
            'c3_par': c3_par,
            'score': info['score'],
            'aut': aut,
            'rigid': rigid,
            'has_self_loop': has_self_loop,
            'deg': deg,
            'blue_deg': blue_deg,
            'black_deg': black_deg,
            'green_deg': green_deg,
            'amber_deg': amber_deg,
            'teal_deg': teal_deg,
            'level_size': level_size,
            'n_parents': n_parents
        }

    H_min = min(H_all); H_max = max(H_all)

    # Assign H-zone
    for mid_val in range(n_merged):
        H = node_cats[mid_val]['H']
        if H == H_min:
            node_cats[mid_val]['H_zone'] = 'FLOOR'
        elif H == H_max:
            node_cats[mid_val]['H_zone'] = 'CEILING'
        else:
            node_cats[mid_val]['H_zone'] = 'BULK'

    # Print all node categories
    print("  %-3s | %-4s | %-3s | %-7s | %-6s | %-5s | %-4s | %-3s | %s | %s | %s | %s | %s | %s | %s | %s" %
          ('mid', 'SC', 'H', 'H_zone', 'c3_par', 'rigid', 'self', 'deg',
           'B', 'K', 'G', 'A', 'T', 'lvl', '#par', 'score'))
    print("  " + "-"*120)
    for mid_val in sorted(range(n_merged), key=lambda m: node_cats[m]['H']):
        nc_v = node_cats[mid_val]
        print("  %-3d | %-4s | %-3d | %-7s | %-6s | %-5s | %-4s | %-3d | %d | %d | %d | %d | %d | %d | %d | %s" %
              (mid_val, 'SC' if nc_v['SC'] else 'NS', nc_v['H'], nc_v['H_zone'],
               nc_v['c3_par'], 'RIGID' if nc_v['rigid'] else 'SYM',
               'LOOP' if nc_v['has_self_loop'] else '-',
               nc_v['deg'], nc_v['blue_deg'], nc_v['black_deg'], nc_v['green_deg'],
               nc_v['amber_deg'], nc_v['teal_deg'], nc_v['level_size'],
               nc_v['n_parents'], nc_v['score']))

    # ============================================================
    # NODE CATEGORY COUNTS
    # ============================================================
    banner("NODE CATEGORY COUNTS (n=%d)" % n)

    cats = [
        ('SC', lambda v: node_cats[v]['SC']),
        ('NS', lambda v: not node_cats[v]['SC']),
        ('RIGID (|Aut|=1)', lambda v: node_cats[v]['rigid']),
        ('SYMMETRIC (|Aut|>1)', lambda v: not node_cats[v]['rigid']),
        ('HAS_SELF_LOOP', lambda v: node_cats[v]['has_self_loop']),
        ('NO_SELF_LOOP', lambda v: not node_cats[v]['has_self_loop']),
        ('c3_EVEN', lambda v: node_cats[v]['c3_par'] == 'EVEN'),
        ('c3_ODD', lambda v: node_cats[v]['c3_par'] == 'ODD'),
        ('FLOOR (H=min)', lambda v: node_cats[v]['H_zone'] == 'FLOOR'),
        ('CEILING (H=max)', lambda v: node_cats[v]['H_zone'] == 'CEILING'),
        ('BULK', lambda v: node_cats[v]['H_zone'] == 'BULK'),
        ('HAS_BLUE_EDGE', lambda v: node_cats[v]['blue_deg'] > 0),
        ('NO_BLUE_EDGE', lambda v: node_cats[v]['blue_deg'] == 0),
        ('HAS_BLACK_EDGE', lambda v: node_cats[v]['black_deg'] > 0),
        ('NO_BLACK_EDGE', lambda v: node_cats[v]['black_deg'] == 0),
        ('HAS_GREEN_EDGE', lambda v: node_cats[v]['green_deg'] > 0),
        ('NO_GREEN_EDGE', lambda v: node_cats[v]['green_deg'] == 0),
        ('HAS_AMBER_EDGE', lambda v: node_cats[v]['amber_deg'] > 0),
        ('NO_AMBER_EDGE', lambda v: node_cats[v]['amber_deg'] == 0),
        ('HAS_TEAL_EDGE', lambda v: node_cats[v]['teal_deg'] > 0),
        ('NO_TEAL_EDGE', lambda v: node_cats[v]['teal_deg'] == 0),
        ('UNIQUE_SCORE', lambda v: Counter(node_cats[m]['score'] for m in range(n_merged))[node_cats[v]['score']] == 1),
        ('SHARED_SCORE', lambda v: Counter(node_cats[m]['score'] for m in range(n_merged))[node_cats[v]['score']] > 1),
        ('SINGLETON_LEVEL (unique H)', lambda v: node_cats[v]['level_size'] == 1),
        ('MULTI_LEVEL (shared H)', lambda v: node_cats[v]['level_size'] > 1),
        ('SINGLE_PARENT', lambda v: node_cats[v]['n_parents'] == 1),
        ('MULTI_PARENT', lambda v: node_cats[v]['n_parents'] > 1),
    ]

    for name, pred in cats:
        count = sum(1 for v in range(n_merged) if pred(v))
        print("  %-30s: %d / %d (%.1f%%)" % (name, count, n_merged, 100*count/n_merged))

    # ============================================================
    # EDGE CATEGORIES
    # ============================================================
    banner("EDGE CATEGORIES (n=%d, E_merged=%d)" % (n, n_merged_edges))

    edge_cats = {}
    for (mi, mj) in merged_edges:
        ni = node_cats[mi]; nj = node_cats[mj]
        ci = merged[mi][0]; cj = merged[mj][0]

        # SC-type
        if ni['SC'] and nj['SC']: sc_type = 'BLUE'
        elif ni['SC'] != nj['SC']: sc_type = 'BLACK'
        else: sc_type = 'GREEN'

        # Score change
        amber = ni['score'] == nj['score']

        # dH
        dH = abs(ni['H'] - nj['H'])
        if dH == 0: dh_type = 'LEVEL'
        elif dH == 2: dh_type = 'STEP'
        else: dh_type = 'LEAP'

        # Fiber (same parent)
        teal = bool(parent_sets.get(ci, set()) & parent_sets.get(cj, set())) if n >= 4 else False

        # c3 parity crossing
        c3_cross = ni['c3_par'] != nj['c3_par']

        # Aut change
        aut_change = ni['aut'] != nj['aut']

        # Degree category
        avg_deg = (ni['deg'] + nj['deg']) / 2

        # H direction (from lower to higher)
        if ni['H'] < nj['H']: h_dir = 'UPHILL'
        elif ni['H'] > nj['H']: h_dir = 'DOWNHILL'
        else: h_dir = 'LEVEL'

        # Thickness from F matrix
        thickness = 0
        for ci_orig in merged[mi]:
            for cj_orig in merged[mj]:
                thickness += F[ci_orig][cj_orig]

        edge_cats[(mi, mj)] = {
            'sc_type': sc_type,
            'amber': amber,
            'teal': teal,
            'dH': dH,
            'dh_type': dh_type,
            'h_dir': h_dir,
            'c3_cross': c3_cross,
            'aut_change': aut_change,
            'thickness': thickness
        }

    # Edge category counts
    edge_cat_list = [
        ('BLUE (SC↔SC)', lambda e: edge_cats[e]['sc_type'] == 'BLUE'),
        ('BLACK (SC↔NS)', lambda e: edge_cats[e]['sc_type'] == 'BLACK'),
        ('GREEN (NS↔NS)', lambda e: edge_cats[e]['sc_type'] == 'GREEN'),
        ('AMBER (same score)', lambda e: edge_cats[e]['amber']),
        ('VIOLET (diff score)', lambda e: not edge_cats[e]['amber']),
        ('TEAL (same parent)', lambda e: edge_cats[e]['teal']),
        ('MAGENTA (diff parent)', lambda e: not edge_cats[e]['teal']),
        ('STEP (|dH|=2)', lambda e: edge_cats[e]['dh_type'] == 'STEP'),
        ('LEAP (|dH|>2)', lambda e: edge_cats[e]['dh_type'] == 'LEAP'),
        ('LEVEL (dH=0)', lambda e: edge_cats[e]['dh_type'] == 'LEVEL'),
        ('c3_CROSSING (parity flip)', lambda e: edge_cats[e]['c3_cross']),
        ('c3_PRESERVING', lambda e: not edge_cats[e]['c3_cross']),
        ('AUT_CHANGE', lambda e: edge_cats[e]['aut_change']),
        ('AUT_PRESERVE', lambda e: not edge_cats[e]['aut_change']),
        ('THIN (thickness ≤ n!)', lambda e: edge_cats[e]['thickness'] <= factorial(n)),
        ('THICK (thickness > n!)', lambda e: edge_cats[e]['thickness'] > factorial(n)),
    ]

    for name, pred in edge_cat_list:
        count = sum(1 for e in merged_edges if pred(e))
        print("  %-35s: %d / %d (%.1f%%)" % (name, count, n_merged_edges, 100*count/n_merged_edges))

    # ============================================================
    # CROSS-TABULATIONS
    # ============================================================
    banner("CROSS-TABULATIONS (n=%d)" % n)

    # SC-type vs AMBER/VIOLET
    print("  SC-type vs Score:")
    for st in ['BLUE', 'BLACK', 'GREEN']:
        amber_c = sum(1 for e in merged_edges if edge_cats[e]['sc_type'] == st and edge_cats[e]['amber'])
        violet_c = sum(1 for e in merged_edges if edge_cats[e]['sc_type'] == st and not edge_cats[e]['amber'])
        total = amber_c + violet_c
        if total > 0:
            print("    %-6s: AMBER=%d, VIOLET=%d, total=%d, AMBER%%=%.1f%%" %
                  (st, amber_c, violet_c, total, 100*amber_c/total))

    # SC-type vs TEAL/MAGENTA
    print("\n  SC-type vs Fiber:")
    for st in ['BLUE', 'BLACK', 'GREEN']:
        teal_c = sum(1 for e in merged_edges if edge_cats[e]['sc_type'] == st and edge_cats[e]['teal'])
        mag_c = sum(1 for e in merged_edges if edge_cats[e]['sc_type'] == st and not edge_cats[e]['teal'])
        total = teal_c + mag_c
        if total > 0:
            print("    %-6s: TEAL=%d, MAGENTA=%d, total=%d, TEAL%%=%.1f%%" %
                  (st, teal_c, mag_c, total, 100*teal_c/total))

    # AMBER vs TEAL
    print("\n  Score vs Fiber:")
    at = sum(1 for e in merged_edges if edge_cats[e]['amber'] and edge_cats[e]['teal'])
    am = sum(1 for e in merged_edges if edge_cats[e]['amber'] and not edge_cats[e]['teal'])
    vt = sum(1 for e in merged_edges if not edge_cats[e]['amber'] and edge_cats[e]['teal'])
    vm = sum(1 for e in merged_edges if not edge_cats[e]['amber'] and not edge_cats[e]['teal'])
    print("    AMBER+TEAL=%d, AMBER+MAGENTA=%d, VIOLET+TEAL=%d, VIOLET+MAGENTA=%d" %
          (at, am, vt, vm))

    # SC-type vs dH
    print("\n  SC-type vs dH:")
    for st in ['BLUE', 'BLACK', 'GREEN']:
        for dht in ['STEP', 'LEAP', 'LEVEL']:
            c = sum(1 for e in merged_edges if edge_cats[e]['sc_type'] == st and edge_cats[e]['dh_type'] == dht)
            if c > 0:
                print("    %s + %s: %d" % (st, dht, c))

    # c3 crossing vs SC-type
    print("\n  c3 crossing vs SC-type:")
    for st in ['BLUE', 'BLACK', 'GREEN']:
        cross = sum(1 for e in merged_edges if edge_cats[e]['sc_type'] == st and edge_cats[e]['c3_cross'])
        pres = sum(1 for e in merged_edges if edge_cats[e]['sc_type'] == st and not edge_cats[e]['c3_cross'])
        total = cross + pres
        if total > 0:
            print("    %-6s: CROSSING=%d, PRESERVING=%d, cross%%=%.1f%%" %
                  (st, cross, pres, 100*cross/total if total > 0 else 0))

    # ============================================================
    # COMPOUND CATEGORIES (intersections)
    # ============================================================
    banner("COMPOUND NODE CATEGORIES (n=%d)" % n)

    compound_nodes = [
        ('SC + RIGID', lambda v: node_cats[v]['SC'] and node_cats[v]['rigid']),
        ('SC + SYMMETRIC', lambda v: node_cats[v]['SC'] and not node_cats[v]['rigid']),
        ('NS + RIGID', lambda v: not node_cats[v]['SC'] and node_cats[v]['rigid']),
        ('NS + SYMMETRIC', lambda v: not node_cats[v]['SC'] and not node_cats[v]['rigid']),
        ('SC + HAS_SELF_LOOP', lambda v: node_cats[v]['SC'] and node_cats[v]['has_self_loop']),
        ('NS + HAS_SELF_LOOP', lambda v: not node_cats[v]['SC'] and node_cats[v]['has_self_loop']),
        ('SC + NO_GREEN', lambda v: node_cats[v]['SC'] and node_cats[v]['green_deg'] == 0),
        ('NS + NO_BLUE', lambda v: not node_cats[v]['SC'] and node_cats[v]['blue_deg'] == 0),
        ('SC + HAS_AMBER', lambda v: node_cats[v]['SC'] and node_cats[v]['amber_deg'] > 0),
        ('NS + HAS_AMBER', lambda v: not node_cats[v]['SC'] and node_cats[v]['amber_deg'] > 0),
        ('SC + c3_EVEN', lambda v: node_cats[v]['SC'] and node_cats[v]['c3_par'] == 'EVEN'),
        ('SC + c3_ODD', lambda v: node_cats[v]['SC'] and node_cats[v]['c3_par'] == 'ODD'),
        ('CEILING + SC', lambda v: node_cats[v]['H_zone'] == 'CEILING' and node_cats[v]['SC']),
        ('FLOOR + SC', lambda v: node_cats[v]['H_zone'] == 'FLOOR' and node_cats[v]['SC']),
        ('SC + NO_BLACK (isolated SC)', lambda v: node_cats[v]['SC'] and node_cats[v]['black_deg'] == 0),
        ('SINGLETON_H + SC', lambda v: node_cats[v]['level_size'] == 1 and node_cats[v]['SC']),
        ('SINGLE_PARENT + SC', lambda v: node_cats[v]['n_parents'] == 1 and node_cats[v]['SC']),
        ('HAS_TEAL + HAS_AMBER', lambda v: node_cats[v]['teal_deg'] > 0 and node_cats[v]['amber_deg'] > 0),
        ('NO_TEAL + NO_AMBER', lambda v: node_cats[v]['teal_deg'] == 0 and node_cats[v]['amber_deg'] == 0),
        ('MAX_DEGREE', lambda v: node_cats[v]['deg'] == max(node_cats[m]['deg'] for m in range(n_merged))),
        ('MIN_DEGREE', lambda v: node_cats[v]['deg'] == min(node_cats[m]['deg'] for m in range(n_merged))),
    ]

    for name, pred in compound_nodes:
        count = sum(1 for v in range(n_merged) if pred(v))
        if count > 0:
            print("  %-35s: %d / %d (%.1f%%)" % (name, count, n_merged, 100*count/n_merged))

    # ============================================================
    # COMPOUND EDGE CATEGORIES
    # ============================================================
    banner("COMPOUND EDGE CATEGORIES (n=%d)" % n)

    compound_edges = [
        ('BLUE + AMBER (SC↔SC, same score)', lambda e: edge_cats[e]['sc_type'] == 'BLUE' and edge_cats[e]['amber']),
        ('BLUE + VIOLET (SC↔SC, diff score)', lambda e: edge_cats[e]['sc_type'] == 'BLUE' and not edge_cats[e]['amber']),
        ('BLUE + TEAL (SC↔SC, same parent)', lambda e: edge_cats[e]['sc_type'] == 'BLUE' and edge_cats[e]['teal']),
        ('BLUE + MAGENTA (SC↔SC, diff parent)', lambda e: edge_cats[e]['sc_type'] == 'BLUE' and not edge_cats[e]['teal']),
        ('BLACK + AMBER', lambda e: edge_cats[e]['sc_type'] == 'BLACK' and edge_cats[e]['amber']),
        ('BLACK + TEAL', lambda e: edge_cats[e]['sc_type'] == 'BLACK' and edge_cats[e]['teal']),
        ('GREEN + AMBER', lambda e: edge_cats[e]['sc_type'] == 'GREEN' and edge_cats[e]['amber']),
        ('GREEN + TEAL', lambda e: edge_cats[e]['sc_type'] == 'GREEN' and edge_cats[e]['teal']),
        ('AMBER + TEAL (same score + same parent)', lambda e: edge_cats[e]['amber'] and edge_cats[e]['teal']),
        ('AMBER + LEVEL (same score + same H)', lambda e: edge_cats[e]['amber'] and edge_cats[e]['dh_type'] == 'LEVEL'),
        ('LEVEL + c3_CROSSING', lambda e: edge_cats[e]['dh_type'] == 'LEVEL' and edge_cats[e]['c3_cross']),
        ('STEP + c3_PRESERVING', lambda e: edge_cats[e]['dh_type'] == 'STEP' and not edge_cats[e]['c3_cross']),
        ('THIN + BLUE', lambda e: edge_cats[e]['thickness'] <= factorial(n) and edge_cats[e]['sc_type'] == 'BLUE'),
        ('THICK + GREEN', lambda e: edge_cats[e]['thickness'] > factorial(n) and edge_cats[e]['sc_type'] == 'GREEN'),
        ('LEAP + c3_CROSSING', lambda e: edge_cats[e]['dh_type'] == 'LEAP' and edge_cats[e]['c3_cross']),
        ('LEAP + AMBER', lambda e: edge_cats[e]['dh_type'] == 'LEAP' and edge_cats[e]['amber']),
    ]

    for name, pred in compound_edges:
        count = sum(1 for e in merged_edges if pred(e))
        if count > 0:
            print("  %-45s: %d / %d (%.1f%%)" % (name, count, n_merged_edges, 100*count/n_merged_edges))

    # ============================================================
    # BOOLEAN NODE PROFILES
    # ============================================================
    banner("BOOLEAN NODE PROFILES (n=%d)" % n)

    # Each node gets a binary profile: (SC, RIGID, SELF_LOOP, HAS_BLUE, HAS_BLACK, HAS_GREEN, HAS_AMBER, HAS_TEAL, c3_EVEN)
    profiles = Counter()
    for v in range(n_merged):
        nc_v = node_cats[v]
        profile = (
            nc_v['SC'],
            nc_v['rigid'],
            nc_v['has_self_loop'],
            nc_v['blue_deg'] > 0,
            nc_v['black_deg'] > 0,
            nc_v['green_deg'] > 0,
            nc_v['amber_deg'] > 0,
            nc_v['teal_deg'] > 0,
            nc_v['c3_par'] == 'EVEN'
        )
        profiles[profile] += 1

    labels = ['SC', 'RIGID', 'LOOP', 'BLUE', 'BLACK', 'GREEN', 'AMBER', 'TEAL', 'c3_EVEN']
    print("  %d distinct boolean profiles out of %d nodes:\n" % (len(profiles), n_merged))
    for profile, count in sorted(profiles.items(), key=lambda x: -x[1]):
        flags = [labels[i] for i in range(9) if profile[i]]
        anti_flags = ['!' + labels[i] for i in range(9) if not profile[i]]
        desc = ', '.join(flags) if flags else '(none set)'
        print("    %d nodes: %s" % (count, desc))

    # ============================================================
    # BOOLEAN EDGE PROFILES
    # ============================================================
    banner("BOOLEAN EDGE PROFILES (n=%d)" % n)

    edge_profiles = Counter()
    for e in merged_edges:
        ec = edge_cats[e]
        profile = (
            ec['sc_type'],
            ec['amber'],
            ec['teal'],
            ec['dh_type'],
            ec['c3_cross'],
            ec['aut_change'],
            ec['thickness'] > factorial(n)
        )
        edge_profiles[profile] += 1

    print("  %d distinct edge profiles out of %d edges:\n" % (len(edge_profiles), n_merged_edges))
    for profile, count in sorted(edge_profiles.items(), key=lambda x: -x[1])[:20]:
        sc_t, amber, teal, dh_t, c3x, autx, thick = profile
        flags = []
        flags.append(sc_t)
        flags.append('AMBER' if amber else 'VIOLET')
        flags.append('TEAL' if teal else 'MAGENTA')
        flags.append(dh_t)
        if c3x: flags.append('c3_CROSS')
        if autx: flags.append('AUT_CHANGE')
        if thick: flags.append('THICK')
        print("    %3d edges: %s" % (count, ' + '.join(flags)))

print("\nDone. opus-2026-03-23-S236")
