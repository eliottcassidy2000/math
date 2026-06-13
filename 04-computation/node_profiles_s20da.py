#!/usr/bin/env python3
"""
node_profiles_s20da.py — Complete node profiling of G_n/Z_2
kind-pasteur-2026-03-23-S20da (overnight)

Every merged class gets a COMPLETE PROFILE of categorical and numerical features.
Then: cluster, correlate, find which features predict which, discover formulas.

NODE FEATURES:
  Intrinsic: H, |Aut|, score_seq, c3, SC/NS, fiber_size (= H/|Aut|)
  Tiling: #grid-sym tilings, grid-sym fraction, #self-loop tiles
  Degree: total_deg, spine_deg (SC-SC), rib_deg (SC-NS), sea_deg (NS-NS)
  Explorer: blue_deg, black_deg, mixed_deg
  H-structure: H_rank, max_neighbor_H, min_neighbor_H, avg_neighbor_H
  Score: score_variance, score_is_palindromic, score_is_regular
  New colors: amber_deg (score-preserving neighbors), teal_deg (same-fiber)
  Parent: parent_class at n-1, #children at n+1 (if available)
  Position: distance_from_transitive, distance_from_regular
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import numpy as np
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  COMPLETE NODE PROFILING OF G_n/Z_2")
print("  kind-pasteur-2026-03-23-S20da")
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

def c3count(a, n):
    c = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if a[i][j] and a[j][k] and a[k][i]: c += 1
                if a[i][k] and a[k][j] and a[j][i]: c += 1
    return c

# ============================================================================
# BUILD FULL PROFILE
# ============================================================================

for n in [5, 6]:
    print(f"\n{'#'*80}")
    print(f"  n = {n}")
    print(f"{'#'*80}")
    t0 = time.time()

    tiles = [(r,c) for r in range(1,n-1) for c in range(1,n-r)]
    m = len(tiles)
    tile_arcs = [(r+c, c-1) for r,c in tiles]

    # Explorer grid symmetry
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

    def build_tournament(tb):
        a = [[0]*n for _ in range(n)]
        for j in range(1, n): a[j][j-1] = 1
        for k in range(m):
            i, j = tile_arcs[k]
            if tb & (1 << k): a[i][j] = 1
            else: a[j][i] = 1
        return a

    # Enumerate
    tiling_data = {}
    class_info = {}
    for tb in range(1 << m):
        T = build_tournament(tb)
        cn = canon(T, n)
        gs = is_grid_sym(tb)
        if cn not in class_info:
            Tc = [[1-T[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
            sc_seq = tuple(sorted(sum(T[i][j] for j in range(n)) for i in range(n)))
            T_del = [[T[i][j] for j in range(n-1)] for i in range(n-1)]
            aut = sum(1 for p in permutations(range(n))
                      if all(T[p[i]][p[j]] == T[i][j] for i in range(n) for j in range(n)))
            class_info[cn] = {
                'H': Hdp(T, n), 'sc': canon(Tc, n) == cn,
                'score': sc_seq, 'c3': c3count(T, n), 'aut': aut,
                'parent': canon(T_del, n-1),
                'tilings': [], 'grid_sym_tilings': 0, 'self_loop_count': 0
            }
        class_info[cn]['tilings'].append(tb)
        if gs: class_info[cn]['grid_sym_tilings'] += 1
        tiling_data[tb] = {'canon': cn, 'grid_sym': gs}

    # Merge
    merge = {}; mid = 0
    for cn in sorted(class_info.keys()):
        if cn in merge: continue
        ci = class_info[cn]
        merge[cn] = mid
        if not ci['sc']:
            T0 = build_tournament(ci['tilings'][0])
            Tc = [[1-T0[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
            cc = canon(Tc, n)
            if cc != cn and cc in class_info: merge[cc] = mid
        mid += 1
    V = mid

    # Merged node info
    nodes = {}
    for cn, ci in class_info.items():
        mi = merge[cn]
        if mi not in nodes:
            fiber = len(ci['tilings'])
            gs_frac = ci['grid_sym_tilings'] / fiber if fiber > 0 else 0
            sc_var = np.var(list(ci['score']))
            is_palin = all(ci['score'][i] + ci['score'][n-1-i] == n-1 for i in range(n//2))
            is_reg = len(set(ci['score'])) == 1
            nodes[mi] = {
                'H': ci['H'], 'sc': ci['sc'], 'score': ci['score'],
                'c3': ci['c3'], 'aut': ci['aut'], 'fiber': fiber,
                'gs_count': ci['grid_sym_tilings'], 'gs_frac': gs_frac,
                'score_var': sc_var, 'palindromic': is_palin, 'regular': is_reg,
                'parent': ci['parent'],
                # Degree counts (to be filled)
                'deg': 0, 'spine_deg': 0, 'rib_deg': 0, 'sea_deg': 0,
                'blue_deg': 0, 'black_deg': 0, 'mixed_deg': 0,
                'amber_deg': 0, 'violet_deg': 0,
                'teal_deg': 0, 'magenta_deg': 0,
                'step_deg': 0, 'leap_deg': 0, 'level_deg': 0,
                'self_loops': 0,
                'max_nbr_H': 0, 'min_nbr_H': 999, 'sum_nbr_H': 0,
            }
        else:
            # Accumulate grid-sym from complement class
            nodes[mi]['gs_count'] += ci['grid_sym_tilings']
            nodes[mi]['fiber'] += len(ci['tilings'])

    # Build edges and count per-node degrees
    edges = {}
    for tb in range(1 << m):
        td = tiling_data[tb]
        src = merge[td['canon']]
        for k in range(m):
            tb2 = tb ^ (1 << k)
            td2 = tiling_data[tb2]
            tgt = merge[td2['canon']]
            if tgt == src:
                nodes[src]['self_loops'] += 1
                continue
            e = (min(src, tgt), max(src, tgt))
            if e not in edges:
                a, b = e
                na, nb = nodes[a], nodes[b]
                dH = abs(na['H'] - nb['H'])
                if na['sc'] and nb['sc']: sc_type = 'SC-SC'
                elif not na['sc'] and not nb['sc']: sc_type = 'NS-NS'
                else: sc_type = 'SC-NS'
                edges[e] = {
                    'sc_type': sc_type, 'dH': dH,
                    'score_same': na['score'] == nb['score'],
                    'same_parent': na['parent'] == nb['parent'],
                    'blue_ct': 0, 'black_ct': 0,
                }
            both_sym = td['grid_sym'] and td2['grid_sym']
            if both_sym: edges[e]['blue_ct'] += 1
            else: edges[e]['black_ct'] += 1

    # Classify edges and accumulate node degrees
    for e, info in edges.items():
        a, b = e
        if info['blue_ct'] > 0 and info['black_ct'] == 0: color = 'blue'
        elif info['blue_ct'] == 0: color = 'black'
        else: color = 'mixed'

        for v in [a, b]:
            nodes[v]['deg'] += 1
            if info['sc_type'] == 'SC-SC': nodes[v]['spine_deg'] += 1
            elif info['sc_type'] == 'SC-NS': nodes[v]['rib_deg'] += 1
            else: nodes[v]['sea_deg'] += 1
            if color == 'blue': nodes[v]['blue_deg'] += 1
            elif color == 'black': nodes[v]['black_deg'] += 1
            else: nodes[v]['mixed_deg'] += 1
            if info['score_same']: nodes[v]['amber_deg'] += 1
            else: nodes[v]['violet_deg'] += 1
            if info['same_parent']: nodes[v]['teal_deg'] += 1
            else: nodes[v]['magenta_deg'] += 1
            if info['dH'] == 0: nodes[v]['level_deg'] += 1
            elif info['dH'] == 2: nodes[v]['step_deg'] += 1
            else: nodes[v]['leap_deg'] += 1

        other_H = [nodes[b]['H'], nodes[a]['H']]
        for i, v in enumerate([a, b]):
            oh = other_H[i]
            nodes[v]['max_nbr_H'] = max(nodes[v]['max_nbr_H'], oh)
            nodes[v]['min_nbr_H'] = min(nodes[v]['min_nbr_H'], oh)
            nodes[v]['sum_nbr_H'] += oh

    # Distance from transitive and regular
    from collections import deque
    adj = defaultdict(set)
    for (a,b) in edges: adj[a].add(b); adj[b].add(a)

    trans_mid = [mi for mi, nd in nodes.items() if nd['H'] == 1][0]
    reg_mid = [mi for mi, nd in nodes.items() if nd['H'] == max(nd2['H'] for nd2 in nodes.values())][0]

    for source, key in [(trans_mid, 'dist_trans'), (reg_mid, 'dist_reg')]:
        dist = {source: 0}; q = deque([source])
        while q:
            v = q.popleft()
            for u in adj[v]:
                if u not in dist: dist[u] = dist[v]+1; q.append(u)
        for mi in nodes: nodes[mi][key] = dist.get(mi, -1)

    print(f"  V={V}, E={len(edges)}, built in {time.time()-t0:.1f}s")

    # ================================================================
    # PRINT FULL NODE PROFILES
    # ================================================================
    print(f"\n  COMPLETE NODE PROFILES:")
    header = (f"  {'mid':>3} {'H':>3} {'SC':>2} {'Aut':>3} {'fib':>4} {'gs':>3} {'c3':>3} "
              f"{'deg':>3} {'spi':>3} {'rib':>3} {'sea':>3} "
              f"{'blu':>3} {'blk':>3} {'mix':>3} "
              f"{'amb':>3} {'vio':>3} {'tea':>3} {'mag':>3} "
              f"{'stp':>3} {'lp':>3} {'lv':>2} "
              f"{'SL':>4} {'dT':>2} {'dR':>2} {'pal':>3} {'svar':>5} {'score'}")
    print(header)

    for mi in sorted(nodes.keys(), key=lambda x: (nodes[x]['H'], x)):
        nd = nodes[mi]
        avg_nh = nd['sum_nbr_H']/nd['deg'] if nd['deg'] > 0 else 0
        print(f"  {mi:3d} {nd['H']:3d} {'S' if nd['sc'] else 'N':>2} {nd['aut']:3d} "
              f"{nd['fiber']:4d} {nd['gs_count']:3d} {nd['c3']:3d} "
              f"{nd['deg']:3d} {nd['spine_deg']:3d} {nd['rib_deg']:3d} {nd['sea_deg']:3d} "
              f"{nd['blue_deg']:3d} {nd['black_deg']:3d} {nd['mixed_deg']:3d} "
              f"{nd['amber_deg']:3d} {nd['violet_deg']:3d} "
              f"{nd['teal_deg']:3d} {nd['magenta_deg']:3d} "
              f"{nd['step_deg']:3d} {nd['leap_deg']:2d} {nd['level_deg']:2d} "
              f"{nd['self_loops']:4d} {nd['dist_trans']:2d} {nd['dist_reg']:2d} "
              f"{'Y' if nd['palindromic'] else 'N':>3} {nd['score_var']:5.2f} {nd['score']}")

    # ================================================================
    # CORRELATIONS
    # ================================================================
    print(f"\n  CORRELATIONS (Pearson r) between node features:")
    features = ['H', 'c3', 'aut', 'fiber', 'gs_count', 'deg', 'spine_deg', 'rib_deg',
                'sea_deg', 'amber_deg', 'teal_deg', 'step_deg', 'self_loops',
                'score_var', 'dist_trans', 'dist_reg']
    vals = {f: np.array([nodes[mi][f] for mi in sorted(nodes.keys())], dtype=float) for f in features}

    strong = []
    for i, f1 in enumerate(features):
        for f2 in features[i+1:]:
            v1, v2 = vals[f1], vals[f2]
            if np.std(v1) > 0.001 and np.std(v2) > 0.001:
                r = np.corrcoef(v1, v2)[0,1]
                if abs(r) > 0.7:
                    strong.append((abs(r), f1, f2, r))

    strong.sort(reverse=True)
    print(f"  Strong correlations (|r| > 0.7):")
    for ar, f1, f2, r in strong[:20]:
        print(f"    {f1:>15} vs {f2:<15}: r = {r:+.3f}")

    # ================================================================
    # FORMULAS: exact relationships
    # ================================================================
    print(f"\n  EXACT RELATIONSHIPS:")

    # fiber = H/|Aut| (already verified)
    all_match = all(abs(nodes[mi]['fiber'] - nodes[mi]['H']/nodes[mi]['aut']) < 0.01 for mi in nodes)
    print(f"    fiber = H/|Aut|: {'EXACT' if all_match else 'FAILS'}")

    # deg = spine + rib + sea
    all_match = all(nodes[mi]['deg'] == nodes[mi]['spine_deg'] + nodes[mi]['rib_deg'] + nodes[mi]['sea_deg'] for mi in nodes)
    print(f"    deg = spine + rib + sea: {'EXACT' if all_match else 'FAILS'}")

    # blue + black + mixed = deg
    all_match = all(nodes[mi]['deg'] == nodes[mi]['blue_deg'] + nodes[mi]['black_deg'] + nodes[mi]['mixed_deg'] for mi in nodes)
    print(f"    deg = blue + black + mixed: {'EXACT' if all_match else 'FAILS'}")

    # amber + violet = deg
    all_match = all(nodes[mi]['deg'] == nodes[mi]['amber_deg'] + nodes[mi]['violet_deg'] for mi in nodes)
    print(f"    deg = amber + violet: {'EXACT' if all_match else 'FAILS'}")

    # teal + magenta = deg
    all_match = all(nodes[mi]['deg'] == nodes[mi]['teal_deg'] + nodes[mi]['magenta_deg'] for mi in nodes)
    print(f"    deg = teal + magenta: {'EXACT' if all_match else 'FAILS'}")

    # step + leap + level = deg
    all_match = all(nodes[mi]['deg'] == nodes[mi]['step_deg'] + nodes[mi]['leap_deg'] + nodes[mi]['level_deg'] for mi in nodes)
    print(f"    deg = step + leap + level: {'EXACT' if all_match else 'FAILS'}")

    # SC nodes: rib_deg > 0? sea_deg = 0?
    sc_nodes = [mi for mi in nodes if nodes[mi]['sc']]
    ns_nodes = [mi for mi in nodes if not nodes[mi]['sc']]
    sc_sea_zero = all(nodes[mi]['sea_deg'] == 0 for mi in sc_nodes)
    ns_spine_zero = all(nodes[mi]['spine_deg'] == 0 for mi in ns_nodes)
    print(f"    SC nodes have sea_deg=0: {'YES' if sc_sea_zero else 'NO'}")
    print(f"    NS nodes have spine_deg=0: {'YES' if ns_spine_zero else 'NO'}")

    # dist_trans + dist_reg >= diameter?
    diam = max(nd['dist_trans'] for nd in nodes.values())
    print(f"    Diameter (from transitive): {diam}")
    print(f"    dist_trans + dist_reg >= diameter for all? "
          f"{'YES' if all(nodes[mi]['dist_trans']+nodes[mi]['dist_reg']>=diam for mi in nodes) else 'NO'}")

    # ================================================================
    # PATTERNS: grid-sym fraction by SC type
    # ================================================================
    print(f"\n  GRID-SYMMETRIC FRACTION:")
    for label, node_set in [('SC', sc_nodes), ('NS', ns_nodes)]:
        fibs = [nodes[mi]['fiber'] for mi in node_set]
        gss = [nodes[mi]['gs_count'] for mi in node_set]
        if sum(fibs) > 0:
            print(f"    {label}: total fiber={sum(fibs)}, grid-sym={sum(gss)}, "
                  f"fraction={sum(gss)/sum(fibs):.4f}")

    # ================================================================
    # SELF-LOOP ANALYSIS
    # ================================================================
    print(f"\n  SELF-LOOP PATTERNS:")
    print(f"    Self-loops / (fiber * m):")
    for mi in sorted(nodes.keys(), key=lambda x: nodes[x]['H']):
        nd = nodes[mi]
        sl_frac = nd['self_loops'] / (nd['fiber'] * m) if nd['fiber'] > 0 else 0
        if nd['sc'] or nd['self_loops'] > 0:
            print(f"      mid={mi:3d} H={nd['H']:3d} {'SC' if nd['sc'] else 'NS'} "
                  f"fiber={nd['fiber']:4d} SL={nd['self_loops']:4d} "
                  f"SL/(fiber*m)={sl_frac:.4f} "
                  f"neutral_arcs={nd['self_loops']/nd['fiber']:.1f}/{m}")

    # ================================================================
    # UNIQUE PROFILES: are any two nodes indistinguishable?
    # ================================================================
    print(f"\n  PROFILE UNIQUENESS:")
    profile_keys = ['H', 'sc', 'c3', 'aut', 'deg', 'spine_deg', 'rib_deg', 'sea_deg']
    profiles = Counter()
    for mi in nodes:
        key = tuple(nodes[mi][k] for k in profile_keys)
        profiles[key] += 1
    dups = {k: v for k, v in profiles.items() if v > 1}
    if dups:
        print(f"    Duplicate profiles (same {profile_keys}):")
        for key, count in dups.items():
            matching = [mi for mi in nodes if tuple(nodes[mi][k] for k in profile_keys) == key]
            print(f"      {dict(zip(profile_keys, key))}: {count} nodes = {matching}")
    else:
        print(f"    All {V} nodes have UNIQUE profiles with {profile_keys}!")

    # What's the MINIMUM feature set that distinguishes all nodes?
    for size in range(1, len(profile_keys)+1):
        from itertools import combinations as comb_iter
        for subset in comb_iter(profile_keys, size):
            profiles_sub = Counter()
            for mi in nodes:
                key = tuple(nodes[mi][k] for k in subset)
                profiles_sub[key] += 1
            if all(v == 1 for v in profiles_sub.values()):
                print(f"    MINIMUM distinguishing set: {subset} (size {size})")
                break
        else:
            continue
        break

print(f"\n  DONE.")
print("=" * 80)
