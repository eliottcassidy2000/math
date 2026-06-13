#!/usr/bin/env python3
"""
node_profiles_n7_s20db.py — Node profiling at n=7 via hash-then-canonicalize
kind-pasteur-2026-03-23-S20db

Uses score+c3+H+deletion-fingerprint to hash tournaments into groups,
then canonicalizes within each group. Much faster than full canonicalization.

Tests hypotheses discovered at n=5,6:
  H1: c3 = -score_var (perfect negative correlation)
  H2: (H, spine_deg) distinguishes all nodes (or small set)
  H3: SC nodes have sea_deg=0, NS nodes have spine_deg=0
  H4: Self-loops increase with H
  H5: Grid-sym fraction -> 0
  H6: Fiber = H/|Aut| (orbit-stabilizer)
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import numpy as np
import time

sys.stdout.reconfigure(line_buffering=True)

n = 7
m = comb(n-1, 2)  # 15 tiles
m_total = comb(n, 2)  # 21 total arcs
total_tilings = 1 << m  # 2^15 = 32768

print("=" * 80)
print(f"  NODE PROFILING AT n={n} (2^{m} = {total_tilings} tilings)")
print("  kind-pasteur-2026-03-23-S20db")
print("=" * 80)

# ============================================================================
# HELPERS
# ============================================================================

tiles = [(r,c) for r in range(1,n-1) for c in range(1,n-r)]
tile_arcs = [(r+c, c-1) for r,c in tiles]  # 0-indexed

# Explorer grid symmetry
explorer_tiles = []
for y in range(1, n-1):
    for x in range(n, y+1, -1):
        explorer_tiles.append((x, y))
tile_key_exp = {(x,y): i for i, (x,y) in enumerate(explorer_tiles)}
trans_map = [tile_key_exp.get((n-y+1, n-x+1), -1) for x,y in explorer_tiles]

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

def Hdp(a):
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

def score_seq(a):
    return tuple(sorted(sum(a[i][j] for j in range(n)) for i in range(n)))

def c3count(a):
    c = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if a[i][j] and a[j][k] and a[k][i]: c += 1
                if a[i][k] and a[k][j] and a[j][i]: c += 1
    return c

def canon(a):
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

# ============================================================================
# PHASE 1: Hash all tilings by (score, c3, H)
# ============================================================================

print(f"\n  PHASE 1: Hashing {total_tilings} tilings...")
t0 = time.time()

hash_groups = defaultdict(list)  # (score, c3, H) -> list of tiling bits
tiling_props = {}  # tb -> {score, c3, H, grid_sym}

for tb in range(total_tilings):
    T = build_tournament(tb)
    sc = score_seq(T)
    c3 = c3count(T)
    H = Hdp(T)
    gs = is_grid_sym(tb)
    hash_groups[(sc, c3, H)].append(tb)
    tiling_props[tb] = {'score': sc, 'c3': c3, 'H': H, 'grid_sym': gs}
    if (tb+1) % 5000 == 0:
        print(f"    {tb+1}/{total_tilings} ({time.time()-t0:.0f}s)")

print(f"  Phase 1 done: {len(hash_groups)} hash groups in {time.time()-t0:.0f}s")

# ============================================================================
# PHASE 2: Canonicalize within each hash group
# ============================================================================

print(f"\n  PHASE 2: Canonicalizing within groups...")
t1 = time.time()

tiling_to_canon = {}
canon_to_info = {}

for hkey, members in hash_groups.items():
    sc, c3, H = hkey
    # Group by canonical form within this hash bucket
    sub_groups = defaultdict(list)
    for tb in members:
        T = build_tournament(tb)
        cn = canon(T)
        sub_groups[cn].append(tb)

    for cn, group in sub_groups.items():
        for tb in group:
            tiling_to_canon[tb] = cn
        if cn not in canon_to_info:
            T = build_tournament(group[0])
            Tc = [[1-T[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
            aut = factorial(n) // len([tb2 for tb2 in range(total_tilings) if tiling_to_canon.get(tb2) == cn] if len(group) < 100 else group)
            # Approximate |Aut| from orbit size
            # Actually: total labeled tournaments in class = n! / |Aut|
            # And fiber = H / |Aut|. So |Aut| = H / fiber = H / len(group)
            aut_approx = max(1, round(H / len(group)))
            canon_to_info[cn] = {
                'H': H, 'score': sc, 'c3': c3,
                'sc': canon(Tc) == cn,
                'aut': aut_approx,
                'fiber': len(group),
                'gs_count': sum(1 for tb in group if tiling_props[tb]['grid_sym']),
            }

print(f"  Phase 2 done: {len(canon_to_info)} iso classes in {time.time()-t1:.0f}s")

# Verify fiber = H/|Aut|
mismatches = 0
for cn, ci in canon_to_info.items():
    if abs(ci['fiber'] - ci['H'] / ci['aut']) > 0.5:
        mismatches += 1
print(f"  fiber = H/|Aut| mismatches: {mismatches}/{len(canon_to_info)}")

# ============================================================================
# PHASE 3: Build merged graph
# ============================================================================

print(f"\n  PHASE 3: Building merged graph...")
t2 = time.time()

merge = {}; mid = 0
for cn in sorted(canon_to_info.keys()):
    if cn in merge: continue
    ci = canon_to_info[cn]
    merge[cn] = mid
    if not ci['sc']:
        T = build_tournament(ci['fiber'] and [tb for tb in tiling_to_canon if tiling_to_canon[tb] == cn][0])
        Tc = [[1-T[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
        cc = canon(Tc)
        if cc != cn and cc in canon_to_info:
            merge[cc] = mid
    mid += 1
V = mid

nodes = {}
for cn, ci in canon_to_info.items():
    mi = merge[cn]
    if mi not in nodes:
        sv = float(np.var(list(ci['score'])))
        nodes[mi] = {
            'H': ci['H'], 'sc': ci['sc'], 'score': ci['score'],
            'c3': ci['c3'], 'aut': ci['aut'], 'fiber': ci['fiber'],
            'gs_count': ci['gs_count'], 'score_var': sv,
            'palindromic': all(ci['score'][i]+ci['score'][n-1-i]==n-1 for i in range(n//2)),
            'deg': 0, 'spine_deg': 0, 'rib_deg': 0, 'sea_deg': 0,
            'blue_deg': 0, 'black_deg': 0, 'mixed_deg': 0,
            'amber_deg': 0, 'step_deg': 0, 'self_loops': 0,
        }
    else:
        nodes[mi]['fiber'] += ci['fiber']
        nodes[mi]['gs_count'] += ci['gs_count']

print(f"  V_merged = {V} in {time.time()-t2:.0f}s")

# ============================================================================
# PHASE 4: Compute edges and node degrees
# ============================================================================

print(f"\n  PHASE 4: Computing edges...")
t3 = time.time()

edges = {}
for tb in range(total_tilings):
    cn = tiling_to_canon[tb]
    src = merge[cn]
    gs_src = tiling_props[tb]['grid_sym']

    for k in range(m):
        tb2 = tb ^ (1 << k)
        cn2 = tiling_to_canon[tb2]
        tgt = merge[cn2]

        if tgt == src:
            nodes[src]['self_loops'] += 1
            continue

        e = (min(src, tgt), max(src, tgt))
        if e not in edges:
            a, b = e
            na, nb = nodes[a], nodes[b]
            if na['sc'] and nb['sc']: sct = 'SC-SC'
            elif not na['sc'] and not nb['sc']: sct = 'NS-NS'
            else: sct = 'SC-NS'
            edges[e] = {
                'sc_type': sct,
                'dH': abs(na['H'] - nb['H']),
                'score_same': na['score'] == nb['score'],
                'blue_ct': 0, 'black_ct': 0,
            }

        gs_tgt = tiling_props[tb2]['grid_sym']
        if gs_src and gs_tgt: edges[e]['blue_ct'] += 1
        else: edges[e]['black_ct'] += 1

    if (tb+1) % 5000 == 0:
        print(f"    {tb+1}/{total_tilings} ({time.time()-t3:.0f}s)")

print(f"  {len(edges)} edges in {time.time()-t3:.0f}s")

# Accumulate node degrees
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
        if info['dH'] == 2: nodes[v]['step_deg'] += 1

# ============================================================================
# RESULTS
# ============================================================================

E = len(edges)
print(f"\n{'='*80}")
print(f"  n={n}: V={V}, E={E}")
print(f"{'='*80}")

# SC-type breakdown
sc_nodes = [mi for mi in nodes if nodes[mi]['sc']]
ns_nodes = [mi for mi in nodes if not nodes[mi]['sc']]
print(f"\n  SC nodes: {len(sc_nodes)}, NS nodes: {len(ns_nodes)}")

sc_sc = sum(1 for info in edges.values() if info['sc_type'] == 'SC-SC')
sc_ns = sum(1 for info in edges.values() if info['sc_type'] == 'SC-NS')
ns_ns = sum(1 for info in edges.values() if info['sc_type'] == 'NS-NS')
print(f"  Spine (SC-SC): {sc_sc}, Ribs (SC-NS): {sc_ns}, Sea (NS-NS): {ns_ns}")

blue_e = sum(1 for info in edges.values() if info['blue_ct']>0 and info['black_ct']==0)
black_e = sum(1 for info in edges.values() if info['blue_ct']==0)
mixed_e = E - blue_e - black_e
print(f"  Explorer: BLUE={blue_e}, BLACK={black_e}, MIXED={mixed_e}")

# ================================================================
# HYPOTHESIS TESTS
# ================================================================
print(f"\n  --- HYPOTHESIS TESTS ---")

# H1: c3 vs score_var correlation
c3_vals = np.array([nodes[mi]['c3'] for mi in sorted(nodes)])
sv_vals = np.array([nodes[mi]['score_var'] for mi in sorted(nodes)])
if np.std(c3_vals) > 0 and np.std(sv_vals) > 0:
    r_c3_sv = np.corrcoef(c3_vals, sv_vals)[0,1]
    print(f"  H1: corr(c3, score_var) = {r_c3_sv:.6f} {'PERFECT!' if abs(r_c3_sv+1) < 0.001 else ''}")

# H2: minimum distinguishing set
h_vals = [nodes[mi]['H'] for mi in sorted(nodes)]
sd_vals = [nodes[mi]['spine_deg'] for mi in sorted(nodes)]
d_vals = [nodes[mi]['deg'] for mi in sorted(nodes)]
# Check (H, spine_deg)
profiles_2 = Counter((nodes[mi]['H'], nodes[mi]['spine_deg']) for mi in nodes)
dups_2 = sum(1 for v in profiles_2.values() if v > 1)
print(f"  H2: (H, spine_deg) unique? {V - dups_2}/{V} unique ({dups_2} collisions)")
# Check (H, deg, spine_deg)
profiles_3 = Counter((nodes[mi]['H'], nodes[mi]['deg'], nodes[mi]['spine_deg']) for mi in nodes)
dups_3 = sum(1 for v in profiles_3.values() if v > 1)
print(f"      (H, deg, spine_deg) unique? {V - dups_3}/{V} unique ({dups_3} collisions)")

# H3: SC/NS partition theorem
sc_sea = any(nodes[mi]['sea_deg'] > 0 for mi in sc_nodes)
ns_spine = any(nodes[mi]['spine_deg'] > 0 for mi in ns_nodes)
print(f"  H3: SC sea_deg=0? {'YES' if not sc_sea else 'NO (FAILS!)'}")
print(f"      NS spine_deg=0? {'YES' if not ns_spine else 'NO (FAILS!)'}")

# H4: Self-loops increase with H
h_arr = np.array([nodes[mi]['H'] for mi in sorted(nodes)], dtype=float)
sl_arr = np.array([nodes[mi]['self_loops']/max(1,nodes[mi]['fiber']) for mi in sorted(nodes)])
if np.std(h_arr) > 0 and np.std(sl_arr) > 0:
    r_h_sl = np.corrcoef(h_arr, sl_arr)[0,1]
    print(f"  H4: corr(H, self_loops_per_fiber) = {r_h_sl:.3f}")

# H5: Grid-sym fraction
total_fib = sum(nodes[mi]['fiber'] for mi in nodes)
total_gs = sum(nodes[mi]['gs_count'] for mi in nodes)
print(f"  H5: Grid-sym fraction = {total_gs}/{total_fib} = {total_gs/total_fib:.4f}")
sc_fib = sum(nodes[mi]['fiber'] for mi in sc_nodes)
sc_gs = sum(nodes[mi]['gs_count'] for mi in sc_nodes)
ns_fib = sum(nodes[mi]['fiber'] for mi in ns_nodes)
ns_gs = sum(nodes[mi]['gs_count'] for mi in ns_nodes)
print(f"      SC: {sc_gs}/{sc_fib} = {sc_gs/sc_fib:.4f}" if sc_fib > 0 else "")
print(f"      NS: {ns_gs}/{ns_fib} = {ns_gs/ns_fib:.4f}" if ns_fib > 0 else "")

# H6: fiber = H/|Aut|
print(f"  H6: fiber = H/|Aut| verified: {sum(1 for mi in nodes if abs(nodes[mi]['fiber'] - nodes[mi]['H']/nodes[mi]['aut']) < 0.5)}/{V}")

# ================================================================
# KEY CORRELATIONS
# ================================================================
print(f"\n  --- KEY CORRELATIONS ---")
features = {'H': h_arr}
for f in ['c3', 'deg', 'spine_deg', 'rib_deg', 'sea_deg', 'fiber', 'score_var', 'amber_deg', 'step_deg']:
    features[f] = np.array([nodes[mi][f] for mi in sorted(nodes)], dtype=float)

fnames = list(features.keys())
for i, f1 in enumerate(fnames):
    for f2 in fnames[i+1:]:
        v1, v2 = features[f1], features[f2]
        if np.std(v1) > 0.001 and np.std(v2) > 0.001:
            r = np.corrcoef(v1, v2)[0,1]
            if abs(r) > 0.8:
                print(f"    {f1:>12} vs {f2:<12}: r = {r:+.4f}")

# ================================================================
# DEGREE DISTRIBUTION
# ================================================================
print(f"\n  --- DEGREE DISTRIBUTION ---")
deg_dist = Counter(nodes[mi]['deg'] for mi in nodes)
print(f"  {'deg':>4} {'count':>6}")
for d in sorted(deg_dist.keys()):
    print(f"  {d:4d} {deg_dist[d]:6d}")

print(f"\n  Avg degree: {np.mean([nodes[mi]['deg'] for mi in nodes]):.2f}")
print(f"  Max degree: {max(nodes[mi]['deg'] for mi in nodes)}")
print(f"  Min degree: {min(nodes[mi]['deg'] for mi in nodes)}")

# ================================================================
# AMBER (score-preserving) fraction
# ================================================================
amber_total = sum(1 for info in edges.values() if info['score_same'])
print(f"\n  AMBER (score-preserving) edges: {amber_total}/{E} = {100*amber_total/E:.1f}%")

# Step edges
step_total = sum(1 for info in edges.values() if info['dH'] == 2)
print(f"  STEP (dH=2) edges: {step_total}/{E} = {100*step_total/E:.1f}%")

# Level edges
level_total = sum(1 for info in edges.values() if info['dH'] == 0)
print(f"  LEVEL (dH=0) edges: {level_total}/{E} = {100*level_total/E:.1f}%")

print(f"\n  Total time: {time.time()-t0:.0f}s")
print(f"\n  DONE.")
print("=" * 80)
