#!/usr/bin/env python3
"""
principal_line_n7_s20cp.py — Principal Line analysis for G_7/Z_2
kind-pasteur-2026-03-23-S20cp

Uses fast hash (score, c3, H, deletion_fingerprint) for 2^21 tournaments.
Builds full G_7/Z_2, then traces the SC blue spine from the transitive class.
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter, deque
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  PRINCIPAL LINE OF G_7/Z_2")
print("  kind-pasteur-2026-03-23-S20cp")
print("=" * 80)

n = 7
m = comb(n, 2)  # 21
total = 1 << m

PAIRS = [(i, j) for i in range(n) for j in range(i+1, n)]

def bits_to_adj(bits):
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(PAIRS):
        if bits & (1 << k): adj[i][j] = 1
        else: adj[j][i] = 1
    return adj

def H_dp_n(adj, nn):
    dp = {}
    for v in range(nn): dp[(1 << v, v)] = 1
    for S in range(1, 1 << nn):
        for v in range(nn):
            if not (S & (1 << v)): continue
            val = dp.get((S, v), 0)
            if val == 0: continue
            for u in range(nn):
                if S & (1 << u): continue
                if adj[v][u]:
                    key = (S | (1 << u), u)
                    dp[key] = dp.get(key, 0) + val
    full = (1 << nn) - 1
    return sum(dp.get((full, v), 0) for v in range(nn))

def score_seq(adj):
    return tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))

def c3_count(adj):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]: c3 += 1
                if adj[i][k] and adj[k][j] and adj[j][i]: c3 += 1
    return c3

def deletion_fingerprint(adj):
    fps = []
    for v in range(n):
        sub = []
        for i in range(n):
            if i == v: continue
            row = []
            for j in range(n):
                if j == v: continue
                row.append(adj[i][j])
            sub.append(row)
        fps.append(H_dp_n(sub, n-1))
    return tuple(sorted(fps))

def canonical_form(adj):
    scores = [sum(adj[i][j] for j in range(n)) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[scores[v]].append(v)
    groups = [sg[s] for s in sorted(set(scores))]
    if all(len(g) == 1 for g in groups):
        perm = [g[0] for g in groups]
        return tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p) + r
    for perm in gp(groups):
        form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best: best = form
    return best

def complement_bits(bits):
    return bits ^ ((1 << m) - 1)

# ============================================================================
# PHASE 1: Build iso classes via hash
# ============================================================================

print(f"\n  PHASE 1: Hashing {total} tournaments...")
t0 = time.time()

hash_to_bits = defaultdict(list)
for bits in range(total):
    adj = bits_to_adj(bits)
    ss = score_seq(adj)
    c3 = c3_count(adj)
    h = H_dp_n(adj, n)
    dfp = deletion_fingerprint(adj)
    hash_to_bits[(ss, c3, h, dfp)].append(bits)
    if (bits+1) % 200000 == 0:
        print(f"    {bits+1}/{total} ({time.time()-t0:.0f}s)")

print(f"  Phase 1: {len(hash_to_bits)} groups in {time.time()-t0:.0f}s")

# ============================================================================
# PHASE 2: Build iso classes with canonical forms
# ============================================================================

print(f"  PHASE 2: Building iso classes...")
t1 = time.time()

class_list = []
bits_to_cid = {}
cid = 0
for hkey, members in hash_to_bits.items():
    canons = {}
    for b in members:
        adj = bits_to_adj(b)
        canon = canonical_form(adj)
        if canon not in canons: canons[canon] = []
        canons[canon].append(b)
    for canon, group in canons.items():
        adj = bits_to_adj(group[0])
        ss, c3, h, dfp = hkey
        class_list.append({
            'cid': cid, 'bits': group[0], 'adj': adj,
            'H': h, 'score': ss, 'c3': c3, 'canon': canon,
            'size': len(group), 'members': group,
            'aut': factorial(n) // len(group)
        })
        for b in group: bits_to_cid[b] = cid
        cid += 1

V = len(class_list)
print(f"  V = {V} (expected 456) in {time.time()-t1:.0f}s")

# ============================================================================
# PHASE 3: SC and complement pairing
# ============================================================================

print(f"  PHASE 3: SC and complement pairing...")
canon_to_cid_map = {d['canon']: d['cid'] for d in class_list}
for data in class_list:
    comp_bits = complement_bits(data['bits'])
    comp_adj = bits_to_adj(comp_bits)
    comp_canon = canonical_form(comp_adj)
    data['comp_cid'] = canon_to_cid_map.get(comp_canon, -1)
    data['sc'] = (data['comp_cid'] == data['cid'])

sc_count = sum(1 for d in class_list if d['sc'])
print(f"  SC = {sc_count}")

# ============================================================================
# PHASE 4: Build edges
# ============================================================================

print(f"  PHASE 4: Computing edges...")
t3 = time.time()

edges = set()
edge_colors = {}
for data in class_list:
    ci = data['cid']
    bits = data['bits']
    for ai in range(m):
        flipped = bits ^ (1 << ai)
        nb = bits_to_cid.get(flipped)
        if nb is not None and nb != ci:
            e = (min(ci, nb), max(ci, nb))
            if e not in edges:
                edges.add(e)
                edge_colors[e] = 'blue' if data['sc'] == class_list[nb]['sc'] else 'black'

E = len(edges)
blue_count = sum(1 for c in edge_colors.values() if c == 'blue')
print(f"  E = {E}, blue = {blue_count}, black = {E - blue_count} ({time.time()-t3:.0f}s)")

# ============================================================================
# PHASE 5: Build merged graph
# ============================================================================

print(f"  PHASE 5: Building merged graph...")
merged_id = {}
mid = 0
for data in class_list:
    ci = data['cid']
    if ci in merged_id: continue
    comp = data['comp_cid']
    merged_id[ci] = mid
    if comp != ci: merged_id[comp] = mid
    mid += 1

V_merged = mid
merged_edges = set()
merged_edge_color = {}
for (a, b) in edges:
    ma, mb = merged_id[a], merged_id[b]
    if ma != mb:
        e = (min(ma, mb), max(ma, mb))
        merged_edges.add(e)
        color = edge_colors.get((min(a,b), max(a,b)), 'blue')
        if e not in merged_edge_color:
            merged_edge_color[e] = color
        elif merged_edge_color[e] != color:
            merged_edge_color[e] = 'mixed'

merged_classes = []
seen = set()
for data in class_list:
    mid_val = merged_id[data['cid']]
    if mid_val not in seen:
        seen.add(mid_val)
        merged_classes.append({
            'mid': mid_val, 'H': data['H'], 'sc': data['sc'],
            'score': data['score'], 'c3': data['c3'], 'aut': data['aut'],
            'size': data['size']
        })
merged_classes.sort(key=lambda x: x['mid'])

blue_adj = defaultdict(set)
black_adj = defaultdict(set)
full_adj = defaultdict(set)
for e, color in merged_edge_color.items():
    a, b = e
    full_adj[a].add(b); full_adj[b].add(a)
    if color == 'blue':
        blue_adj[a].add(b); blue_adj[b].add(a)
    elif color == 'black':
        black_adj[a].add(b); black_adj[b].add(a)
    else:
        blue_adj[a].add(b); blue_adj[b].add(a)
        black_adj[a].add(b); black_adj[b].add(a)

E_merged = len(merged_edges)
blue_m = sum(1 for c in merged_edge_color.values() if c == 'blue')
black_m = sum(1 for c in merged_edge_color.values() if c == 'black')
print(f"  V_merged={V_merged}, E_merged={E_merged}, blue={blue_m}, black={black_m}")

# ============================================================================
# PHASE 6: PRINCIPAL LINE ANALYSIS
# ============================================================================

print(f"\n{'='*70}")
print(f"  PRINCIPAL LINE OF G_7/Z_2")
print(f"{'='*70}")

# Find transitive class
transitive = None
for c in merged_classes:
    if c['H'] == 1:
        transitive = c['mid']
        break

tc = merged_classes[transitive]
print(f"\n  TRANSITIVE: mid={transitive}, H=1, score={tc['score']}, sc={tc['sc']}, |Aut|={tc['aut']}")
print(f"  Blue neighbors: {sorted(blue_adj[transitive])}")
print(f"  Black neighbors: {sorted(black_adj[transitive])}")

# SC blue spine BFS
sc_mids = set(c['mid'] for c in merged_classes if c['sc'])
spine_parent = {transitive: None}
spine_depth = {transitive: 0}
queue = deque([transitive])
spine_order = [transitive]
while queue:
    v = queue.popleft()
    for u in blue_adj[v]:
        if u in spine_parent or u not in sc_mids: continue
        spine_parent[u] = v
        spine_depth[u] = spine_depth[v] + 1
        queue.append(u)
        spine_order.append(u)

print(f"\n  SC BLUE SPINE: {len(spine_order)} vertices, max depth {max(spine_depth.values())}")

# Print spine tree structure (compact)
print(f"\n  SPINE TREE (BFS order):")
for v in spine_order:
    c = merged_classes[v]
    d = spine_depth[v]
    children = [u for u in blue_adj[v] if u in spine_parent and spine_parent[u] == v]
    indent = "  " * d
    print(f"    {indent}[{v}] H={c['H']:>3d} score={c['score']} c3={c['c3']:>2d} "
          f"|Aut|={c['aut']:>2d} ch={children}")

# Find longest path in blue SC subgraph from transitive
def find_all_paths(start, adj_d, valid):
    all_p = []
    stack = [(start, [start])]
    while stack:
        v, path = stack.pop()
        ext = False
        for u in adj_d[v]:
            if u in valid and u not in set(path):
                stack.append((u, path + [u]))
                ext = True
        if not ext:
            all_p.append(path)
    return all_p

all_paths = find_all_paths(transitive, blue_adj, sc_mids)
longest = max(all_paths, key=len)
print(f"\n  PRINCIPAL AXIS: {len(longest)} vertices, {len(longest)-1} edges")
axis_H = [merged_classes[v]['H'] for v in longest]
print(f"  Path: {longest}")
print(f"  H: {axis_H}")

# Score and c3 along axis
print(f"\n  AXIS DETAIL:")
for i, v in enumerate(longest):
    c = merged_classes[v]
    score = c['score']
    mean_s = sum(score) / len(score)
    var_s = sum((s - mean_s)**2 for s in score) / len(score)
    print(f"    [{i:>2d}] mid={v:>3d} H={c['H']:>3d} score={score} var={var_s:.2f} c3={c['c3']:>2d} |Aut|={c['aut']:>2d}")

axis_set = set(longest)
axis_pos = {v: i for i, v in enumerate(longest)}

# Perpendicular structure: BFS from axis into full graph
perp_dist = {}
nearest_axis = {}
for v in longest:
    perp_dist[v] = 0
    nearest_axis[v] = v

bfs_q = deque(longest)
visited = set(longest)
while bfs_q:
    v = bfs_q.popleft()
    for u in full_adj[v]:
        if u not in visited:
            visited.add(u)
            perp_dist[u] = perp_dist[v] + 1
            nearest_axis[u] = nearest_axis[v]
            bfs_q.append(u)

# Group by (axis_position, perp_dist)
perp_groups = defaultdict(list)
for v in range(V_merged):
    ax = nearest_axis.get(v)
    if ax is not None:
        perp_groups[(axis_pos.get(ax, -1), perp_dist.get(v, 0))].append(v)

# Summary: how many vertices at each perpendicular distance from each axis point
print(f"\n  PERPENDICULAR BRANCHING SUMMARY:")
max_perp = max(perp_dist.values()) if perp_dist else 0
print(f"  Max perpendicular distance: {max_perp}")
for d in range(max_perp + 1):
    count = sum(1 for v in range(V_merged) if perp_dist.get(v) == d and v not in axis_set) if d > 0 else len(axis_set)
    if d == 0:
        print(f"  d={d}: {count} (axis itself)")
    else:
        sc_at_d = sum(1 for v in range(V_merged) if perp_dist.get(v) == d and v not in axis_set and merged_classes[v]['sc'])
        ns_at_d = count - sc_at_d
        print(f"  d={d}: {count} (SC={sc_at_d}, NS={ns_at_d})")

# Bilateral analysis: for each axis vertex, count branches and their properties
print(f"\n  BILATERAL ANALYSIS PER AXIS VERTEX:")
total_above = 0
total_below = 0
total_level = 0
for i, v in enumerate(longest):
    ax_H = merged_classes[v]['H']
    branches = []
    for d in range(1, 10):
        branches.extend(perp_groups.get((i, d), []))

    if not branches: continue

    above = [u for u in branches if merged_classes[u]['H'] > ax_H]
    below = [u for u in branches if merged_classes[u]['H'] < ax_H]
    at_level = [u for u in branches if merged_classes[u]['H'] == ax_H]
    total_above += len(above)
    total_below += len(below)
    total_level += len(at_level)

    sc_b = sum(1 for u in branches if merged_classes[u]['sc'])
    ns_b = len(branches) - sc_b

    print(f"    Axis[{i:>2d}] H={ax_H:>3d}: {len(branches)} branches "
          f"(SC={sc_b},NS={ns_b}) above={len(above)} below={len(below)} level={len(at_level)}")

print(f"\n  GLOBAL BILATERAL BALANCE: above={total_above}, below={total_below}, level={total_level}")
print(f"  Balance ratio: {total_above}/{total_below+total_level+total_above} = "
      f"{total_above/(total_above+total_below+total_level):.3f}")

# Black edges from axis
print(f"\n  BLACK CONNECTIVITY FROM AXIS:")
total_black_from_axis = 0
for i, v in enumerate(longest):
    bn = black_adj[v]
    total_black_from_axis += len(bn)
    if i < 5 or i >= len(longest) - 5:  # show first and last 5
        on_ax = sum(1 for u in bn if u in axis_set)
        off_ax = len(bn) - on_ax
        print(f"    Axis[{i:>2d}] H={merged_classes[v]['H']:>3d}: {len(bn)} black "
              f"(on_axis={on_ax}, off={off_ax})")

# Score variance profile
print(f"\n  SCORE VARIANCE ALONG AXIS:")
variances = []
for v in longest:
    score = merged_classes[v]['score']
    mean_s = sum(score) / len(score)
    var_s = sum((s - mean_s)**2 for s in score) / len(score)
    variances.append(var_s)
print(f"  Variances: {[f'{v:.2f}' for v in variances]}")
print(f"  Monotone decreasing? {all(variances[i] >= variances[i+1] for i in range(len(variances)-1))}")

# Check: does the axis visit ALL SC vertices?
sc_on_axis = len(axis_set & sc_mids)
sc_off_axis = len(sc_mids - axis_set)
print(f"\n  SC COVERAGE: {sc_on_axis} on axis, {sc_off_axis} off axis (total SC = {len(sc_mids)})")
if sc_off_axis > 0:
    for v in sc_mids - axis_set:
        c = merged_classes[v]
        d = perp_dist.get(v, -1)
        print(f"    Off-axis SC: mid={v} H={c['H']} score={c['score']} perp_dist={d}")

# H-differences along axis
diffs = [axis_H[i+1] - axis_H[i] for i in range(len(axis_H)-1)]
print(f"\n  H-DIFFERENCES ALONG AXIS: {diffs}")
print(f"  Max jump: {max(abs(d) for d in diffs)}")
print(f"  Uphill steps: {sum(1 for d in diffs if d > 0)}")
print(f"  Downhill steps: {sum(1 for d in diffs if d < 0)}")
print(f"  Level steps: {sum(1 for d in diffs if d == 0)}")

# Degree along axis (in full merged graph)
print(f"\n  DEGREE ALONG AXIS:")
for i, v in enumerate(longest):
    d_total = len(full_adj[v])
    d_blue = len(blue_adj[v])
    d_black = len(black_adj[v])
    print(f"    [{i:>2d}] H={merged_classes[v]['H']:>3d} deg={d_total:>3d} (blue={d_blue:>3d}, black={d_black:>3d})")

# ============================================================================
# PHASE 7: THE TWO BRANCHES FROM TRANSITIVE
# ============================================================================

print(f"\n{'='*70}")
print(f"  THE TWO BRANCHES FROM TRANSITIVE")
print(f"{'='*70}")

trans_blue_nbrs = sorted(blue_adj[transitive])
print(f"  Transitive has {len(trans_blue_nbrs)} blue neighbors: {trans_blue_nbrs}")

for bi, branch_root in enumerate(trans_blue_nbrs):
    print(f"\n  BRANCH {bi}: root = {branch_root}")
    # BFS from branch root through blue SC edges (excluding transitive)
    branch_visited = {transitive, branch_root}
    branch_order = [branch_root]
    bq = deque([branch_root])
    while bq:
        v = bq.popleft()
        for u in blue_adj[v]:
            if u not in branch_visited and u in sc_mids:
                branch_visited.add(u)
                branch_order.append(u)
                bq.append(u)

    for v in branch_order:
        c = merged_classes[v]
        d_from_trans = spine_depth.get(v, -1)
        # Count NS classes attached via black
        ns_attached = sum(1 for u in black_adj[v] if not merged_classes[u]['sc'])
        print(f"    [{v:>3d}] H={c['H']:>3d} depth={d_from_trans} score={c['score']} "
              f"c3={c['c3']:>2d} |Aut|={c['aut']:>2d} NS_attached={ns_attached}")

    print(f"  Branch {bi} size: {len(branch_order)} SC vertices")
    print(f"  Branch {bi} H range: [{min(merged_classes[v]['H'] for v in branch_order)}, "
          f"{max(merged_classes[v]['H'] for v in branch_order)}]")

print(f"\n  Total time: {time.time()-t0:.0f}s")
print("=" * 80)
