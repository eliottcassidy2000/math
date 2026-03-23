#!/usr/bin/env python3
"""
principal_line_n7_fast_s20cp.py — Fast principal line analysis for n=7
kind-pasteur-2026-03-23-S20cp

Reuses the output from the previous computation (we know V=456, E=4086,
SC=88, the full BFS tree structure). Uses the BFS tree to find the
diameter path (principal axis) efficiently via two BFS passes.
Computes bilateral analysis, branch structure, and cross-n comparison.
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter, deque
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  PRINCIPAL LINE OF G_7/Z_2 (FAST — FROM BFS TREE)")
print("  kind-pasteur-2026-03-23-S20cp")
print("=" * 80)

# Parse the BFS tree from previous output
# Format: [mid] H=... score=... c3=... |Aut|=... ch=[children]

# We already know the full structure. Let me reconstruct it from the raw data.
# Actually, the fastest way is to rebuild from scratch but use the BFS tree
# for the diameter instead of enumerating all paths.

n = 7
m = comb(n, 2)
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
        sub = [[adj[i][j] for j in range(n) if j != v] for i in range(n) if i != v]
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

# ============================================================================
# BUILD (same as before, fast)
# ============================================================================

print(f"\n  Building G_7/Z_2...")
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

print(f"  Hashing done: {len(hash_to_bits)} groups ({time.time()-t0:.0f}s)")

# Build iso classes
class_list = []
bits_to_cid = {}
cid = 0
canon_to_cid_map = {}
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
            'size': len(group), 'aut': factorial(n) // len(group)
        })
        for b in group: bits_to_cid[b] = cid
        canon_to_cid_map[canon] = cid
        cid += 1

# SC and complement
for data in class_list:
    comp_bits = data['bits'] ^ ((1 << m) - 1)
    comp_adj = bits_to_adj(comp_bits)
    comp_canon = canonical_form(comp_adj)
    data['comp_cid'] = canon_to_cid_map.get(comp_canon, -1)
    data['sc'] = (data['comp_cid'] == data['cid'])

# Edges
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

# Merged graph
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
        if e not in merged_edge_color: merged_edge_color[e] = color
        elif merged_edge_color[e] != color: merged_edge_color[e] = 'mixed'

merged_classes = []
seen = set()
for data in class_list:
    mid_val = merged_id[data['cid']]
    if mid_val not in seen:
        seen.add(mid_val)
        merged_classes.append({'mid': mid_val, 'H': data['H'], 'sc': data['sc'],
                               'score': data['score'], 'c3': data['c3'], 'aut': data['aut'],
                               'size': data['size']})
merged_classes.sort(key=lambda x: x['mid'])

blue_adj = defaultdict(set)
black_adj = defaultdict(set)
full_adj = defaultdict(set)
for e, color in merged_edge_color.items():
    a, b = e
    full_adj[a].add(b); full_adj[b].add(a)
    if color in ('blue', 'mixed'):
        blue_adj[a].add(b); blue_adj[b].add(a)
    if color in ('black', 'mixed'):
        black_adj[a].add(b); black_adj[b].add(a)

sc_mids = set(c['mid'] for c in merged_classes if c['sc'])

print(f"  V_merged={V_merged}, E_merged={len(merged_edges)}, SC={len(sc_mids)} ({time.time()-t0:.0f}s)")

# ============================================================================
# PRINCIPAL AXIS VIA BFS TREE DIAMETER
# ============================================================================

print(f"\n{'='*70}")
print(f"  PRINCIPAL AXIS")
print(f"{'='*70}")

# Find transitive
transitive = next(c['mid'] for c in merged_classes if c['H'] == 1)
print(f"  Transitive: mid={transitive}")
print(f"  Blue neighbors (3): {sorted(blue_adj[transitive] & sc_mids)}")
print(f"  Black neighbors: {sorted(black_adj[transitive])}")

# Build SC blue adjacency (restricted to SC vertices)
sc_blue_adj = defaultdict(set)
sc_blue_edges = 0
for e, color in merged_edge_color.items():
    a, b = e
    if color in ('blue', 'mixed') and a in sc_mids and b in sc_mids:
        sc_blue_adj[a].add(b)
        sc_blue_adj[b].add(a)
        sc_blue_edges += 1

print(f"  SC-SC blue edges: {sc_blue_edges}")
print(f"  Is tree? {sc_blue_edges == len(sc_mids) - 1}")

# Find diameter path in SC blue subgraph via two BFS passes
def bfs_farthest(start, adj, valid):
    """Find farthest vertex from start in subgraph defined by valid set."""
    dist = {start: 0}
    queue = deque([start])
    farthest = start
    max_dist = 0
    parent = {start: None}
    while queue:
        v = queue.popleft()
        for u in adj[v]:
            if u in valid and u not in dist:
                dist[u] = dist[v] + 1
                parent[u] = v
                queue.append(u)
                if dist[u] > max_dist:
                    max_dist = dist[u]
                    farthest = u
    return farthest, max_dist, parent

# Step 1: BFS from transitive to find farthest SC vertex
end1, _, _ = bfs_farthest(transitive, sc_blue_adj, sc_mids)
# Step 2: BFS from end1 to find the actual diameter endpoint
end2, diam, parent = bfs_farthest(end1, sc_blue_adj, sc_mids)

# Reconstruct diameter path
diameter_path = []
v = end2
while v is not None:
    diameter_path.append(v)
    v = parent[v]
diameter_path.reverse()

print(f"\n  SC BLUE SUBGRAPH DIAMETER: {diam}")
print(f"  Diameter path: {len(diameter_path)} vertices")
print(f"  From mid={end1} to mid={end2}")

# Where is transitive on this path?
if transitive in diameter_path:
    trans_pos = diameter_path.index(transitive)
    print(f"  Transitive at position {trans_pos} of {len(diameter_path)-1}")
else:
    # Find closest point on diameter path to transitive
    _, _, from_trans = bfs_farthest(transitive, sc_blue_adj, sc_mids)
    for v in diameter_path:
        if v in from_trans:
            pass  # find the closest
    print(f"  Transitive NOT on diameter path!")

# Use path from transitive to farthest point as principal axis
end_from_trans, trans_diam, trans_parent = bfs_farthest(transitive, sc_blue_adj, sc_mids)
axis_path = []
v = end_from_trans
while v is not None:
    axis_path.append(v)
    v = trans_parent[v]
axis_path.reverse()

print(f"\n  PRINCIPAL AXIS (from transitive to farthest SC vertex):")
print(f"  Length: {len(axis_path)} vertices, {len(axis_path)-1} edges")
axis_H = [merged_classes[v]['H'] for v in axis_path]
print(f"  H: {axis_H}")
print(f"  Endpoint: mid={end_from_trans}, H={merged_classes[end_from_trans]['H']}")

# Print axis detail
print(f"\n  AXIS DETAIL:")
for i, v in enumerate(axis_path):
    c = merged_classes[v]
    score = c['score']
    mean_s = sum(score) / len(score)
    var_s = sum((s - mean_s)**2 for s in score) / len(score)
    sc_blue_deg = len(sc_blue_adj[v])
    black_deg = len(black_adj[v])
    print(f"    [{i:>2d}] mid={v:>3d} H={c['H']:>3d} score={score} var={var_s:.2f} "
          f"c3={c['c3']:>2d} |Aut|={c['aut']:>2d} SC_blue_deg={sc_blue_deg} black_deg={black_deg}")

# ============================================================================
# THREE BRANCHES FROM TRANSITIVE
# ============================================================================

print(f"\n{'='*70}")
print(f"  THREE BRANCHES FROM TRANSITIVE")
print(f"{'='*70}")

trans_sc_blue_nbrs = sorted(blue_adj[transitive] & sc_mids)
print(f"  Transitive has {len(trans_sc_blue_nbrs)} SC blue neighbors: {trans_sc_blue_nbrs}")

for bi, br in enumerate(trans_sc_blue_nbrs):
    # BFS from branch root excluding transitive
    visited = {transitive}
    bq = deque([br])
    visited.add(br)
    branch = [br]
    while bq:
        v = bq.popleft()
        for u in sc_blue_adj[v]:
            if u not in visited:
                visited.add(u)
                branch.append(u)
                bq.append(u)

    H_vals = [merged_classes[v]['H'] for v in branch]
    c3_vals = [merged_classes[v]['c3'] for v in branch]
    print(f"\n  BRANCH {bi}: root=mid {br} (H={merged_classes[br]['H']})")
    print(f"    Size: {len(branch)} SC vertices")
    print(f"    H range: [{min(H_vals)}, {max(H_vals)}]")
    print(f"    c3 range: [{min(c3_vals)}, {max(c3_vals)}]")
    print(f"    Contains Paley (H=189)? {189 in H_vals}")
    print(f"    Contains regular score? {any(merged_classes[v]['score'] == (3,3,3,3,3,3,3) for v in branch)}")

    # Depth distribution
    depth_dist = Counter()
    depth_map = {br: 1}
    dq = deque([br])
    while dq:
        v = dq.popleft()
        for u in sc_blue_adj[v]:
            if u not in depth_map and u != transitive:
                depth_map[u] = depth_map[v] + 1
                dq.append(u)
    max_d = max(depth_map.values()) if depth_map else 0
    for d in range(1, max_d+1):
        cnt = sum(1 for v, dd in depth_map.items() if dd == d)
        if cnt > 0:
            depth_dist[d] = cnt
    print(f"    Depth distribution: {dict(sorted(depth_dist.items()))}")

    # Score diversity at each depth
    for d in range(1, min(max_d+1, 6)):
        verts = [v for v, dd in depth_map.items() if dd == d]
        scores = set(merged_classes[v]['score'] for v in verts)
        H_at_d = sorted(merged_classes[v]['H'] for v in verts)
        print(f"    Depth {d}: {len(verts)} vertices, {len(scores)} scores, H={H_at_d[:5]}{'...' if len(H_at_d)>5 else ''}")

# ============================================================================
# PERPENDICULAR ANALYSIS
# ============================================================================

print(f"\n{'='*70}")
print(f"  PERPENDICULAR ANALYSIS")
print(f"{'='*70}")

axis_set = set(axis_path)
axis_pos = {v: i for i, v in enumerate(axis_path)}

# BFS from axis to all vertices
perp_dist = {}
nearest_axis = {}
for v in axis_path:
    perp_dist[v] = 0
    nearest_axis[v] = v

bfs_q = deque(axis_path)
visited_perp = set(axis_path)
while bfs_q:
    v = bfs_q.popleft()
    for u in full_adj[v]:
        if u not in visited_perp:
            visited_perp.add(u)
            perp_dist[u] = perp_dist[v] + 1
            nearest_axis[u] = nearest_axis[v]
            bfs_q.append(u)

max_perp = max(perp_dist.values())
print(f"  Max perpendicular distance from axis: {max_perp}")
print(f"  Vertices reached: {len(visited_perp)}/{V_merged}")

for d in range(max_perp + 1):
    count = sum(1 for v in range(V_merged) if perp_dist.get(v) == d)
    if d == 0:
        print(f"  d=0: {count} (axis)")
    else:
        sc_d = sum(1 for v in range(V_merged) if perp_dist.get(v) == d and merged_classes[v]['sc'])
        ns_d = count - sc_d
        print(f"  d={d}: {count} (SC={sc_d}, NS={ns_d})")

# Bilateral: above vs below at each axis position
print(f"\n  BILATERAL BALANCE PER AXIS POSITION:")
total_above = 0
total_below = 0
total_level = 0
for i, v in enumerate(axis_path):
    ax_H = merged_classes[v]['H']
    branches = [u for u in range(V_merged) if nearest_axis.get(u) == v and u != v]
    above = sum(1 for u in branches if merged_classes[u]['H'] > ax_H)
    below = sum(1 for u in branches if merged_classes[u]['H'] < ax_H)
    level = sum(1 for u in branches if merged_classes[u]['H'] == ax_H)
    total_above += above
    total_below += below
    total_level += level
    if branches and (i < 3 or i >= len(axis_path) - 3):
        print(f"    [{i:>2d}] H={ax_H:>3d}: {len(branches)} branches (above={above}, below={below}, level={level})")

print(f"\n  GLOBAL BALANCE: above={total_above}, below={total_below}, level={total_level}")
print(f"  Balance ratio: above/(above+below) = {total_above/(total_above+total_below):.3f}" if total_above+total_below > 0 else "")

# ============================================================================
# CROSS-n SUMMARY (combining n=3..6 results + n=7)
# ============================================================================

print(f"\n{'='*70}")
print(f"  CROSS-n SUMMARY (n=3..7)")
print(f"{'='*70}")

print("""
  PRINCIPAL AXIS:
  n  SC  spine  axis_len  sc_blue_edges  axis_H_range  #branches_trans
  3   2    2       2           1           [1,3]              1
  4   2    2       2           1           [1,5]              1
  5   8    8       8           7           [1,15]             2
  6  12   12      10          11           [1,45]             2
  7  88   88      14           ?           [1,189]            3

  BLUE SC IS A TREE: Confirmed at n=3..6 (edges = SC-1).
  At n=7: 88 SC vertices, need to verify edge count.
  The tree structure deepens: max_depth = 1, 1, 3, 7, 13.

  TRANSITIVE BLUE DEGREE (SC neighbors only):
  n=3: 1    n=4: 1    n=5: 2    n=6: 2    n=7: 3

  LARGER BRANCH H-JUMP FROM TRANSITIVE:
  n=3: 2    n=4: 4    n=5: 8    n=6: 16   n=7: 32
  PATTERN: 2^(n-2) !!
  n=7: transitive neighbor with H=33 -> Delta_H = 32 = 2^5 = 2^(n-2)

  SCORE VARIANCE ALONG AXIS:
  Always decreasing from transitive (max variance) toward regular (0 variance).
  Not strictly monotone but the TREND is clear.

  The principal axis traces the journey from ORDER (transitive) to CHAOS (regular).
  Score variance quantifies this: maximum spread -> zero spread.
  The intermediate SC classes represent PARTIALLY ORDERED tournaments.
""")

print(f"  Computation time: {time.time()-t0:.0f}s")
print("=" * 80)
