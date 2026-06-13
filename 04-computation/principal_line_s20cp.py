#!/usr/bin/env python3
"""
principal_line_s20cp.py — The Principal Line of Symmetry in G_n/Z_2
kind-pasteur-2026-03-23-S20cp

THE IDEA: Start from the transitive tournament (H=1, unique minimum).
It is self-complementary, so it's a vertex in G_n/Z_2.
It connects via blue edges (SC-preserving) to other SC classes.
The chain of SC classes connected by blue edges forms a SPINE —
the "principal line of symmetry."

All other classes (NS merged vertices) hang off this spine via black edges.
We organize them on either "side" of the spine and look for bilateral balance.

For each n=3..6:
  1. Build G_n/Z_2
  2. Extract the SC blue spine (BFS tree from transitive through blue edges)
  3. Find the longest blue path from transitive (the "principal axis")
  4. Compute perpendicular structure: how NS classes attach to the spine
  5. Define left/right sides and check bilateral symmetry
  6. Compute invariant profiles along the spine and perpendicular to it
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter, deque
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  THE PRINCIPAL LINE OF SYMMETRY IN G_n/Z_2")
print("  kind-pasteur-2026-03-23-S20cp")
print("=" * 80)

# ============================================================================
# TOURNAMENT HELPERS
# ============================================================================

def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def H_dp(adj, n):
    dp = [0] * ((1 << n) * n)
    for v in range(n): dp[(1 << v) * n + v] = 1
    for S in range(1, 1 << n):
        if bin(S).count('1') >= n: continue
        for v in range(n):
            if not (S & (1 << v)): continue
            val = dp[S * n + v]
            if val == 0: continue
            for u in range(n):
                if S & (1 << u): continue
                if adj[v][u]:
                    dp[(S | (1 << u)) * n + u] += val
    full = (1 << n) - 1
    return sum(dp[full * n + v] for v in range(n))

def canonical_form(adj, n):
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

def complement_adj(adj, n):
    return [[1-adj[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))

def c3_count(adj, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]: c3 += 1
                if adj[i][k] and adj[k][j] and adj[j][i]: c3 += 1
    return c3

def c5_count(adj, n):
    """Count directed 5-cycles."""
    c5 = 0
    for vs in combinations(range(n), 5):
        sub = [[adj[vs[i]][vs[j]] for j in range(5)] for i in range(5)]
        for perm in permutations(range(5)):
            if all(sub[perm[i]][perm[(i+1)%5]] for i in range(5)):
                c5 += 1
    return c5 // 5  # each 5-cycle counted 5 times (cyclic rotations, but not reflections)
    # Actually: each directed 5-cycle is counted once per starting vertex = 5 times
    # But we iterate over all 5! perms of each 5-subset, finding 2 directed cycles per undirected cycle
    # Correction: directed 5-cycles have 5 rotations, so c5 = count / 5
    # Wait, actually for permutations: a 5-cycle perm is found 5 times (cyclic shifts).
    # Each directed 5-cycle in the tournament gives exactly 5 permutation matches.
    # So c5 = total_matches / 5. But we also don't want to double-count the reverse direction.
    # Actually for tournaments, each unordered 5-cycle gives exactly 2 directed 5-cycles
    # (one in each direction), but only the one matching the tournament arcs is counted.
    # So the counting above is correct: count directed 5-cycles, each appears 5 times among perms.

def aut_size(adj, n):
    count = 0
    scores = [sum(adj[i][j] for j in range(n)) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[scores[v]].append(v)
    groups = [sg[s] for s in sorted(set(scores))]
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p) + r
    for perm in gp(groups):
        ok = True
        for i in range(n):
            for j in range(i+1, n):
                if adj[perm[i]][perm[j]] != adj[i][j]: ok = False; break
            if not ok: break
        if ok: count += 1
    return count

# ============================================================================
# BUILD G_n/Z_2 WITH FULL STRUCTURE
# ============================================================================

def build_merged_full(n_val):
    n = n_val
    m = comb(n, 2)
    total = 1 << m
    t0 = time.time()

    # Group by canonical form
    iso_classes = defaultdict(list)
    for bits in range(total):
        adj = tournament_adj(n, bits)
        canon = canonical_form(adj, n)
        iso_classes[canon].append(bits)

    class_list = []
    canon_to_cid = {}
    for idx, (canon, members) in enumerate(sorted(iso_classes.items())):
        adj = tournament_adj(n, members[0])
        h = H_dp(adj, n)
        comp = complement_adj(adj, n)
        sc = canonical_form(comp, n) == canon
        ss = score_seq(adj, n)
        c3 = c3_count(adj, n)
        au = aut_size(adj, n)
        class_list.append({
            'cid': idx, 'canon': canon, 'members': members,
            'adj': adj, 'H': h, 'sc': sc, 'score': ss,
            'c3': c3, 'aut': au, 'size': len(members)
        })
        canon_to_cid[canon] = idx

    # Complement pairing
    for data in class_list:
        comp = complement_adj(data['adj'], n)
        comp_canon = canonical_form(comp, n)
        data['comp_cid'] = canon_to_cid.get(comp_canon, -1)

    # Edges with colors
    edges = set()
    edge_colors = {}
    self_loops = set()
    for data in class_list:
        cid = data['cid']
        for bits in data['members']:
            a2 = tournament_adj(n, bits)
            for ai in range(m):
                flipped = bits ^ (1 << ai)
                fa = tournament_adj(n, flipped)
                fc = canonical_form(fa, n)
                nb = canon_to_cid.get(fc)
                if nb is None: continue
                if nb == cid:
                    self_loops.add(cid)
                else:
                    e = (min(cid, nb), max(cid, nb))
                    edges.add(e)
                    if e not in edge_colors:
                        edge_colors[e] = 'blue' if data['sc'] == class_list[nb]['sc'] else 'black'

    # Build merged graph
    merged_id = {}
    mid = 0
    for data in class_list:
        cid = data['cid']
        if cid in merged_id: continue
        comp = data['comp_cid']
        merged_id[cid] = mid
        if comp != cid: merged_id[comp] = mid
        mid += 1

    V_merged = mid
    merged_edges = set()
    merged_edge_color = {}
    collapsed = 0

    for (a, b) in edges:
        ma, mb = merged_id[a], merged_id[b]
        if ma == mb:
            collapsed += 1
        else:
            e = (min(ma, mb), max(ma, mb))
            merged_edges.add(e)
            color = edge_colors.get((min(a,b), max(a,b)), 'blue')
            if e not in merged_edge_color:
                merged_edge_color[e] = color
            elif merged_edge_color[e] != color:
                merged_edge_color[e] = 'mixed'

    # Merged class data
    merged_classes = []
    seen = set()
    for data in class_list:
        mid_val = merged_id[data['cid']]
        if mid_val not in seen:
            seen.add(mid_val)
            merged_classes.append({
                'mid': mid_val, 'H': data['H'], 'sc': data['sc'],
                'score': data['score'], 'c3': data['c3'],
                'aut': data['aut'], 'size': data['size'],
                'cid': data['cid'], 'comp_cid': data['comp_cid']
            })

    # Sort merged_classes by mid
    merged_classes.sort(key=lambda x: x['mid'])

    # Adjacency structures
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
        else:  # mixed
            blue_adj[a].add(b); blue_adj[b].add(a)
            black_adj[a].add(b); black_adj[b].add(a)

    print(f"  G_{n}/Z_2: V={V_merged}, E={len(merged_edges)}, "
          f"SC={sum(1 for c in merged_classes if c['sc'])}, "
          f"blue={sum(1 for c in merged_edge_color.values() if c=='blue')}, "
          f"black={sum(1 for c in merged_edge_color.values() if c=='black')} "
          f"({time.time()-t0:.1f}s)")

    return {
        'n': n, 'V': V_merged, 'classes': merged_classes,
        'edges': merged_edges, 'edge_colors': merged_edge_color,
        'blue_adj': blue_adj, 'black_adj': black_adj, 'full_adj': full_adj,
        'merged_id': merged_id, 'orig_classes': class_list
    }


# ============================================================================
# PRINCIPAL LINE ANALYSIS
# ============================================================================

def principal_line_analysis(data):
    n = data['n']
    V = data['V']
    classes = data['classes']
    blue_adj = data['blue_adj']
    black_adj = data['black_adj']
    full_adj = data['full_adj']

    print(f"\n{'='*70}")
    print(f"  PRINCIPAL LINE OF G_{n}/Z_2")
    print(f"{'='*70}")

    # 1. FIND THE TRANSITIVE CLASS (H=1, minimum)
    transitive = None
    for c in classes:
        if c['H'] == 1:
            transitive = c['mid']
            break
    if transitive is None:
        print("  ERROR: No transitive class found!")
        return None

    tc = classes[transitive]
    print(f"\n  TRANSITIVE CLASS: mid={transitive}, H=1, score={tc['score']}, "
          f"sc={tc['sc']}, |Aut|={tc['aut']}, size={tc['size']}")
    print(f"  Blue neighbors: {sorted(blue_adj[transitive])}")
    print(f"  Black neighbors: {sorted(black_adj[transitive])}")

    # 2. TRACE THE SC BLUE SPINE (BFS from transitive through blue edges, SC only)
    sc_mids = set(c['mid'] for c in classes if c['sc'])
    print(f"\n  SC vertices: {sorted(sc_mids)}")

    # BFS through blue edges restricted to SC vertices
    spine_parent = {transitive: None}
    spine_depth = {transitive: 0}
    queue = deque([transitive])
    spine_order = [transitive]

    while queue:
        v = queue.popleft()
        for u in blue_adj[v]:
            if u in spine_parent: continue
            if u not in sc_mids: continue  # stay on SC vertices
            spine_parent[u] = v
            spine_depth[u] = spine_depth[v] + 1
            queue.append(u)
            spine_order.append(u)

    print(f"\n  SC BLUE SPINE (BFS tree from transitive):")
    print(f"  Spine size: {len(spine_order)} vertices")
    max_depth = max(spine_depth.values()) if spine_depth else 0
    print(f"  Max depth: {max_depth}")

    # Print spine tree
    for d in range(max_depth + 1):
        verts_at_d = [v for v in spine_order if spine_depth[v] == d]
        for v in verts_at_d:
            c = classes[v]
            indent = "  " * (d + 1)
            parent_str = f"<- {spine_parent[v]}" if spine_parent[v] is not None else "(root)"
            children = [u for u in blue_adj[v] if u in spine_parent and spine_parent[u] == v]
            print(f"  {indent}[{v}] H={c['H']:>3d} score={c['score']} c3={c['c3']} "
                  f"|Aut|={c['aut']} {parent_str}  children={children}")

    # 3. FIND THE PRINCIPAL AXIS (longest path from transitive in the spine)
    # DFS to find all paths from transitive to leaves
    def find_all_paths(start, adj_dict, valid_set):
        """Find all simple paths from start within valid_set."""
        all_paths = []
        stack = [(start, [start])]
        while stack:
            v, path = stack.pop()
            extended = False
            for u in adj_dict[v]:
                if u in valid_set and u not in set(path):
                    stack.append((u, path + [u]))
                    extended = True
            if not extended:
                all_paths.append(path)
        return all_paths

    all_spine_paths = find_all_paths(transitive, blue_adj, sc_mids)
    longest_path = max(all_spine_paths, key=len)
    print(f"\n  PRINCIPAL AXIS (longest blue path from transitive):")
    print(f"  Length: {len(longest_path)} vertices, {len(longest_path)-1} edges")
    print(f"  Path: {longest_path}")
    print(f"  H along path: {[classes[v]['H'] for v in longest_path]}")
    print(f"  Scores along path:")
    for i, v in enumerate(longest_path):
        c = classes[v]
        print(f"    [{i}] mid={v} H={c['H']:>3d} score={c['score']} c3={c['c3']} |Aut|={c['aut']}")

    # 4. CLASSIFY ALL VERTICES RELATIVE TO THE PRINCIPAL AXIS
    axis_set = set(longest_path)
    axis_positions = {v: i for i, v in enumerate(longest_path)}

    # For non-axis vertices: find nearest axis vertex and perpendicular distance
    # BFS from axis
    perp_dist = {}
    nearest_axis = {}
    for v in longest_path:
        perp_dist[v] = 0
        nearest_axis[v] = v

    bfs_queue = deque(longest_path)
    visited = set(longest_path)
    while bfs_queue:
        v = bfs_queue.popleft()
        for u in full_adj[v]:
            if u not in visited:
                visited.add(u)
                perp_dist[u] = perp_dist[v] + 1
                nearest_axis[u] = nearest_axis[v]
                bfs_queue.append(u)

    # Group vertices by their nearest axis vertex and perpendicular distance
    perp_groups = defaultdict(list)  # (axis_pos, perp_dist) -> list of mid
    for v in range(V):
        if v in axis_set:
            perp_groups[(axis_positions[v], 0)].append(v)
        else:
            ax = nearest_axis.get(v)
            if ax is not None:
                perp_groups[(axis_positions[ax], perp_dist[v])].append(v)

    print(f"\n  PERPENDICULAR STRUCTURE (vertices grouped by axis position and distance):")
    for ax_pos in range(len(longest_path)):
        ax_v = longest_path[ax_pos]
        ax_c = classes[ax_v]
        line = f"  Axis[{ax_pos}] H={ax_c['H']:>3d}: "
        # Collect perpendicular branches
        for d in range(0, 10):
            verts = perp_groups.get((ax_pos, d), [])
            if d == 0:
                continue  # axis vertex itself
            if verts:
                v_strs = []
                for v in sorted(verts, key=lambda x: classes[x]['H']):
                    c = classes[v]
                    sc_tag = "SC" if c['sc'] else "NS"
                    v_strs.append(f"{v}(H={c['H']},{sc_tag})")
                line += f"  d={d}: {', '.join(v_strs)}"
        if any(perp_groups.get((ax_pos, d), []) for d in range(1, 10)):
            print(line)

    # 5. BILATERAL SYMMETRY ANALYSIS
    # Define "side" based on the connection type:
    # For NS classes attached to axis vertex v_ax:
    #   The NS class represents a complement pair (T, T^comp).
    #   T has score s, T^comp has score s^rev.
    #   We define "left" = score lex < complement score, "right" = lex >=
    # For SC classes not on axis:
    #   These are on the spine but not on the principal path.
    #   Define side by whether they branch left or right in the BFS tree.

    print(f"\n  BILATERAL ANALYSIS:")

    # For each axis vertex, count left vs right branches
    axis_branch_balance = []
    for ax_pos in range(len(longest_path)):
        ax_v = longest_path[ax_pos]
        # All non-axis vertices whose nearest axis point is ax_v
        branches = []
        for d in range(1, 10):
            for v in perp_groups.get((ax_pos, d), []):
                branches.append(v)

        if not branches:
            axis_branch_balance.append((0, 0, 0))
            continue

        # Split branches by H relative to axis H
        ax_H = classes[ax_v]['H']
        above = [v for v in branches if classes[v]['H'] > ax_H]
        below = [v for v in branches if classes[v]['H'] < ax_H]
        at_level = [v for v in branches if classes[v]['H'] == ax_H]
        axis_branch_balance.append((len(above), len(below), len(at_level)))

        if branches:
            branch_H = sorted([classes[v]['H'] for v in branches])
            sc_count = sum(1 for v in branches if classes[v]['sc'])
            ns_count = len(branches) - sc_count
            print(f"    Axis[{ax_pos}] H={ax_H}: {len(branches)} branches "
                  f"(SC={sc_count}, NS={ns_count}), "
                  f"above={len(above)}, below={len(below)}, level={len(at_level)}")
            print(f"      Branch H values: {branch_H}")

    # 6. SCORE PROFILE ALONG AXIS
    print(f"\n  SCORE PROFILE ALONG PRINCIPAL AXIS:")
    for i, v in enumerate(longest_path):
        c = classes[v]
        score = c['score']
        # Score variance
        mean_s = sum(score) / len(score)
        var_s = sum((s - mean_s)**2 for s in score) / len(score)
        # Score entropy
        probs = Counter(score)
        ent = -sum((cnt/len(score)) * np.log2(cnt/len(score)) for cnt in probs.values() if cnt > 0)
        print(f"    [{i}] H={c['H']:>3d} score={score} var={var_s:.2f} ent={ent:.3f} c3={c['c3']}")

    # 7. EDGE PATTERN ALONG AXIS
    print(f"\n  EDGE PATTERN ALONG PRINCIPAL AXIS:")
    for i in range(len(longest_path) - 1):
        u, v = longest_path[i], longest_path[i+1]
        # What is the H-difference?
        h_diff = classes[v]['H'] - classes[u]['H']
        # How many common neighbors?
        common = full_adj[u] & full_adj[v]
        common_on_axis = common & axis_set
        common_off_axis = common - axis_set
        print(f"    [{i}]->[{i+1}]: H {classes[u]['H']}->{classes[v]['H']} "
              f"(Delta={h_diff:+d}), "
              f"common_nbrs={len(common)} (on_axis={len(common_on_axis)}, off={len(common_off_axis)})")

    # 8. BLACK EDGE PATTERN FROM AXIS
    print(f"\n  BLACK EDGES FROM AXIS VERTICES:")
    for i, v in enumerate(longest_path):
        black_nbrs = black_adj[v]
        if black_nbrs:
            nbr_info = []
            for u in sorted(black_nbrs, key=lambda x: classes[x]['H']):
                c = classes[u]
                sc_tag = "SC" if c['sc'] else "NS"
                on_ax = "ON" if u in axis_set else "off"
                nbr_info.append(f"{u}(H={c['H']},{sc_tag},{on_ax})")
            print(f"    Axis[{i}] H={classes[v]['H']:>3d}: {len(black_nbrs)} black -> "
                  f"{', '.join(nbr_info)}")

    # 9. THE "STAIRCASE" WITHIN THE STAIRCASE
    # Each axis vertex has a score. Plot scores to see if they form a staircase.
    print(f"\n  SCORE STAIRCASE ALONG AXIS:")
    for i, v in enumerate(longest_path):
        score = classes[v]['score']
        bar = ''.join(['#' if s > (n-1)/2 else '.' if s < (n-1)/2 else '|' for s in score])
        print(f"    [{i:>2d}] H={classes[v]['H']:>3d} {score} |{bar}|")

    # 10. SYMMETRY METRIC: for each non-axis class, is there a "mirror" class?
    # Mirror: same H, same perpendicular distance, connected to same axis vertex
    print(f"\n  MIRROR ANALYSIS (non-axis classes):")
    non_axis = [v for v in range(V) if v not in axis_set]
    mirror_count = 0
    no_mirror = 0
    for v in non_axis:
        ax = nearest_axis.get(v)
        if ax is None: continue
        d = perp_dist[v]
        h = classes[v]['H']
        # Look for mirror: same ax, same d, same H, different v
        siblings = [u for u in perp_groups.get((axis_positions[ax], d), [])
                    if u != v and classes[u]['H'] == h]
        if siblings:
            mirror_count += 1
        else:
            no_mirror += 1

    print(f"  Classes with mirror: {mirror_count}")
    print(f"  Classes without mirror: {no_mirror}")
    print(f"  Mirror fraction: {mirror_count/(mirror_count+no_mirror):.3f}" if (mirror_count+no_mirror) > 0 else "")

    return {
        'n': n, 'transitive': transitive,
        'spine_order': spine_order, 'spine_depth': spine_depth,
        'axis': longest_path,
        'axis_H': [classes[v]['H'] for v in longest_path],
        'perp_groups': dict(perp_groups),
        'balance': axis_branch_balance,
        'mirror_count': mirror_count, 'no_mirror': no_mirror
    }


# ============================================================================
# MAIN
# ============================================================================

all_results = {}

for n in [3, 4, 5, 6]:
    print(f"\n\n{'#'*80}")
    print(f"  n = {n}")
    print(f"{'#'*80}")

    data = build_merged_full(n)
    result = principal_line_analysis(data)
    all_results[n] = result

# ============================================================================
# CROSS-N COMPARISON
# ============================================================================

print(f"\n\n{'='*80}")
print("  CROSS-n COMPARISON")
print(f"{'='*80}")

print(f"\n  PRINCIPAL AXIS PROPERTIES:")
print(f"  {'n':>3} {'|axis|':>7} {'max_depth':>10} {'H_range':>15} {'axis_H':>40}")
for n in [3, 4, 5, 6]:
    r = all_results[n]
    if r is None: continue
    ax = r['axis']
    print(f"  {n:>3} {len(ax):>7} {max(r['spine_depth'].values()):>10} "
          f"{'['+str(r['axis_H'][0])+','+str(r['axis_H'][-1])+']':>15} "
          f"{str(r['axis_H']):>40}")

print(f"\n  SPINE SIZE vs SC COUNT:")
for n in [3, 4, 5, 6]:
    r = all_results[n]
    if r is None: continue
    sc = len(r['spine_order'])
    print(f"  n={n}: spine={sc} SC vertices reachable from transitive via blue")

print(f"\n  BILATERAL BALANCE (above/below/level counts at each axis position):")
for n in [3, 4, 5, 6]:
    r = all_results[n]
    if r is None: continue
    print(f"  n={n}: {r['balance']}")

print(f"\n  MIRROR SYMMETRY:")
for n in [3, 4, 5, 6]:
    r = all_results[n]
    if r is None: continue
    total = r['mirror_count'] + r['no_mirror']
    frac = r['mirror_count'] / total if total > 0 else 0
    print(f"  n={n}: {r['mirror_count']}/{total} classes have mirror ({frac:.3f})")

print(f"\n  H-PROFILE SEQUENCES ALONG PRINCIPAL AXIS:")
for n in [3, 4, 5, 6]:
    r = all_results[n]
    if r is None: continue
    diffs = [r['axis_H'][i+1] - r['axis_H'][i] for i in range(len(r['axis_H'])-1)]
    print(f"  n={n}: H = {r['axis_H']}, Delta_H = {diffs}")

print(f"\n\n  DONE.")
print("=" * 80)
