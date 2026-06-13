#!/usr/bin/env python3
"""
merged_n7_deep_s20co.py — Deep analysis of G_7 and G_7/Z_2
kind-pasteur-2026-03-23-S20co

Uses fast hash approach: (score, c3, H_deletion_fingerprint)
to classify 2^21 = 2,097,152 tournaments into 456 iso classes,
then computes 30+ invariants of G_7 and G_7/Z_2.

Expected runtime: ~10-15 minutes.
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  DEEP ANALYSIS OF G_7 AND G_7/Z_2")
print("  kind-pasteur-2026-03-23-S20co")
print("=" * 80)

n = 7
m = comb(n, 2)  # 21
total = 1 << m   # 2097152

PAIRS = [(i, j) for i in range(n) for j in range(i+1, n)]

print(f"  n={n}, m={m}, |tournaments|={total}")

# ============================================================================
# HELPERS
# ============================================================================

def bits_to_adj(bits):
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(PAIRS):
        if bits & (1 << k):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def H_dp_n(adj, nn):
    dp = {}
    for v in range(nn):
        dp[(1 << v, v)] = 1
    for S in range(1, 1 << nn):
        for v in range(nn):
            if not (S & (1 << v)):
                continue
            val = dp.get((S, v), 0)
            if val == 0:
                continue
            for u in range(nn):
                if S & (1 << u):
                    continue
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
    """Sorted tuple of H(T-v) for each v."""
    fps = []
    for v in range(n):
        sub_adj = []
        for i in range(n):
            if i == v:
                continue
            row = []
            for j in range(n):
                if j == v:
                    continue
                row.append(adj[i][j])
            sub_adj.append(row)
        fps.append(H_dp_n(sub_adj, n-1))
    return tuple(sorted(fps))

def canonical_form(adj):
    scores = [sum(adj[i][j] for j in range(n)) for i in range(n)]
    score_groups = defaultdict(list)
    for v in range(n):
        score_groups[scores[v]].append(v)
    sorted_scores = sorted(set(scores))
    groups = [score_groups[s] for s in sorted_scores]

    if all(len(g) == 1 for g in groups):
        perm = [g[0] for g in groups]
        return tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))

    best = None
    def gen_perms(groups):
        if not groups:
            yield []
            return
        for p in permutations(groups[0]):
            for rest in gen_perms(groups[1:]):
                yield list(p) + rest
    for perm in gen_perms(groups):
        form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

def complement_bits(bits):
    return bits ^ ((1 << m) - 1)

# ============================================================================
# PHASE 1: Hash all tournaments
# ============================================================================

print(f"\n  PHASE 1: Hashing all {total} tournaments...")
t0 = time.time()

hash_to_bits = defaultdict(list)  # hash -> list of bits

for bits in range(total):
    adj = bits_to_adj(bits)
    ss = score_seq(adj)
    c3 = c3_count(adj)
    h = H_dp_n(adj, n)
    # Use deletion fingerprint for disambiguation
    dfp = deletion_fingerprint(adj)

    hkey = (ss, c3, h, dfp)
    hash_to_bits[hkey].append(bits)

    if (bits + 1) % 100000 == 0:
        elapsed = time.time() - t0
        rate = (bits + 1) / elapsed
        eta = (total - bits - 1) / rate
        print(f"    {bits+1}/{total} ({elapsed:.0f}s, ETA {eta:.0f}s, {len(hash_to_bits)} groups)")

print(f"  Phase 1 done: {len(hash_to_bits)} hash groups in {time.time()-t0:.0f}s")

# Check for ambiguity
ambiguous = sum(1 for members in hash_to_bits.values()
                if len(set(canonical_form(bits_to_adj(b)) for b in members)) > 1)
print(f"  Ambiguous groups (multiple iso classes): {ambiguous}")

# ============================================================================
# PHASE 2: Build iso classes
# ============================================================================

print(f"\n  PHASE 2: Building iso classes...")
t1 = time.time()

# If no ambiguity, each hash group = one iso class
class_list = []
bits_to_cid = {}
cid = 0

for hkey, members in hash_to_bits.items():
    # Check if this group has multiple canonical forms
    canons = {}
    for b in members:
        adj = bits_to_adj(b)
        canon = canonical_form(adj)
        if canon not in canons:
            canons[canon] = []
        canons[canon].append(b)

    for canon, group_members in canons.items():
        adj = bits_to_adj(group_members[0])
        ss, c3, h, dfp = hkey
        class_list.append({
            'cid': cid, 'bits': group_members[0], 'adj': adj,
            'H': h, 'score': ss, 'c3': c3, 'canon': canon,
            'size': len(group_members), 'members': group_members,
            'dfp': dfp
        })
        for b in group_members:
            bits_to_cid[b] = cid
        cid += 1

V = len(class_list)
print(f"  |V(G_7)| = {V} (expected 456)")
print(f"  Phase 2 done in {time.time()-t1:.1f}s")

# ============================================================================
# PHASE 3: SC status and complement pairing
# ============================================================================

print(f"\n  PHASE 3: SC and complement pairing...")
t2 = time.time()

canon_to_cid_map = {d['canon']: d['cid'] for d in class_list}

for data in class_list:
    comp_bits = complement_bits(data['bits'])
    comp_adj = bits_to_adj(comp_bits)
    comp_canon = canonical_form(comp_adj)
    comp_cid = canon_to_cid_map.get(comp_canon, -1)
    data['comp_cid'] = comp_cid
    data['sc'] = (comp_cid == data['cid'])

sc_count = sum(1 for d in class_list if d['sc'])
print(f"  SC classes: {sc_count} (expected 88)")
v_merged = (V + sc_count) // 2
print(f"  V_merged = {v_merged} (expected 272)")

# Compute |Aut| = n! / class_size
for data in class_list:
    data['aut'] = factorial(n) // data['size']

print(f"  Phase 3 done in {time.time()-t2:.1f}s")

# ============================================================================
# PHASE 4: Build edges (sample representative + flip)
# ============================================================================

print(f"\n  PHASE 4: Computing edges...")
t3 = time.time()

edges = set()
self_loops = set()
edge_colors = {}

# For each class, flip each of m arcs on the representative
# and look up the resulting class
for data in class_list:
    cid_i = data['cid']
    bits = data['bits']

    for arc_idx in range(m):
        flipped = bits ^ (1 << arc_idx)
        nb_cid = bits_to_cid.get(flipped)
        if nb_cid is None:
            continue
        if nb_cid == cid_i:
            self_loops.add(cid_i)
        else:
            e = (min(cid_i, nb_cid), max(cid_i, nb_cid))
            if e not in edges:
                edges.add(e)
                sc_i = data['sc']
                sc_j = class_list[nb_cid]['sc']
                edge_colors[e] = 'blue' if sc_i == sc_j else 'black'

    if (data['cid'] + 1) % 100 == 0:
        elapsed = time.time() - t3
        print(f"    {data['cid']+1}/{V} ({elapsed:.1f}s, {len(edges)} edges)")

# IMPORTANT: flipping just the representative gives only a SUBSET of edges!
# A class with size k has k representatives, each giving potentially different neighbors.
# We need to also flip arcs on other members, or check all members.

# Actually, re-reading: the representative flip gives ALL neighbors because
# if two iso classes are connected, then SOME representative of one class
# flips to SOME member of the other. Since we stored bits_to_cid for ALL
# tournaments, looking up any flipped tournament works.
# BUT: we only flip the representative's m arcs. The representative may not
# reach all neighboring classes via its m flips.
# Need to check: for representative bits, flipping arc k gives tournament T'.
# T' might be isomorphic to a class that the representative can't reach,
# but another member can.

# Solution: use ALL members of each class (but that's 2M total flips = expensive)
# Alternative: check if current edge count matches expected

print(f"\n  Representative-only edges: {len(edges)}")
print(f"  Expected: 4086")

if len(edges) < 4086:
    print(f"  Missing edges! Using full enumeration for remaining...")
    # For classes with missing edges, try more members
    # Build adj list from what we have
    adj_list = defaultdict(set)
    for (a, b) in edges:
        adj_list[a].add(b)
        adj_list[b].add(a)

    # Check degree vs expected
    # Actually let's just enumerate all members for small classes
    # and sample for large ones
    for data in class_list:
        cid_i = data['cid']
        for bits in data['members']:
            for arc_idx in range(m):
                flipped = bits ^ (1 << arc_idx)
                nb_cid = bits_to_cid.get(flipped)
                if nb_cid is None:
                    continue
                if nb_cid != cid_i:
                    e = (min(cid_i, nb_cid), max(cid_i, nb_cid))
                    if e not in edges:
                        edges.add(e)
                        sc_i = data['sc']
                        sc_j = class_list[nb_cid]['sc']
                        edge_colors[e] = 'blue' if sc_i == sc_j else 'black'

        if (data['cid'] + 1) % 50 == 0:
            elapsed = time.time() - t3
            print(f"    Full pass: {data['cid']+1}/{V} ({elapsed:.1f}s, {len(edges)} edges)")

E = len(edges)
blue = sum(1 for c in edge_colors.values() if c == 'blue')
black = sum(1 for c in edge_colors.values() if c == 'black')
print(f"\n  |E(G_7)| = {E} (expected 4086, blue={blue}, black={black})")
print(f"  Self-loops: {len(self_loops)}")
print(f"  Phase 4 done in {time.time()-t3:.1f}s")

# ============================================================================
# PHASE 5: Build merged graph
# ============================================================================

print(f"\n  PHASE 5: Building merged graph G_7/Z_2...")

merged_id = {}
mid = 0
for data in class_list:
    cid_i = data['cid']
    if cid_i in merged_id:
        continue
    comp = data['comp_cid']
    merged_id[cid_i] = mid
    if comp != cid_i:
        merged_id[comp] = mid
    mid += 1

V_merged = mid
print(f"  V_merged = {V_merged}")

merged_edges = set()
collapsed = 0
merged_edge_color = {}

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

E_merged = len(merged_edges)
twin = E - E_merged - collapsed
blue_m = sum(1 for c in merged_edge_color.values() if c == 'blue')
black_m = sum(1 for c in merged_edge_color.values() if c == 'black')
mixed_m = sum(1 for c in merged_edge_color.values() if c == 'mixed')

print(f"  E_merged = {E_merged}, collapsed = {collapsed}, twin = {twin}")
print(f"  Check: {E_merged}+{collapsed}+{twin} = {E_merged+collapsed+twin} = {E}")
print(f"  Colors: blue={blue_m}, black={black_m}, mixed={mixed_m}")

# Build adjacency matrix
A_merged = np.zeros((V_merged, V_merged), dtype=int)
for (a, b) in merged_edges:
    A_merged[a][b] = 1
    A_merged[b][a] = 1

A_orig = np.zeros((V, V), dtype=int)
for (a, b) in edges:
    A_orig[a][b] = 1
    A_orig[b][a] = 1

# Merged class data
merged_classes = []
seen = set()
for data in class_list:
    mid_val = merged_id[data['cid']]
    if mid_val not in seen:
        seen.add(mid_val)
        merged_classes.append({
            'mid': mid_val, 'H': data['H'], 'sc': data['sc'],
            'score': data['score'], 'c3': data['c3'], 'aut': data['aut']
        })

# ============================================================================
# PHASE 6: INVARIANTS
# ============================================================================

def compute_invariants(A, V, name, classes_data=None):
    print(f"\n  {'='*60}")
    print(f"  INVARIANTS OF {name} (V={V}, E={int(np.sum(A)//2)})")
    print(f"  {'='*60}")

    A_f = A.astype(float)
    E_count = int(np.sum(A) // 2)
    deg = np.sum(A, axis=1).astype(int)

    print(f"  Degree: min={np.min(deg)}, max={np.max(deg)}, avg={np.mean(deg):.3f}")
    print(f"  Density: {2*E_count/(V*(V-1)):.6f}" if V > 1 else "  Density: N/A")

    # Adjacency spectrum
    eigenvalues = np.sort(np.linalg.eigvalsh(A_f))[::-1]
    print(f"  Spectral radius: {eigenvalues[0]:.6f}")
    print(f"  Spectral gap: {eigenvalues[0]-eigenvalues[1]:.6f}" if V > 1 else "")
    print(f"  Graph energy: {np.sum(np.abs(eigenvalues)):.6f}")
    print(f"  Estrada index: {np.sum(np.exp(eigenvalues)):.6f}")
    print(f"  Top 10 eigenvalues: [{', '.join(f'{x:.4f}' for x in eigenvalues[:10])}]")

    # Laplacian
    D = np.diag(deg.astype(float))
    L = D - A_f
    lap_eigs = np.sort(np.linalg.eigvalsh(L))
    alg_conn = float(lap_eigs[1]) if V > 1 else 0
    print(f"  Algebraic connectivity: {alg_conn:.6f}")
    print(f"  Laplacian spectral radius: {float(lap_eigs[-1]):.6f}")

    # Kirchhoff index
    nonzero_lap = lap_eigs[lap_eigs > 1e-10]
    if len(nonzero_lap) > 0:
        kirchhoff = float(V * np.sum(1.0 / nonzero_lap))
        print(f"  Kirchhoff index: {kirchhoff:.4f}")

    # Spanning trees (may overflow)
    if len(nonzero_lap) == V - 1:
        log_st = np.sum(np.log(nonzero_lap)) - np.log(V)
        print(f"  log(spanning trees): {log_st:.4f}")
        if log_st < 100:
            st = round(np.exp(log_st))
            print(f"  Spanning trees: {st}")

    # Connectivity and distances
    dist_matrix = np.full((V, V), V + 1)
    for start in range(V):
        dist_matrix[start][start] = 0
        queue = [start]
        visited = {start}
        d = 0
        while queue:
            next_q = []
            d += 1
            for v in queue:
                for u in range(V):
                    if A[v][u] > 0 and u not in visited:
                        visited.add(u)
                        dist_matrix[start][u] = d
                        next_q.append(u)
            queue = next_q

    reachable = dist_matrix[dist_matrix <= V]
    diameter = int(np.max(reachable))
    wiener = int(np.sum(dist_matrix[dist_matrix <= V]) // 2)
    avg_dist = float(np.mean(dist_matrix[(dist_matrix > 0) & (dist_matrix <= V)]))
    print(f"  Diameter: {diameter}")
    print(f"  Wiener index: {wiener}")
    print(f"  Avg distance: {avg_dist:.4f}")

    # Triangles
    A2 = A_f @ A_f
    triangles = int(np.trace(A_f @ A2) / 6)
    print(f"  Triangles: {triangles}")

    # Clustering coefficient
    cc_sum = 0
    cc_count = 0
    for v in range(V):
        nbrs = [u for u in range(V) if A[v][u] > 0]
        k = len(nbrs)
        if k >= 2:
            t = sum(1 for i in range(len(nbrs)) for j in range(i+1, len(nbrs))
                    if A[nbrs[i]][nbrs[j]] > 0)
            cc_sum += 2 * t / (k * (k - 1))
            cc_count += 1
    avg_cc = cc_sum / cc_count if cc_count > 0 else 0
    print(f"  Avg clustering: {avg_cc:.6f}")

    # Forman-Ricci curvature
    forman_curvs = []
    for i in range(V):
        for j in range(i+1, V):
            if A[i][j] > 0:
                di, dj = int(deg[i]), int(deg[j])
                tri_ij = sum(1 for k in range(V) if k != i and k != j
                             and A[i][k] > 0 and A[j][k] > 0)
                forman_curvs.append(4 - di - dj + 3 * tri_ij)
    if forman_curvs:
        print(f"  Forman-Ricci: min={min(forman_curvs)}, max={max(forman_curvs)}, "
              f"avg={np.mean(forman_curvs):.4f}, total={sum(forman_curvs)}")

    # Degree entropy
    deg_counts = Counter(deg)
    probs = np.array(list(deg_counts.values()), dtype=float) / V
    entropy = -np.sum(probs * np.log2(probs + 1e-30))
    print(f"  Degree entropy: {entropy:.4f} bits")

    # Betweenness centrality
    betweenness = np.zeros(V)
    for s in range(V):
        stack = []
        pred = [[] for _ in range(V)]
        sigma = np.zeros(V)
        sigma[s] = 1
        d_arr = np.full(V, -1)
        d_arr[s] = 0
        queue = [s]
        while queue:
            v = queue.pop(0)
            stack.append(v)
            for w in range(V):
                if A[v][w] > 0:
                    if d_arr[w] < 0:
                        d_arr[w] = d_arr[v] + 1
                        queue.append(w)
                    if d_arr[w] == d_arr[v] + 1:
                        sigma[w] += sigma[v]
                        pred[w].append(v)
        delta = np.zeros(V)
        while stack:
            w = stack.pop()
            for v in pred[w]:
                delta[v] += (sigma[v] / sigma[w]) * (1 + delta[w])
            if w != s:
                betweenness[w] += delta[w]
    betweenness /= 2
    print(f"  Betweenness: max={np.max(betweenness):.4f}, avg={np.mean(betweenness):.4f}")

    # Ramanujan
    d_avg = np.mean(deg)
    threshold = 2 * np.sqrt(d_avg - 1) if d_avg > 1 else 0
    second_eig = max(abs(eigenvalues[1]), abs(eigenvalues[-1]))
    ramanujan = second_eig <= threshold
    ratio = float(second_eig / threshold) if threshold > 0 else float('inf')
    print(f"  Ramanujan: {ramanujan} (ratio={ratio:.4f})")

    # H-gradient analysis (if we have H data)
    if classes_data:
        H_vals = np.array([c['H'] for c in classes_data])
        uphill = 0
        downhill = 0
        level_count = 0
        for i in range(V):
            for j in range(i+1, V):
                if A[i][j] > 0:
                    if H_vals[i] < H_vals[j]:
                        uphill += 1
                    elif H_vals[i] > H_vals[j]:
                        downhill += 1
                    else:
                        level_count += 1
        print(f"  H-gradient: uphill={uphill}, downhill={downhill}, level={level_count}")
        print(f"  DAG: {level_count == 0}")
        print(f"  H levels: {len(set(H_vals))}")
        print(f"  H range: [{min(H_vals)}, {max(H_vals)}]")

    return {
        'V': V, 'E': E_count, 'diameter': diameter, 'wiener': wiener,
        'triangles': triangles, 'spectral_radius': float(eigenvalues[0]),
        'alg_conn': alg_conn, 'graph_energy': float(np.sum(np.abs(eigenvalues))),
        'kirchhoff': kirchhoff if 'kirchhoff' in dir() else None,
        'avg_cc': avg_cc, 'avg_dist': avg_dist, 'entropy': entropy,
        'ramanujan': ramanujan, 'eigenvalues': eigenvalues
    }


# Compute invariants for G_7
print(f"\n  Computing G_7 invariants...")
res_g7 = compute_invariants(A_orig, V, "G_7", class_list)

# Compute invariants for G_7/Z_2
print(f"\n  Computing G_7/Z_2 invariants...")
res_g7m = compute_invariants(A_merged, V_merged, "G_7/Z_2", merged_classes)

# ============================================================================
# PHASE 7: Blue/black subgraph analysis
# ============================================================================

print(f"\n  BLUE/BLACK SUBGRAPH ANALYSIS OF G_7/Z_2")
print(f"  {'-'*50}")

A_blue = np.zeros((V_merged, V_merged), dtype=int)
A_black = np.zeros((V_merged, V_merged), dtype=int)
for e, color in merged_edge_color.items():
    a, b = e
    if color == 'blue':
        A_blue[a][b] = 1; A_blue[b][a] = 1
    elif color == 'black':
        A_black[a][b] = 1; A_black[b][a] = 1
    else:
        A_blue[a][b] = 1; A_blue[b][a] = 1  # mixed -> count as both

# Connected components
def cc(A_sub, V):
    visited = [False]*V
    comps = []
    for s in range(V):
        if visited[s]: continue
        comp = []
        q = [s]
        visited[s] = True
        while q:
            v = q.pop(0)
            comp.append(v)
            for u in range(V):
                if A_sub[v][u] > 0 and not visited[u]:
                    visited[u] = True
                    q.append(u)
        comps.append(comp)
    return comps

blue_comps = cc(A_blue, V_merged)
black_comps = cc(A_black, V_merged)
print(f"  Blue subgraph: {len(blue_comps)} components, sizes {sorted([len(c) for c in blue_comps], reverse=True)[:10]}")
print(f"  Black subgraph: {len(black_comps)} components, sizes {sorted([len(c) for c in black_comps], reverse=True)[:10]}")

blue_degs = np.sum(A_blue, axis=1)
black_degs = np.sum(A_black, axis=1)
print(f"  Blue degree: min={np.min(blue_degs)}, max={np.max(blue_degs)}, avg={np.mean(blue_degs):.2f}")
print(f"  Black degree: min={np.min(black_degs)}, max={np.max(black_degs)}, avg={np.mean(black_degs):.2f}")

# SC vs NS breakdown
sc_mids = [c['mid'] for c in merged_classes if c['sc']]
ns_mids = [c['mid'] for c in merged_classes if not c['sc']]
print(f"  SC merged vertices: {len(sc_mids)}, NS merged vertices: {len(ns_mids)}")

if sc_mids:
    sc_blue = [blue_degs[m] for m in sc_mids]
    sc_black = [black_degs[m] for m in sc_mids]
    print(f"  SC blue degree: mean={np.mean(sc_blue):.2f}, range=[{min(sc_blue)},{max(sc_blue)}]")
    print(f"  SC black degree: mean={np.mean(sc_black):.2f}, range=[{min(sc_black)},{max(sc_black)}]")
if ns_mids:
    ns_blue = [blue_degs[m] for m in ns_mids]
    ns_black = [black_degs[m] for m in ns_mids]
    print(f"  NS blue degree: mean={np.mean(ns_blue):.2f}, range=[{min(ns_blue)},{max(ns_blue)}]")
    print(f"  NS black degree: mean={np.mean(ns_black):.2f}, range=[{min(ns_black)},{max(ns_black)}]")

# ============================================================================
# PHASE 8: DEGREE-H CORRELATION
# ============================================================================

print(f"\n  DEGREE-H CORRELATION ON G_7/Z_2")
print(f"  {'-'*50}")

H_vals = np.array([c['H'] for c in merged_classes])
c3_vals = np.array([c['c3'] for c in merged_classes])
aut_vals = np.array([c['aut'] for c in merged_classes])
deg_vals = np.sum(A_merged, axis=1).astype(float)

if np.std(H_vals) > 0:
    print(f"  corr(H, degree) = {np.corrcoef(H_vals, deg_vals)[0,1]:.4f}")
if np.std(c3_vals) > 0:
    print(f"  corr(c3, degree) = {np.corrcoef(c3_vals, deg_vals)[0,1]:.4f}")
    print(f"  corr(c3, H) = {np.corrcoef(c3_vals, H_vals)[0,1]:.4f}")
if np.std(aut_vals) > 0:
    print(f"  corr(|Aut|, degree) = {np.corrcoef(aut_vals, deg_vals)[0,1]:.4f}")

# Top vertices by deg*H
importance = deg_vals * H_vals
ranked = sorted(range(V_merged), key=lambda i: importance[i], reverse=True)
print(f"  Top 10 by deg*H:")
for rank, idx in enumerate(ranked[:10]):
    c = merged_classes[idx]
    print(f"    {rank+1}. mid={idx} H={c['H']} deg={int(deg_vals[idx])} c3={c['c3']} sc={c['sc']} aut={c['aut']}")

# ============================================================================
# SUMMARY
# ============================================================================

print(f"\n\n  {'='*60}")
print(f"  SUMMARY FOR n=7")
print(f"  {'='*60}")
print(f"  G_7: V={V}, E={E}, diameter={res_g7['diameter']}, triangles={res_g7['triangles']}")
print(f"  G_7/Z_2: V={V_merged}, E={E_merged}, diameter={res_g7m['diameter']}, triangles={res_g7m['triangles']}")
print(f"  SC={sc_count}, collapsed={collapsed}, twin={twin}")
print(f"  Blue merged={blue_m}, Black merged={black_m}, Mixed={mixed_m}")
print(f"  Spectral radius G_7: {res_g7['spectral_radius']:.4f}")
print(f"  Spectral radius G_7/Z_2: {res_g7m['spectral_radius']:.4f}")
print(f"  Algebraic connectivity G_7: {res_g7['alg_conn']:.4f}")
print(f"  Algebraic connectivity G_7/Z_2: {res_g7m['alg_conn']:.4f}")
print(f"  Graph energy G_7: {res_g7['graph_energy']:.4f}")
print(f"  Graph energy G_7/Z_2: {res_g7m['graph_energy']:.4f}")

print(f"\n  Total runtime: {time.time()-t0:.0f}s")
print("=" * 80)
