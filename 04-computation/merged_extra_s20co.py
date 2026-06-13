#!/usr/bin/env python3
"""
merged_extra_s20co.py — Extra computations for G_n and G_n/Z_2
kind-pasteur-2026-03-23-S20co

Computes:
  1. Independence polynomial of G_6 and G_6/Z_2 (V=56 and V=34, small ind. number)
  2. Clique complex Betti numbers of G_6 and G_6/Z_2
  3. Full sequence catalog with OEIS checks
  4. H-gradient DAG analysis on ORIGINAL G_n (not merged)
  5. Degree-H correlation analysis
  6. Edge weight distribution analysis on merged graph
  7. Connected components of blue/black subgraphs in merged graph
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  EXTRA COMPUTATIONS FOR MERGED META-GRAPH")
print("  kind-pasteur-2026-03-23-S20co")
print("=" * 80)

# ============================================================================
# TOURNAMENT HELPERS (copied from main script)
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
    for v in range(n):
        dp[(1 << v) * n + v] = 1
    for S in range(1, 1 << n):
        pc = bin(S).count('1')
        if pc >= n:
            continue
        for v in range(n):
            if not (S & (1 << v)):
                continue
            val = dp[S * n + v]
            if val == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if adj[v][u]:
                    dp[(S | (1 << u)) * n + u] += val
    full = (1 << n) - 1
    return sum(dp[full * n + v] for v in range(n))

def canonical_form(adj, n):
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

def complement_adj(adj, n):
    return [[1 - adj[i][j] if i != j else 0 for j in range(n)] for i in range(n)]

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

def aut_size(adj, n):
    count = 0
    scores = [sum(adj[i][j] for j in range(n)) for i in range(n)]
    score_groups = defaultdict(list)
    for v in range(n):
        score_groups[scores[v]].append(v)
    sorted_scores = sorted(set(scores))
    groups = [score_groups[s] for s in sorted_scores]
    def gen_perms(groups):
        if not groups:
            yield []
            return
        for p in permutations(groups[0]):
            for rest in gen_perms(groups[1:]):
                yield list(p) + rest
    for perm in gen_perms(groups):
        ok = True
        for i in range(n):
            for j in range(i+1, n):
                if adj[perm[i]][perm[j]] != adj[i][j]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            count += 1
    return count

# ============================================================================
# BUILD G_n for n=3..6
# ============================================================================

def build_all(n_val):
    n = n_val
    m = comb(n, 2)
    total = 1 << m

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
        aut = aut_size(adj, n)
        class_list.append({
            'cid': idx, 'canon': canon, 'members': members,
            'adj': adj, 'H': h, 'sc': sc, 'score': ss,
            'c3': c3, 'aut': aut, 'size': len(members), 'bits': members[0]
        })
        canon_to_cid[canon] = idx

    for data in class_list:
        comp = complement_adj(data['adj'], n)
        comp_canon = canonical_form(comp, n)
        data['comp_cid'] = canon_to_cid.get(comp_canon, -1)

    # Build edges
    edges = set()
    self_loops = set()
    edge_colors = {}

    for data in class_list:
        cid = data['cid']
        for bits in data['members']:
            adj2 = tournament_adj(n, bits)
            for arc_idx in range(m):
                flipped = bits ^ (1 << arc_idx)
                fadj = tournament_adj(n, flipped)
                fcanon = canonical_form(fadj, n)
                nb = canon_to_cid.get(fcanon)
                if nb is None:
                    continue
                if nb == cid:
                    self_loops.add(cid)
                else:
                    e = (min(cid, nb), max(cid, nb))
                    edges.add(e)
                    if e not in edge_colors:
                        sc_i = class_list[cid]['sc']
                        sc_j = class_list[nb]['sc']
                        edge_colors[e] = 'blue' if sc_i == sc_j else 'black'

    V = len(class_list)
    sc_count = sum(1 for d in class_list if d['sc'])

    # Build adjacency matrix
    A = np.zeros((V, V), dtype=int)
    for (a, b) in edges:
        A[a][b] = 1
        A[b][a] = 1

    # Build merged
    merged_id = {}
    mid = 0
    for data in class_list:
        cid = data['cid']
        if cid in merged_id:
            continue
        comp = data['comp_cid']
        merged_id[cid] = mid
        if comp != cid:
            merged_id[comp] = mid
        mid += 1

    V_merged = mid
    merged_edges = set()
    collapsed = 0
    for (a, b) in edges:
        ma, mb = merged_id[a], merged_id[b]
        if ma == mb:
            collapsed += 1
        else:
            merged_edges.add((min(ma, mb), max(ma, mb)))

    A_merged = np.zeros((V_merged, V_merged), dtype=int)
    for (a, b) in merged_edges:
        A_merged[a][b] = 1
        A_merged[b][a] = 1

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

    return {
        'n': n, 'V': V, 'E': len(edges), 'classes': class_list,
        'edges': edges, 'edge_colors': edge_colors, 'self_loops': self_loops,
        'sc_count': sc_count, 'A': A,
        'V_merged': V_merged, 'E_merged': len(merged_edges),
        'merged_edges': merged_edges, 'A_merged': A_merged,
        'merged_classes': merged_classes, 'merged_id': merged_id,
        'collapsed': collapsed
    }


# ============================================================================
# PART 1: INDEPENDENCE POLYNOMIAL (optimized for small ind. number)
# ============================================================================

def independence_polynomial(A, V, max_check=None):
    """Compute full independence polynomial. Stop at alpha."""
    alpha = [0] * (V + 1)
    alpha[0] = 1
    for size in range(1, V + 1):
        count = 0
        total_combos = comb(V, size)
        if max_check and total_combos > max_check:
            # Use backtracking instead
            count = count_independent_sets_backtrack(A, V, size)
        else:
            for subset in combinations(range(V), size):
                if all(A[subset[i]][subset[j]] == 0
                       for i in range(len(subset)) for j in range(i+1, len(subset))):
                    count += 1
        alpha[size] = count
        if count == 0:
            break
    while alpha and alpha[-1] == 0:
        alpha.pop()
    return alpha


def count_independent_sets_backtrack(A, V, target_size):
    """Count independent sets of exact size target_size using backtracking."""
    count = 0
    def backtrack(start, current, excluded):
        nonlocal count
        if len(current) == target_size:
            count += 1
            return
        remaining = target_size - len(current)
        for v in range(start, V):
            if v in excluded:
                continue
            if V - v < remaining:
                break
            new_excluded = excluded | {u for u in range(V) if A[v][u] > 0}
            backtrack(v + 1, current + [v], new_excluded)
    backtrack(0, [], set())
    return count


# ============================================================================
# PART 2: CLIQUE COMPLEX BETTI NUMBERS
# ============================================================================

def clique_complex_betti(A, V, max_dim=None):
    """Compute Betti numbers of the clique complex."""
    cliques_by_dim = defaultdict(list)
    cliques_by_dim[0] = [(v,) for v in range(V)]

    for size in range(2, V + 1):
        if max_dim and size - 1 > max_dim:
            break
        found_any = False
        if comb(V, size) > 500000:
            break
        for subset in combinations(range(V), size):
            if all(A[subset[i]][subset[j]] > 0
                   for i in range(len(subset)) for j in range(i+1, len(subset))):
                cliques_by_dim[size - 1].append(subset)
                found_any = True
        if not found_any:
            break

    max_d = max(cliques_by_dim.keys()) if cliques_by_dim else 0
    f_vector = [len(cliques_by_dim[d]) for d in range(max_d + 1)]
    print(f"    Clique complex f-vector: {f_vector}")

    betti = []
    for d in range(max_d + 1):
        if d == 0:
            if 1 in cliques_by_dim and len(cliques_by_dim[1]) > 0:
                edges_list = cliques_by_dim[1]
                vertices_list = cliques_by_dim[0]
                n_v = len(vertices_list)
                n_e = len(edges_list)
                bd = np.zeros((n_v, n_e), dtype=float)
                v_index = {v: i for i, v in enumerate(vertices_list)}
                for j, (u, w) in enumerate(edges_list):
                    bd[v_index[(u,)]][j] = 1
                    bd[v_index[(w,)]][j] = -1
                rank_d1 = np.linalg.matrix_rank(bd)
                betti.append(n_v - rank_d1)
            else:
                betti.append(len(cliques_by_dim[0]))
        else:
            simplices_d = cliques_by_dim.get(d, [])
            simplices_dm1 = cliques_by_dim.get(d-1, [])
            simplices_dp1 = cliques_by_dim.get(d+1, [])

            if len(simplices_d) == 0:
                betti.append(0)
                continue

            s_dm1_index = {s: i for i, s in enumerate(simplices_dm1)}
            bd_d = np.zeros((len(simplices_dm1), len(simplices_d)), dtype=float)
            for j, sigma in enumerate(simplices_d):
                for k in range(len(sigma)):
                    face = sigma[:k] + sigma[k+1:]
                    if face in s_dm1_index:
                        bd_d[s_dm1_index[face]][j] = (-1)**k

            ker_d = len(simplices_d) - np.linalg.matrix_rank(bd_d)

            if len(simplices_dp1) > 0:
                s_d_index = {s: i for i, s in enumerate(simplices_d)}
                bd_dp1 = np.zeros((len(simplices_d), len(simplices_dp1)), dtype=float)
                for j, sigma in enumerate(simplices_dp1):
                    for k in range(len(sigma)):
                        face = sigma[:k] + sigma[k+1:]
                        if face in s_d_index:
                            bd_dp1[s_d_index[face]][j] = (-1)**k
                im_dp1 = np.linalg.matrix_rank(bd_dp1)
            else:
                im_dp1 = 0

            betti.append(int(ker_d - im_dp1))

    return betti, f_vector


# ============================================================================
# PART 3: H-GRADIENT DAG ANALYSIS ON ORIGINAL G_n
# ============================================================================

def h_gradient_analysis(data):
    """Detailed H-gradient analysis on G_n (original, not merged)."""
    n = data['n']
    classes = data['classes']
    edges = data['edges']

    print(f"\n  H-GRADIENT ANALYSIS OF G_{n} (ORIGINAL)")
    print(f"  {'-'*50}")

    uphill = 0
    downhill = 0
    level = 0
    level_edges = []

    for (a, b) in edges:
        ha = classes[a]['H']
        hb = classes[b]['H']
        if ha < hb:
            uphill += 1
        elif ha > hb:
            downhill += 1
        else:
            level += 1
            level_edges.append((a, b, ha))

    # By symmetry, if (T, T') is an edge, (T^op, T'^op) is also an edge
    # So every uphill edge has a "mirror" downhill edge via complement
    print(f"  Uphill: {uphill}, Downhill: {downhill}, Level: {level}")
    print(f"  DAG: {level == 0}")
    if level_edges:
        print(f"  Level edges:")
        for (a, b, h) in level_edges:
            sa = classes[a]['sc']
            sb = classes[b]['sc']
            print(f"    ({a},{b}) at H={h} | SC=({sa},{sb}) | scores=({classes[a]['score']},{classes[b]['score']})")

    # Check if DAG ignoring complements (orient by H, break ties by SC)
    # Number of distinct H levels
    h_vals = sorted(set(c['H'] for c in classes))
    print(f"  H levels: {len(h_vals)}: {h_vals}")

    # Longest monotone path
    # Build DAG from uphill edges
    adj_dag = defaultdict(set)
    for (a, b) in edges:
        ha = classes[a]['H']
        hb = classes[b]['H']
        if ha < hb:
            adj_dag[a].add(b)
        elif ha > hb:
            adj_dag[b].add(a)

    # Longest path in DAG via DP
    topo_order = sorted(range(len(classes)), key=lambda i: classes[i]['H'])
    longest = [1] * len(classes)
    for v in topo_order:
        for u in adj_dag[v]:
            longest[u] = max(longest[u], longest[v] + 1)
    max_chain = max(longest)
    print(f"  Longest H-increasing chain: {max_chain}")

    # Width (max antichain at any H level in DAG)
    max_width = max(Counter(c['H'] for c in classes).values())
    print(f"  Max width (vertices at same H): {max_width}")

    # Count H-monotone Hamiltonian paths (how many ways to visit all classes in H-order?)
    # This is expensive, skip for large V
    V = len(classes)
    if V <= 15:
        mono_ham_count = count_monotone_ham_paths(adj_dag, V, topo_order, classes)
        print(f"  H-monotone Hamiltonian paths: {mono_ham_count}")

    return {
        'uphill': uphill, 'downhill': downhill, 'level': level,
        'is_dag': level == 0, 'max_chain': max_chain, 'max_width': max_width,
        'h_levels': len(h_vals)
    }


def count_monotone_ham_paths(adj_dag, V, topo_order, classes):
    """Count Hamiltonian paths that strictly increase H (allowing level edges)."""
    # Build full adjacency (both uphill and level edges)
    adj_mono = defaultdict(set)
    for v in range(V):
        for u in adj_dag[v]:
            adj_mono[v].add(u)

    # DP over subsets
    dp = {}
    for v in range(V):
        dp[(1 << v, v)] = 1
    for S in range(1, 1 << V):
        for v in range(V):
            if not (S & (1 << v)):
                continue
            val = dp.get((S, v), 0)
            if val == 0:
                continue
            for u in adj_mono[v]:
                if S & (1 << u):
                    continue
                key = (S | (1 << u), u)
                dp[key] = dp.get(key, 0) + val
    full = (1 << V) - 1
    return sum(dp.get((full, v), 0) for v in range(V))


# ============================================================================
# PART 4: BLUE/BLACK SUBGRAPH ANALYSIS
# ============================================================================

def blue_black_subgraph_analysis(data):
    """Analyze blue and black subgraphs of the merged graph separately."""
    n = data['n']
    V = data['V_merged']
    edges = data['merged_edges']

    # Determine edge color in merged graph
    merged_edge_colors = {}
    merged_id = data['merged_id']
    for (a, b) in data['edges']:
        ma, mb = merged_id[a], merged_id[b]
        if ma == mb:
            continue
        e = (min(ma, mb), max(ma, mb))
        color = data['edge_colors'].get((min(a,b), max(a,b)), 'blue')
        if e not in merged_edge_colors:
            merged_edge_colors[e] = color
        elif merged_edge_colors[e] != color:
            merged_edge_colors[e] = 'mixed'

    blue_edges = [e for e, c in merged_edge_colors.items() if c == 'blue']
    black_edges = [e for e, c in merged_edge_colors.items() if c == 'black']
    mixed_edges = [e for e, c in merged_edge_colors.items() if c == 'mixed']

    print(f"\n  BLUE/BLACK SUBGRAPH ANALYSIS OF G_{n}/Z_2")
    print(f"  {'-'*50}")
    print(f"  Blue edges: {len(blue_edges)}, Black edges: {len(black_edges)}, Mixed: {len(mixed_edges)}")

    # Blue subgraph
    A_blue = np.zeros((V, V), dtype=int)
    for (a, b) in blue_edges:
        A_blue[a][b] = 1
        A_blue[b][a] = 1

    # Black subgraph
    A_black = np.zeros((V, V), dtype=int)
    for (a, b) in black_edges:
        A_black[a][b] = 1
        A_black[b][a] = 1

    # Connected components
    def connected_components(A_sub, V):
        visited = [False] * V
        components = []
        for start in range(V):
            if visited[start]:
                continue
            comp = []
            queue = [start]
            visited[start] = True
            while queue:
                v = queue.pop(0)
                comp.append(v)
                for u in range(V):
                    if A_sub[v][u] > 0 and not visited[u]:
                        visited[u] = True
                        queue.append(u)
            components.append(comp)
        return components

    blue_comps = connected_components(A_blue, V)
    black_comps = connected_components(A_black, V)

    print(f"  Blue subgraph: {len(blue_comps)} components, sizes {sorted([len(c) for c in blue_comps], reverse=True)}")
    print(f"  Black subgraph: {len(black_comps)} components, sizes {sorted([len(c) for c in black_comps], reverse=True)}")

    # Spectral radius of each subgraph
    if len(blue_edges) > 0:
        eigs_blue = np.sort(np.linalg.eigvalsh(A_blue.astype(float)))[::-1]
        print(f"  Blue spectral radius: {eigs_blue[0]:.4f}")
    if len(black_edges) > 0:
        eigs_black = np.sort(np.linalg.eigvalsh(A_black.astype(float)))[::-1]
        print(f"  Black spectral radius: {eigs_black[0]:.4f}")

    # Blue density, black density
    blue_density = 2 * len(blue_edges) / (V * (V-1)) if V > 1 else 0
    black_density = 2 * len(black_edges) / (V * (V-1)) if V > 1 else 0
    print(f"  Blue density: {blue_density:.6f}, Black density: {black_density:.6f}")

    # SC-type analysis of merged vertices
    merged_classes = data['merged_classes']
    sc_verts = [c['mid'] for c in merged_classes if c['sc']]
    ns_verts = [c['mid'] for c in merged_classes if not c['sc']]
    print(f"  SC merged vertices: {len(sc_verts)}, NS merged vertices: {len(ns_verts)}")

    # Blue-degree distribution for SC vs NS
    blue_deg_sc = [sum(A_blue[v]) for v in sc_verts]
    blue_deg_ns = [sum(A_blue[v]) for v in ns_verts]
    black_deg_sc = [sum(A_black[v]) for v in sc_verts]
    black_deg_ns = [sum(A_black[v]) for v in ns_verts]

    if blue_deg_sc:
        print(f"  SC blue degree: mean={np.mean(blue_deg_sc):.2f}, range=[{min(blue_deg_sc)},{max(blue_deg_sc)}]")
    if blue_deg_ns:
        print(f"  NS blue degree: mean={np.mean(blue_deg_ns):.2f}, range=[{min(blue_deg_ns)},{max(blue_deg_ns)}]")
    if black_deg_sc:
        print(f"  SC black degree: mean={np.mean(black_deg_sc):.2f}, range=[{min(black_deg_sc)},{max(black_deg_sc)}]")
    if black_deg_ns:
        print(f"  NS black degree: mean={np.mean(black_deg_ns):.2f}, range=[{min(black_deg_ns)},{max(black_deg_ns)}]")


# ============================================================================
# PART 5: DEGREE-H CORRELATION AND CENTRALITY ANALYSIS
# ============================================================================

def centrality_analysis(data):
    """Analyze correlation between tournament invariants and graph-theoretic centrality."""
    n = data['n']
    V = data['V_merged']
    A = data['A_merged']
    classes = data['merged_classes']

    print(f"\n  CENTRALITY-INVARIANT CORRELATION ON G_{n}/Z_2")
    print(f"  {'-'*50}")

    H_vals = np.array([c['H'] for c in classes])
    c3_vals = np.array([c['c3'] for c in classes])
    aut_vals = np.array([c['aut'] for c in classes])
    deg_vals = np.sum(A, axis=1).astype(float)

    # Correlations
    if V > 2 and np.std(H_vals) > 0:
        corr_H_deg = np.corrcoef(H_vals, deg_vals)[0, 1]
        print(f"  corr(H, degree) = {corr_H_deg:.4f}")
    if V > 2 and np.std(c3_vals) > 0:
        corr_c3_deg = np.corrcoef(c3_vals, deg_vals)[0, 1]
        corr_c3_H = np.corrcoef(c3_vals, H_vals)[0, 1]
        print(f"  corr(c3, degree) = {corr_c3_deg:.4f}")
        print(f"  corr(c3, H) = {corr_c3_H:.4f}")
    if V > 2 and np.std(aut_vals) > 0:
        corr_aut_deg = np.corrcoef(aut_vals, deg_vals)[0, 1]
        print(f"  corr(|Aut|, degree) = {corr_aut_deg:.4f}")

    # Eigenvector centrality
    eigenvalues, eigenvectors = np.linalg.eigh(A.astype(float))
    idx = np.argmax(eigenvalues)
    evec = np.abs(eigenvectors[:, idx])
    evec /= np.max(evec)

    if V > 2 and np.std(evec) > 0:
        corr_evec_H = np.corrcoef(evec, H_vals)[0, 1]
        print(f"  corr(eigenvector_centrality, H) = {corr_evec_H:.4f}")

    # PageRank-style: which vertices are most important?
    # (simplified: use degree * H as importance)
    importance = deg_vals * H_vals
    ranked = sorted(range(V), key=lambda i: importance[i], reverse=True)
    print(f"  Top 5 by deg*H:")
    for rank, idx in enumerate(ranked[:5]):
        c = classes[idx]
        print(f"    {rank+1}. mid={idx} H={c['H']} deg={int(deg_vals[idx])} c3={c['c3']} sc={c['sc']} score={c['score']}")


# ============================================================================
# MAIN
# ============================================================================

all_data = {}
for n in [3, 4, 5, 6]:
    print(f"\n{'#'*80}")
    print(f"  n = {n}")
    print(f"{'#'*80}")
    t0 = time.time()
    data = build_all(n)
    print(f"  Built in {time.time()-t0:.1f}s")
    all_data[n] = data

    # Independence polynomial of G_n
    print(f"\n  --- Independence Polynomial of G_{n} ---")
    ip_gn = independence_polynomial(data['A'], data['V'], max_check=500000)
    I2_gn = sum(a * (2**k) for k, a in enumerate(ip_gn))
    Ineg1_gn = sum(a * ((-1)**k) for k, a in enumerate(ip_gn))
    print(f"  I(G_{n}, x) = {ip_gn}")
    print(f"  I(G_{n}, 2) = {I2_gn} (meta-H)")
    print(f"  I(G_{n}, -1) = {Ineg1_gn} (Euler char)")

    # Independence polynomial of G_n/Z_2
    print(f"\n  --- Independence Polynomial of G_{n}/Z_2 ---")
    ip_merged = independence_polynomial(data['A_merged'], data['V_merged'], max_check=500000)
    I2_m = sum(a * (2**k) for k, a in enumerate(ip_merged))
    Ineg1_m = sum(a * ((-1)**k) for k, a in enumerate(ip_merged))
    print(f"  I(G_{n}/Z_2, x) = {ip_merged}")
    print(f"  I(G_{n}/Z_2, 2) = {I2_m} (meta-H)")
    print(f"  I(G_{n}/Z_2, -1) = {Ineg1_m} (Euler char)")

    # Betti numbers of G_n
    if data['V'] <= 56:
        print(f"\n  --- Clique Complex Betti Numbers of G_{n} ---")
        betti_gn, fvec_gn = clique_complex_betti(data['A'], data['V'], max_dim=5)
        euler_gn = sum((-1)**k * b for k, b in enumerate(betti_gn))
        print(f"  beta(G_{n}) = {betti_gn}")
        print(f"  chi(G_{n}) = {euler_gn}")

    # Betti numbers of G_n/Z_2
    print(f"\n  --- Clique Complex Betti Numbers of G_{n}/Z_2 ---")
    betti_m, fvec_m = clique_complex_betti(data['A_merged'], data['V_merged'], max_dim=5)
    euler_m = sum((-1)**k * b for k, b in enumerate(betti_m))
    print(f"  beta(G_{n}/Z_2) = {betti_m}")
    print(f"  chi(G_{n}/Z_2) = {euler_m}")

    # H-gradient analysis
    h_gradient_analysis(data)

    # Blue/black subgraph analysis
    blue_black_subgraph_analysis(data)

    # Centrality analysis
    centrality_analysis(data)


# ============================================================================
# SEQUENCE CATALOG
# ============================================================================

print(f"\n\n{'='*80}")
print("  COMPLETE SEQUENCE CATALOG")
print(f"{'='*80}")

print("""
NEW SEQUENCES (first computed here, check against OEIS):

=== G_n sequences (n=3,4,5,6) ===
V(G_n) = 2, 4, 12, 56                  [A000568 - known]
E(G_n) = 1, 5, 30, 290                 [check OEIS]
self_loops(G_n) = 1, 2, 7, 30          [check OEIS]
triangles(G_n) = 0, 2, 21, 248         [check OEIS]
diameter(G_n) = 1, 2, 3, 4             [check OEIS - matches n-2!]
girth(G_n) = inf, 3, 3, 3              [trivial]
spanning_trees(G_n) = 1, 8, 2347680, ? [check OEIS]
Wiener(G_n) = 1, 7, 108, 3314          [check OEIS]
I(G_n, 2) = 5, 13, 793, ?              [NEW meta-H]
I(G_n, -1) = -1, -2, 1, ?              [NEW Euler char]

=== G_n/Z_2 sequences (n=3,4,5,6) ===
V(G_n/Z_2) = 2, 3, 10, 34             [= (A000568 + SC)/2]
E(G_n/Z_2) = 1, 3, 21, 143            [check OEIS]
collapsed = 0, 0, 0, 5                 [check OEIS]
twin = 0, 2, 9, 142                    [check OEIS]
triangles(G_n/Z_2) = 0, 1, 12, 139    [check OEIS]
diameter(G_n/Z_2) = 1, 1, 3, 4        [check OEIS]
spanning_trees(G_n/Z_2) = 1, 3, 32159, ?  [check OEIS]
Wiener(G_n/Z_2) = 1, 3, 73, 1138      [check OEIS]
I(G_n/Z_2, 2) = 5, 7, 301, ?          [NEW merged meta-H]
I(G_n/Z_2, -1) = -1, -2, 1, ?         [NEW Euler char]
blue_merged = 1, 1, 13, 98             [check OEIS]
black_merged = 0, 2, 8, 45             [check OEIS]

=== Spectral sequences ===
spectral_radius(G_n) = 1.00, 2.56, 5.58, 11.67
spectral_radius(G_n/Z_2) = 1.00, 2.00, 4.64, 9.78
alg_conn(G_n) = 2.00, 2.00, 1.60, 1.96
alg_conn(G_n/Z_2) = 2.00, 3.00, 1.43, 1.47
graph_energy(G_n) = 2.0, 5.1, 20.2, 137.4
graph_energy(G_n/Z_2) = 2.0, 4.0, 15.9, 73.2
Estrada(G_n) = 3.1, 14.5, 282.8, 118439
Estrada(G_n/Z_2) = 3.1, 8.1, 116.7, 18037

=== Topological sequences ===
beta_0 always = 1 (connected)
beta_1(G_n) = 0, 0, 2, ?
beta_1(G_n/Z_2) = 0, 0, 2, ?
chi(clique_cmplx G_n) = 1, 1, -1, ?
chi(clique_cmplx G_n/Z_2) = 1, 1, -1, ?

=== Curvature sequences ===
Forman_avg(G_n) = 2.0, 2.4, -0.97, -11.18
Forman_avg(G_n/Z_2) = 2.0, 3.0, -0.19, -6.45
Forman_total(G_n) = 2, 12, -29, -3242
Forman_total(G_n/Z_2) = 2, 9, -4, -923

=== Key findings ===
1. diameter(G_n) = n-2 CONFIRMED for n=3..6 (CONJECTURE from earlier)
2. G_n/Z_2 is NOT a DAG at n>=5 (downhill edges appear)
3. G_n is NOT a DAG at n>=5 (but uphill >> downhill)
4. beta_1 = 2 appears at n=5 for BOTH G_n and G_n/Z_2
5. Forman curvature goes negative at n=5 (graph becomes hyperbolic)
6. Ramanujan property: G_4 YES, G_5 YES, G_6 NO (transition at n=6)
""")

print("\n  DONE.")
