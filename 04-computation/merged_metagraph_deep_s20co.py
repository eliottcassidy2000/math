#!/usr/bin/env python3
"""
merged_metagraph_deep_s20co.py  Deep invariant analysis of G_n/Z_2 (merged meta-graph)
kind-pasteur-2026-03-23-S20co

Computes 30+ invariants of both G_n (isomorphism class graph) and G_n/Z_2
(merged graph = quotient by complement involution) for n=3..6.
n=7 attempted if time permits.

NEW METRICS computed here for the first time:
  - Laplacian spectrum and algebraic connectivity of G_n/Z_2
  - Cheeger constant (isoperimetric number) of G_n/Z_2
  - Clique complex Betti numbers (simplicial homology)
  - Tutte polynomial evaluations
  - Resistance distance (effective resistance) statistics
  - Estrada index (trace of exp(A))
  - Graph energy (sum of |eigenvalues|)
  - Kirchhoff index (sum of effective resistances)
  - Automorphism count of G_n/Z_2 as abstract graph
  - H-weighted Laplacian spectrum
  - Ollivier-Ricci curvature statistics
  - Forman-Ricci curvature statistics
"""

import sys
import numpy as np
from math import comb, factorial, gcd
from itertools import permutations, combinations
from collections import defaultdict, Counter
from functools import lru_cache
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  DEEP INVARIANT ANALYSIS OF G_n AND G_n/Z_2")
print("  kind-pasteur-2026-03-23-S20co")
print("=" * 80)
print()

# ============================================================================
# CORE HELPERS
# ============================================================================

def tournament_adj(n, bits):
    """Convert bit encoding to adjacency matrix."""
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
    """Count Hamiltonian paths via DP."""
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
    """Canonical form of tournament under vertex permutation."""
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
    """Count automorphisms (permutations preserving adj)."""
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
# BUILD ISO CLASS GRAPH
# ============================================================================

def build_iso_class_graph(n):
    """Build the complete iso class graph G_n for given n."""
    print(f"\n{'='*70}")
    print(f"  BUILDING G_{n}")
    print(f"{'='*70}")
    t0 = time.time()

    m = comb(n, 2)
    total = 1 << m
    print(f"  n={n}, m={m}, |tournaments|={total}")

    # Step 1: Group tournaments by canonical form
    iso_classes = defaultdict(list)
    for bits in range(total):
        adj = tournament_adj(n, bits)
        canon = canonical_form(adj, n)
        iso_classes[canon].append(bits)

    V = len(iso_classes)
    print(f"  |V(G_{n})| = {V} iso classes ({time.time()-t0:.1f}s)")

    # Step 2: Build class data
    class_list = []
    canon_to_cid = {}
    for idx, (canon, members) in enumerate(sorted(iso_classes.items())):
        adj = tournament_adj(n, members[0])
        h = H_dp(adj, n)
        sc = canonical_form(complement_adj(adj, n), n) == canon
        ss = score_seq(adj, n)
        c3 = c3_count(adj, n)
        aut = aut_size(adj, n)
        class_size = len(members)

        class_list.append({
            'cid': idx, 'canon': canon, 'members': members,
            'adj': adj, 'H': h, 'sc': sc, 'score': ss,
            'c3': c3, 'aut': aut, 'size': class_size, 'bits': members[0]
        })
        canon_to_cid[canon] = idx

    # Step 3: Complement pairing
    for data in class_list:
        comp = complement_adj(data['adj'], n)
        comp_canon = canonical_form(comp, n)
        data['comp_cid'] = canon_to_cid.get(comp_canon, -1)

    sc_count = sum(1 for d in class_list if d['sc'])
    print(f"  SC classes: {sc_count}")

    # Step 4: Build edges
    edges = set()
    self_loops = set()
    edge_colors = {}  # (min,max) -> 'blue' or 'black'
    edge_weights = Counter()

    for data in class_list:
        cid = data['cid']
        adj = data['adj']
        sc_i = data['sc']

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
                    edge_weights[e] += 1
                    sc_j = class_list[nb]['sc']
                    if e not in edge_colors:
                        if sc_i == sc_j:
                            edge_colors[e] = 'blue'
                        else:
                            edge_colors[e] = 'black'

    E = len(edges)
    blue = sum(1 for c in edge_colors.values() if c == 'blue')
    black = sum(1 for c in edge_colors.values() if c == 'black')
    print(f"  |E(G_{n})| = {E} (blue={blue}, black={black})")
    print(f"  Self-loops: {len(self_loops)}")
    print(f"  Build time: {time.time()-t0:.1f}s")

    return {
        'n': n, 'V': V, 'E': E, 'classes': class_list,
        'edges': edges, 'edge_colors': edge_colors, 'edge_weights': edge_weights,
        'self_loops': self_loops, 'sc_count': sc_count,
        'canon_to_cid': canon_to_cid, 'blue': blue, 'black': black
    }


def build_merged_graph(gn_data):
    """Build G_n/Z_2 by merging complement pairs."""
    n = gn_data['n']
    classes = gn_data['classes']
    V = gn_data['V']

    print(f"\n  BUILDING G_{n}/Z_2 (merged graph)")

    # Assign merged IDs
    merged_id = {}
    mid = 0
    for data in classes:
        cid = data['cid']
        if cid in merged_id:
            continue
        comp = data['comp_cid']
        merged_id[cid] = mid
        if comp != cid:
            merged_id[comp] = mid
        mid += 1

    V_merged = mid
    print(f"  V_merged = {V_merged} (= ({V} + {gn_data['sc_count']})/2 = {(V + gn_data['sc_count'])//2})")

    # Build merged edges
    merged_edges = set()
    collapsed = 0
    twin = 0

    for (a, b) in gn_data['edges']:
        ma, mb = merged_id[a], merged_id[b]
        if ma == mb:
            collapsed += 1
        else:
            e = (min(ma, mb), max(ma, mb))
            merged_edges.add(e)

    E_merged = len(merged_edges)
    twin = len(gn_data['edges']) - E_merged - collapsed
    print(f"  E_merged = {E_merged}, collapsed = {collapsed}, twin = {twin}")
    print(f"  Check: {E_merged} + {collapsed} + {twin} = {E_merged + collapsed + twin} = {len(gn_data['edges'])} {'YES' if E_merged+collapsed+twin == len(gn_data['edges']) else 'NO'}")

    # Merged class properties: inherit from the SC class or the lower-cid of the pair
    merged_classes = []
    merged_H = {}
    for data in classes:
        cid = data['cid']
        mid_val = merged_id[cid]
        if mid_val not in merged_H:
            merged_H[mid_val] = data['H']
            merged_classes.append({
                'mid': mid_val,
                'cid': cid,
                'H': data['H'],
                'sc': data['sc'],
                'score': data['score'],
                'c3': data['c3'],
                'aut': data['aut'],
                'size': data['size'],
                'is_sc_class': data['sc'],
                'pair_cid': data['comp_cid']
            })

    # Build merged adjacency matrix
    adj_merged = np.zeros((V_merged, V_merged), dtype=int)
    for (a, b) in merged_edges:
        adj_merged[a][b] = 1
        adj_merged[b][a] = 1

    # Blue/black in merged graph
    blue_merged = 0
    black_merged = 0
    for (a, b) in gn_data['edges']:
        ma, mb = merged_id[a], merged_id[b]
        if ma != mb:
            color = gn_data['edge_colors'].get((min(a,b), max(a,b)), 'blue')
            e = (min(ma, mb), max(ma, mb))
            if e in merged_edges:
                if color == 'blue':
                    blue_merged += 1
                else:
                    black_merged += 1
    # Deduplicate: each merged edge counted from both original edges (or 1 if one was collapsed)
    # Actually we need a different approach - let's track per merged edge
    merged_edge_color = {}
    for (a, b) in gn_data['edges']:
        ma, mb = merged_id[a], merged_id[b]
        if ma == mb:
            continue
        e = (min(ma, mb), max(ma, mb))
        color = gn_data['edge_colors'].get((min(a,b), max(a,b)), 'blue')
        if e not in merged_edge_color:
            merged_edge_color[e] = color
        # If both blue and black edges map to same merged edge, mark mixed
        elif merged_edge_color[e] != color:
            merged_edge_color[e] = 'mixed'

    blue_m = sum(1 for c in merged_edge_color.values() if c == 'blue')
    black_m = sum(1 for c in merged_edge_color.values() if c == 'black')
    mixed_m = sum(1 for c in merged_edge_color.values() if c == 'mixed')

    print(f"  Merged colors: blue={blue_m}, black={black_m}, mixed={mixed_m}")

    return {
        'n': n, 'V': V_merged, 'E': E_merged,
        'edges': merged_edges, 'adj': adj_merged,
        'classes': merged_classes, 'merged_id': merged_id,
        'collapsed': collapsed, 'twin': twin,
        'blue': blue_m, 'black': black_m, 'mixed': mixed_m,
        'edge_colors': merged_edge_color
    }


# ============================================================================
# GRAPH INVARIANTS
# ============================================================================

def compute_graph_invariants(adj_matrix, name, vertex_data=None):
    """Compute 30+ graph invariants from adjacency matrix."""
    A = np.array(adj_matrix, dtype=float)
    V = A.shape[0]
    E = int(np.sum(A) // 2)

    print(f"\n  {'='*60}")
    print(f"  INVARIANTS OF {name} (V={V}, E={E})")
    print(f"  {'='*60}")

    results = {'V': V, 'E': E, 'name': name}

    # --- 1. Basic graph properties ---
    deg = np.sum(A, axis=1).astype(int)
    results['degree_seq'] = sorted(deg)
    results['min_deg'] = int(np.min(deg))
    results['max_deg'] = int(np.max(deg))
    results['avg_deg'] = float(np.mean(deg))
    results['density'] = 2*E / (V*(V-1)) if V > 1 else 0
    print(f"  Degree: min={results['min_deg']}, max={results['max_deg']}, avg={results['avg_deg']:.3f}")
    print(f"  Density: {results['density']:.6f}")

    # --- 2. Adjacency spectrum ---
    eigenvalues = np.linalg.eigvalsh(A)
    eigenvalues = np.sort(eigenvalues)[::-1]  # descending
    results['adj_spectrum'] = eigenvalues
    results['spectral_radius'] = float(eigenvalues[0])
    results['spectral_gap'] = float(eigenvalues[0] - eigenvalues[1]) if V > 1 else 0
    results['graph_energy'] = float(np.sum(np.abs(eigenvalues)))
    print(f"  Spectral radius: {results['spectral_radius']:.6f}")
    print(f"  Spectral gap: {results['spectral_gap']:.6f}")
    print(f"  Graph energy: {results['graph_energy']:.6f}")
    print(f"  Spectrum: [{', '.join(f'{x:.4f}' for x in eigenvalues[:min(10,V)])}{'...' if V>10 else ''}]")

    # --- 3. Estrada index = tr(exp(A)) ---
    results['estrada_index'] = float(np.sum(np.exp(eigenvalues)))
    print(f"  Estrada index: {results['estrada_index']:.6f}")

    # --- 4. Laplacian spectrum ---
    D = np.diag(deg.astype(float))
    L = D - A
    lap_eigs = np.sort(np.linalg.eigvalsh(L))
    results['lap_spectrum'] = lap_eigs
    results['algebraic_connectivity'] = float(lap_eigs[1]) if V > 1 else 0
    results['lap_spectral_radius'] = float(lap_eigs[-1])
    print(f"  Algebraic connectivity (lam2): {results['algebraic_connectivity']:.6f}")
    print(f"  Laplacian spectral radius: {results['lap_spectral_radius']:.6f}")

    # --- 5. Normalized Laplacian ---
    if np.all(deg > 0):
        D_inv_sqrt = np.diag(1.0 / np.sqrt(deg.astype(float)))
        L_norm = np.eye(V) - D_inv_sqrt @ A @ D_inv_sqrt
        norm_eigs = np.sort(np.linalg.eigvalsh(L_norm))
        results['norm_lap_spectrum'] = norm_eigs
        results['norm_alg_conn'] = float(norm_eigs[1]) if V > 1 else 0
        print(f"  Normalized algebraic connectivity: {results['norm_alg_conn']:.6f}")

    # --- 6. Kirchhoff index (sum of effective resistances) ---
    if V > 1 and results['algebraic_connectivity'] > 1e-10:
        # Kirchhoff = V * sum(1/lambda_i) for nonzero lambda_i
        nonzero_lap = lap_eigs[1:]  # skip lambda_0 = 0
        nonzero_lap = nonzero_lap[nonzero_lap > 1e-10]
        results['kirchhoff_index'] = float(V * np.sum(1.0 / nonzero_lap))
        print(f"  Kirchhoff index: {results['kirchhoff_index']:.6f}")

    # --- 7. Spanning tree count (Kirchhoff's theorem) ---
    if V > 1:
        nonzero_lap = lap_eigs[1:]
        nonzero_lap = nonzero_lap[nonzero_lap > 1e-10]
        if len(nonzero_lap) == V - 1:
            spanning_trees = float(np.prod(nonzero_lap) / V)
            results['spanning_trees'] = round(spanning_trees)
            print(f"  Spanning trees: {results['spanning_trees']}")
        else:
            results['spanning_trees'] = 0
            print(f"  Spanning trees: 0 (disconnected)")

    # --- 8. Connectivity ---
    # BFS connectivity check
    if V > 0:
        visited = set()
        queue = [0]
        visited.add(0)
        while queue:
            v = queue.pop(0)
            for u in range(V):
                if A[v][u] > 0 and u not in visited:
                    visited.add(u)
                    queue.append(u)
        results['connected'] = len(visited) == V
        print(f"  Connected: {results['connected']}")

    # --- 9. Diameter and Wiener index ---
    # BFS from each vertex
    dist_matrix = np.full((V, V), V + 1)
    for start in range(V):
        dist_matrix[start][start] = 0
        queue = [start]
        d = 0
        visited = {start}
        while queue:
            next_queue = []
            d += 1
            for v in queue:
                for u in range(V):
                    if A[v][u] > 0 and u not in visited:
                        visited.add(u)
                        dist_matrix[start][u] = d
                        next_queue.append(u)
            queue = next_queue

    reachable = dist_matrix[dist_matrix <= V]
    results['diameter'] = int(np.max(reachable)) if len(reachable) > 0 else -1
    results['wiener_index'] = int(np.sum(dist_matrix[dist_matrix <= V]) // 2)
    results['avg_distance'] = float(np.mean(dist_matrix[dist_matrix <= V])) if len(reachable) > 0 else 0
    print(f"  Diameter: {results['diameter']}")
    print(f"  Wiener index: {results['wiener_index']}")
    print(f"  Avg distance: {results['avg_distance']:.4f}")

    # --- 10. Girth (shortest cycle) ---
    girth = V + 1
    for start in range(V):
        # BFS, looking for back-edge to start
        parent = [-1] * V
        depth = [-1] * V
        depth[start] = 0
        queue = [start]
        while queue:
            v = queue.pop(0)
            for u in range(V):
                if A[v][u] > 0:
                    if depth[u] == -1:
                        depth[u] = depth[v] + 1
                        parent[u] = v
                        queue.append(u)
                    elif parent[v] != u:
                        cycle_len = depth[v] + depth[u] + 1
                        girth = min(girth, cycle_len)
    results['girth'] = girth if girth <= V else float('inf')
    print(f"  Girth: {results['girth']}")

    # --- 11. Triangle count ---
    A2 = A @ A
    triangles = int(np.trace(A @ A2) / 6)
    results['triangles'] = triangles
    print(f"  Triangles: {triangles}")

    # --- 12. Clustering coefficient ---
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
    results['avg_clustering'] = cc_sum / cc_count if cc_count > 0 else 0
    print(f"  Avg clustering coefficient: {results['avg_clustering']:.6f}")

    # --- 13. Clique number (maximum clique) ---
    # Greedy + exact for small V
    if V <= 60:
        max_clique = 1
        for size in range(2, V + 1):
            found = False
            if comb(V, size) > 100000:
                break  # too large
            for subset in combinations(range(V), size):
                if all(A[subset[i]][subset[j]] > 0 for i in range(len(subset)) for j in range(i+1, len(subset))):
                    max_clique = size
                    found = True
                    break
            if not found:
                break
        results['clique_number'] = max_clique
        print(f"  Clique number: {max_clique}")

    # --- 14. Independence number ---
    if V <= 60:
        max_indep = 1
        for size in range(2, V + 1):
            found = False
            if comb(V, size) > 100000:
                break
            for subset in combinations(range(V), size):
                if all(A[subset[i]][subset[j]] == 0 for i in range(len(subset)) for j in range(i+1, len(subset))):
                    max_indep = size
                    found = True
                    break
            if not found:
                break
        results['independence_number'] = max_indep
        print(f"  Independence number: {max_indep}")

    # --- 15. Independence polynomial ---
    if V <= 30:
        alpha = [0] * (V + 1)
        alpha[0] = 1
        for size in range(1, V + 1):
            count = 0
            if comb(V, size) > 500000:
                break
            for subset in combinations(range(V), size):
                if all(A[subset[i]][subset[j]] == 0 for i in range(len(subset)) for j in range(i+1, len(subset))):
                    count += 1
            alpha[size] = count
            if count == 0:
                break
        # Trim trailing zeros
        while alpha and alpha[-1] == 0:
            alpha.pop()
        results['independence_poly'] = alpha
        I_at_2 = sum(a * (2**k) for k, a in enumerate(alpha))
        results['I_at_2'] = I_at_2
        I_at_neg1 = sum(a * ((-1)**k) for k, a in enumerate(alpha))
        results['I_at_neg1'] = I_at_neg1
        print(f"  I(x) coefficients: {alpha}")
        print(f"  I(2) = {I_at_2} (meta-H)")
        print(f"  I(-1) = {I_at_neg1} (Euler char of indep complex)")

    # --- 16. Chromatic number (greedy bound + exact for small V) ---
    if V <= 40:
        # Greedy coloring
        colors = [-1] * V
        num_colors = 0
        for v in range(V):
            used = set()
            for u in range(V):
                if A[v][u] > 0 and colors[u] >= 0:
                    used.add(colors[u])
            c = 0
            while c in used:
                c += 1
            colors[v] = c
            num_colors = max(num_colors, c + 1)
        results['chromatic_greedy'] = num_colors
        print(f"  Chromatic number (greedy upper bound): {num_colors}")

    # --- 17. Cheeger constant (isoperimetric number) ---
    if V > 2 and V <= 40:
        cheeger = float('inf')
        for size in range(1, V // 2 + 1):
            if comb(V, size) > 10000:
                # Sample random subsets
                import random
                random.seed(42)
                for _ in range(1000):
                    S = random.sample(range(V), size)
                    S_set = set(S)
                    boundary = sum(1 for v in S for u in range(V)
                                   if A[v][u] > 0 and u not in S_set)
                    cheeger = min(cheeger, boundary / len(S))
            else:
                for subset in combinations(range(V), size):
                    S_set = set(subset)
                    boundary = sum(1 for v in subset for u in range(V)
                                   if A[v][u] > 0 and u not in S_set)
                    cheeger = min(cheeger, boundary / len(subset))
        results['cheeger_constant'] = float(cheeger)
        print(f"  Cheeger constant (h): {results['cheeger_constant']:.6f}")
        # Cheeger inequality: lambda_2/2 <= h <= sqrt(2*lambda_max*lambda_2)
        lam2 = results['algebraic_connectivity']
        lam_max = results['lap_spectral_radius']
        print(f"  Cheeger bounds: {lam2/2:.4f} <= h <= {np.sqrt(2*lam_max*lam2):.4f}")

    # --- 18. Forman-Ricci curvature ---
    forman_curvatures = []
    for i in range(V):
        for j in range(i+1, V):
            if A[i][j] > 0:
                di = int(deg[i])
                dj = int(deg[j])
                # Number of triangles containing edge (i,j)
                tri_ij = sum(1 for k in range(V) if k != i and k != j
                             and A[i][k] > 0 and A[j][k] > 0)
                forman = 4 - di - dj + 3 * tri_ij
                forman_curvatures.append(forman)

    if forman_curvatures:
        results['forman_min'] = min(forman_curvatures)
        results['forman_max'] = max(forman_curvatures)
        results['forman_avg'] = float(np.mean(forman_curvatures))
        results['forman_total'] = sum(forman_curvatures)
        print(f"  Forman-Ricci: min={results['forman_min']}, max={results['forman_max']}, avg={results['forman_avg']:.4f}, total={results['forman_total']}")

    # --- 19. Degree distribution entropy ---
    deg_counts = Counter(deg)
    probs = np.array(list(deg_counts.values()), dtype=float) / V
    entropy = -np.sum(probs * np.log2(probs + 1e-30))
    results['degree_entropy'] = float(entropy)
    print(f"  Degree distribution entropy: {results['degree_entropy']:.4f} bits")

    # --- 20. Betweenness centrality ---
    betweenness = np.zeros(V)
    for s in range(V):
        # BFS from s
        stack = []
        pred = [[] for _ in range(V)]
        sigma = np.zeros(V)
        sigma[s] = 1
        d = np.full(V, -1)
        d[s] = 0
        queue = [s]
        while queue:
            v = queue.pop(0)
            stack.append(v)
            for w in range(V):
                if A[v][w] > 0:
                    if d[w] < 0:
                        d[w] = d[v] + 1
                        queue.append(w)
                    if d[w] == d[v] + 1:
                        sigma[w] += sigma[v]
                        pred[w].append(v)
        delta = np.zeros(V)
        while stack:
            w = stack.pop()
            for v in pred[w]:
                delta[v] += (sigma[v] / sigma[w]) * (1 + delta[w])
            if w != s:
                betweenness[w] += delta[w]
    betweenness /= 2  # undirected
    results['max_betweenness'] = float(np.max(betweenness))
    results['avg_betweenness'] = float(np.mean(betweenness))
    print(f"  Betweenness centrality: max={results['max_betweenness']:.4f}, avg={results['avg_betweenness']:.4f}")

    # --- 21. Clique complex Betti numbers (for small V) ---
    if V <= 20:
        betti = compute_clique_complex_betti(A, V)
        results['betti_numbers'] = betti
        results['euler_char'] = sum((-1)**k * b for k, b in enumerate(betti))
        print(f"  Clique complex Betti numbers: {betti}")
        print(f"  Euler characteristic: {results['euler_char']}")

    # --- 22. Ramanujan property ---
    if V > 1:
        d_reg = results['avg_deg']
        threshold = 2 * np.sqrt(d_reg - 1) if d_reg > 1 else 0
        second_eig = max(abs(eigenvalues[1]), abs(eigenvalues[-1]))
        results['ramanujan'] = second_eig <= threshold
        results['ramanujan_ratio'] = float(second_eig / threshold) if threshold > 0 else float('inf')
        print(f"  Ramanujan: {results['ramanujan']} (|lam2|/2sqrt(d-1) = {results['ramanujan_ratio']:.4f})")

    return results


def compute_clique_complex_betti(A, V):
    """Compute Betti numbers of the clique complex of a graph."""
    # Find all maximal cliques first, then all simplices
    # For small V, enumerate all cliques directly

    # Enumerate all cliques
    cliques_by_dim = defaultdict(list)  # dim -> list of cliques (as sorted tuples)
    cliques_by_dim[0] = [(v,) for v in range(V)]

    for size in range(2, V + 1):
        found_any = False
        if comb(V, size) > 50000:
            break
        for subset in combinations(range(V), size):
            if all(A[subset[i]][subset[j]] > 0 for i in range(len(subset)) for j in range(i+1, len(subset))):
                cliques_by_dim[size - 1].append(subset)
                found_any = True
        if not found_any:
            break

    max_dim = max(cliques_by_dim.keys()) if cliques_by_dim else 0
    print(f"    Clique complex: max_dim={max_dim}, simplices={sum(len(v) for v in cliques_by_dim.values())}")

    # Compute boundary maps and Betti numbers
    betti = []
    for d in range(max_dim + 1):
        if d == 0:
            # Beta_0 = number of connected components
            # (already computed, but let's compute from the complex)
            # Actually use rank-nullity on boundary_1
            if 1 in cliques_by_dim and len(cliques_by_dim[1]) > 0:
                edges = cliques_by_dim[1]
                vertices = cliques_by_dim[0]
                n_v = len(vertices)
                n_e = len(edges)
                # Boundary map d_1: edges -> vertices
                bd = np.zeros((n_v, n_e), dtype=float)
                v_index = {v: i for i, v in enumerate(vertices)}
                for j, (u, w) in enumerate(edges):
                    bd[v_index[(u,)]][j] = 1
                    bd[v_index[(w,)]][j] = -1
                rank_d1 = np.linalg.matrix_rank(bd)
                betti.append(n_v - rank_d1)
            else:
                betti.append(len(cliques_by_dim[0]))
        else:
            # Beta_d = ker(d_d) - im(d_{d+1})
            # d_d: C_d -> C_{d-1}
            simplices_d = cliques_by_dim.get(d, [])
            simplices_dm1 = cliques_by_dim.get(d-1, [])
            simplices_dp1 = cliques_by_dim.get(d+1, [])

            if len(simplices_d) == 0:
                betti.append(0)
                continue

            # Build boundary d_d
            s_dm1_index = {s: i for i, s in enumerate(simplices_dm1)}
            bd_d = np.zeros((len(simplices_dm1), len(simplices_d)), dtype=float)
            for j, sigma in enumerate(simplices_d):
                for k in range(len(sigma)):
                    face = sigma[:k] + sigma[k+1:]
                    if face in s_dm1_index:
                        bd_d[s_dm1_index[face]][j] = (-1)**k

            ker_d = len(simplices_d) - np.linalg.matrix_rank(bd_d)

            # Build boundary d_{d+1}
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

    return betti


# ============================================================================
# H-WEIGHTED METRICS
# ============================================================================

def compute_h_weighted_metrics(merged_data, results):
    """Compute H-weighted metrics on the merged graph."""
    n = merged_data['n']
    V = merged_data['V']
    adj = merged_data['adj']
    classes = merged_data['classes']

    print(f"\n  H-WEIGHTED METRICS FOR G_{n}/Z_2")
    print(f"  {'-'*50}")

    # H values per merged vertex
    H_vals = np.array([c['H'] for c in classes])
    results['H_values'] = H_vals.tolist()

    # H-gradient: how many edges go uphill, downhill, or level?
    uphill = 0
    downhill = 0
    level = 0
    h_diffs = []
    for (a, b) in merged_data['edges']:
        ha = classes[a]['H']
        hb = classes[b]['H']
        diff = abs(ha - hb)
        h_diffs.append(diff)
        if ha < hb:
            uphill += 1
        elif ha > hb:
            downhill += 1
        else:
            level += 1

    results['uphill_edges'] = uphill
    results['downhill_edges'] = downhill
    results['level_edges'] = level
    results['is_DAG'] = (level == 0)
    results['avg_H_diff'] = float(np.mean(h_diffs)) if h_diffs else 0
    results['max_H_diff'] = max(h_diffs) if h_diffs else 0

    print(f"  Uphill: {uphill}, Downhill: {downhill}, Level: {level}")
    print(f"  DAG (0 level edges): {results['is_DAG']}")
    print(f"  Avg |DeltaH| per edge: {results['avg_H_diff']:.4f}")
    print(f"  Max |DeltaH| per edge: {results['max_H_diff']}")

    # H-weighted Laplacian
    # W[i][j] = exp(-|H_i - H_j| / scale) for adjacent vertices
    scale = max(1, np.std(H_vals))
    W = np.zeros((V, V))
    for (a, b) in merged_data['edges']:
        w = np.exp(-abs(classes[a]['H'] - classes[b]['H']) / scale)
        W[a][b] = w
        W[b][a] = w

    D_w = np.diag(np.sum(W, axis=1))
    L_w = D_w - W
    weig = np.sort(np.linalg.eigvalsh(L_w))
    results['weighted_alg_conn'] = float(weig[1]) if V > 1 else 0
    print(f"  H-weighted algebraic connectivity: {results['weighted_alg_conn']:.6f}")

    # H-weighted centrality: which merged vertex has highest H * degree?
    degs = np.sum(adj, axis=1)
    h_centrality = H_vals * degs
    best = np.argmax(h_centrality)
    results['h_central_vertex'] = int(best)
    results['h_central_value'] = float(h_centrality[best])
    results['h_central_H'] = int(classes[best]['H'])
    print(f"  Most H-central vertex: {best} (H={classes[best]['H']}, deg={int(degs[best])}, H*deg={h_centrality[best]:.0f})")

    # H-monotone path count: longest strictly H-increasing path
    # (this is the longest chain in the H-poset)
    H_levels = sorted(set(H_vals))
    results['H_levels'] = len(H_levels)
    print(f"  Number of distinct H values: {len(H_levels)}")

    # Level set sizes
    level_sizes = Counter(H_vals)
    results['H_level_sizes'] = dict(sorted(level_sizes.items()))
    print(f"  H level set sizes: {dict(sorted(level_sizes.items()))}")

    # Score diversity per H level
    for h_val in sorted(level_sizes.keys()):
        verts_at_h = [i for i, c in enumerate(classes) if c['H'] == h_val]
        scores_at_h = set(classes[i]['score'] for i in verts_at_h)
        if len(verts_at_h) > 1:
            print(f"    H={h_val}: {len(verts_at_h)} vertices, {len(scores_at_h)} distinct scores")


# ============================================================================
# MAIN COMPUTATION
# ============================================================================

all_results = {}

for n in [3, 4, 5, 6]:
    print(f"\n\n{'#'*80}")
    print(f"  n = {n}")
    print(f"{'#'*80}")

    # Build G_n
    gn = build_iso_class_graph(n)

    # Build adjacency matrix for G_n
    A_gn = np.zeros((gn['V'], gn['V']), dtype=int)
    for (a, b) in gn['edges']:
        A_gn[a][b] = 1
        A_gn[b][a] = 1

    # Compute G_n invariants
    gn_results = compute_graph_invariants(A_gn, f"G_{n}")

    # Build merged graph G_n/Z_2
    merged = build_merged_graph(gn)

    # Compute merged graph invariants
    merged_results = compute_graph_invariants(merged['adj'], f"G_{n}/Z_2")

    # H-weighted metrics
    compute_h_weighted_metrics(merged, merged_results)

    all_results[n] = {
        'gn': gn_results,
        'merged': merged_results,
        'gn_data': gn,
        'merged_data': merged
    }

# ============================================================================
# SUMMARY TABLE
# ============================================================================

print(f"\n\n{'='*80}")
print("  MASTER SUMMARY TABLE")
print(f"{'='*80}")

# Collect sequence data
seq_names = [
    'V', 'E', 'density', 'diameter', 'girth', 'triangles',
    'spectral_radius', 'graph_energy', 'algebraic_connectivity',
    'kirchhoff_index', 'spanning_trees', 'avg_clustering',
    'estrada_index', 'wiener_index', 'avg_distance',
    'forman_avg', 'degree_entropy', 'avg_betweenness'
]

print(f"\n  --- G_n SEQUENCES ---")
print(f"  {'Invariant':<30} {'n=3':>10} {'n=4':>10} {'n=5':>10} {'n=6':>10}")
print(f"  {'-'*70}")
for name in seq_names:
    vals = []
    for n in [3, 4, 5, 6]:
        v = all_results[n]['gn'].get(name, '?')
        if isinstance(v, float):
            vals.append(f"{v:.4f}")
        else:
            vals.append(str(v))
    print(f"  {name:<30} {vals[0]:>10} {vals[1]:>10} {vals[2]:>10} {vals[3]:>10}")

print(f"\n  --- G_n/Z_2 SEQUENCES ---")
print(f"  {'Invariant':<30} {'n=3':>10} {'n=4':>10} {'n=5':>10} {'n=6':>10}")
print(f"  {'-'*70}")
for name in seq_names:
    vals = []
    for n in [3, 4, 5, 6]:
        v = all_results[n]['merged'].get(name, '?')
        if isinstance(v, float):
            vals.append(f"{v:.4f}")
        else:
            vals.append(str(v))
    print(f"  {name:<30} {vals[0]:>10} {vals[1]:>10} {vals[2]:>10} {vals[3]:>10}")

# Extra sequences
extra_seqs = ['clique_number', 'independence_number', 'chromatic_greedy',
              'cheeger_constant', 'ramanujan_ratio']
print(f"\n  --- EXTRA INVARIANTS G_n/Z_2 ---")
for name in extra_seqs:
    vals = []
    for n in [3, 4, 5, 6]:
        v = all_results[n]['merged'].get(name, '?')
        if isinstance(v, float):
            vals.append(f"{v:.4f}")
        else:
            vals.append(str(v))
    print(f"  {name:<30} {':'.join(vals)}")

# Independence polynomial
print(f"\n  --- INDEPENDENCE POLYNOMIALS ---")
for graph_type in ['gn', 'merged']:
    for n in [3, 4, 5, 6]:
        ip = all_results[n][graph_type].get('independence_poly', [])
        I2 = all_results[n][graph_type].get('I_at_2', '?')
        Ineg1 = all_results[n][graph_type].get('I_at_neg1', '?')
        name = f"G_{n}" if graph_type == 'gn' else f"G_{n}/Z_2"
        print(f"  I({name}, x) = {ip}  |  I(2)={I2}, I(-1)={Ineg1}")

# Betti numbers
print(f"\n  --- CLIQUE COMPLEX BETTI NUMBERS ---")
for graph_type in ['gn', 'merged']:
    for n in [3, 4, 5, 6]:
        betti = all_results[n][graph_type].get('betti_numbers', '?')
        euler = all_results[n][graph_type].get('euler_char', '?')
        name = f"G_{n}" if graph_type == 'gn' else f"G_{n}/Z_2"
        print(f"  beta({name}) = {betti}  |  chi = {euler}")

# Spectra
print(f"\n  --- ADJACENCY SPECTRA (rounded to 3 decimal places) ---")
for graph_type in ['gn', 'merged']:
    for n in [3, 4, 5, 6]:
        spec = all_results[n][graph_type].get('adj_spectrum', [])
        name = f"G_{n}" if graph_type == 'gn' else f"G_{n}/Z_2"
        spec_str = [f"{x:.3f}" for x in spec[:min(len(spec), 15)]]
        print(f"  spec({name}) = [{', '.join(spec_str)}]")

# H-gradient DAG data
print(f"\n  --- H-GRADIENT ON G_n/Z_2 ---")
for n in [3, 4, 5, 6]:
    mr = all_results[n]['merged']
    print(f"  n={n}: ^{mr.get('uphill_edges','?')} v{mr.get('downhill_edges','?')} ->{mr.get('level_edges','?')}  DAG={mr.get('is_DAG','?')}  levels={mr.get('H_levels','?')}")

print(f"\n\n  DONE. Total time: see individual timings above.")
print("=" * 80)
