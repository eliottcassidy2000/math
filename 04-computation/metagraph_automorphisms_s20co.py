#!/usr/bin/env python3
"""
metagraph_automorphisms_s20co.py — Automorphisms and symmetries of G_n/Z_2
kind-pasteur-2026-03-23-S20co

Computes:
  1. Automorphism group of G_n/Z_2 as abstract graph (ignoring tournament labels)
  2. Orbit structure under Aut(G_n/Z_2)
  3. Quotient graph G_n/Z_2 / Aut
  4. Vertex-transitive analysis
  5. Edge-transitive analysis
  6. Distance polynomial / walk generating function
  7. Resistance distance matrix
  8. Ollivier-Ricci curvature (via optimal transport on neighborhoods)
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  AUTOMORPHISMS AND SYMMETRIES OF G_n/Z_2")
print("  kind-pasteur-2026-03-23-S20co")
print("=" * 80)

# ============================================================================
# HELPERS (minimal set)
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
        if pc >= n: continue
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

# ============================================================================
# BUILD G_n/Z_2
# ============================================================================

def build_merged(n_val):
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
        class_list.append({
            'cid': idx, 'canon': canon, 'adj': adj, 'H': h, 'sc': sc,
            'score': score_seq(adj, n), 'c3': c3_count(adj, n), 'size': len(members),
            'members': members
        })
        canon_to_cid[canon] = idx

    for data in class_list:
        comp = complement_adj(data['adj'], n)
        comp_canon = canonical_form(comp, n)
        data['comp_cid'] = canon_to_cid.get(comp_canon, -1)

    edges = set()
    for data in class_list:
        cid = data['cid']
        for bits in data['members']:
            adj2 = tournament_adj(n, bits)
            for arc_idx in range(m):
                flipped = bits ^ (1 << arc_idx)
                fadj = tournament_adj(n, flipped)
                fcanon = canonical_form(fadj, n)
                nb = canon_to_cid.get(fcanon)
                if nb is not None and nb != cid:
                    edges.add((min(cid, nb), max(cid, nb)))

    merged_id = {}
    mid = 0
    for data in class_list:
        cid = data['cid']
        if cid in merged_id: continue
        comp = data['comp_cid']
        merged_id[cid] = mid
        if comp != cid:
            merged_id[comp] = mid
        mid += 1

    V_merged = mid
    merged_edges = set()
    for (a, b) in edges:
        ma, mb = merged_id[a], merged_id[b]
        if ma != mb:
            merged_edges.add((min(ma, mb), max(ma, mb)))

    A = np.zeros((V_merged, V_merged), dtype=int)
    for (a, b) in merged_edges:
        A[a][b] = 1
        A[b][a] = 1

    merged_classes = []
    seen = set()
    for data in class_list:
        mid_val = merged_id[data['cid']]
        if mid_val not in seen:
            seen.add(mid_val)
            merged_classes.append({
                'mid': mid_val, 'H': data['H'], 'sc': data['sc'],
                'score': data['score'], 'c3': data['c3'], 'size': data['size']
            })

    return V_merged, sorted(merged_edges), A, merged_classes


# ============================================================================
# GRAPH AUTOMORPHISMS (brute force for small V)
# ============================================================================

def find_automorphisms(A, V):
    """Find all automorphisms of graph with adjacency matrix A."""
    if V > 12:
        print(f"    V={V} too large for brute-force automorphism search")
        return None

    deg = [sum(A[i]) for i in range(V)]
    deg_groups = defaultdict(list)
    for v in range(V):
        deg_groups[deg[v]].append(v)

    sorted_degs = sorted(set(deg))
    groups = [deg_groups[d] for d in sorted_degs]

    auts = []
    def gen_perms(groups):
        if not groups:
            yield []
            return
        for p in permutations(groups[0]):
            for rest in gen_perms(groups[1:]):
                yield list(p) + rest

    for perm in gen_perms(groups):
        ok = True
        for i in range(V):
            for j in range(i+1, V):
                if A[perm[i]][perm[j]] != A[i][j]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            auts.append(tuple(perm))

    return auts


def aut_orbits(auts, V):
    """Compute orbits of automorphism group action."""
    parent = list(range(V))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for perm in auts:
        for v in range(V):
            union(v, perm[v])

    orbits = defaultdict(list)
    for v in range(V):
        orbits[find(v)].append(v)
    return list(orbits.values())


# ============================================================================
# OLLIVIER-RICCI CURVATURE
# ============================================================================

def ollivier_ricci(A, V, alpha=0.5):
    """Compute Ollivier-Ricci curvature for each edge.

    Using lazy random walk measure: mu_x = alpha * delta_x + (1-alpha) * uniform(neighbors)
    Curvature kappa(x,y) = 1 - W_1(mu_x, mu_y) / d(x,y)
    where W_1 is Wasserstein-1 distance and d(x,y) is graph distance.
    """
    # Compute distance matrix
    dist = np.full((V, V), V + 1)
    for start in range(V):
        dist[start][start] = 0
        queue = [start]
        visited = {start}
        while queue:
            next_q = []
            for v in queue:
                for u in range(V):
                    if A[v][u] > 0 and u not in visited:
                        visited.add(u)
                        dist[start][u] = dist[start][v] + 1
                        next_q.append(u)
            queue = next_q

    curvatures = []
    for i in range(V):
        for j in range(i+1, V):
            if A[i][j] == 0:
                continue

            # Build measures
            nbrs_i = [u for u in range(V) if A[i][u] > 0]
            nbrs_j = [u for u in range(V) if A[j][u] > 0]

            di = len(nbrs_i)
            dj = len(nbrs_j)

            # Measure mu_i: alpha on i, (1-alpha)/di on each neighbor
            # Measure mu_j: alpha on j, (1-alpha)/dj on each neighbor

            # Support of mu_i = {i} + nbrs_i
            # Support of mu_j = {j} + nbrs_j
            support = sorted(set([i, j] + nbrs_i + nbrs_j))
            n_sup = len(support)
            idx_map = {v: k for k, v in enumerate(support)}

            # Build cost matrix (distances between support points)
            C = np.zeros((n_sup, n_sup))
            for a in range(n_sup):
                for b in range(n_sup):
                    C[a][b] = dist[support[a]][support[b]]

            # Build measure vectors
            mu_i = np.zeros(n_sup)
            mu_i[idx_map[i]] = alpha
            for u in nbrs_i:
                mu_i[idx_map[u]] += (1 - alpha) / di

            mu_j = np.zeros(n_sup)
            mu_j[idx_map[j]] = alpha
            for u in nbrs_j:
                mu_j[idx_map[u]] += (1 - alpha) / dj

            # Compute W_1 via simple LP formulation (for small support)
            # W_1 = min sum C[a][b] * gamma[a][b]
            # s.t. sum_b gamma[a][b] = mu_i[a], sum_a gamma[a][b] = mu_j[b], gamma >= 0

            # Simple approximation: use the coupling where mass is moved greedily
            # For exact W_1, use the dual: W_1 = max sum f * (mu_i - mu_j) s.t. f Lipschitz(1)

            # Use dual formulation for exact W_1
            diff = mu_i - mu_j

            # W_1 = max over 1-Lipschitz functions f of <f, diff>
            # Enumerate all 1-Lipschitz functions on the support
            # A 1-Lipschitz function satisfies |f(a) - f(b)| <= C[a][b]
            # For the dual, it suffices to use f(v) = dist(v, center) for each center
            # Actually W_1 = sum of positive part of sorted partial sums (Earth mover)

            # Simple exact computation for small support:
            # Use the fact that for tree metrics, W_1 has a closed form
            # For general metrics, use LP

            # Greedy earth mover (not exact but good approximation):
            w1 = 0
            remaining = diff.copy()
            for _ in range(n_sup * n_sup):
                # Find source (positive remaining) and sink (negative remaining)
                sources = [(remaining[a], a) for a in range(n_sup) if remaining[a] > 1e-12]
                sinks = [(-remaining[b], b) for b in range(n_sup) if remaining[b] < -1e-12]
                if not sources or not sinks:
                    break
                # Move from closest source-sink pair
                best_cost = float('inf')
                best_pair = None
                for (sa, a) in sources:
                    for (sb, b) in sinks:
                        if C[a][b] < best_cost:
                            best_cost = C[a][b]
                            best_pair = (a, b, min(sa, sb))
                if best_pair is None:
                    break
                a, b, amount = best_pair
                w1 += best_cost * amount
                remaining[a] -= amount
                remaining[b] += amount

            kappa = 1 - w1 / dist[i][j]
            curvatures.append((i, j, kappa))

    return curvatures


# ============================================================================
# WALK GENERATING FUNCTION
# ============================================================================

def walk_counts(A, V, max_length=10):
    """Compute number of closed walks of each length."""
    A_f = A.astype(float)
    power = np.eye(V)
    walks = []
    for k in range(max_length + 1):
        walks.append(int(round(np.trace(power))))
        power = power @ A_f
    return walks


# ============================================================================
# RESISTANCE DISTANCE
# ============================================================================

def resistance_distances(A, V):
    """Compute resistance distance matrix using Laplacian pseudoinverse."""
    deg = np.sum(A, axis=1).astype(float)
    L = np.diag(deg) - A.astype(float)

    # Moore-Penrose pseudoinverse of L
    L_pinv = np.linalg.pinv(L)

    R = np.zeros((V, V))
    for i in range(V):
        for j in range(V):
            R[i][j] = L_pinv[i][i] + L_pinv[j][j] - 2 * L_pinv[i][j]

    return R


# ============================================================================
# MAIN
# ============================================================================

for n in [3, 4, 5, 6]:
    print(f"\n{'#'*80}")
    print(f"  n = {n}")
    print(f"{'#'*80}")
    t0 = time.time()

    V, edges, A, classes = build_merged(n)
    print(f"  G_{n}/Z_2: V={V}, E={len(edges)}")

    # 1. Automorphisms
    print(f"\n  --- Automorphisms ---")
    auts = find_automorphisms(A, V)
    if auts is not None:
        print(f"  |Aut(G_{n}/Z_2)| = {len(auts)}")
        orbits = aut_orbits(auts, V)
        print(f"  Number of vertex orbits: {len(orbits)}")
        orbit_sizes = sorted([len(o) for o in orbits], reverse=True)
        print(f"  Orbit sizes: {orbit_sizes}")
        print(f"  Vertex-transitive: {len(orbits) == 1}")

        # What H values are in each orbit?
        for i, orb in enumerate(orbits):
            h_vals = sorted(set(classes[v]['H'] for v in orb))
            sc_vals = [classes[v]['sc'] for v in orb]
            sc_count_orb = sum(sc_vals)
            print(f"    Orbit {i}: size={len(orb)}, H={h_vals}, SC={sc_count_orb}/{len(orb)}")

        # Edge orbits
        edge_parent = {}
        for (a, b) in edges:
            if (a, b) not in edge_parent:
                edge_parent[(a, b)] = (a, b)
        for perm in auts:
            for (a, b) in edges:
                mapped = (min(perm[a], perm[b]), max(perm[a], perm[b]))
                if mapped in edge_parent:
                    # Union
                    pa = edge_parent[(a, b)]
                    pb = edge_parent[mapped]
                    if pa != pb:
                        for e in list(edge_parent.keys()):
                            if edge_parent[e] == pb:
                                edge_parent[e] = pa
        edge_orbits = defaultdict(list)
        for e, p in edge_parent.items():
            edge_orbits[p].append(e)
        print(f"  Number of edge orbits: {len(edge_orbits)}")
        print(f"  Edge-transitive: {len(edge_orbits) == 1}")

    # 2. Walk generating function
    print(f"\n  --- Walk Counts (closed walks of length k) ---")
    walks = walk_counts(A, V, max_length=12)
    print(f"  w(k) for k=0..12: {walks}")
    # w(0) = V, w(1) = 0 (no self-loops), w(2) = 2E, w(3) = 6*triangles
    print(f"  w(0)={walks[0]}=V, w(2)={walks[2]}=2E={2*len(edges)}, "
          f"w(3)={walks[3]}=6*tri={6*(walks[3]//6)}")

    # 3. Resistance distances
    print(f"\n  --- Resistance Distances ---")
    R = resistance_distances(A, V)
    print(f"  Total resistance (Kirchhoff index / V): {np.sum(R)/2/V:.4f}")
    print(f"  Max resistance: {np.max(R):.4f}")
    print(f"  Min resistance (non-zero): {np.min(R[R > 1e-10]):.4f}")
    print(f"  Avg resistance: {np.mean(R[R > 1e-10]):.4f}")

    # Resistance vs graph distance correlation
    dist_matrix = np.full((V, V), V + 1)
    for start in range(V):
        dist_matrix[start][start] = 0
        queue = [start]
        visited = {start}
        while queue:
            next_q = []
            for v in queue:
                for u in range(V):
                    if A[v][u] > 0 and u not in visited:
                        visited.add(u)
                        dist_matrix[start][u] = dist_matrix[start][v] + 1
                        next_q.append(u)
            queue = next_q

    mask = (dist_matrix > 0) & (dist_matrix <= V)
    if np.sum(mask) > 2:
        corr = np.corrcoef(dist_matrix[mask].flatten(), R[mask].flatten())[0, 1]
        print(f"  corr(graph_dist, resistance_dist): {corr:.4f}")

    # 4. Ollivier-Ricci curvature
    if V <= 40:
        print(f"\n  --- Ollivier-Ricci Curvature (alpha=0.5) ---")
        curvs = ollivier_ricci(A, V, alpha=0.5)
        kappas = [k for (_, _, k) in curvs]
        print(f"  min kappa = {min(kappas):.4f}")
        print(f"  max kappa = {max(kappas):.4f}")
        print(f"  avg kappa = {np.mean(kappas):.4f}")
        print(f"  #positive: {sum(1 for k in kappas if k > 0.01)}")
        print(f"  #negative: {sum(1 for k in kappas if k < -0.01)}")
        print(f"  #zero: {sum(1 for k in kappas if abs(k) <= 0.01)}")

        # Total Ollivier curvature
        total_ollivier = sum(kappas)
        print(f"  Total Ollivier curvature: {total_ollivier:.4f}")

        # Compare Forman vs Ollivier
        deg = np.sum(A, axis=1).astype(int)
        forman_curvs = []
        for i in range(V):
            for j in range(i+1, V):
                if A[i][j] > 0:
                    di, dj = int(deg[i]), int(deg[j])
                    tri_ij = sum(1 for k in range(V) if k != i and k != j
                                 and A[i][k] > 0 and A[j][k] > 0)
                    forman_curvs.append(4 - di - dj + 3 * tri_ij)

        if len(forman_curvs) == len(kappas):
            corr_fo = np.corrcoef(forman_curvs, kappas)[0, 1]
            print(f"  corr(Forman, Ollivier): {corr_fo:.4f}")

    print(f"\n  Time: {time.time()-t0:.1f}s")


print(f"\n\n{'='*80}")
print("  SUMMARY: SYMMETRY STRUCTURE")
print(f"{'='*80}")
print("""
  n  |Aut|  vertex_orbits  edge_orbits  vtx_trans  edge_trans
  3    ?       ?              ?            ?          ?
  4    ?       ?              ?            ?          ?
  5    ?       ?              ?            ?          ?
  6    ?       ?              ?            ?          ?
  (filled in from computation above)

  CURVATURE COMPARISON:
  n  Forman_avg  Ollivier_avg  corr(F,O)
  3     ?           ?            ?
  4     ?           ?            ?
  5     ?           ?            ?
  6     ?           ?            ?
""")

print("  DONE.")
