#!/usr/bin/env python3
"""
tutte_metagraph_s20co.py — Tutte polynomial of G_n/Z_2
kind-pasteur-2026-03-23-S20co

Computes the Tutte polynomial T(x,y) of the merged meta-graph for n=3,4,5.
n=6 has 34 vertices and 143 edges -- too large for naive deletion-contraction
but we can extract special evaluations.

The Tutte polynomial encodes:
  T(1,1) = number of spanning trees
  T(2,1) = number of spanning forests
  T(1,2) = number of connected spanning subgraphs
  T(2,2) = 2^|E|
  T(0,2) = number of spanning subgraphs with same number of components as G
  T(2,0) = number of acyclic orientations
  T(0,-2) = (-1)^|V|-1 * flow polynomial at k=2 (if connected)

Also computes:
  - Chromatic polynomial P(k) from Tutte
  - Reliability polynomial R(p) from Tutte
  - Flow polynomial F(k) from Tutte
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  TUTTE POLYNOMIAL OF G_n/Z_2")
print("  kind-pasteur-2026-03-23-S20co")
print("=" * 80)

# ============================================================================
# TOURNAMENT HELPERS (same as before)
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
            'cid': idx, 'canon': canon, 'members': members,
            'adj': adj, 'H': h, 'sc': sc, 'score': score_seq(adj, n),
            'c3': c3_count(adj, n), 'size': len(members)
        })
        canon_to_cid[canon] = idx

    for data in class_list:
        comp = complement_adj(data['adj'], n)
        comp_canon = canonical_form(comp, n)
        data['comp_cid'] = canon_to_cid.get(comp_canon, -1)

    edges = set()
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
                if nb is not None and nb != cid:
                    e = (min(cid, nb), max(cid, nb))
                    edges.add(e)
                    if e not in edge_colors:
                        edge_colors[e] = 'blue' if data['sc'] == class_list[nb]['sc'] else 'black'

    V = len(class_list)
    sc_count = sum(1 for d in class_list if d['sc'])

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

    merged_classes = []
    seen = set()
    for data in class_list:
        mid_val = merged_id[data['cid']]
        if mid_val not in seen:
            seen.add(mid_val)
            merged_classes.append({'mid': mid_val, 'H': data['H'], 'sc': data['sc']})

    return V_merged, sorted(merged_edges), merged_classes


# ============================================================================
# TUTTE POLYNOMIAL via deletion-contraction
# ============================================================================

def connected_components_count(V, edge_set):
    """Count connected components of graph with V vertices and given edge set."""
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
    for (a, b) in edge_set:
        union(a, b)
    return len(set(find(x) for x in range(V)))


def tutte_polynomial_subset_expansion(V, edges):
    """Compute Tutte polynomial via subset expansion (Sokal's formula).

    T(x,y) = sum over A subset of E: (x-1)^(k(A)-k(E)) * (y-1)^(k(A)+|A|-|V(A)|)

    where k(A) = number of connected components of (V, A).
    Actually, the correct formula uses the full vertex set V for all subsets.
    """
    E = len(edges)
    print(f"    Tutte via subset expansion: V={V}, E={E}")

    if E > 25:
        print(f"    Too many edges for subset expansion ({E} > 25). Skipping.")
        return None

    k_E = connected_components_count(V, edges)  # k(E)

    # Store as polynomial in (x-1) and (y-1)
    # T(x,y) = sum c_{i,j} (x-1)^i (y-1)^j
    max_deg = V + E
    coeffs = defaultdict(int)  # (i, j) -> coefficient

    for mask in range(1 << E):
        subset = [edges[k] for k in range(E) if mask & (1 << k)]
        size = bin(mask).count('1')
        k_A = connected_components_count(V, subset)

        # nullity = |A| + k(A) - |V| (for the subgraph)
        # Actually for Tutte: rank r(A) = |V| - k(A), nullity n(A) = |A| - r(A) = |A| - |V| + k(A)
        r_A = V - k_A
        r_E = V - k_E
        n_A = size - r_A

        i = r_E - r_A  # power of (x-1)
        j = n_A         # power of (y-1)

        if i >= 0 and j >= 0:
            coeffs[(i, j)] += 1

    return dict(coeffs)


def evaluate_tutte(coeffs, x, y):
    """Evaluate T(x,y) from coefficient dict."""
    if coeffs is None:
        return None
    result = 0
    for (i, j), c in coeffs.items():
        result += c * (x-1)**i * (y-1)**j
    return result


def chromatic_from_tutte(coeffs, V, k):
    """Chromatic polynomial P(G, k) = (-1)^(V-k(E)) * k^(k(E)) * T(1-k, 0).
    For connected graphs: P(G, k) = (-1)^(V-1) * k * T(1-k, 0)."""
    if coeffs is None:
        return None
    return (-1)**(V-1) * k * evaluate_tutte(coeffs, 1-k, 0)


# ============================================================================
# MAIN
# ============================================================================

for n in [3, 4, 5]:
    print(f"\n{'#'*80}")
    print(f"  n = {n}")
    print(f"{'#'*80}")
    t0 = time.time()

    V, edges, classes = build_merged(n)
    E = len(edges)
    print(f"  G_{n}/Z_2: V={V}, E={E}")

    # Compute Tutte polynomial
    coeffs = tutte_polynomial_subset_expansion(V, edges)

    if coeffs is not None:
        print(f"\n  Tutte polynomial T(x,y) coefficients (in basis (x-1)^i * (y-1)^j):")
        max_i = max(i for (i, j) in coeffs.keys())
        max_j = max(j for (i, j) in coeffs.keys())
        for i in range(max_i + 1):
            row = []
            for j in range(max_j + 1):
                c = coeffs.get((i, j), 0)
                row.append(str(c))
            print(f"    i={i}: {', '.join(row)}")

        # Special evaluations
        print(f"\n  Special evaluations:")
        t11 = evaluate_tutte(coeffs, 1, 1)
        t21 = evaluate_tutte(coeffs, 2, 1)
        t12 = evaluate_tutte(coeffs, 1, 2)
        t22 = evaluate_tutte(coeffs, 2, 2)
        t20 = evaluate_tutte(coeffs, 2, 0)
        t02 = evaluate_tutte(coeffs, 0, 2)
        print(f"    T(1,1) = {t11} (spanning trees)")
        print(f"    T(2,1) = {t21} (spanning forests)")
        print(f"    T(1,2) = {t12} (connected spanning subgraphs)")
        print(f"    T(2,2) = {t22} (= 2^E = {2**E})")
        print(f"    T(2,0) = {t20} (acyclic orientations)")
        print(f"    T(0,2) = {t02} (spanning subgraphs same #components)")

        # Chromatic polynomial evaluations
        print(f"\n  Chromatic polynomial P(G, k):")
        for k in range(1, 8):
            pk = chromatic_from_tutte(coeffs, V, k)
            print(f"    P({k}) = {pk}")

        # Chromatic number = smallest k with P(k) > 0
        for k in range(1, 20):
            pk = chromatic_from_tutte(coeffs, V, k)
            if pk > 0:
                print(f"    Chromatic number: {k}")
                break

        # Reliability polynomial R(p) at various p
        print(f"\n  Reliability polynomial R(p) (probability graph stays connected):")
        for p10 in [1, 2, 3, 5, 7, 9, 10]:
            p = p10 / 10.0
            # R(p) = sum over connected spanning subgraphs S:
            # R(p) = (1-p)^(V-k(E)) * p^(k(E)) * sum c_{i,j} (p/(1-p))^j * ((1-p)/p * something)
            # Actually: R(p) = p^(V-1) * (1-p)^(E-V+1) * T(1, 1/(1-p))
            if p < 1 and p > 0:
                rp = p**(V-1) * (1-p)**(E-V+1) * evaluate_tutte(coeffs, 1, 1/(1-p))
                print(f"    R({p:.1f}) = {rp:.6f}")

    print(f"  Time: {time.time()-t0:.1f}s")


# ============================================================================
# n=6 SPECIAL EVALUATIONS (too large for full Tutte, use subset formulas)
# ============================================================================

print(f"\n{'#'*80}")
print(f"  n = 6 (SPECIAL EVALUATIONS ONLY)")
print(f"{'#'*80}")

t0 = time.time()
V6, edges6, classes6 = build_merged(6)
E6 = len(edges6)
print(f"  G_6/Z_2: V={V6}, E={E6}")
print(f"  2^E = {2**E6} -- too large for subset expansion")

# Chromatic polynomial via deletion-contraction would be O(phi^E) ~ huge
# But we can compute P(k) directly by counting proper colorings
print(f"\n  Computing chromatic polynomial via direct coloring count...")

# Build adjacency list
adj_list = defaultdict(set)
for (a, b) in edges6:
    adj_list[a].add(b)
    adj_list[b].add(a)

def count_proper_colorings(V, adj_list, k, max_count=10**8):
    """Count proper k-colorings via backtracking."""
    count = 0
    colors = [-1] * V

    def backtrack(v):
        nonlocal count
        if v == V:
            count += 1
            if count >= max_count:
                return
            return
        used = set()
        for u in adj_list[v]:
            if colors[u] >= 0:
                used.add(colors[u])
        for c in range(k):
            if c not in used:
                colors[v] = c
                backtrack(v + 1)
                if count >= max_count:
                    return
                colors[v] = -1

    backtrack(0)
    return count

print(f"  P(1) = 0 (trivially)")
for k in range(2, 8):
    pk = count_proper_colorings(V6, adj_list, k, max_count=10**9)
    print(f"  P({k}) = {pk}")
    if pk > 0 and k <= 6:
        print(f"    -> Chromatic number <= {k}")

print(f"\n  Time: {time.time()-t0:.1f}s")

# ============================================================================
# SUMMARY
# ============================================================================

print(f"\n\n{'='*80}")
print("  TUTTE POLYNOMIAL SUMMARY")
print(f"{'='*80}")
print("""
  G_3/Z_2: K_2 (complete graph on 2 vertices)
    T(x,y) = x
    P(k) = k(k-1)
    chi = 2

  G_4/Z_2: K_3 (complete graph on 3 vertices)
    T(x,y) = x^2 + x + y
    P(k) = k(k-1)(k-2)
    chi = 3

  G_5/Z_2: 10 vertices, 21 edges
    T(x,y) = (see above)
    P(k) = (see above)
    chi = (see above)

  G_6/Z_2: 34 vertices, 143 edges
    chi = (see above, from direct coloring)
""")

print("  DONE.")
