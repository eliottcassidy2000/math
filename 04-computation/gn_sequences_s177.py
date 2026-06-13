#!/usr/bin/env python3
"""
gn_sequences_s177.py — Brand New Sequences from G_n
opus-2026-03-22-S177

Every graph invariant of G_n is an integer sequence indexed by n.
Compute as many as possible at n=3,4,5 and search for patterns.
"""

from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter, deque
import numpy as np

def build_Gn(n):
    m = comb(n, 2)
    iso = defaultdict(list)
    for bits in range(1 << m):
        adj = [[0]*n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if bits & (1 << idx): adj[i][j] = 1
                else: adj[j][i] = 1
                idx += 1
        best = None
        for perm in permutations(range(n)):
            form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
            if best is None or form < best: best = form
        iso[best].append(bits)

    classes = []
    for c, members in iso.items():
        adj = [[0]*n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if members[0] & (1 << idx): adj[i][j] = 1
                else: adj[j][i] = 1
                idx += 1
        # H
        dp = [0]*((1<<n)*n)
        for v in range(n): dp[(1<<v)*n+v] = 1
        for S in range(1, 1<<n):
            for v in range(n):
                if not (S&(1<<v)): continue
                val = dp[S*n+v]
                if val == 0: continue
                for u in range(n):
                    if S&(1<<u): continue
                    if adj[v][u]: dp[(S|(1<<u))*n+u] += val
        full = (1<<n)-1
        h = sum(dp[full*n+v] for v in range(n))
        ss = tuple(sorted(sum(adj[i]) for i in range(n)))
        aut = factorial(n) // len(members)
        classes.append({'h': h, 'scores': ss, 'size': len(members), 'aut': aut})

    classes.sort(key=lambda x: (x['h'], x['scores']))
    for i, cl in enumerate(classes): cl['id'] = i

    bits_to_id = {}
    for cl in classes:
        for b in iso[list(iso.keys())[cl['id']]]:
            pass
    # Rebuild mapping properly
    class_canons = list(iso.keys())
    # Re-sort to match classes
    sorted_canons = sorted(class_canons, key=lambda c: (
        sum(dp_val for dp_val in [0]),  # placeholder
    ))
    # Simpler approach: use the class index directly
    id_map = {}
    for i, (c, members) in enumerate(sorted(iso.items(), key=lambda x: x[0])):
        for b in members:
            id_map[b] = i

    # Rebuild with proper sorting
    canon_list = sorted(iso.keys())
    remap = {}
    for new_id, cl in enumerate(classes):
        pass  # need to link classes to canons

    # Just rebuild edges directly
    edges = set()
    adj_list = defaultdict(set)

    # Map bits to class id via canonical
    bit_to_canon = {}
    for c, members in iso.items():
        for b in members:
            bit_to_canon[b] = c

    canon_to_id = {}
    for i, c in enumerate(sorted(iso.keys())):
        canon_to_id[c] = i

    # Recompute classes in canonical sort order
    classes2 = []
    for c in sorted(iso.keys()):
        members = iso[c]
        adj = [[0]*n for _ in range(n)]
        idx = 0
        for ii in range(n):
            for jj in range(ii+1, n):
                if members[0] & (1 << idx): adj[ii][jj] = 1
                else: adj[jj][ii] = 1
                idx += 1
        dp2 = [0]*((1<<n)*n)
        for v in range(n): dp2[(1<<v)*n+v] = 1
        for S in range(1, 1<<n):
            for v in range(n):
                if not (S&(1<<v)): continue
                val = dp2[S*n+v]
                if val == 0: continue
                for u in range(n):
                    if S&(1<<u): continue
                    if adj[v][u]: dp2[(S|(1<<u))*n+u] += val
        full = (1<<n)-1
        h = sum(dp2[full*n+v] for v in range(n))
        classes2.append({'h': h, 'size': len(members), 'aut': factorial(n)//len(members),
                        'scores': tuple(sorted(sum(adj[i]) for i in range(n)))})

    N = len(classes2)

    for bits in range(1 << m):
        c1 = bit_to_canon[bits]
        id1 = canon_to_id[c1]
        for bit in range(m):
            c2 = bit_to_canon[bits ^ (1 << bit)]
            id2 = canon_to_id[c2]
            if id1 != id2:
                edges.add((min(id1,id2), max(id1,id2)))
                adj_list[id1].add(id2)
                adj_list[id2].add(id1)

    return N, edges, adj_list, classes2

print("=" * 72)
print("  BRAND NEW SEQUENCES FROM G_n")
print("  opus-2026-03-22-S177")
print("=" * 72)
print()

# Build G_n for n=3,4,5
results = {}
for n in [3, 4, 5]:
    print(f"Building G_{n}...")
    N, edges, adj_list, classes = build_Gn(n)
    results[n] = (N, edges, adj_list, classes)
    print(f"  G_{n}: {N} nodes, {len(edges)} edges")

print()

# ==========================================================================
# COMPUTE EVERY GRAPH INVARIANT
# ==========================================================================

all_seqs = {}

for n in [3, 4, 5]:
    N, edges, adj_list, classes = results[n]
    degrees = [len(adj_list[i]) for i in range(N)]

    # 1. Nodes (=iso classes, A000568)
    all_seqs.setdefault('nodes', []).append(N)

    # 2. Edges
    all_seqs.setdefault('edges', []).append(len(edges))

    # 3. Max degree
    all_seqs.setdefault('max_degree', []).append(max(degrees))

    # 4. Min degree
    all_seqs.setdefault('min_degree', []).append(min(degrees))

    # 5. Sum of degrees = 2|E|
    all_seqs.setdefault('sum_degrees', []).append(sum(degrees))

    # 6. Number of self-complementary classes
    n_sc = sum(1 for cl in classes if cl['scores'] == tuple(sorted(n-1-s for s in cl['scores'])))
    all_seqs.setdefault('sc_classes', []).append(n_sc)

    # 7. Diameter
    diam = 0
    for i in range(N):
        dist = [-1]*N
        dist[i] = 0
        q = deque([i])
        while q:
            u = q.popleft()
            for v in adj_list[u]:
                if dist[v] == -1:
                    dist[v] = dist[u] + 1
                    q.append(v)
        diam = max(diam, max(d for d in dist if d >= 0))
    all_seqs.setdefault('diameter', []).append(diam)

    # 8. Number of triangles in G_n
    n_triangles = 0
    for (i, j) in edges:
        common = adj_list[i] & adj_list[j]
        n_triangles += len(common)
    n_triangles //= 3
    all_seqs.setdefault('triangles', []).append(n_triangles)

    # 9. Independence number (max independent set)
    alpha = 0
    for mask in range(1 << N):
        sel = [i for i in range(N) if mask & (1 << i)]
        indep = True
        for a in range(len(sel)):
            for b in range(a+1, len(sel)):
                if sel[b] in adj_list[sel[a]]:
                    indep = False
                    break
            if not indep: break
        if indep:
            alpha = max(alpha, len(sel))
    all_seqs.setdefault('independence_number', []).append(alpha)

    # 10. Chromatic number (greedy upper bound)
    # Use greedy coloring
    colors = [-1]*N
    order = sorted(range(N), key=lambda i: -degrees[i])
    for v in order:
        used = {colors[u] for u in adj_list[v] if colors[u] >= 0}
        c = 0
        while c in used:
            c += 1
        colors[v] = c
    chromatic_ub = max(colors) + 1
    all_seqs.setdefault('chromatic_ub', []).append(chromatic_ub)

    # 11. Number of H-level edges (same H)
    level = sum(1 for (i,j) in edges if classes[i]['h'] == classes[j]['h'])
    all_seqs.setdefault('level_edges', []).append(level)

    # 12. Number of H values (distinct H)
    n_h = len(set(cl['h'] for cl in classes))
    all_seqs.setdefault('h_values', []).append(n_h)

    # 13. Meta-H = I(G_n, 2)
    meta_coeffs = [0]*(N+1)
    for mask in range(1 << N):
        sel = [i for i in range(N) if mask & (1 << i)]
        indep = True
        for a in range(len(sel)):
            for b in range(a+1, len(sel)):
                if sel[b] in adj_list[sel[a]]:
                    indep = False
                    break
            if not indep: break
        if indep:
            meta_coeffs[len(sel)] += 1
    meta_H = sum(c * 2**k for k, c in enumerate(meta_coeffs))
    all_seqs.setdefault('meta_H', []).append(meta_H)

    # 14. Total |Aut| across all classes
    total_aut = sum(cl['aut'] for cl in classes)
    all_seqs.setdefault('total_aut', []).append(total_aut)

    # 15. Product of class sizes
    prod_sizes = 1
    for cl in classes:
        prod_sizes *= cl['size']
    all_seqs.setdefault('prod_sizes', []).append(prod_sizes)

    # 16. Number of score-ambiguous classes
    score_groups = defaultdict(list)
    for i, cl in enumerate(classes):
        score_groups[cl['scores']].append(i)
    ambig = sum(len(ids) for ids in score_groups.values() if len(ids) > 1)
    all_seqs.setdefault('score_ambiguous', []).append(ambig)

    # 17. Girth of G_n
    girth = float('inf')
    for start in range(N):
        dist = [-1]*N
        dist[start] = 0
        parent = [-1]*N
        q = deque([start])
        while q:
            u = q.popleft()
            for v in adj_list[u]:
                if dist[v] == -1:
                    dist[v] = dist[u] + 1
                    parent[v] = u
                    q.append(v)
                elif parent[u] != v:
                    girth = min(girth, dist[u] + dist[v] + 1)
    all_seqs.setdefault('girth', []).append(girth if girth < float('inf') else 0)

    # 18. Number of maximal cliques
    # Simple: count cliques
    max_clique = 0
    for mask in range(1 << N):
        sel = [i for i in range(N) if mask & (1 << i)]
        is_clique = True
        for a in range(len(sel)):
            for b in range(a+1, len(sel)):
                if sel[b] not in adj_list[sel[a]]:
                    is_clique = False
                    break
            if not is_clique: break
        if is_clique:
            max_clique = max(max_clique, len(sel))
    all_seqs.setdefault('clique_number', []).append(max_clique)

    # 19. Wiener index (sum of all pairwise distances)
    wiener = 0
    for i in range(N):
        dist = [-1]*N
        dist[i] = 0
        q = deque([i])
        while q:
            u = q.popleft()
            for v in adj_list[u]:
                if dist[v] == -1:
                    dist[v] = dist[u] + 1
                    q.append(v)
        wiener += sum(d for d in dist if d > 0)
    wiener //= 2
    all_seqs.setdefault('wiener_index', []).append(wiener)

    # 20. Skip perfect matchings (too expensive for brute force at n=5)
    all_seqs.setdefault('perfect_matchings', []).append(-1)

print()
print("=" * 72)
print("  ALL SEQUENCES FROM G_n (n=3,4,5)")
print("=" * 72)
print()

print(f"  {'Sequence':25s}  {'n=3':>6s}  {'n=4':>6s}  {'n=5':>6s}  notes")
print(f"  {'─'*25}  {'─'*6}  {'─'*6}  {'─'*6}  {'─'*30}")

for name, vals in sorted(all_seqs.items()):
    v3 = str(vals[0]) if len(vals) > 0 else '?'
    v4 = str(vals[1]) if len(vals) > 1 else '?'
    v5 = str(vals[2]) if len(vals) > 2 else '?'

    note = ""
    if name == 'nodes': note = "A000568 (known)"
    elif name == 'h_values': note = "from S165"
    elif name == 'meta_H': note = "I(G_n, 2)"
    elif name == 'sc_classes': note = "A000171 related"

    print(f"  {name:25s}  {v3:>6s}  {v4:>6s}  {v5:>6s}  {note}")

print()

# Flag potentially new sequences
print(f"  POTENTIALLY NEW SEQUENCES (3 terms, check OEIS):")
for name, vals in sorted(all_seqs.items()):
    if name in ['nodes', 'h_values', 'sc_classes']:
        continue  # known
    print(f"    {name}: {vals}")

print()
print("opus-2026-03-22-S177")
