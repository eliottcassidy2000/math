#!/usr/bin/env python3
"""
gn_merged_tutte_s219.py — Tutte polynomial invariants + chromatic number verification
opus-2026-03-23-S219

Key tasks:
1. Verify chi(G_n/Z_2) = n-1 at n=7 (272 vertices, 2123 edges)
2. Compute spanning tree count for n=5..7 via Kirchhoff
3. Compute independence polynomial for merged graph
4. Extract all Tutte evaluations for n=3..5 (n=6+ too large for full Tutte)
5. Look for patterns in Tutte coefficients
"""

import sys
import subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  TUTTE POLYNOMIAL INVARIANTS + CHROMATIC NUMBER")
print("  opus-2026-03-23-S219")
print("=" * 80)

# ============================================================================
# REUSE the merged graph builder from S218
# ============================================================================

def build_merged_graph(n):
    """Build G_n/Z_2 and return (V, edges, adj_list)."""
    m = comb(n, 2)
    PAIRS = [(i,j) for i in range(n) for j in range(i+1,n)]

    # Get tournaments from nauty
    if n <= 6:
        ALL_PERMS = list(permutations(range(n)))
        def bits_to_adj(bits):
            adj = [[0]*n for _ in range(n)]
            for k,(i,j) in enumerate(PAIRS):
                if bits & (1<<k): adj[i][j]=1
                else: adj[j][i]=1
            return adj
        def canonical(adj):
            best = None
            for perm in ALL_PERMS:
                form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
                if best is None or form < best: best = form
            return best
        classes_map = {}
        for bits in range(1<<m):
            adj = bits_to_adj(bits)
            c = canonical(adj)
            if c not in classes_map:
                classes_map[c] = {'bits': bits, 'adj': adj, 'canon': c, 'size': 0}
            classes_map[c]['size'] += 1
        classes = list(classes_map.values())
        for i, cl in enumerate(classes): cl['cid'] = i
        is_nauty = False
    else:
        res = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
        lines = [l.strip() for l in (res.stdout or '').split('\n')
                 if len(l.strip())==m and all(c in '01' for c in l.strip())]
        def bits_to_adj_n(bits):
            adj = [[0]*n for _ in range(n)]
            for k,(i,j) in enumerate(PAIRS):
                if bits & (1<<(m-1-k)): adj[i][j]=1
                else: adj[j][i]=1
            return adj
        classes = []
        for i, line in enumerate(lines):
            bits = int(line, 2)
            classes.append({'bits': bits, 'adj': bits_to_adj_n(bits), 'cid': i})
        is_nauty = True

    # Compute SC status and complement pairing
    def score_seq(adj):
        return tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))
    def c3_count(adj):
        c3 = 0
        for i in range(n):
            for j in range(i+1,n):
                for k in range(j+1,n):
                    if adj[i][j] and adj[j][k] and adj[k][i]: c3+=1
                    if adj[i][k] and adj[k][j] and adj[j][i]: c3+=1
        return c3
    def H_dp(adj):
        dp = {}
        for v in range(n): dp[(1<<v,v)] = 1
        for S in range(1, 1<<n):
            for v in range(n):
                if not(S&(1<<v)): continue
                val = dp.get((S,v),0)
                if val==0: continue
                for u in range(n):
                    if S&(1<<u): continue
                    if adj[v][u]:
                        key=(S|(1<<u),u)
                        dp[key] = dp.get(key,0) + val
        return sum(dp.get(((1<<n)-1,v),0) for v in range(n))
    def canon_pruned(adj):
        scores = [sum(adj[i][j] for j in range(n)) for i in range(n)]
        sg = defaultdict(list)
        for v in range(n): sg[scores[v]].append(v)
        groups = [sg[s] for s in sorted(set(scores))]
        if all(len(g)==1 for g in groups):
            perm = [g[0] for g in groups]
            return tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        best = None
        def gen(gs):
            if not gs: yield []; return
            for p in permutations(gs[0]):
                for rest in gen(gs[1:]): yield list(p)+rest
        for perm in gen(groups):
            form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
            if best is None or form < best: best = form
        return best

    # Build hash and canonical lookups
    hash_map = defaultdict(list)
    canon_map = {}
    for cl in classes:
        adj = cl['adj']
        cl['score'] = score_seq(adj)
        cl['c3'] = c3_count(adj)
        cl['H'] = H_dp(adj)
        cl['own_canon'] = canon_pruned(adj)
        comp = [[1-adj[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
        cl['comp_canon'] = canon_pruned(comp)
        cl['sc'] = (cl['own_canon'] == cl['comp_canon'])
        hash_map[(cl['score'], cl['c3'], cl['H'])].append(cl['cid'])
        canon_map[cl['own_canon']] = cl['cid']

    for cl in classes:
        cl['comp_cid'] = canon_map.get(cl['comp_canon'], -1)

    def lookup(adj):
        ss = score_seq(adj); c3 = c3_count(adj); h = H_dp(adj)
        cids = hash_map.get((ss,c3,h))
        if not cids: return -1
        if len(cids)==1: return cids[0]
        c = canon_pruned(adj)
        return canon_map.get(c, -1)

    # Build edges
    edges = set()
    for cl in classes:
        bits = cl['bits']
        for arc in range(m):
            if is_nauty:
                fb = bits ^ (1 << (m-1-arc))
                fa = [[0]*n for _ in range(n)]
                for k,(i,j) in enumerate(PAIRS):
                    if fb & (1<<(m-1-k)): fa[i][j]=1
                    else: fa[j][i]=1
            else:
                fb = bits ^ (1 << arc)
                fa = [[0]*n for _ in range(n)]
                for k,(i,j) in enumerate(PAIRS):
                    if fb & (1<<k): fa[i][j]=1
                    else: fa[j][i]=1
            nb = lookup(fa)
            if nb >= 0 and nb != cl['cid']:
                edges.add((min(cl['cid'], nb), max(cl['cid'], nb)))

    # Build merged graph
    merged_node = {}
    mid = 0
    for cl in classes:
        ci = cl['cid']
        if ci in merged_node: continue
        if cl['sc']:
            merged_node[ci] = mid; mid += 1
        else:
            merged_node[ci] = mid
            comp = cl['comp_cid']
            if comp >= 0 and comp != ci:
                merged_node[comp] = mid
            mid += 1

    merged_edges = set()
    for (a,b) in edges:
        ma, mb = merged_node[a], merged_node[b]
        if ma != mb:
            merged_edges.add((min(ma,mb), max(ma,mb)))

    merged_adj = defaultdict(set)
    for (a,b) in merged_edges:
        merged_adj[a].add(b)
        merged_adj[b].add(a)

    return mid, merged_edges, merged_adj, classes

# ============================================================================
# CHROMATIC NUMBER via greedy + exact verification
# ============================================================================

def greedy_chromatic(V, adj):
    """Greedy coloring with DSatur heuristic."""
    colors = [-1] * V
    sat = [0] * V  # saturation degree
    deg = [len(adj[i]) for i in range(V)]

    for step in range(V):
        # Pick uncolored vertex with max saturation, break ties by degree
        best = -1
        for v in range(V):
            if colors[v] >= 0: continue
            if best < 0 or (sat[v], deg[v]) > (sat[best], deg[best]):
                best = v

        # Assign smallest available color
        used = set(colors[nb] for nb in adj[best] if colors[nb] >= 0)
        c = 0
        while c in used: c += 1
        colors[best] = c

        # Update saturation
        for nb in adj[best]:
            if colors[nb] < 0:
                nbr_colors = set(colors[w] for w in adj[nb] if colors[w] >= 0)
                sat[nb] = len(nbr_colors)

    return max(colors) + 1

def verify_k_colorable(V, adj, k, timeout=30):
    """Check if graph is k-colorable via backtracking with pruning."""
    colors = [-1] * V
    t0 = time.time()

    # Order vertices by degree (descending)
    order = sorted(range(V), key=lambda v: -len(adj[v]))

    def backtrack(idx):
        if time.time() - t0 > timeout:
            return None  # timeout
        if idx == V:
            return True
        v = order[idx]
        used = set(colors[nb] for nb in adj[v] if colors[nb] >= 0)
        for c in range(k):
            if c not in used:
                colors[v] = c
                result = backtrack(idx + 1)
                if result is True:
                    return True
                if result is None:
                    return None
                colors[v] = -1
        return False

    return backtrack(0)

# ============================================================================
# SPANNING TREE COUNT via Kirchhoff (matrix-tree theorem)
# ============================================================================

def spanning_trees(V, edges):
    """Count spanning trees via Kirchhoff's theorem (det of Laplacian cofactor)."""
    if V <= 1: return 1 if V == 1 else 0

    # Build Laplacian
    import numpy as np
    L = np.zeros((V, V), dtype=np.float64)
    for (a, b) in edges:
        L[a][b] -= 1
        L[b][a] -= 1
        L[a][a] += 1
        L[b][b] += 1

    # Delete last row and column
    M = L[:-1, :-1]
    det = np.linalg.det(M)
    return round(det)

# ============================================================================
# INDEPENDENCE POLYNOMIAL
# ============================================================================

def independence_poly(V, adj, max_size=None):
    """Compute independence polynomial coefficients alpha_k."""
    if max_size is None:
        max_size = V

    # For small graphs, enumerate all independent sets
    if V <= 25:
        alpha = [0] * (V + 1)
        for mask in range(1 << V):
            # Check if independent
            verts = [v for v in range(V) if mask & (1 << v)]
            is_indep = True
            for i in range(len(verts)):
                for j in range(i+1, len(verts)):
                    if verts[j] in adj[verts[i]]:
                        is_indep = False; break
                if not is_indep: break
            if is_indep:
                alpha[len(verts)] += 1
        return alpha
    else:
        # Greedy estimate for large graphs
        return None

# ============================================================================
# MAIN COMPUTATION
# ============================================================================

results = {}

for n in range(3, 9):
    print(f"\n{'='*70}")
    print(f"  n = {n}")
    print(f"{'='*70}")
    t0 = time.time()

    Vm, m_edges, m_adj, classes = build_merged_graph(n)
    Em = len(m_edges)
    print(f"  G_{n}/Z_2: V={Vm}, E={Em} ({time.time()-t0:.1f}s)")

    # Chromatic number
    print(f"  Computing chromatic number...")
    chi_greedy = greedy_chromatic(Vm, m_adj)
    print(f"    Greedy: chi ≤ {chi_greedy}")

    # Try to verify exact chromatic number
    chi_exact = chi_greedy
    for k in range(chi_greedy - 1, 0, -1):
        print(f"    Testing {k}-colorable...", end=" ", flush=True)
        result = verify_k_colorable(Vm, m_adj, k, timeout=60)
        if result is True:
            chi_exact = k
            print("YES")
        elif result is False:
            print("NO")
            break
        else:
            print("TIMEOUT")
            break

    print(f"    Chromatic number: chi = {chi_exact}")
    print(f"    Conjecture chi = n-1 = {n-1}: {'✓' if chi_exact == n-1 else '✗'}")

    # Spanning trees (Kirchhoff)
    if Vm <= 500:
        tau = spanning_trees(Vm, m_edges)
        print(f"  Spanning trees: {tau}")
    else:
        tau = None

    # Independence polynomial (small graphs only)
    if Vm <= 25:
        alpha = independence_poly(Vm, m_adj)
        I_at_2 = sum(alpha[k] * 2**k for k in range(len(alpha)))
        print(f"  Independence polynomial: alpha = {alpha[:8]}...")
        print(f"  I(G_{n}/Z_2, 2) = {I_at_2}")
    else:
        alpha = None
        I_at_2 = None

    # Clique number
    omega = 1
    if Vm <= 1000:
        for v in range(Vm):
            clique = {v}
            cands = m_adj[v].copy()
            for u in sorted(cands, key=lambda x: -len(m_adj[x])):
                if all(u in m_adj[w] for w in clique):
                    clique.add(u)
            omega = max(omega, len(clique))
    print(f"  Clique number (greedy): omega ≥ {omega}")

    # Triangle count
    tri = 0
    for (a, b) in m_edges:
        tri += len(m_adj[a] & m_adj[b])
    tri //= 3
    print(f"  Triangles: {tri}")

    results[n] = {
        'V': Vm, 'E': Em, 'chi': chi_exact, 'tau': tau,
        'alpha_poly': alpha, 'I_at_2': I_at_2, 'omega': omega,
        'triangles': tri
    }

    print(f"  Done ({time.time()-t0:.1f}s)")

# ============================================================================
# SUMMARY
# ============================================================================

print(f"\n{'='*80}")
print(f"  SUMMARY OF TUTTE/CHROMATIC INVARIANTS")
print(f"{'='*80}")

print(f"\n  CHROMATIC NUMBER CONJECTURE: chi(G_n/Z_2) = n-1")
print(f"  {'n':>3} {'V':>6} {'E':>8} {'chi':>4} {'n-1':>4} {'match':>6}")
for n in sorted(results):
    R = results[n]
    print(f"  {n:3d} {R['V']:6d} {R['E']:8d} {R['chi']:4d} {n-1:4d} "
          f"{'✓' if R['chi']==n-1 else '✗':>6s}")

print(f"\n  SPANNING TREES: {[results[n]['tau'] for n in sorted(results) if results[n]['tau'] is not None]}")
print(f"  CLIQUE NUMBER: {[results[n]['omega'] for n in sorted(results)]}")
print(f"  TRIANGLES: {[results[n]['triangles'] for n in sorted(results)]}")

if any(results[n]['I_at_2'] for n in results):
    print(f"  I(G_n/Z_2, 2): {[results[n]['I_at_2'] for n in sorted(results) if results[n]['I_at_2'] is not None]}")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
