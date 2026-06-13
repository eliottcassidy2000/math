#!/usr/bin/env python3
"""
even_graph_metagraph_s315.py — The Even Graph Metagraph as First-Class Object
opus-2026-03-24-S315

THE SETUP:
- Fix K_n and a spanning tree T_0 (the base path P_0 = n-1 → n-2 → ... → 0)
- The cycle space of K_n has dimension m = C(n-1,2) over GF(2)
- Fundamental cycles (one per non-tree edge) form a basis
- Each tiling = binary vector in {0,1}^m = which fundamental cycles to XOR
- Result: an EVEN SUBGRAPH of K_n (all degrees even)

THE METAGRAPH E_n:
- Vertices = isomorphism classes of even graphs on n vertices
- Edges = connected by single fundamental-cycle XOR (= single tile flip)

THIS IS THE PARALLEL UNIVERSE to G_n (tournament metagraph).
Same labeled structure Q_m, different quotient by S_n.

INVARIANTS TO COMPUTE:
- V, E of E_n
- Degree sequence
- Chromatic number χ(E_n)
- Clique number ω(E_n)  
- Spectrum
- Even-graph invariants per class: |E|, degree seq, components, triangles, |Aut|
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

def build_even_graph_metagraph(n):
    """Build the even graph metagraph E_n."""
    # All edges of K_n
    edges = [(i,j) for i in range(n) for j in range(i+1, n)]
    m_full = len(edges)  # C(n,2)
    edge_idx = {e: i for i, e in enumerate(edges)}
    
    # Base path (spanning tree): 0-1, 1-2, ..., (n-2)-(n-1)
    base_edges = [(i, i+1) for i in range(n-1)]
    non_base = [e for e in edges if e not in base_edges]
    m = len(non_base)  # C(n-1,2)
    
    # Fundamental cycles: for non-base edge (i,j), cycle = path i→i+1→...→j plus (i,j)
    fund_cycles = []
    for e in non_base:
        i, j = e
        cycle_bits = 0
        for k in range(i, j):
            cycle_bits |= (1 << edge_idx[(k, k+1)])
        cycle_bits |= (1 << edge_idx[(i, j)])
        fund_cycles.append(cycle_bits)
    
    # Map tiling → even graph (as bit vector on edges of K_n)
    def tiling_to_even(tile_bits):
        result = 0
        for k in range(m):
            if tile_bits & (1 << k):
                result ^= fund_cycles[k]
        return result
    
    # Graph properties from bit vector
    def graph_degree_seq(g_bits):
        deg = [0] * n
        for k, (i,j) in enumerate(edges):
            if g_bits & (1 << k):
                deg[i] += 1; deg[j] += 1
        return tuple(sorted(deg))
    
    def graph_edge_count(g_bits):
        return bin(g_bits).count('1')
    
    def graph_triangles(g_bits):
        count = 0
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if (g_bits & (1 << edge_idx[(i,j)])) and \
                       (g_bits & (1 << edge_idx[(i,k)])) and \
                       (g_bits & (1 << edge_idx[(j,k)])):
                        count += 1
        return count
    
    def graph_components(g_bits):
        adj = defaultdict(set)
        for k, (i,j) in enumerate(edges):
            if g_bits & (1 << k):
                adj[i].add(j); adj[j].add(i)
        visited = set()
        comps = 0
        for v in range(n):
            if v not in visited:
                comps += 1
                stack = [v]
                while stack:
                    u = stack.pop()
                    if u in visited: continue
                    visited.add(u)
                    for w in adj[u]:
                        if w not in visited:
                            stack.append(w)
        return comps
    
    # Canonical form under S_n (for even graphs)
    all_perms = list(permutations(range(n)))
    def graph_canon(g_bits):
        best = None
        for p in all_perms:
            nb = 0
            for k, (i,j) in enumerate(edges):
                pi, pj = min(p[i], p[j]), max(p[i], p[j])
                nk = edge_idx[(pi, pj)]
                if g_bits & (1 << k):
                    nb |= (1 << nk)
            if best is None or nb < best:
                best = nb
        return best
    
    # Enumerate all 2^m even graphs via tilings
    class_map = {}  # canon → class_id
    class_info = {}  # class_id → properties
    tiling_class = {}  # tiling → class_id
    class_id = 0
    
    for tb in range(1 << m):
        eg = tiling_to_even(tb)
        cn = graph_canon(eg)
        if cn not in class_map:
            class_map[cn] = class_id
            ds = graph_degree_seq(eg)
            ec = graph_edge_count(eg)
            tri = graph_triangles(eg)
            comp = graph_components(eg)
            class_info[class_id] = {
                'canon': cn, 'deg_seq': ds, 'edges': ec,
                'triangles': tri, 'components': comp,
                'fiber': 0
            }
            class_id += 1
        cid = class_map[cn]
        class_info[cid]['fiber'] += 1
        tiling_class[tb] = cid
    
    V = class_id
    
    # Compute |Aut| from orbit-stabilizer: |orbit| = n!/|Aut|
    for cid in range(V):
        fiber = class_info[cid]['fiber']
        # fiber = number of tilings mapping to this class
        # but orbit size = number of labeled even graphs in this class
        # orbit = fiber (since each tiling gives a distinct labeled graph)
        # Actually: orbit = fiber? Let me verify...
        # The tiling→even graph map is injective on labeled graphs.
        # But not all labeled graphs in the iso class need be in the image
        # (only those containing the base path as cycle-space generators).
        # Actually, the cycle space map IS surjective on even graphs (proven above).
        # So fiber = number of labeled even graphs in this iso class = orbit size.
        orbit = fiber
        aut = factorial(n) // orbit if orbit > 0 else factorial(n)
        class_info[cid]['aut'] = aut
        class_info[cid]['orbit'] = orbit
    
    # Build metagraph edges: tile flip connects classes
    meta_edges = set()
    for tb in range(1 << m):
        cid_a = tiling_class[tb]
        for tile in range(m):
            tb2 = tb ^ (1 << tile)
            cid_b = tiling_class[tb2]
            if cid_a != cid_b:
                meta_edges.add((min(cid_a, cid_b), max(cid_a, cid_b)))
    
    adj_list = defaultdict(set)
    for (a, b) in meta_edges:
        adj_list[a].add(b); adj_list[b].add(a)
    
    return V, class_info, meta_edges, adj_list, tiling_class, m

def dsatur_color(V, adj):
    """DSATUR greedy coloring."""
    colors = [-1] * V
    sat = [0] * V
    adj_colors = [set() for _ in range(V)]
    degrees = [len(adj[v]) for v in range(V)]
    v0 = max(range(V), key=lambda v: degrees[v])
    colors[v0] = 0
    for u in adj[v0]:
        adj_colors[u].add(0); sat[u] = len(adj_colors[u])
    for step in range(1, V):
        best = -1; bs = -1; bd = -1
        for v in range(V):
            if colors[v] >= 0: continue
            if sat[v] > bs or (sat[v] == bs and degrees[v] > bd):
                best = v; bs = sat[v]; bd = degrees[v]
        v = best; used = adj_colors[v]
        c = 0
        while c in used: c += 1
        colors[v] = c
        for u in adj[v]:
            if colors[u] < 0:
                adj_colors[u].add(c); sat[u] = len(adj_colors[u])
    return max(colors) + 1, colors

def find_max_clique(V, adj, edges):
    """Bron-Kerbosch with pivoting."""
    best = []
    count = [0]
    def bk(R, P, X):
        nonlocal best
        count[0] += 1
        if count[0] > 5000000: return
        if not P and not X:
            if len(R) > len(best): best = list(R)
            return
        if len(R) + len(P) <= len(best): return
        pivot_cands = P | X
        if not pivot_cands: return
        pivot = max(pivot_cands, key=lambda u: len(adj[u] & P))
        for v in list(P - adj[pivot]):
            bk(R | {v}, P & adj[v], X & adj[v])
            P.remove(v); X.add(v)
    bk(set(), set(range(V)), set())
    return best

print("=" * 72)
print("  THE EVEN GRAPH METAGRAPH E_n")
print("  opus-2026-03-24-S315")
print("=" * 72)

for n in range(3, 8):
    print(f"\n{'#'*72}")
    print(f"  n = {n}")
    print(f"{'#'*72}")
    
    if n >= 8:
        print("  Skipping (too large for brute-force canonicalization)")
        continue
    
    V, ci, edges, adj, tc, m = build_even_graph_metagraph(n)
    E = len(edges)
    
    print(f"  V(E_{n}) = {V} even graph iso classes")
    print(f"  E(E_{n}) = {E} edges")
    print(f"  m = {m} tiles (cycle-space dimension)")
    
    # Degree stats
    degrees = [len(adj[v]) for v in range(V)]
    if V > 0:
        print(f"  Degree: min={min(degrees)}, max={max(degrees)}, mean={np.mean(degrees):.1f}")
    
    # Class details
    print(f"\n  EVEN GRAPH CLASSES:")
    print(f"  {'ID':>3} {'|E|':>4} {'deg_seq':>20} {'tri':>4} {'comp':>5} {'|Aut|':>6} {'fiber':>6}")
    for cid in range(V):
        c = ci[cid]
        print(f"  {cid:>3} {c['edges']:>4} {str(c['deg_seq']):>20} {c['triangles']:>4} {c['components']:>5} {c['aut']:>6} {c['fiber']:>6}")
    
    # Chromatic number and clique number
    if V > 1:
        clique = find_max_clique(V, adj, edges)
        omega = len(clique)
        chi_upper, dsatur_cols = dsatur_color(V, adj)
        print(f"\n  ω(E_{n}) = {omega} (max clique)")
        print(f"  χ(E_{n}) ≤ {chi_upper} (DSATUR upper bound)")
        
        # Try to find exact χ via backtracking for small V
        if V <= 50:
            for k in range(omega, chi_upper):
                cols = [-1] * V
                order = sorted(range(V), key=lambda v: -len(adj[v]))
                count = [0]
                def bt(idx):
                    if count[0] > 5000000: return False
                    if idx == V: return True
                    count[0] += 1
                    v = order[idx]
                    used = {cols[u] for u in adj[v] if cols[u] >= 0}
                    for c in range(k):
                        if c not in used:
                            cols[v] = c
                            if bt(idx + 1): return True
                            cols[v] = -1
                    return False
                if bt(0):
                    chi_upper = k
                    print(f"  χ(E_{n}) = {k} (exact)")
                    break
        
        print(f"  n-1 = {n-1}")
        if chi_upper == n-1:
            print(f"  *** χ(E_{n}) = n-1 = {n-1}! Same as tournament metagraph! ***")
        elif chi_upper < n-1:
            print(f"  *** χ(E_{n}) = {chi_upper} < n-1 = {n-1}! DIFFERENT from tournaments! ***")
        else:
            print(f"  χ(E_{n}) = {chi_upper} > n-1 = {n-1}")
    
    # Spectrum
    if V > 1:
        A = np.zeros((V, V))
        for (a, b) in edges: A[a][b] = 1; A[b][a] = 1
        evals = sorted(np.linalg.eigvalsh(A))
        print(f"\n  SPECTRUM:")
        print(f"  λ_max = {evals[-1]:.4f}, λ_min = {evals[0]:.4f}")
        if abs(evals[0]) > 1e-10:
            hoffman = 1 + evals[-1] / abs(evals[0])
            print(f"  Hoffman bound: {hoffman:.4f}")
    
    # Fiber distribution (how many tilings per class)
    fibers = sorted([ci[cid]['fiber'] for cid in range(V)])
    print(f"\n  FIBER DISTRIBUTION:")
    print(f"  Total tilings: {sum(fibers)} = 2^{m}")
    print(f"  Fibers: {fibers}")

# COMPARISON TABLE
print(f"\n{'='*72}")
print(f"  COMPARISON: TOURNAMENT G_n vs EVEN GRAPH E_n")
print(f"{'='*72}")
print(f"  {'n':>3} {'V(G_n)':>8} {'V(E_n)':>8} {'E(G_n)':>8} {'E(E_n)':>8} {'χ(G_n)':>7} {'χ(E_n)':>7}")

# We know tournament metagraph from previous work
tourn_data = {
    3: (1, 0),  # V, E for merged
    4: (3, 3),
    5: (10, 21),  # merged = (A000568(5)+SC)/2 but using unmerged here
    6: (34, 143),
    7: (272, 2123),
}

print("DONE.")
