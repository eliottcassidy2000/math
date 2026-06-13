#!/usr/bin/env python3
"""
overnight_round2_s290.py — Graph tournaments at n=8 (Q_3) and n=9 (Paley9)
opus-2026-03-24-S290

Construct T_Q3 (3-cube graph tournament, n=8) and T_Paley9 (n=9).
"""

import sys
from math import comb, factorial
from collections import defaultdict

sys.stdout.reconfigure(line_buffering=True)

def find_ham_paths_graph(adj, n, max_paths=500):
    """Find Hamiltonian paths in a graph given by adjacency dict."""
    paths = []
    def dfs(path, visited):
        if len(paths) >= max_paths: return
        if len(path) == n:
            paths.append(tuple(path))
            return
        last = path[-1]
        for nb in sorted(adj[last]):
            if nb not in visited:
                visited.add(nb)
                path.append(nb)
                dfs(path, visited)
                path.pop()
                visited.remove(nb)
    for start in range(n):
        if len(paths) >= max_paths: break
        dfs([start], {start})
    return paths

def build_graph_tournament(ham_path, graph_edges, n):
    """Build tournament from graph + Hamiltonian path."""
    label = {}
    for k, v in enumerate(ham_path):
        label[v] = n - k

    adj = [[0]*(n+1) for _ in range(n+1)]
    for i in range(n):
        for j in range(i+1, n):
            is_edge = (min(i,j), max(i,j)) in graph_edges
            li = label[i]; lj = label[j]
            a = max(li, lj); b = min(li, lj)
            if is_edge:
                adj[a][b] = 1  # unflipped
            else:
                adj[b][a] = 1  # flipped
    return adj, label

def H_count(adj, n):
    dp = {}
    for v in range(1, n+1): dp[(1 << (v-1), v)] = 1
    full = (1 << n) - 1
    for S in range(1, 1 << n):
        for v in range(1, n+1):
            if not (S & (1 << (v-1))): continue
            val = dp.get((S, v), 0)
            if val == 0: continue
            for u in range(1, n+1):
                if S & (1 << (u-1)): continue
                if adj[v][u]:
                    dp[(S | (1 << (u-1)), u)] = dp.get((S | (1 << (u-1)), u), 0) + val
    return sum(dp.get((full, v), 0) for v in range(1, n+1))

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i][j] for j in range(1, n+1) if j != i) for i in range(1, n+1)))

# ============================================================
# Q_3: 3-cube graph on 8 vertices
# ============================================================

print("=" * 72)
print("  GRAPH TOURNAMENTS: Q_3 (n=8) and Paley(9)")
print("=" * 72)

# Q_3: vertices = binary strings of length 3
# Edge iff Hamming distance 1
q3_edges = set()
for i in range(8):
    for j in range(i+1, 8):
        if bin(i ^ j).count('1') == 1:
            q3_edges.add((i, j))

print(f"\n  Q_3: {len(q3_edges)} edges on 8 vertices")
print(f"  Edges: {sorted(q3_edges)}")

# Hamiltonian paths in Q_3
print(f"  Finding Hamiltonian paths...")
q3_adj = defaultdict(set)
for (a,b) in q3_edges: q3_adj[a].add(b); q3_adj[b].add(a)
q3_paths = find_ham_paths_graph(q3_adj, 8, max_paths=500)
print(f"  Found {len(q3_paths)} directed paths (capped at 500)")

# Build tournaments from different paths
unique_q3 = {}
for path in q3_paths:
    adj, label = build_graph_tournament(path, q3_edges, 8)
    h = H_count(adj, 8)
    ss = score_seq(adj, 8)
    key = (ss, h)
    if key not in unique_q3:
        unique_q3[key] = {'count': 0, 'path': path, 'H': h, 'scores': ss}
    unique_q3[key]['count'] += 1

print(f"\n  {len(unique_q3)} DISTINCT tournament types from Q_3")
for key, data in sorted(unique_q3.items(), key=lambda x: x[1]['H']):
    ss = data['scores']
    palindromic = (ss == tuple(7-s for s in reversed(ss)))
    szele = factorial(8) / (2**7)
    print(f"    H={data['H']}, scores={ss}, "
          f"#paths={data['count']}, palindromic={palindromic}, "
          f"H/Szele={data['H']/szele:.4f}")

# ============================================================
# Paley graph of order 9
# ============================================================

# Paley(9): vertices = GF(9) = {0,1,2,α,α+1,α+2,2α,2α+1,2α+2}
# where α² = α+1 (primitive element of GF(9))
# Edge iff difference is a nonzero square in GF(9)
# GF(9) squares: {0, 1, 2, α+1, 2α+2} (the nonzero squares are {1, 2, α+1, 2α+2})
# Actually easier: Paley(9) = strongly regular (9,4,1,2)

# Paley(9) edges: x~y iff x-y is a nonzero square in GF(9)
# GF(9) = F_3[x]/(x²+1) where x² = -1 = 2
# Elements: 0, 1, 2, x, x+1, x+2, 2x, 2x+1, 2x+2
# Label: 0→0, 1→1, 2→2, x→3, x+1→4, x+2→5, 2x→6, 2x+1→7, 2x+2→8

# Squares in GF(9)^*:
# 1²=1, 2²=1 (in F_3), x²=2, (x+1)²=x²+2x+1=2+2x+1=2x,
# (x+2)²=x²+x+1=2+x+1=x, (2x)²=x²=2, (2x+1)²=(2x)²+2·2x+1=2+x+1=x+2... hmm
# Let me just use the known fact: Paley(9) is the unique srg(9,4,1,2).
# It's also the Cayley graph of Z_3 × Z_3 with connection set = QR.

# Paley(9) adjacency from standard construction:
# I'll use the known edge list for srg(9,4,1,2)
# Each vertex has degree 4, so 9×4/2 = 18 edges

# The Paley graph of order 9 can be constructed as:
# V = {0,...,8}, edge iff (i-j) mod 9 is a quadratic residue
# Wait, 9 is not prime. Paley(9) uses GF(9), not Z_9.
# For GF(9): QR = {1, 2, x+1, 2x+2} (4 nonzero squares)
# This is 4-regular: each vertex connects to 4 others.

# Let me just construct it explicitly.
# Using the known Paley(9) = L(K_{3,3}) (line graph of K_{3,3})
# Actually Paley(9) is NOT L(K_{3,3}). Let me use a direct construction.

# GF(9) = F_3[α]/(α²+1), so α² = -1 = 2.
# Elements: 0, 1, 2, α, α+1, α+2, 2α, 2α+1, 2α+2
# Mapping: 0→0, 1→1, 2→2, α→3, α+1→4, α+2→5, 2α→6, 2α+1→7, 2α+2→8

# Multiplication table needed for squares...
# Actually, let me just use a simpler construction.
# For Paley(q) with q=p^k≡1(mod 4): edge iff difference is a nonzero square.
# q=9, 9≡1(mod 4). GF(9) squares = {x² : x ∈ GF(9)*}.

# GF(9) = F_3[α] with α²=2 (or α²=-1 since -1=2 in F_3).
# The nonzero elements: 1, 2, α, α+1, α+2, 2α, 2α+1, 2α+2
# Squares: 1²=1, 2²=1, α²=2, (α+1)²=α²+2α+1=2+2α+1=2α,
#   (α+2)²=α²+α+1=2+α+1=α, (2α)²=4α²=α²=2,
#   (2α+1)²=4α²+4α+1=α²+α+1=2+α+1=α, wait that's same as (α+2)²...
# Let me be more careful:
# (2α+1)² = (2α)² + 2·(2α)·1 + 1² = α² + 4α + 1 = 2 + α + 1 = α (in F_3)
# Hmm, α = (α+2)² and α = (2α+1)². So squares of {1,2,α,α+1,α+2,2α,2α+1,2α+2} are:
# {1,1,2,2α,α,2,α,...} = {1, 2, α, 2α}
# So QR = {1, 2, α, 2α} = {1, 2, 3, 6} in our labeling.

paley9_QR = {1, 2, 3, 6}  # quadratic residues of GF(9)

# GF(9) addition: we need to compute i-j in GF(9) for all pairs
# Elements as (a,b) in F_3²: 0=(0,0), 1=(1,0), 2=(2,0), 3=(0,1), 4=(1,1), 5=(2,1), 6=(0,2), 7=(1,2), 8=(2,2)

gf9 = [(i % 3, i // 3) for i in range(9)]

def gf9_sub(i, j):
    """Compute i - j in GF(9) = F_3²."""
    ai, bi = gf9[i]; aj, bj = gf9[j]
    return ((ai - aj) % 3) + ((bi - bj) % 3) * 3

paley9_edges = set()
for i in range(9):
    for j in range(i+1, 9):
        diff = gf9_sub(i, j)
        if diff in paley9_QR:
            paley9_edges.add((i, j))
        # Also check reverse
        diff_rev = gf9_sub(j, i)
        if diff_rev in paley9_QR:
            paley9_edges.add((i, j))

print(f"\n  Paley(9): {len(paley9_edges)} edges on 9 vertices")
# Should be 18 edges for srg(9,4,1,2)
degs = [0]*9
for (a,b) in paley9_edges: degs[a] += 1; degs[b] += 1
print(f"  Degrees: {degs}")
print(f"  4-regular? {all(d == 4 for d in degs)}")

if len(paley9_edges) == 18:
    print(f"  Finding Hamiltonian paths...")
    p9_adj = defaultdict(set)
    for (a,b) in paley9_edges: p9_adj[a].add(b); p9_adj[b].add(a)
    p9_paths = find_ham_paths_graph(p9_adj, 9, max_paths=500)
    print(f"  Found {len(p9_paths)} directed paths (capped at 500)")

    unique_p9 = {}
    for path in p9_paths[:100]:  # compute H for first 100
        adj, label = build_graph_tournament(path, paley9_edges, 9)
        h = H_count(adj, 9)
        ss = score_seq(adj, 9)
        key = (ss, h)
        if key not in unique_p9:
            unique_p9[key] = {'count': 0, 'path': path, 'H': h, 'scores': ss}
        unique_p9[key]['count'] += 1

    print(f"\n  {len(unique_p9)} DISTINCT tournament types from Paley(9) (from {min(100,len(p9_paths))} paths)")
    for key, data in sorted(unique_p9.items(), key=lambda x: x[1]['H']):
        ss = data['scores']
        palindromic = (ss == tuple(8-s for s in reversed(ss)))
        szele = factorial(9) / (2**8)
        print(f"    H={data['H']}, scores={ss}, "
              f"#paths={data['count']}, palindromic={palindromic}, "
              f"H/Szele={data['H']/szele:.4f}")

print("\nDONE.")
