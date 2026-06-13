#!/usr/bin/env python3
"""
petersen_tournament_s288.py — The Petersen graph AS a tournament on 10 vertices
opus-2026-03-24-S288

CONSTRUCTION:
1. Start with the Petersen graph (10 vertices, 15 edges)
2. Find a Hamiltonian path through it (9 of the 15 edges)
3. The 9 path edges = base path arcs (unflipped)
4. The remaining 6 Petersen edges = unflipped tiles
5. The 30 non-Petersen edges = flipped tiles
6. Orient: unflipped = high→low, flipped = low→high

This defines a specific tournament T_Pet on 10 vertices.
Study its H, score sequence, |Aut|, iso class.
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict

sys.stdout.reconfigure(line_buffering=True)

# ============================================================
# THE PETERSEN GRAPH
# ============================================================

def petersen_edges():
    """Petersen graph: outer pentagon + inner pentagram + spokes."""
    edges = set()
    for i in range(5): edges.add((i, (i+1)%5))          # outer pentagon
    for i in range(5): edges.add((5+i, 5+(i+2)%5))      # inner pentagram
    for i in range(5): edges.add((i, 5+i))               # spokes
    return {(min(a,b), max(a,b)) for (a,b) in edges}

# ============================================================
# FIND ALL HAMILTONIAN PATHS IN THE PETERSEN GRAPH
# ============================================================

def find_ham_paths(edges, n=10):
    """Find all Hamiltonian paths in the Petersen graph."""
    adj = defaultdict(set)
    for (a,b) in edges:
        adj[a].add(b); adj[b].add(a)

    paths = []
    def dfs(path, visited):
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
        dfs([start], {start})

    return paths

# ============================================================
# BUILD THE PETERSEN TOURNAMENT
# ============================================================

def build_petersen_tournament(ham_path, pet_edges, n=10):
    """Build tournament from Petersen graph + Hamiltonian path.

    Labels: vertex visited k-th in path gets label (n-k+1).
    So first vertex = label n, last = label 1.
    Base path: n → n-1 → ... → 1.

    Petersen edges (not in path) = unflipped tiles (high→low).
    Non-Petersen edges = flipped tiles (low→high).
    """
    # Create labeling
    label = {}
    for k, v in enumerate(ham_path):
        label[v] = n - k  # first vertex = n, last = 1

    # Build adjacency matrix (1-indexed labels 1..n)
    adj = [[0]*(n+1) for _ in range(n+1)]

    for i in range(n):
        for j in range(i+1, n):
            # Original Petersen vertices
            pet_edge = (min(i,j), max(i,j)) in pet_edges
            li = label[i]; lj = label[j]
            a = max(li, lj); b = min(li, lj)

            if pet_edge:
                # Petersen edge → unflipped → high label beats low label
                adj[a][b] = 1
            else:
                # Non-Petersen → flipped → low label beats high label
                adj[b][a] = 1

    return adj, label

def tournament_properties(adj, n):
    """Compute properties of a tournament (1-indexed)."""
    # Score sequence
    scores = []
    for i in range(1, n+1):
        s = sum(adj[i][j] for j in range(1, n+1) if j != i)
        scores.append(s)
    score_seq = tuple(sorted(scores))

    # Hamiltonian path count
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
    H = sum(dp.get((full, v), 0) for v in range(1, n+1))

    # 3-cycle count
    c3 = 0
    for i in range(1, n+1):
        for j in range(i+1, n+1):
            for k in range(j+1, n+1):
                if adj[i][j] and adj[j][k] and adj[k][i]: c3 += 1
                if adj[i][k] and adj[k][j] and adj[j][i]: c3 += 1

    # Number of tiles set to 1 (flipped)
    n_flipped = 0; n_unflipped = 0
    for a in range(1, n+1):
        for b in range(1, a):
            if a - b >= 2:  # tile
                if adj[a][b]:
                    n_unflipped += 1
                else:
                    n_flipped += 1

    return {
        'H': H, 'scores': score_seq, 'c3': c3,
        'n_flipped': n_flipped, 'n_unflipped': n_unflipped,
        'score_detail': [sum(adj[i][j] for j in range(1, n+1) if j != i) for i in range(1, n+1)],
    }

# ============================================================
# MAIN
# ============================================================

print("=" * 72)
print("  THE PETERSEN TOURNAMENT")
print("  opus-2026-03-24-S288")
print("=" * 72)

n = 10
pet = petersen_edges()
print(f"\n  Petersen graph: {len(pet)} edges on {n} vertices")
print(f"  Edges: {sorted(pet)}")

# Find Hamiltonian paths
print(f"\n  Finding Hamiltonian paths in Petersen graph...")
ham_paths = find_ham_paths(pet, n)
print(f"  Found {len(ham_paths)} directed Hamiltonian paths")
print(f"  = {len(ham_paths)//2} undirected paths (each traversable both ways)")

# Analyze a few Hamiltonian paths
unique_tournaments = {}

for idx, path in enumerate(ham_paths):
    adj, label = build_petersen_tournament(path, pet, n)
    props = tournament_properties(adj, n)

    key = (props['scores'], props['c3'], props['H'])
    if key not in unique_tournaments:
        unique_tournaments[key] = {
            'path': path, 'props': props, 'count': 0, 'label': label
        }
    unique_tournaments[key]['count'] += 1

print(f"\n  {len(unique_tournaments)} DISTINCT tournament types from {len(ham_paths)} paths")

print(f"\n  PETERSEN TOURNAMENTS:")
print(f"  {'#':>3} {'H':>6} {'c3':>4} {'scores':>30} {'#paths':>7} {'flipped':>8} {'unflipped':>10}")

for i, (key, data) in enumerate(sorted(unique_tournaments.items(), key=lambda x: x[0][2])):
    props = data['props']
    print(f"  {i+1:>3} {props['H']:>6} {props['c3']:>4} {str(props['scores']):>30} "
          f"{data['count']:>7} {props['n_flipped']:>8} {props['n_unflipped']:>10}")

# Detailed analysis of each type
for i, (key, data) in enumerate(sorted(unique_tournaments.items(), key=lambda x: x[0][2])):
    props = data['props']
    path = data['path']
    label = data['label']

    print(f"\n  TYPE {i+1}: H={props['H']}, scores={props['scores']}")
    print(f"    Hamiltonian path: {' → '.join(str(v) for v in path)}")
    print(f"    Labels: {dict(sorted(label.items()))}")
    print(f"    Score per vertex (label order): {props['score_detail']}")
    print(f"    Flipped tiles: {props['n_flipped']}/{props['n_flipped']+props['n_unflipped']}")

    # Is it SC?
    # Check if score sequence is palindromic
    ss = props['scores']
    is_palindrome = (ss == tuple(n-1-s for s in reversed(ss)))
    print(f"    Palindromic scores (SC candidate): {is_palindrome}")

    # Szele comparison
    szele = factorial(n) / (2 ** (n-1))
    print(f"    H / Szele expected ({szele:.0f}): {props['H']/szele:.4f}")

print(f"\n{'='*72}")
print("  SYNTHESIS")
print(f"{'='*72}")
print(f"""
  The Petersen graph, viewed as a tournament on 10 vertices:
  - 15 edges → oriented high→low (unflipped = "Petersen arcs")
  - 30 non-edges → oriented low→high (flipped)
  - Base path chosen from a Hamiltonian path in the Petersen graph
  - 6 remaining Petersen edges become unflipped tiles
  - 30 non-Petersen pairs become flipped tiles

  The Petersen graph IS Hamiltonian-path-rich (many paths exist)
  but NOT Hamiltonian-cycle (no cycle — a famous result).
  Each Hamiltonian path choice gives a potentially different tournament.

  The number of tiles: C(9,2) = 36
  Unflipped tiles: 6 (= 15 - 9 = Petersen edges minus path edges)
  Flipped tiles: 30
  Flipped fraction: 30/36 = 83.3%
""")

print("DONE.")
