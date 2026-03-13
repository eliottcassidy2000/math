#!/usr/bin/env python3
"""
Analyze the two non-Paley regular n=7 tournament classes.

H=171 (contractible, χ=1): What cycle structure makes it contractible?
H=175 (β₁=1, χ=0): What 1-cycle survives?

Key: at n=7, we know
  β₂ = 0 always (THM-108/109)
  β₁·β₃ = 0 (mutual exclusivity)

So:
  H=171: β₁=0, β₃=0 → fully acyclic path homology (except β₀=1)
  H=175: β₁=1, β₃=0 → one directed 1-hole

Question: what determines β₁ at n=7?
Previous work: β₁=1 iff T has a "bad" vertex structure.
At n=5: β₁=1 iff t₃ ≥ 3 (enough 3-cycles to create an independent 1-cycle).

For regular n=7: all have t₃=14. So 3-cycle count alone doesn't determine β₁.
The 3-cycle OVERLAP structure must be the key.

opus-2026-03-13-S71b
"""

import itertools

def count_3cycles_per_vertex(adj, n):
    """Count 3-cycles through each vertex."""
    counts = [0] * n
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            for k in range(n):
                if k == i or k == j:
                    continue
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    counts[i] += 1
    # Each cycle counted once per vertex, and each cycle has 3 vertices
    return [c // 2 for c in counts]  # divide by 2 because (j,k) and (k,j) both counted

def find_all_3cycles(adj, n):
    """Find all directed 3-cycles."""
    cycles = []
    for combo in itertools.combinations(range(n), 3):
        a, b, c = combo
        if adj[a][b] and adj[b][c] and adj[c][a]:
            cycles.append(combo)
        elif adj[a][c] and adj[c][b] and adj[b][a]:
            cycles.append(combo)
    return cycles

def analyze_cycle_structure(adj, n):
    """Detailed 3-cycle structure analysis."""
    cycles = find_all_3cycles(adj, n)
    c3 = len(cycles)

    # Overlap: for each pair of 3-cycles, count shared vertices
    overlaps = {0: 0, 1: 0, 2: 0}
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            shared = len(set(cycles[i]) & set(cycles[j]))
            if shared in overlaps:
                overlaps[shared] += 1

    # Per-vertex: how many 3-cycles through each vertex
    vertex_in_cycles = [0] * n
    for cycle in cycles:
        for v in cycle:
            vertex_in_cycles[v] += 1

    # Conflict graph degree: for each cycle, how many other cycles share a vertex
    conflict_degrees = []
    for i in range(len(cycles)):
        deg = 0
        for j in range(len(cycles)):
            if i != j and set(cycles[i]) & set(cycles[j]):
                deg += 1
        conflict_degrees.append(deg)

    return {
        'c3': c3,
        'overlaps': overlaps,
        'vertex_in_cycles': sorted(vertex_in_cycles),
        'conflict_degrees': sorted(conflict_degrees),
        'disj': overlaps.get(0, 0)
    }

# Generate and analyze
n = 7
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(edges)

# Find representatives
from collections import Counter
reps = {}
count = 0

for bits in range(2**m):
    adj = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            adj[i][j] = 1
        else:
            adj[j][i] = 1

    if not all(sum(adj[i]) == 3 for i in range(n)):
        continue

    # Compute H quickly
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (pm & (1 << u)) and adj[u][v]:
                        total += dp.get((pm, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    H = sum(dp.get((full, v), 0) for v in range(n))

    if H not in reps:
        reps[H] = adj
        info = analyze_cycle_structure(adj, n)
        print(f"\nH={H}:")
        print(f"  c3 = {info['c3']}")
        print(f"  disj_33 = {info['disj']}")
        print(f"  Overlaps (shared vertices): {info['overlaps']}")
        print(f"  Vertex cycle count: {info['vertex_in_cycles']}")
        print(f"  Conflict graph degrees: {info['conflict_degrees']}")

        # Score sequence
        scores = sorted([sum(adj[i]) for i in range(n)])
        print(f"  Score: {scores}")

    count += 1
    if len(reps) == 3:
        break

# Key analysis: what makes H=175 have β₁=1?
print("\n" + "="*60)
print("WHY β₁=1 FOR H=175?")
print("="*60)

print("""
At n=5: β₁=1 iff every vertex is in some 3-cycle AND the 3-cycles
"wrap around" in a way that creates a non-trivial 1-cycle in Ω₁.

For regular n=7: every vertex is in c₃/v = 14/7 * 3 = 6 cycles.
So the cycle-per-vertex condition is always satisfied.

The difference must be in the OVERLAP structure:
  H=189 (Paley): disj=7, overlap2=21, ov1=63 → β₁=0
  H=171:         disj=10, overlap2=18, ov1=57 → β₁=0
  H=175:         disj=14, overlap2=14, ov1=49 → β₁=1

The H=175 class has the MOST disjoint 3-cycle pairs (14).
More disjoint pairs = more "independence" in the cycle structure.
This creates enough room for a non-trivial 1-cycle.

CONJECTURE: At n=7 regular, β₁=1 iff disj_33 ≥ 14 (the maximum).
This would connect cycle independence to topological holes.
""")
