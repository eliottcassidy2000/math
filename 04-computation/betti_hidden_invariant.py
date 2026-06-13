#!/usr/bin/env python3
"""
betti_hidden_invariant.py - What distinguishes β₁=0 from β₁=1 within same isomorphism class?

KEY FINDING: At n=5, (t₃=3, score=(1,1,2,3,3)) has 120 with β₁=0 and 120 with β₁=1.
These tournaments have IDENTICAL:
  - F-polynomial [9,30,42,30,9]
  - 3-cycle count (t₃=3)
  - Score sequence (1,1,2,3,3)
  - Strong connectivity (all SC=True)

But DIFFERENT path homology!

QUESTION: Are these isomorphic (same tournament up to relabeling)?
If YES: path homology is NOT a tournament isomorphism invariant (unlikely).
If NO: what distinguishes them? It must be the CYCLE GRAPH structure.

Let me compute:
1. The odd-cycle graph Ω(T) for each
2. The independence number α(Ω(T))
3. Whether the 3-cycles overlap in specific patterns

Author: opus-2026-03-07-S46e
"""
import sys
sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from path_homology_v2 import path_betti_numbers
from itertools import combinations, permutations
from collections import Counter, defaultdict

def count_t3(A, n):
    return sum(1 for i,j,k in combinations(range(n), 3)
               if (A[i][j] and A[j][k] and A[k][i]) or
                  (A[i][k] and A[k][j] and A[j][i]))

def score_seq(A, n):
    return tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))

def f_polynomial(A, n):
    counts = [0] * n
    for perm in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[perm[i]][perm[i+1]])
        counts[fwd] += 1
    return tuple(counts)

def find_3cycles(A, n):
    """Return list of directed 3-cycle vertex sets."""
    cycles = []
    for combo in combinations(range(n), 3):
        i, j, k = combo
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            cycles.append(frozenset(combo))
    return cycles

def cycle_graph_structure(cycles):
    """Compute the structure of the 3-cycle overlap graph.
    Vertices = 3-cycles, edges = pairs sharing a vertex."""
    n_cycles = len(cycles)
    adj = [[0]*n_cycles for _ in range(n_cycles)]
    for i in range(n_cycles):
        for j in range(i+1, n_cycles):
            if len(cycles[i] & cycles[j]) > 0:
                adj[i][j] = 1
                adj[j][i] = 1

    # Compute: number of edges, max independent set size
    n_edges = sum(adj[i][j] for i in range(n_cycles) for j in range(i+1, n_cycles))

    # Max independent set (brute force for small sizes)
    max_indep = 0
    for size in range(n_cycles, 0, -1):
        found = False
        for subset in combinations(range(n_cycles), size):
            if all(adj[subset[a]][subset[b]] == 0 for a in range(len(subset)) for b in range(a+1, len(subset))):
                max_indep = size
                found = True
                break
        if found:
            break

    # Overlap pattern: how many vertices does each pair of cycles share?
    overlap_sizes = Counter()
    for i in range(n_cycles):
        for j in range(i+1, n_cycles):
            overlap_sizes[len(cycles[i] & cycles[j])] += 1

    return {
        'n_cycles': n_cycles,
        'n_edges': n_edges,
        'max_indep': max_indep,
        'overlap_sizes': dict(overlap_sizes),
    }

def isomorphism_type(A, n):
    """Compute a canonical form for the tournament (up to vertex relabeling)."""
    # Use the sorted row-sum profile as a coarse invariant
    # For finer: compute the sorted adjacency matrix under all permutations
    # This is expensive but n=5 has only 120 permutations
    best = None
    for perm in permutations(range(n)):
        mat = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if best is None or mat < best:
            best = mat
    return best

# === n=5: Focus on (t₃=3, score=(1,1,2,3,3)) class ===
print("=" * 60)
print("n=5: t₃=3, score=(1,1,2,3,3) — hidden invariant")
print("=" * 60)

n = 5
m = n*(n-1)//2

b1_0_tournaments = []
b1_1_tournaments = []

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
            idx += 1

    t3 = count_t3(A, n)
    score = score_seq(A, n)
    if t3 != 3 or score != (1, 1, 2, 3, 3):
        continue

    beta = path_betti_numbers(A, n)
    b1 = int(beta[1]) if len(beta) > 1 else 0

    cycles = find_3cycles(A, n)
    cg = cycle_graph_structure(cycles)

    if b1 == 0:
        b1_0_tournaments.append({'A': A, 'bits': bits, 'cycles': cycles, 'cg': cg})
    else:
        b1_1_tournaments.append({'A': A, 'bits': bits, 'cycles': cycles, 'cg': cg})

print(f"\nβ₁=0: {len(b1_0_tournaments)} tournaments")
print(f"β₁=1: {len(b1_1_tournaments)} tournaments")

# Compare cycle graph structures
cg0_types = Counter()
cg1_types = Counter()
for d in b1_0_tournaments:
    cg = d['cg']
    key = (cg['n_cycles'], cg['n_edges'], cg['max_indep'], tuple(sorted(cg['overlap_sizes'].items())))
    cg0_types[key] += 1
for d in b1_1_tournaments:
    cg = d['cg']
    key = (cg['n_cycles'], cg['n_edges'], cg['max_indep'], tuple(sorted(cg['overlap_sizes'].items())))
    cg1_types[key] += 1

print(f"\nβ₁=0 cycle graph types: {dict(cg0_types)}")
print(f"β₁=1 cycle graph types: {dict(cg1_types)}")

# Show a specific example from each
print("\nExample β₁=0:")
d = b1_0_tournaments[0]
print(f"  bits={d['bits']:010b}")
for i in range(n):
    out = [j for j in range(n) if d['A'][i][j]]
    print(f"  {i} → {out}")
print(f"  3-cycles: {[set(c) for c in d['cycles']]}")

print("\nExample β₁=1:")
d = b1_1_tournaments[0]
print(f"  bits={d['bits']:010b}")
for i in range(n):
    out = [j for j in range(n) if d['A'][i][j]]
    print(f"  {i} → {out}")
print(f"  3-cycles: {[set(c) for c in d['cycles']]}")

# Check isomorphism types
print("\n" + "=" * 60)
print("ISOMORPHISM TYPES")
print("=" * 60)

# This is expensive but doable for 240 tournaments at n=5
iso_by_beta = {0: Counter(), 1: Counter()}
for b1_val, tournaments in [(0, b1_0_tournaments), (1, b1_1_tournaments)]:
    for d in tournaments:
        iso = isomorphism_type(d['A'], n)
        iso_by_beta[b1_val][iso] += 1

n_iso_0 = len(iso_by_beta[0])
n_iso_1 = len(iso_by_beta[1])
print(f"\nβ₁=0: {n_iso_0} isomorphism types (120 tournaments)")
print(f"β₁=1: {n_iso_1} isomorphism types (120 tournaments)")

# Do any isomorphism types appear in BOTH classes?
overlap = set(iso_by_beta[0].keys()) & set(iso_by_beta[1].keys())
print(f"\nOverlapping isomorphism types: {len(overlap)}")
if overlap:
    print("  → Path homology is NOT an isomorphism invariant! Different labelings give different β₁!")
else:
    print("  → Path homology IS an isomorphism invariant (as expected)")
    print(f"  → There are {n_iso_0} iso types with β₁=0 and {n_iso_1} with β₁=1")
    print(f"  → These are genuinely different tournament structures")
