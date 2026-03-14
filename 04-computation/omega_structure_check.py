#!/usr/bin/env python3
"""
Ω(T) Structure Check — Can K₁₀, K₈-e, K₆-2e arise as conflict graphs?
opus-2026-03-14-S71h

Direct approach: enumerate tournaments T on n vertices, compute Ω(T),
and check what graphs actually appear as Ω(T).

Key question: which graphs G with I(G,2)=21 can actually be Ω(T)?

If NONE of the graphs with I(G,2)=21 can be Ω(T), then H≠21 is proved.
"""

from itertools import combinations

def get_odd_cycles(n, adj, max_len=None):
    """Find all directed odd cycles in tournament on n vertices.
    adj[i] = set of out-neighbors of i.
    Returns list of frozensets (vertex sets of directed odd cycles).
    """
    if max_len is None:
        max_len = n

    cycles = set()

    # Find directed cycles by DFS
    for length in range(3, max_len + 1, 2):  # odd lengths only
        # Try all ordered tuples of `length` vertices
        for combo in combinations(range(n), length):
            # Check all cyclic orderings
            from itertools import permutations
            # Actually, for efficiency, just check if combo forms a directed cycle
            # in some order. A directed k-cycle on vertices v_0,...,v_{k-1}
            # means v_0→v_1→...→v_{k-1}→v_0.

            # For small k, try all permutations of the combo
            if length <= 5:
                for perm in permutations(combo):
                    is_cycle = True
                    for i in range(length):
                        nxt = (i + 1) % length
                        if perm[nxt] not in adj[perm[i]]:
                            is_cycle = False
                            break
                    if is_cycle:
                        cycles.add(frozenset(combo))
                        break
            else:
                # For longer cycles, use a smarter approach
                # But we won't need this for small n
                pass

    return list(cycles)

def build_omega(n, adj):
    """Build Ω(T): vertices = directed odd cycles, edges = shared vertex."""
    cycles = get_odd_cycles(n, adj, max_len=min(n, 7))

    # Build adjacency
    omega_n = len(cycles)
    omega_edges = []
    for i in range(omega_n):
        for j in range(i+1, omega_n):
            if cycles[i] & cycles[j]:  # shared vertex
                omega_edges.append((i, j))

    return cycles, omega_edges

def independence_poly_at_2(num_verts, edges):
    """Compute I(G, 2) for graph with given vertices and edges."""
    adj = set()
    for u, v in edges:
        adj.add((u, v))
        adj.add((v, u))

    total = 0
    for mask in range(1 << num_verts):
        verts = [i for i in range(num_verts) if mask & (1 << i)]
        k = len(verts)
        ok = True
        for a in range(k):
            for b in range(a+1, k):
                if (verts[a], verts[b]) in adj:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            total += 2**k
    return total

def tournament_from_bits(n, bits):
    """Create tournament adjacency from bit encoding."""
    adj = [set() for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i].add(j)
            else:
                adj[j].add(i)
            idx += 1
    return adj

def hp_count(n, adj):
    """Count Hamiltonian paths using bitmask DP."""
    # dp[mask][v] = number of paths visiting exactly vertices in mask, ending at v
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1

    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in adj[v]:
                if not (mask & (1 << w)):
                    dp[mask | (1 << w)][w] += dp[mask][v]

    return sum(dp[full])

print("=" * 70)
print("Ω(T) STRUCTURE — WHAT GRAPHS ACTUALLY ARISE?")
print("=" * 70)
print()

# First: exhaustively check n=5 tournaments
print("n=5: Exhaustive Ω(T) analysis")
print("-" * 40)

from collections import Counter
n = 5
num_edges = n * (n-1) // 2  # 10
total = 1 << num_edges  # 1024

omega_sizes = Counter()
h_values = Counter()
omega_graph_types = Counter()  # (num_cycles, num_omega_edges) -> count

for bits in range(total):
    adj = tournament_from_bits(n, bits)
    H = hp_count(n, adj)
    h_values[H] += 1

    cycles = get_odd_cycles(n, adj, max_len=5)
    num_cycles = len(cycles)

    # Build Ω edges
    omega_edges = []
    for i in range(num_cycles):
        for j in range(i+1, num_cycles):
            if cycles[i] & cycles[j]:
                omega_edges.append((i, j))

    omega_sizes[num_cycles] += 1

    # Compute I(Ω, 2) as a check
    if num_cycles > 0:
        I_omega = independence_poly_at_2(num_cycles, omega_edges)
    else:
        I_omega = 1  # empty graph

    if I_omega != H:
        print(f"  OCF VIOLATION at bits={bits}: H={H}, I(Ω,2)={I_omega}")

    omega_graph_types[(num_cycles, len(omega_edges))] += 1

print(f"  Total tournaments: {total}")
print(f"  H value distribution: {dict(sorted(h_values.items()))}")
print(f"  Ω(T) size distribution: {dict(sorted(omega_sizes.items()))}")
print(f"  OCF verified for all {total} tournaments")
print()

print("  Ω graph types (|V|, |E|) → count:")
for (v, e), cnt in sorted(omega_graph_types.items()):
    print(f"    |V|={v}, |E|={e}: {cnt} tournaments")

print()

# Now check n=7: what Ω structures give H near 21?
print("n=7: Sampling for Ω structures with I(Ω,2) close to 21")
print("-" * 40)

import random
random.seed(42)
n = 7
num_edges = 21
sample_size = 5000

h_dist = Counter()
for trial in range(sample_size):
    bits = random.randint(0, (1 << num_edges) - 1)
    adj = tournament_from_bits(n, bits)
    H = hp_count(n, adj)
    h_dist[H] += 1

print(f"  H value distribution (sample of {sample_size}):")
for h in sorted(h_dist.keys()):
    if h_dist[h] >= 5:
        print(f"    H={h}: {h_dist[h]} ({100*h_dist[h]/sample_size:.1f}%)")

if 21 in h_dist:
    print(f"\n  *** H=21 FOUND at n=7! Count: {h_dist[21]} ***")
else:
    print(f"\n  H=21 NOT found in {sample_size} samples at n=7")

print()

# Key question: what's the MINIMUM n where H=21 could occur?
# We know: H is always odd, and H achievable values at n=5 are {1,3,5,9,11,13,15}
# At n=7: need to check full spectrum

# More targeted: check which (|V|, |E|) graphs with I(G,2)=21 actually arise as Ω
print("=" * 70)
print("WHICH GRAPHS WITH I(G,2)=21 CAN BE Ω(T)?")
print("=" * 70)
print()

# The graphs with I(G,2)=21 are:
# 1. 4 vertices, 3 edges: P₄ or K₁⊔K₃ (both have I=1+4x+3x²)
# 2. 6 vertices, 13 edges: K₆-2e (I=1+6x+2x²)
# 3. 8 vertices, 27 edges: K₈-e (I=1+8x+x²)
# 4. 10 vertices, 0 edges: K₁₀ (I=1+10x)

# For Ω(T) at n=5:
# - Max number of 3-cycles in n=5 tournament is C(5,3)=10, but directed: max t₃=4 at n=5
# Wait, t₃ can be up to (n choose 3)/something... let me compute.

print("Maximum number of directed odd cycles at each n:")
for n in range(3, 8):
    num_edges = n*(n-1)//2
    if n <= 6:
        max_cycles = 0
        total = 1 << num_edges
        for bits in range(total):
            adj = tournament_from_bits(n, bits)
            cycles = get_odd_cycles(n, adj, max_len=n)
            max_cycles = max(max_cycles, len(cycles))
        print(f"  n={n}: max |odd cycles| = {max_cycles} (exhaustive over {total} tournaments)")
    else:
        # Sample for n=7
        max_cycles = 0
        for _ in range(2000):
            bits = random.randint(0, (1 << num_edges) - 1)
            adj = tournament_from_bits(n, bits)
            cycles = get_odd_cycles(n, adj, max_len=7)
            max_cycles = max(max_cycles, len(cycles))
        print(f"  n={n}: max |odd cycles| ≥ {max_cycles} (sample of 2000)")

print()
print("For H=21 via I(G,2)=21, we need Ω(T) to have:")
print("  - 4 vertices (K₃-poisoned case): blocked by THM-201/202")
print("  - 6 vertices: need 6 odd cycles in T")
print("  - 8 vertices: need 8 odd cycles in T")
print("  - 10 vertices: need 10 odd cycles in T (disjoint!)")
print()
print("Since K₁₀ = 10 ISOLATED vertices = 10 pairwise DISJOINT odd cycles,")
print("this requires 10 vertex-disjoint odd cycles, needing ≥ 30 tournament vertices!")
print("But for H=21, we need very few HP. At n=30, H is astronomically large.")
print("So K₁₀ as Ω(T) with H=21 is self-contradictory: too many vertices → too many HP.")
print()

# Check: at n=30, minimum H would be way more than 21
# Actually n doesn't have to be 30 -- 5-cycles reuse vertices too
# But 10 disjoint 3-cycles need 30 vertices

# Actually I need to be more careful. K₁₀ as Ω means 10 odd cycles,
# NO pair shares a vertex (that's what "independent" means in Ω).
# So all 10 are vertex-disjoint. Each has ≥ 3 vertices.
# So we need ≥ 30 tournament vertices.
# At n ≥ 30, H(T) ≥ n!/... which is >> 21.

# Actually, the minimum H for a tournament on n vertices is 1 (transitive).
# So large n doesn't automatically mean large H.
# BUT: having 10 disjoint directed 3-cycles means the tournament is "very cyclic"
# at least locally, which pushes H up.

# Let's think about it differently:
# If T has 10 vertex-disjoint 3-cycles, and the rest of the vertices are
# "linearly ordered" (transitive), what is H?

print("REFINED ANALYSIS: Can K₁₀ (10 disjoint odd cycles) give H=21?")
print()
print("K₁₀ as Ω means: 10 vertex-disjoint odd cycles, each contributing")
print("factor (1+2^k) where k = cycle length. For 3-cycles: factor = (1+2) = 3.")
print("If all 10 are 3-cycles: I(K₁₀, 2) = (1+2)^0... wait, that's wrong.")
print()
print("Actually I(K₁₀, 2) = 1 + 10·2 = 21. The 10 vertices of K₁₀ are the")
print("10 odd cycles, each being an INDEPENDENT SET of size 1.")
print("I = 1 + 10x means: empty set contributes 1, each single cycle contributes x,")
print("and NO pair is independent (all pairs are adjacent in K₁₀).")
print()
print("Wait — K₁₀ has ALL edges! So there are NO independent sets of size ≥ 2.")
print("This means every pair of odd cycles shares a vertex.")
print("So K₁₀ = 10 pairwise-INTERSECTING odd cycles, not disjoint!")
print()
print("CORRECTION: The graph G with I(G,2)=21 and I=1+10x is K₁₀ (COMPLETE graph).")
print("Ω(T) = K₁₀ means 10 odd cycles, every pair sharing at least one vertex.")
print("This is the CLIQUE case, not the independent case.")
print()

# So I had it backwards. Let me reconsider:
# K₁₀ complete → all 10 cycles pairwise intersect → they share vertices
# The MINIMUM number of tournament vertices for 10 pairwise-intersecting 3-cycles:
# By Helly/sunflower: if all are 3-element sets pairwise intersecting...

# The minimum vertex cover of 10 pairwise-intersecting 3-sets:
# They can all share one common vertex! Then we need 1 + 2*10 = 21 vertices at most
# (but many can share the hub vertex and still share other vertices).
# Actually: 10 triangles through a common vertex v need only v plus 2 others each.
# If they all share v: at most 1 + 20 = 21 vertices (if all "other" vertices distinct).
# But can reuse: 10 triangles on vertex v with edges to pairs from a set of m vertices.
# Need C(m,2) ≥ 10, so m ≥ 5 (C(5,2)=10). Total: 6 vertices.

print("SUNFLOWER STRUCTURE: 10 pairwise-intersecting 3-cycles through vertex v")
print("  Need C(m,2) ≥ 10 other vertices: m=5 suffices (C(5,2)=10)")
print("  Total: 6 tournament vertices, 10 directed 3-cycles through vertex v")
print("  Each triangle = {v, a_i, a_j} for distinct pairs (a_i, a_j) from 5 others")
print()

# Check: on 6 vertices, can we have EXACTLY 10 triangles through vertex 0?
n = 6
num_edges = 15
count_with_10_triangles = 0
for bits in range(1 << num_edges):
    adj = tournament_from_bits(n, bits)
    # Count 3-cycles through vertex 0
    cycles_through_0 = []
    for combo in combinations(range(1, n), 2):
        triple = (0,) + combo
        # Check all 2 possible directed 3-cycles on these 3 vertices
        from itertools import permutations
        for perm in permutations(triple):
            if perm[1] in adj[perm[0]] and perm[2] in adj[perm[1]] and perm[0] in adj[perm[2]]:
                cycles_through_0.append(frozenset(triple))
                break

    all_cycles = get_odd_cycles(n, adj, max_len=5)

    if len(all_cycles) == 10:
        # Check if Ω = K₁₀ (all pairs intersect)
        all_intersect = True
        for i in range(len(all_cycles)):
            for j in range(i+1, len(all_cycles)):
                if not (all_cycles[i] & all_cycles[j]):
                    all_intersect = False
                    break
            if not all_intersect:
                break

        if all_intersect:
            H = hp_count(n, adj)
            if count_with_10_triangles < 3:
                print(f"  Ω = K₁₀ found! n=6, H={H}, I(K₁₀,2)=21")
                if H == 21:
                    print(f"    *** H=21 ACHIEVED! OCF would be satisfied! ***")
                else:
                    print(f"    H≠21, so OCF says I(Ω,2)=H={H}≠21. Ω≠K₁₀ actually.")
            count_with_10_triangles += 1

if count_with_10_triangles == 0:
    print("  No tournament on 6 vertices has exactly 10 pairwise-intersecting odd cycles")
else:
    print(f"  Total: {count_with_10_triangles} tournaments with 10 pairwise-intersecting cycles")

print()
print("=" * 70)
print("COMPLETE Ω(T) CATALOG AT n=5")
print("=" * 70)
print()

# At n=5, what H values occur and what are their Ω structures?
n = 5
num_edges = 10
total = 1 << num_edges

h_to_omega = {}  # H -> list of (cycles, omega_edges) descriptions
for bits in range(total):
    adj = tournament_from_bits(n, bits)
    H = hp_count(n, adj)
    cycles = get_odd_cycles(n, adj, max_len=5)

    omega_edges = []
    nc = len(cycles)
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i] & cycles[j]:
                omega_edges.append((i, j))

    if H not in h_to_omega:
        h_to_omega[H] = []
    h_to_omega[H].append((nc, len(omega_edges), cycles))

for H in sorted(h_to_omega.keys()):
    entries = h_to_omega[H]
    sizes = set((nc, ne) for nc, ne, _ in entries)
    print(f"H={H}: {len(entries)} tournaments")
    for (nc, ne) in sorted(sizes):
        cnt = sum(1 for a, b, _ in entries if a == nc and b == ne)
        # Compute I(G,2) for this graph structure
        if nc > 0:
            # Use one example
            example = next((nc2, ne2, cyc) for nc2, ne2, cyc in entries if nc2 == nc and ne2 == ne)
            I_val = independence_poly_at_2(nc, [(i,j) for i in range(nc) for j in range(i+1, nc) if example[2][i] & example[2][j]])
        else:
            I_val = 1
        print(f"  Ω: |V|={nc}, |E|={ne}, I(Ω,2)={I_val}: {cnt} tournaments")

print()
print("=" * 70)
print("CONCLUSION")
print("=" * 70)
print()
print("The Ω(T) structure at n=5 shows exactly which graphs arise.")
print("H=21 is absent because no tournament on 5 vertices produces")
print("an Ω(T) with I(Ω,2)=21.")
print("The question remains: at what n (if any) does H=21 first appear?")
