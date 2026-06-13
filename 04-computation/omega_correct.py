#!/usr/bin/env python3
"""
CORRECTED Ω(T) computation — count ALL distinct directed odd cycles.
opus-2026-03-14-S71h

Bug fix: the previous script collapsed multiple directed cycles on the same
vertex set into one. A directed k-cycle is specified by a cyclic ordering
(up to rotation). On n=5 vertices, there can be multiple directed 5-cycles.
"""

from itertools import combinations, permutations
from collections import Counter

def find_all_directed_odd_cycles(n, adj, max_len=None):
    """Find ALL distinct directed odd cycles in a tournament.

    Each cycle is represented by its canonical form:
    the lexicographically smallest rotation of its vertex sequence.

    Returns list of (canonical_tuple, vertex_frozenset) pairs.
    """
    if max_len is None:
        max_len = n

    cycles = []
    seen = set()

    for length in range(3, max_len + 1, 2):  # odd lengths only
        for combo in combinations(range(n), length):
            # Try all permutations of this vertex set as potential cycles
            for perm in permutations(combo):
                # Check if perm forms a directed cycle: p[0]→p[1]→...→p[k-1]→p[0]
                is_cycle = True
                for i in range(length):
                    nxt = (i + 1) % length
                    if perm[nxt] not in adj[perm[i]]:
                        is_cycle = False
                        break

                if is_cycle:
                    # Normalize: find canonical rotation (smallest starting vertex)
                    rotations = [perm[i:] + perm[:i] for i in range(length)]
                    canonical = min(rotations)

                    if canonical not in seen:
                        seen.add(canonical)
                        cycles.append((canonical, frozenset(combo)))

    return cycles

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

print("=" * 70)
print("CORRECTED Ω(T) COMPUTATION")
print("=" * 70)
print()

# First: check how many directed 5-cycles a tournament on 5 vertices can have
print("Directed 5-cycle counts at n=5:")
print("-" * 40)
n = 5
num_edges = 10
total = 1 << num_edges

cycle5_counts = Counter()
for bits in range(total):
    adj = tournament_from_bits(n, bits)
    cycles = find_all_directed_odd_cycles(n, adj, max_len=5)
    # Count 5-cycles
    c5 = sum(1 for c, vs in cycles if len(c) == 5)
    cycle5_counts[c5] += 1

print("  Distribution of #(directed 5-cycles):")
for k in sorted(cycle5_counts.keys()):
    print(f"    {k} five-cycles: {cycle5_counts[k]} tournaments")

print()

# Now verify OCF with corrected Ω
print("=" * 70)
print("OCF VERIFICATION WITH CORRECTED Ω AT n=5")
print("=" * 70)
print()

violations = 0
h_to_omega_info = {}

for bits in range(total):
    adj = tournament_from_bits(n, bits)
    H = hp_count(n, adj)

    cycles = find_all_directed_odd_cycles(n, adj, max_len=5)
    nc = len(cycles)

    # Build Ω edges: two cycles are adjacent iff they share a vertex
    omega_edges = []
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i][1] & cycles[j][1]:  # frozensets intersect
                omega_edges.append((i, j))

    I_omega = independence_poly_at_2(nc, omega_edges) if nc > 0 else 1

    if I_omega != H:
        violations += 1
        if violations <= 5:
            print(f"  VIOLATION bits={bits}: H={H}, I(Ω,2)={I_omega}, |Ω|={nc}")
            # Show cycle details
            for k, (cyc, vs) in enumerate(cycles):
                print(f"    cycle {k}: {cyc} (vertices {sorted(vs)})")

    # Record for catalog
    key = (nc, len(omega_edges))
    if H not in h_to_omega_info:
        h_to_omega_info[H] = Counter()
    h_to_omega_info[H][key] += 1

if violations == 0:
    print("  *** OCF VERIFIED FOR ALL 1024 TOURNAMENTS AT n=5 ***")
else:
    print(f"\n  VIOLATIONS: {violations}/{total}")

print()
print("Ω(T) catalog at n=5:")
for H in sorted(h_to_omega_info.keys()):
    print(f"  H={H}:")
    for (nc, ne), cnt in sorted(h_to_omega_info[H].items()):
        I_val = f"I={independence_poly_at_2(nc, list(combinations(range(nc), 2))[:ne])}" if nc > 0 else "I=1"
        print(f"    |V|={nc}, |E|={ne}: {cnt} tournaments")

print()

# Now check n=4 for sanity
print("=" * 70)
print("OCF VERIFICATION AT n=4")
print("=" * 70)
print()

n = 4
num_edges = 6
total = 1 << num_edges
violations = 0

for bits in range(total):
    adj = tournament_from_bits(n, bits)
    H = hp_count(n, adj)

    cycles = find_all_directed_odd_cycles(n, adj, max_len=3)  # only 3-cycles at n=4
    nc = len(cycles)

    omega_edges = []
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i][1] & cycles[j][1]:
                omega_edges.append((i, j))

    I_omega = independence_poly_at_2(nc, omega_edges) if nc > 0 else 1

    if I_omega != H:
        violations += 1
        if violations <= 3:
            print(f"  VIOLATION bits={bits}: H={H}, I(Ω,2)={I_omega}, |Ω|={nc}")
            for k, (cyc, vs) in enumerate(cycles):
                print(f"    cycle {k}: {cyc} (vertices {sorted(vs)})")

if violations == 0:
    print("  *** OCF VERIFIED FOR ALL 64 TOURNAMENTS AT n=4 ***")
else:
    print(f"  VIOLATIONS: {violations}/{total}")

print()

# Now n=3
print("OCF at n=3:")
n = 3
for bits in range(8):
    adj = tournament_from_bits(n, bits)
    H = hp_count(n, adj)
    cycles = find_all_directed_odd_cycles(n, adj, max_len=3)
    nc = len(cycles)
    omega_edges = []
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i][1] & cycles[j][1]:
                omega_edges.append((i, j))
    I_omega = independence_poly_at_2(nc, omega_edges) if nc > 0 else 1
    status = "✓" if I_omega == H else "✗"
    print(f"  bits={bits}: H={H}, |Ω|={nc}, I(Ω,2)={I_omega} {status}")
