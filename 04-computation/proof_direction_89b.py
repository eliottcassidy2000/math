#!/usr/bin/env python3
"""
PROOF DIRECTION FOR THM-209: Why H = IP(G, 2)?
opus-2026-03-14-S89b

The formula H(T) = Σ_{S ⊆ OddCyc(T), pairwise disjoint} 2^|S| says:
each collection of k pairwise disjoint odd cycles contributes 2^k
Hamiltonian paths.

KEY INSIGHT from the OCF (Odd Cycle Collection Formula):
  H(T) ≡ 1 (mod 2) because H = Σ_Ω (-1)^{n-comp(Ω)} where Ω ranges
  over odd-cycle covers. The 2^k comes from inclusion-exclusion on the
  cycle orientations.

APPROACH: Consider the "cycle reversal" operation. For each odd cycle C,
reversing all arcs of C gives a different tournament T' with the SAME
set of Hamiltonian paths EXCEPT those that traverse C.

More precisely: if C is a directed 3-cycle i→j→k→i, then reversing C
gives i←j←k←i (= i→k→j→i). Any Hamiltonian path that uses two
consecutive arcs from C will be "flipped" — some paths are created,
others destroyed.

The coefficient 2 for a single cycle C means: the DIFFERENCE in H
between T and T' (where T' reverses C) is ±2. This is because:
  - Exactly 2 Hamiltonian paths are affected by the cycle reversal
  - These 2 paths correspond to the 2 ways of traversing the 3-cycle
    within a Hamiltonian path

Let's verify this "cycle reversal" interpretation.
"""

from itertools import combinations, permutations
from collections import defaultdict

def compute_H(n, adj):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if (mask, v) not in dp: continue
            val = dp[(mask, v)]
            for u in range(n):
                if mask & (1 << u): continue
                if adj[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + val
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def get_ham_paths(n, adj):
    """Return the set of all Hamiltonian paths (as tuples)."""
    paths = []
    for perm in permutations(range(n)):
        ok = True
        for i in range(n-1):
            if not adj[perm[i]][perm[i+1]]:
                ok = False
                break
        if ok:
            paths.append(perm)
    return paths

def tournament_adj(n, bits):
    adj = [[False]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
            idx += 1
    return adj

def reverse_cycle(adj, cycle_verts, cycle_direction):
    """Reverse a directed cycle in the tournament.
    cycle_direction is a list of arcs [(a,b), (b,c), (c,a), ...]"""
    new_adj = [row[:] for row in adj]
    for (a, b) in cycle_direction:
        new_adj[a][b] = False
        new_adj[b][a] = True
    return new_adj

def get_directed_3cycles(n, adj):
    """Return list of (vertex_set, directed_arcs) for each 3-cycle."""
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    cycles.append((frozenset([i,j,k]), [(i,j),(j,k),(k,i)]))
                elif adj[i][k] and adj[k][j] and adj[j][i]:
                    cycles.append((frozenset([i,j,k]), [(i,k),(k,j),(j,i)]))
    return cycles

print("="*70)
print("PROOF DIRECTION: Why does each disjoint cycle set contribute 2^k?")
print("="*70)

# ===== PART 1: Cycle reversal and H difference =====
print("\n" + "="*70)
print("PART 1: CYCLE REVERSAL — ΔH for reversing a single 3-cycle")
print("="*70)

n = 5
m = n*(n-1)//2

delta_h_dist = defaultdict(int)
for bits in range(1 << m):
    adj = tournament_adj(n, bits)
    H_orig = compute_H(n, adj)
    cycles = get_directed_3cycles(n, adj)

    for vs, arcs in cycles:
        adj_rev = reverse_cycle(adj, vs, arcs)
        H_rev = compute_H(n, adj_rev)
        delta = H_rev - H_orig
        delta_h_dist[delta] += 1

print(f"  n={n}: ΔH distribution when reversing a 3-cycle:")
for d in sorted(delta_h_dist.keys()):
    print(f"    ΔH = {d}: {delta_h_dist[d]} cases")

print(f"\n  Note: ΔH when reversing a 3-cycle is NOT always ±2.")
print(f"  The reversal can change MULTIPLE cycles simultaneously.")

# ===== PART 2: Path decomposition — which paths use which cycles? =====
print("\n" + "="*70)
print("PART 2: WHICH PATHS USE WHICH CYCLE ARCS?")
print("="*70)

# For a specific tournament, classify each Hamiltonian path by which cycles it traverses
bits = 0b0101010101  # arbitrary
adj = tournament_adj(n, bits)
H = compute_H(n, adj)
paths = get_ham_paths(n, adj)
cycles = get_directed_3cycles(n, adj)

print(f"  Tournament (bits={bits}): H={H}, {len(cycles)} 3-cycles")
for i, (vs, arcs) in enumerate(cycles):
    arc_set = set(arcs)
    print(f"    C{i}: {sorted(vs)}, arcs={arcs}")

    # Which paths use at least 2 consecutive arcs from this cycle?
    using_cycle = 0
    for path in paths:
        path_arcs = set((path[j], path[j+1]) for j in range(n-1))
        shared = path_arcs & arc_set
        if len(shared) >= 2:
            using_cycle += 1
    print(f"      Paths using ≥2 arcs from C{i}: {using_cycle}")

# ===== PART 3: Inclusion-exclusion interpretation =====
print("\n" + "="*70)
print("PART 3: INCLUSION-EXCLUSION / MÖBIUS FUNCTION")
print("="*70)

print("""
The OCF says: H(T) = Σ_Ω (-1)^{n-comp(Ω)} where Ω ranges over
"odd-cycle covers" of T. The 2^k pattern suggests a BINARY choice
at each cycle: include it or not.

For the Möbius function on the Boolean lattice of cycle subsets:
  μ(∅, S) = (-1)^|S|

And the inclusion-exclusion gives:
  f(S) = Σ_{U ⊆ S} (-1)^{|S|-|U|} g(U)

The formula H = Σ_S 2^|S| where S ranges over disjoint cycle sets
is equivalent to:
  H = Σ_S 2^|S| = evaluation of IP at x=2

This is the PARTITION FUNCTION of the hard-core lattice gas at λ=2.
""")

# ===== PART 4: Direct path analysis at n=4 =====
print("\n" + "="*70)
print("PART 4: DIRECT PATH ANALYSIS AT n=4")
print("="*70)

n = 4
m = n*(n-1)//2

for bits in range(1 << m):
    adj = tournament_adj(n, bits)
    H = compute_H(n, adj)
    cycles = get_directed_3cycles(n, adj)
    t3 = len(cycles)
    paths = get_ham_paths(n, adj)

    if t3 == 1 and H == 3:
        print(f"\n  Tournament with t₃=1, H=3 (bits={bits}):")
        vs, arcs = cycles[0]
        arc_set = set(arcs)
        print(f"    3-cycle: {arcs}")
        print(f"    Hamiltonian paths:")
        for path in paths:
            path_arcs = [(path[j], path[j+1]) for j in range(n-1)]
            shared = [a for a in path_arcs if a in arc_set]
            marker = " ← uses cycle arcs" if shared else ""
            print(f"      {'→'.join(map(str, path))} (cycle arcs: {shared}){marker}")
        break

# ===== PART 5: n=4, the 2-path contribution =====
print("\n" + "="*70)
print("PART 5: UNDERSTANDING THE +2 CONTRIBUTION")
print("="*70)

print("""
For n=4 with one 3-cycle on {a,b,c} and a fourth vertex d:
  - The transitive part gives H=1 (the unique path through all vertices in score order)
  - The 3-cycle adds 2 more paths — the "cycle-traversing" paths

How does this work? If the cycle is a→b→c→a and d beats a,b,c:
  Transitive path: d→a→b→c (or some order)
  Cycle paths: d→c→a→b and d→b→c→a (the two ways to traverse the cycle)

Wait, that's 3 paths total? Yes: H = 1 + 2 = 3.

But the "1" is the path that doesn't use the cycle at all?
No — in a tournament, ALL paths use arcs. The question is whether
the path TRAVERSES the cycle (uses 2 consecutive cycle arcs).

Let me check: for the cycle a→b→c→a:
  - Path d→a→b→c: uses arcs (a,b) and (b,c) from the cycle. TWO cycle arcs.
  - Path d→c→a→b: uses arc (c,a) and (a,b). TWO cycle arcs.
  - Path d→b→c→a: uses arcs (b,c) and (c,a). TWO cycle arcs.

All three paths use exactly 2 cycle arcs! So the +2 is not "2 cycle paths",
it's something more subtle.
""")

# Detailed analysis for t3=1 tournament at n=4
n = 4
for bits in range(1 << (n*(n-1)//2)):
    adj = tournament_adj(n, bits)
    cycles = get_directed_3cycles(n, adj)
    H = compute_H(n, adj)
    if len(cycles) == 1 and H == 3:
        paths = get_ham_paths(n, adj)
        vs, arcs = cycles[0]
        arc_set = set(arcs)
        print(f"  Example: bits={bits}")
        for v in range(n):
            out_neighbors = [u for u in range(n) if u != v and adj[v][u]]
            print(f"    Vertex {v} beats: {out_neighbors}")
        print(f"  Cycle: {arcs}")
        print(f"  All {H} Hamiltonian paths:")
        for path in paths:
            path_arcs = [(path[j], path[j+1]) for j in range(n-1)]
            cycle_arcs_used = [a for a in path_arcs if a in arc_set]
            non_cycle = [a for a in path_arcs if a not in arc_set]
            print(f"    {'→'.join(map(str, path))}: cycle arcs={cycle_arcs_used}, other={non_cycle}")
        break

# ===== PART 6: Cycle reversal at n=4 — the pairing =====
print("\n" + "="*70)
print("PART 6: THE PAIRING — Paths in T vs paths in T' (reversed cycle)")
print("="*70)

n = 4
for bits in range(1 << (n*(n-1)//2)):
    adj = tournament_adj(n, bits)
    cycles = get_directed_3cycles(n, adj)
    H = compute_H(n, adj)
    if len(cycles) == 1 and H == 3:
        vs, arcs = cycles[0]
        adj_rev = reverse_cycle(adj, vs, arcs)
        H_rev = compute_H(n, adj_rev)
        paths_orig = set(get_ham_paths(n, adj))
        paths_rev = set(get_ham_paths(n, adj_rev))

        common = paths_orig & paths_rev
        only_orig = paths_orig - paths_rev
        only_rev = paths_rev - paths_orig

        print(f"  Original (1 cycle, H={H}):")
        for p in sorted(paths_orig):
            print(f"    {'→'.join(map(str, p))}")
        print(f"  Reversed (T', H={H_rev}):")
        cycles_rev = get_directed_3cycles(n, adj_rev)
        for p in sorted(paths_rev):
            print(f"    {'→'.join(map(str, p))}")
        print(f"\n  Common paths: {len(common)}")
        for p in sorted(common):
            print(f"    {'→'.join(map(str, p))}")
        print(f"  Only in original: {len(only_orig)}")
        for p in sorted(only_orig):
            print(f"    {'→'.join(map(str, p))}")
        print(f"  Only in reversed: {len(only_rev)}")
        for p in sorted(only_rev):
            print(f"    {'→'.join(map(str, p))}")

        print(f"\n  H(T) - H(T') = {H - H_rev}")
        print(f"  |only_orig| - |only_rev| = {len(only_orig) - len(only_rev)}")
        print(f"  Cycles in T': {[(sorted(vs), arcs) for vs, arcs in cycles_rev]}")
        break

# ===== PART 7: General pattern — transitive base + cycle contributions =====
print("\n" + "="*70)
print("PART 7: TRANSITIVE BASE + CYCLE CONTRIBUTION DECOMPOSITION")
print("="*70)

n = 5
m = n*(n-1)//2

print(f"  For each n=5 tournament, decompose H into contributions:")
print(f"  H = 1 (base) + Σ_C 2 (each cycle) + Σ_{{C,C'}} 4 (disjoint pairs)")
print()

for bits in [0b0000000000, 0b0000000001, 0b0000000011, 0b1010101010, 0b1111111111]:
    adj = tournament_adj(n, bits)
    H = compute_H(n, adj)
    cycles = get_directed_3cycles(n, adj)

    # Also get 5-cycles
    t5 = 0
    for combo in combinations(range(n), 5):
        for perm in permutations(combo):
            ok = True
            for idx in range(5):
                if not adj[perm[idx]][perm[(idx+1)%5]]:
                    ok = False
                    break
            if ok:
                t5 += 1
    t5 //= 5

    t3 = len(cycles)
    expected = 1 + 2*(t3 + t5)
    print(f"  bits={bits}: H={H}, t₃={t3}, t₅={t5}, 1+2(t₃+t₅)={expected}, match={H==expected}")

print("\n" + "="*70)
print("DONE — Proof direction analysis")
print("="*70)

print("""
PROOF STRATEGY SUMMARY:

The formula H = IP(G, 2) can be understood through the OCF
(Odd Cycle Collection Formula):

H(T) = Σ_{Ω ∈ OddCycCovers(T)} (-1)^{n - comp(Ω)}

The key insight is that each INDEPENDENT SET of pairwise
vertex-disjoint odd cycles S = {C₁, ..., C_k} contributes
exactly 2^k to this sum.

This is because:
1. The cycles in S partition some of the vertices into disjoint cycles
2. Each cycle C_i of length ℓ_i contributes a factor of 2 to the path
   count (the two ways to "unwind" the cycle into a path segment)
3. The disjointness ensures these factors multiply independently

The factor 2 per cycle is related to the fact that a directed
odd cycle of length ℓ has exactly 2 Hamiltonian paths (the two
directions of traversal), and ℓ is odd, so reversing contributes
with the same sign.

FORMAL PROOF DIRECTION:
Use the OCF expansion and group terms by the maximal independent
set they contain. Show that the contribution of each independent
set S is exactly 2^|S| by:
(a) Identifying a bijection between Hamiltonian paths and
    (independent set, remaining path) pairs
(b) Showing the "remaining path" factor contributes 1 (the base)
(c) The "independent set" factor contributes 2^|S| (binary choice per cycle)
""")
