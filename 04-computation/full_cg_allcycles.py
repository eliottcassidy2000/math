#!/usr/bin/env python3
"""
full_cg_allcycles.py — opus-2026-03-14-S76

DEFINITIVE TEST: The OCF uses ALL directed odd cycles (not chordless).

Previous script showed:
- "All cycles" matches H at n=5 for T0, T1, T2, T4 but NOT T3
- Need to check the cycle counting for T3

The issue: "find_all_odd_cycles_with_chords" counts each VERTEX SET once,
but there can be MULTIPLE directed cycles on the same vertex set.
For a 5-cycle on {a,b,c,d,e}, there could be two different cyclic orderings
that both form directed cycles.

CORRECTION: We should count DIRECTED cycles, not vertex sets.
A directed odd cycle is a cyclic permutation class of vertices.
On k vertices, there are (k-1)!/2 possible cyclic orderings
(since a→b→c→a and a→c→b→a are different cycles).

Actually, the conflict graph has vertices = directed odd cycles,
where "directed odd cycle" means an oriented cycle.
But two orientations of the same set of vertices use the same vertices,
so they are adjacent in the conflict graph.

Wait: in a tournament, on a set of 3 vertices {a,b,c}, there is
at most ONE directed 3-cycle (either a→b→c→a or a→c→b→a, not both).

On 5 vertices, there can be multiple directed 5-cycles with different
orderings. But we only care about the VERTEX SET for disjointness.

The question: does the OCF count each vertex set once (as one "cycle node")
or does it count each directed cycle separately?

Let me check by comparing H = I(Ω, 2).
"""

from itertools import combinations, permutations
import random

def random_tournament(n):
    adj = [0] * n
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i] |= (1 << j)
            else:
                adj[j] |= (1 << i)
    return adj

def ham_paths(adj, n):
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for idx in range(n-1):
            if not (adj[perm[idx]] & (1 << perm[idx+1])):
                ok = False
                break
        if ok:
            count += 1
    return count

def find_directed_odd_cycles_by_vertex_set(adj, n):
    """For each subset of odd size, check if it supports a directed cycle.
    Count each vertex set at most once."""
    cycles = []
    for length in range(3, n + 1, 2):
        for combo in combinations(range(n), length):
            verts = list(combo)
            fs = frozenset(combo)
            for perm in permutations(verts):
                is_cycle = True
                for idx in range(length):
                    if not (adj[perm[idx]] & (1 << perm[(idx+1) % length])):
                        is_cycle = False
                        break
                if is_cycle:
                    cycles.append(fs)
                    break  # one per vertex set
    return cycles

def count_directed_cycles_on_subset(adj, verts):
    """Count how many distinct directed cycles exist on the given vertex set."""
    length = len(verts)
    verts = list(verts)
    count = 0
    seen = set()
    for perm in permutations(verts):
        is_cycle = True
        for idx in range(length):
            if not (adj[perm[idx]] & (1 << perm[(idx+1) % length])):
                is_cycle = False
                break
        if is_cycle:
            # Normalize: smallest rotation
            canon = min(perm[i:] + perm[:i] for i in range(length))
            canon = tuple(canon)
            if canon not in seen:
                seen.add(canon)
                count += 1
    return count

def compute_independence_poly(cycles_list):
    """Independence polynomial of conflict graph."""
    m = len(cycles_list)
    if m == 0:
        return [1]

    conflict = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if len(cycles_list[i] & cycles_list[j]) > 0:
                conflict[i][j] = True
                conflict[j][i] = True

    alpha = {}
    alpha[0] = 1
    max_k = 0

    for mask in range(1, 1 << m):
        bits = []
        temp = mask
        idx = 0
        while temp:
            if temp & 1:
                bits.append(idx)
            temp >>= 1
            idx += 1

        is_indep = True
        for i in range(len(bits)):
            for j in range(i+1, len(bits)):
                if conflict[bits[i]][bits[j]]:
                    is_indep = False
                    break
            if not is_indep:
                break

        if is_indep:
            k = len(bits)
            alpha[k] = alpha.get(k, 0) + 1
            max_k = max(max_k, k)

    return [alpha.get(k, 0) for k in range(max_k + 1)]

# ====================================================================
print("=" * 70)
print("PART 1: DETAILED n=5 ANALYSIS")
print("=" * 70)
print()

random.seed(99)
n = 5
for trial in range(10):
    adj = random_tournament(n)
    H = ham_paths(adj, n)

    # Method A: vertex-set cycles (each vertex set counted once)
    cycles_vs = find_directed_odd_cycles_by_vertex_set(adj, n)
    alpha_vs = compute_independence_poly(cycles_vs)
    i2_vs = sum(2**k * a for k, a in enumerate(alpha_vs))
    im1_vs = sum((-1)**k * a for k, a in enumerate(alpha_vs))

    # Count individual directed cycles per vertex set
    cycle_details = {}
    for fs in cycles_vs:
        nc = count_directed_cycles_on_subset(adj, fs)
        cycle_details[len(fs)] = cycle_details.get(len(fs), [])
        cycle_details[len(fs)].append((fs, nc))

    print(f"T{trial}: H={H}")
    print(f"  Vertex-set method: I(2)={i2_vs}, match={H==i2_vs}")
    for length in sorted(cycle_details.keys()):
        for fs, nc in cycle_details[length]:
            if nc > 1:
                print(f"    {length}-cycle on {set(fs)}: {nc} directed cycles!")
    print(f"  Alpha: {alpha_vs}")

    # If mismatch, check what's missing
    if H != i2_vs:
        gap = H - i2_vs
        print(f"  GAP = {gap}")
        # The gap might be explained by multiple directed cycles on same vertex set
        total_extra = 0
        for length in sorted(cycle_details.keys()):
            for fs, nc in cycle_details[length]:
                if nc > 1:
                    total_extra += nc - 1
        print(f"  Extra directed cycles (beyond 1 per vertex set): {total_extra}")
    print()

# ====================================================================
print("=" * 70)
print("PART 2: COUNTING DIRECTED CYCLES CORRECTLY")
print("=" * 70)
print()
print("On a set of k vertices, a directed cycle visits each vertex once")
print("in cyclic order. The number of distinct directed k-cycles on k")
print("labeled vertices is (k-1)! (fixing one vertex, permuting the rest).")
print("But in a tournament, at most SOME of these are actually directed.")
print()
print("For 3 vertices: (3-1)! = 2 possible directed 3-cycles.")
print("A tournament on 3 vertices has exactly 1 or 0 directed 3-cycles.")
print("(It's 1 iff the tournament is NOT transitive.)")
print()
print("For 5 vertices: (5-1)! = 24 possible directed 5-cycles.")
print("But as directed CYCLES: we identify cyclic rotations,")
print("so (5-1)! / 5... no, (5-1)! already fixes one vertex.")
print("Actually: there are (5-1)!/2 = 12 distinct unoriented cycles,")
print("each giving 1 directed cycle. So up to 12 directed 5-cycles.")
print()

# Count directed 5-cycles at n=5
n = 5
random.seed(42)
for trial in range(5):
    adj = random_tournament(n)
    H = ham_paths(adj, n)

    # Count 3-cycles
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i]&(1<<j)) and (adj[j]&(1<<k)) and (adj[k]&(1<<i)):
                    c3 += 1
                elif (adj[i]&(1<<k)) and (adj[k]&(1<<j)) and (adj[j]&(1<<i)):
                    c3 += 1

    # Count ALL directed 5-cycles (as canonical rotations)
    verts = list(range(n))
    c5_directed = set()
    for perm in permutations(verts):
        is_cycle = True
        for idx in range(5):
            if not (adj[perm[idx]] & (1 << perm[(idx+1) % 5])):
                is_cycle = False
                break
        if is_cycle:
            # Canonical: smallest rotation
            canon = min(tuple(perm[i:] + perm[:i]) for i in range(5))
            c5_directed.add(canon)

    print(f"T{trial}: H={H}, #3-cycles={c3}, #directed-5-cycles={len(c5_directed)}")
    # All 5-cycles are on the SAME vertex set {0,1,2,3,4}
    # So they're all mutually adjacent in CG(T)
    # α₁ = c3 + len(c5_directed)
    # All 5-cycles share vertices with each other AND with all 3-cycles
    # that use any of {0,1,2,3,4} = all 3-cycles!
    # So all cycles are mutually adjacent → CG is complete → α₂ = 0
    a1 = c3 + len(c5_directed)
    print(f"  All cycles are on overlapping vertices at n=5")
    print(f"  α₁={a1}, α₂=0 (complete conflict graph)")
    print(f"  I(2) = 1 + 2·{a1} = {1+2*a1}")
    print(f"  Match? {H == 1+2*a1}")
    print()

# ====================================================================
print()
print("=" * 70)
print("PART 3: THE CORRECT OCF — WHAT ARE THE VERTICES OF Ω(T)?")
print("=" * 70)
print()
print("HYPOTHESIS 1: Vertices = vertex sets supporting a directed odd cycle")
print("HYPOTHESIS 2: Vertices = individual directed odd cycles")
print()
print("Testing at n=5:")
print("At n=5, every 3-cycle is on a unique vertex set (C(5,3)=10 possible),")
print("and there's only one 5-vertex set.")
print("So Hyp 1 and Hyp 2 agree for 3-cycles.")
print("For 5-cycles: Hyp 1 has 1 node, Hyp 2 has multiple nodes.")
print()

# Build specific tournament
n = 5
# All arcs: 0→1, 0→2, 0→3, 0→4, 1→2, 2→3, 3→4, 4→1, 1→3, 2→4
# This is QR_5
adj_qr5 = [0]*5
adj_qr5[0] = (1<<1)|(1<<2)|(1<<3)|(1<<4)  # 0 beats all
# 1→2, 1→3
adj_qr5[1] = (1<<2)|(1<<3)
# 2→3, 2→4
adj_qr5[2] = (1<<3)|(1<<4)
# 3→4, 3→... nah let me build QR_5 properly

# QR_5: quadratic residues mod 5 are {1,4}
# i→j iff j-i mod 5 ∈ {1,4}
adj_qr5 = [0]*5
for i in range(5):
    for j in range(5):
        if i != j and (j-i) % 5 in {1, 4}:
            adj_qr5[i] |= (1 << j)

H_qr5 = ham_paths(adj_qr5, 5)
print(f"QR_5: H = {H_qr5}")

# Score sequence
scores = [bin(a).count('1') for a in adj_qr5]
print(f"  Scores: {scores} (should be all 2)")

# Count all cycle types
c3_list = []
for i in range(5):
    for j in range(i+1, 5):
        for k in range(j+1, 5):
            if (adj_qr5[i]&(1<<j)) and (adj_qr5[j]&(1<<k)) and (adj_qr5[k]&(1<<i)):
                c3_list.append(frozenset([i,j,k]))
            elif (adj_qr5[i]&(1<<k)) and (adj_qr5[k]&(1<<j)) and (adj_qr5[j]&(1<<i)):
                c3_list.append(frozenset([i,j,k]))

# 5-cycles
c5_set = set()
for perm in permutations(range(5)):
    ok = True
    for idx in range(5):
        if not (adj_qr5[perm[idx]] & (1 << perm[(idx+1)%5])):
            ok = False
            break
    if ok:
        canon = min(tuple(perm[i:]+perm[:i]) for i in range(5))
        c5_set.add(canon)

print(f"  3-cycles: {len(c3_list)}")
print(f"  Directed 5-cycles: {len(c5_set)}")
for c in sorted(c5_set):
    print(f"    {c}")

# Now: Hyp 2 (directed cycles as vertices)
all_nodes = list(c3_list) + [frozenset(range(5))] * len(c5_set)
# Wait, this overcounts. Each directed 5-cycle has the SAME vertex set.
# In the conflict graph: all 5-cycle nodes are adjacent to each other
# and to all 3-cycle nodes (since they all share vertices).
# So they form a clique with everything else.

# Hyp 1: vertex sets
vs_nodes = list(c3_list) + ([frozenset(range(5))] if c5_set else [])
alpha_vs = compute_independence_poly(vs_nodes)
i2_vs = sum(2**k * a for k, a in enumerate(alpha_vs))

# Hyp 2: directed cycles (each 5-cycle is a separate node, all on {0,1,2,3,4})
# All 5-cycle nodes share all vertices with each other and with all 3-cycle nodes
dc_nodes = list(c3_list) + [frozenset(range(5))] * len(c5_set)
alpha_dc = compute_independence_poly(dc_nodes)
i2_dc = sum(2**k * a for k, a in enumerate(alpha_dc))

print(f"\n  Hyp 1 (vertex sets): α = {alpha_vs}, I(2) = {i2_vs}")
print(f"  Hyp 2 (directed cycles): α = {alpha_dc}, I(2) = {i2_dc}")
print(f"  H = {H_qr5}")
print(f"  Hyp 1 matches? {H_qr5 == i2_vs}")
print(f"  Hyp 2 matches? {H_qr5 == i2_dc}")

# ====================================================================
print()
print("=" * 70)
print("PART 4: RE-READING THE DEFINITION MORE CAREFULLY")
print("=" * 70)
print()
print("From definitions.md:")
print("  'Conflict graph Ω(T): vertices are directed odd cycles of T'")
print()
print("'Directed odd cycles' most naturally means:")
print("Each cyclic sequence a→b→c→a is ONE directed cycle.")
print("Different orderings of the same vertex set that give different")
print("directed cycles are different vertices of Ω(T).")
print()
print("But in a tournament on 3 vertices {a,b,c}:")
print("If a→b→c→a is a directed cycle, then a→c→b→a is NOT")
print("(because we'd need c→b, b→a, a→c, which is the reverse).")
print("In a tournament, exactly one of the two orientations works.")
print("So for 3-cycles: 1 directed cycle per vertex set.")
print()
print("For 5-cycles on {a,b,c,d,e} in a tournament:")
print("Multiple distinct directed 5-cycles can exist!")
print("E.g., a→b→c→d→e→a AND a→c→e→b→d→a could both be directed cycles.")
print()

# ====================================================================
print()
print("=" * 70)
print("PART 5: EXHAUSTIVE n=5 WITH DIRECTED 5-CYCLES AS SEPARATE NODES")
print("=" * 70)
print()

n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
num_edges = len(edges)
total = 2 ** num_edges

h_match_vs = 0
h_match_dc = 0
im1_max_dc = -999
im1_violations_dc = 0

for bits in range(total):
    adj = [0] * n
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            adj[i] |= (1 << j)
        else:
            adj[j] |= (1 << i)

    H = ham_paths(adj, n)

    # 3-cycles
    c3 = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i]&(1<<j)) and (adj[j]&(1<<k)) and (adj[k]&(1<<i)):
                    c3.append(frozenset([i,j,k]))
                elif (adj[i]&(1<<k)) and (adj[k]&(1<<j)) and (adj[j]&(1<<i)):
                    c3.append(frozenset([i,j,k]))

    # Directed 5-cycles (canonical rotations)
    c5_canonical = set()
    for perm in permutations(range(n)):
        ok = True
        for idx in range(5):
            if not (adj[perm[idx]] & (1 << perm[(idx+1)%5])):
                ok = False
                break
        if ok:
            canon = min(tuple(perm[i:]+perm[:i]) for i in range(5))
            c5_canonical.add(canon)

    num_c5 = len(c5_canonical)

    # Method: vertex sets (1 node per 5-cycle vertex set)
    vs_nodes = list(c3) + ([frozenset(range(n))] if num_c5 > 0 else [])
    alpha_vs = compute_independence_poly(vs_nodes)
    i2_vs = sum(2**k * a for k, a in enumerate(alpha_vs))

    # Method: directed cycles (each 5-cycle = separate node)
    # All on same vertex set, so all adjacent to each other and all 3-cycles
    dc_nodes = list(c3) + [frozenset(range(n))] * num_c5
    alpha_dc = compute_independence_poly(dc_nodes)
    i2_dc = sum(2**k * a for k, a in enumerate(alpha_dc))
    im1_dc = sum((-1)**k * a for k, a in enumerate(alpha_dc))

    if H == i2_vs:
        h_match_vs += 1
    if H == i2_dc:
        h_match_dc += 1
    if im1_dc > 1:
        im1_violations_dc += 1
    im1_max_dc = max(im1_max_dc, im1_dc)

print(f"n=5 exhaustive ({total} tournaments):")
print(f"  H = I(2) matches (vertex-set method): {h_match_vs}/{total}")
print(f"  H = I(2) matches (directed-cycle method): {h_match_dc}/{total}")
print(f"  I(-1) max (directed-cycle): {im1_max_dc}")
print(f"  I(-1) > 1 violations (directed-cycle): {im1_violations_dc}")

# ====================================================================
print()
print("=" * 70)
print("PART 6: THE REAL FORMULA — MAYBE NOT INDEPENDENCE POLYNOMIAL?")
print("=" * 70)
print()
print("The OCF formula is: H(T) = Σ_{S independent in Ω(T)} 2^|S|")
print("This IS the independence polynomial at x=2.")
print()
print("But WAIT: the original OCF from THM-070 says:")
print("H = I(Ω(T), 2)")
print("And this was VERIFIED at small n.")
print()
print("Let me check if the verification used 3-cycles only!")
print()

# Check: at n=5, using ONLY 3-cycles
n = 5
edges_n5 = [(i,j) for i in range(n) for j in range(i+1,n)]
total_n5 = 2**len(edges_n5)

match_3only = 0
for bits in range(total_n5):
    adj = [0] * n
    for idx, (i,j) in enumerate(edges_n5):
        if bits & (1 << idx):
            adj[i] |= (1 << j)
        else:
            adj[j] |= (1 << i)

    H = ham_paths(adj, n)
    c3 = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i]&(1<<j)) and (adj[j]&(1<<k)) and (adj[k]&(1<<i)):
                    c3.append(frozenset([i,j,k]))
                elif (adj[i]&(1<<k)) and (adj[k]&(1<<j)) and (adj[j]&(1<<i)):
                    c3.append(frozenset([i,j,k]))

    alpha = compute_independence_poly(c3)
    i2 = sum(2**k * a for k, a in enumerate(alpha))

    if H == i2:
        match_3only += 1

print(f"Using ONLY 3-cycles at n=5: {match_3only}/{total_n5} match H=I(2)")
print()
print("If this is <100%, then 3-cycles alone don't suffice.")
print("If 100%, then the OCF might only need 3-cycles!")
