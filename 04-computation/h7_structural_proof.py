#!/usr/bin/env python3
"""
h7_structural_proof.py — opus-2026-03-14-S71g

GOAL: Prove H=7 impossible for ALL tournaments.

H=7 requires:
  I(Ω, 2) = 7
  ⟺ α₀=1, α₁=3, α_k=0 for k≥2
  ⟺ Ω = K₃ (complete graph on 3 vertices)
  ⟺ Exactly 3 odd cycles, all pairwise sharing a vertex

Strategy: Show that in any tournament, 3 pairwise-conflicting odd cycles
force the existence of additional odd cycles (contradicting α₁=3).

Key sub-questions:
1. Can a tournament have exactly 3 odd cycles (t₃+t₅+t₇+...=3)?
2. If yes, can all 3 pairwise share vertices?
3. What additional cycles does the "3 pairwise conflicting" structure create?
"""

from itertools import permutations, combinations
from collections import defaultdict
from math import comb

def make_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_hp(A, n):
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if A[v][u]:
                    dp[S | (1 << u)][u] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

def count_directed_cycles(A, n, length):
    """Count directed cycles of given length in tournament A."""
    count = 0
    for combo in combinations(range(n), length):
        for perm in permutations(combo):
            is_cycle = True
            for i in range(length):
                if not A[perm[i]][perm[(i+1) % length]]:
                    is_cycle = False
                    break
            if is_cycle:
                count += 1
    return count // length  # each cycle counted 'length' times

def find_all_odd_cycles(A, n, max_length=None):
    """Find all directed odd cycles. Return list of (frozenset of vertices)."""
    if max_length is None:
        max_length = n
    cycles = []
    for length in range(3, max_length+1, 2):
        for combo in combinations(range(n), length):
            for perm in permutations(combo):
                # Check if it's a valid cycle starting at min vertex
                if perm[0] != min(combo):
                    continue
                is_cycle = True
                for i in range(length):
                    if not A[perm[i]][perm[(i+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    cycles.append(frozenset(combo))
    # Remove duplicates (same vertex set, different starting point)
    return list(set(cycles))

def find_all_odd_cycles_detailed(A, n, max_length=None):
    """Find all directed odd cycles with vertex lists."""
    if max_length is None:
        max_length = n
    cycles = []
    seen = set()
    for length in range(3, max_length+1, 2):
        for combo in combinations(range(n), length):
            for perm in permutations(combo):
                is_cycle = True
                for i in range(length):
                    if not A[perm[i]][perm[(i+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    # Canonical form: start at min vertex
                    min_idx = perm.index(min(combo))
                    canon = tuple(perm[min_idx:] + perm[:min_idx])
                    if canon not in seen:
                        seen.add(canon)
                        cycles.append(canon)
    return cycles

def build_conflict_graph(cycles):
    """Build conflict graph: cycles sharing a vertex are adjacent."""
    n_cyc = len(cycles)
    adj = [[False]*n_cyc for _ in range(n_cyc)]
    for i in range(n_cyc):
        si = set(cycles[i]) if isinstance(cycles[i], tuple) else cycles[i]
        for j in range(i+1, n_cyc):
            sj = set(cycles[j]) if isinstance(cycles[j], tuple) else cycles[j]
            if si & sj:  # share a vertex
                adj[i][j] = adj[j][i] = True
    return adj

def independence_poly(adj, n):
    """Compute independence polynomial of graph with adjacency matrix adj."""
    alpha = [0] * (n + 1)
    for mask in range(1 << n):
        vertices = [i for i in range(n) if mask & (1 << i)]
        independent = True
        for i in range(len(vertices)):
            for j in range(i+1, len(vertices)):
                if adj[vertices[i]][vertices[j]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            alpha[len(vertices)] += 1
    return alpha

# ============================================================
# Part 1: At each n, find tournaments with exactly t₃+t₅+...=3
# ============================================================
print("=" * 70)
print("TOURNAMENTS WITH EXACTLY 3 ODD CYCLES")
print("=" * 70)

for n in range(3, 8):
    total_edges = n * (n - 1) // 2
    num_t = 2 ** total_edges

    if num_t > 2**21:
        print(f"\n  n={n}: too large for exhaustive, sampling...")
        import random
        random.seed(42)
        sample = 200000
        found = 0
        for _ in range(sample):
            bits = random.randint(0, num_t - 1)
            A = make_tournament(bits, n)
            cycles = find_all_odd_cycles_detailed(A, n)
            if len(cycles) == 3:
                found += 1
                adj = build_conflict_graph(cycles)
                alpha = independence_poly(adj, 3)
                H = count_hp(A, n)
                vsets = [set(c) for c in cycles]
                all_conflict = all(vsets[i] & vsets[j] for i in range(3) for j in range(i+1, 3))
                if found <= 5 or all_conflict:
                    print(f"    Found! H={H}, cycles={cycles}")
                    print(f"      All pairwise conflicting? {all_conflict}")
                    print(f"      alpha={alpha}, I(Omega,2)={sum(alpha[k]*2**k for k in range(len(alpha)))}")
        print(f"    Total with exactly 3 cycles: {found}/{sample}")
        continue

    found_3 = []
    for bits in range(num_t):
        A = make_tournament(bits, n)
        cycles = find_all_odd_cycles_detailed(A, n)
        if len(cycles) == 3:
            H = count_hp(A, n)
            adj = build_conflict_graph(cycles)
            alpha = independence_poly(adj, 3)
            vsets = [set(c) for c in cycles]
            all_conflict = all(vsets[i] & vsets[j] for i in range(3) for j in range(i+1, 3))
            found_3.append((bits, H, cycles, all_conflict, alpha))

    print(f"\n  n={n}: {len(found_3)} tournaments with exactly 3 odd cycles")
    if found_3:
        all_conflicting = [x for x in found_3 if x[3]]
        not_conflicting = [x for x in found_3 if not x[3]]
        print(f"    All pairwise conflicting: {len(all_conflicting)}")
        print(f"    Not all conflicting: {len(not_conflicting)}")
        for bits, H, cycles, ac, alpha in found_3[:10]:
            Ival = sum(alpha[k]*2**k for k in range(len(alpha)))
            print(f"    H={H}, all_conf={ac}, alpha={alpha}, I(Ω,2)={Ival}")
            print(f"      cycles: {cycles}")

# ============================================================
# Part 2: Structural analysis of 3 pairwise-conflicting triangles
# ============================================================
print(f"\n{'='*70}")
print("CAN 3 TRIANGLES ALL PAIRWISE SHARE VERTICES?")
print(f"{'='*70}")

# 3 triangles C1, C2, C3 with Ci ∩ Cj ≠ ∅ for all i≠j.
# Case 1: All share a common vertex v.
#   C1 = {v, a1, b1}, C2 = {v, a2, b2}, C3 = {v, a3, b3}
#   with {a_i, b_i} pairwise disjoint → need 7 vertices.
#   But at n=7, having vertex v in 3 triangles and 6 other vertices
#   creates MANY more triangles.
# Case 2: Pairwise sharing but no common vertex.
#   C1∩C2={a}, C1∩C3={b}, C2∩C3={c} with a≠b≠c.
#   C1={a,b,x}, C2={a,c,y}, C3={b,c,z}, need 6 vertices min.

print("\nCase 1: Common vertex v in all 3 triangles")
print("  C1={v,a,b}, C2={v,c,d}, C3={v,e,f}, n≥7")
print("  At n=7: vertices {v,a,b,c,d,e,f}")
print("  The tournament on {a,b,c,d,e,f} has C(6,3)=20 triples")
print("  Each triple forms 0 or 2 directed 3-cycles")
print("  Need: v is in EXACTLY 3 triangles, and {a,b,c,d,e,f} has 0 triangles")
print("  i.e., the sub-tournament on {a,b,c,d,e,f} is transitive")

# Verify: if the 6 non-v vertices form a transitive tournament,
# how many triangles does v participate in?
# v→a,v→b means {v,a,b} is a triangle iff b→a→v... no.
# Triangle on {v,x,y}: need v→x→y→v (or v→y→x→v).
# If v→x and v→y: triangle iff x→y→v or y→x→v. But v→y, so y→v impossible.
# Wait: v→x means the arc goes from v to x. For a triangle v→x→y→v,
# we need y→v. So the vertices y that beat v are the in-neighbors of v.

# In a tournament on {v}∪S where |S|=6:
# Let v have out-degree d in S (d in-neighbors, 6-d out-neighbors).
# Triangle through v: v→x→y→v requires x∈out(v), y∈in(v), x→y.
# Number of such triangles = # edges from out(v) to in(v) in the sub-tournament on S.

# If out(v) = {a,b,c}, in(v) = {d,e,f}:
# Sub-tournament on S has some edges from {a,b,c} to {d,e,f}.
# Among the 9 pairs (a,d),(a,e),...,(c,f), each has one direction.
# Triangle count through v = # pairs where out-vertex → in-vertex = t_v.
# We need t_v = 3.

# Can {a,b,c,d,e,f} be arranged transitively while having exactly 3 edges
# from {a,b,c} to {d,e,f}? The transitive tournament on 6 vertices has a
# specific structure.

print("\n  Checking: transitive sub-tournament + controlled v-degree")
for d_out in range(7):
    d_in = 6 - d_out
    # In transitive tournament on S, the number of edges from any
    # set of size d_out to its complement of size d_in depends on
    # which vertices are in which set.
    # Min edges out→in: pick top d_out vertices as out.
    # Max edges out→in: pick bottom d_out vertices.
    # Since transitive has total order 0<1<...<5 where i→j iff i<j,
    # edges from A to B = #{(a,b): a∈A, b∈B, a<b}

    if d_out == 0 or d_in == 0:
        continue

    # Try all subsets of size d_out from {0,1,2,3,4,5}
    possible_tv = set()
    for subset in combinations(range(6), d_out):
        complement = [x for x in range(6) if x not in subset]
        edges = sum(1 for a in subset for b in complement if a < b)
        possible_tv.add(edges)

    if 3 in possible_tv:
        print(f"  d_out={d_out}: possible t_v through v = {sorted(possible_tv)}, includes 3 ✓")
    else:
        print(f"  d_out={d_out}: possible t_v through v = {sorted(possible_tv)}, NO 3")

# ============================================================
# Part 3: The critical question — does transitive-on-6 + 3 triangles
# through v FORCE additional non-v triangles?
# ============================================================
print(f"\n{'='*70}")
print("CRITICAL TEST: 7-VERTEX TOURNAMENTS WITH EXACTLY 3 TRIANGLES")
print(f"{'='*70}")

# At n=7, enumerate tournaments where the non-v sub (6 vertices) is transitive
# and v has exactly 3 triangles through it.
# Then check: does this tournament have exactly 3 odd cycles total?
# (Must also count 5-cycles and 7-cycles!)

n = 7
# Fix vertex 0 as v. Vertices 1-6 form transitive tournament: i→j iff i<j.
# v's arcs: for each vertex i in 1..6, either v→i or i→v.
# 2^6 = 64 choices.

print("\n  Fixing v=0, vertices 1-6 in transitive order")
results = []

for v_bits in range(64):
    # Build tournament
    A = [[0]*7 for _ in range(7)]
    # Transitive on 1-6: i→j iff i<j
    for i in range(1, 7):
        for j in range(i+1, 7):
            A[i][j] = 1
    # v=0's arcs
    for i in range(6):
        if v_bits & (1 << i):
            A[0][i+1] = 1  # v→(i+1)
        else:
            A[i+1][0] = 1  # (i+1)→v

    # Count triangles through v
    out_v = [i+1 for i in range(6) if v_bits & (1 << i)]
    in_v = [i+1 for i in range(6) if not (v_bits & (1 << i))]
    # Triangle v→x→y→v: x∈out_v, y∈in_v, x→y (i.e., x<y in transitive)
    t_v = sum(1 for x in out_v for y in in_v if A[x][y])

    # Non-v triangles: in transitive tournament, there are 0.
    # So total t₃ = t_v.

    if t_v != 3:
        continue

    # Count ALL odd cycles (3, 5, 7)
    all_cycles = find_all_odd_cycles_detailed(A, 7)
    t3 = sum(1 for c in all_cycles if len(c) == 3)
    t5 = sum(1 for c in all_cycles if len(c) == 5)
    t7 = sum(1 for c in all_cycles if len(c) == 7)
    total = len(all_cycles)

    H = count_hp(A, 7)

    # Check if all cycles pairwise share vertices
    if total == 3:
        adj = build_conflict_graph(all_cycles)
        vsets = [set(c) for c in all_cycles]
        all_conf = all(vsets[i] & vsets[j] for i in range(3) for j in range(i+1, 3))
    else:
        all_conf = None

    results.append((v_bits, t3, t5, t7, total, H, all_cycles if total <= 10 else f"{total} cycles", all_conf))

print(f"  Tournaments with t₃=3 through v: {len(results)}")
for vb, t3, t5, t7, total, H, cycles, ac in results:
    out_v = [i+1 for i in range(6) if vb & (1 << i)]
    in_v = [i+1 for i in range(6) if not (vb & (1 << i))]
    print(f"    out(v)={out_v}, in(v)={in_v}")
    print(f"    t₃={t3}, t₅={t5}, t₇={t7}, total={total}, H={H}", end="")
    if ac is not None:
        print(f", all_conflicting={ac}")
    else:
        print()
    if total <= 10 and isinstance(cycles, list):
        for c in cycles:
            print(f"      cycle: {c}")

# ============================================================
# Part 4: General structural argument
# ============================================================
print(f"\n{'='*70}")
print("STRUCTURAL ARGUMENT: WHY 3 PAIRWISE-CONFLICTING CYCLES FORCE MORE")
print(f"{'='*70}")

# Key insight: if C1, C2, C3 are 3-cycles all sharing vertex v,
# then v has out-degree 3 and in-degree 3 (in the sub-tournament
# restricted to the cycle vertices). The 3 out-neighbors of v
# form a tournament, and the 3 in-neighbors form a tournament.
# The inter-connections create 5-cycles!

# Let's check: with 3 triangles through v, all on disjoint other vertices,
# do 5-cycles necessarily arise?

print("\n  Case: C1={v,a,b}, C2={v,c,d}, C3={v,e,f} with")
print("  v→a→b→v, v→c→d→v, v→e→f→v")
print("  So out(v)={a,c,e}, in(v)={b,d,f}")
print("  The sub-tournament on {a,b,c,d,e,f} must have a→b→v as part of cycle,")
print("  and also arcs between {a,c,e} and {b,d,f}.")
print()
print("  A 5-cycle through v: v→a→?→?→?→v uses 4 other vertices.")
print("  E.g., v→a→c→d→b→v (if a→c, c→d, d→b, b→v all exist).")
print("  We know: b→v (yes), c→d (yes, part of C2: v→c→d→v).")
print("  So we need: a→c and d→b.")
print()
print("  Testing: which 5-cycles can exist?")

# Enumerate at n=7: fix v=0, C1={0,1,2}, C2={0,3,4}, C3={0,5,6}
# Cycle arcs: 0→1→2→0, 0→3→4→0, 0→5→6→0
# So: A[0][1]=A[1][2]=A[2][0]=1, A[0][3]=A[3][4]=A[4][0]=1, A[0][5]=A[5][6]=A[6][0]=1
# Remaining arcs to decide: {1,2} vs {3,4}, {1,2} vs {5,6}, {3,4} vs {5,6}
# That's 4+4+4 = 12 arcs → 2^12 = 4096 tournaments

n = 7
base_arcs = {(0,1), (1,2), (2,0), (0,3), (3,4), (4,0), (0,5), (5,6), (6,0)}

# The remaining undecided pairs:
pairs_12_34 = [(i,j) for i in [1,2] for j in [3,4]]
pairs_12_56 = [(i,j) for i in [1,2] for j in [5,6]]
pairs_34_56 = [(i,j) for i in [3,4] for j in [5,6]]
all_remaining = pairs_12_34 + pairs_12_56 + pairs_34_56

print(f"\n  Fixing 3 triangles through v=0: {{0,1,2}}, {{0,3,4}}, {{0,5,6}}")
print(f"  Remaining undecided pairs: {len(all_remaining)}")
print(f"  Total completions: {2**len(all_remaining)}")

exactly_3_total = 0
h7_count = 0

for bits in range(2**len(all_remaining)):
    A = [[0]*7 for _ in range(7)]
    # Set base arcs
    for (i,j) in base_arcs:
        A[i][j] = 1
    # Set remaining arcs
    for idx, (i,j) in enumerate(all_remaining):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    # Find all odd cycles
    all_cycles = find_all_odd_cycles_detailed(A, 7)
    total = len(all_cycles)

    if total == 3:
        exactly_3_total += 1
        adj = build_conflict_graph(all_cycles)
        vsets = [set(c) for c in all_cycles]
        all_conf = all(vsets[i] & vsets[j] for i in range(3) for j in range(i+1, 3))
        H = count_hp(A, 7)
        if H == 7:
            h7_count += 1
        print(f"    FOUND: total=3, H={H}, all_conf={all_conf}")
        for c in all_cycles:
            print(f"      {c}")

if exactly_3_total == 0:
    print(f"\n  RESULT: No completion has exactly 3 odd cycles!")
    print(f"  The 3 base triangles ALWAYS force additional odd cycles.")

print(f"\n  H=7 count: {h7_count}")

# Count by number of odd cycles
cycle_counts = defaultdict(int)
for bits in range(2**len(all_remaining)):
    A = [[0]*7 for _ in range(7)]
    for (i,j) in base_arcs:
        A[i][j] = 1
    for idx, (i,j) in enumerate(all_remaining):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    all_cycles = find_all_odd_cycles_detailed(A, 7)
    cycle_counts[len(all_cycles)] += 1

print(f"\n  Distribution of total odd cycle count:")
for k in sorted(cycle_counts.keys()):
    print(f"    {k} cycles: {cycle_counts[k]} tournaments")

print(f"\n  MINIMUM odd cycles with 3 forced triangles: {min(cycle_counts.keys())}")

# ============================================================
# Part 5: What creates the extra cycles?
# ============================================================
print(f"\n{'='*70}")
print("MECHANISM: WHY 3 TRIANGLES FORCE EXTRA CYCLES")
print(f"{'='*70}")

# Find the tournament(s) achieving the minimum cycle count
min_cycles = min(cycle_counts.keys())
print(f"\n  Minimum cycle count: {min_cycles}")
print(f"  Examples:")

count_shown = 0
for bits in range(2**len(all_remaining)):
    A = [[0]*7 for _ in range(7)]
    for (i,j) in base_arcs:
        A[i][j] = 1
    for idx, (i,j) in enumerate(all_remaining):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    all_cycles = find_all_odd_cycles_detailed(A, 7)
    if len(all_cycles) == min_cycles and count_shown < 3:
        count_shown += 1
        H = count_hp(A, 7)
        t3 = sum(1 for c in all_cycles if len(c) == 3)
        t5 = sum(1 for c in all_cycles if len(c) == 5)
        t7 = sum(1 for c in all_cycles if len(c) == 7)
        print(f"\n    bits={bits}: t₃={t3}, t₅={t5}, t₇={t7}, H={H}")
        for c in all_cycles:
            print(f"      {c}")

        # Check conflict structure
        adj = build_conflict_graph(all_cycles)
        vsets = [set(c) for c in all_cycles]
        print(f"      Pairwise conflicts:")
        for i in range(len(all_cycles)):
            for j in range(i+1, len(all_cycles)):
                shared = vsets[i] & vsets[j]
                print(f"        C{i}∩C{j} = {shared} {'CONFLICT' if shared else 'INDEPENDENT'}")

        alpha = independence_poly(adj, len(all_cycles))
        Ival = sum(alpha[k]*2**k for k in range(len(alpha)))
        print(f"      α = {alpha[:6]}, I(Ω,2) = {Ival}")

print(f"\n{'='*70}")
print("ANALYSIS COMPLETE")
print(f"{'='*70}")
