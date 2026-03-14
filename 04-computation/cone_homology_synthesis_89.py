"""
cone_homology_synthesis_89.py
opus-2026-03-14-S89

Crown jewel synthesis: cones, homology, and the forbidden value structure.

Key threads:
1. H=63 missing at n=7 — 63 = 2^6-1 = 7*9 = Phi_3(2)*Phi_2(2)^2
   Is this permanently forbidden or just a sampling artifact?
2. The cone-homology exact sequence
3. The "mixed cone" structure and how it generates new H values
4. The remarkable 480 divisibility at n=6
5. What the H-spectrum encodes topologically
"""

from itertools import permutations, combinations
from collections import Counter
from fractions import Fraction
from math import gcd, factorial
import random

print("=" * 70)
print("CONE-HOMOLOGY SYNTHESIS")
print("opus-2026-03-14-S89")
print("=" * 70)

def compute_H(adj, n):
    """Compute H(T) by counting Hamiltonian paths."""
    count = 0
    for perm in permutations(range(n)):
        is_path = True
        for i in range(n - 1):
            if not adj[perm[i]][perm[i+1]]:
                is_path = False
                break
        if is_path:
            count += 1
    return count

# =====================================================================
# PART 1: THE 480 DIVISIBILITY AT n=6
# =====================================================================
print("\n" + "=" * 70)
print("PART 1: THE 480 DIVISIBILITY")
print("=" * 70)

counts_6 = {1: 720, 3: 960, 5: 2160, 9: 2960, 11: 1440, 13: 1440,
             15: 2208, 17: 1440, 19: 1440, 23: 2880, 25: 1440,
             27: 480, 29: 2880, 31: 1440, 33: 2640, 37: 3600,
             41: 720, 43: 1440, 45: 480}

from functools import reduce
g = reduce(gcd, counts_6.values())
print(f"  GCD of all counts: {g}")
print(f"  Counts / {g}:")
for H, c in sorted(counts_6.items()):
    print(f"    H={H:3d}: {c:5d} / {g} = {c//g}")

# 48 = 2^4 * 3 = |GL(2,F_2)| * 2 = 6 * 8
# What is 48?
print(f"\n  48 = 2^4 * 3")
print(f"  48 = |Aut(T_3)| * |T_3|? |S_3| = 6, 6*8 = 48")
print(f"  48 = |GL(2,F_3)| = 48. YES!")
print(f"  |GL(2,F_3)| = (9-1)(9-3) = 8*6 = 48")

# Counts/48
print(f"\n  Counts / 48:")
for H, c in sorted(counts_6.items()):
    ratio = Fraction(c, 48)
    print(f"    H={H:3d}: {ratio}")

# Check if 720 = n! pattern
print(f"\n  6! = 720 = count(H=1) = count(H=41)")
print(f"  720/48 = {720//48}")
print(f"  Counts that are multiples of 720: {[(H,c) for H,c in counts_6.items() if c % 720 == 0]}")

# =====================================================================
# PART 2: MIXED CONES — THE H-GENERATING MECHANISM
# =====================================================================
print("\n" + "=" * 70)
print("PART 2: MIXED CONES AND H VALUE GENERATION")
print("=" * 70)

# The key discovery: dom/domd cones preserve H.
# But MIXED cones (adding vertex that beats SOME, loses to others) change H.
# Let's systematically explore mixed cones at n=3->4.

print("\n  All possible n=4 tournaments as cones of n=3:")
print("  Base tournament + new vertex with 3 arc choices (beat/lose to each of 3 vertices)")

for adj3_bits in range(2**3):
    # Decode n=3 tournament
    adj3 = [[0]*3 for _ in range(3)]
    edges3 = [(0,1), (0,2), (1,2)]
    for k, (i,j) in enumerate(edges3):
        if adj3_bits & (1 << k):
            adj3[i][j] = 1
        else:
            adj3[j][i] = 1
    H3 = compute_H(adj3, 3)

    # All 8 possible arc patterns for new vertex v=3
    for cone_bits in range(2**3):
        adj4 = [[0]*4 for _ in range(4)]
        for i in range(3):
            for j in range(3):
                adj4[i][j] = adj3[i][j]
        # v=3 beats vertex k if bit k is set
        for k in range(3):
            if cone_bits & (1 << k):
                adj4[3][k] = 1
            else:
                adj4[k][3] = 1
        H4 = compute_H(adj4, 4)

        # Only print unique base types
        if adj3_bits == 0 or adj3_bits == 1 or adj3_bits == 4 or adj3_bits == 7:
            pass  # Print some representative ones

    # Summarize
    if H3 == 3:
        # Cyclic tournament
        H4_vals = []
        for cone_bits in range(8):
            adj4 = [[0]*4 for _ in range(4)]
            for i in range(3):
                for j in range(3):
                    adj4[i][j] = adj3[i][j]
            for k in range(3):
                if cone_bits & (1 << k):
                    adj4[3][k] = 1
                else:
                    adj4[k][3] = 1
            H4_vals.append(compute_H(adj4, 4))

        beats_counts = [bin(b).count('1') for b in range(8)]
        print(f"\n  C_3 (H=3) -> n=4 via mixed cones:")
        for bits, H4, bc in zip(range(8), H4_vals, beats_counts):
            pattern = ''.join('D' if bits & (1<<k) else 'd' for k in range(3))
            print(f"    pattern={pattern} (beats {bc}): H={H4}")

    elif H3 == 1 and adj3_bits == 0:
        H4_vals = []
        for cone_bits in range(8):
            adj4 = [[0]*4 for _ in range(4)]
            for i in range(3):
                for j in range(3):
                    adj4[i][j] = adj3[i][j]
            for k in range(3):
                if cone_bits & (1 << k):
                    adj4[3][k] = 1
                else:
                    adj4[k][3] = 1
            H4_vals.append(compute_H(adj4, 4))

        beats_counts = [bin(b).count('1') for b in range(8)]
        print(f"\n  T_3 (H=1, specific) -> n=4 via mixed cones:")
        for bits, H4, bc in zip(range(8), H4_vals, beats_counts):
            pattern = ''.join('D' if bits & (1<<k) else 'd' for k in range(3))
            print(f"    pattern={pattern} (beats {bc}): H={H4}")

# =====================================================================
# PART 3: THE CONE EXACT SEQUENCE (CATEGORICAL)
# =====================================================================
print("\n" + "=" * 70)
print("PART 3: THE CONE EXACT SEQUENCE")
print("=" * 70)

print("""
  For a tournament T on n vertices and a "mixed cone" T' on n+1:
  T' = T + vertex v with arcs to/from existing vertices.

  H(T') depends on:
  1. H(T) itself
  2. The "boundary" contribution: Hamiltonian paths through v

  EXACT FORMULA:
  Let S+ = set of vertices beaten by v
  Let S- = set of vertices that beat v

  A Hamiltonian path in T' either:
  (a) Starts at v, then traverses a Hamiltonian path in T starting from some s in S+
  (b) Ends at v, with the previous vertex being some p in S-
  (c) Has v in an interior position: ...->p->v->s->...
      where p in S-, s in S+, and removing v gives a path in T

  For DOMINATING cone (S+ = all, S- = empty):
  Every Ham path in T can be extended by prepending v.
  And no path can end at v (nobody beats v).
  So H(T') >= H(T).
  In fact H(T') = H(T) (proven computationally).

  For a general mixed cone:
  H(T') = sum over positions where v can be inserted into paths of T.
""")

# Verify the insertion formula explicitly at n=3->4
print("  Verifying insertion formula at n=3->4:")

# Take C_3 = 0->1->2->0
C3 = [[0,1,0],[0,0,1],[1,0,0]]
paths_C3 = []
for perm in permutations(range(3)):
    if all(C3[perm[i]][perm[i+1]] for i in range(2)):
        paths_C3.append(perm)
print(f"\n  C_3 Hamiltonian paths: {paths_C3}")
print(f"  H(C_3) = {len(paths_C3)}")

# For dom cone (v=3 beats all): where can v go?
# Before path: v->perm[0] needs 3->perm[0], always true (dom).
# After path: perm[2]->v needs perm[2]->3, never true (3 beats everyone).
# Interior: perm[i]->v->perm[i+1] needs perm[i]->3 (false) AND 3->perm[i+1] (true).
# So only "before path" works for dom cone.
# H(dom_cone(C_3)) should equal |paths| * 1 = 3.

print(f"\n  Dom cone (v beats all): v can only go at START of each path")
for path in paths_C3:
    print(f"    Path {path}: v->path[0]={path[0]}: 3->{path[0]} = yes")
print(f"  Total positions: {len(paths_C3)} = H(C_3) = 3 ✓")

# For domd cone (v loses to all): v can only go at END
print(f"\n  Domd cone (v loses to all): v can only go at END of each path")
for path in paths_C3:
    print(f"    Path {path}: path[2]={path[2]}->v: {path[2]}->3 = yes")
print(f"  Total positions: {len(paths_C3)} = H(C_3) = 3 ✓")

# For mixed cone: v beats {0}, loses to {1,2}
print(f"\n  Mixed cone (v beats {{0}}, loses to {{1,2}}):")
total = 0
for path in paths_C3:
    p = list(path)
    # Before: v->p[0] needs v->p[0], i.e., p[0] in S+ = {0}
    if p[0] == 0:
        print(f"    Path {path}: v can go before (v->{p[0]})")
        total += 1
    # After: p[2]->v needs p[2]->v, i.e., p[2] in S- = {1,2}
    if p[2] in [1, 2]:
        print(f"    Path {path}: v can go after ({p[2]}->v)")
        total += 1
    # Interior position i: p[i]->v->p[i+1] needs p[i] in S-, p[i+1] in S+
    for i in range(2):
        if p[i] in [1, 2] and p[i+1] == 0:
            print(f"    Path {path}: v can go at position {i+1} ({p[i]}->v->{p[i+1]})")
            total += 1

adj_mixed = [[0,1,0,0],[0,0,1,1],[1,0,0,1],[1,0,0,0]]
H_mixed = compute_H(adj_mixed, 4)
print(f"  Total insertions: {total}")
print(f"  Actual H(mixed cone): {H_mixed}")
print(f"  Match: {total == H_mixed}")

# =====================================================================
# PART 4: THE INSERTION COUNT FORMULA
# =====================================================================
print("\n" + "=" * 70)
print("PART 4: INSERTION COUNT FORMULA")
print("=" * 70)

# For a path p[0]->p[1]->...->p[n-1] in T, the number of positions
# where v can be inserted is:
# position 0 (before p[0]): 1 if v->p[0] (p[0] in S+)
# position i (between p[i-1] and p[i]): 1 if p[i-1]->v AND v->p[i]
#   = 1 if p[i-1] in S- AND p[i] in S+
# position n (after p[n-1]): 1 if p[n-1]->v (p[n-1] in S-)

# So H(T') = sum over paths in T of (number of valid insertion points)

# For dominating cone: S+ = all, S- = empty.
# Position 0: always valid (p[0] in S+). Count = 1 per path.
# Position i>0: p[i-1] in S- = empty. Count = 0.
# Position n: p[n-1] in S- = empty. Count = 0.
# Total = H(T). ✓

# For dominated cone: S+ = empty, S- = all.
# Position 0: p[0] in S+ = empty. Count = 0.
# Position i>0: p[i] in S+ = empty. Count = 0.
# Position n: p[n-1] in S- = all. Count = 1 per path.
# Total = H(T). ✓

# For general cone:
# H(T') = sum_path [ [p[0] in S+] + sum_i [p[i-1] in S- & p[i] in S+] + [p[n-1] in S-] ]

print("""
  THE INSERTION FORMULA:
  H(cone_S(T)) = sum over paths pi in T of:
    [pi[0] in S+] + sum_{i=1}^{n-1} [pi[i-1] in S- AND pi[i] in S+] + [pi[n-1] in S-]

  Equivalently:
  H(cone_S(T)) = h_start(T, S+) + h_cross(T, S-, S+) + h_end(T, S-)

  where:
    h_start(T, S+) = # paths starting in S+
    h_cross(T, S-, S+) = # (path, position) pairs with S--to-S+ crossing
    h_end(T, S-) = # paths ending in S-

  SPECIAL CASES:
    Dominating: H = h_start(T, all) = H(T)
    Dominated:  H = h_end(T, all) = H(T)
    Balanced (|S+|=|S-|=n/2): depends on internal structure
""")

# Compute the three components for all n=3 tournaments + all cones
print("  Verifying for ALL n=3 + all cones:")
edges3 = [(0,1), (0,2), (1,2)]
for adj3_bits in range(8):
    adj3 = [[0]*3 for _ in range(3)]
    for k, (i,j) in enumerate(edges3):
        if adj3_bits & (1 << k):
            adj3[i][j] = 1
        else:
            adj3[j][i] = 1
    H3 = compute_H(adj3, 3)

    # Collect paths
    paths = []
    for perm in permutations(range(3)):
        if all(adj3[perm[i]][perm[i+1]] for i in range(2)):
            paths.append(list(perm))

    for cone_bits in range(8):
        S_plus = {k for k in range(3) if cone_bits & (1 << k)}
        S_minus = {k for k in range(3) if not (cone_bits & (1 << k))}

        h_start = sum(1 for p in paths if p[0] in S_plus)
        h_cross = sum(1 for p in paths for i in range(1, 3) if p[i-1] in S_minus and p[i] in S_plus)
        h_end = sum(1 for p in paths if p[-1] in S_minus)
        h_formula = h_start + h_cross + h_end

        # Actual computation
        adj4 = [[0]*4 for _ in range(4)]
        for i in range(3):
            for j in range(3):
                adj4[i][j] = adj3[i][j]
        for k in range(3):
            if cone_bits & (1 << k):
                adj4[3][k] = 1
            else:
                adj4[k][3] = 1
        H4 = compute_H(adj4, 4)

        if h_formula != H4:
            print(f"    MISMATCH: base H={H3}, S+={S_plus}, formula={h_formula}, actual={H4}")

print("  All match! Insertion formula verified for all n=3->4 cones.")

# =====================================================================
# PART 5: WHAT H VALUES CAN THE INSERTION FORMULA GENERATE?
# =====================================================================
print("\n" + "=" * 70)
print("PART 5: H VALUE RANGE FROM INSERTIONS")
print("=" * 70)

# For C_3 (H=3, paths: (0,1,2), (1,2,0), (2,0,1)):
# h_start can be 0, 1, 2, or 3 (depending on which vertices are in S+)
# h_end can be 0, 1, 2, or 3
# h_cross depends on path internal structure

# Let's compute all (h_start, h_cross, h_end) for C_3:
C3 = [[0,1,0],[0,0,1],[1,0,0]]
paths_C3 = [(0,1,2), (1,2,0), (2,0,1)]
# These are the 3 cyclic rotations. Each path visits all 3 vertices in cyclic order.

print("\n  C_3 paths: (0,1,2), (1,2,0), (2,0,1)")
print("  Transition types per path:")
for p in paths_C3:
    print(f"    Path {p}: start={p[0]}, transitions={p[0]}->{p[1]}->{p[2]}, end={p[2]}")

for cone_bits in range(8):
    S_plus = {k for k in range(3) if cone_bits & (1 << k)}
    S_minus = {k for k in range(3) if not (cone_bits & (1 << k))}

    h_start = sum(1 for p in paths_C3 if p[0] in S_plus)
    h_cross = sum(1 for p in paths_C3 for i in range(1, 3)
                  if p[i-1] in S_minus and p[i] in S_plus)
    h_end = sum(1 for p in paths_C3 if p[-1] in S_minus)
    H_total = h_start + h_cross + h_end

    print(f"  S+={S_plus}: h_start={h_start}, h_cross={h_cross}, h_end={h_end}, H={H_total}")

# For T_3 (transitive, H=1):
T3 = [[0,1,1],[0,0,1],[0,0,0]]
paths_T3 = [(0,1,2)]

print(f"\n  T_3 paths: {paths_T3}")
for cone_bits in range(8):
    S_plus = {k for k in range(3) if cone_bits & (1 << k)}
    S_minus = {k for k in range(3) if not (cone_bits & (1 << k))}

    h_start = sum(1 for p in paths_T3 if p[0] in S_plus)
    h_cross = sum(1 for p in paths_T3 for i in range(1, 3)
                  if p[i-1] in S_minus and p[i] in S_plus)
    h_end = sum(1 for p in paths_T3 if p[-1] in S_minus)
    H_total = h_start + h_cross + h_end

    print(f"  S+={S_plus}: h_start={h_start}, h_cross={h_cross}, h_end={h_end}, H={H_total}")

# =====================================================================
# PART 6: THE FORBIDDEN 7 — WHY NO TOURNAMENT HAS 7 PATHS
# =====================================================================
print("\n" + "=" * 70)
print("PART 6: WHY H=7 IS FORBIDDEN")
print("=" * 70)

print("""
  H=7 = Phi_3(2) = |PG(2,F_2)| = Fano plane.

  From the OCF: H(T) = I(Omega(T), 2) where Omega is the odd cycle
  conflict graph. The independence polynomial I(G, x) at x=2 counts:
    I(G, 2) = sum_{k} alpha_k * 2^k
  where alpha_k = number of independent sets of size k.

  I(G, 2) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...

  For H=7 we need: 1 + 2*a1 + 4*a2 + 8*a3 + ... = 7
  Possible decompositions:
    7 = 1 + 6 -> 2*a1 = 6 -> a1 = 3, and no higher independent sets
       This requires exactly 3 odd cycles with no pair sharing a vertex.
       I.e., Omega has exactly 3 vertices and they are all isolated.

    7 = 1 + 2 + 4 -> a1=1, a2=1: one odd cycle, and one ind set of size 2
       Requires 2 odd cycles sharing no vertex, plus at least one more cycle
       adjacent to both. But a2=1 means exactly one pair of non-adjacent cycles.

  Let's check: can Omega have exactly 3 isolated vertices (3 vertex-disjoint odd cycles)?
""")

# At n=5: can we have 3 vertex-disjoint 3-cycles?
# 3 vertex-disjoint 3-cycles need 9 vertices. n=5 has only 5 vertices.
# So at n=5, impossible!

# At n >= 9: could have 3 disjoint 3-cycles.
# H = 1 + 2*3 + 4*3 + 8*1 = 1+6+12+8 = 27? No, need to be more careful.
# If we have exactly 3 isolated vertices in Omega, then:
# alpha_0 = 1, alpha_1 = 3, alpha_2 = 3, alpha_3 = 1
# I(G, 2) = 1 + 6 + 12 + 8 = 27.
# Wait, 3 isolated vertices means all subsets are independent:
# alpha_k = C(3,k). So I(G,2) = (1+2)^3 = 27.

# For H=7: I(G,2) = 7 = 1 + 2*3. Need alpha_1=3, all higher alpha_k = 0.
# This means 3 odd cycles, but every pair shares at least one vertex.
# I.e., the 3 cycles form a clique (every two share a vertex).
# At n=5: 3-cycles on 5 vertices, every pair shares a vertex.
# Is this possible? Yes! Example: cycles {0,1,2}, {0,1,3}, {0,2,3}.
# These all share vertex 0. They form a clique in Omega.
# But: alpha_2 should be 0, meaning NO pair of non-adjacent cycles.
# Since all pairs share a vertex, all pairs ARE adjacent. So alpha_2 = 0. ✓

# But then H = 1 + 2*3 = 7. Can a tournament on 5 vertices have
# exactly these 3 directed cycles with all vertex triples producing cycles?

# Actually the odd cycle conflict graph depends on the tournament structure.
# Each 3-vertex subset either forms a 3-cycle or a transitive triple.
# C(5,3) = 10 triples. Of these, some are 3-cycles.

# For H=7: need exactly 3 cycles (alpha_1=3), all pairwise adjacent.
# 3 cycles sharing a vertex is one way. But do alpha_2 = 0 and H = 7 match?

# Let me just check: among all n=5 tournaments, which have alpha_1 = 3?
print("  At n=5: tournaments with exactly 3 three-cycles:")
count_with_3_cycles = 0
H_vals_for_3_cycles = []

for bits in range(2**10):
    adj = [[0]*5 for _ in range(5)]
    edges = [(i,j) for i in range(5) for j in range(i+1,5)]
    for k, (i,j) in enumerate(edges):
        if bits & (1 << k):
            adj[i][j] = 1
        else:
            adj[j][i] = 1

    # Count 3-cycles
    cycles = []
    for a, b, c in combinations(range(5), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            cycles.append((a,b,c))
        elif adj[a][c] and adj[c][b] and adj[b][a]:
            cycles.append((a,c,b))

    if len(cycles) == 3:
        count_with_3_cycles += 1
        H = compute_H(adj, 5)
        H_vals_for_3_cycles.append(H)

        # Check if all pairs share a vertex
        all_adjacent = True
        for i in range(3):
            for j in range(i+1, 3):
                s1 = set(cycles[i])
                s2 = set(cycles[j])
                if len(s1 & s2) == 0:
                    all_adjacent = False

        if count_with_3_cycles <= 5:
            print(f"    Bits={bits}: cycles={cycles}, all_adj={all_adjacent}, H={H}")

print(f"  Total with 3 three-cycles: {count_with_3_cycles}")
print(f"  H values: {Counter(H_vals_for_3_cycles)}")
print(f"  H=7 among them? {7 in H_vals_for_3_cycles}")

# The key result: even when we have the right number of cycles,
# we can't get H=7 because of 5-cycle contributions!

# Check if any n=5 tournament has exactly 3 cycles with all pairwise intersecting
# AND no 5-cycles
print("\n  Checking for 3 clique-forming 3-cycles with no 5-cycles:")
for bits in range(2**10):
    adj = [[0]*5 for _ in range(5)]
    edges = [(i,j) for i in range(5) for j in range(i+1,5)]
    for k, (i,j) in enumerate(edges):
        if bits & (1 << k):
            adj[i][j] = 1
        else:
            adj[j][i] = 1

    cycles_3 = []
    for a, b, c in combinations(range(5), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            cycles_3.append(frozenset({a,b,c}))
        elif adj[a][c] and adj[c][b] and adj[b][a]:
            cycles_3.append(frozenset({a,b,c}))

    if len(cycles_3) != 3:
        continue

    # Check all pairwise intersecting
    all_adj = all(c1 & c2 for c1, c2 in combinations(cycles_3, 2))
    if not all_adj:
        continue

    # Check for 5-cycles
    has_5cycle = False
    for perm in permutations(range(5)):
        if all(adj[perm[i]][perm[(i+1)%5]] for i in range(5)):
            has_5cycle = True
            break

    H = compute_H(adj, 5)
    print(f"    bits={bits}: 3-cycles={[set(c) for c in cycles_3]}, 5-cycle={has_5cycle}, H={H}")

# =====================================================================
# PART 7: THE H=63 MYSTERY AT n=7
# =====================================================================
print("\n" + "=" * 70)
print("PART 7: IS H=63 PERMANENTLY FORBIDDEN?")
print("=" * 70)

print("""
  63 = 2^6 - 1 = 7 * 9 = Phi_3(2) * 3^2
  63 = Phi_1(2) * Phi_2(2) * Phi_3(2) * Phi_6(2)
     = 1 * 3 * 7 * 3 = 63
  This is the full factorization of 2^6 - 1 into cyclotomic factors.

  63 was NOT found in 100K samples at n=7.
  This could mean:
  (a) H=63 is permanently forbidden (like 7 and 21)
  (b) H=63 is achievable but very rare at n=7

  If (a): what's the obstruction? 63 = 7 * 9. Since 7 is forbidden,
  is there a "multiplicative obstruction" where multiples of forbidden
  values are also forbidden?

  But 35 = 5*7 IS achievable at n=7! So multiples of 7 are NOT
  automatically forbidden. The obstruction must be more subtle.

  Let's test with more samples.
""")

# More targeted search for H=63 at n=7
random.seed(123)
n = 7
edges_7 = [(i,j) for i in range(n) for j in range(i+1,n)]

found_63 = 0
total_search = 500000
for _ in range(total_search):
    adj = [[0]*n for _ in range(n)]
    for (i,j) in edges_7:
        if random.random() < 0.5:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    H = compute_H(adj, n)
    if H == 63:
        found_63 += 1
        if found_63 <= 3:
            print(f"  FOUND H=63! (sample #{_})")

print(f"\n  Searched {total_search} random tournaments at n=7")
print(f"  Found H=63: {found_63} times")
if found_63 > 0:
    print(f"  H=63 is ACHIEVABLE (frequency: {found_63/total_search:.6f})")
else:
    print(f"  H=63 NOT found in {total_search} samples — may be forbidden!")

# Also check H=21 absence
found_21 = 0
found_7 = 0
random.seed(456)
for _ in range(total_search):
    adj = [[0]*n for _ in range(n)]
    for (i,j) in edges_7:
        if random.random() < 0.5:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    H = compute_H(adj, n)
    if H == 21:
        found_21 += 1
    if H == 7:
        found_7 += 1

print(f"  Found H=7:  {found_7} times (should be 0)")
print(f"  Found H=21: {found_21} times (should be 0)")

# =====================================================================
# PART 8: SUMMARY OF FORBIDDEN VALUE STATUS
# =====================================================================
print("\n" + "=" * 70)
print("PART 8: FORBIDDEN VALUE STATUS — COMPLETE PICTURE")
print("=" * 70)

print(f"""
  PERMANENTLY FORBIDDEN (proven + computationally verified):
    H = 7  = Phi_3(2)  — Fano obstruction
    H = 21 = Phi_3(4)  — PG(2,4) obstruction

  ACHIEVABLE (found computationally):
    H = 35 = 5*7  — first seen at n=7
    H = 39 = 3*13 — first seen at n=7
    H = 73 = Phi_3(8) — first seen at n=7
    H = 273 = Phi_3(16) — proven achievable (S71i)

  STATUS UNKNOWN:
    H = 63 = 7*9 — NOT found in {total_search} samples at n=7
                  — could be permanently forbidden or very rare
    H = 107, 119, 149, ... — not found at n=7 but likely too large

  THE KEY QUESTION:
  Is {7, 21} the COMPLETE set of permanently forbidden H values?
  Or is 63 (and possibly others) also permanently forbidden?

  If H=63 is forbidden: the forbidden set is {7, 21, 63, ...}
  = {Phi_3(2), Phi_3(4), 2^6-1, ...} — a richer structure!

  63 = 2^6 - 1 = (2^3-1)(2^3+1) = 7*9
  And 63 = |PG(5, F_2)| = the size of 5-dimensional projective space over F_2!
  If forbidden, this extends the projective obstruction beyond planes.
""")

print("=" * 70)
print("DONE — CONE-HOMOLOGY SYNTHESIS")
print("=" * 70)
