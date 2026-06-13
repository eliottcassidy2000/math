"""
vitali_geometric_structure.py -- kind-pasteur-2026-03-13-S61

We now know: dc7 depends on the FULL local data:
  (internal S-arcs) + (S-outside projection) + (outside-outside arcs)

The outside-outside arcs form a TOURNAMENT on 3 vertices.
A tournament on 3 vertices is either a 3-cycle (C3) or a transitive triple (T3).

Key question: does the XX tournament type (C3 vs T3) determine dc7?
If so, this connects to the CIRCULAR structure of the complement.

Deeper: the Vitali atom S sits inside the full tournament T.
The "fiber" over S consists of:
  (a) S-internal: the (1,1,2,2) tournament on S (up to iso: just 1 type)
  (b) S-X arcs: 4×3 = 12 arcs, each directed S→X or X→S
  (c) X-X arcs: 3 arcs forming a tournament on the outside vertices

The Vitali ATOM reverses (a) while keeping (b) and (c) fixed.
Since dc7 depends on all three components, the "hidden dimension" is the
joint structure of the S→X projection AND the X-X tournament.

Hypothesis: dc7 = f(internal) ⊗ g(projection) ⊗ h(XX_type)
where the three factors interact multiplicatively.

Let's also explore: what are the "2" and "1" in the score (1,1,2,2)?
These vertices behave asymmetrically under the Vitali operation.
The winners (score 2) and losers (score 1) may have different "coupling"
to the outside, which determines dc7.
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def lambda_graph(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v:
                    continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

def reverse_subtournament(A, n, subset):
    B = A.copy()
    for i in subset:
        for j in subset:
            if i != j:
                B[i][j] = A[j][i]
    return B

def sub_scores(A, n, subset):
    k = len(subset)
    return tuple(sorted([sum(A[subset[i]][subset[j]] for j in range(k) if i != j) for i in range(k)]))

def count_directed_k_cycles(A, n, k):
    Ak = np.linalg.matrix_power(A, k)
    return int(np.trace(Ak)) // k

def xx_type(A, outside):
    """Classify outside tournament as C3 or T3"""
    x0, x1, x2 = outside
    # Check for 3-cycle
    c = (A[x0][x1]*A[x1][x2]*A[x2][x0] + A[x0][x2]*A[x2][x1]*A[x1][x0])
    return "C3" if c > 0 else "T3"

n = 7
total_bits = n * (n-1) // 2

print("=" * 60)
print("GEOMETRIC STRUCTURE OF VITALI dc7")
print("=" * 60)

np.random.seed(42)
examples = []

for trial in range(8000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        if not np.array_equal(L, lambda_graph(B, n)):
            continue

        S = set(subset)
        S_sorted = sorted(S)
        outside = sorted(k for k in range(n) if k not in S)

        c7_A = count_directed_k_cycles(A, n, 7)
        c7_B = count_directed_k_cycles(B, n, 7)
        dc7 = c7_B - c7_A

        # Internal tournament classification
        losers = [v for v in S_sorted if sum(A[v][w] for w in S_sorted if w != v) == 1]
        winners = [v for v in S_sorted if sum(A[v][w] for w in S_sorted if w != v) == 2]

        # S→X projection scores
        s_to_x_scores = {}
        for s in S_sorted:
            s_to_x_scores[s] = sum(A[s][x] for x in outside)

        # X→S projection scores
        x_to_s_scores = {}
        for x in outside:
            x_to_s_scores[x] = sum(A[x][s] for s in S_sorted)

        # XX type
        xxt = xx_type(A, outside)

        # "Coupling vector": for each S-vertex, its outdegree to outside
        loser_out = tuple(sorted(s_to_x_scores[l] for l in losers))
        winner_out = tuple(sorted(s_to_x_scores[w] for w in winners))

        examples.append({
            'bits': bits,
            'subset': subset,
            'dc7': dc7,
            'losers': losers,
            'winners': winners,
            'outside': outside,
            'xxt': xxt,
            'loser_out': loser_out,
            'winner_out': winner_out,
            's_to_x': s_to_x_scores,
            'x_to_s': x_to_s_scores,
            'A': A.copy(),
        })
        break

    if len(examples) >= 500:
        break

print(f"Collected {len(examples)} examples")

# Analysis 1: XX type vs dc7
print(f"\n--- XX Tournament Type vs dc7 ---")
xx_dc7 = defaultdict(list)
for ex in examples:
    xx_dc7[ex['xxt']].append(ex['dc7'])

for xxt, dc7s in sorted(xx_dc7.items()):
    print(f"  {xxt}: dc7 dist = {dict(sorted(Counter(dc7s).items()))}, n={len(dc7s)}")

# Analysis 2: Loser/Winner coupling vs dc7
print(f"\n--- Coupling Structure vs dc7 ---")
coupling_dc7 = defaultdict(list)
for ex in examples:
    key = (ex['loser_out'], ex['winner_out'], ex['xxt'])
    coupling_dc7[key].append(ex['dc7'])

print(f"Coupling groups: {len(coupling_dc7)}")
ambig = sum(1 for dc7s in coupling_dc7.values() if len(set(dc7s)) > 1)
print(f"Ambiguous for dc7: {ambig}")

for key, dc7s in sorted(coupling_dc7.items()):
    dc7_dist = Counter(dc7s)
    if len(dc7s) >= 2:
        print(f"  L_out={key[0]}, W_out={key[1]}, XX={key[2]}: dc7={dict(sorted(dc7_dist.items()))} ({len(dc7s)})")

# Analysis 3: Detailed coupling — which S-vertex connects to which outside vertex
print(f"\n--- Detailed Arc Pattern vs dc7 ---")

# For each S-vertex s, its arc pattern to outside: (A[s][x0], A[s][x1], A[s][x2])
# But we need to account for the internal S-structure too.
# Key insight: in (1,1,2,2), the two losers each beat exactly 1 other S-vertex,
# and the two winners each beat exactly 2 other S-vertices.
# The "missing" arc is which S-vertex the loser beats and which the winner loses to.

# The internal structure: losers L1, L2; winners W1, W2
# Each loser beats exactly 1 of {L2, W1, W2} and loses to the other 2.
# Each winner beats exactly 2 of {L1, L2, W_other} and loses to 1.

# In fact, in a (1,1,2,2) tournament:
# - Both winners beat both losers (4 arcs: W1→L1, W1→L2, W2→L1, W2→L2)
# - Between the two winners: one beats the other (1 arc)
# - Between the two losers: one beats the other (1 arc)
# Total: 6 arcs. Let's verify.

print(f"\n  Verifying (1,1,2,2) structure:")
for ex in examples[:5]:
    A = ex['A']
    L = ex['losers']
    W = ex['winners']
    print(f"    S={ex['subset']}, L={L}, W={W}")
    # Check W→L arcs
    for w in W:
        for l in L:
            print(f"      A[{w}->{l}]={A[w][l]}", end="  ")
    print()
    # W-W arc
    print(f"      A[{W[0]}->{W[1]}]={A[W[0]][W[1]]}", end="  ")
    # L-L arc
    print(f"      A[{L[0]}->{L[1]}]={A[L[0]][L[1]]}")

# So the (1,1,2,2) is NOT just W→L. Let me check.
wl_pattern = Counter()
for ex in examples:
    A = ex['A']
    L, W = ex['losers'], ex['winners']
    pattern = (A[W[0]][L[0]], A[W[0]][L[1]], A[W[1]][L[0]], A[W[1]][L[1]])
    wl_pattern[pattern] += 1

print(f"\n  W->L patterns: {dict(wl_pattern)}")

# The key arc structure:
# Not all 4 W→L arcs point the same way. Let me enumerate the actual structure.
print(f"\n  Full internal structure types:")
internal_types = Counter()
for ex in examples:
    A = ex['A']
    S_sorted = sorted(ex['subset'])
    # Represent as adjacency within S
    key = tuple(A[i][j] for i in S_sorted for j in S_sorted if i != j)
    internal_types[key] += 1

print(f"  Distinct internal structures: {len(internal_types)}")
for key, cnt in internal_types.most_common(5):
    print(f"    count={cnt}: {key}")

# Analysis 4: The "7-cycle counting" interpretation
print(f"\n{'='*60}")
print("7-CYCLE DECOMPOSITION UNDER VITALI FLIP")
print(f"{'='*60}")

# A 7-cycle on 7 vertices visits each vertex exactly once.
# The 4 S-vertices and 3 outside vertices are visited in some interleaved order.
# The cycle uses some arcs within S and some between S and outside.

# Under the Vitali flip, only S-internal arcs change.
# So dc7 = #{7-cycles using reversed S-arcs favorably} - #{using them unfavorably}

# How many S-internal arcs can a 7-cycle use?
# A 7-cycle visits 7 vertices. The 4 S-vertices appear in some positions.
# Between consecutive S-vertices in the cycle (with possible outside vertices between),
# if two S-vertices are adjacent in the cycle, that's an S-S arc.

# Let's count directly: for each 7-cycle, how many S-internal arcs does it use?

# First, enumerate all directed 7-cycles in T and T'
cycle_analysis = []
for ex_idx in range(min(20, len(examples))):
    ex = examples[ex_idx]
    A = ex['A']
    B = reverse_subtournament(A, n, list(ex['subset']))
    S = set(ex['subset'])

    # Count 7-cycles and classify by number of S-internal arcs
    # A directed Hamiltonian cycle is a cyclic permutation.
    # Enumerate all permutations of {0,...,6} starting from 0
    c7_by_internal = Counter()
    c7_by_internal_B = Counter()

    from itertools import permutations as perms
    for perm in perms(range(1, n)):
        cycle = (0,) + perm
        # Check if this is a directed cycle in A
        valid_A = all(A[cycle[i]][cycle[(i+1) % n]] == 1 for i in range(n))
        valid_B = all(B[cycle[i]][cycle[(i+1) % n]] == 1 for i in range(n))

        if valid_A or valid_B:
            # Count S-internal arcs
            s_arcs = sum(1 for i in range(n)
                        if cycle[i] in S and cycle[(i+1) % n] in S)
            if valid_A:
                c7_by_internal[s_arcs] += 1
            if valid_B:
                c7_by_internal_B[s_arcs] += 1

    # Normalize by cycle count (each cycle counted 7 times)
    c7_A_check = sum(v for v in c7_by_internal.values()) // 7
    c7_B_check = sum(v for v in c7_by_internal_B.values()) // 7

    cycle_analysis.append({
        'dc7': ex['dc7'],
        'c7_A': c7_A_check,
        'c7_B': c7_B_check,
        'A_by_internal': dict(c7_by_internal),
        'B_by_internal': dict(c7_by_internal_B),
        'xxt': ex['xxt'],
    })

    if ex_idx < 5:
        print(f"\nExample {ex_idx}: dc7={ex['dc7']}, XX={ex['xxt']}")
        print(f"  A: c7={c7_A_check}, by S-arcs: {dict(c7_by_internal)}")
        print(f"  B: c7={c7_B_check}, by S-arcs: {dict(c7_by_internal_B)}")

# Summary: which S-arc counts contribute to dc7?
print(f"\n--- S-arc count distribution in 7-cycles ---")
all_A_arcs = Counter()
all_B_arcs = Counter()
for ca in cycle_analysis:
    for k, v in ca['A_by_internal'].items():
        all_A_arcs[k] += v
    for k, v in ca['B_by_internal'].items():
        all_B_arcs[k] += v

print(f"  In A: {dict(sorted(all_A_arcs.items()))}")
print(f"  In B: {dict(sorted(all_B_arcs.items()))}")

# Key: which S-arc counts appear?
# A 7-cycle through 4 S-vertices can have 0, 1, 2, 3, or 4 consecutive S-S arcs.
# The maximum is 4 (if all S-vertices are consecutive: ...S1-S2-S3-S4-X1-X2-X3-...)
# The minimum is 0 (if no two S-vertices are adjacent in the cycle)

# Under reversal, cycles with k S-internal arcs:
# - k=0: no change (no reversed arcs used) => same cycle in both A and B
# - k>0: some arcs reversed, cycle may exist in one but not other

# So dc7 = sum_k (#{7-cycles in B with k S-arcs} - #{7-cycles in A with k S-arcs})
# The k=0 terms cancel! Only k>=1 terms contribute.
# This means dc7 depends only on cycles that TRAVERSE the Vitali atom.

print(f"\n  k=0 cycles equal? ", end="")
k0_equal = all(
    ca['A_by_internal'].get(0, 0) == ca['B_by_internal'].get(0, 0)
    for ca in cycle_analysis
)
print(f"{'YES' if k0_equal else 'NO'}")

# For cycles with k S-arcs: under reversal, k arcs flip.
# If k is EVEN: the same number of forward/backward arcs => cycle may survive
# If k is ODD: parity changes => cycle direction reverses

# Actually: reversing ALL S-arcs means if a cycle traverses S-arcs
# a1->a2, a2->a3, ..., it becomes a2->a1, a3->a2, ...
# This reverses the DIRECTION of the S-internal sub-path.

# Hmm, this is more nuanced. Let me think about the cycle structure.
# A 7-cycle visits 4 S-vertices. In the cycle, these appear at positions i1 < i2 < i3 < i4.
# Between consecutive S-vertices, there may be 0 or more outside vertices.

# The S-internal arcs are those where cycle[j] ∈ S AND cycle[j+1] ∈ S.
# Under reversal: A'[s1][s2] = A[s2][s1], so arc s1->s2 becomes s2->s1.

# For a cycle to survive the reversal at positions where both endpoints are in S:
# the arc needs to be reversed (s1->s2 becomes s2->s1), but the cycle needs s2->s1.
# So the cycle survives iff at every S-S adjacency, the original had s2->s1
# (i.e., the cycle was already going "against" the original arc direction).
# After reversal, s2->s1 stays as s1->s2 which matches the cycle direction s1->s2? No...

# Let me be more precise:
# Original cycle: ... -> s1 -> s2 -> ... requires A[s1][s2] = 1
# After reversal: A'[s1][s2] = A[s2][s1]
# So the cycle requires A[s2][s1] = 1 after reversal, i.e., s2 beat s1 originally.
# This means: a cycle that used arc s1->s2 in A survives in B only if s2 also beats s1,
# which is impossible (tournament). So NO cycle with an S-S arc survives the flip!

# Instead: cycles in B with S-S arcs were NOT present in A (they needed the reversed direction).

# Therefore:
# Let f(k) = #{7-cycles in A with exactly k S-S arcs}
# Let g(k) = #{7-cycles in B with exactly k S-S arcs}
# Cycles with k=0: f(0) = g(0) (not affected by reversal)
# Cycles with k>0: f(k) and g(k) are DISJOINT (no cycle with S-S arcs exists in both)

# So dc7 = c7_B - c7_A = sum_k (g(k) - f(k)) = sum_{k>0} (g(k) - f(k))

# Moreover, for each k>0, g(k) = #{cycles that USE the B-arcs at exactly k S-S positions}
# These are cycles that used the REVERSED arcs in their k S-S positions.

# Can we say more? The relationship between f(k) and g(k)?
# For a specific cycle visiting S-vertices in order s_{i1}, s_{i2}, s_{i3}, s_{i4}:
# If exactly k of the 4 "gaps" between consecutive S-vertices in the cycle are direct S-S arcs,
# then under reversal, those k arcs flip direction.
# The cycle with REVERSED S-S arcs and SAME other arcs may or may not be valid.

# Actually, let me just verify the "no shared cycle" claim:
print(f"\n  Verifying: no 7-cycle with S-S arcs exists in both A and B:")
shared_with_ss = 0
for ca in cycle_analysis:
    for k in ca['A_by_internal']:
        if k > 0 and k in ca['B_by_internal']:
            # Check if same SPECIFIC cycles appear
            pass  # Can't tell from counts alone
print(f"  (Need per-cycle tracking for exact verification)")

# Instead, count the parity structure
print(f"\n--- Parity of S-arc count in 7-cycles ---")
for ca in cycle_analysis[:5]:
    print(f"  dc7={ca['dc7']}: A arcs={ca['A_by_internal']}, B arcs={ca['B_by_internal']}")

# Analysis 5: The "permutation pattern" of S-vertices in the cycle
print(f"\n{'='*60}")
print("PERMUTATION PATTERN OF S IN 7-CYCLES")
print(f"{'='*60}")

# In a 7-cycle on vertices {0,...,6}, the 4 S-vertices appear in some
# relative order (up to cyclic rotation and reflection).
# The number of distinct patterns: S(4,7) = ?
# Actually: place 4 vertices among 7 cycle positions.
# The gaps between consecutive S-vertices (cyclically) sum to 7.
# Gaps can be 1 (S-S adjacent) or 2+ (with outside vertices between).
# With 4 S-vertices and 3 outside, gaps sum to 7.
# Gap sizes are positive integers summing to 7 with 4 terms.
# Gap = 1 means S-S adjacent (counts as S-arc).
# Gap >= 2 means at least one outside vertex between.

# Partitions of 7 into 4 positive parts:
# 1+1+1+4, 1+1+2+3, 1+2+2+2, 2+2+2+1 (same as 1+2+2+2)
# Actually: 1+1+1+4, 1+1+2+3, 1+2+2+2
# Cyclic patterns: these are necklaces.

# With 4 gaps summing to 7:
# Compositions (ordered) with multiplicities:
# (1,1,1,4): 4 arrangements (which gap is 4)
# (1,1,2,3): 4!/2! = 12 arrangements
# (1,2,2,2): 4 arrangements
# (2,2,2,1): same as (1,2,2,2)
# But modulo cyclic rotation (4 rotations):
# (1,1,1,4): 1 necklace
# (1,1,2,3): 3 necklaces (1,1,2,3), (1,2,1,3), (1,2,3,1) up to rotation...
# Actually: (1,1,2,3), (1,2,3,1), (2,3,1,1), (3,1,1,2) are rotations of each other.
# So just 12/4 = 3 necklaces.
# (1,2,2,2): 4/4 = 1 necklace (if all rotations distinct)
# Actually (1,2,2,2) rotations: (1,2,2,2), (2,2,2,1), (2,2,1,2), (2,1,2,2) — all distinct.
# So 1 necklace.

# Total cyclic patterns: 1 + 3 + 1 = 5 necklaces.
# With S-arc counts: 3, 2, 1, 0 S-arcs (gap=1 means S-arc).

# S-arcs = #{gaps of size 1}:
# (1,1,1,4): 3 S-arcs
# (1,1,2,3): 2 S-arcs
# (1,1,3,2): same as above
# (1,2,1,3): 2 S-arcs but different arrangement
# (1,2,2,2): 1 S-arc
# None with 0 S-arcs from above compositions...
# Wait: (2,2,2,1) has 1 S-arc. What about (2,2,3,0)? No, gaps must be >= 1.
# To have 0 S-arcs, all gaps >= 2. Four gaps >= 2 sum to >= 8 > 7. Impossible!

# So EVERY 7-cycle on 7 vertices with 4 S-vertices has at LEAST 1 S-internal arc.
# This means k=0 is IMPOSSIBLE when all 7 vertices are visited!

# Wait, but our data shows k=0 cycles. Let me recheck...
# Oh! k=0 means the cycle visits 4 S-vertices but no two are adjacent in the cycle.
# But we just showed this requires gaps summing to 7 with each >= 2,
# and 4 gaps >= 2 sum to >= 8. Contradiction!

print("THEOREM: In a Hamiltonian 7-cycle on 7 vertices,")
print("if 4 vertices are in S, then at least 1 S-S arc must be used.")
print("Proof: 4 gaps between consecutive S-vertices sum to 7.")
print("Each gap >= 1. If all gaps >= 2, sum >= 8 > 7. Contradiction.")
print("So at least one gap = 1, meaning an S-S arc.")
print()

# Verify against our data
print("Verification against data:")
k0_total = 0
for ca in cycle_analysis:
    k0_A = ca['A_by_internal'].get(0, 0)
    k0_B = ca['B_by_internal'].get(0, 0)
    k0_total += k0_A + k0_B

print(f"  Total k=0 cycles across all examples: {k0_total}")

if k0_total == 0:
    print("  CONFIRMED: No 7-cycle avoids all S-S arcs!")
    print("  This means EVERY Hamiltonian cycle passes through the Vitali atom.")
    print("  Therefore: c7_A and c7_B are completely DISJOINT sets of cycles!")
    print("  dc7 = c7_B - c7_A (no cancellation, pure difference)")
else:
    print(f"  UNEXPECTED: Found k=0 cycles! Investigating...")

# Final synthesis
print(f"\n{'='*60}")
print("SYNTHESIS: VITALI ATOM AS TOPOLOGICAL DEFECT")
print(f"{'='*60}")

print("""
The Vitali atom is a TOPOLOGICAL DEFECT in the tournament:

1. EVERY Hamiltonian cycle MUST traverse the atom (at least 1 S-S arc).
2. Reversing the atom creates a COMPLETELY DIFFERENT set of Hamiltonian cycles.
3. The sets c7(A) and c7(B) are DISJOINT — no cycle exists in both.
4. dc7 = |c7(B)| - |c7(A)| is the "topological charge" of the defect.

The "hidden higher-dimensional structure" is:
- The atom sits at the intersection of ALL Hamiltonian cycles.
- Its "charge" (dc7) depends on how the outside tournament couples
  to the atom through the S-X arcs AND the X-X arcs.
- The X-X arcs determine the "winding" of outside vertices around the atom.
- This is a purely TOPOLOGICAL invariant — not captured by sigma (Level 1.5).

The hierarchy:
  Level 1 (lambda): sees the atom exists, can't see its effect on c7
  Level 1.5 (sigma): sees the ±1 perturbation pattern, but can't resolve dc7
  Level 2 (A): sees the full coupling including X-X, determines dc7

The sigma perturbation (±1 on 4*(n-4) pairs) is the "shadow" of the atom
on the 2-point correlation function. The full dc7 requires the 3-point
(and higher) correlation — exactly the X-X arcs.
""")

# Quantitative: distribution of S-arc counts in cycles
print(f"Distribution of S-arc counts in Hamiltonian cycles:")
for k in sorted(all_A_arcs.keys()):
    a_count = all_A_arcs.get(k, 0)
    b_count = all_B_arcs.get(k, 0)
    print(f"  k={k}: in A: {a_count} (×7), in B: {b_count} (×7), diff: {b_count - a_count}")

print("\nDone.")
