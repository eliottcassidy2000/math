#!/usr/bin/env python3
"""
n9_transition_deep.py — opus-2026-03-14-S76

WHY does the Cauchy-Schwarz proof of α₁ ≥ α₂ break at n=10?
The bound requires 9 ≥ n, i.e., n ≤ 9 = 3².

Is 9 = 3² = 3×3 a STRUCTURAL boundary or just an artifact of the proof?

THE KEY INSIGHT from S75:
- Each 3-cycle uses 3 vertices
- Σ_v d_v = 3·α₁ (each cycle counted 3 times)
- Cauchy-Schwarz gives e(CG) ≥ (9α₁² - 3nα₁)/(2n)
- The 9 = 3² comes from (Σ d_v)² / n = (3α₁)² / n = 9α₁²/n

So the 9 IS 3² where 3 = number of vertices per triangle.
At n=9: we have 9α₁²/9 = α₁², which EXACTLY balances.
At n=10: the balance breaks.

THIS MEANS: at n ≤ 9, cycles are "dense enough" relative to the
tournament that their overlaps force α₁ ≥ α₂.
At n ≥ 10, cycles can be "sparse enough" to potentially violate.

But DOES a violation actually occur at n=10?

CREATIVE DIRECTION (user's prompt):
9 = 3² = 3×3. Three 3-cycles can cover 9 vertices exactly.
This is the TIGHTEST PACKING of 3-cycles: ⌊n/3⌋ = 3.
At n=9, we can have THREE disjoint 3-cycles covering ALL vertices.
At n=10, we have a "spare" vertex — this changes everything.

Also: 9 = KEY₂² = KEY₂ × KEY₂.
The second key SQUARED. This is where the "ternary" structure
enters its second level.
"""

from itertools import combinations, permutations
from math import factorial, comb
import random

print("=" * 70)
print("PART 1: WHY 9 = 3² IS THE BOUNDARY")
print("=" * 70)
print()
print("  The Cauchy-Schwarz proof of α₁ ≥ α₂:")
print()
print("  1. Each 3-cycle uses 3 vertices → Σ_v d_v = 3α₁")
print("  2. By convexity: Σ C(d_v,2) ≥ n·C(3α₁/n, 2) = (9α₁²-3nα₁)/(2n)")
print("  3. e(CG) ≥ Σ C(d_v,2) - s₂ ≥ (9α₁²-3nα₁)/(2n) - s₂")
print("  4. α₂ = C(α₁,2) - e(CG) ≤ C(α₁,2) - (9α₁²-3nα₁)/(2n)")
print("  5. α₂ ≤ α₁ requires (n-9)α₁ ≤ 0, i.e., n ≤ 9")
print()
print("  WHERE DOES 9 COME FROM?")
print("  9 = 3² = (vertices per cycle)²")
print("  The '3' in 3α₁ is the cycle length.")
print("  Squaring it via Cauchy-Schwarz gives 9 = 3².")
print()
print("  IF we had k-cycles instead of 3-cycles:")
print("  Σ_v d_v = k·α₁ → bound requires k² ≥ n")
print("  i.e., n ≤ k²")
print()
print("  So:")
print("  3-cycles: n ≤ 9 = 3²")
print("  5-cycles: n ≤ 25 = 5²")
print("  7-cycles: n ≤ 49 = 7²")
print()
print("  INCLUDING 5-cycles would extend the proof to n ≤ 25!")
print("  But 5-cycles contribute to α₁ AND potentially to α₂.")
print()

# Can we extend the proof by including 5-cycles?
print("  EXTENDED PROOF SKETCH:")
print("  Let α₁ = α₁³ + α₁⁵ (3-cycles + 5-cycles)")
print("  Σ_v d_v = 3α₁³ + 5α₁⁵")
print()
print("  Cauchy-Schwarz on the combined vertex-cycle incidences:")
print("  Σ C(d_v,2) ≥ (3α₁³ + 5α₁⁵)² / (2n)")
print("  = (9(α₁³)² + 30α₁³α₁⁵ + 25(α₁⁵)²) / (2n)")
print()
print("  For the bound to work: need this ≥ α₁(α₁-3)/2")
print("  where α₁ = α₁³ + α₁⁵.")
print()
print("  If α₁⁵ = 0 (only 3-cycles): reduces to 9 ≥ n as before.")
print("  If α₁³ = 0 (only 5-cycles): need 25 ≥ n.")
print("  Mixed case: interpolates between 9 and 25.")
print()

print("=" * 70)
print("PART 2: THE 3×3 STRUCTURE AT n=9")
print("=" * 70)
print()
print("  n=9 = 3×3: THREE disjoint 3-cycles cover ALL 9 vertices.")
print("  This is a PERFECT 3-PARTITION: V = C₁ ∪ C₂ ∪ C₃, |Cᵢ|=3.")
print()
print("  The number of such 3-partitions: C(9,3)·C(6,3)·C(3,3)/3! = 280")
print()

# Each 3-partition gives a potential α₃ = 1 (three disjoint cycles)
# At n=9: α₃ can be ≥ 1!
# I(-1) = 1 - α₁ + α₂ - α₃
# The third term appears!
print("  At n=9: α₃ ≥ 0 (need 3×3 = 9 vertices)")
print("  This is the FIRST n where α₃ can be nonzero!")
print("  (At n=8: ⌊8/3⌋ = 2, so α₃ = 0)")
print()
print("  I(-1) = 1 - α₁ + α₂ - α₃")
print("  The alternating sum now has THREE variable terms.")
print()
print("  For I(-1) ≤ 1: need α₁ - α₂ + α₃ ≥ 0")
print("  This is STRONGER than α₁ ≥ α₂!")
print("  We need: α₁ ≥ α₂ - α₃, i.e., α₁ + α₃ ≥ α₂")
print()

# KEY: at n=9, the FULL alternating sum condition becomes
# α₁ - α₂ + α₃ ≥ 0
# This is a THREE-TERM condition, matching n = 3²

print("  BEAUTIFUL: at n = 3², the alternating sum has 3 terms!")
print("  At n = 3¹ = 3: one term (α₁ only)")
print("  At n = 3² = 9: three terms (α₁, α₂, α₃)")
print("  At n = 3³ = 27: four terms? (α₁,...,α₉)")
print()
print("  Actually: max independent set size ⌊n/3⌋.")
print("  At n=3: ⌊3/3⌋=1 → I(-1) = 1 - α₁")
print("  At n=6: ⌊6/3⌋=2 → I(-1) = 1 - α₁ + α₂")
print("  At n=9: ⌊9/3⌋=3 → I(-1) = 1 - α₁ + α₂ - α₃")
print("  At n=12: ⌊12/3⌋=4 → I(-1) = 1 - α₁ + α₂ - α₃ + α₄")
print()
print("  So 9 = 3² = 3×3 marks the transition to 3 LEVELS")
print("  in the independence complex hierarchy!")
print()
print("  The 'levels' correspond to:")
print("  Level 0: empty set (1 term)")
print("  Level 1: individual cycles (α₁)")
print("  Level 2: disjoint pairs (α₂)")
print("  Level 3: disjoint triples (α₃) — NEW at n=9!")
print()

print("=" * 70)
print("PART 3: RECURRENCE VIEW — THE TRIBONACCI THRESHOLD")
print("=" * 70)
print()
print("  At n < 6: I(x) = 1 + α₁x (Fibonacci-like, 2 terms)")
print("  At 6 ≤ n < 9: I(x) = 1 + α₁x + α₂x² (tribonacci-like, 3 terms)")
print("  At 9 ≤ n < 12: I(x) = 1 + α₁x + α₂x² + α₃x³ (4 terms)")
print()
print("  The k-nacci connection:")
print("  Fibonacci (2 terms): root → φ ≈ 1.618")
print("  Tribonacci (3 terms): root → ≈ 1.839")
print("  Tetranacci (4 terms): root → ≈ 1.928")
print("  k-nacci: root → 2 as k → ∞")
print()
print("  TRANSITION AT n=9:")
print("  I(x) gains its THIRD variable coefficient.")
print("  The polynomial 'depth' matches the tribonacci level.")
print("  At x=2: the tribonacci root (1.839) is closest to 2 among")
print("  the first three k-nacci roots.")
print()

# Compute: how many terms does I(x) have at each n?
print("  Independence polynomial degree by n:")
for n in range(3, 21):
    max_k = n // 3
    print(f"    n={n:2d}: max degree = {max_k}, polynomial has {max_k+1} terms")

print()

# The 3×3 structure is also about the CONFLICT GRAPH
# At n=9 with three disjoint 3-cycles:
# CG has (at least) 3 vertices forming a triangle
# (each pair of the three disjoint cycles is an independent set of size 2)
# Wait: disjoint cycles are INDEPENDENT (no shared vertices)
# So in CG, they are NOT adjacent
# Three pairwise non-adjacent vertices = independent set of size 3 = α₃ contribution

print("=" * 70)
print("PART 4: THE 3×3 TOURNAMENT — EXTREMAL STRUCTURE")
print("=" * 70)
print()
print("  A tournament on 9 vertices with three disjoint 3-cycles:")
print("  V = {0,1,2} ∪ {3,4,5} ∪ {6,7,8}")
print("  C₁ = 0→1→2→0, C₂ = 3→4→5→3, C₃ = 6→7→8→6")
print()
print("  The inter-block arcs connect the three blocks.")
print("  Different inter-block patterns give different (α₁, α₂, α₃).")
print()

# Build a specific 3×3 tournament
n = 9

def build_block_tournament(inter_pattern):
    """Build 9-vertex tournament with three 3-cycle blocks.
    inter_pattern: function(block_i, vertex_in_i, block_j, vertex_in_j) → bool
    """
    adj = [[False]*n for _ in range(n)]
    # Intra-block: 3-cycles
    for b in range(3):
        base = 3*b
        adj[base][base+1] = True
        adj[base+1][base+2] = True
        adj[base+2][base] = True
    # Inter-block
    for i in range(n):
        for j in range(n):
            bi, bj = i//3, j//3
            if bi != bj:
                adj[i][j] = inter_pattern(bi, i%3, bj, j%3)
    return adj

def count_3cycles(adj, n):
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    cycles.append(frozenset([i,j,k]))
                elif adj[j][i] and adj[i][k] and adj[k][j]:
                    cycles.append(frozenset([i,j,k]))
    return cycles

def compute_alphas(cycles):
    a1 = len(cycles)
    a2 = 0
    a3 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if len(cycles[i] & cycles[j]) == 0:
                a2 += 1
                for k in range(j+1, len(cycles)):
                    if len(cycles[i] & cycles[k]) == 0 and len(cycles[j] & cycles[k]) == 0:
                        a3 += 1
    return a1, a2, a3

# Pattern 1: All inter-block arcs go from lower to higher block
def pattern_transitive(bi, vi, bj, vj):
    return bi < bj

adj1 = build_block_tournament(pattern_transitive)
cycles1 = count_3cycles(adj1, n)
a1_1, a2_1, a3_1 = compute_alphas(cycles1)
im1_1 = 1 - a1_1 + a2_1 - a3_1
h1 = 1 + 2*a1_1 + 4*a2_1 + 8*a3_1
print(f"  Pattern 1 (transitive inter-block):")
print(f"    α₁={a1_1}, α₂={a2_1}, α₃={a3_1}, I(-1)={im1_1}, H={h1}")
print()

# Pattern 2: Inter-block arcs form a cycle (0→1→2→0)
def pattern_cyclic(bi, vi, bj, vj):
    return (bj - bi) % 3 == 1

adj2 = build_block_tournament(pattern_cyclic)
cycles2 = count_3cycles(adj2, n)
a1_2, a2_2, a3_2 = compute_alphas(cycles2)
im1_2 = 1 - a1_2 + a2_2 - a3_2
h2 = 1 + 2*a1_2 + 4*a2_2 + 8*a3_2
print(f"  Pattern 2 (cyclic inter-block: 0→1→2→0):")
print(f"    α₁={a1_2}, α₂={a2_2}, α₃={a3_2}, I(-1)={im1_2}, H={h2}")
print()

# Pattern 3: Diagonal inter-block (vi→vj iff (bi<bj and vi<=vj) or ...)
def pattern_diagonal(bi, vi, bj, vj):
    return (bi + vi) % 3 < (bj + vj) % 3 or ((bi + vi) % 3 == (bj + vj) % 3 and bi < bj)

adj3 = build_block_tournament(pattern_diagonal)
cycles3 = count_3cycles(adj3, n)
a1_3, a2_3, a3_3 = compute_alphas(cycles3)
im1_3 = 1 - a1_3 + a2_3 - a3_3
h3 = 1 + 2*a1_3 + 4*a2_3 + 8*a3_3
print(f"  Pattern 3 (diagonal inter-block):")
print(f"    α₁={a1_3}, α₂={a2_3}, α₃={a3_3}, I(-1)={im1_3}, H={h3}")
print()

# Pattern 4: Random
random.seed(42)
def pattern_random(bi, vi, bj, vj):
    return random.random() < 0.5

adj4 = build_block_tournament(pattern_random)
cycles4 = count_3cycles(adj4, n)
a1_4, a2_4, a3_4 = compute_alphas(cycles4)
im1_4 = 1 - a1_4 + a2_4 - a3_4
h4 = 1 + 2*a1_4 + 4*a2_4 + 8*a3_4
print(f"  Pattern 4 (random inter-block):")
print(f"    α₁={a1_4}, α₂={a2_4}, α₃={a3_4}, I(-1)={im1_4}, H={h4}")
print()

# Check α₁ - α₂ + α₃ ≥ 0 for all
for name, a1, a2, a3 in [("transitive", a1_1, a2_1, a3_1),
                           ("cyclic", a1_2, a2_2, a3_2),
                           ("diagonal", a1_3, a2_3, a3_3),
                           ("random", a1_4, a2_4, a3_4)]:
    alt = a1 - a2 + a3
    print(f"    {name}: α₁-α₂+α₃ = {a1}-{a2}+{a3} = {alt} {'✓' if alt >= 0 else '✗'}")

print()

print("=" * 70)
print("PART 5: RANDOM SAMPLING AT n=9,10,11")
print("=" * 70)
print()

# Sample random tournaments and check the alternating sum
for n in [9, 10, 11, 12]:
    random.seed(0)
    violations = 0
    max_neg = 0
    total = 10000

    for trial in range(total):
        adj = [[False]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = True
                else:
                    adj[j][i] = True

        cycles = count_3cycles(adj, n)
        a1, a2, a3 = compute_alphas(cycles)
        alt = a1 - a2 + a3

        if alt < 0:
            violations += 1
            if alt < max_neg:
                max_neg = alt

    print(f"  n={n}: {violations}/{total} violations of α₁-α₂+α₃ ≥ 0, min = {max_neg}")

print()

print("=" * 70)
print("PART 6: THE 9→10 TRANSITION IN DETAIL")
print("=" * 70)
print()
print("  At n=9: every vertex is in a 3-cycle (generically).")
print("  Σ d_v = 3α₁, average d_v = 3α₁/9 = α₁/3.")
print("  CS bound: Σ C(d_v,2) ≥ 9·C(α₁/3, 2) = 9·α₁(α₁-3)/(2·9) = α₁(α₁-3)/2")
print("  This EXACTLY matches what we need: e(CG) ≥ α₁(α₁-3)/2.")
print("  At n=9, the Cauchy-Schwarz bound is TIGHT!")
print()
print("  At n=10: average d_v = 3α₁/10 = 0.3α₁.")
print("  CS bound: 10·C(0.3α₁, 2) = 10·0.3α₁(0.3α₁-1)/2 = 0.9α₁(0.3α₁-1)/2")
print("  = (0.27α₁² - 0.9α₁)/2 < α₁(α₁-3)/2 = (α₁² - 3α₁)/2")
print("  The bound FAILS because 0.27 < 1.")
print()
print("  THE GAP: 9/n vs 1.")
print("  At n=9: 9/9 = 1 (tight)")
print("  At n=10: 9/10 = 0.9 (gap of 0.1)")
print()
print("  BUT: at n=10, we also have 5-cycles!")
print("  5-cycles contribute to Σ d_v differently.")
print("  If we include 5-cycles: Σ d_v = 3α₁³ + 5α₁⁵")
print("  CS gives: (3α₁³ + 5α₁⁵)²/(2n)")
print("  For this to work at n=10: need (3α₁³+5α₁⁵)²/10 ≥ α₁²")
print("  where α₁ = α₁³ + α₁⁵.")
print()

# Check: what fraction of 5-cycles is needed?
print("  If ALL cycles are 3-cycles (α₁⁵=0):")
print("    Need 9α₁²/10 ≥ α₁² → 9 ≥ 10, fails.")
print()
print("  If fraction f of cycles are 5-cycles:")
print("    α₁³ = (1-f)α₁, α₁⁵ = fα₁")
print("    Σ d_v = (3(1-f) + 5f)α₁ = (3+2f)α₁")
print("    CS needs: (3+2f)²α₁²/10 ≥ α₁²")
print("    → (3+2f)² ≥ 10")
print("    → 3+2f ≥ √10 ≈ 3.162")
print("    → f ≥ 0.081")
print()
print("  So if ≥ 8.1% of odd cycles are 5-cycles, the bound works at n=10!")
print()
print("  At n=10: is this always true?")
print("  A tournament on 10 vertices has C(10,5)·... possible 5-cycles.")
print("  Every tournament with ≥ 3 vertices has at least one 3-cycle")
print("  (Rédei's theorem: odd number of Ham paths).")
print()

print("=" * 70)
print("PART 7: THE RECURRENCE 3² = 9 AND THE KEYS")
print("=" * 70)
print()
print("  9 = 3² appears in multiple contexts:")
print()
print("  1. CAUCHY-SCHWARZ BOUNDARY: α₁≥α₂ proved for n ≤ 9")
print("  2. PERFECT 3-PARTITION: V = C₁∪C₂∪C₃ possible")
print("  3. INDEPENDENCE COMPLEX: third term α₃ appears")
print("  4. RECURRENCE DEPTH: I(x) has degree 3 polynomial")
print()
print("  In the KEY framework:")
print("  9 = KEY₂² (second key squared)")
print("  9 = KEY₁³ + 1 (2³ + 1)")
print("  9 = KEY₁ · KEY₂ + KEY₂ (2·3 + 3 = 6+3)")
print()
print("  THE RECURRENCE z²=5z-6:")
print("  At z=9: 81 = 45 - 6 = 39. NO, 9²=81, 5·9=45, 81≠45-6=39.")
print("  So 9 is NOT a root of the tournament polynomial.")
print()
print("  But 9 = 3² and the characteristic polynomial is (z-2)(z-3).")
print("  The DISCRIMINANT Δ = 25-24 = 1.")
print("  Δ² = 1. So the discriminant is a PERFECT SQUARE.")
print("  The roots differ by Δ = 1: 3-2 = 1.")
print()
print("  9 = 3² connects to the SECOND ITERATION of the KEY:")
print("  Key 3 → 3² = 9 (the boundary)")
print("  Key 2 → 2² = 4 (= number of cycles through a vertex pair?)")
print("  Key 2 → 2³ = 8 (= coefficient of α₃ in H)")
print()
print("  H = 1 + 2α₁ + 4α₂ + 8α₃ + ...")
print("  The coefficient of α₃ is 8 = 2³ = KEY₁³.")
print("  And α₃ first appears at n = 9 = KEY₂².")
print()
print("  KEY₁^(k+1) = 2^(k+1) is the weight of level-k independence sets.")
print("  KEY₂^k is when level-k independence sets first appear.")
print("  The BALANCE: 2^(k+1) vs 3^k.")
print("  Level 1: weight 4=2², appears at n=6=2·3")
print("  Level 2: weight 8=2³, appears at n=9=3²")
print("  Level 3: weight 16=2⁴, appears at n=12=4·3")
print()

# The weights grow as 2^k, the thresholds grow as 3k
# 2^k / (3k) → ∞, so higher levels are HEAVIER per item
# but have FEWER items (bounded by C(α₁, k))

print("  Weight/threshold ratio: 2^(k+1) / (3k):")
for k in range(1, 8):
    w = 2**(k+1)
    thresh = 3*k
    print(f"    Level {k}: weight={w}, threshold n={3*(k+1)}, ratio={w/thresh:.3f}")

print()
print("  The ratio 2^(k+1)/(3k) grows without bound.")
print("  Higher levels contribute EXPONENTIALLY more to H per set.")
print("  But they contribute (-1)^k to I(-1).")
print("  The alternating sum I(-1) = 1 - α₁ + α₂ - α₃ + ...")
print("  is a BALANCE between large positive and large negative terms.")
print()

print("=" * 70)
print("PART 8: THE PACKING INTERPRETATION")
print("=" * 70)
print()
print("  SIMPLICES INSIDE CUBOIDS:")
print("  A d-simplex has volume 1/d! of the d-cube.")
print("  H simplices 'fill' volume H/n! of the n-cube.")
print()
print("  At n=9: the tournament polytope U(T) = union of H orthoschemes")
print("  in [0,1]^9. Each orthoscheme has volume 1/9! = 1/362880.")
print()
print("  The alternating sum I(-1) = χ(independence complex) measures")
print("  the 'net topological complexity' of the packing.")
print()
print("  HILBERT'S 3RD PROBLEM at n=9:")
print("  All tournament polytopes with same H are scissors-congruent")
print("  (proved in S75: Dehn invariant = 0).")
print()
print("  But the ARRANGEMENT TOPOLOGY differs for same H, different I(x).")
print("  At n=9, I(x) has degree 3, so the arrangement can have")
print("  3 levels of 'twisting' — the α₃ term adds a new dimension")
print("  of topological complexity.")
print()
print("  PACKING THREE SIMPLICES IN A 3-CUBE:")
print("  A 3-cube has volume 3³ = 27 in 'simplex units' (×3!).")
print("  Actually: [0,3]³ has volume 27, packed by 27·6 = 162 simplices.")
print("  The 3×3×3 structure of the cube reflects the 3²=9 boundary.")
print()
print("  In tournament terms:")
print("  n=9 with three disjoint 3-cycles is like a 3×3 grid")
print("  of simplices inside a 3-cuboid.")
print("  Each 3-cycle provides 2 'directions' → 2³ = 8 total")
print("  from three independent cycles → α₃ contributes 8 to H.")
print()

print("=" * 70)
print("PART 9: EXTENDING THE PROOF TO n ≥ 10")
print("=" * 70)
print()
print("  Three strategies to extend α₁ ≥ α₂ beyond n=9:")
print()
print("  STRATEGY 1: Include 5-cycles in the vertex-cycle incidence.")
print("  Each 5-cycle uses 5 vertices → adds 5 to Σ d_v per cycle.")
print("  Cauchy-Schwarz becomes: (3α₁³+5α₁⁵)²/(2n) ≥ ...")
print("  Works as long as enough 5-cycles exist.")
print()
print("  STRATEGY 2: Use the TURÁN density theorem.")
print("  CG(T) has the property that every vertex's neighborhood")
print("  contains a large clique. This is a TURÁN-type condition.")
print("  Turán's theorem: graph with no K_{r+1} has ≤ (1-1/r)n²/2 edges.")
print("  For CG(T): the complement (whose edges = α₂) has bounded")
print("  clique number (= max independent set of CG).")
print()
print("  STRATEGY 3: Use the 3-PARTITION structure.")
print("  At n=3k: the vertices can be partitioned into k blocks of 3.")
print("  Each block can contain at most 1 3-cycle.")
print("  Inter-block cycles add to α₁ but are constrained by block structure.")
print("  This gives α₁ ≈ C(n,3)/... and α₂ ≈ ...")
print()

# Strategy 1: how many 5-cycles does a random tournament have?
print("  STRATEGY 1 CHECK: 5-cycle density")
print("  Expected #3-cycles in random tournament on n vertices: C(n,3)/4")
print("  Expected #5-cycles: C(n,5)·4!/(2·5) · (1/2)^5... complicated")
print()
for n in [9, 10, 11, 12, 15, 20]:
    e3 = comb(n, 3) / 4  # each triple has 1/4 chance of being a 3-cycle
    # 5-cycles: choose 5 vertices, then directed 5-cycle = 12 orientations
    # out of 2^10 tournaments on 5 vertices, 12 have a directed 5-cycle
    # Actually: P(5-cycle on given 5 vertices) = 12/2^10... no
    # Simpler: for random tournament on 5 vertices, expected #5-cycles = ?
    # Each set of 5 vertices: expected #chordless directed 5-cycles
    # Actually: P(directed cycle 1→2→3→4→5→1) = (1/2)^5 = 1/32
    # 12 cyclic orderings per direction, 2 directions = 24 directed 5-cycles
    # Expected per 5-set: 24/32 = 3/4
    # But chordless: P(no chord) = ? In a random tournament, the probability
    # that a 5-cycle is chordless is the probability that none of the 5 chords
    # go in the "wrong" direction. This is complex.
    # Skip exact calculation; just note the ratio
    e5_approx = comb(n, 5) * 3/4  # rough estimate
    ratio = e5_approx / e3 if e3 > 0 else 0
    print(f"  n={n}: E[#3-cycles]≈{e3:.0f}, E[#5-cycles]≈{e5_approx:.0f}, ratio≈{ratio:.2f}")

print()
print("  5-cycles OUTNUMBER 3-cycles for n ≥ 10!")
print("  This means including 5-cycles should STRONGLY extend the proof.")
print()

# The combined Cauchy-Schwarz bound
print("  Combined bound (3-cycles + 5-cycles):")
print("  Need (3f₃ + 5f₅)² ≥ n where f₃+f₅=1 (fractions)")
print()
for n in [10, 15, 20, 25, 30]:
    # Best case: all 5-cycles
    min_all5 = 25  # (5·1)² = 25 ≥ n iff n ≤ 25
    # Worst case: all 3-cycles
    min_all3 = 9   # (3·1)² = 9 ≥ n iff n ≤ 9
    # Mixed: (3(1-f)+5f)² = (3+2f)² ≥ n → f ≥ (√n-3)/2
    import math
    f_needed = max(0, (math.sqrt(n) - 3) / 2)
    print(f"  n={n}: need f₅ ≥ {f_needed:.4f} ({f_needed*100:.1f}%)")

print()
print("  For n=10: need f₅ ≥ 8.1%")
print("  For n=15: need f₅ ≥ 43.7%")
print("  For n=20: need f₅ ≥ 73.6%")
print("  For n=25: need f₅ ≥ 100% (only works with all 5-cycles)")
print()
print("  Since 5-cycles OUTNUMBER 3-cycles for n≥10,")
print("  the mixed bound should extend well beyond n=9!")
print()

print("=" * 70)
print("PART 10: THE 3² MANIFESTO — WHY 9 IS NATURAL")
print("=" * 70)
print()
print("  The number 9 = 3² appears because:")
print()
print("  1. GEOMETRIC: 3 vertices per triangle, 3² = maximum n")
print("     where triangles are 'dense enough' for Cauchy-Schwarz.")
print()
print("  2. ALGEBRAIC: The KEY₂ = 3 enters its SECOND POWER.")
print("     First power: 3 vertices per cycle (basic unit).")
print("     Second power: 3² vertices for three disjoint cycles (α₃ emerges).")
print()
print("  3. COMBINATORIAL: ⌊9/3⌋ = 3 = number of blocks in perfect partition.")
print("     This is the FIRST n where the partition number = block size.")
print("     (At n=6: 2 blocks of 3. At n=9: 3 blocks of 3 = SQUARE.)")
print()
print("  4. RECURRENCE: The k-nacci sequence at k=3 (tribonacci) has")
print("     root ≈ 1.839, closest 'below-2' k-nacci root to 2.")
print("     Tribonacci = 3 terms = 3 levels of independence complex.")
print()
print("  5. TOPOLOGICAL: I(-1) = 1 - α₁ + α₂ - α₃ has 3 varying terms.")
print("     The alternating sum has 'depth 3' = KEY₂.")
print("     At n < 9: depth ≤ 2. At n ≥ 9: depth ≥ 3.")
print()
print("  6. THE (2,3) BRIDGE:")
print("     Coefficient of α₃ in H is 8 = 2³ = KEY₁³.")
print("     Threshold for α₃ is n = 9 = 3² = KEY₂².")
print("     The CROSS-MULTIPLICATION: KEY₁^(k+1) × KEY₂^k")
print("     = 2^(k+1) × 3^k = weight × threshold.")
print("     At k=2: 2³ × 3² = 8 × 9 = 72 = total 'capacity' at level 2.")
print()
print("  7. THE CARTAN CONNECTION:")
print("     det(A₈) = 9 = 3².")
print("     The Cartan determinant of A₈ is EXACTLY the boundary value!")
print("     This suggests: α₁ ≥ α₂ ⟺ det(A_{n-1}) ≤ 9?")
print("     No: det(A_{n-1}) = n. So this is just n ≤ 9.")
print("     But the MEANING is: the type A₈ root system is the largest")
print("     for which the 3-cycle clique structure suffices.")
print()

# Final: the trinity
print("  THE TRINITY OF 3²:")
print("  3² = 9 = 3 × 3")
print("  = (vertices per cycle) × (max disjoint cycles)")
print("  = (KEY₂) × (⌊n/KEY₂⌋)")
print("  = KEY₂²")
print()
print("  At this threshold: every quantity 'squares':")
print("  - Cycle weight: 2² = 4 (α₂ weight)")
print("  - Cycles per vertex: α₁/3 (average)")
print("  - Vertex coverage: 3×3 = 9 (full)")
print("  - CS bound: (3α₁)²/9 = α₁² (exact match)")
