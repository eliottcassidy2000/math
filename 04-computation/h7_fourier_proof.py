#!/usr/bin/env python3
"""
h7_fourier_proof.py — opus-2026-03-14-S71g

Attempt to prove H=7 impossible using Fourier analysis.

Key facts:
1. Ĥ[∅] = n!/2^{n-1} (mean H)
2. |Ĥ[S]| = (n-2k)!/2^{n-2} for |S|=2k (even Walsh degree)
3. Ĥ[S] = 0 for odd |S|
4. H(T) = Σ_S Ĥ[S]·χ_S(T) where χ_S(T) ∈ {±1}

For H(T) = 7:
  7 = n!/2^{n-1} + Σ_{S≠∅} Ĥ[S]·χ_S(T)

The sum has magnitude bounded by Σ|Ĥ[S]|.
If 7 < mean - Σ|Ĥ[S≠∅]|, then H=7 is impossible.
(Minimum possible H = mean - Σ|Ĥ[S≠∅]|)

But at small n, the minimum H is 1 (transitive tournament), which IS close to 7.
So the bound won't work directly for small n.

Instead: use the OCF connection. H = I(Ω, 2) where Ω is the conflict graph.
H=7 requires Ω=K₃ (unique). Show K₃ can't be Ω(T) for any tournament T.
"""

from math import factorial, comb
from itertools import combinations

print("=" * 60)
print("FOURIER BOUNDS ON H")
print("=" * 60)

for n in range(3, 12):
    mean_H = factorial(n) / 2**(n-1)

    # Number of nonzero Fourier coefficients at each even degree 2k
    # From THM-077: only S corresponding to unions of even-length paths contribute
    # The number of such S at degree 2k is related to the number of ways to
    # partition 2k arcs into even-length paths in K_n

    # Simpler bound: at degree 2k, there are at most C(m, 2k) coefficients
    # where m = C(n,2) = number of arcs
    m = comb(n, 2)

    total_fourier_energy = 0
    for k in range(1, n//2 + 1):
        deg = 2*k
        if n - deg < 0:
            break
        coeff_mag = factorial(n - deg) / 2**(n-2)
        # Upper bound on number of nonzero coefficients at this degree
        # Actually, the exact count matters. From THM-077:
        # nonzero at degree 2k iff S is a union of edge-disjoint even-length paths
        # For degree 2 (k=1): # nonzero = # adjacent arc pairs = n*(n-1)*(n-2)/2
        # (each pair of arcs sharing a vertex, minus some)
        # Actually at degree 2: each S is a 2-path (length-2 path in the line graph)
        # Number of 2-paths in K_n = n * C(n-1, 2) = n*(n-1)*(n-2)/2
        # But we use UNdirected arcs as variables, so...

        # Let's just compute the total Fourier energy bound
        num_at_deg = comb(m, deg)  # UPPER bound
        contribution = num_at_deg * coeff_mag
        total_fourier_energy += contribution

    min_H_bound = mean_H - total_fourier_energy
    max_H_bound = mean_H + total_fourier_energy

    print(f"  n={n}: mean H = {mean_H:.1f}")
    print(f"         Fourier bound: H ∈ [{min_H_bound:.1f}, {max_H_bound:.1f}]")
    print(f"         Actual min H = 1 (transitive)")
    print(f"         H=7 possible by bound? {min_H_bound <= 7 <= max_H_bound}")

# The bound is too loose because C(m, 2k) vastly overcounts nonzero coefficients.
# Let's use the exact count of nonzero coefficients.

print(f"\n{'='*60}")
print("EXACT NONZERO FOURIER COUNT (from THM-077)")
print(f"{'='*60}")

# From THM-077: Ĥ[S] ≠ 0 iff S is a union of even-length paths in K_n
# (considering arcs as undirected edges for the path structure)
#
# At degree 2 (|S|=2): S = {e1, e2} where e1, e2 share a vertex (2-path)
# Count: for each vertex v, choose 2 of its n-1 neighbors = C(n-1, 2)
# Total: n * C(n-1, 2) (but each 2-path counted once since it has unique center)
# Actually: n * C(n-1, 2) counts ordered? No.
# A 2-path a-v-b has center v and endpoints a,b.
# For each v: C(n-1, 2) such paths. Total: n * C(n-1, 2).
# For arcs: each 2-path involves 2 arcs. The arc pair {(a,v),(v,b)} or similar.
# In undirected terms: 2-paths in K_n. Count = n * C(n-1, 2).
# But in directed (Walsh) terms: each undirected 2-path a-v-b gives TWO arcs.
# The Walsh subset S = {arc1, arc2} where arc1 connects a-v and arc2 connects v-b.
# For the undirected pair {a,v} we can have x_{av} (a<v) as variable.
# So S = {x_{min(a,v),max(a,v)}, x_{min(v,b),max(v,b)}}.

# Number of undirected 2-paths in K_n = n * C(n-1, 2)
# At degree 4: unions of two disjoint 2-paths, or single 4-path
# This gets complicated. Let me just check computationally for small n.

for n in range(3, 8):
    m = comb(n, 2)

    # Degree 2: count 2-paths
    count_deg2 = n * comb(n-1, 2)

    # Degree 4: count of S with |S|=4 that are unions of even-length paths
    # Case A: two disjoint 2-paths
    # Case B: single 4-path (4 edges forming a path in K_n)
    # 4-paths in K_n: choose 5 vertices, order them as path = 5!/2 (divide by reversal)
    # Actually no: 4-path uses 5 vertices in a specific order.
    # Number of labeled paths of length 4 in K_n = n*(n-1)*(n-2)*(n-3)*(n-4)/2
    if n >= 5:
        count_4path = n * (n-1) * (n-2) * (n-3) * (n-4) // 2
    else:
        count_4path = 0

    # Two disjoint 2-paths: choose 2 non-overlapping 2-paths
    # Each 2-path uses 3 vertices. Two disjoint use 6 vertices (if no overlap).
    # But they can share vertices as long as edge-disjoint? No: "disjoint 2-paths"
    # means vertex-disjoint for the path decomposition.
    # Actually, in the Fourier formula, S needs to decompose into EDGE-disjoint
    # even-length paths. Vertex sharing is allowed.

    # For now, let's just compute the total Fourier deviation bound
    coeff_deg2 = factorial(n-2) / 2**(n-2)
    if n >= 5:
        coeff_deg4 = factorial(n-4) / 2**(n-2)
    else:
        coeff_deg4 = 0

    max_dev = count_deg2 * coeff_deg2
    if n >= 5:
        # Upper bound for degree 4 contribution
        count_deg4_upper = comb(m, 4)
        max_dev += count_deg4_upper * coeff_deg4
    if n >= 7:
        coeff_deg6 = factorial(n-6) / 2**(n-2)
        count_deg6_upper = comb(m, 6)
        max_dev += count_deg6_upper * coeff_deg6

    mean_H = factorial(n) / 2**(n-1)
    min_H = mean_H - max_dev

    print(f"  n={n}: deg-2 terms={count_deg2}, |coeff|={coeff_deg2:.4f}")
    print(f"         max deviation from mean ≥ {max_dev:.1f}")
    print(f"         mean={mean_H:.1f}, min_H_bound={min_H:.1f}")

# ============================================================
# The real proof: structural, not Fourier
# ============================================================
print(f"\n{'='*60}")
print("THE STRUCTURAL PROOF: Ω=K₃ IMPOSSIBLE FOR TOURNAMENTS")
print(f"{'='*60}")

print("""
THEOREM: No tournament T has Ω(T) = K₃.

PROOF:
  Ω(T) = K₃ requires exactly 3 odd cycles C₁, C₂, C₃ in T,
  with every pair sharing at least one vertex.

  CASE 1: All three are 3-cycles.
    Let C₁ = (a→b→c→a), C₂ = (d→e→f→d), C₃ = (g→h→i→g).
    Each pair shares a vertex, so the vertex sets have pairwise
    non-empty intersection.

    Sub-case 1a: All share a common vertex v.
      WLOG v=a=d=g. Then C₁=(v→b→c→v), C₂=(v→e→f→v), C₃=(v→h→i→v).
      If {b,c}∩{e,f}=∅ and {b,c}∩{h,i}=∅ and {e,f}∩{h,i}=∅:
        T has ≥7 vertices: v,b,c,e,f,h,i.
        The sub-tournament on {b,c,e,f,h,i} must have 0 directed 3-cycles
        (else we'd have >3 total). This means it's transitive.
        But in a 7-vertex tournament, the arcs between the 3 "cycle pairs"
        {b,c}, {e,f}, {h,i} create 5-cycles through v.
        Specifically: v→b→?→?→?→v using vertices from different cycles.
        Since the sub-tournament on 6 vertices is transitive,
        there exists a Hamiltonian path b₁→b₂→...→b₆ through {b,c,e,f,h,i}.
        The 5-cycle v→b₁→b₂→b₃→b₄→v exists iff b₄→v.
        Since v has 3 in-neighbors {c,f,i} among the 6 vertices,
        we need b₄ ∈ {c,f,i}. In the transitive order, this happens
        for specific configurations.

        *** NEED TO VERIFY: does every such configuration create
            at least one additional odd cycle (5-cycle or 7-cycle)? ***

    Sub-case 1b: Some cycle pairs share 2 vertices.
      If C₁ and C₂ share 2 vertices, say C₁=(a→b→c→a) and C₂=(a→b→d→a).
      This means a→b is an arc in both. For C₂: a→b→d→a, so b→d→a.
      For C₁: a→b→c→a, so b→c→a.
      Now c and d are both beat by b and both beat a.
      The arc between c and d creates either c→d or d→c.
      If c→d: then a→b→c→d→a is a 4-cycle (EVEN, not counted in Ω).
      If d→c: then a→b→d→c→a is a 4-cycle (EVEN, not counted).
      Either way, no new ODD cycle from {a,b,c,d} alone.
      But the sub-tournament on {a,b,c,d} has 4 vertices.
      It has t₃ ≥ 2 (C₁ and C₂) and possibly more.

      Actually: {a,b,c,d}. Arcs: a→b, b→c, c→a (from C₁), b→d, d→a (from C₂).
      Remaining: c vs d. If c→d: {a,b,c,d} has arcs a→b, b→c, b→d, c→a, c→d, d→a.
      3-cycles: (a→b→c→a)✓, (a→b→d→a)✓, (b→c→d→? : need d→b, but b→d)✗.
      So t₃=2 on these 4 vertices. c→d case: no new 3-cycle.
      If d→c: arcs a→b, b→c, b→d, c→a, d→a, d→c.
      3-cycles: (a→b→c→a)✓, (a→b→d→a)✓, (b→d→c→? : need c→b, but b→c)✗.
      So t₃=2 again. Good.

      Now C₃ must share a vertex with both C₁ and C₂.
      C₃ shares with C₁ ⟹ C₃ ∩ {a,b,c} ≠ ∅.
      C₃ shares with C₂ ⟹ C₃ ∩ {a,b,d} ≠ ∅.
      C₃ is a 3-cycle using (some of {a,b,c,d}) ∪ (new vertices).

      If C₃ entirely within {a,b,c,d}: we showed t₃=2 on these, so C₃
      must be one of C₁ or C₂. But C₃ is distinct. Contradiction.

      So C₃ uses at least one new vertex e.
      C₃ = (x→y→e→x) or similar, with x,y ∈ {a,b,c,d}.
      This creates a 5-vertex sub-tournament {a,b,c,d,e}.
      Need: no OTHER 3-cycle or 5-cycle is created.

      *** This requires very specific arc patterns. ***
      *** Computational verification needed. ***

  CASE 2: One or more cycles have length ≥ 5.
    A 5-cycle uses 5 vertices. Three 5-cycles all pairwise sharing vertices
    use at most 15 vertices. But each 5-cycle creates MANY sub-3-cycles
    (a 5-vertex tournament typically has 1-5 directed 3-cycles).
    So having a 5-cycle almost always means t₃ ≥ 1, adding to the count.

    Exception: if T restricted to the 5-cycle vertices has t₃=0
    (transitive on those 5), but then the 5-cycle doesn't exist
    in a transitive tournament. Contradiction.

    Actually: a directed 5-cycle on vertices {a,b,c,d,e} requires
    a→b→c→d→e→a. The sub-tournament on these 5 has this 5-cycle
    plus ALL other arcs. These other arcs create 3-cycles.
    For the 5-cycle a→b→c→d→e→a, the remaining arcs are:
    a-c, a-d, b-d, b-e, c-e. If ALL go "forward" (a→c, a→d, b→d, b→e, c→e),
    the tournament has score sequence (4,3,2,1,0) — TRANSITIVE! But then
    there's no 5-cycle. Contradiction.

    So a 5-cycle always coexists with at least one 3-cycle.
    Therefore if we have a 5-cycle as one of our 3 cycles,
    we automatically have a 3-cycle as another, but then the
    total cycle count is ≥ 3 only if we have exactly
    {one 5-cycle, two 3-cycles, no other cycles}. Unlikely.
""")

# Verify: does a directed 5-cycle always create 3-cycles?
print("Verification: 5-cycle subtournament 3-cycle count")
print("  5-cycle: 0→1→2→3→4→0")
print("  Remaining 5 arcs: (0,2),(0,3),(1,3),(1,4),(2,4)")
print("  Each can go either way: 2^5 = 32 completions")

n = 5
base_5cycle = [(0,1),(1,2),(2,3),(3,4),(4,0)]
remaining = [(0,2),(0,3),(1,3),(1,4),(2,4)]

for bits in range(32):
    A = [[0]*5 for _ in range(5)]
    for (i,j) in base_5cycle:
        A[i][j] = 1
    for idx, (i,j) in enumerate(remaining):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    # Count 3-cycles
    t3 = 0
    for combo in combinations(range(5), 3):
        for p in [(combo[0],combo[1],combo[2]),
                  (combo[0],combo[2],combo[1])]:
            if A[p[0]][p[1]] and A[p[1]][p[2]] and A[p[2]][p[0]]:
                t3 += 1
    t3 //= 1  # each 3-cycle found once with canonical orientation

    if t3 == 0:
        print(f"    bits={bits}: t₃=0 (NO 3-cycles alongside 5-cycle)")
        # This shouldn't happen — let's see
        print(f"      Arcs: {[(i,j) for i in range(5) for j in range(5) if A[i][j]]}")

min_t3 = float('inf')
max_t3 = 0
for bits in range(32):
    A = [[0]*5 for _ in range(5)]
    for (i,j) in base_5cycle:
        A[i][j] = 1
    for idx, (i,j) in enumerate(remaining):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1
    t3 = 0
    for a in range(5):
        for b in range(5):
            for c in range(5):
                if a!=b and b!=c and a!=c:
                    if A[a][b] and A[b][c] and A[c][a]:
                        t3 += 1
    t3 //= 3  # each cycle counted 3 times
    min_t3 = min(min_t3, t3)
    max_t3 = max(max_t3, t3)

print(f"\n  5-cycle completions: t₃ ∈ [{min_t3}, {max_t3}]")
print(f"  t₃ ≥ {min_t3} ALWAYS when a 5-cycle exists on 5 vertices")

if min_t3 >= 1:
    print(f"\n  PROVED: A directed 5-cycle always creates ≥{min_t3} directed 3-cycles!")
    print(f"  Therefore: if Ω has a 5-cycle vertex, it also has ≥{min_t3} 3-cycle vertices.")
    print(f"  This means α₁ ≥ 1 + {min_t3} = {1+min_t3} just from one 5-cycle.")

# Now check: 5-cycle + its forced 3-cycles, do the 3-cycles conflict with the 5-cycle?
print(f"\n  Checking conflict between 5-cycle and its forced 3-cycles:")
for bits in range(32):
    A = [[0]*5 for _ in range(5)]
    for (i,j) in base_5cycle:
        A[i][j] = 1
    for idx, (i,j) in enumerate(remaining):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    # Find all cycles
    five_cycle_verts = frozenset(range(5))  # the 5-cycle uses all 5 vertices
    three_cycles = []
    for a in range(5):
        for b in range(a+1, 5):
            for c in range(b+1, 5):
                if (A[a][b] and A[b][c] and A[c][a]) or \
                   (A[a][c] and A[c][b] and A[b][a]):
                    three_cycles.append(frozenset([a,b,c]))

    if three_cycles:
        # All 3-cycles share vertices with 5-cycle (they're subsets of same 5 vertices!)
        all_conflict = all(five_cycle_verts & tc for tc in three_cycles)
        # This is trivially true since all vertices are in {0,1,2,3,4}

# The key point: when we have a 5-cycle on vertices {0,1,2,3,4},
# any 3-cycle on a subset of these vertices trivially shares vertices with the 5-cycle.
# So they're automatically in conflict!
# The question is: do the 3-cycles conflict with EACH OTHER?
# Since they're on subsets of {0,1,2,3,4}, they may or may not share vertices.

print(f"\n  Key: 3-cycles forced by a 5-cycle are subsets of the SAME 5 vertices.")
print(f"  So they automatically share vertices with the 5-cycle (conflict).")
print(f"  And two 3-cycles on 5 vertices can share 0, 1, or 2 vertices.")
print(f"  Independent pair possible only if they use disjoint 3-subsets of {0,...,4}.")
print(f"  But C(5,3)=10 triples, and disjoint pairs of triples from 5 vertices")
print(f"  require 6 vertices — IMPOSSIBLE! So on 5 vertices, ALL pairs of 3-cycles")
print(f"  share at least one vertex!")

print(f"\n  THEREFORE: On 5 vertices, any 5-cycle + its forced 3-cycles gives")
print(f"  a conflict graph that is COMPLETE (all pairwise conflicting).")
print(f"  The independence number of this Ω restricted to 5 vertices is 1.")
print(f"  So α₂ = 0 is satisfied, but α₁ ≥ 1 + min_t3 = {1+min_t3} > 3")
print(f"  whenever min_t3 ≥ 3, which... let me check.")

if min_t3 >= 3:
    print(f"  YES! min_t3={min_t3} ≥ 3, so α₁ ≥ 4 > 3 and I(Ω,2) ≥ 1+4·2=9 > 7.")
else:
    print(f"  min_t3={min_t3} < 3, need more analysis.")
    print(f"  If min_t3=1 or 2, we could have 5-cycle + 1 or 2 three-cycles = 2 or 3 total.")
    print(f"  But we need EXACTLY 3 total for Ω=K₃.")

print(f"\n{'='*60}")
print("DONE")
print(f"{'='*60}")
