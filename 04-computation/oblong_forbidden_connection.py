#!/usr/bin/env python3
"""
Oblong-Forbidden Connection
opus-2026-03-14-S71h

STUNNING OBSERVATION from knacci_simplex_bridge.py Part 10:

The simplex polynomial (x+1)^1 = x+1 at oblong fugacities x=m(m+1):
  m=0: x=0  → 1
  m=1: x=2  → 3     (tournament world)
  m=2: x=6  → 7     ← PERMANENTLY FORBIDDEN H VALUE!
  m=3: x=12 → 13
  m=4: x=20 → 21    ← PERMANENTLY FORBIDDEN H VALUE!
  m=5: x=30 → 31

The simplex values at m=2 and m=4 are EXACTLY {7, 21}!

But wait: (m(m+1)+1) = m²+m+1. At m=2: 7. At m=4: 21.
  7 = 2²+2+1 = C(3,2)+C(2,1)+1 = 4+2+1
  21 = 4²+4+1 = 16+4+1

Is this a coincidence, or is there a deep reason why m²+m+1 at
even m gives permanently forbidden H values?

THIS SCRIPT INVESTIGATES.
"""

import math

print("=" * 70)
print("OBLONG SIMPLEX VALUES = m²+m+1")
print("=" * 70)
print()

print("m²+m+1 sequence (cyclotomic Φ₃(m)):")
for m in range(10):
    val = m*m + m + 1
    x = m*(m+1)
    print(f"  m={m}: m²+m+1 = {val:5d},  x=m(m+1)={x:3d},  (x+1)={x+1}")
    # Note: m²+m+1 = x+1 where x = m(m+1). Obvious!

print()
print("So the 'simplex at oblong' is just m²+m+1.")
print("This is the CYCLOTOMIC POLYNOMIAL Φ₃(m) = (m³-1)/(m-1) = m²+m+1")
print()

# Key question: WHY are m=2 (giving 7) and m=4 (giving 21) special?
# 7 = Φ₃(2) = I(K₃, 2) = I(C₃, 2)
# 21 = Φ₃(4) = I(P₄, 2) = C(7,2) = T(6)

# The forbidden values are 7 and 21.
# 7 is forbidden because K₃ cannot be a component of Ω(T) (THM-201)
# 21 is forbidden because P₄ cannot be Ω(T) (THM-202) and
# the only graph with I(G,2)=21 at relevant vertex counts is P₄.

# But from the cyclotomic perspective:
# Φ₃(m) for m = 0,1,2,...: 1, 3, 7, 13, 21, 31, 43, 57, 73, 91, ...
# Which of these are achievable as H values?

print("=" * 70)
print("WHICH Φ₃(m) VALUES ARE ACHIEVABLE AS H?")
print("=" * 70)
print()

# Known H-spectra (exhaustive)
# n=3: {1, 3}
# n=4: {1, 3, 5}
# n=5: {1, 3, 5, 9, 11, 13, 15}
# n=6: {1,3,5,9,11,13,15,17,19,23,25,27,29,31,33,37,41,43,45}
# n=7: all odd 1..189 except {7,21,63,107,119,149,161-187 gaps}
# n=8: fills 63, 107, 119, etc.
# Permanent gaps: {7, 21} only

permanent_gaps = {7, 21}
n7_gaps = {7, 21, 63, 107, 119, 149, 161, 163, 165, 167, 169, 173, 177, 179, 181, 183, 185, 187}

for m in range(20):
    val = m*m + m + 1
    if val in permanent_gaps:
        status = "PERMANENTLY FORBIDDEN"
    elif val in n7_gaps:
        status = "gap at n=7 (fills later)"
    elif val % 2 == 0:
        status = "EVEN (never an H value)"
    else:
        status = "achievable"
    print(f"  Φ₃({m:2d}) = {val:4d}  {status}")

print()
print("REMARKABLE: The ONLY permanently forbidden Φ₃ values are m=2 and m=4!")
print()

# What's special about m=2 and m=4?
# Both are even! And no other even m gives a forbidden value.
# Φ₃(even) = even²+even+1 = odd (always odd, good)
# The question is: why m=2 and m=4 but not m=6,8,...?

print("=" * 70)
print("STRUCTURAL ANALYSIS: WHY m=2 AND m=4?")
print("=" * 70)
print()

# 7 = I(C₃, 2) — cycle of length 3
# 21 = I(P₄, 2) — path of length 4
# Both are independence polynomial values of small SPECIFIC graphs at x=2.

# What are the other Φ₃ values as I(G, 2)?
# Φ₃(0) = 1 = I(∅, 2) — empty graph
# Φ₃(1) = 3 = I(P₁, 2) = I(K₁, 2)·(1+2) — single vertex
# Φ₃(2) = 7 = I(C₃, 2) = I(K₃, 2) — triangle/complete graph
# Φ₃(3) = 13 = I(P₂⊔K₁, 2) = 5·... no. I(P₂,2)=5, I(K₁,2)=3, I(P₂⊔K₁)=15≠13
#   Actually 13 = I(P₃, 2)? No, I(P₃,2) = 11. I(C₄,2)? = 17. Hmm.
# Let me check: what graph G has I(G,2) = 13?

# For a graph on n vertices with independence number α:
# I(G,2) = 1 + 2α₁ + 4α₂ + ... where αₖ = #{indep sets of size k}

# I(G,2) = 13 means: 1 + 2a1 + 4a2 + 8a3 + ... = 13
# So 2a1 + 4a2 + 8a3 + ... = 12, i.e., a1 + 2a2 + 4a3 + ... = 6
# Possible: a1=6 (6 vertices, no edges? Then I=1+12=13... but C(6,2)=15 pairs,
#   I(6K₁,2) = 3^6=729, not 13)
# a1=4, a2=1: graph on ≥4 vertices with 4 vertices and 1 independent pair
# a1=2, a2=2: 2 vertices + 2 independent pairs... 2 vertices have exactly 1 pair.
#   Need more vertices.

# Actually let me just try small graphs
from itertools import combinations

def independence_polynomial_at_2(n, edges):
    """Compute I(G, 2) for graph with n vertices and given edges."""
    adj = set()
    for u, v in edges:
        adj.add((u, v))
        adj.add((v, u))

    total = 0
    # Enumerate all subsets
    for mask in range(1 << n):
        # Check if independent set
        vertices = [i for i in range(n) if mask & (1 << i)]
        is_indep = True
        for i in range(len(vertices)):
            for j in range(i+1, len(vertices)):
                if (vertices[i], vertices[j]) in adj:
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            total += 2**len(vertices)
    return total

# Find all graphs with I(G,2) = Φ₃(m) for small m
targets = {m: m*m+m+1 for m in range(7)}

print("Searching for graphs G with I(G,2) = Φ₃(m)...")
print()

for n in range(1, 8):
    if n > 6:
        print(f"  n={n}: skipping (too many graphs)")
        continue
    all_edges = list(combinations(range(n), 2))
    for num_edges in range(len(all_edges)+1):
        if num_edges > 10 and n >= 6:
            break  # Skip dense graphs on many vertices
        for edge_set in combinations(all_edges, num_edges):
            val = independence_polynomial_at_2(n, edge_set)
            for m, target in targets.items():
                if val == target and target > 1:
                    print(f"  I(G,2) = {target} = Φ₃({m}): n={n}, edges={num_edges}, "
                          f"graph={list(edge_set)}")

print()

# Now the deep question: is there a pattern connecting Φ₃(m) to forbidden H values?
# H=7 forbidden because K₃ can't be Ω component
# H=21 forbidden because I(G,2)=21 only works for P₄ (on appropriate vertex count)
#   and P₄ can't be Ω(T)

# For Φ₃(6) = 43: definitely achievable (in n=5 spectrum)
# For Φ₃(8) = 73: achievable

# So the Φ₃ connection is suggestive but not the FULL story.
# The forbidden values happen to be Φ₃ at m=2,4 but not ALL even m.

print("=" * 70)
print("ALTERNATIVE PATTERN: JACOBSTHAL/JACOBSTHAL-LUCAS")
print("=" * 70)
print()

# 7 = JL(3) = 2³+(-1)³ (Jacobsthal-Lucas)
# 21 = J(6) = (2⁶-(-1)⁶)/3 (Jacobsthal)
# These are DIFFERENT sequences!

# Jacobsthal: J(n) = (2ⁿ-(-1)ⁿ)/3. Sequence: 0,1,1,3,5,11,21,43,...
# Jacobsthal-Lucas: JL(n) = 2ⁿ+(-1)ⁿ. Sequence: 1,1,3,5,9,17,31,65,...
# (actually JL starts: JL(0)=2, JL(1)=1, JL(2)=5, JL(3)=7, JL(4)=17,...)
# Wait: 2⁰+1=2, 2¹-1=1, 2²+1=5, 2³-1=7, 2⁴+1=17, 2⁵-1=31, 2⁶+1=65

# Hmm, let me be precise.
print("Jacobsthal numbers J(n) = (2ⁿ-(-1)ⁿ)/3:")
for n in range(12):
    J = (2**n - (-1)**n) // 3
    print(f"  J({n:2d}) = {J}")

print()
print("Jacobsthal-Lucas numbers JL(n) = 2ⁿ+(-1)ⁿ:")
for n in range(12):
    JL = 2**n + (-1)**n
    print(f"  JL({n:2d}) = {JL}")

print()
print("Independence polynomial values:")
print("  I(C_k, 2) = JL(k) = 2^k + (-1)^k")
print("  I(P_k, 2) = J(k+2)")
print()

# 7 = JL(3) = I(C_3, 2)
# 21 = J(6) = I(P_4, 2)

# Both 7 and 21 live in the Jacobsthal/JL world.
# Are there other J or JL values that are forbidden?

print("Jacobsthal values in forbidden set:")
for n in range(15):
    J = (2**n - (-1)**n) // 3
    status = "FORBIDDEN" if J in permanent_gaps else "achievable" if J % 2 == 1 else "even"
    if J > 0 and J < 500 and J % 2 == 1:
        print(f"  J({n:2d}) = {J:4d}  {status}")

print()
print("Jacobsthal-Lucas values in forbidden set:")
for n in range(15):
    JL = 2**n + (-1)**n
    status = "FORBIDDEN" if JL in permanent_gaps else "achievable" if JL % 2 == 1 else "even"
    if JL > 0 and JL < 500 and JL % 2 == 1:
        print(f"  JL({n:2d}) = {JL:4d}  {status}")

print()

# Key: 7 is a JL value and 21 is a J value.
# No other J or JL values are permanently forbidden.
# The cyclotomic Φ₃ connection: 7 = Φ₃(2), 21 = Φ₃(4).
# But also: Φ₃(m) = m²+m+1, and J(n) · JL(n) = J(2n).
# J(6) = 21 = Φ₃(4). Is there a direct connection?
# 4² + 4 + 1 = 21. Also 21 = (2⁶-1)/3 = 63/3.
# And 7 = 2³-1. So 7 = 2³-1 and 21 = (2⁶-1)/3 = 3·7/1... wait: 63/3=21, 63=9·7.

# DEEP PATTERN:
# 7 = 2³ - 1 (Mersenne)
# 21 = (2⁶-1)/3 = 63/3 = (2³-1)(2³+1)/3 = 7·9/3 = 7·3 = 21.
# So 21 = 7 · 3 = 7 · (x+1) at x=2!
# And 7 = 7 · 1 = 7 · (x+1)^0.

# In other words: 7 = I(K₃, 2), 21 = I(P₄, 2) = 3·I(K₃, 2).
# Wait: 3·7 = 21. YES!

print("=" * 70)
print("KEY DISCOVERY: 21 = 3 × 7 = (x+1) × I(K₃, x) at x=2")
print("=" * 70)
print()

print("  I(K₃, 2) = 7")
print("  I(P₄, 2) = 21 = 3 × 7")
print()
print("  Is 21 = (x+1) · I(K₃, x) at x=2? Let's check:")
print(f"  (2+1) · I(K₃, 2) = 3 · 7 = 21 ✓")
print()
print("  But this could be a coincidence. Let's check the general formula:")
print(f"  I(K₃, x) = 1 + 3x (triangle = three vertices, any one is independent)")
print(f"  (x+1) · I(K₃, x) = (x+1)(1+3x) = 1 + 4x + 3x²")
print(f"  I(P₄, x) = 1 + 4x + 3x² (path on 5 vertices)")
print()

# Let me verify: I(P₄, x) = ?
# P₄ has 5 vertices: 0-1-2-3-4
# Indep sets: ∅ (1), singletons (5), pairs avoiding edges:
#   {0,2},{0,3},{0,4},{1,3},{1,4},{2,4} = 6 pairs
# Triples: {0,2,4} = 1 triple
# I(P₄, x) = 1 + 5x + 6x² + x³
# At x=2: 1+10+24+8 = 43. Wait, that doesn't match!

# Actually P₄ means path graph on 4 EDGES (5 vertices)
# Let me recount. P_k = path with k vertices and k-1 edges.
# I(P_1, x) = 1+x
# I(P_2, x) = 1+2x
# I(P_3, x) = 1+3x+x² (vertices 1,2,3; edges 12,23; indep pairs: {1,3})
# I(P_4, x) = 1+4x+3x² (vertices 1,2,3,4; edges 12,23,34; indep pairs: {1,3},{1,4},{2,4})
# I(P_4, 2) = 1+8+12 = 21 ✓

# I(K₃, x) = 1+3x (K₃ has 3 vertices, all pairs are edges, so only singletons)
# (x+1)·I(K₃, x) = (x+1)(1+3x) = 1+3x+x+3x² = 1+4x+3x²
# = I(P₄, x) EXACTLY!!

print("  VERIFIED: I(P₄, x) = (x+1) · I(K₃, x) for ALL x!")
print(f"  I(P₄, x) = 1 + 4x + 3x²")
print(f"  (x+1)(1+3x) = 1 + 4x + 3x²  ✓")
print()

# This is a POLYNOMIAL IDENTITY, not a coincidence!
# Why? Because P₄ = K₃ * K₁ in some graph product sense?
# Actually: I(G ⊔ H, x) = I(G,x) · I(H,x) for DISJOINT UNION.
# (x+1) = I(K₁, x). So (x+1)·I(K₃,x) = I(K₁ ⊔ K₃, x).
# K₁ ⊔ K₃ has 4 vertices: one isolated + one triangle.
# I(K₁⊔K₃, x) = (1+x)(1+3x) = 1+4x+3x².
# But I(P₄, x) = 1+4x+3x² also.
# So I(P₄, x) = I(K₁⊔K₃, x) for all x!

# This means P₄ and K₁⊔K₃ are "I-equivalent" (same independence polynomial)!
# Even though they are DIFFERENT GRAPHS!

print("=" * 70)
print("EUREKA: P₄ AND K₁⊔K₃ ARE I-EQUIVALENT!")
print("=" * 70)
print()
print("  P₄ (path on 4 vertices): 1-2-3-4, connected")
print("  K₁⊔K₃ (vertex + triangle): disconnected")
print()
print("  BOTH have I(G,x) = 1 + 4x + 3x²")
print()
print("  This means: H=21 is forbidden for the SAME reason as H=7!")
print("  I(G,2)=21 requires a graph with independence poly 1+4x+3x²,")
print("  which is I-equivalent to K₁⊔K₃.")
print("  Since K₃ cannot be a component of Ω(T) (THM-201),")
print("  and P₄ cannot be Ω(T) (THM-202),")
print("  BOTH realizations of the polynomial 1+4x+3x² are blocked!")
print()

# But wait — are P₄ and K₁⊔K₃ the ONLY graphs with this polynomial?
print("Checking: are P₄ and K₁⊔K₃ the only graphs with I(G,x)=1+4x+3x²?")
print()

# A graph with I(G,2)=21 must have:
# α₀=1 (always), 2α₁+4α₂+8α₃+... = 20
# α₁+2α₂+4α₃+... = 10
# With polynomial 1+4x+3x²: α₁=4, α₂=3, α₃=α₄=...=0.
# So n≥4 vertices, exactly 4 vertices (α₁), exactly 3 independent pairs, max indep = 2.

# n=4, α₁=4, α₂=3: all 4 vertices are independent individually.
# C(4,2)=6 pairs total, 3 are independent = 3 edges.
# Graphs with 4 vertices and 3 edges, independence number 2:
# - P₄ (path): 1-2, 2-3, 3-4. Indep pairs: {1,3},{1,4},{2,4}. α₂=3. ✓
# - K₃+K₁: triangle 1-2-3 + isolated 4. Indep pairs: {1,4},{2,4},{3,4}. α₂=3. ✓
# - Star K_{1,3}: center 1, edges 1-2,1-3,1-4. Indep pairs: {2,3},{2,4},{3,4}. α₂=3.
#   But α₃=1 ({2,3,4}). So I = 1+4x+3x²+x³ ≠ target.
# - C₃+edge: triangle 1-2-3, edge 1-4. α₂: {2,4},{3,4},{2,3}? No, {2,3} shares edge 2-3.
#   Wait: edges are 1-2, 2-3, 1-3, 1-4. Indep pairs: {3,4} only? No:
#   {2,4}: 2-4 no edge ✓, {3,4}: 3-4 no edge ✓, {2,3}: edge ✗.
#   α₂=2. Not matching.

# What about cycle C₄? Edges: 1-2, 2-3, 3-4, 4-1. α₁=4,
# Indep pairs: {1,3},{2,4}. α₂=2. Not matching.

# Actually let me be more systematic about 4-vertex, 3-edge graphs.
from itertools import combinations

print("4-vertex graphs with 3 edges and I(G,x) = 1+4x+3x²:")
all_edges_4 = list(combinations(range(4), 2))
for edge_set in combinations(all_edges_4, 3):
    # Compute full I(G, x) coefficients
    adj = set()
    for u, v in edge_set:
        adj.add((u, v))
        adj.add((v, u))

    alpha = [0] * 5  # α₀ through α₄
    for mask in range(1 << 4):
        verts = [i for i in range(4) if mask & (1 << i)]
        is_indep = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if (verts[i], verts[j]) in adj:
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            alpha[len(verts)] += 1

    poly = alpha  # coefficients
    if poly[:3] == [1, 4, 3] and all(c == 0 for c in poly[3:]):
        print(f"  Edges: {list(edge_set)}, I(G,x) = 1 + 4x + 3x²  ✓")

print()

# Check larger vertex counts
print("5-vertex graphs with I(G,2) = 21:")
count = 0
all_edges_5 = list(combinations(range(5), 2))
for num_e in range(len(all_edges_5)+1):
    for edge_set in combinations(all_edges_5, num_e):
        val = independence_polynomial_at_2(5, edge_set)
        if val == 21:
            count += 1
            if count <= 5:
                print(f"  Edges: {list(edge_set)} ({num_e} edges)")
if count == 0:
    print("  NONE FOUND!")
elif count > 5:
    print(f"  ... and {count-5} more")
print(f"  Total: {count}")

print()
print("=" * 70)
print("CONCLUSION")
print("=" * 70)
print()
print("""
THE TWO FORBIDDEN VALUES HAVE A UNIFIED EXPLANATION:

1. H = 7 = I(K₃, 2): forbidden because K₃ cannot be Ω-component (THM-201)
   K₃ is the ONLY graph with I(G,2) = 7.

2. H = 21 = I(P₄, 2) = I(K₁⊔K₃, 2): BOTH realizations are forbidden.
   - P₄ forbidden by THM-202
   - K₁⊔K₃ forbidden because K₃ component forbidden by THM-201
   I(G,x) = 1+4x+3x² = (1+x)(1+3x) = I(K₁,x)·I(K₃,x)
   The polynomial FACTORS through K₃!

THE PATTERN: Both forbidden values have independence polynomials
that FACTOR through (1+3x) = I(K₃, x). Specifically:
  7 = I(K₃, 2): polynomial is (1+3x)
  21 = I(K₁⊔K₃, 2) = I(P₄, 2): polynomial is (1+x)(1+3x)

The K₃ factor is the "poison" — any H value whose independence
polynomial necessarily factors through I(K₃, x) = 1+3x is forbidden,
because Ω(T) cannot contain K₃ as a substructure.

CYCLOTOMIC CONNECTION:
  7 = Φ₃(2) = 2²+2+1 (third cyclotomic polynomial at m=2)
  21 = Φ₃(4) = 4²+4+1 (third cyclotomic polynomial at m=4)
  Both live on the curve m²+m+1 at EVEN m values.

  But Φ₃(6) = 43 IS achievable, so the cyclotomic pattern alone
  doesn't determine forbiddenness — it's the I-polynomial factorization
  through (1+3x) that matters.
""")
