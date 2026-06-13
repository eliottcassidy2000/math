"""
five_fold_petersen_map.py -- kind-pasteur-2026-03-14-S67

The Petersen graph has 10 vertices = C(5,2) = 2-element subsets of {1,2,3,4,5}.
The five exceptional Lie groups correspond to five elements of {1,2,3,4,5}.
Each Petersen vertex {i,j} maps to a PAIR of exceptional groups.

This script explores what graph-theoretic properties of the Petersen graph
correspond to algebraic properties of the exceptional Lie groups under
various assignments.

KEY IDEA: assign {1,2,3,4,5} -> {G2, F4, E6, E7, E8} and study which
assignment makes the Petersen structure most meaningful.
"""

from itertools import permutations, combinations
import numpy as np

# Exceptional Lie groups
GROUPS = {
    'G2': {'rank': 2, 'dim': 14, 'h': 6, 'roots_pos': 6, 'det': 1},
    'F4': {'rank': 4, 'dim': 52, 'h': 12, 'roots_pos': 24, 'det': 1},
    'E6': {'rank': 6, 'dim': 78, 'h': 12, 'roots_pos': 36, 'det': 3},
    'E7': {'rank': 7, 'dim': 133, 'h': 18, 'roots_pos': 63, 'det': 2},
    'E8': {'rank': 8, 'dim': 248, 'h': 30, 'roots_pos': 120, 'det': 1},
}
GROUP_NAMES = ['G2', 'F4', 'E6', 'E7', 'E8']

# Petersen graph: K(5,2)
# Vertices: 2-element subsets of {0,1,2,3,4}
# Edge: disjoint subsets
pet_vertices = list(combinations(range(5), 2))
pet_edges = []
for i in range(len(pet_vertices)):
    for j in range(i+1, len(pet_vertices)):
        if len(set(pet_vertices[i]) & set(pet_vertices[j])) == 0:
            pet_edges.append((i, j))

print("=" * 70)
print("FIVE-FOLD PETERSEN MAP: EXCEPTIONAL GROUPS ON K(5,2)")
print("=" * 70)

print(f"\nPetersen K(5,2) has {len(pet_vertices)} vertices and {len(pet_edges)} edges")
print(f"Each vertex = {{i,j}} subset of {{0,1,2,3,4}}")
print(f"Edge connects {{i,j}} and {{k,l}} iff disjoint (i.e., {{i,j,k,l}} = 4 elements)")
print()

# Under assignment sigma: {0,1,2,3,4} -> {G2,F4,E6,E7,E8}
# vertex {i,j} -> pair (sigma(i), sigma(j))
# edge {i,j}--{k,l} -> the two PAIRS are disjoint, using 4 of 5 groups

# For each Petersen vertex (pair of groups), compute interaction invariants
def pair_invariant(g1_name, g2_name, invariant_fn):
    """Compute some invariant from a pair of groups."""
    g1, g2 = GROUPS[g1_name], GROUPS[g2_name]
    return invariant_fn(g1, g2)

# Natural assignment: ordered by rank
# G2(r=2), F4(r=4), E6(r=6), E7(r=7), E8(r=8)
natural_assignment = ['G2', 'F4', 'E6', 'E7', 'E8']

print("NATURAL ASSIGNMENT (ordered by rank):")
print("  0->G2, 1->F4, 2->E6, 3->E7, 4->E8")
print()

# For each vertex, show the pair and its invariants
print("  Vertex  | Pair       | Sum_r | Prod_r | Sum_h | Prod_h | Sum_dim | Sum_roots+")
print("  " + "-"*85)

for idx, (i, j) in enumerate(pet_vertices):
    g1, g2 = natural_assignment[i], natural_assignment[j]
    d1, d2 = GROUPS[g1], GROUPS[g2]
    sr = d1['rank'] + d2['rank']
    pr = d1['rank'] * d2['rank']
    sh = d1['h'] + d2['h']
    ph = d1['h'] * d2['h']
    sd = d1['dim'] + d2['dim']
    srp = d1['roots_pos'] + d2['roots_pos']
    print(f"  {{{i},{j}}}={idx:2d} | ({g1},{g2}){'':<{5-len(g1)-len(g2)}} | {sr:5d} | {pr:6d} | {sh:5d} | {ph:6d} | {sd:7d} | {srp:10d}")

# Check which edges (pairs of disjoint pairs) have nice properties
print(f"\nEDGES (disjoint pairs of pairs = using 4 of 5 groups):")
print(f"  Edge            | Missing | Sum_4r | Prod_4r | Sum_4h")
print("  " + "-"*65)

for ei, ej in pet_edges:
    i1, j1 = pet_vertices[ei]
    i2, j2 = pet_vertices[ej]
    groups_used = set([i1, j1, i2, j2])
    missing = (set(range(5)) - groups_used).pop()
    g_missing = natural_assignment[missing]

    g_used = [natural_assignment[k] for k in [i1, j1, i2, j2]]
    sr = sum(GROUPS[g]['rank'] for g in g_used)
    pr = 1
    for g in g_used:
        pr *= GROUPS[g]['rank']
    sh = sum(GROUPS[g]['h'] for g in g_used)

    # This vertex pair plus that vertex pair = 4 groups
    # Missing group is the "excluded" one
    # In Petersen, each vertex has degree 3 and non-edge has degree 6
    # Non-edge = overlapping pairs

    print(f"  {{{i1},{j1}}}--{{{i2},{j2}}} | miss={missing}({g_missing}) | {sr:6d} | {pr:7d} | {sh:6d}")

# Now the KEY question: does the Petersen adjacency encode anything about Lie pairs?
print(f"\n{'='*70}")
print("PETERSEN ADJACENCY AS LIE INTERACTION")
print("=" * 70)

print("""
In Petersen K(5,2), vertices {i,j} and {k,l} are:
  ADJACENT (edge) iff {i,j} and {k,l} are DISJOINT
  NON-ADJACENT iff they SHARE an element

For exceptional Lie groups:
  Edge (disjoint pairs): the 4 groups are "independent" (no shared element)
  Non-edge (overlapping): the groups share a "common root" element

What makes two exceptional groups "interact"?
  - F4/E6 share h=12 (same Coxeter number)
  - E6/E7/E8 are the E-series (share Dynkin structure)
  - G2/F4 are non-simply-laced
  - G2/E8 are the "bookends" (h ratio = 5)
""")

# Check: which assignment makes the maximum independent set
# correspond to a meaningful Lie sub-family?
print("MAXIMUM INDEPENDENT SETS (size 4 = alpha(Petersen)):")
print("Each is a 'star' — all pairs containing a fixed element")
for center in range(5):
    star = [pet_vertices[idx] for idx in range(len(pet_vertices))
            if center in pet_vertices[idx]]
    g_center = natural_assignment[center]
    others = [natural_assignment[v[0] if v[1] == center else v[1]] for v in star]
    print(f"  Center={center}({g_center}): pairs with {others}")
    # These 4 pairs all contain g_center
    # They form an independent set because any two share g_center

print(f"""
Under the natural assignment, the 5 max independent sets are:
  Star(G2): {{G2,F4}}, {{G2,E6}}, {{G2,E7}}, {{G2,E8}} — G2 paired with everything
  Star(F4): {{F4,G2}}, {{F4,E6}}, {{F4,E7}}, {{F4,E8}} — F4 paired with everything
  Star(E6): {{E6,G2}}, {{E6,F4}}, {{E6,E7}}, {{E6,E8}} — E6 paired with everything
  Star(E7): {{E7,G2}}, {{E7,F4}}, {{E7,E6}}, {{E7,E8}} — E7 paired with everything
  Star(E8): {{E8,G2}}, {{E8,F4}}, {{E8,E6}}, {{E8,E7}} — E8 paired with everything

The "dependent core" (6 vertices not in the star):
  When center=E8: the core is all pairs NOT involving E8
  = {{G2,F4}}, {{G2,E6}}, {{G2,E7}}, {{F4,E6}}, {{F4,E7}}, {{E6,E7}}
  = all pairs from {{G2,F4,E6,E7}} = K_4 on the 4 non-E8 groups
  This K_4 has 6 vertices and 3 edges in Petersen (= 3 disjoint pairs among 4 groups)
""")

# The h-sum and h-product for Petersen edges vs non-edges
print("=" * 70)
print("h-SUMS ON PETERSEN EDGES vs NON-EDGES")
print("=" * 70)

h_vals = [GROUPS[g]['h'] for g in natural_assignment]  # [6, 12, 12, 18, 30]

edge_h_sums = []
nonedge_h_sums = []

for i in range(10):
    for j in range(i+1, 10):
        vi, vj = pet_vertices[i], pet_vertices[j]
        h_sum = h_vals[vi[0]] + h_vals[vi[1]] + h_vals[vj[0]] + h_vals[vj[1]]
        is_edge = (i, j) in pet_edges or (j, i) in pet_edges
        h_pair = h_vals[vi[0]] + h_vals[vi[1]]
        h_pair2 = h_vals[vj[0]] + h_vals[vj[1]]

        if is_edge:
            edge_h_sums.append(h_pair + h_pair2)
        else:
            nonedge_h_sums.append(h_pair + h_pair2)

print(f"\nh-pair sums for each vertex:")
for idx, (i, j) in enumerate(pet_vertices):
    h_sum = h_vals[i] + h_vals[j]
    g1, g2 = natural_assignment[i], natural_assignment[j]
    lie_interp = ""
    if h_sum == 18: lie_interp = "= h(E7)"
    elif h_sum == 24: lie_interp = "= h(F4) + h(E6)"
    elif h_sum == 30: lie_interp = "= h(E8)"
    elif h_sum == 36: lie_interp = "= 2*h(E7)"
    elif h_sum == 42: lie_interp = "= h(G2)*7 = h(G2)*rank(E7)"
    elif h_sum == 48: lie_interp = "= 4*h(F4)"
    elif h_sum == 12: lie_interp = "= h(F4)"
    print(f"  {{{i},{j}}} ({g1},{g2}): h_sum = {h_sum} {lie_interp}")

print(f"\nSum of ALL h-pair sums = {sum(h_vals[i]+h_vals[j] for i,j in pet_vertices)}")
print(f"  = sum over each element of h_i * (number of pairs containing i)")
print(f"  = sum h_i * 4 = 4 * {sum(h_vals)} = {4 * sum(h_vals)}")
total_h = sum(h_vals)
print(f"  Total h = {total_h} = {'+'.join(str(h) for h in h_vals)} = {total_h}")
print(f"  4 * total_h = {4*total_h}")

# Now check: what invariant is CONSTANT on Petersen edges?
print(f"\n{'='*70}")
print("INVARIANT CONSTANCY ON PETERSEN EDGES")
print("=" * 70)

print("\nFor each edge {{i,j}}--{{k,l}}, the missing element m has h_m:")
for ei, ej in pet_edges:
    i1, j1 = pet_vertices[ei]
    i2, j2 = pet_vertices[ej]
    used = {i1, j1, i2, j2}
    missing = (set(range(5)) - used).pop()
    h_missing = h_vals[missing]
    h_used_sum = sum(h_vals[k] for k in used)
    h_used_prod = 1
    for k in used:
        h_used_prod *= h_vals[k]
    r_used_sum = sum(GROUPS[natural_assignment[k]]['rank'] for k in used)
    print(f"  miss={missing}({natural_assignment[missing]}): "
          f"h_miss={h_missing}, h_used_sum={h_used_sum}, "
          f"r_used_sum={r_used_sum}, h_prod={h_used_prod}")

print(f"""
KEY OBSERVATION:
  Each Petersen EDGE excludes exactly one of the 5 groups.
  The edges partition into 5 classes by which group is excluded:
    miss=0 (G2):  3 edges, h_miss=6
    miss=1 (F4):  3 edges, h_miss=12
    miss=2 (E6):  3 edges, h_miss=12
    miss=3 (E7):  3 edges, h_miss=18
    miss=4 (E8):  3 edges, h_miss=30

  Each group is "missed" by exactly 3 edges = degree of K(5,2).
  h_used_sum = total_h - h_miss = {total_h} - h_miss
""")

# The complementary view: Johnson graph J(5,2) has 30 edges
# Each J(5,2) edge connects overlapping pairs
print("=" * 70)
print("RANK PRODUCTS ON PETERSEN COMPLEMENT (JOHNSON J(5,2))")
print("=" * 70)

r_vals = [GROUPS[g]['rank'] for g in natural_assignment]  # [2, 4, 6, 7, 8]

print(f"\nRanks: {list(zip(natural_assignment, r_vals))}")
print(f"\nFor J(5,2) edges (overlapping pairs, share element k):")
print(f"  Shared element gives 'interaction rank' = rank of shared group")

# J(5,2) edges: vertices that share an element
j_edges = []
for i in range(10):
    for j in range(i+1, 10):
        if len(set(pet_vertices[i]) & set(pet_vertices[j])) > 0:
            j_edges.append((i, j))

# For each shared element, count the edges
for shared in range(5):
    edges_via = [(i, j) for i, j in j_edges
                 if shared in set(pet_vertices[i]) & set(pet_vertices[j])]
    g_shared = natural_assignment[shared]
    r_shared = r_vals[shared]
    print(f"  Shared={shared}({g_shared}, r={r_shared}): "
          f"{len(edges_via)} edges (C(4,1)·... = 6 edges each)")

print(f"\nTotal J(5,2) edges: {len(j_edges)} = h(E8) = 30")

# The rank polynomial: sum of rank products over J(5,2) edges
rank_prod_sum = sum(r_vals[list(set(pet_vertices[i]) & set(pet_vertices[j]))[0]]
                    for i, j in j_edges
                    if len(set(pet_vertices[i]) & set(pet_vertices[j])) == 1)
print(f"Sum of interaction ranks: {rank_prod_sum}")
# This should equal sum_k r_k * C(4,1) = ... no
# Each element k is shared by C(4,1) = 4 pairs of vertices
# Actually: element k appears in 4 vertices. These form K_4 in J(5,2) = 6 edges each.
# So each k contributes 6 J(5,2) edges.
print(f"Expected: sum_k r_k * 6 = 6 * {sum(r_vals)} = {6*sum(r_vals)}")
# That's not right either — it's not r_k per edge, it's just counting.

# More interesting: product of the two pair-ranks for each J(5,2) edge
print(f"\nFor each J(5,2) edge (overlapping pairs {{a,b}}--{{a,c}}):")
print(f"  The interaction is mediated by group a")
print(f"  Pair ranks: r_a+r_b and r_a+r_c")
print(f"  Product: (r_a+r_b)(r_a+r_c)")

# Actually the most natural invariant is the independence polynomial of Petersen
# evaluated at x = some Lie-derived value
print(f"\n{'='*70}")
print("I(Petersen, x) AT LIE-DERIVED VALUES")
print("=" * 70)

# I(P,x) = 1 + 10x + 30x^2 + 30x^3 + 5x^4
def I_pet(x):
    return 1 + 10*x + 30*x**2 + 30*x**3 + 5*x**4

# At each h value and rank
for name, data in GROUPS.items():
    r = data['rank']
    h = data['h']
    print(f"  I(P, h({name})) = I(P, {h}) = {I_pet(h)}")

print()
for name, data in GROUPS.items():
    r = data['rank']
    print(f"  I(P, r({name})) = I(P, {r}) = {I_pet(r)}")

print()
for name, data in GROUPS.items():
    r = data['rank']
    h = data['h']
    h1 = h + 1
    print(f"  I(P, h+1({name})) = I(P, {h1}) = {I_pet(h1)}")

print()
# I(P, 2) = 461 — this would be H if Petersen were a tournament CG
print(f"  I(P, KEY_1=2) = {I_pet(2)} (= H if Petersen were a CG)")
print(f"  I(P, KEY_2=3) = {I_pet(3)}")
print(f"  I(P, -1) = {I_pet(-1)} (Euler char of indep complex)")

# Check: I(P, -KEY_1) and I(P, -KEY_2)
print(f"  I(P, -2) = {I_pet(-2)}")
print(f"  I(P, -3) = {I_pet(-3)}")

# The most remarkable evaluations
print(f"\nSPECIAL VALUES:")
print(f"  I(P, 0) = 1")
print(f"  I(P, 1) = 76 = 4 * 19 = rank(F4) * (h(E7)+1)")
print(f"  I(P, -1) = -4 = -rank(F4)")
print(f"  I(P, 2) = 461 = prime")
print(f"  I(P, 3) = 1516 = 4 * 379 = rank(F4) * 379")

# Is I(P,1) = 76 meaningful? 76 = sum of first 4 triangular numbers... no
# 76 = T(12) - T(4) ... no
# 76 = 2^2 * 19 = KEY_1^2 * (h(E7)+1)
print(f"  76 = KEY_1^2 * (h(E7)+1) = 4 * 19")

# Summary
print(f"\n{'='*70}")
print("FIVE-FOLD SYNTHESIS")
print("=" * 70)
print(f"""
The Petersen graph K(5,2) encodes the exceptional Lie groups through:

1. VERTEX = PAIR of groups
   10 vertices = C(5,2) = all pairs from {{G2,F4,E6,E7,E8}}

2. EDGE = DISJOINT PAIRS (complementary, using 4 of 5 groups)
   15 edges, each missing exactly 1 group
   3 edges miss each group (= degree of Petersen)

3. NON-EDGE = OVERLAPPING PAIRS (interacting, sharing 1 group)
   30 non-edges = h(E8) = Coxeter number of E8
   6 non-edges share each group (= degree of J(5,2))

4. INDEPENDENT SETS = families of pairwise non-interacting pairs
   Maximum independent set = 4 pairs all containing the same group
   = the "star" of one exceptional group
   alpha = 4 = rank(F4)

5. EIGENVALUE SPECTRUM
   {{3, 1^5, -2^4}} = {{KEY_2, 1, -KEY_1}}
   These are the roots of z^3-2z^2-5z+6 = 0
   whose coefficients are {{1, -KEY_1, -(KEY_1+KEY_2), h(G2)}}

6. PRODUCTS WITH TOURNAMENT POLYNOMIAL
   (z-3)(z^2-5z+6) has coefficients {{1, -rank(E8), 21, -h(E7)}}
   (z+2)(z^2-5z+6) has coefficients {{1, -KEY_2, -rank(F4), h(F4)}}
   Both encode Lie data through elementary symmetric polynomials

7. CYCLOTOMIC BRIDGE
   h+1 = Phi_3(key) for 4 of 5 exceptional groups:
     Phi_3(KEY_1) = 7 = h(G2)+1
     Phi_3(KEY_2) = 13 = h(F4)+1 = h(E6)+1
     Phi_3(KEY_1+KEY_2) = 31 = h(E8)+1
   E7 is the outlier: h(E7)+1 = 19 = 3^3-2^3 (recurrence formula)

8. FORBIDDEN H VALUES
   H=7 = Phi_3(KEY_1) and H=21 = Phi_3(KEY_1^2)
   are exactly the H values where the CG would need to be K_3 or K_3+K_1
   — both have exactly KEY_2 = 3 edges, which is FORBIDDEN
   by the tournament structure's forcing theorems.
""")

if __name__ == "__main__":
    pass
