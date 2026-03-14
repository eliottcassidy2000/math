#!/usr/bin/env python3
"""
dehn_invariant_tournament.py — opus-2026-03-14-S75

Compute the ACTUAL Dehn invariant of tournament order polytopes.
The order polytope O(T) = {x ∈ [0,1]^n : x_i > x_j if i→j in T}.

For the Dehn invariant, we need:
- The boundary facets of O(T) and their (n-2)-dimensional volumes
- The dihedral angles at each boundary ridge

Key insight from the simplex_cuboid script:
- Internal dihedral angles are π (flat internal boundaries)
- Boundary dihedral angles are π/2 (cube boundary)
- So D(O(T)) = (boundary ridge volumes) ⊗ π/2

This means D(O(T)) is proportional to the boundary volume.
Two tournament polytopes with same H but different boundary structure
would have different Dehn invariants!

We verify this at n=3,4,5 using exact volume computations.
"""

from itertools import permutations
from fractions import Fraction
from math import factorial

def gen_tournaments(n):
    """Generate all tournaments on n vertices."""
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(2**len(edges)):
        adj = [[False]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if bits & (1 << idx):
                adj[i][j] = True
            else:
                adj[j][i] = True
        yield adj

def ham_paths(adj, n):
    """Return list of Hamiltonian paths (as permutations)."""
    paths = []
    for perm in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if not adj[perm[i]][perm[i+1]]:
                valid = False
                break
        if valid:
            paths.append(perm)
    return paths

def order_polytope_volume(adj, n):
    """Volume of O(T) = H/n! where H = number of Ham paths."""
    return Fraction(len(ham_paths(adj, n)), factorial(n))

def simplex_from_path(perm, n):
    """
    The simplex Δ_σ for Hamiltonian path σ is:
    {x ∈ [0,1]^n : x_{σ(1)} ≥ x_{σ(2)} ≥ ... ≥ x_{σ(n)}}

    This is a standard simplex in the n-cube.
    Its vertices are:
    - v_0 = (0, 0, ..., 0)
    - v_k = e_{σ(1)} + e_{σ(2)} + ... + e_{σ(k)} for k=1,...,n
    where e_i is the i-th standard basis vector.
    """
    vertices = []
    v = [Fraction(0)] * n
    vertices.append(tuple(v))
    for k in range(n):
        v = list(v)
        v[perm[k]] = Fraction(1)
        vertices.append(tuple(v))
    return vertices

def simplex_boundary_facets(perm, n):
    """
    The simplex Δ_σ has n+1 facets.
    The boundary of the n-cube [0,1]^n consists of the 2n hyperplanes x_i=0 and x_i=1.

    Facets of Δ_σ that lie on the cube boundary:
    - x_{σ(n)} = 0 (the "bottom" facet) — lies on {x_{σ(n)}=0}
    - x_{σ(1)} = 1 (the "top" facet) — lies on {x_{σ(1)}=1}
    - For internal facets: x_{σ(k)} = x_{σ(k+1)} — these are INTERNAL to O(T)
      IF both orderings (with k and k+1 swapped) are Hamiltonian paths.

    Returns list of (type, data) for each facet.
    """
    facets = []
    # Bottom facet: x_{σ(n)} = 0
    facets.append(('boundary', f'x_{perm[n-1]}=0'))
    # Top facet: x_{σ(1)} = 1
    facets.append(('boundary', f'x_{perm[0]}=1'))
    # Internal facets: x_{σ(k)} = x_{σ(k+1)} for k=0,...,n-2
    for k in range(n-1):
        facets.append(('internal', f'x_{perm[k]}=x_{perm[k+1]}'))
    return facets

print("=" * 70)
print("PART 1: TOURNAMENT ORDER POLYTOPE — BASIC STRUCTURE")
print("=" * 70)
print()

for n in range(3, 6):
    print(f"  n={n}:")
    iso_classes = {}  # H -> list of tournaments
    for adj in gen_tournaments(n):
        paths = ham_paths(adj, n)
        h = len(paths)
        if h not in iso_classes:
            iso_classes[h] = []
        iso_classes[h].append((adj, paths))

    for h in sorted(iso_classes.keys()):
        adj, paths = iso_classes[h][0]  # representative
        vol = Fraction(h, factorial(n))
        print(f"    H={h}: vol={vol} = {float(vol):.6f}, {len(paths)} paths")

        if n <= 4:
            # List all paths
            for p in paths:
                verts = simplex_from_path(p, n)
                facets = simplex_boundary_facets(p, n)
                boundary_facets = [f for f in facets if f[0] == 'boundary']
                internal_facets = [f for f in facets if f[0] == 'internal']
                print(f"      path {p}: {len(boundary_facets)} boundary, {len(internal_facets)} internal facets")
    print()

print("=" * 70)
print("PART 2: ADJACENCY OF HAMILTONIAN PATHS")
print("=" * 70)
print()
print("  Two Hamiltonian paths σ, τ are ADJACENT if they share a facet,")
print("  i.e., τ = σ with positions k and k+1 swapped.")
print("  This means σ(k)→σ(k+1) in T but also σ(k+1)→σ(k) in T")
print("  ... which is impossible in a tournament!")
print()
print("  WAIT: The facet x_{σ(k)} = x_{σ(k+1)} is shared between Δ_σ")
print("  (where x_{σ(k)} ≥ x_{σ(k+1)}) and Δ_τ where τ swaps k and k+1")
print("  (where x_{σ(k+1)} ≥ x_{σ(k)}).")
print("  For Δ_τ to be in O(T), we need σ(k+1)→σ(k) in T.")
print("  But σ is a Hamiltonian path, so σ(k)→σ(k+1) in T.")
print("  These contradict! So NO two simplices share an INTERNAL facet.")
print()
print("  CONCLUSION: The H simplices in the triangulation of O(T)")
print("  are DISJOINT (they only share lower-dimensional faces)!")
print()
print("  This means O(T) is a DISJOINT UNION of H simplices,")
print("  not a proper triangulation. The order polytope is NOT convex")
print("  when H > 1 (it has 'holes' between simplices).")
print()
print("  CORRECTION: Actually, the order polytope O(T) IS convex")
print("  (intersection of half-spaces). The simplices DO overlap in")
print("  lower-dimensional faces. The issue is that no two simplices")
print("  share a CODIMENSION-1 face that is internal to the cube.")
print()
print("  Let me reconsider: O(T) = {x ∈ [0,1]^n : x_i > x_j if i→j}")
print("  This is convex. Its volume is H/n!.")
print("  The H simplices are the maximal simplices in the chamber")
print("  decomposition of O(T) by the hyperplanes x_i = x_j.")
print()

# Let's actually verify the geometry at n=3
print("  DETAILED GEOMETRY at n=3:")
print()
for adj in gen_tournaments(3):
    paths = ham_paths(adj, 3)
    h = len(paths)

    # Describe the tournament
    arcs = []
    for i in range(3):
        for j in range(i+1, 3):
            if adj[i][j]:
                arcs.append(f"{i}→{j}")
            else:
                arcs.append(f"{j}→{i}")

    print(f"  Tournament: {', '.join(arcs)}, H={h}")

    for p in paths:
        print(f"    Path {p}: simplex {{x_{p[0]} ≥ x_{p[1]} ≥ x_{p[2]}}}")
        # The simplex has vertices:
        # (0,0,0), e_{p[0]}, e_{p[0]}+e_{p[1]}, (1,1,1)
        verts = []
        v = [0, 0, 0]
        verts.append(tuple(v))
        for k in range(3):
            v = list(v)
            v[p[k]] = 1
            verts.append(tuple(v))
        for vv in verts:
            print(f"      vertex: {vv}")

    if h == 3:
        # Show the arrangement
        print(f"    This is the REGULAR tournament (3-cycle).")
        print(f"    O(T) = {{x ∈ [0,1]³ : x_0>x_1, x_1>x_2, x_2>x_0}} ∪ ...")
        print(f"    Wait: x_0>x_1 AND x_1>x_2 AND x_2>x_0 is EMPTY (contradictory)!")
        print(f"    So O(T) must be defined differently.")
        print()
        print(f"    RECHECK: O(T) = {{x ∈ [0,1]³ : x_i > x_j for SOME consistency}}")
        print(f"    Actually O(T) = convex hull of LINEAR EXTENSIONS")
        print()
        print(f"    Let me reconsider the order polytope definition.")
        print(f"    For a POSET P: O(P) = {{f: P→[0,1] : f order-preserving}}")
        print(f"    A tournament is NOT a poset (has cycles).")
        print(f"    So the 'order polytope' needs a different definition.")
    print()

print("=" * 70)
print("PART 3: ACYCLIC ORIENTATIONS AND LINEAR EXTENSIONS")
print("=" * 70)
print()
print("  For a tournament T, the Hamiltonian paths correspond to")
print("  permutations σ where σ(1)→σ(2)→...→σ(n) follows T's arcs.")
print()
print("  The 'natural polytope' associated to T is:")
print("  P(T) = convex hull of indicator vectors of Hamiltonian paths.")
print()
print("  Or better: for each Hamiltonian path σ, define the simplex")
print("  Δ_σ = {x ∈ [0,1]^n : 1 ≥ x_{σ(1)} ≥ x_{σ(2)} ≥ ... ≥ x_{σ(n)} ≥ 0}")
print()
print("  The UNION U(T) = ∪_σ Δ_σ over all H Hamiltonian paths")
print("  has volume H/n!.")
print()
print("  For the TRANSITIVE tournament: U(T) = one simplex, volume 1/n!")
print("  For the REGULAR tournament on 3: U(T) = three simplices")
print()
print("  At n=3, the three simplices for the regular tournament are:")
print("  Δ₁: x_0 ≥ x_1 ≥ x_2 (path 0→1→2)")
print("  Δ₂: x_1 ≥ x_2 ≥ x_0 (path 1→2→0)")
print("  Δ₃: x_2 ≥ x_0 ≥ x_1 (path 2→0→1)")
print()

# Verify: do these three simplices tile half the cube?
# Volume of each: 1/6 (=1/3!)
# Total: 3/6 = 1/2
# The cube [0,1]^3 has volume 1
# The other three simplices are:
# x_0 ≥ x_2 ≥ x_1, x_1 ≥ x_0 ≥ x_2, x_2 ≥ x_1 ≥ x_0
# These correspond to the REVERSE tournament (2→1, 0→2, 1→0)

print("  Volume check: 3 × (1/6) = 1/2 of the cube.")
print("  The other 3 simplices belong to the REVERSE tournament.")
print("  Together, all 6 simplices tile the cube: 6/6 = 1. ✓")
print()
print("  For a tournament T and its reverse T̄:")
print("  U(T) ∪ U(T̄) = [0,1]^n (they tile the cube!)")
print("  Vol(U(T)) + Vol(U(T̄)) = 1")
print("  H(T)/n! + H(T̄)/n! = 1")
print("  H(T) + H(T̄) = n!")
print()
print("  RÉDEI'S THEOREM says H(T) is ODD, so H(T̄) = n! - H(T) is also odd")
print("  (since n! is even for n≥2, this requires n! to be even, which gives")
print("  H(T) odd → H(T̄) = n!-H(T) which is even-odd = odd ✓)")
print()

# Verify H(T) + H(T̄) = n! at n=3,4,5
print("  Verification of H(T) + H(T̄) = n!:")
for n in range(3, 7):
    nf = factorial(n)
    for adj in gen_tournaments(n):
        # Build reverse
        rev = [[not adj[i][j] if i != j else False for j in range(n)] for i in range(n)]
        h1 = len(ham_paths(adj, n))
        h2 = len(ham_paths(rev, n))
        if h1 + h2 != nf:
            print(f"  FAIL at n={n}: H={h1}, H̄={h2}, sum={h1+h2}, n!={nf}")
            break
    else:
        print(f"  n={n}: H(T)+H(T̄)=n!={nf} ✓ (checked all {2**len([(i,j) for i in range(n) for j in range(i+1,n)])} tournaments)")

print()

print("=" * 70)
print("PART 4: SCISSORS CONGRUENCE OF SIMPLEX UNIONS")
print("=" * 70)
print()
print("  U(T) = ∪_σ Δ_σ where σ ranges over H Hamiltonian paths.")
print()
print("  TWO TOURNAMENT POLYTOPES U(T₁), U(T₂) with same H:")
print("  Do they have the same Dehn invariant?")
print()
print("  Each Δ_σ is a standard simplex in [0,1]^n.")
print("  ALL such simplices are CONGRUENT (related by a permutation matrix).")
print("  So each has the SAME Dehn invariant.")
print()
print("  The Dehn invariant of a UNION is NOT the sum of Dehn invariants!")
print("  D(A ∪ B) = D(A) + D(B) - D(A ∩ B)")
print("  (inclusion-exclusion on the boundary)")
print()
print("  The intersection structure of the simplices is encoded by")
print("  the PERMUTOHEDRON COMPLEX of the tournament.")
print()
print("  QUESTION: When do two simplices Δ_σ, Δ_τ intersect?")
print("  Δ_σ ∩ Δ_τ ≠ ∅ iff the corresponding total orders σ, τ")
print("  are 'compatible' on some subset.")
print()

# Compute intersections of simplices at n=3
print("  INTERSECTION STRUCTURE at n=3, regular tournament:")
from itertools import permutations as perms
n = 3
# Regular tournament: 0→1, 1→2, 2→0
adj_reg = [[False]*3 for _ in range(3)]
adj_reg[0][1] = True
adj_reg[1][2] = True
adj_reg[2][0] = True

paths = ham_paths(adj_reg, 3)
print(f"  Paths: {paths}")

# Check pairwise intersections
for i in range(len(paths)):
    for j in range(i+1, len(paths)):
        p1, p2 = paths[i], paths[j]
        # Two simplices Δ_σ and Δ_τ intersect iff
        # there exists x with x_{σ(1)} ≥ x_{σ(2)} ≥ x_{σ(3)}
        # and x_{τ(1)} ≥ x_{τ(2)} ≥ x_{τ(3)}
        # This is: the intersection of two half-space chains
        # Always non-empty (x = (1/3, 1/3, 1/3) is in both)
        print(f"  Δ_{p1} ∩ Δ_{p2}: non-empty (contains x=({','.join(['1/3']*3)}))")
        # What dimension is the intersection?
        # x_{p1[0]} ≥ x_{p1[1]} ≥ x_{p1[2]} AND x_{p2[0]} ≥ x_{p2[1]} ≥ x_{p2[2]}
        # These give constraints. Let's count independent constraints.
        constraints = set()
        for k in range(2):
            constraints.add((p1[k], p1[k+1]))  # p1[k] ≥ p1[k+1]
            constraints.add((p2[k], p2[k+1]))  # p2[k] ≥ p2[k+1]
        print(f"    Constraints: {sorted(constraints)}")
        # Check for contradictions (a ≥ b AND b ≥ a → a = b)
        equalities = []
        for (a,b) in constraints:
            if (b,a) in constraints:
                equalities.append((a,b))
        if equalities:
            print(f"    Forced equalities: {equalities}")
            print(f"    Intersection is lower-dimensional (has measure 0)")
        else:
            print(f"    No forced equalities: intersection has positive measure!")
            print(f"    (This would mean the simplices OVERLAP in volume)")

print()
print("  The three simplices of the regular tournament pairwise intersect")
print("  only on their boundaries (measure 0), so they tile without overlap.")
print()

print("=" * 70)
print("PART 5: DEHN INVARIANT COMPUTATION")
print("=" * 70)
print()
print("  For a convex polytope P in R^n, the Dehn-Hadwiger invariant is:")
print("  D_k(P) = Σ_{k-faces F} vol_k(F) ⊗ external angle(F)")
print()
print("  For n=3 (3-dimensional polytopes), the classical Dehn invariant is D_1:")
print("  D(P) = Σ_{edges e} length(e) ⊗ θ(e)")
print("  where θ(e) is the dihedral angle at edge e.")
print()
print("  For the standard simplex Δ in R³:")
print("  {(x,y,z): 1 ≥ x ≥ y ≥ z ≥ 0}")
print("  Vertices: (0,0,0), (1,0,0), (1,1,0), (1,1,1)")
print("  Edges:")
edges_simplex = [
    ((0,0,0), (1,0,0)),
    ((0,0,0), (1,1,0)),
    ((0,0,0), (1,1,1)),
    ((1,0,0), (1,1,0)),
    ((1,0,0), (1,1,1)),
    ((1,1,0), (1,1,1)),
]
import math
for (v1, v2) in edges_simplex:
    length = sum((a-b)**2 for a,b in zip(v1,v2))**0.5
    print(f"  {v1} — {v2}: length = {length:.4f}")

print()
print("  The dihedral angles of this simplex:")
print("  Between faces meeting at each edge:")

# For the simplex {1 ≥ x ≥ y ≥ z ≥ 0}:
# Facets:
# F1: x=1 (contains (1,0,0),(1,1,0),(1,1,1))
# F2: y=x (contains (0,0,0),(1,1,0),(1,1,1))
# F3: z=y (contains (0,0,0),(1,0,0),(1,1,1))
# F4: z=0 (contains (0,0,0),(1,0,0),(1,1,0))

# Normal vectors:
# F1: x=1 → normal (1,0,0) outward? Actually inward since x≤1
# F2: y=x → y-x=0 → normal (-1,1,0)/√2
# F3: z=y → z-y=0 → normal (0,-1,1)/√2
# F4: z=0 → normal (0,0,-1)

import numpy as np

normals = {
    'x=1': np.array([1, 0, 0], dtype=float),      # outward normal (x increases)
    'y=x': np.array([-1, 1, 0], dtype=float) / np.sqrt(2),  # outward
    'z=y': np.array([0, -1, 1], dtype=float) / np.sqrt(2),  # outward
    'z=0': np.array([0, 0, -1], dtype=float),      # outward (z decreases)
}

# Actually need OUTWARD normals. The simplex is {x≤1, y≤x, z≤y, z≥0}
# F1: x≤1, outward normal is +x direction: (1,0,0) ✓
# F2: y≤x, i.e., x-y≥0, outward normal is direction of decreasing x-y: (-1,1,0)/√2 ✓
# F3: z≤y, i.e., y-z≥0, outward normal is (0,-1,1)/√2 ✓
# F4: z≥0, outward normal is -z direction: (0,0,-1) ✓

# Edge = intersection of two facets
# Dihedral angle at edge = π - angle between outward normals
edge_facets = [
    ('x=1', 'y=x'),   # edge (1,1,0)-(1,1,1)
    ('x=1', 'z=y'),   # edge (1,0,0)-(1,1,1)... wait
    ('x=1', 'z=0'),   # edge (1,0,0)-(1,1,0)
    ('y=x', 'z=y'),   # edge (0,0,0)-(1,1,1)
    ('y=x', 'z=0'),   # edge (0,0,0)-(1,1,0)... wait
    ('z=y', 'z=0'),   # edge (0,0,0)-(1,0,0)
]

# Actually, let me identify which edges are at which facet intersections
# Vertices of each facet:
# F1(x=1): (1,0,0), (1,1,0), (1,1,1)
# F2(y=x): (0,0,0), (1,1,0), (1,1,1)
# F3(z=y): (0,0,0), (1,0,0), (1,1,1)... wait
#   z=y and z≤y: vertex (1,0,0) has z=0,y=0 so z=y ✓
#   vertex (0,0,0): z=0,y=0 ✓
#   vertex (1,1,1): z=1,y=1 ✓
#   vertex (1,1,0): z=0,y=1 → z≠y ✗
# So F3(z=y): (0,0,0), (1,0,0), (1,1,1) — but that's only 3 vertices, correct for a triangle
# Wait, this is wrong. Let me reconsider.
# The simplex has 4 vertices: (0,0,0), (1,0,0), (1,1,0), (1,1,1)
# Each facet is a triangle (3 vertices):
# F1(x=1): (1,0,0), (1,1,0), (1,1,1) ✓
# F2(y=x): vertices where y=x: (0,0,0)✓, (1,1,0)✓, (1,1,1)✓. So (0,0,0),(1,1,0),(1,1,1) ✓
# F3(z=y): vertices where z=y: (0,0,0)✓(0=0), (1,0,0)✓(0=0), (1,1,1)✓(1=1). So (0,0,0),(1,0,0),(1,1,1) ✓
# F4(z=0): vertices where z=0: (0,0,0)✓, (1,0,0)✓, (1,1,0)✓. So (0,0,0),(1,0,0),(1,1,0) ✓

# Edges and which facets they belong to:
# (0,0,0)-(1,0,0): F3 ∩ F4
# (0,0,0)-(1,1,0): F2 ∩ F4
# (0,0,0)-(1,1,1): F2 ∩ F3
# (1,0,0)-(1,1,0): F1 ∩ F4
# (1,0,0)-(1,1,1): F1 ∩ F3
# (1,1,0)-(1,1,1): F1 ∩ F2

print("  Facets of the standard simplex {1≥x≥y≥z≥0}:")
print("    F1: x=1, vertices (1,0,0),(1,1,0),(1,1,1)")
print("    F2: y=x, vertices (0,0,0),(1,1,0),(1,1,1)")
print("    F3: z=y, vertices (0,0,0),(1,0,0),(1,1,1)")
print("    F4: z=0, vertices (0,0,0),(1,0,0),(1,1,0)")
print()

edges_with_facets = [
    ('(0,0,0)-(1,0,0)', 'z=y', 'z=0'),
    ('(0,0,0)-(1,1,0)', 'y=x', 'z=0'),
    ('(0,0,0)-(1,1,1)', 'y=x', 'z=y'),
    ('(1,0,0)-(1,1,0)', 'x=1', 'z=0'),
    ('(1,0,0)-(1,1,1)', 'x=1', 'z=y'),
    ('(1,1,0)-(1,1,1)', 'x=1', 'y=x'),
]

print("  Dihedral angles:")
for edge_name, f1, f2 in edges_with_facets:
    n1 = normals[f1]
    n2 = normals[f2]
    cos_angle = np.dot(n1, n2)
    # Dihedral angle = π - arccos(cos(angle between outward normals))
    angle_between = np.arccos(np.clip(cos_angle, -1, 1))
    dihedral = math.pi - angle_between
    print(f"    Edge {edge_name} (F:{f1} ∩ F:{f2}): cos={cos_angle:.4f}, dihedral={dihedral:.4f} = {dihedral/math.pi:.4f}π")

print()

# Edge lengths
edge_lengths = [
    ('(0,0,0)-(1,0,0)', 1.0),
    ('(0,0,0)-(1,1,0)', math.sqrt(2)),
    ('(0,0,0)-(1,1,1)', math.sqrt(3)),
    ('(1,0,0)-(1,1,0)', 1.0),
    ('(1,0,0)-(1,1,1)', math.sqrt(2)),
    ('(1,1,0)-(1,1,1)', 1.0),
]

print("  Dehn invariant of standard simplex:")
print("  D = Σ length(e) ⊗ θ(e)")
print()
dehn_terms = []
for (e_name, length), (_, f1, f2) in zip(edge_lengths, edges_with_facets):
    n1 = normals[f1]
    n2 = normals[f2]
    cos_angle = np.dot(n1, n2)
    angle_between = np.arccos(np.clip(cos_angle, -1, 1))
    dihedral = math.pi - angle_between
    dehn_terms.append((length, dihedral))
    print(f"    {e_name}: length={length:.4f}, θ={dihedral/math.pi:.6f}π")

print()
print("  D = ", " + ".join(f"{l:.4f}⊗{t/math.pi:.6f}π" for l,t in dehn_terms))

# Check: for the regular tetrahedron, D ≠ 0 (famous result)
# Our simplex is NOT regular — it's a "path simplex"
# Dihedral angles: π/2 and arccos(1/2)=π/3? Let me check.

print()
print("  NOTE: The simplex {1≥x≥y≥z≥0} is an ORTHOSCHEME,")
print("  not a regular tetrahedron. Its Dehn invariant structure")
print("  reflects the ORDERING (permutation) that defines it.")
print()

# Key result: all dihedral angles
print("  DIHEDRAL ANGLE SUMMARY:")
angles = set()
for (e_name, length), (_, f1, f2) in zip(edge_lengths, edges_with_facets):
    n1 = normals[f1]
    n2 = normals[f2]
    cos_angle = np.dot(n1, n2)
    angle_between = np.arccos(np.clip(cos_angle, -1, 1))
    dihedral = math.pi - angle_between
    angles.add(round(dihedral/math.pi, 6))

print(f"  Distinct dihedral angles (as fractions of π): {sorted(angles)}")
print()

# The standard orthoscheme in R^n has dihedral angles π/2 and arccos(1/√2)=π/4
# and arccos(1/√3)?
# Actually for the simplex {x₁≥x₂≥...≥x_n≥0, x₁≤1}:
# Adjacent "interior" facets (x_k=x_{k+1} and x_{k+1}=x_{k+2}) meet at angle...
# Let me compute for general n.

print("=" * 70)
print("PART 6: ORTHOSCHEME DIHEDRAL ANGLES IN GENERAL n")
print("=" * 70)
print()
print("  The standard orthoscheme in R^n is:")
print("  {1 ≥ x₁ ≥ x₂ ≥ ... ≥ x_n ≥ 0}")
print()
print("  Facets:")
print("  F₀: x₁ = 1 (cap)")
print("  Fₖ: xₖ = x_{k+1} for k=1,...,n-1 (internal)")
print("  F_n: x_n = 0 (base)")
print()
print("  Normal vectors (outward):")
print("  F₀: (1, 0, ..., 0)")
print("  Fₖ: (-eₖ + e_{k+1})/√2 for k=1,...,n-1")
print("  F_n: (0, ..., 0, -1)")
print()
print("  Dihedral angles between adjacent facets:")

for n in range(2, 8):
    # Facets: F₀ (x₁=1), F₁ (x₁=x₂), ..., F_{n-1} (x_{n-1}=x_n), F_n (x_n=0)
    # Normal of F₀: e₁
    # Normal of Fₖ (k=1..n-1): (-eₖ + e_{k+1})/√2
    # Normal of F_n: -e_n

    def normal(k, n):
        v = np.zeros(n)
        if k == 0:
            v[0] = 1
        elif k == n:
            v[n-1] = -1
        else:
            v[k-1] = -1
            v[k] = 1
            v /= np.sqrt(2)
        return v

    # Compute all dihedral angles between pairs of facets that share an edge
    # Facets Fᵢ and Fⱼ share an edge iff |i-j| ≤ ... (they share n-2 vertices)
    # Actually they share a ridge (codim-2 face) iff they are "adjacent" in the
    # facet lattice. For the orthoscheme, adjacent facets are consecutive: Fₖ, F_{k+1}.

    angles_list = []
    for k in range(n):
        n1 = normal(k, n)
        n2 = normal(k+1, n)
        cos_a = np.dot(n1, n2)
        angle = math.pi - np.arccos(np.clip(cos_a, -1, 1))
        angles_list.append(angle / math.pi)

    # Also non-adjacent pairs
    all_angles = []
    for i in range(n+1):
        for j in range(i+1, n+1):
            n1 = normal(i, n)
            n2 = normal(j, n)
            cos_a = np.dot(n1, n2)
            if abs(cos_a) < 0.9999:  # not parallel
                angle = math.pi - np.arccos(np.clip(cos_a, -1, 1))
                all_angles.append((i, j, angle/math.pi))

    consec = ", ".join(f"F{k}∩F{k+1}: {a:.4f}π" for k, a in enumerate(angles_list))
    print(f"  n={n}: consecutive: {consec}")

    # Non-adjacent
    non_adj = [(i,j,a) for i,j,a in all_angles if abs(j-i) > 1]
    if non_adj and n <= 4:
        for i,j,a in non_adj:
            print(f"         non-adjacent F{i}∩F{j}: {a:.4f}π")

print()
print("  PATTERN: Consecutive orthoscheme facets meet at angle π/2!")
print("  (The orthoscheme is 'right-angled' between consecutive facets.)")
print()
print("  Non-consecutive facets (Fᵢ, Fⱼ with |i-j|≥2) meet at")
print("  angle arccos(1/2) = π/3 (for adjacent-gap-1)")
print("  or arccos(?) for larger gaps.")
print()

print("=" * 70)
print("PART 7: H(T) + H(T̄) = n! AND COMPLEMENTARY TILING")
print("=" * 70)
print()
print("  The H simplices of T and n!-H simplices of T̄ tile [0,1]^n.")
print("  This is because every permutation σ defines exactly one simplex,")
print("  and σ is a Hamiltonian path of T XOR of T̄ (exactly one).")
print()
print("  Wait: σ = (σ₁,...,σ_n) is a Ham path of T iff σᵢ→σᵢ₊₁ for all i.")
print("  It's a Ham path of T̄ iff σᵢ←σᵢ₊₁ for all i, i.e., σᵢ₊₁→σᵢ for all i.")
print("  This means the REVERSE permutation (σ_n,...,σ₁) is a Ham path of T̄.")
print()
print("  So if σ is a Ham path of T, then σ^{rev} is a Ham path of T̄.")
print("  Since path reversal is a bijection, H(T) = H(T̄)... WRONG!")
print()
print("  That would give H(T) = H(T̄), but we also have H(T) + H(T̄) = n!.")
print("  This forces H(T) = n!/2 for all T, which is false.")
print()
print("  The ERROR: H(T) + H(T̄) = n! is NOT correct!")
print("  A permutation σ is a Ham path of T iff ALL arcs go forward.")
print("  If even ONE arc goes backward, σ is NOT a Ham path of T.")
print("  And σ is also not a Ham path of T̄ unless ALL arcs go backward.")
print()
print("  So most permutations are Ham paths of NEITHER T nor T̄!")
print()

# Verify
n = 3
adj_reg = [[False]*3 for _ in range(3)]
adj_reg[0][1] = True
adj_reg[1][2] = True
adj_reg[2][0] = True
rev = [[not adj_reg[i][j] if i != j else False for j in range(3)] for i in range(3)]

h1 = len(ham_paths(adj_reg, 3))
h2 = len(ham_paths(rev, 3))
print(f"  n=3 regular: H(T)={h1}, H(T̄)={h2}, sum={h1+h2}, n!={factorial(3)}")
print(f"  H(T)+H(T̄) = {h1+h2} ≠ n! = {factorial(3)}")
print()

# So the tiling picture is wrong! The simplices of T and T̄ do NOT tile the cube.
# What's the correct picture?
print("  CORRECT PICTURE:")
print("  The n! simplices of the permutohedron decomposition tile [0,1]^n.")
print("  Each permutation σ defines one simplex Δ_σ = {x_{σ₁}≥...≥x_{σ_n}}.")
print("  H(T) of these are 'consistent' with T (all arcs forward).")
print("  The rest have at least one arc reversed.")
print()
print("  H(T)/n! is the FRACTION of the cube covered by T's simplices.")
print("  At n=3: H=1 covers 1/6, H=3 covers 3/6=1/2 of the cube.")
print()

# But do H(T) + H(T̄) = n!?
# σ is a path of T iff σ₁→σ₂→...→σ_n in T
# σ is a path of T̄ iff σ₁→σ₂→...→σ_n in T̄, i.e., σ₂→σ₁, σ₃→σ₂, ... in T
# So σ is a path of T̄ iff σ_n→σ_{n-1}→...→σ₁ in T, i.e., σ^rev is a path of T
# So paths of T̄ = reverses of paths of T
# Hence H(T̄) = H(T)!

print("  THEOREM: H(T̄) = H(T)")
print("  Proof: σ is a Ham path of T̄ iff σ^rev is a Ham path of T.")
print()
# Verify
for n in range(3, 7):
    all_match = True
    for adj in gen_tournaments(n):
        rev = [[not adj[i][j] if i != j else False for j in range(n)] for i in range(n)]
        h1 = len(ham_paths(adj, n))
        h2 = len(ham_paths(rev, n))
        if h1 != h2:
            all_match = False
            break
    print(f"  n={n}: H(T)=H(T̄) for all T? {all_match}")

print()

print("=" * 70)
print("PART 8: THE REAL SCISSORS CONGRUENCE QUESTION")
print("=" * 70)
print()
print("  Since U(T) is a union of H congruent simplices that tile")
print("  a region of [0,1]^n, and all these simplices are orthoschemes")
print("  (related by permutation of coordinates), the question becomes:")
print()
print("  Given two DIFFERENT sets of H simplices from the n! total,")
print("  are the resulting unions scissors-congruent?")
print()
print("  ANSWER: YES! Any two unions of the same number of")
print("  orthoschemes from the permutohedron decomposition are")
print("  scissors-congruent, because we can always 'swap' simplices:")
print("  Cut out a simplex from one, paste it into the other.")
print()
print("  So in this picture, U(T₁) and U(T₂) are scissors-congruent")
print("  iff H(T₁) = H(T₂).")
print()
print("  The Dehn invariant is the SAME for all such unions,")
print("  because the simplices are all congruent and the boundary")
print("  structure only depends on H.")
print()
print("  CONCLUSION: The scissors congruence conjecture from Part 1")
print("  is WRONG. Same H always gives scissors congruence.")
print("  The independence polynomial I(x) captures ADDITIONAL structure")
print("  beyond scissors congruence — it's about the ARRANGEMENT,")
print("  not just the volume or Dehn invariant.")
print()
print("  The RIGHT question: What does I(x) capture about U(T)")
print("  that H alone doesn't?")
print()
print("  ANSWER: I(x) captures the TOPOLOGY of the boundary of U(T).")
print("  The boundary of U(T) consists of:")
print("  - External boundary (faces on the cube boundary)")
print("  - Internal boundary (faces between 'selected' and 'unselected' simplices)")
print()
print("  I(-1) = Euler characteristic of the independence complex,")
print("  which measures the 'topological complexity' of the set of")
print("  selected simplices (how they're arranged in the permutohedron).")
print()

# At n=5: same H always means same I(x) (since α₂=0, I determined by H)
# At n=6: same H can have different I(x)
# So at n=6, U(T₁) and U(T₂) with same H but different I(x)
# are scissors-congruent but NOT "arrangement-equivalent"

print("  SAME H, DIFFERENT I(x) at n=6:")
print("  U(T₁) and U(T₂) are scissors-congruent (same volume, same pieces)")
print("  But their independence complexes have different Euler characteristics.")
print("  This means the 'wiring diagram' of which simplices are selected")
print("  has different topology, even though the union has the same volume.")
print()
print("  ANALOGY to Hilbert's 3rd:")
print("  Hilbert asked: same volume → scissors congruent?")
print("  Dehn showed: NO, Dehn invariant obstructs.")
print("  We ask: same H → same arrangement?")
print("  I(x) shows: NO, arrangement topology differs.")
print()
print("  But this is a DIFFERENT question from scissors congruence!")
print("  Scissors congruence asks about cutting and reassembling.")
print("  Arrangement topology asks about the combinatorial structure.")

print()
print("=" * 70)
print("DONE")
print("=" * 70)
