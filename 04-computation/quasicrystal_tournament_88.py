#!/usr/bin/env python3
"""
quasicrystal_tournament_88.py — opus-2026-03-14-S88

The quasicrystal-tournament bridge.

The Fibonacci tiling (LSLLS...) with L=3, S=2 is a 1D quasicrystal.
Its 2D analog is the Penrose tiling, which has:
  - Icosahedral symmetry group (A₅ ≅ PSL(2,4) ≅ PSL(2,5))
  - Golden ratio scaling
  - 5-fold rotational symmetry

Tournament connections:
  - A₅ ≅ PSL(2,4) acts on PG(1,4) = 5 points
  - A₅ = rotation group of icosahedron
  - |A₅| = 60 = |PGL(2,4)|
  - The icosahedron has φ (golden ratio) proportions
  - 5 = number of vertices in smallest non-trivial regular tournament

Also: PSL(2,4) ≅ PSL(2,5) ≅ A₅ is the UNIQUE simple group of order 60.
This "exceptional isomorphism" connects:
  - F₄ (4-element field, Baer subplanes)
  - F₅ (5-element field, regular pentagons)
  - A₅ (alternating group, icosahedral symmetry)
"""

from math import comb, factorial, gcd, sqrt
from collections import Counter

# ══════════════════════════════════════════════════════════════════
# PART 1: The exceptional isomorphism PSL(2,4) ≅ PSL(2,5) ≅ A₅
# ══════════════════════════════════════════════════════════════════

print("=" * 70)
print("PART 1: THE EXCEPTIONAL ISOMORPHISM")
print("=" * 70)

# PSL(2,q) for small q:
# PSL(2,2) ≅ S₃ (order 6)
# PSL(2,3) ≅ A₄ (order 12)
# PSL(2,4) ≅ A₅ (order 60) ← exceptional!
# PSL(2,5) ≅ A₅ (order 60) ← exceptional!
# PSL(2,7) ≅ GL(3,2) (order 168) ← our Fano automorphism group!

for q in [2, 3, 4, 5, 7, 8, 9, 11, 13]:
    order = q * (q*q - 1) // gcd(2, q-1)
    print(f"  |PSL(2, F_{q})| = {order}")

print(f"""
EXCEPTIONAL ISOMORPHISMS:
  PSL(2,2) = S_3             (order 6)
  PSL(2,3) = A_4             (order 12)
  PSL(2,4) = PSL(2,5) = A_5  (order 60)  ← THE KEY ONE
  PSL(2,7) = GL(3,2)         (order 168)

WHY PSL(2,4) = PSL(2,5) = A_5:
  PSL(2,4) acts on PG(1,4) = 5 points (projective line over F_4).
  PSL(2,5) acts on PG(1,5) = 6 points.
  A_5 acts on 5 elements.
  All three have order 60 = 2^2 * 3 * 5.

  The isomorphism PSL(2,4) -> A_5 comes from the action on 5 points
  of PG(1,4) = F_4 union {{infinity}} = {{0, 1, alpha, alpha+1, infinity}}.

  This connects F_4 (the field behind Baer subplanes) to A_5
  (the icosahedral group behind the golden ratio)!
""")

# ══════════════════════════════════════════════════════════════════
# PART 2: The icosahedral structure
# ══════════════════════════════════════════════════════════════════

print("=" * 70)
print("PART 2: ICOSAHEDRAL STRUCTURE AND TOURNAMENTS")
print("=" * 70)

phi = (1 + sqrt(5)) / 2

print(f"""
THE ICOSAHEDRON:
  12 vertices, 30 edges, 20 faces (triangular)
  Rotation group: A_5 (order 60)
  Full symmetry: A_5 x Z_2 (order 120 = |S_5|)

  Vertex coordinates (up to scaling):
  (0, +/-1, +/-phi), cyclic permutations

  phi = {phi:.6f} (golden ratio)

  Key numbers:
  12 vertices = C(4,2) = number of edges in K_4
  30 edges = C(6,2) - 6 = 15 - ... no
  30 = 2 * 15 = 2 * C(6,2) / ...
  20 faces = C(6,3) = 20 3-subsets of 6 elements

  Wait: 20 = C(6,3)? Yes! 20 = 6!/(3!*3!)
  And 20 triangular faces of icosahedron.
  This connects to the 20 3-cycles in some tournament?

  Actually for n=5: number of directed 3-cycles in a regular
  tournament on 5 vertices:
  Total 3-subsets: C(5,3) = 10
  Each gives 2 orientations (CW, CCW).
  In any tournament, each triple has exactly 1 of 2 cycle orientations
  OR is a path (if 0 or 1 of the 2 3-cycles exist).

  For a regular tournament on 5 (score sequence 2,2,2,2,2):
  alpha_1 = n(n^2-1)/24 = 5*24/24 = 5
  So 5 directed 3-cycles.
  These 5 correspond to... 5 faces of a polytope?
  5 = number of vertices of icosahedron / ... no.

ICOSAHEDRAL TOURNAMENT:
  The icosahedron graph has 12 vertices, 30 edges.
  An orientation of the icosahedron graph gives a partial tournament.
  But it's not a complete tournament (only 30 of C(12,2)=66 edges).

  However: the DUAL of the icosahedron is the DODECAHEDRON:
  20 vertices, 30 edges, 12 faces (pentagonal)
  And the Petersen graph IS a subgraph of the dodecahedron graph!

  Petersen graph: 10 vertices, 15 edges, girth 5
  = K(5,2): vertices = 2-element subsets of {{1,...,5}}
  = the Kneser graph

  10 = Petersen vertices = BIBD points at n=6
  15 = Petersen edges = C(6,2) = BIBD blocks at n=6 = duads of S_6
""")

# ══════════════════════════════════════════════════════════════════
# PART 3: Regular tournaments on 5 vertices
# ══════════════════════════════════════════════════════════════════

print("=" * 70)
print("PART 3: REGULAR TOURNAMENTS ON 5 VERTICES")
print("=" * 70)

# Build all tournaments on 5 vertices
import numpy as np

n = 5
arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
n_arcs = len(arcs)  # C(5,2) = 10

regular_count = 0
regular_examples = []

for mask in range(2**n_arcs):
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(arcs):
        if mask & (1 << k):
            adj[i][j] = 1
        else:
            adj[j][i] = 1

    # Check regularity (all scores = 2)
    scores = [sum(adj[i]) for i in range(n)]
    if scores == [2]*n:
        # Count 3-cycles
        cycles = 0
        for a in range(n):
            for b in range(a+1, n):
                for c in range(b+1, n):
                    if adj[a][b] and adj[b][c] and adj[c][a]:
                        cycles += 1
                    elif adj[a][c] and adj[c][b] and adj[b][a]:
                        cycles += 1

        # Compute H
        from itertools import combinations as combs
        odd_cycles = []
        for size in [3, 5]:
            for subset in combs(range(n), size):
                sub_adj = [[adj[i][j] for j in subset] for i in subset]
                s = len(subset)
                # Check if subset forms a directed cycle
                # A cycle visits each vertex exactly once
                if size == 3:
                    a, b, c = 0, 1, 2
                    if sub_adj[a][b] and sub_adj[b][c] and sub_adj[c][a]:
                        odd_cycles.append(frozenset(subset))
                    elif sub_adj[a][c] and sub_adj[c][b] and sub_adj[b][a]:
                        odd_cycles.append(frozenset(subset))
                elif size == 5:
                    # Check all Hamiltonian cycles of K_5
                    from itertools import permutations
                    found = False
                    for perm in permutations(range(5)):
                        is_cycle = True
                        for k in range(5):
                            if not sub_adj[perm[k]][perm[(k+1) % 5]]:
                                is_cycle = False
                                break
                        if is_cycle:
                            found = True
                            break
                    if found:
                        odd_cycles.append(frozenset(subset))

        # Build conflict graph
        oc_list = list(set(odd_cycles))
        nc = len(oc_list)

        # Independence polynomial at x=2
        # I(Omega, 2) = sum_{independent sets S} 2^|S|
        # For small nc, enumerate all subsets
        H = 0
        for sub_mask in range(2**nc):
            is_indep = True
            selected = [i for i in range(nc) if sub_mask & (1 << i)]
            for i in range(len(selected)):
                for j in range(i+1, len(selected)):
                    if oc_list[selected[i]] & oc_list[selected[j]]:
                        is_indep = False
                        break
                if not is_indep:
                    break
            if is_indep:
                H += 2**len(selected)

        regular_count += 1
        regular_examples.append({
            'mask': mask,
            'adj': [row[:] for row in adj],
            '3cycles': cycles,
            'odd_cycles': nc,
            'H': H
        })

print(f"Regular tournaments on n=5: {regular_count}")
print(f"Distinct H values: {sorted(set(t['H'] for t in regular_examples))}")

# Group by H
h_groups = Counter(t['H'] for t in regular_examples)
for h, count in sorted(h_groups.items()):
    example = next(t for t in regular_examples if t['H'] == h)
    print(f"  H={h}: {count} tournaments, 3-cycles={example['3cycles']}, odd_cycles={example['odd_cycles']}")

# ══════════════════════════════════════════════════════════════════
# PART 4: A₅ and the 5 Platonic solids
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 4: THE PLATONIC-TOURNAMENT CORRESPONDENCE")
print("=" * 70)

print(f"""
THE 5 PLATONIC SOLIDS AND TOURNAMENT NUMBERS:

  Solid        V   E   F  |Aut|  f-vector    Tournament link
  Tetrahedron  4   6   4    24   (4,6,4)     |S_4|=24=Golay length
  Cube         8  12   6    48   (8,12,6)    8=T(3), 12=C(4,3)*3
  Octahedron   6  12   8    48   (6,12,8)    6=LCM(2,3)=period
  Dodecahedron 20  30  12   120  (20,30,12)  120=|S_5|=|Aut(Petersen)|
  Icosahedron  12  30  20   120  (12,30,20)  12=2*6, 30=C(6,2)*2

  V - E + F = 2 (Euler's formula) for all five.

  THE EXCEPTIONAL SOLID: The icosahedron has A_5 rotation symmetry.
  A_5 = PSL(2,4) = PSL(2,5) = the exceptional isomorphism.

  THE DUAL PAIR: Icosahedron (12V,20F) <-> Dodecahedron (20V,12F)
  This duality swaps vertices and faces.
  20 = C(6,3): face count of icosahedron = 3-subsets of 6
  12 = 2*6: vertex count = twice the period

  THE CUBE-OCTAHEDRON PAIR: (8V,6F) <-> (6V,8F)
  8 = number of tournaments at n=3
  6 = period of tournament parity

  THE SELF-DUAL: Tetrahedron (4V,4F) = K_4
  4 = |F_4| = the field
  |Aut(Tet)| = 24 = Golay code length!

ICOSAHEDRAL NUMBERS:
  12 vertices: 12 = 2 * 6 = 2 * period
  30 edges: 30 = 2 * 15 = 2 * C(6,2) = 2 * (duads of S_6)
  20 faces: 20 = C(6,3) = 3-subsets of 6
  60 rotations = |A_5| = |PSL(2,4)|
  120 symmetries = |S_5| = |Aut(Petersen)|

  The icosahedron's face structure C(6,3) connects it to
  the 3-SUBSETS of a 6-element set — exactly the context
  where tournament 3-cycles live at n=6!
""")

# ══════════════════════════════════════════════════════════════════
# PART 5: The golden ratio in tournament theory
# ══════════════════════════════════════════════════════════════════

print("=" * 70)
print("PART 5: GOLDEN RATIO APPEARANCES IN TOURNAMENTS")
print("=" * 70)

print(f"""
WHERE φ = {phi:.6f} APPEARS:

1. FIBONACCI WORD: The tournament word L=3, S=2 with L/S ratio -> phi

2. PROJECTIVE PLANE RATIO: |PG(2,F_3)|/|PG(2,F_4)| = 13/21 -> 1/phi

3. CYCLOTOMIC: Phi_6(phi) = phi^2 - phi + 1 = 2 (the generator!)

4. ICOSAHEDRAL COORDINATES: vertices at (0, +/-1, +/-phi)
   and cyclic permutations

5. CONTINUED FRACTION: phi = 1 + 1/(1 + 1/(1 + ...)) = [1;1,1,1,...]
   The SIMPLEST continued fraction = the "most irrational" number.
   Tournament parity is fundamentally about the interplay of 2 and 3,
   and phi arises from the continued fraction [1;1,1,...] where
   each step alternates between adding a 2-tile and a 3-tile.

6. EIGENVALUE CONNECTION: The Fibonacci recurrence matrix
   [[1,1],[1,0]] has eigenvalues phi and psi = -1/phi.
   The tournament parity matrix M(-1) = [[1,-1],[1,0]] has
   eigenvalues that are 6th roots of unity.
   M(2) = [[1,2],[1,0]] has eigenvalues phi*sqrt(2)... no.
   Actually [[1,2],[1,0]] has char poly x^2 - x - 2 = (x-2)(x+1).
   Eigenvalues: 2 and -1.
   Hmm, these are the tournament generators, not phi!

   So: M(2) has eigenvalues 2 (arc generator) and -1 (sign flip).
   M(-1) has eigenvalues that are 6th roots of unity (period 6).
   The Fibonacci matrix M(1) has eigenvalues phi and psi.

   These three matrices share the same form [[1,a],[1,0]]:
   a = 1 -> Fibonacci (golden ratio)
   a = 2 -> Tournament generators (2 and -1)
   a = -1 -> Tournament parity (period 6)

   THE MATRIX FAMILY [[1,a],[1,0]]:
   Eigenvalues = (1 +/- sqrt(1+4a)) / 2
   a = 1: (1+/-sqrt(5))/2 = phi, psi
   a = 2: (1+/-3)/2 = 2, -1
   a = -1: (1+/-sqrt(-3))/2 = primitive 6th roots of unity

   All three are governed by discriminant D = 1 + 4a:
   a=1: D=5 (sqrt(5) -> golden ratio)
   a=2: D=9=3^2 (rational eigenvalues)
   a=-1: D=-3 (complex eigenvalues -> periodicity)
""")

# Verify eigenvalues
for a, name in [(1, "Fibonacci"), (2, "Tournament gen"), (-1, "Tournament parity")]:
    M = np.array([[1, a], [1, 0]], dtype=complex)
    evals = np.linalg.eigvals(M)
    D = 1 + 4*a
    print(f"  M({a:+d}) [{name:20s}]: eigenvalues = {evals[0]:.4f}, {evals[1]:.4f}  (D={D})")

# ══════════════════════════════════════════════════════════════════
# PART 6: The Penrose-Tournament connection
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 6: PENROSE TILINGS AND TOURNAMENT STRUCTURE")
print("=" * 70)

print(f"""
PENROSE TILINGS:
  Use two tiles: thick rhombus (angle 2pi/5) and thin rhombus (angle pi/5)
  Tile ratio -> phi (golden ratio)
  5-fold rotational symmetry

  The de Bruijn construction: a Penrose tiling is the projection of a
  5D lattice Z^5 through a 2D plane with irrational slope.

  The 1D analog: Fibonacci word = projection of Z^2 through a line
  with slope 1/phi.

TOURNAMENT ANALOG:
  A tournament on n vertices = a complete directed graph.
  The SCORE SEQUENCE is the degree vector.
  Regular tournaments have constant score (n-1)/2.

  The regular tournament lattice:
  At n=5: 24 regular tournaments (on labeled vertices)
  At n=7: 2640 regular tournaments

  QUASICRYSTALLINE TOURNAMENTS?
  A quasicrystalline tournament would have:
  - Long-range order (structured H-spectrum)
  - No translational symmetry (non-circulant)
  - Self-similar scaling (Fibonacci-like structure)

  The QR_7 tournament IS crystalline (circulant, periodic).
  But the AP_7 has less symmetry (|Aut|=7 vs 21).
  The "Middle" class has |Aut|=3, the least symmetric.

  Conjecture: as n grows, "quasicrystalline" tournaments
  (those with self-similar but non-periodic structure)
  might dominate the H-spectrum.

THE DEEP CONNECTION:
  Penrose tilings are governed by the golden ratio phi.
  Tournament tilings are governed by the generators 2 and 3.
  The bridge: phi = (1+sqrt(5))/2, and sqrt(5) connects to
  the discriminant of x^2-x-1.

  But 5 = the smallest number of vertices for a non-trivial
  regular tournament! (n must be odd, n=3 is trivial: all isomorphic.)

  So: phi (golden ratio) ↔ n=5 (smallest non-trivial tournament)
  And A_5 = PSL(2,4) = PSL(2,5) connects:
  - The icosahedral symmetry (phi)
  - The F_4 structure (Baer subplanes)
  - The tournament on 5 vertices
  All through one exceptional isomorphism!
""")

# ══════════════════════════════════════════════════════════════════
# PART 7: Summary
# ══════════════════════════════════════════════════════════════════

print("=" * 70)
print("SUMMARY: THE QUASICRYSTAL-TOURNAMENT DICTIONARY")
print("=" * 70)

print(f"""
  QUASICRYSTAL CONCEPT    TOURNAMENT ANALOG
  ─────────────────────   ──────────────────
  Two tile types (L,S)    Two generators (3,2)
  Golden ratio phi        Period 6 / Fibonacci F_8=21
  Penrose rhombi          Odd cycle types (3,5,7,...)
  5-fold symmetry         A_5 = PSL(2,4) = PSL(2,5)
  Z^5 lattice             Score sequence space
  Irrational projection   Cyclotomic evaluation (Phi_3(2)=7)
  Quasiperiodicity        Non-circulant tournaments
  Inflation rule L->LS    H recursion via OCF
  Diffraction pattern     Fourier spectrum of H(T)

  THE MATRIX FAMILY M(a) = [[1,a],[1,0]]:
  a = -1: Period 6 (tournament parity, 6th roots)     D = -3
  a =  1: Golden ratio (Fibonacci, quasicrystals)      D = 5
  a =  2: Generators (tournament arc/cycle counts)      D = 9

  The discriminants: -3, 5, 9
  Sum: -3 + 5 + 9 = 11 = |QR_11| (next Paley prime!)
  Product: -3 * 5 * 9 = -135 = -5 * 27
""")
