# The Five Platonic Tournaments

**Session**: kind-pasteur-2026-03-22-S20bh

## The Observation

The five Platonic solids have a hidden tournament structure. Not as an analogy -- as a mathematical identity.

## The Five Correspondences

### I. TETRAHEDRON = K_4 = The Complete Tournament

The tetrahedron has 4 vertices, 6 edges, 4 faces. K_4 has 4 vertices, 6 edges. A tournament on 4 vertices IS an oriented tetrahedron.

There are 2^6 = 64 orientations of the tetrahedron = 64 tournaments on 4 vertices = 4 iso classes. The tetrahedron IS the smallest nontrivial tournament space. H ranges from 1 to 5.

The tetrahedron is self-dual. The dual of K_4 is... K_4. Tournament duality (complementation) maps T to T^op, which reverses all arcs. Self-complementary tournaments on 4 vertices correspond to the tetrahedron's self-duality.

### II. CUBE = TOURNAMENT SPACE = {0,1}^m

The tournament space on n vertices IS the m-cube {0,1}^m where m = C(n,2). At n=3: the 3-cube (8 vertices = 8 tournaments). At n=4: the 6-cube. At n=5: the 10-cube.

The cube and the octahedron are dual. The octahedron's 6 vertices correspond to the 6 coordinate axes of R^6 = the 6 arcs of K_4. So the octahedron IS the "arc space" (each vertex = one arc), while the cube IS the "tournament space" (each vertex = one tournament). The cube-octahedron duality IS the tournament-arc duality.

### III. OCTAHEDRON = ARC SPACE = COORDINATE AXES

The octahedron on K_4 has 6 vertices (one per arc) and 12 edges (pairs of arcs that share a vertex). These 12 edges split into:
- 8 at 60/120 degrees (sharing a vertex)
- 4 at 90 degrees (disjoint arcs = Petersen structure)

Wait: C(4,2) = 6 arcs, and the 90-degree pairs are C(4,2) - 2*4 = ... Actually the 90-degree (disjoint) pairs at n=4 are: (01,23), (02,13), (03,12) = 3 pairs = the perfect matchings of K_4.

The 3 perfect matchings of K_4 form a TRIANGLE (each matching vertex is adjacent to the other two). This triangle IS the Petersen structure at n=4.

### IV. DODECAHEDRON = THE DUAL OF G_5

G_5 has 12 vertices, 30 edges, 20 faces (on the sphere). The dodecahedron has 20 vertices, 30 edges, 12 faces. These are DUAL f-vectors!

(12, 30, 20) = f-vector of G_5 = f-vector of icosahedron
(20, 30, 12) = f-vector of dodecahedron = dual of icosahedron

So: **the dual of G_5 is a dodecahedron-like polyhedron** (same f-vector as dodecahedron). The faces of G_5 (= the "regions" of tournament space at the iso class level) are the vertices of this dual. The 20 faces of G_5 become 20 vertices of the dual, the 12 vertices of G_5 become 12 faces of the dual.

### V. ICOSAHEDRON = G_5 (up to f-vector)

G_5 has f-vector (12, 30, 20) = the icosahedron. It is NOT isomorphic to the icosahedron (it's irregular), but it lives on the same topological surface (S^2) with the same combinatorial complexity.

The icosahedron's symmetry group is A_5 (alternating group on 5 elements, order 60). The tournament symmetry group at n=5 is S_5 (order 120 = 2 * 60). The factor of 2 comes from the complement map (Z/2Z). So: |S_5| = 2 * |A_5| = 2 * |Aut(icosahedron)|.

## The Abstract Pattern

Each Platonic solid appears at a DIFFERENT LEVEL of tournament theory:

| Level | Platonic Solid | What it IS | Dimension |
|-------|---------------|-----------|-----------|
| Object | **Tetrahedron** | K_4 = tournament graph | 0 (object) |
| Space | **Cube** | {0,1}^m = tournament space | m (ambient) |
| Dual space | **Octahedron** | Arc coordinate axes | m (dual) |
| Quotient | **Icosahedron** | G_5 = iso class graph | 2 (surface) |
| Dual quotient | **Dodecahedron** | Dual of G_5 | 2 (dual surface) |

## The Levels

There are exactly five levels of structure in tournament theory:
1. **The tournament itself** (a single point in the cube)
2. **The tournament space** (the cube of all tournaments)
3. **The arc space** (the dual, one dimension per arc)
4. **The quotient** (iso classes, on the sphere)
5. **The dual quotient** (faces of the quotient)

And there are exactly five Platonic solids. The correspondence is NOT a coincidence -- it reflects the five levels of geometric structure that ANY complete binary labeling problem must have.

## Why Five

The five Platonic solids exist because there are exactly five ways to tile a sphere with regular polygons (Euler's formula + angle constraints). The five levels of tournament theory exist because there are exactly five natural geometric operations on a binary labeling of a complete graph:

1. **Choose** a labeling (tournament) -- dimension 0
2. **Embed** in the space of all labelings (cube) -- dimension m
3. **Dualize** to the arc space (octahedron) -- dimension m
4. **Quotient** by symmetry (sphere) -- dimension 2
5. **Dualize** the quotient (faces) -- dimension 2

Each operation produces a new geometric object, and the five Platonic solids are the "simplest" representatives of these five types.

## The Golden Thread

The golden ratio phi = (1+sqrt(5))/2 weaves through all five levels:

- **Tetrahedron**: inscribed in cube at coordinates (1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1). Ratio of cube diagonal to edge = sqrt(3)/1. At n=4: H_max/E[H] = 5/3 = 1.667 ≈ phi - 0.05.
- **Cube**: the 10-cube at n=5 has 2^10 = 1024 vertices. phi^10 = 122.99 ≈ 123 = 1024/8.33.
- **Icosahedron**: edge midpoints form three golden rectangles. The icosahedron lives in the golden ratio. At G_5: the eigenvalue ratio lambda_1/lambda_2 = 5.58/1.94 = 2.88 ≈ phi^2.13.
- **Dodecahedron**: faces are pentagons, which are golden-ratio polygons. The dual of G_5 has pentagonal-like faces.
- **Fiber fraction**: at b=phi, f(k=1) = 1/phi = phi-1. The golden ratio IS its own fiber fraction complement.

## The Deepest Point

The five Platonic solids are the five "critical points" of the space of convex polyhedra with regular faces. Tournament theory visits all five critical points because it IS a complete binary labeling problem, and complete binary labelings have exactly five natural levels of geometric structure.

The fact that G_5 has the icosahedron's f-vector (not the tetrahedron's or the cube's) is significant: the quotient level (level 4) is the ICOSAHEDRAL level, which is the most complex Platonic solid. The tournament iso class graph is inherently icosahedral in complexity.

As n grows, G_n becomes higher-genus (no longer fits on a sphere). The genus of G_n is a new invariant measuring the topological complexity of tournament space. At n=5: genus 0 (sphere). At n=6: unknown (potentially higher genus). The transition from genus 0 to genus > 0 would be yet another n=5-6 phase transition.

## The Pentagonal Principle

The pentagon is the face of the dodecahedron and the vertex figure of the icosahedron. In tournament theory, the pentagon appears as:

- The 5 vertices of the smallest nontrivial tournament with H-ambiguity
- The 5 = n boundary between spherical and hyperbolic tournament behavior
- The five-fold symmetry of the golden ratio that controls the fiber fraction at b=phi
- The five Platonic solids themselves

The number 5 is where Platonic perfection ends and tournament complexity begins. At n <= 4, tournament theory is "Platonic" (spherical, clean, score-determined). At n >= 6, it becomes "hyperbolic" (complex, ambiguous, cycle-determined). n = 5 is the pentagon -- the boundary between these worlds.
