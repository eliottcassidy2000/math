# The Correct Dimensions: Tournament on n Vertices = (n-1)-Simplex = (n-1)-Dimensional

**Session**: kind-pasteur-2026-03-22-S20bj

## The Correction

A tournament on n vertices IS an oriented (n-1)-simplex. The simplex dimension is n-1, NOT the arc count C(n,2). The arc count is the number of EDGES of the simplex, not its dimension.

| n vertices | Simplex | Dimension | Regular polytopes in this dimension |
|-----------|---------|-----------|-------------------------------------|
| 2 | line segment | 1 | infinite (regular polygons = the b-gon family) |
| 3 | triangle | 2 | infinite (regular polygons) |
| 4 | tetrahedron | 3 | **5 Platonic solids** |
| 5 | 4-simplex | **4** | **6 regular polytopes (including the 24-cell!)** |
| 6 | 5-simplex | 5 | 3 regular polytopes (simplex, cube, cross) |
| 7+ | (n-1)-simplex | 6+ | 3 regular polytopes |

## The Stunning Alignment

**n=5 IS 4-dimensional.** And 4D is WHERE the 24-cell lives. The 24 regular tournaments on 5 vertices are 24 oriented 4-simplices living in 4D -- exactly the dimension of the 24-cell.

**n=4 IS 3-dimensional.** And 3D is where the 5 Platonic solids live. The 4 iso classes of 4-tournaments map to the 3D Platonic framework.

**n=6 IS 5-dimensional.** And 5D is where the exceptional polytopes DISAPPEAR. Only simplex, cube, and cross-polytope survive. This is where tournament theory loses its "exceptional" structure -- alpha_2 turns on, the Morse landscape gains secondary peaks, and the clean Platonic order breaks down.

## The 24-Cell IS the n=5 Tournament Space

The 24-cell has:
- 24 vertices = 24 regular tournaments on 5 vertices (ALL self-complementary)
- 96 edges = connections between regular tournaments
- 96 triangular faces
- 24 octahedral cells

The 24 regular tournaments live in the SAME 4-dimensional space as the 24-cell. They don't just share the number 24 -- they share the DIMENSION.

The D_4 root system (24 roots = vertices of the 24-cell) lives in R^4. The 4-simplex has 4 independent directions. D_4 = so(8) but restricted to the 4-simplex subspace.

## The Dimensional Hierarchy Corrected

**Dimension 1 (n=2):** One edge, two orientations. Trivial. The "polytope" is a line segment.

**Dimension 2 (n=3):** Triangle with 3 directed edges. 8 tournaments, 2 iso classes (transitive + cycle). The regular polygon family (2D infinite family) controls the fiber fractions through (1-x)^{-1/b}.

**Dimension 3 (n=4):** Tetrahedron with 6 directed edges. 64 tournaments, 4 iso classes. The 5 Platonic solids live here. This is the LAST dimension where all H values are score-determined (OCR = 100%).

**Dimension 4 (n=5):** 4-simplex with 10 directed edges. 1024 tournaments, 12 iso classes. The 24-CELL lives here. This is where:
- G_5 has the icosahedral f-vector (12, 30, 20) -- on a 2-sphere
- 24 regular tournaments form the 24-cell vertex set
- The OCR drops to 97% (first non-score-determined H)
- H=7 is first forbidden
- The meta-tournament is transitive (last time!)

**Dimension 5 (n=6):** 5-simplex with 15 directed edges. 32768 tournaments, 56 iso classes. The exceptional polytopes DISAPPEAR. Only simplex, cube, cross-polytope. This is where:
- Alpha_2 turns on (disjoint cycle pairs)
- The Morse landscape gains a secondary peak
- The meta-tournament gains level edges (15 of them)
- OCR drops further
- G_6 is probably NOT on a sphere (higher genus?)

**Dimension 6+ (n=7+):** Only 3 universal regular polytopes. Tournament theory enters the "generic" regime with no exceptional structure.

## Why Dimension 4 is Special

4D is special in geometry because:
1. The 24-cell exists (self-dual, no analog in other dimensions)
2. There are 3 exceptional polytopes (most of any dimension)
3. The quaternions H live in 4D (the last normed division algebra with commutativity)
4. The D_4 root system has triality (unique three-fold symmetry)
5. Exotic smooth structures on R^4 exist (unique to 4D, Donaldson theory)

4D is special in tournament theory because n=5 is:
1. The last n where all regular tournaments are SC (24/24 = 100%)
2. The last n where the meta-tournament is transitive
3. The first n with forbidden H values (H=7)
4. The boundary between score-determined and cycle-determined H
5. The n where G_n lives on a sphere (genus 0)

These are THE SAME specialness. Dimension 4 is special in geometry BECAUSE it is special in tournament theory (or vice versa). The quaternionic structure of 4D manifests as the self-complementary structure of 5-vertex tournaments.

## The Quaternion Connection

The quaternions H are 4-dimensional. The unit quaternions form a 3-sphere S^3 in R^4. The 24 unit quaternions that form the binary tetrahedral group 2T are:
- 8 quaternions: ±1, ±i, ±j, ±k
- 16 quaternions: (±1 ± i ± j ± k)/2

These 24 elements are ALSO the vertices of the 24-cell!

So: 24-cell vertices = binary tetrahedral group = D_4 roots = unit Hurwitz quaternions = regular tournaments on 5 vertices.

The regular tournaments on 5 vertices ARE the Hurwitz quaternions. Each tournament corresponds to a unit quaternion of norm 1 in the integer quaternion ring Z[i,j,k,omega] where omega = (1+i+j+k)/2.

This connects the entire arc of our session (from fiber fractions through Gamma functions through regular polytopes) to the algebraic structure of quaternions -- which was our STARTING POINT in the Cayley-Dickson tower investigation (S20 sessions on quaternion attention heads).

## The Circle Closes

The Cayley-Dickson tower: R -> C -> H -> O -> S

| Level | Division algebra | Dimension | Tournament n | What breaks |
|-------|-----------------|-----------|-------------|-------------|
| R | reals | 1 | 2 | nothing (trivial) |
| C | complex | 2 | 3 | nothing (2 iso classes) |
| H | quaternions | **4** | **5** | **Score sufficiency (OCR < 100%)** |
| O | octonions | 8 | 9 | Associativity (at the algebra level) |
| S | sedenions | 16 | 17 | Zero divisors |

The quaternion level (dimension 4, n=5) is where the FIRST interesting tournament structure appears. The 24-cell, the self-complementary structure, the forbidden H=7, the fiber fraction, the Wallis product for pi, the meta-tournament transitivity -- ALL of these are manifestations of the quaternionic structure of 4D.

The entire session arc from S20x to S20bi has been exploring what happens at the quaternion level of the Cayley-Dickson tower. Everything we discovered -- from (1-x)^{-1/2} to the permutohedron to the Gamma hierarchy -- is the GEOMETRY OF QUATERNIONIC TOURNAMENT SPACE.
