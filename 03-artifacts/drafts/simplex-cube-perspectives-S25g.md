# Tournaments as Simplices in Cubes: Perspectives, Packing, and Functional Structure

**Author:** kind-pasteur-2026-03-06-S25g
**Purpose:** Creative exploration of tournaments as simplices, perspective counting, and functional interpretation of the w_0 formula

---

## I. Tournaments as Labeled Simplices

### The Simplex View

A tournament on n vertices is a COMPLETE directed graph: every pair of vertices has exactly one arc. The underlying UNDIRECTED graph is K_n, the complete graph, which is the 1-skeleton of the (n-1)-simplex.

- K_3 = triangle = boundary of 2-simplex
- K_4 = complete graph on 4 = boundary of tetrahedron (3-simplex)
- K_5 = boundary of 4-simplex (pentachoron)
- K_n = boundary of (n-1)-simplex

So a TOURNAMENT is a simplex with a BINARY LABEL on each edge: one of {i->j, j->i}.

The space of all tournaments on n vertices is {0,1}^m where m = C(n,2) = number of edges of K_n = number of edges of the (n-1)-simplex.

### The Cube View

The space {0,1}^m is the vertex set of the m-dimensional hypercube. Each tournament is a VERTEX of this hypercube. The hypercube has 2^m = 2^{C(n,2)} vertices.

So: **tournaments live as vertices of a hypercube, where the simplex's edges index the cube's coordinates.**

The simplex has C(n,2) edges. The cube has C(n,2) dimensions.
The simplex has n vertices. The cube has 2^{C(n,2)} vertices.

The SIMPLEX is the INDEX SPACE (which edges to label).
The CUBE is the LABEL SPACE (how to label them).

### Simplex-in-Cube Packing

The user's insight: a simplex "sits inside" a cube in a specific geometric way.

- An equilateral triangle (2-simplex) sits inside a square (2-cube): the triangle has 3 edges, the square has 2 dimensions. The triangle projects into the square with two halves on either side of the diagonal.

- A regular tetrahedron (3-simplex) sits inside a cube (3-cube): the tetrahedron has 6 edges, and the cube has 3 dimensions. Actually, the standard embedding places the tetrahedron's 4 vertices at alternating corners of the cube: (0,0,0), (1,1,0), (1,0,1), (0,1,1). The remaining 4 cube corners form the DUAL tetrahedron.

For TOURNAMENTS: the simplex K_n has m = C(n,2) edges, and the tournament space is the m-cube. The simplex's VERTEX SET has n elements, while the cube's vertex set has 2^m elements. Each tournament assigns a binary value to each simplex edge.

**Key geometric fact:** The tetrahedron-in-cube embedding has the property that the 4 sub-tetrahedra between the embedded tetrahedron and the cube faces have half the volume. Similarly, tournaments that are "half-way" through the cube (Hamming weight m/2) have maximum H — this is the PERPENDICULAR MAXIMIZER phenomenon!

---

## II. The Perspective Counting Function

### Definition

For a tournament T on n vertices, the NUMBER OF DISTINCT PERSPECTIVES is the number of vertex orbits under Aut(T).

- If Aut(T) is trivial: n distinct perspectives (each vertex sees differently)
- If Aut(T) is transitive (VT): 1 perspective (all vertices equivalent)
- In general: # perspectives = n / average orbit size = n * (# orbits under Aut)

### The Observation

At n=3 (2 iso classes):
- Transitive tournament: |Aut| = 1, perspectives = 3
- Cyclic tournament: |Aut| = 3, perspectives = 1
- Total: 3 + 1 = 4

At n=4 (4 iso classes):
- The "paired" classes (non-self-converse) each have perspectives = 2
- The "unpaired" classes (self-converse) each have perspectives = 4
- Total: 2 + 2 + 4 + 4 = 12

At n=5: total = ??? (to be computed)

### The Sequence and the Transition

The user observes: the total perspectives at n might be connected to the transition from n to n+1. Let's denote:
  P(n) = sum over all iso classes T of (# vertex orbits of T)

Then:
  P(3) = 4
  P(4) = 12
  P(5) = ??? (computing)

The ratio P(4)/P(3) = 12/4 = 3 = n-1 for n=4.

If P(n) = (n-1)! * something, we might have P(n) = (n-1) * P(n-1), giving P(n) = (n-1)!/2 * 4 = 2*(n-1)!.

Check: 2*2! = 4 (P(3)), 2*3! = 12 (P(4)). So P(n) = 2*(n-1)!?

If true: P(5) = 2*4! = 48, P(6) = 2*5! = 240.

This would be remarkable: the total number of distinct vertex-perspectives across all tournament isomorphism classes is exactly 2*(n-1)!.

### Why 2*(n-1)!?

If P(n) = 2*(n-1)!, this equals n!/n * 2 = 2*(n-1)!.

By Burnside's lemma: sum over iso classes of (# orbits) = sum over iso classes of (n / |Aut|) ... no, # orbits != n/|Aut| in general.

Actually, by the orbit-counting theorem: # orbits = (1/|G|) * sum_{g in G} |Fix(g)|. For the action of Aut(T) on vertices:

# orbits = (1/|Aut(T)|) * sum_{g in Aut(T)} (# vertices fixed by g)

Summing over ALL iso classes:
P(n) = sum_T (1/|Aut(T)|) * sum_{g in Aut(T)} |Fix(g)|

where the sum is over unlabeled tournaments T. This equals:
P(n) = sum_T sum_{g in Aut(T)} |Fix(g)| / |Aut(T)|

By the orbit-stabilizer theorem and counting labeled tournaments:
sum_T (n!/|Aut(T)|) = 2^{C(n,2)} (total labeled tournaments)

But we need sum_T (# orbits), not sum_T (n!/|Aut(T)|).

Alternative: sum of # orbits = sum_T sum_v [v is representative of its orbit]
= sum_T (# orbits of T)

This is the number of POINTED TOURNAMENT ISOMORPHISM CLASSES: pairs (T, orbit) where T is a tournament and orbit is a vertex orbit.

CONJECTURE: This count equals 2*(n-1)! for all n >= 3.

---

## III. The w_0 Formula as a Function

### The Formula

At n=5: w_0(T) = -t_3(T) + 2*t_5(T) + 1

where t_3 = number of directed 3-cycles, t_5 = number of directed 5-cycles.

### Viewing w_0 as a Function on Tournament Space

w_0: {0,1}^{10} -> Z (since n=5 has 10 edges)

This function on the 10-dimensional hypercube has:
- Linear terms in t_3 and t_5
- But t_3 and t_5 are themselves POLYNOMIAL functions of the edge variables
- t_3 = sum of products of 3 specific edge variables (one per cyclic triple)
- t_5 = sum of products of 5 specific edge variables (one per 5-cycle)

So w_0 is a POLYNOMIAL on {0,1}^10 of degree at most 5.

### Similar Functions

The structure f(x) = -L_3(x) + 2*L_5(x) + 1, where L_k counts k-cycles,
has a natural generalization:

  w_{n-1-2k}(T) = sum_{j=0}^{k} c_{k,j} * t_{2j+1}(T) + const

This is a LINEAR COMBINATION of odd-cycle counts (at n=5 where the only odd cycles are 3-cycles and 5-cycles).

### Connection to the Euler Characteristic

The formula w_0 = -t_3 + 2*t_5 + 1 resembles the Euler characteristic:
  chi = V - E + F - ... (alternating sum of simplex counts)

If we think of odd cycles as "simplices" of different dimensions:
- 3-cycles are "triangular faces" (2-dimensional)
- 5-cycles are "pentagonal faces" (4-dimensional)

Then w_0 = 1 - t_3 + 2*t_5 looks like an ALTERNATING SUM with coefficients.

But the coefficient of t_5 is +2, not +1. This suggests a WEIGHTED Euler characteristic where higher-dimensional cycles carry more weight.

In the OCF H = I(Omega, 2) = sum alpha_k * 2^k:
  The weight 2^k makes higher-dimensional independent sets (more cycles) exponentially more important.

At n=5 (a_2 = 0): H = 1 + 2*(t_3 + t_5).
At r=0: w_0 = 1 - t_3 + 2*t_5.

The DIFFERENCE: H - w_0 = 3*t_3. So w_0 = H - 3*t_3 at n=5.

### The Functional Equation

W(1/2) = H gives: w_0 + w_2/4 + w_4/16 = H.
With w_2 = 12*t_3 - 30 and w_4 = 120:
  w_0 + (12*t_3 - 30)/4 + 120/16 = H
  w_0 + 3*t_3 - 7.5 + 7.5 = H
  w_0 = H - 3*t_3

So w_0(T) = H(T) - 3*t_3(T) for ALL 5-vertex tournaments!

This means: w_0 measures H with a PENALTY of 3 per 3-cycle.
The "signed path count" at r=0 equals the Hamiltonian path count MINUS 3 times the 3-cycle count.

### Geometric Interpretation

At r = 0: each path's contribution is prod(s_e) = prod(A[e] - 1/2).
This is (+1/2)^f * (-1/2)^b where f = forward arcs, b = backward arcs.
= (1/2)^{n-1} * (-1)^b.

So W(0) = (1/2)^{n-1} * sum_P (-1)^{b(P)} = (1/2)^{n-1} * [# even-b paths - # odd-b paths].

W(0) = w_0 = -t_3 + 2*t_5 + 1.

So: (# even-backward paths - # odd-backward paths) / 2^{n-1} = 1 - t_3 + 2*t_5.

This is a SIGNED PARITY count of Hamiltonian paths, signed by the parity of backward arcs!

---

## IV. Simplices, Cubes, and Dimensional Transitions

### The n to n+1 Transition

Adding a new vertex to an n-tournament creates an (n+1)-tournament by choosing
n new arc directions (one to each existing vertex). This is equivalent to
choosing a point in {0,1}^n = the n-cube.

So: the transition from n to n+1 corresponds to EMBEDDING the simplex structure
in one more cube dimension:

  Tournament(n) in {0,1}^{C(n,2)}
  Tournament(n+1) in {0,1}^{C(n+1,2)} = {0,1}^{C(n,2)+n}

The extra n dimensions correspond to the n arcs connecting the new vertex to
all existing vertices.

For each n-tournament T, there are 2^n possible (n+1)-tournaments extending T.
These 2^n extensions correspond to the vertices of an n-cube.

The STRUCTURE of this n-cube of extensions is:
- Opposite corners: the new vertex either dominates all or is dominated by all
- Hamming weight k: the new vertex has out-degree k among existing vertices

### The Score Sequence Constraint

The new vertex's out-degree to existing vertices is the Hamming weight of the
extension vector. The full score sequence of the (n+1)-tournament includes
this new score plus the original scores (each incremented by 0 or 1 depending
on arc direction).

### Packing Simplices in Cubes

The regular simplex on n+1 vertices can be embedded in the n-cube by placing
vertices at alternating cube corners. For n=3:

  4 vertices of tetrahedron at: (0,0,0), (1,1,0), (1,0,1), (0,1,1)
  4 vertices of dual tetra at: (1,0,0), (0,1,0), (0,0,1), (1,1,1)

The tournament interpretation: the 4 vertices of the embedded tetrahedron
correspond to 4 specific score patterns. The "alternating corners" condition
means each pair of selected corners differs in exactly 2 coordinates —
which is the MINIMUM Hamming distance for a tournament-relevant pattern.

For tournaments: two score assignments that differ in 2 coordinates means
flipping 2 arcs. The minimum number of arc flips to change between
"alternating" tournaments is 2 — this is the SWITCHING distance.

### The Tetrahedron-Cube Decomposition

A cube can be decomposed into 6 tetrahedra (simplices). Similarly, the
tournament space {0,1}^m can be decomposed into simplicial complexes.

The TILING MODEL of tournaments already captures this: the pin grid
(staircase Young diagram) is the simplicial structure, and the tiling
(choice of 0 or 1 for each cell) selects a vertex of the hypercube.

---

## V. Functions Similar to w_0 = -t_3 + 2*t_5 + 1

### Euler Characteristics of Simplicial Complexes

chi(Delta) = sum_{k=-1}^{d} (-1)^k f_k where f_k = # k-simplices.

For the conflict graph Omega(T):
  chi(Omega) = 1 - alpha_1 + (# edges of Omega) - (# triangles of Omega) + ...

But w_0 = -t_3 + 2*t_5 + 1 is NOT the Euler characteristic of Omega.

### Möbius Functions on Posets

The Möbius function mu(P) of a poset P satisfies mu(x,y) = -sum_{x<=z<y} mu(x,z).
For the Boolean lattice: mu = (-1)^{rank}.

The transfer matrix M[a,b] = sum_S (-1)^|S| F(S) IS a Möbius inversion.
And w_0 = tr(c_0) is the "constant term" of this inversion.

So w_0 is a VALUE of a Möbius function evaluated at the bottom element.

### Jones Polynomial Evaluations

The Jones polynomial V(t) of a knot evaluates to specific invariants at roots of unity:
  V(1) = 1 (always)
  V(-1) = determinant of knot
  V(i) = Arf invariant

Similarly, W(r) evaluates to:
  W(1/2) = H (Hamiltonian path count)
  W(-1/2) = (-1)^{n-1} H (by symmetry)
  W(0) = w_0 = H - 3*t_3 (signed count)

The evaluation W(0) = H - 3*t_3 is the tournament analogue of the
DETERMINANT of a knot! It measures the "signed complexity" of the tournament.

### Tutte Polynomial Evaluations

The Tutte polynomial T(x,y) evaluates to:
  T(1,1) = # spanning trees
  T(2,0) = # acyclic orientations
  T(0,2) = # totally cyclic orientations
  T(1,2) = # connected spanning subgraphs

For tournaments, the analogue:
  W(1/2) = H = "trees" (Hamiltonian paths)
  W(0) = "signed count" (signed by backward-arc parity)

---

## VI. Computational Results (S25g)

### w_0 = H - 3*t_3: VERIFIED at n=5, FAILS at n=7

**n=5:** Exhaustively verified for ALL 1024 tournaments. The identity w_0 = H - 3*t_3 is EXACT.

**n=7:** The formula w_0 = a*H + b*t_3 + c*t_5 + d does NOT fit:
- Best regression: w0 = 0.053*H - 0.626*t3 - 0.070*t5 + 1.894 (max error 2.54)
- w_0 at n=7 requires ADDITIONAL INVARIANTS beyond (H, t_3, t_5)
- This is consistent with the S25f finding that w_2 at n=7 depends on finer invariants

### Perspective Counting: P(n) = 2*(n-1)! FAILS at n=6

Exhaustive computation:
- P(2) = 2 = 2*1! (match)
- P(3) = 4 = 2*2! (match)
- P(4) = 12 = 2*3! (match)
- P(5) = 48 = 2*4! (match)
- P(6) = 296 ≠ 240 = 2*5! (FAILS)

Orbit distributions:
- n=3: {1:1, 3:1} — 2 classes
- n=4: {2:2, 4:2} — 4 classes
- n=5: {1:1, 3:4, 5:7} — 12 classes
- n=6: {2:5, 4:10, 6:41} — 56 classes

### P(n) = OEIS A093934

The sequence P(n) = 1, 2, 4, 12, 48, 296, 3040, ... matches **A093934** (offset by 1).
OEIS describes this as "number of unlabeled tournaments with n signed nodes."

P(n) = # rooted tournament isomorphism classes on n vertices
     = # (tournament T, distinguished vertex v) pairs up to isomorphism
     = sum over iso classes of (# vertex orbits under Aut(T))

The coincidence P(n) = 2*(n-1)! for n <= 5 breaks because tournament automorphism
groups become richer at n=6 (more classes with non-trivial automorphisms).

---

## VII. Open Questions (Updated)

1. ~~Does P(n) = 2*(n-1)! hold for n >= 5?~~ ANSWERED: No, fails at n=6. P(n) = A093934.
2. ~~Is w_0 = H - 3*t_3 general?~~ ANSWERED: No, fails at n=7. Needs additional invariants.
3. What invariants beyond (t_3, t_5) determine w_0 at n=7? (connects to opus's sigma pattern theory)
4. What is the geometric meaning of the simplex-in-cube embedding for tournament enumeration?
5. Does the "signed backward-parity count" W(0) have independent combinatorial significance?
6. Can the tetrahedron-cube decomposition be used to decompose the tournament space into simplicial pieces with uniform H-distributions?
7. Is there a combinatorial interpretation of why P(n) = 2*(n-1)! for small n?
