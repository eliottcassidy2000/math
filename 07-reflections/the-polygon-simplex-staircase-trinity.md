# The Polygon-Simplex-Staircase Trinity

**Session**: kind-pasteur-2026-03-22-S20be

## The Three Objects

A tournament on n vertices is simultaneously:

1. **A binary labeling of the 1-skeleton of the (n-1)-simplex** (each edge of K_n gets a direction)
2. **A binary tiling of the staircase Young diagram** delta_{n-2} (each cell is 0 or 1)
3. **A point inside the permutohedron** (which is the zonotope dual to the A_{n-1} root system, which IS the regular polygon/simplex structure)

These three descriptions are not independent -- they are different faces of the same geometric object. The connections run through the Gamma function hierarchy Gamma(1/b)^b.

## The Triangle of Equivalences

```
         SIMPLEX (K_n)
        /            \
       /              \
  tournament arc = edge = root alpha_{ij} of A_{n-1}
     /                    \
    /                      \
STAIRCASE (delta_{n-2}) === PERMUTOHEDRON (Pi_n)
   |                              |
   tile (i,j) with range d    score class = face
   |                              |
   H = 1 + 2^d (single tile)  score = projection to base
```

## The Simplex Face

A tournament on n vertices IS a directed graph on K_n. The complete graph K_n is the 1-skeleton of the (n-1)-simplex. The k-faces of the simplex are the k-subsets of vertices. The tournament structure assigns a DIRECTION to every edge (1-face).

The key invariants map to simplex faces:
- **0-faces** (vertices): n tournament players
- **1-faces** (edges): C(n,2) arcs = tournament arcs
- **2-faces** (triangles): C(n,3) potential 3-cycles
- **k-faces**: C(n,k+1) potential (k+1)-cycles

The Coxeter angle structure on the A_{n-1} root system gives the arc interactions:
- 60 degrees: cooperative arcs (share a vertex in the same role)
- 90 degrees: independent arcs (disjoint = Petersen/Kneser structure)
- 120 degrees: conflicting arcs (share a vertex in opposite roles)

## The Staircase Face

The staircase Young diagram delta_{n-2} = (n-2, n-3, ..., 2, 1) encodes the tournament as a binary tiling. Each cell (i,j) of the staircase corresponds to the non-base-path arc connecting vertices i and j+i+1.

The "range" d = j (the column index in the staircase) determines the information content: tile at range d carries 2^d bits. This is the H = 1 + 2^d formula.

The staircase rows have lengths n-2, n-3, ..., 1 -- this is the triangle number T_{n-2}. The total number of cells is C(n-1,2) = the cycle space dimension.

## The Permutohedron Face

The permutohedron Pi_n is the convex hull of all permutations of (1,2,...,n), embedded in R^n. It is a zonotope -- the Minkowski sum of the line segments corresponding to the positive roots of A_{n-1}.

The score map sends a tournament to a point in Pi_n (the score sequence). The fibers of this map are the score classes. The permutohedron face lattice encodes the possible score sequences (Landau's theorem gives the feasibility conditions).

From Kolesnik-Sanchez (DCG 2024): tournament enumeration = mixed volume computation on the permutohedron.

## The Gamma Function Unifies Them

The three faces are unified by the Gamma function hierarchy:

**Gamma(1/b)^b** governs the fiber fraction for base-b comparison systems.

At integer b, this connects to:
- **Simplex**: b-ary labelings of the (n-1)-simplex edges
- **Staircase**: b-ary tilings with tile weights b^d
- **Permutohedron**: b-ary score sequences

The reflection formula **Gamma(1/b) * Gamma(1-1/b) = pi/sin(pi/b)** introduces the REGULAR b-GON:

sin(pi/b) is the half-chord of the regular b-gon inscribed in a unit circle. This means:
- **The regular b-gon controls the fiber fraction** through sin(pi/b)
- **The regular polygon and the simplex are dual**: the b-gon lives in 2D (the complex plane), while the simplex lives in (n-1)D
- **The staircase mediates**: it converts the 2D polygon structure (chord lengths) into the (n-1)D simplex structure (edge labelings)

## The Regular b-gon -> Simplex Bridge

For the regular b-gon inscribed in a unit circle:
- b = 3 (equilateral triangle): the simplex in 2D
- b = 4 (square): the hypercube face
- b = 5 (pentagon): the icosahedral structure (Platonic solids end here)
- b = 6 (hexagon): the tiling of the Euclidean plane
- b = infinity (circle): the continuous limit

The simplex Delta_{n-1} has the symmetry group S_n (permutations of vertices). The regular b-gon has the dihedral group D_b (rotations and reflections). These are connected through the Coxeter-Dynkin diagram:

A_{n-1} <-> simplex <-> S_n
I_2(b) <-> regular b-gon <-> D_b

The staircase delta_{n-2} encodes how the A_{n-1} (simplex) structure is "viewed" through the I_2(b) (polygon) lens at the fiber level.

## The Egg Drop Connection

From S20ao: k eggs and n floors correspond to k-dimensional simplices. The egg drop capacity C(x,k) is the number of k-faces of a simplex.

Now: the regular b-gon gives the base of the comparison system. The simplex gives the dimensionality.

So: **a b-ary comparison system with k "destructive resources" lives on a k-simplex inscribed in a b-gon environment**. The capacity is:

Capacity(k, b, n) ~ n^k / (k! * Gamma(1/b)^b)

This generalizes both:
- The egg drop formula: C(x,k) = x^k/k! for k eggs
- The fiber fraction: f ~ 1/Gamma(1/b)^{b/2} for b-ary comparisons

## The Deep Geometric Picture

Tournament space is a **fiber bundle** where:
- **Base**: the permutohedron Pi_n (score sequences)
- **Fiber**: the cycle space at each score (even graphs)
- **Total space**: the (n-1)-simplex's edge labelings

The fiber fraction f(n) = (1/2)_k/k! measures the fiber thickness. Its generating function (1-x)^{-1/2} is the two-sheeted cover of the Riemann sphere branched at x=1.

The **regular b-gon** enters through the Gamma function: the chord length sin(pi/b) determines Gamma(1/b) via the reflection formula, which determines the fiber fraction, which determines the bundle geometry.

The **staircase** is the interface between base and fiber: its rows (indexed by range d) have lengths that decrease linearly, creating the triangular shape. The tile weights 2^d are the place values of the binary adder that computes H from the tiling.

## The Summary

POLYGON <-> SIMPLEX <-> STAIRCASE

are three faces of the same geometric object:
- **Polygon**: the 2D cross-section (circle/chord structure, Gamma function)
- **Simplex**: the high-dimensional structure (vertex set, edge labelings, S_n symmetry)
- **Staircase**: the intermediary (Young diagram, tiling model, transfer matrix)

The number pi = Gamma(1/2)^2 appears because tournaments are BINARY (b=2), and the regular 2-gon has chord length sin(pi/2) = 1. The generalized pi Gamma(1/b)^b for other bases comes from the chord length sin(pi/b) of the regular b-gon.

All three constants -- pi, e, and gamma (Euler-Mascheroni) -- appear in the asymptotic regime: Gamma(1/b)^b ~ b^b * e^{-gamma}, and pi = Gamma(1/2)^2 is the b=2 special case.

**The polygon-simplex-staircase trinity IS tournament theory, viewed from three geometric angles.**
