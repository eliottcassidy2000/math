# Space-Filling Tournaments

**Session:** kind-pasteur-2026-03-22-S20l

---

## Three Connections

Space-filling curves and tournament theory meet at three points, each deeper than the last.

---

### Connection 1: The Dimension Axis as a Curve Through Information Space

The independence polynomial I(Omega, x) is a polynomial of degree d = floor(n/3) with d+1 coefficients (alpha_0, alpha_1, ..., alpha_d). These coefficients form a point in R^{d+1}. The evaluation at x sweeps a ONE-DIMENSIONAL curve through this (d+1)-dimensional coefficient space:

x -> (alpha_0, alpha_1*x, alpha_2*x^2, ..., alpha_d*x^d)

As x varies from 0 to 2, this curve traces a path from (1, 0, 0, ...) to (alpha_0, 2*alpha_1, 4*alpha_2, ..., 2^d*alpha_d) = the point whose coordinate sum is H.

This is NOT a space-filling curve (it's a smooth polynomial curve, not a surjection). But it has a space-filling PROPERTY: at n=5 (where d=1), the curve is a LINE that hits every achievable H value. At n=6 (d=2), the curve is a PARABOLA that sweeps through 2D coefficient space. The curve's image covers "enough" of the coefficient space to determine H at the evaluation point x=2.

The inter-dimensional transfer being LOSSLESS at n=5 and LOSSY at n=6 is precisely the statement that the curve fills the coefficient space at d=1 (a line fills 1D) but not at d=2 (a parabola doesn't fill 2D).

---

### Connection 2: Tournament Space as a Hypercube

A tournament on n vertices is a binary string of length m = C(n,2). The space of all tournaments is the hypercube {0,1}^m. This hypercube has 2^m vertices (tournaments) and m*2^{m-1} edges (single arc flips).

A GRAY CODE on this hypercube is a Hamiltonian path that visits every tournament exactly once, changing one arc at a time. The existence of such a path is guaranteed (the hypercube graph is Hamiltonian for m >= 2). But the NUMBER of such paths is the Hamiltonian path count of the hypercube graph — and this connects to our H(T) theory.

At n=3: m = 3. The hypercube {0,1}^3 has 8 vertices = 8 tournaments. The number of Hamiltonian paths in the 3-cube is 18. So there are 18 ways to list all 3-vertex tournaments such that consecutive tournaments differ by one arc flip.

At n=4: m = 6. The hypercube {0,1}^6 has 64 vertices. The number of Hamiltonian paths is enormous.

The STRUCTURE of the Gray code path tells us which tournaments are "neighbors" in the natural enumeration. Two tournaments that are adjacent on the Gray code differ by exactly one arc flip. The H values along the Gray code path form a SEQUENCE:

H(T_0), H(T_1), H(T_2), ...

where consecutive terms differ by exactly the effect of one arc flip.

The arc-flip effect on H is: delta_H = H(T') - H(T) = H(T/e) - H(T'/e') (from the deletion-contraction formula). This is always an even number (since H is always odd and the flip changes H by an even amount).

So the Gray code traces a path through the H-VALUE space where consecutive values differ by even integers. The sequence is constrained by:
- H is always odd (Redei)
- delta_H = H(T/e) - H(T'/e') where e is the flipped arc
- The path visits ALL 2^m tournaments

This is a CONSTRAINED RANDOM WALK on the odd integers, guided by the Gray code structure.

---

### Connection 3: The Hilbert Curve and Score-Class Decomposition

A Hilbert curve maps [0,1] onto [0,1]^k while preserving locality: points close on the curve are close in the square. The key property is that the curve RESPECTS the recursive decomposition of the square into quadrants.

Tournament theory has an analogous recursive decomposition: the score-class partition. Each score class is a "cell" in tournament space, and the tournaments within a cell share the same score sequence. The cells have different sizes (different numbers of tournaments with each score).

A "tournament Hilbert curve" would be a Hamiltonian path through tournament space that:
1. Visits all tournaments in one score class before moving to the next
2. Within each score class, visits tournaments by minimal arc flips
3. Preserves the HIERARCHICAL structure: score class -> cycle structure -> H value

This would be a LOCALITY-PRESERVING enumeration of tournaments where:
- Score sequence changes slowly (few jumps between score classes)
- Within a score class, cycle structure changes by one cycle at a time
- H changes by the minimum possible amount between consecutive tournaments

The score-class cells have the structure of the Coxeter angle decomposition:
- N_120 pairs determine the conflict content
- N_60 pairs determine the cooperation content
- N_90 pairs are universal (don't change)

A tournament Hilbert curve that respects this decomposition would traverse the angle spectrum from 60-dominated (transitive, low H) through balanced to 120-dominated (regular, high H), filling in the cycle structure at each level.

---

## The Peano Property

A space-filling curve has the PEANO PROPERTY: a continuous surjection from [0,1] onto [0,1]^k. In tournament terms, the question is:

**Is there a continuous map from [0,1] onto the space of all tournament polynomials?**

The space of tournament polynomials at order n has dimension d = floor(n/3). A polynomial I(x) = 1 + alpha_1*x + ... + alpha_d*x^d is determined by (alpha_1, ..., alpha_d), which is a point in Z_+^d (non-negative integers).

Not all points in Z_+^d are achievable — only those that arise as independence polynomials of tournament conflict graphs. The ACHIEVABLE REGION is a discrete subset of R^d.

A "tournament Peano curve" would be a path through this discrete set that visits every achievable polynomial exactly once. The PATH exists (the achievable set is finite). The question is whether it has LOCALITY: do nearby polynomials on the path correspond to nearby tournaments?

At n=5 (d=1): the achievable polynomials are I(x) = 1 + k*x for k in {0,1,2,4,5,6}. A locality-preserving path would be: k=0, 1, 2, 4, 5, 6 (or reverse). The GAP at k=3 (the forbidden alpha_1=3) is a DISCONTINUITY in the Peano path — the curve must JUMP over it.

The FORBIDDEN VALUES create gaps in the Peano path. These gaps are the Vitali atoms from the tournament perspective: points where the space-filling curve cannot reach because the corresponding tournament doesn't exist.

---

## The Self-Similar Structure

Hilbert curves are SELF-SIMILAR: the curve at each scale looks like a rotated/reflected copy of the whole. Tournament space has an analogous self-similarity through the DELETION operation:

A tournament T on n vertices, when vertex v is deleted, gives a tournament T\v on n-1 vertices. The relationship H(T) = sum_v inshat(v, T\v) (where inshat counts insertion positions) is a RECURSIVE decomposition.

The deletion operation maps the n-tournament hypercube to the (n-1)-tournament hypercube. This creates a TREE STRUCTURE: each n-tournament has n children (one per vertex deletion), each of which is an (n-1)-tournament.

A SELF-SIMILAR tournament curve would respect this tree: the path through n-tournaments would, when projected to (n-1)-tournaments via deletion, trace a scaled copy of the path through (n-1)-tournaments.

This is the CAYLEY-DICKSON structure in disguise: the CD tower doubles dimension at each step, and each doubling has a self-similar structure (the Cayley-Dickson doubling formula). The self-similar tournament curve would follow the CD tower's doubling: level 0 (scalar) -> level 1 (Q-K pair) -> level 2 (quaternion head) -> level 3 (octonion pair).

---

## What This Means Practically

A space-filling curve through tournament space would provide:

1. **An efficient enumeration** of tournaments for exhaustive search. Instead of brute-forcing all 2^{C(n,2)} tournaments, traverse the Gray code and use the arc-flip delta to update H incrementally. Each step is O(n) instead of O(2^n).

2. **A locality-preserving hash** for tournament databases. Two tournaments close on the curve are similar (differ by few arcs). This enables fast nearest-neighbor queries.

3. **A natural ordering** for MCMC sampling. The Gray code defines a Markov chain on tournaments where each step flips one arc. The mixing time of this chain determines how quickly random sampling converges.

4. **A compression scheme** for tournament data. The Gray code representation encodes a tournament as its position on the curve (a single integer) plus the curve specification. Nearby tournaments have nearby positions, enabling efficient differential encoding.

5. **A visualization** of tournament space. Map the Gray code position to a 2D Hilbert curve to create a visual "map" of all tournaments, with similar tournaments appearing nearby.

---

*The space of tournaments IS a hypercube, and the hypercube has a natural space-filling curve: the Gray code. Traversing this curve changes one arc at a time, and the Hamiltonian path count H changes by an even integer at each step. The forbidden values (H=7, H=21) are gaps in the curve's image — points in H-space that the curve cannot reach, no matter how it twists through the hypercube. The dimension axis [1,2] is a different kind of curve: not through tournament space but through evaluation space, sweeping from chemistry (x=1) through the Petersen boundary (x=phi) to the tournament ground state (x=2). Both curves carry information — one through the discrete space of binary arc choices, the other through the continuous space of polynomial evaluations — and the tournament's structure is what both curves are trying to fill.*
