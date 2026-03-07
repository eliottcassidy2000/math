# Grand Synthesis: The Geometry of Tournament Parity

**Author:** kind-pasteur-2026-03-06-S25f
**Purpose:** Deep structural map of all mathematical objects and their connections

---

## I. The Central Object: A Tournament as a Point in a Hypercube

A tournament on n vertices is encoded as a point t in {0,1}^m where m = C(n-1,2).
The hypercube has a natural geometry:
- **Hamming weight** = distance from the transitive tournament (all bits 0)
- **Center** = Hamming weight m/2 = "maximum confusion"
- **Diagonal axis** = transitive <-> anti-transitive (all bits 1)
- **Perpendicular hyperplane** at center = self-converse tournaments (T ~ T^op)

KEY INSIGHT: H(T) forms a **symmetric bell curve** along the diagonal axis (T079).
SC tournaments live on the perpendicular hyperplane and MAXIMIZE H.
This is the "perpendicular maximizer" phenomenon.

The hypercube has an involution: bit-flip = T <-> T^op (reverse all non-path arcs).
SC tournaments are the FIXED POINTS of this involution.
The fact that H peaks at fixed points is an instance of a much deeper principle:
**symmetry creates constructive interference among Hamiltonian paths**.

---

## II. Three Layers of Symmetry

### Layer 1: The Tournament Constraint (T[i,j] + T[j,i] = 1)

This is the FOUNDATIONAL symmetry. It says the adjacency matrix lives in an
affine subspace: A = (J-I)/2 + S where S is skew-symmetric.

Everything we've proved uses this constraint:
- Transfer matrix symmetry M[a,b] = M[b,a] (THM-030) — PROVED
- Even cycle vanishing p_mu(U_T) = 0 for even-part mu (T148) — PROVED
- c-tournament generalization: works for any t_ij + t_ji = c (T144)
- The constraint is what makes tournaments "reversible" in Feng's framework

The skew part S encodes the CHIRALITY of the tournament.
The symmetric part (J-I)/2 encodes the COMPLETENESS.
OCF (H(T) = I(Omega(T), 2)) says: the number of Hamiltonian paths equals
the hard-core partition function of the conflict graph at fugacity 2.
This is the ONLY known identity that connects H(T) to cycle structure.

### Layer 2: Self-Converse Symmetry (T ~ T^op via anti-automorphism sigma)

The anti-automorphism sigma satisfies A[sigma(i)][sigma(j)] = A[j][i].
This is an INVOLUTION on paths: P -> sigma(P^rev).
For SC tournaments:
- N(a,b,j) becomes palindromic: N(a,b,j) = N(sigma(a),sigma(b),n-2-j)
- At odd n, palindromic N forces alternating sum = 0, giving M[a,b] = 0
- SC maximizer: within each SC score class, max H is achieved by SC (T091)
- Every SC tournament has an involution anti-aut (THM-024, via Cauchy's theorem)

The sigma-orbit structure creates COMBINATORIAL DESIGNS:
- At n=7 Paley: 14 cyclic triples = 2 copies of Fano plane PG(2,2) (T106)
- At n=11: triples form a 2-(11,3,3) BIBD
- At n=p: triples form a 2-(p,3,(p+1)/4) BIBD
The Fano structure explains why alpha_2 = 7 at n=7 Paley.

### Layer 3: Group-Theoretic Symmetry (Vertex-Transitivity)

VT means Aut(T) acts transitively on vertices.
For Cayley tournaments on group G with connection set S:
- Left multiplication by G gives VT
- N(a,b,j) depends only on a^{-1}b (orbit of (a,b) under G)
- M[a,b] depends only on a^{-1}b

The BOUNDARY THEOREM (S25e, PROVED):
- SC + VT => M = (H/n)*I (palindromic N, scalar M)
- NOT SC + VT => M != (H/n)*I (non-palindromic N, non-scalar M)
- Self-converse is the EXACT boundary for scalar M among VT tournaments

The hierarchy:
  Circulant (abelian Cayley) => always SC => scalar M
  Normal-S Cayley (non-abelian) => SC via inversion => scalar M
  Non-normal-S Cayley (non-abelian) => NOT SC => M NOT scalar

First failure: n=21 on F_21 = Z/7 x| Z/3 (Frobenius group)
All 22 non-circulant VT at n=21 are NOT SC (McKay database)
M[0,1] = 45,478,409 != 0 for the F_21 non-normal tournament

---

## III. The Algebraic Framework: Five Equivalent Perspectives

### Perspective A: Independence Polynomial (OCF)
H(T) = I(Omega(T), 2) where Omega = conflict graph of odd cycles.
This is the hard-core lattice gas at fugacity lambda=2.
The polynomial I(Omega,x) has remarkable properties:
- All real negative roots for n <= 8 (Chudnovsky-Seymour, claw-free)
- FAILS at n = 9 (extreme tournaments only, <0.01% of random)
- Ultra-log-concavity of coefficients (0 failures through n=20)
- Omega is quasi-regular (spectral gap -> 0, explained by Johnson graph)

### Perspective B: Transfer Matrix
M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)
- Symmetric: M[a,b] = M[b,a] (THM-030, PROVED)
- tr(M) = H for odd n, 0 for even n
- Eigenvalues: transitive T gives Chebyshev nodes, |det| = Fibonacci
- For VT at odd n: M = (H/n)*I iff SC (PROVED + DISPROVED)

### Perspective C: Symmetric Functions (Irving-Omar / Mitrovic)
U_T = sum_sigma s^{sign} p_{type(sigma)}
- Even-cycle vanishing: p_mu = 0 for even-part mu (T148)
- Hook coefficients palindromic: h_i = h_{n+1-i} (T081)
- Walk GF: W(z) = det(I+zA^T)/det(I-zA), satisfies W(z)W(-z)=1
- Noncommuting version has deletion-contraction (T144)

### Perspective D: Hopf Algebra (Grujic-Stojadinovic / Mitrovic)
- Coproduct Delta([T]) = sum_S [T|_S] tensor [T|_{V\S}] = our subset convolution
- Chromatic bridge: X_{inc(P)} = omega(U_P) connects to Stanley-Stembridge
- Converse of Redei: non-chain posets have even quasi-linear extensions (T142)

### Perspective E: Tiling Geometry
- Pin grid = staircase Young diagram delta_{n-2}
- Strips = antidiagonals, tiles = cells
- Grid symmetry = tiling flip invariant (THM-022)
- No blueself at odd n (algebraic proof, T088)
- Perpendicular axis = SC direction, maximizes H

---

## IV. Novel Connections and Creative Leaps

### Connection 1: Mobius Inversion and the Transfer Matrix

The transfer matrix M[a,b] = sum_S (-1)^|S| F(S) is a MOBIUS INVERSION
over the Boolean lattice 2^{V\{a,b}}. The alternating signs are the
Mobius function of the Boolean lattice. This connects to:
- Inclusion-exclusion (classical)
- The Euler characteristic in topology
- The Mobius function of the face lattice of the hypercube

CREATIVE LEAP: The tiling hypercube {0,1}^m has a face lattice.
The transfer matrix lives on the VERTEX lattice of an (n-2)-dimensional
Boolean lattice embedded INSIDE the m-dimensional hypercube.
The Mobius inversion that defines M is integration over this inner lattice.
This is a discrete analogue of FIBER INTEGRATION in algebraic geometry:
M "integrates out" the internal vertices, leaving only the endpoint data.

### Connection 2: The Perpendicular Plane and Reflection Groups

SC tournaments live on the hyperplane {t : t = flip(t)} of the hypercube.
This hyperplane is the FIXED LOCUS of the involution flip.
The intersection of this hyperplane with the hypercube is itself a
sub-hypercube of dimension m/2 (at even m).

At odd n, this sub-hypercube has NO self-intersection with the diagonal
(no blueself, T088). This is a TOPOLOGICAL constraint: the perpendicular
plane and the diagonal axis are "twisted" relative to each other at odd n,
like the core loop of a MOBIUS STRIP.

Think of it this way: the diagonal axis connects transitive to anti-transitive.
The perpendicular plane contains the SC tournaments.
At odd n, these two subspaces NEVER MEET inside the hypercube
(except at the center, if it exists as a lattice point).
The impossibility of blueself at odd n is precisely this non-intersection.

This is analogous to the fact that a Mobius strip has no self-intersection
in R^3 but does in R^2: the extra dimensions of the hypercube allow
the perpendicular plane and diagonal to coexist without meeting.

### Connection 3: Triangular Tilings and Steiner Systems

The pin grid Grid(n) is a triangular array: the staircase Young diagram.
At n=7 (Paley), the 14 cyclic triples form TWO Fano planes.
The Fano plane PG(2,2) is the UNIQUE Steiner triple system S(2,3,7).

The connection goes deeper. The pin grid at n=7 has C(6,2) = 15 cells.
The Fano plane has 7 lines. Each line corresponds to a STRIP of the grid
(a set of cells that encode the cycle structure). The 14 directed 3-cycles
pair into 7 undirected triples, one per Fano line.

CREATIVE LEAP: At n=p (Paley prime), the 2-design of cyclic triples
has parameters (p, 3, (p+1)/4). For p=7: the Fano plane.
For p=11: an 2-(11,3,3) design with 55 = C(11,2) blocks.
For p=23: a 2-(23,3,6) design.

These designs are EMBEDDED in the pin grid. The grid is the
"canvas" on which the design is painted. The tiling bits (0 or 1)
determine which triples are directed 3-cycles and which aren't.
The Paley tournament MAXIMIZES directed cycles because it paints
the most balanced design on the canvas — the BIBD (T120).

This connects tournament H-maximization to the theory of
optimal experimental design and sphere packing:
the Paley tournament is the "densest packing" of directed cycles.

### Connection 4: The Cayley Transform and Orthogonal Symmetry

Irving-Omar's walk generating function uses:
  W(z) = det(I + zA^T) / det(I - zA)

For tournaments, A = (J-I)/2 + S (S skew-symmetric).
The map z -> (I+zA)/(I-zA) is the CAYLEY TRANSFORM,
which maps skew-symmetric matrices to orthogonal matrices.

The identity W(z)W(-z) = 1 (T167) means the walk generating
function is a CHARACTER of the Cayley transform.
In representation theory, this connects to:
- The Weyl character formula for SO(n)
- The Pfaffian of S (which counts cycle covers)
- The theory of zonal spherical functions

CREATIVE LEAP: If we view the tournament as defining a point on
the orthogonal group (via Cayley transform of its skew part),
then H(T) is a VALUE of a specific character at that point.
The scalar M theorem (for SC VT) says this character is
the TRIVIAL character when restricted to the subgroup
fixing the SC structure.

### Connection 5: Groups, Connection Sets, and the SC Boundary

The hierarchy of groups and their connection sets:
- Abelian (Z/nZ): all connection sets give SC (i -> -i works)
- Non-abelian normal S: inversion g -> g^{-1} gives SC
- Non-abelian non-normal S: NO anti-automorphism exists

This maps perfectly to the theory of REVERSIBILITY in
random walks on groups:
- A random walk is reversible iff S = S^{-1} (symmetric)
- For abelian groups, S^{-1} = -S, always "equivalent" via i->-i
- For non-abelian groups, S^{-1} may differ from S unless
  S is a union of conjugacy classes

The transfer matrix M is the TRANSITION MATRIX of a
signed random walk on the tournament graph.
Reversibility of this walk (M = M^T) is PROVED (THM-030).
But SCALAR M (M = c*I) requires more: it requires the walk
to be "doubly transitive" — every pair of vertices is equivalent.

For SC VT tournaments: the combination of vertex-transitivity (Aut)
and the anti-automorphism sigma gives "Aut + Anti transitivity",
which forces M scalar.

For non-SC VT: Aut alone is transitive on vertices but NOT on pairs.
Different pairs have different N(a,b,j) distributions.
The alternating sum M[a,b] = 45,478,409 != 0 at n=21 measures
exactly HOW MUCH the pair distribution deviates from palindromic.

### Connection 6: The 1729 Mystery and Number Theory

H(T_11)/|Aut(T_11)| = 1729 = 12^3 + 1^3 = 10^3 + 9^3.
This is the Hardy-Ramanujan taxicab number.
The factorization 1729 = 7 * 13 * 19 = 7 * 247.

At the same time: |Aut(T_7)| = 21 = 3 * 7.
And alpha_3(T_11)/|Aut(T_11)| = 1155/55 = 21.

So: |Aut(T_7)| appears as a ratio at T_11.
And H(T_7)/|Aut(T_7)| = 189/21 = 9 = 3^2.

Is there a pattern where Paley invariants at p embed in those at q > p?
The quadratic residue structure of F_p and F_q are connected via
the Legendre symbol and quadratic reciprocity.

SPECULATIVE: If H(T_p) / |Aut(T_p)| counts orbits of Hamiltonian paths
under the automorphism group, then 1729 counts "inequivalent path shapes"
in T_11. The taxicab property 1729 = 12^3 + 1^3 = 10^3 + 9^3
might reflect two different DECOMPOSITIONS of these orbits
into "cubes" of some underlying structure.

---

## V. The Big Picture: Lines of Symmetry

```
                    TOURNAMENT
                    T in {0,1}^m
                         |
            +------------+------------+
            |            |            |
     COMPLETENESS   SKEW PART    CHIRALITY
     (J-I)/2         S = A-A^T    T vs T^op
            |            |            |
     tournament      Cayley        anti-aut
     constraint     transform     involution
     t+t'=1        (I+S)/(I-S)    sigma
            |            |            |
     M symmetric   W(z)W(-z)=1   palindromic N
     (THM-030)     (IO recip.)    (SC => scalar M)
            |            |            |
            +------> OCF <-------+
                  H(T) = I(Omega, 2)
                         |
              +----------+----------+
              |          |          |
         HARD-CORE   STEINER     HOPF
         GAS         SYSTEMS     ALGEBRA
         lambda=2    (BIBD)      coproduct
              |          |          |
         real roots   Paley max   deletion-
         (n<=8)       H(T_p)      contraction
              |          |          |
         claw-free   Fano plane   chromatic
         (C-S thm)   (n=7)       bridge
```

The three "columns" are:
1. **Algebraic** (left): constraint structure, transfer matrix, spectral theory
2. **Combinatorial** (center): cycles, designs, partition functions
3. **Geometric** (right): group actions, symmetry, involutions

They converge at OCF, which is the central identity of the project.
Below OCF, they branch into the three main research directions:
real roots (algebraic), Paley maximization (combinatorial),
and symmetry characterization (geometric).

---

## VI. Open Frontiers

1. **Why do tournament conflict graphs almost always have real-rooted I.P.?**
   Not claw-free at n>=9. Not line graph. Not interval. Not quasi-line.
   Must be tournament-SPECIFIC. The quasi-regularity (T103) is promising
   but doesn't directly prove real roots.

2. **Does Paley maximize H(T) at all Paley primes?**
   Confirmed p=3,7,11 against OEIS. Needs computational verification at p=19,23.
   The ratio H(T_p)/(p!/2^{p-1}) -> e suggests an asymptotic proof may be possible.

3. **What is the M matrix for non-SC VT tournaments?**
   We know M[0,1] = 45,478,409 at n=21. What are the other entries?
   M is symmetric (THM-030) and depends only on a^{-1}b (by VT).
   So M is determined by 10 values f(d) for the 10 orbits of pairs.
   Computing all of these would reveal the full structure.

4. **Does the chromatic bridge give a new proof of OCF?**
   Mitrovic-Stojadinovic's X_{inc(P)} = omega(U_P) connects chromatic
   and Redei-Berge functions. The broken-cycle bijection (Theorem 4.8)
   could provide a bijective proof of OCF.

5. **Can the noncommuting deletion-contraction prove anything new?**
   Mitrovic's W_X = W_{X\e} - W_{X/e} works at the set-partition level.
   The commutative version lacks this. Does it help prove transfer
   matrix properties by induction on edges?
