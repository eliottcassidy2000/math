# Deep Connections: W(r), Cayley Transform, and the Spectral Theory of Tournaments

**Author:** kind-pasteur-2026-03-06-S25f
**Purpose:** Explore the deepest connections between the weighted path polynomial, transfer matrix, and orthogonal group theory

---

## I. The W(r) = tr(M(r)) Identity and Its Meaning

Opus S27 proved: at odd n, the weighted path polynomial

$$W(r) = \sum_P \prod_{e \in P} (r + s_e)$$

equals the transfer matrix trace tr(M(r)), where s_e = A[u,v] - 1/2 is the "signed edge weight."

This is profound because it connects TWO different decompositions of path enumeration:

### Decomposition 1 (W): Product over edges
Each path P contributes a product of (n-1) factors (r + s_e).
This is a "multiplicative" decomposition — the path's contribution is the product of its local edge decisions.

### Decomposition 2 (tr(M)): Inclusion-exclusion over subsets
The transfer matrix M[a,b] = sum_S (-1)^|S| E_a(S) B_b(R) decomposes each path into a "front half" ending at a and a "back half" starting at a.
The trace sums over all "midpoints" a.

The equality W(r) = tr(M(r)) says: **the multiplicative edge-by-edge perspective equals the additive inclusion-exclusion perspective**.

This is analogous to the identity in linear algebra: det(A) = sum_sigma sgn(sigma) prod a_{i,sigma(i)}.
The determinant (a global object) equals a sum of products (local objects).

### Coefficient Consequences

Expanding W(r) = sum_k w_k * r^k where w_k = sum_P e_{n-1-k}(s_P):

- w_{n-1} = H (the Hamiltonian path count — coefficient of r^{n-1})
- w_{n-2} = 0 (odd symmetry from s -> -s exchange)
- w_{n-3} = sum_P e_2(s_P) = 2*(n-2)! * t_3 - (n-2)!*C(n,3)/2

The last is opus's algebraic proof: **the third coefficient of W(r) is a linear function of the 3-cycle count t_3**. This means:
- W(r) at degree n-1: counts paths (graph-level)
- W(r) at degree n-3: counts 3-cycles (cycle-level)
- W(r) at degree n-5: should count 5-cycles (conjecture)

So W(r) is a "generating function in r that stratifies tournament invariants by odd-cycle complexity."

---

## II. Connection to the Cayley Transform

### The Walk Generating Function

Irving-Omar define:
  W_IO(z) = det(I + zA^T) / det(I - zA)

where A is the tournament adjacency matrix. They prove W_IO(z)W_IO(-z) = 1 (reciprocity).

For tournaments, A = (J-I)/2 + S where S is skew-symmetric.
The Cayley transform maps S -> Q = (I+S)(I-S)^{-1}, producing an orthogonal matrix Q.

### The Key Bridge

Our W(r) is NOT the same as Irving-Omar's W_IO(z), but they are related:
- W_IO(z) is a RATIONAL function (ratio of determinants)
- Our W(r) is a POLYNOMIAL (sum of products over paths)

However, both encode the same tournament data through different lenses.

The connection is: the eigenvalues of A (or equivalently S) determine both:
1. The roots of det(I - zA) and det(I + zA^T), hence W_IO
2. The coefficients of W(r) via the transfer matrix eigenvalues

For the Cayley transform Q = (I+S)(I-S)^{-1}:
- The eigenvalues of S are ±i*lambda_k (purely imaginary, since S is real skew-symmetric)
- The eigenvalues of Q are (1+i*lambda_k)/(1-i*lambda_k) = e^{2i*arctan(lambda_k)}
- These lie on the unit circle — Q is orthogonal

So: **the tournament T defines a point on the orthogonal group O(n) via the Cayley transform of its skew part**.

### The Character Interpretation

W_IO(z) = det(I + zA^T) / det(I - zA)
         = prod_k (1 + z*mu_k) / (1 - z*mu_k)

where mu_k are eigenvalues of A.

For tournament: mu_k = 1/2 + i*lambda_k/2 where lambda_k are eigenvalues of S.

The reciprocity W_IO(z)*W_IO(-z) = 1 follows from:
  prod_k [(1+z*mu_k)(1-z*mu_k)] / [(1-z*mu_k)(1+z*mu_k)] = 1

This is a TAUTOLOGY — it holds for any matrix. But the fact that W_IO specializes to tournament-specific invariants at particular z-values is NOT tautological.

CREATIVE INSIGHT: The representation-theoretic content is in the FACTORING of W_IO into characters of irreducible representations of O(n). The tournament constraint (A + A^T = J - I) means the eigenvalues of A lie on a specific AFFINE SUBSPACE of C^n, and this constraint is what produces the rich structure.

---

## III. Pfaffian Structure and Tournament Determinants

For a tournament T on n vertices with skew-symmetric part S:
- If n is even: det(S) = Pf(S)^2 (Cayley's theorem)
- If n is odd: det(S) = 0 (skew-symmetric matrices of odd order are singular)

The Pfaffian Pf(S) counts CYCLE COVERS of the tournament (with signs).
A cycle cover of T is a set of vertex-disjoint directed cycles covering all vertices.

### Connection to H(T)

The number of Hamiltonian paths H(T) and the Pfaffian Pf(S) are both polynomial functions of the tournament entries. But they count fundamentally different objects:
- H(T) counts paths (open chains)
- Pf(S) counts cycle covers (closed loops)

The OCF H(T) = I(Omega(T), 2) bridges this gap by expressing paths in terms of ODD CYCLES (via the conflict graph). The independence polynomial I(Omega,x) decomposes H into contributions from collections of odd cycles.

### A Pfaffian-Path Duality?

SPECULATIVE CONNECTION: Is there a DUALITY between:
- Pfaffian(S) = "cycle covers count" (even n)
- H(T) = "path count" (all n)

Such a duality would say: the number of ways to CLOSE the tournament into cycles is related to the number of ways to OPEN it into paths.

For even n: the transfer matrix M is singular (tr(M) = 0 at odd n... wait, tr(M) = H at odd n, 0 at even n). So at even n, paths "cancel" in the transfer matrix while cycles "survive" in the Pfaffian.

At odd n: paths "survive" (tr(M) = H) while cycles "cancel" (det(S) = 0).

This is a PATH-CYCLE DUALITY governed by the parity of n:
- Odd n: paths dominate, cycles cancel (det(S) = 0, tr(M) = H)
- Even n: cycles dominate, paths cancel (det(S) = Pf(S)^2, tr(M) = 0)

---

## IV. The BIBD Embedding and Finite Geometry

### Recap: Paley Tournaments and Designs

For Paley tournament T_p on prime p:
- Vertices = F_p (field of p elements)
- Arc i->j iff j-i is a quadratic residue mod p
- The cyclic 3-cycles form a 2-(p, 3, (p+1)/4) BIBD

At p=7: this BIBD is the Fano plane PG(2,2), the unique S(2,3,7).
At p=11: this is a 2-(11,3,3) BIBD with 55 blocks.
At p=23: this is a 2-(23,3,6) BIBD.

### The Pin Grid as a Canvas

The pin grid Grid(n) is the staircase Young diagram delta_{n-2} with m = C(n-1,2) cells.
Each cell (i,j) with i < j encodes the arc direction: T[i,j] or T[j,i].

The 3-cycles of T correspond to TRIANGLES in the pin grid:
- Three cells (i,j), (j,k), (i,k) form a triangle
- A directed 3-cycle (i,j,k) exists iff the three cells have the "cyclic" pattern

For the Paley tournament, the BIBD structure means the triangles are MAXIMALLY SPREAD:
every pair of vertices appears in exactly (p+1)/4 cyclic triples.

### Connection to Tiling Maximization

The Hamiltonian path count H(T) = I(Omega(T), 2) depends on:
1. The number of odd cycles (alpha_k coefficients of I.P.)
2. How these cycles INTERSECT (the conflict graph structure)

The Paley tournament maximizes H because:
1. It has the most 3-cycles (regular score sequence)
2. The 3-cycles form a BIBD, meaning they are optimally spread
3. This optimal spreading MINIMIZES conflicts in Omega
4. Fewer conflicts = larger independence polynomial at x=2

This is a PACKING argument: the Paley tournament achieves the "densest non-conflicting packing" of odd cycles, maximizing I(Omega, 2).

### Finite Geometry Perspective

The Fano plane PG(2,2) has:
- 7 points, 7 lines
- Each line contains 3 points
- Each point lies on 3 lines
- Any 2 points determine a unique line

At n=7 Paley, the 14 directed 3-cycles pair into 7 undirected triples,
one per Fano line. The SEVEN independence polynomial coefficient alpha_2 = 7
equals the number of Fano lines! This cannot be coincidence.

CREATIVE INSIGHT: The independence number alpha_k of Omega(T_p) might equal
a GEOMETRIC INVARIANT of the BIBD:
- alpha_1 = number of 3-cycles = p*(p-1)/6 * ... (need to check)
- alpha_2 = number of "compatible pairs" of 3-cycles = ???

For the Fano plane: alpha_2 = 7 = number of lines.
Two 3-cycles are independent (non-conflicting) iff they share no vertex,
which in the Fano plane means... actually the 7 Fano lines DO share vertices.

Let me reconsider. alpha_2 counts pairs of vertex-disjoint 3-cycles.
At n=7 with 7 vertices and 3-vertex cycles: a pair of disjoint cycles uses 6 vertices.
So alpha_2 = number of pairs of disjoint 3-cycles = 7.
These 7 pairs correspond to the 7 ways to partition 6 of the 7 vertices into two triples.

---

## V. The Transfer Matrix as a Random Walk

### M as Transition Matrix

The transfer matrix M[a,b] = sum_S (-1)^|S| E_a(S) B_b(R) can be interpreted as:
- Start a random walk at vertex a
- At each step, choose the next vertex from the "front" or "back" half
- Accumulate a sign (-1)^|S| based on which halves are used
- End at vertex b

The symmetry M[a,b] = M[b,a] (THM-030) means this walk is REVERSIBLE.
This is NOT the same as the usual random walk on the tournament (which is irreversible).

### Connection to Reversibility in Group Theory

For a Cayley tournament on group G with connection set S:
- The random walk is irreversible: Prob(g -> gs) != Prob(gs -> g) unless S = S^{-1}
- But the TRANSFER MATRIX walk is always reversible (M = M^T)

The transfer matrix "symmetrizes" the tournament by the inclusion-exclusion over subsets.
This is analogous to how the Laplacian L = D - A symmetrizes the adjacency matrix.

For SC VT tournaments: M = (H/n)*I means the transfer walk is not just reversible but ISOTROPIC — it goes nowhere on average, distributing equally to all vertices.

For non-SC VT at n=21: M[0,1] = 45,478,409 means the walk has a PREFERRED DIRECTION — it preferentially connects vertices that are "close" in the group structure.

### The Spectrum of M

For the transitive tournament T_n (all arcs i->j for i < j):
- M has eigenvalues at Chebyshev nodes
- |det(M)| = Fibonacci number F_n
- The spectral gap of M measures "how far T is from circulant"

QUESTION: What are the eigenvalues of M for the Paley tournament?
Since Paley is SC + VT, M = (H/n)*I, so all eigenvalues equal H/n.
This means the Paley tournament has MAXIMALLY DEGENERATE M spectrum.

For the transitive tournament: M has maximally SPREAD spectrum.
For the Paley tournament: M has maximally COLLAPSED spectrum.

This is another instance of the "perpendicular maximizer" phenomenon:
maximum H corresponds to maximum spectral degeneracy of M.

---

## VI. The Hopf Algebra and Deletion-Contraction

### Grujic-Stojadinovic Coproduct

The coproduct Delta([T]) = sum_S [T|_S] tensor [T|_{V\S}] is our subset convolution.
This gives the space of tournament isomorphism classes the structure of a HOPF ALGEBRA.

The key property: the coproduct is COCOMMUTATIVE (since tensor product is commutative).
This means the dual algebra is COMMUTATIVE — a polynomial ring.

The GENERATORS of this polynomial ring are the PRIMITIVE elements:
tournaments T such that Delta(T) = T tensor 1 + 1 tensor T.
These are the tournaments that cannot be decomposed as "direct sums."

### Mitrovic's Deletion-Contraction

The noncommutative deletion-contraction W_X = W_{X\e} - W_{X/e} works at the set-partition level. In the commutative specialization, this becomes:

H(T) = H(T\e) - H(T/e)

where T\e deletes an edge and T/e contracts an edge.

For tournaments: "deleting an edge" means removing the arc between two vertices i,j and summing over both possible orientations. "Contracting an edge" means identifying i and j into a single vertex.

This gives a RECURSIVE structure: H(T) can be computed by deletion-contraction on any edge, decomposing the problem into smaller tournaments.

### Connection to the Tutte Polynomial

The Tutte polynomial T_G(x,y) satisfies deletion-contraction for graphs.
For tournaments, the analogue is the Redei-Berge symmetric function U_T.
The chromatic bridge (Mitrovic): X_{inc(P)} = omega(U_P) connects the chromatic symmetric function to U_T.

The Stanley-Stembridge conjecture (PROVED by Hikita 2024) states that the chromatic symmetric function of (3+1)-free posets is e-positive.

CREATIVE CONNECTION: If tournaments are viewed as "total orders with cycles,"
then the Redei-Berge function U_T is the "tournament analogue" of the chromatic function.
The even-cycle vanishing (T148) — p_mu(U_T) = 0 for even-part mu — is the tournament
analogue of the acyclicity condition for chromatic functions.

---

## VII. Number-Theoretic Patterns

### The 1729 Mystery Deepened

H(T_p) / |Aut(T_p)| for Paley tournaments:
- p=3: H=3, |Aut|=3, ratio = 1
- p=7: H=189, |Aut|=21, ratio = 9 = 3^2
- p=11: H=95095, |Aut|=55, ratio = 1729 = 7*13*19

Observations:
- |Aut(T_p)| = p*(p-1)/2 (the affine group AGL(1,p) restricted to QR)
- The ratio counts ORBITS of Hamiltonian paths under Aut
- 1729 = 12^3 + 1^3 = 10^3 + 9^3 (Hardy-Ramanujan taxicab)

### Quadratic Reciprocity Connection

The Paley tournament T_p is defined by the Legendre symbol: i->j iff (j-i|p) = 1.
Quadratic reciprocity says: (p|q) * (q|p) = (-1)^{(p-1)(q-1)/4}

Could quadratic reciprocity connect H(T_p) at different primes?
If (p|q) = 1, then the "p-structure" embeds in F_q via the QR inclusion.
This might explain why |Aut(T_7)| = 21 appears as a ratio at T_11:
  alpha_3(T_11)/|Aut(T_11)| = 1155/55 = 21 = |Aut(T_7)|

### Factorization Pattern

1729 = 7 * 13 * 19
Note: 7 = p, 13 = 2p-1, 19 = ???
Actually: 7, 13, 19 are primes in arithmetic progression with common difference 6.
And 6 = p-1 where p = 7. So: the orbit count at p=11 factors into primes
spaced by (p-1) for the PREVIOUS Paley prime p=7.

This is HIGHLY SPECULATIVE but tantalizing. At p=19 or p=23, we could test
whether the orbit count factors similarly.

### Predicted: H(T_19) and H(T_23)

If the Paley maximizer conjecture holds, H(T_19) should be the maximum
over all 19-vertex tournaments. |Aut(T_19)| = 19*9 = 171.
The orbit count H(T_19)/171 would reveal whether the prime-factorization
pattern continues.

---

## VIII. Summary: The Five Lines of Symmetry

1. **Tournament constraint line** (t+t'=1): Creates M symmetry, even-cycle vanishing, OCF
2. **Self-converse line** (T ~ T^op): Creates palindromic N, scalar M for VT, perpendicular maximizer
3. **Vertex-transitive line** (Aut transitive): Creates orbit structure, Cayley connection
4. **Paley/QR line** (Legendre symbol): Creates BIBD, maximum H, number-theoretic patterns
5. **Spectral line** (eigenvalues of A): Creates Cayley transform point, W(z) reciprocity, real roots

These five lines intersect at the PALEY TOURNAMENT, which lies at the maximum of H
and at the center of the symmetry structure. The Paley tournament is the unique point
where all five symmetries are simultaneously maximized.

---

## IX. Open Questions from This Analysis

1. **PATH-CYCLE DUALITY**: Is the duality between H and Pf(S) at odd/even n a formal mathematical duality (e.g., Poincare duality on some space)?

2. **SPECTRAL DEGENERACY**: Does maximum H always correspond to maximally degenerate M spectrum? (True for Paley; check for non-SC maximizers.)

3. **BIBD PACKING**: Can we formalize the "Paley maximizes H because BIBD maximizes non-conflicting cycle packing" argument? This would require showing I(Omega, 2) is maximized when cycles form a BIBD.

4. **W(r) STRATIFICATION**: Do the coefficients w_{n-2k-1} of W(r) correspond to (2k+1)-cycle counts? If so, W(r) is a "cycle-stratified generating function."

5. **QUADRATIC RECIPROCITY AND H**: Is there a formal connection between quadratic reciprocity and the embedding of H(T_p) invariants in H(T_q)?

6. **HOPF PRIMITIVE ELEMENTS**: What are the primitive tournaments in the Grujic-Stojadinovic Hopf algebra? Are they the "indecomposable" tournaments? How does this connect to the OCF?
