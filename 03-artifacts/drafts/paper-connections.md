# Connections to External Papers

## CRITICAL NEW DISCOVERY (from this analysis session)

**The Transfer Matrix M[a,b] = sum_S (-1)^|S| E_a(S)*B_b(M\S) is SYMMETRIC.**

Verified: n=4 (3000/3000), n=5 (3000/3000), n=6 (1200/1200), n=7 (300/300), n=8 (60/60).

This means: sum_S (-1)^|S| E_i(S)*B_j(R) = sum_S (-1)^|S| E_j(S)*B_i(R).

The Even-Odd Split Lemma is an IMMEDIATE COROLLARY (off-diagonal difference = 0).
But symmetry is STRONGER: it constrains the diagonal entries M[i][i] and M[j][j] too.

Connection to Feng's Dual Burnside: Q = AB is symmetric when the chain satisfies
detailed balance (reversibility). Our transfer matrix being symmetric suggests a hidden
reversibility in the Hamiltonian path decomposition — the "forward leg" (paths ending at
a vertex) and "backward leg" (paths starting at a vertex) are dual in the Burnside sense.

---

## Paper 1: "A Theory of Tournament Representations" (Rajkumar, Veerathu, Mir 2021)
arXiv: 2110.05188

### Summary
Develops a structural theory for embedding tournaments in R^d via skew-symmetric matrices.
Key concepts: flip classes, forbidden configurations, locally transitive tournaments, doubly
regular tournaments, feedback arc sets, sign rank.

### Direct Connections to Our Work

**CONNECTION 1: Flip Classes and Arc-Flip Induction**

Their "flip class" F(T) = {phi_S(T) : S subset [n]} is exactly the orbit of T under the
operation of reversing all arcs across a cut (S, S^c). Our arc-flip induction proof strategy
(THM-015) flips ONE arc at a time, walking from the transitive tournament to any target.
Their flip reverses ALL arcs across a cut simultaneously.

KEY INSIGHT: Every tournament is in the flip class of some R-cone (their Prop 1). An R-cone
is a tournament with a "universal" vertex (beats or loses to everyone). This is exactly the
structure we exploit in vertex deletion: T with vertex v is an R-cone over T-v when v beats
everyone (or loses to everyone). Our H(T) - H(T-v) = 2*sum mu(C) formula is about the
relationship between the cone and its base.

**CONNECTION 2: Locally Transitive Tournaments = Rank 2**

Their Theorem 6: A tournament is locally transitive iff it's rank 2 (representable by a
rank 2 skew-symmetric matrix). Locally transitive means: for every vertex i, both T_i^+
(vertices i beats) and T_i^- (vertices beating i) are transitive sub-tournaments.

RELEVANCE: For locally transitive tournaments, the conflict graph Omega(T) has a very
constrained structure. If T_i^+ is always transitive, then any odd cycle must "mix" the
in-neighborhood and out-neighborhood of every vertex it contains. This severely limits
which odd cycles can exist and should make our OCF formula I(Omega(T), 2) = H(T) easier
to prove for this class.

CONJECTURE: For locally transitive tournaments, every odd cycle is a 3-cycle, and
Omega(T) is the "triangle graph" of T. If true, OCF reduces to a statement about 3-cycles
only, which is our THM-005 setup.

**CONNECTION 3: Doubly Regular Tournaments and the Hadamard Conjecture**

Their Conjecture 1: Rank 2d tournaments forbid (4d-1)-doubly-regular-cones. They show
this implies the Hadamard conjecture. Doubly regular tournaments on n vertices have
|T_i^+| = |T_j^+| and |T_i^+ ∩ T_j^+| = k for all i,j (constant).

RELEVANCE: Doubly regular tournaments are extremal for many tournament properties.
Since H(T) is always odd (Redei's theorem), and our formula H(T) = I(Omega(T), 2) holds,
this constrains I(Omega, 2) mod 2 for ALL tournaments. For doubly regular tournaments,
the cycle structure of Omega is highly symmetric, potentially giving closed-form expressions
for I(Omega, 2).

**CONNECTION 4: Feedback Arc Sets and Arc-Flip Distance**

Their mu(T) = min_sigma min_{T' in F(T)} theta(sigma, T') measures the minimum "flip
feedback node set" — the fewest vertices involved in feedback arcs across all tournaments
in the flip class. Our arc-flip induction traverses a path in the "arc-flip graph" from the
transitive tournament to T, proving delta_H = delta_I at each step.

POTENTIAL: Their Theorem 11 says any tournament can be represented in dimension
2(mu(T)+1). If mu(T) is small, the tournament is "close to transitive" in a precise sense.
Our induction proof might be structurally simpler for low-mu tournaments, with the general
case requiring mu(T) induction steps. This could give a proof by induction on mu(T).

**CONNECTION 5: Sign Rank and Polynomial Identities**

Their Theorem 12 connects tournament representation to sign-rank of matrices. Our Signed
Position Identity sum_{P: i->j} (-1)^{pos(i)} = sum_{P': j->i} (-1)^{pos(j)} is a
polynomial identity in the arc variables. The sign pattern of the tournament matrix
determines which monomials appear. Sign-rank bounds constrain which cancellation patterns
are possible, potentially providing a route to proving our identity.

SPECULATIVE: The tournament matrix T can be written as T = (J + M)/2 where J is all-ones
and M is a {+1,-1} skew-symmetric matrix. Our polynomial identity lives in the ring
Z[{T[a][b]}] / (T[a][b]+T[b][a]-1). The sign-rank of M constrains the algebraic structure
of this quotient ring.

---

## Paper 2: "The Dual Burnside Process" (Feng 2025)
arXiv: 2510.25202

### Summary
Studies a dual Markov chain to the classical Burnside process for sampling from group orbits.
Key result: the dual has stationary law pi(g) proportional to |X_g| (fixed-point set size),
is reversible, and shares nonzero spectrum with the primal chain.

### Direct Connections to Our Work

**CONNECTION 6: S_n Acting on Tournaments via Coordinate Permutation**

The paper's coordinate-permutation model has G = S_n acting on X = [k]^n by permuting
coordinates. For tournaments, the relevant action is S_n acting on the set of all tournaments
on n vertices by relabeling. The orbits are isomorphism classes of tournaments.

Our OCF identity H(T) = I(Omega(T), 2) is isomorphism-invariant: if T' is isomorphic to T,
both sides are equal. This means the identity is constant on S_n-orbits.

The Burnside process on tournament isomorphism classes would have:
- X = set of all tournaments on [n]
- G = S_n (vertex permutations)
- Orbits = non-isomorphic tournaments
- X_g = {T : g(T) = T} = tournaments with automorphism g
- G_T = Aut(T) = automorphism group

The dual chain walks on permutations g with |X_g| > 0 (permutations that fix at least
one tournament).

**CONNECTION 7: Block-Flip Matrix and Our Arc-Flip Structure**

Their block-flip matrix M = [[0,A],[B,0]] with M^2 = [[Q,0],[0,K]] decomposes the
dynamics into two alternating legs. This is structurally parallel to our arc-flip induction:

Our setup: Start from tournament T. Flip arc (i,j) to get T'. We have:
  H(T') - H(T) = delta_H,  I(Omega(T'),2) - I(Omega(T),2) = delta_I

The "forward leg" A maps tournaments to arc-flip changes, the "backward leg" B maps
changes back to tournaments. Two steps of the block-flip matrix reproduce the dynamics
of the individual chains. Compare to our Even-Odd Split Lemma: the alternating sum
sum (-1)^|S| Delta(S,R) = 0 is a vanishing condition on a "block-flip" type structure
where even-S and odd-S terms alternate.

**CONNECTION 8: Primal-Dual Spectral Correspondence and Signed Position Identity**

Their Theorem 3.9: Spec_{!=0}(Q) = Spec_{!=0}(K) — the dual and primal chains share
nonzero eigenvalues. The eigenvector correspondence (Theorem 3.10) shows the forward
and backward legs A, B are mutual inverses on eigenspaces.

ANALOGY: In our Signed Position Identity, we have two "sides":
  LHS = sum_{P: i->j} (-1)^{pos(i)}   (paths where i precedes j)
  RHS = sum_{P': j->i} (-1)^{pos(j)}   (paths where j precedes i)

These are "dual" sums: one fixes i's role, the other fixes j's role. The identity LHS = RHS
says these dual viewpoints agree, just as the primal and dual Burnside chains share spectra.

DEEPER: The factorization Q = AB, K = BA where A and B are stochastic matrices is
exactly the structure we'd need. Define:
  A(P, S) = "contribution of path P to the S-term in the adj decomposition"
  B(S, P') = "contribution of subset S to path P'"
Then sum_S A(P,S)*B(S,P') would give a "transfer kernel" and the identity AB = BA
(spectrally) could encode our alternating sum vanishing.

**CONNECTION 9: Lumping, Orbits, and the Independence Polynomial**

Their lumping theory (Section 3.5) shows when a chain on a fine state space projects
cleanly to a chain on orbits. The TV-preservation criterion (Theorem 3.30) requires
constant sign on blocks.

For our problem: The "fine" state space is individual directed Hamiltonian paths.
The "lumped" state space could be:
(a) Paths grouped by which odd cycles they "cover"
(b) Paths grouped by insertion signature

Our sum H(T) = sum_P 1 (counting all paths) needs to equal I(Omega(T), 2) = sum over
independent sets of odd cycles. This is exactly a lumping: we need to show the path count
distributes correctly onto independent sets of Omega.

The independence polynomial I(G, x) = sum_{k} alpha_k x^k counts independent sets
weighted by x^k. At x=2, each independent set of k cycles contributes 2^k. This is the
weight a set of k vertex-disjoint odd cycles should receive from the path count.

**CONNECTION 10: Fixed-Point Counting and Tournament Automorphisms**

The dual Burnside stationary law pi(g) proportional to |X_g| counts how many tournaments
are fixed by permutation g. For the identity permutation, |X_e| = |X| = all tournaments.
For a transposition (i j), |X_{(ij)}| counts tournaments where swapping vertices i,j is
an automorphism (T has the same structure whether you call a vertex i or j).

RELEVANCE: Our arc-flip identity concerns what happens when we reverse the arc between
i and j. A tournament where (i j) is an automorphism satisfies T[i][k] = T[j][k] for all
k != i,j, meaning the arc between i and j is the ONLY difference between viewing from i
vs j. For such tournaments, our delta_H should have a particularly clean form, potentially
giving a base case for the general identity.

---

## Creative Synthesis: A Unified View

### The "Spectral OCF" Conjecture

Combine the ideas: consider the Markov chain on tournaments defined by random arc flips.
The stationary distribution is uniform on all 2^(n choose 2) tournaments on n labeled vertices.

Define f(T) = H(T) and g(T) = I(Omega(T), 2). OCF claims f = g.

Under a single arc flip (i,j):
  delta f = H(T') - H(T)    [our delta_H]
  delta g = I(Omega(T'),2) - I(Omega(T),2)    [our delta_I]

If we could show that delta f = delta g for ALL flips (our THM-015 approach), and both
agree at the transitive tournament (f = g = 1), then f = g everywhere.

The Dual Burnside framework suggests studying this "dual" to the path-counting chain:
instead of walking on tournaments, walk on the arc-flips themselves, with the dual
stationary law weighted by how many tournaments each flip "fixes" (preserves H = I).

### The "Flip Class + OCF" Approach

From Rajkumar et al.: every tournament is in the flip class of an R-cone. An R-cone has
a universal vertex v. For an R-cone, H(T) and I(Omega(T), 2) can both be expressed in
terms of T-v (the base tournament).

Specifically: if v beats everyone, every Hamiltonian path either starts or ends at v (since
v has no in-arcs from the rest). So H(T) = H_start(v) + H_end(v) where these count paths
starting/ending at v. These relate directly to H(T-v) via insertion counts.

For the conflict graph: every odd cycle through v uses exactly one in-arc and one out-arc
of v. Since v beats everyone, the in-arc to v must come from the cycle closing back.
This constrains Omega(T) for R-cones.

STRATEGY: Prove OCF for R-cones (simpler structure), then use the flip-class machinery
to extend to all tournaments. The flip phi_S reverses arcs across a cut — if we can track
how both H and I(Omega, 2) change under phi_S, we get OCF for all T from OCF for R-cones.

### The "Alternating Sum as Burnside Orbit Count"

Our Even-Odd Split Lemma: sum (-1)^|S| Delta(S, R) = 0.

This alternating sum looks like an inclusion-exclusion / Mobius inversion. In the Burnside
framework, alternating sums arise naturally from character theory:

  sum_{g in G} chi(g) |X_g| = |G| * (number of orbits carrying representation chi)

For chi = sign character of S_m (where m = |M| = n-2), this gives:
  sum_{sigma in S_m} sgn(sigma) |X_sigma| = |S_m| * (# orbits with sign = -1)

If we identify our subsets S of M with permutations (or with elements of a group acting on
the "middle vertices"), the alternating sum sum (-1)^|S| Delta(S, M\S) could be an instance
of such a character-weighted Burnside sum.

SPECULATIVE CONNECTION: The group Z_2^m acts on M = V\{i,j} by "flipping" which side
of the path each vertex is on (left part S vs right part R). The character (-1)^|S| is the
determinant character of this group action. Our lemma says the determinant-weighted count
vanishes — meaning the "signed orbit count" is zero.

This would connect our Even-Odd Split to the representation theory of Z_2^m acting on
Hamiltonian path decompositions, with the Dual Burnside process providing the sampling
framework.

---

Source: opus-2026-03-05-S4b
Papers: arXiv:2110.05188, arXiv:2510.25202
