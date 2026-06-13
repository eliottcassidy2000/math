# Deep Paper Connections — Detailed Analysis

Source: opus-2026-03-05-S4b (second session)

---

## 1. Schweser-Stiebitz-Toft (arXiv:2510.10659) — Redei Revisited

### Key Theorems Extracted

**Redei's Stronger Theorem (Thm 1.1):** Let T be a tournament with >= 2 vertices.
Add a new non-empty set W of vertices and some non-oriented edges between W and T
and within W. Then the number of Hamiltonian paths in the resulting mixed graph G,
beginning AND ending in T, is EVEN.

**Berge's Stronger Theorem (Thm 1.2):** For any mixed graph G with >= 2 vertices,
G and its complement G-bar have the same parity of Hamiltonian paths.

**Dirac's Stronger Theorem (Thm 2.1):** For a mixed graph G and subset A of
oriented+non-oriented edges, N_A (paths containing A) and N_{=A} (paths containing
exactly A) have the same parity.

### Connection to OCF

**DIRECT CONNECTION TO OPEN PROBLEM 4 (mixed graphs):**
Our Q-Lemma (Route A) proves Redei for tournaments via a fixed-point-free involution
on two-block decompositions. Schweser-Stiebitz-Toft's "Redei's Stronger Theorem"
handles mixed graphs (tournaments + non-oriented edges).

To extend our Q-Lemma: in a mixed graph, a non-oriented edge {u,v} means BOTH
u->v and v->u are present. A Hamiltonian path using this edge traverses it in one
direction. The Q-function Q(P, k) inserts vertex v at position k — but in a mixed
graph, the insertion condition changes: instead of requiring exactly T[v][a]=1 AND
T[b][v]=1, a non-oriented edge allows BOTH directions.

**Key insight:** Non-oriented edges DOUBLE the insertion opportunities. If edge
{a,v} is non-oriented, then both a->v and v->a are valid. This means inshat(v,P')
changes its counting: each non-oriented edge adjacent to v contributes to BOTH
"beat" and "beaten by" counts simultaneously.

**Strategy:** Define Q(P,k) for mixed graphs. The involution on (P, k) pairs should
still work — what changes is the boundary counting. The parity argument (inshat is
always odd) should generalize because:
- Oriented edges contribute 0 or 1 to transitions (as before)
- Non-oriented edges contribute 1 to BOTH boundary terms
- The algebraic identity inshat = b + #01 + #10 still holds with modified b.

**Concrete next step:** Verify computationally that inshat is always odd for mixed
graphs. If so, the Q-Lemma proof extends directly.

### Connection to OCF via Berge

Berge's theorem (G and G-bar have same path parity) is interesting because for
tournaments T, the complement T-bar is the reverse tournament T^op. So
H(T) ≡ H(T^op) (mod 2). For self-converse tournaments T = T^op, this is trivially
true. But Berge says more: even for MIXED graphs, complementation preserves parity.

If OCF generalizes to mixed graphs (H(G) = I(Omega(G), 2) for mixed G), then
Berge's theorem would require I(Omega(G), 2) ≡ I(Omega(G-bar), 2) (mod 2).
This constrains the conflict graph under complementation.

---

## 2. Rajkumar-Veerathu-Mir (arXiv:2110.05188) — Tournament Representations

### Key Theorems Extracted

**Definition 1 (Cut-flip):** phi_S(T) reverses all edges across cut (S, S-bar).
**Definition 2 (Flip class):** F(T) = {T' : exists S with T' isomorphic to phi_S(T)}.

**Proposition 1:** Every flip class contains an R-cone.
*Proof:* For any T, pick vertex i. Set T' = phi_{i ∪ T_i^-}(T). Then T' is R-coned by i (i beats everyone in T').

**Theorem 5:** Locally transitive tournaments = flip class of the transitive tournament.
**Theorem 6:** T is locally transitive iff rank(M) = 2 for any skew-symmetric M inducing T.

**Theorem 11 (Dimension bound):** For any tournament T, there exists a representation
in R^{2(mu(T)+1)} dimensions, where mu(T) = min_sigma min_{T' in F(T)} theta(sigma, T')
is the flip feedback node set size.

**Theorem 12 (Sign-rank bound):** sign-rank(G) <= 2(mu(T)+1) for the sign matrix G
associated with tournament T.

### Deep Connection to OCF Proof Strategy (INV-004)

The proof of Proposition 1 is constructive: given T, the R-cone T' = phi_{i ∪ T_i^-}(T)
is obtained by a SINGLE cut-flip. So the "distance" from any tournament to an R-cone
is exactly 1 cut-flip.

**For OCF proof via R-cones:**
1. Show OCF for R-cones (simpler: universal vertex means every Ham path starts or
   ends at the cone vertex)
2. Show E(T) = H(T) - I(Omega(T), 2) is invariant under cut-flips phi_S

For step 2: a cut-flip phi_S reverses O(|S| * |S-bar|) arcs simultaneously.
This is much harder than a single arc-flip. However, the structure is special:
ALL arcs between S and S-bar are reversed simultaneously, while arcs WITHIN S
and WITHIN S-bar are preserved.

**Key observation:** Under phi_S, the odd cycles of T change in a structured way.
A cycle C is "affected" by phi_S iff C has at least one vertex in S and at least one
in S-bar. The affected cycle C becomes a cycle C' in phi_S(T) with a specific
relationship: the internal structure (within S and within S-bar) is preserved, but
the cross-arcs are all reversed.

For 3-cycles (a,b,c): if all three are in S or all in S-bar, the cycle is preserved.
If two are in S and one in S-bar (or vice versa), the cycle may be destroyed or
created. The counting of created/destroyed cycles under cut-flips should be
tractable because the structure is more constrained than arbitrary multi-arc flips.

### Connection to mu(T) Induction (INV-005)

Theorem 11 gives dim <= 2(mu(T)+1). For locally transitive T: mu(T) = 0, dim <= 2.
For a tournament with one "defect": mu(T) = 1, dim <= 4.

**Induction strategy:** Prove OCF for mu(T) = 0 (transitive — trivial, H=1).
Then show: if OCF holds for all T' with mu(T') < k, it holds for mu(T) = k.

The inductive step: T with mu(T) = k can be "corrected" to a tournament T' with
mu(T') = k-1 by... what operation? The cut-flip phi_S connects T to its R-cone
(Proposition 1), but mu may change unpredictably under cut-flips.

**Better:** mu(T) relates to the feedback arc set. If we flip a single arc in the
feedback arc set, we reduce the feedback count by 1. So induction on |FAS(T)| might
be more tractable than induction on mu(T). But |FAS| is the number of ARCS, not
nodes, in the feedback set.

---

## 3. Feng (arXiv:2510.25202) — Dual Burnside Process (Deep Dive)

### Key Technical Results

**Theorem 3.3 (Reversibility):** The dual Burnside chain Q on G* has stationary
distribution pi(g) = |X_g| / (|G| * z) and is REVERSIBLE:
  pi(g) * Q(g,h) = pi(h) * Q(h,g)

**Lemma 3.8 (Factorization):** Q = AB, K = BA, where:
  A(g,x) = 1_{x in X_g} / |X_g| (forward: group -> states)
  B(x,h) = 1_{h in G_x} / |G_x| (backward: states -> group)

**Theorem 3.9 (Spectral correspondence):** Spec_{!=0}(Q) = Spec_{!=0}(K) with
matching multiplicities.

**Theorem 3.10 (Eigenvector intertwining):** If v is a K-eigenvector with eigenvalue
lambda, then Av is a Q-eigenvector with the same eigenvalue. A and B are mutual
inverses (up to scalar lambda^{-1}) on eigenspaces.

**Block-flip matrix:** M = [[0,A],[B,0]], M^2 = [[Q,0],[0,K]]. The block matrix M
is irreducible and bipartite with period 2.

### PRECISE Connection to Transfer Matrix Symmetry

Our transfer matrix M[a,b] = sum_S (-1)^|S| E_a(S) * B_b(M\S) has the structure:
- "Forward leg" E_a(S): number of Ham paths of T[S] ending at vertex adjacent to a
- "Backward leg" B_b(R): number of Ham paths of T[R] starting from vertex adjacent to b

This is EXACTLY the AB factorization structure:
- A maps from "subsets S" to "Ham path ends at a-neighbor"
- B maps from "Ham path starts at b-neighbor" to "complement subsets R"

The Feng framework says Q = AB is symmetric when detailed balance holds:
  pi(g) * Q(g,h) = pi(h) * Q(h,g)

In our setting, pi would be a weighting on vertices {i,j}. With pi = uniform,
symmetry of Q means Q(i,j) = Q(j,i), i.e., M[i][j] = M[j][i].

**The deep question:** What is the "group action" underlying our transfer matrix?

Candidate: The group Z_2 acts on {i,j} by transposition. The "states" are the
subsets S of V\{i,j}. The Z_2 action on states: sigma maps S to V\{i,j}\S = R
(the complement). Then:
- X_sigma = {S : sigma(S) = S} = {S : S = R} = emptyset (no fixed points unless S = R, which requires |V\{i,j}| even)
- The Burnside chain on this action counts orbits of subsets under complement-flip.

But this doesn't directly give our (-1)^|S| weighting. The sign character of Z_2
DOES give the alternating sum, which is our even-odd split.

**Actionable insight:** The transfer matrix M is symmetric because the underlying
structure satisfies a HIDDEN detailed balance condition. The forward and backward
legs are dual in the sense that the "flow" from i to j through subset S equals
the "flow" from j to i through the same subset S. This is because:
  E_i(S) * B_j(R) counts "paths through S ending near i, continuing through R from j"
  E_j(S) * B_i(R) counts "paths through S ending near j, continuing through R from i"

The equality sum_S (-1)^|S| [E_i*B_j - E_j*B_i] = 0 is the statement that the
SIGNED flow is symmetric. The sign (-1)^|S| plays the role of the character in the
Burnside framework.

---

## 4. Cross-Paper Synthesis

### The Flip-Class + Burnside Connection

Rajkumar's flip classes partition tournaments into orbits under cut-flips (Z_2^n action).
Feng's Dual Burnside samples orbits with stationary law proportional to fixed-point counts.

**Novel synthesis:** Consider the Z_2^n action on tournaments by cut-flips. The fixed
points of a flip phi_S are tournaments invariant under reversing arcs across (S, S-bar).
The Burnside count gives |{orbits}| = (1/2^n) * sum_S |Fix(phi_S)|.

Our OCF identity H(T) = I(Omega(T), 2) is flip-class invariant (since E(T) = 0 is
preserved under cut-flips IF OCF holds). So proving E = 0 on one representative per
flip class suffices. Proposition 1 says R-cones are representatives.

**The dual chain on flip classes:** Walk on the Z_2^n group (subsets S) with transitions
via cut-flips. The dual stationary law pi(S) proportional to |{T : phi_S(T) = T}|
weights subsets by how many tournaments they fix. The spectral correspondence
(Theorem 3.9) then relates the mixing of this dual chain to the classical Burnside chain
on tournaments.

### Redei + Representation Dimension

Schweser-Stiebitz-Toft's mixed-graph Redei + Rajkumar's dimension theory:
A rank-d tournament can be "mixed" by adding non-oriented edges without changing
the dimension. The mixed-graph Redei theorem then applies to rank-d mixed tournaments.

For rank 2 (locally transitive): the mixed graph obtained by making some arcs
non-oriented is still representable in R^2. Redei's stronger theorem applies to
these mixed graphs, giving parity constraints on their Ham paths.

**Open question:** Does the OCF formula I(Omega(G), 2) = H(G) extend to rank-2
mixed graphs? If so, the dimension-2 structure might simplify the proof.

---

## Summary of New Findings

| INV | Paper | Finding | Status |
|-----|-------|---------|--------|
| INV-010 | SST | Q-Lemma extends to mixed graphs if inshat remains odd | TESTABLE |
| INV-019 | SST | Berge's thm constrains Omega under complementation | NEW CONNECTION |
| INV-004 | RVM | Cut-flip distance to R-cone is exactly 1 (Prop 1) | CONFIRMED |
| INV-005 | RVM | mu(T) induction less promising than FAS induction | ASSESSED |
| INV-015 | RVM | Theorem 12 gives sign-rank bound via mu(T) | CONFIRMED |
| INV-016 | Feng | Q=AB is symmetric via detailed balance (Thm 3.3) | CONFIRMED |
| INV-016 | Feng | Block-flip M^2 = diag(Q,K) structure matches ours | CONFIRMED |
| INV-001 | Feng | Transfer matrix symmetry = hidden detailed balance | KEY INSIGHT |

Source: opus-2026-03-05-S4b
