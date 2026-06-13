# Transfer Matrix Symmetry: Web Research Compilation

**Session:** kind-pasteur-2026-03-06-S23 (extended research)
**Goal:** Exhaustive web research on all connections to the transfer matrix symmetry problem

---

## 1. The Grinberg-Stanley Det/Per Formula (arXiv:2412.10572)

The "Revisiting" paper by Grinberg-Stanley establishes a **fundamental decomposition**:

**Proposition 2:**
```
ham(D) = sum_{S subset [n]} det(A_bar[S]) * per(A[S^c])
```

This decomposes Hamiltonian path enumeration into products of determinants of complement submatrices and permanents of original submatrices. This is DIRECTLY related to our transfer matrix:

Our M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)

The E_a and B_b are Hamiltonian path counts in subgraphs — they are PERMANENTS (sums over permutations with all positive signs) of restricted adjacency submatrices. The (-1)^|S| factor introduces the determinantal sign.

**Key connection:** Their walk generating function
```
W_D(z) = det(I + z*X*A_bar) / det(I - z*X*A)
```
encodes ALL path information. The coefficient extraction L_n W_D(1) yields ham(D).

**Why this matters for us:** Our transfer matrix M is essentially the BILINEAR form obtained from this det/per decomposition, where we fix endpoints a,b and split the internal vertices into two groups S, R. The symmetry M[a,b] = M[b,a] is asking whether this bilinear form is symmetric — a question about the interplay between det and per.

## 2. The Rédei-Berge Hopf Algebra (arXiv:2402.07606)

The comultiplication in the Rédei-Berge Hopf algebra is:

**Delta([X]) = sum_{S subset V} [X|_S] tensor [X|_{V\S}]**

This IS our subset convolution! The coalgebra structure directly encodes the decomposition we use for M[a,b].

**Antipode formula (Theorem 4.9):**
```
S(U_X) = (-1)^|V| * U_{X_bar}
```

This connects the Hopf antipode to the complementary digraph, with a sign depending on vertex count. For tournaments, X_bar = X^op (converse), so:
```
S(U_T) = (-1)^n * U_{T^op}
```

**Implication:** The antipode is an INVOLUTION (S^2 = id in a cocommutative Hopf algebra). So S(S(U_T)) = U_T, which gives:
```
(-1)^n U_{T^op^op} = U_T  =>  U_{T^op} = (-1)^n S(U_T)
```
This is closely related to our T^op equivalence M_{T^op} = (-1)^{n-2} M_T.

## 3. The Noncommuting Rédei-Berge Function (arXiv:2504.20968, Mitrovic)

**DELETION-CONTRACTION (Theorem 3.7):**
```
W_X = W_{X\e} - W_{X/e}
```

This property is ABSENT in the commutative case but present in the noncommuting version. This is potentially very powerful for inductive proofs.

**Tournament specialization (Corollary 3.12):**
```
W_X = sum_{sigma in S_V(X), odd cycles} 2^{psi(sigma)} p_{Type(sigma)}
```

All cycles must have odd length; psi(sigma) counts nontrivial cycles.

**Key insight for us:** If we could express M[a,b] as a specialization of W_X evaluated at specific noncommuting variables, the deletion-contraction property might give an inductive proof of symmetry.

## 4. Feng's Dual Burnside Process (arXiv:2510.25202)

**Q = AB, K = BA factorization:**
- A(g,x) = 1{x in X_g}/|X_g| (forward leg)
- B(x,h) = 1{h in G_x}/|G_x| (backward leg)

Our transfer matrix has a similar factorization:
- M = E^T * Lambda * B where Lambda = diag((-1)^|S|)

The Feng paper proves:
- **Spec_nonzero(Q) = Spec_nonzero(K)** (shared eigenvalues)
- **Reversibility:** pi(g)Q(g,h) = pi(h)Q(h,g) (detailed balance)
- **Mixing time equivalence:** |t_mix(Q) - t_mix(K)| <= 1

**Crucial parallel:** Feng's Q is reversible (symmetric under detailed balance) because the AB factorization inherits symmetry from the group action. Our M is symmetric because the E^T Lambda B factorization inherits symmetry from the tournament constraint t_ij + t_ji = c.

The tournament constraint plays the role of "detailed balance" in the Markov chain analogy!

## 5. Björklund's Parity of Directed Hamiltonian Cycles (arXiv:1301.7250)

**Key formula:** A combinatorial formula for ham_cycles(D) mod any positive integer, based on inclusion-exclusion over cycle covers.

**Relevant insight:** The relationship det(A) ≡ per(A) (mod 2) for 0-1 matrices means that over F_2, the determinantal signs become invisible. This is why Rédei's parity theorem works — the (-1) signs don't matter mod 2.

For our problem: at c=0 (pure skew), the "even r-powers" conjecture says that the interplay between det-like ((-1)^|S|) and per-like (path products) contributions has a specific parity structure. Björklund's approach of working modulo primes could potentially resolve this.

## 6. Skew-Adjacency Matrices of Tournaments

For a tournament T, the skew-adjacency matrix S has:
- S_ij = 1 if i -> j, S_ji = -1 if i -> j
- S is skew-symmetric: S + S^T = 0

Properties:
- det(S) = Pf(S)^2 for even n
- det(S) = 0 for odd n
- Principal minors encode subgraph structure

**Connection to our c=0 case:** When c=0, our arc weights ARE the entries of the skew-adjacency matrix (up to scaling). The c=0 proof via path reversal is essentially using the skew-symmetry S = -S^T.

**Bounded principal minors (2023 paper):** Tournaments where all principal minors of S are bounded. "Local orders" have all principal minors at most 1. This connects to the structure of the transfer matrix entries.

## 7. Doubly Regular Tournaments and Skew Hadamard Matrices

DRTs are equivalent to skew Hadamard matrices of order n+1. A DRT on n=4t+3 vertices has:
- Every vertex: in-degree = out-degree = 2t+1
- For any two vertices v,w: exactly t common predecessors

The skew Hadamard matrix H satisfies H + H^T = 2I, i.e., h_ij + h_ji = 2 for off-diagonal entries. This is EXACTLY our c-tournament condition with c=2!

**Spectral characterization:** Tournament eigenvalues are {(n-1)/2, (-1 ± sqrt(n))/2} with specific multiplicities.

**Connection to Paley:** Paley tournaments T_p are DRTs (for prime p ≡ 3 mod 4). We already know T_p maximizes H(T). The DRT = skew Hadamard connection means the Paley maximizer lives in a very structured algebraic class.

## 8. The Power-Sum Expansion and p-Positivity

**Grinberg-Stanley (Theorem 1.38):** For tournaments:
```
U_D = polynomial in p_1, 2p_3, 2p_5, 2p_7, ...
```
with NONNEGATIVE integer coefficients.

This p-positivity is remarkable. It means the symmetric function U_T can be written using only ODD power sums, with even multiplicity (the factor of 2).

**Connection to our even r-powers:** The power-sum expansion uses only odd-indexed p_k. In our r=c/2 expansion, only even powers of r appear. These feel related — both are "parity selection rules" for different bases.

## 9. The Reciprocity Formula

**From the Hopf algebra:**
```
u_X(-m) = (-1)^|V| * u_{X_bar}(m)
```

For tournaments (X_bar = X^op):
```
u_T(-m) = (-1)^n * u_{T^op}(m) = (-1)^n * u_T(m)
```

The last equality uses U_{T^op} = U_T (path reversal symmetry of the symmetric function).

This means u_T(m) has DEFINITE PARITY in m: u_T(-m) = (-1)^n u_T(m). So u_T is even/odd depending on n.

**This is EXACTLY our parity finding!** Our M[a,b] has definite s-parity (-1)^{n-2}, which is the same statement with a shift of 2 (because M involves n-2 internal vertices).

## 10. The Converse Invariance Polynomial (arXiv:2407.17051)

**Digraph polynomial:**
```
P_D(x) = sum_v (1+x)^{d+(v)} (1-x)^{d-(v)}
```

A digraph D is converse invariant iff P_D(x) = P_{-D}(x), which happens iff **all odd coefficients vanish**.

**Theorem 3.2:** P_D(x) = P_{-D}(x) iff all odd coefficients vanish up to max odd degree.

This is ANOTHER instance of the "only even powers" pattern! Converse invariance (D ≅ D^op as unlabeled graphs) corresponds to vanishing of odd polynomial coefficients. Compare with our transfer matrix symmetry corresponding to vanishing of odd r-powers.

---

## Synthesis: The Grand Pattern

Multiple independent threads converge on the same structure:

1. **Transfer matrix symmetry** (our problem): M[a,b]=M[b,a] ⟺ only even r-powers
2. **Rédei-Berge reciprocity**: u_T(-m) = (-1)^n u_T(m) ⟺ definite parity in m
3. **Converse invariance polynomial**: P_D = P_{-D} ⟺ only even x-powers
4. **Hopf antipode**: S(U_T) = (-1)^n U_{T^op} = (-1)^n U_T
5. **p-positivity**: U_T uses only p_1, 2p_3, 2p_5, ... (odd power sums)
6. **DRT/skew Hadamard**: H + H^T = 2I (the c=2 case of our c-tournament)

All of these are manifestations of the same underlying principle: **tournaments have a "parity symmetry" inherited from the constraint t_ij + t_ji = constant.** The skew-symmetric part s_ij = -s_ji provides a natural Z/2 grading (parity), and this grading propagates through all algebraic structures built from the tournament.

The proof of transfer matrix symmetry should exploit this parity structure. The most promising approach appears to be:

1. Express M[a,b] using the Grinberg-Stanley det/per formula
2. Use the tournament constraint to relate det and per contributions
3. Show the odd-parity terms cancel via the skew symmetry

Or alternatively, use the Mitrovic deletion-contraction in the noncommuting setting to prove symmetry inductively.

---

## 11. The Cover Polynomial Reciprocity (Chung-Graham)

For an n-vertex digraph D, the cover polynomial satisfies:
```
C(D; 1, 0) = ham(D)  (Hamiltonian path count)
C(D; 0, 1) = per(A)  (permanent)
C(D; 0, -1) = (-1)^n det(A)  (determinant)
```

**Reciprocity formula:**
```
C(D; x, y) = (-1)^n C(D'; -x-y, y)
```
where D' is the complement digraph.

For tournaments, D' = D^op, so:
```
C(T; x, y) = (-1)^n C(T^op; -x-y, y)
```

At (x,y)=(1,0): C(T;1,0) = (-1)^n C(T^op; -1, 0). But C(T^op; 1, 0) = ham(T^op) = ham(T) by path reversal. So this gives:
```
ham(T) = (-1)^n C(T; -1, 0)
```

This is EXACTLY a specialization of the Rédei-Berge reciprocity!

**Connection to our problem:** The cover polynomial's reciprocity is a global statement. Our transfer matrix M[a,b] is an endpoint-conditioned version. Can the cover polynomial be "refined" to track endpoints? The answer is YES — Chow's path-cycle symmetric function does exactly this, and our M[a,b] corresponds to extracting the coefficient of x_a x_b in a specific specialization.

## 12. The Even Cycle Vanishing Connection (opus-S10)

Opus PROVED: p_mu(U_T) = 0 whenever mu has an even part. The proof uses cycle reversal involution (reverse an even-length cycle, sign changes by (-1)^{k-1} = -1).

**The key bridge question:** Our even-r-powers conjecture says M(r,s) has only even powers of r=c/2. The Even Cycle Vanishing says U_T uses only odd-part cycle types in its power-sum expansion.

These ARE related but not identical:
- U_T is a GLOBAL symmetric function encoding ALL Hamiltonian path information
- M[a,b] is an ENDPOINT-CONDITIONED quantity

**Hypothesis:** M[a,b] can be expressed as a specific coefficient extraction from U_T or a closely related object. If so, the even-cycle-vanishing propagates to give the even-r-powers property.

**Evidence:** The Grinberg-Stanley formula ham(D) = sum_S det(A_bar[S]) per(A[S^c]) is exactly the "total" version of our M[a,b] = sum_S (-1)^|S| E_a(S) B_b(R). The total ham(D) = sum_{a,b} M[a,b] (or a related sum). If M has definite parity in s, the sum also has definite parity, consistent with even-cycle-vanishing.

**Possible proof strategy:** Express M[a,b] using the endpoint-refined Rédei-Berge function, then use the involution argument at the refined level.

## 13. El Sahili's Converse Invariance of Oriented Path Types

El Sahili and Ghazo Hanna (2023): A tournament T and T^op contain the same number of oriented Hamiltonian paths of any given TYPE.

Here "type" = the forward/backward arc sequence along the path. An arc is "forward" if it goes in the direction v_i -> v_{i+1}, "backward" otherwise.

**Connection:** This is a STRONGER version of our U_{T^op} = U_T. It says not just the total count but the TYPE-specific counts are equal for T and T^op.

**For us:** If the forward/backward type is preserved, then endpoint-conditioned counts (paths from a to b of each type) should also be related between T and T^op. This could give the endpoint-level parity structure we need.

## 14. Converse Invariance Polynomial (Ai et al. 2025)

The digraph polynomial P_D(x) = sum_v (1+x)^{d+(v)} (1-x)^{d-(v)} satisfies:
P_D(x) = P_{-D}(x) iff D is converse invariant iff ALL ODD COEFFICIENTS VANISH.

This is another instance of the "only even powers" pattern:
- Converse invariance ↔ vanishing of odd coefficients
- Transfer matrix symmetry ↔ vanishing of odd r-powers
- Even Cycle Vanishing ↔ only odd-part cycle types in p-expansion

ALL of these are manifestations of the T ↔ T^op involution acting on polynomial coefficients.

## 15. The Björklund Determinant Sums

Björklund (2014): Introduces "Labeled Hamiltonicity" and relates it to a "Labeled Cycle Cover Sum" over F_2. The cycle cover sum is evaluated via DETERMINANTS.

**Key technical insight:** Over characteristic 2, per(A) = det(A). This is why Rédei's theorem (ham(T) is odd) can be proved algebraically — the permanent and determinant coincide mod 2.

**For us:** Our transfer matrix M = sum_S (-1)^|S| E_a B_b mixes permanent-like (path products) and determinant-like ((-1)^|S| signs) structures. Over characteristic 2, these collapse, and symmetry becomes trivial. The even-r-powers conjecture is about what happens BEYOND mod-2 — it's a refinement of the mod-2 structure.

## 16. Unimodular Tournaments and Skew-Adjacency Determinants

A tournament is UNIMODULAR if det(S) = 1 where S is the skew-adjacency matrix. Equivalent to all eigenvalues being algebraic units.

Unimodular tournaments include:
- Transitive tournaments (at even n)
- Their switches

A tournament has no "diamonds" (vertex dominating/dominated by 3-cycle) iff it's switching-equivalent to a transitive tournament.

**Connection to transfer matrix:** The skew-adjacency matrix S has det(S) = Pf(S)^2 (even n) or 0 (odd n). Our c=0 proof uses the skew-symmetry S = -S^T. The unimodular condition puts a specific constraint on the "total weight" of the tournament that might constrain M[a,b].

---

## Grand Synthesis: Three Pillars of the Proof

The web research reveals THREE complementary pillars for a proof of transfer matrix symmetry:

### Pillar 1: Involution (proven for global, conjectured for endpoint)
- **Global:** Even Cycle Vanishing (opus-S10, PROVED). Cycle reversal involution cancels even-part cycle types.
- **Endpoint:** Need to show the same involution works at the endpoint-conditioned level.
- **Key tool:** Mitrovic's deletion-contraction for the noncommuting Rédei-Berge function could give an inductive handle.

### Pillar 2: Algebraic Structure
- **Hopf algebra:** Comultiplication = subset convolution = our M[a,b] decomposition.
- **Antipode:** S(U_T) = (-1)^n U_{T^op} = (-1)^n U_T for tournaments.
- **Det/Per formula:** ham(D) = sum_S det(A_bar[S]) per(A[S^c]).
- **Key insight:** M[a,b] is an endpoint-conditioned version of this det/per formula.

### Pillar 3: Parity Structure
- **Cover polynomial reciprocity:** C(T; x,y) = (-1)^n C(T^op; -x-y, y).
- **U_T uses only odd power sums** (p-positivity with only odd-part types).
- **Even-r-powers:** M(r,s) uses only even powers of r.
- **Converse invariance polynomial:** only even coefficients.
- **All are consequences of the T ↔ T^op involution combined with the tournament constraint t_ij + t_ji = c.**

### The Remaining Gap

The gap between "global parity" (Even Cycle Vanishing, proved) and "endpoint parity" (even-r-powers, verified but unproved) is exactly the gap between:
- The symmetric function U_T (which sees all paths at once)
- The transfer matrix entry M[a,b] (which conditions on endpoints)

Bridging this gap requires either:
1. Expressing M[a,b] as a coefficient of U_T or a related object, then propagating the parity
2. Adapting the involution proof directly to the endpoint-conditioned setting
3. Using deletion-contraction (Mitrovic) to prove even-r-powers inductively

## 17. ALGEBRAIC PROOF: M[a,b](-r) = M[b,a](r) (NEW, proved in S23)

**This is a clean algebraic proof of the equivalence between the two symmetry conjectures.**

**Theorem:** For any c-tournament with weights t_ij = r + s_ij (where s_ij = -s_ji, r = c/2):
```
M[a,b](-r, s) = M[b,a](r, s)
```

**Proof (5 steps):**
1. **T(-r) = -T(r)^T:** Since T(-r)[i,j] = -r + s_ij and -T(r)^T[i,j] = -(r + s_ji) = -(r - s_ij) = -r + s_ij. ✓
2. **Path reversal under negated transpose:** E_a(V; -T^T) = (-1)^{|V|-1} B_a(V; T) (reversing a path of k arcs under -T^T gives (-1)^k times the reversed path weight under T).
3. **Similarly:** B_b(V; -T^T) = (-1)^{|V|-1} E_b(V; T).
4. **Substitution:** M[a,b](-r) = sum_S (-1)^|S| · [(-1)^|S| B_a(S∪{a}; T)] · [(-1)^|R| E_b(R∪{b}; T)] = sum_S (-1)^{|R|} B_a(S∪{a}) E_b(R∪{b}).
5. **Relabel S↔R:** This equals sum_{S'} (-1)^{|S'|} E_b(S'∪{b}) B_a(R'∪{a}) = M[b,a](r). ✓

**Corollary:** The following are equivalent:
- (i) M[a,b] = M[b,a] (transfer matrix symmetry)
- (ii) M[a,b](r) has only even powers of r
- (iii) M[a,b](r,-s) = (-1)^{n-2} M[a,b](r,s) (s-parity)

Proof: (i)⟺(ii) from M[a,b](-r) = M[b,a](r): M[a,b] = M[b,a] iff M[a,b](r) = M[a,b](-r) iff only even r-powers.

**Verified computationally:** n=4,5,6 in `04-computation/m_negr_equals_m_swap.py`.

## 18. Noncommutative Rédei-Berge: Key Properties (arXiv:2504.20968)

Mitrovic's noncommutative W_X has:
- **W_X = W_{X^op}** (Theorem 3.4) — the same T↔T^op symmetry
- **Deletion-contraction:** W_X = W_{X\e} - W_{X/e}↑ (powerful induction tool)
- **p-expansion:** W_X = sum_{σ} (-1)^{φ(σ)} p_{Type(σ)} (same as commutative)
- **Multiplicativity:** W_{X·Y} = W_X · W_Y for certain graph operations

The noncommutative version tracks vertex positions but does NOT directly encode endpoint-specific path counts. Extracting M[a,b] from W_X remains an open problem.

## 19. Converse Invariance Polynomial (arXiv:2407.17051, Ai et al. 2025)

For an oriented graph D, define the "converse invariance polynomial":
```
P_D(x) = sum_u (1+x)^{d+(u)} * (1-x)^{d-(u)}
```
**Theorem (Ai):** If D is converse invariant (same count in T and T^op for all tournaments T), then P_D(x) = P_{-D}(x), which requires ALL ODD COEFFICIENTS to vanish.

**Parallel to our problem:** Our even-r-powers conjecture also requires odd coefficients to vanish. The Ai polynomial measures degree-sequence compatibility; our polynomial measures path-weight compatibility. Both capture T↔T^op symmetry through parity constraints.

El Sahili proved all oriented graphs with max degree ≤ 2 (including Hamiltonian paths of any type) are converse invariant. But this is about UNWEIGHTED counts and doesn't directly give our weighted M[a,b] symmetry.

## 20. Walk Generating Function and Endpoints (arXiv:2412.10572)

The walk generating function W_D(z) = det(I + zXA^T) / det(I - zXA) for tournaments.

The (a,b) entry of (I - zXA)^{-1} generates walks from a to b. Extracting the multilinear coefficient gives Hamiltonian walk counts. For c-tournaments:
- A = rJ' + S (where J' = J-I, S skew-symmetric)
- The walk GF becomes det(I + zX(rJ'-S)) / det(I - zX(rJ'+S))

Under r → -r: this transforms via the T(-r) = -T^T identity, potentially connecting to the M(-r) = M(b,a) result.

## Open Connections Still to Explore

- [ ] The LGV lemma for non-intersecting paths — can M[a,b] be expressed as a determinant of a path matrix?
- [ ] Björklund's modular formula for ham_cycles — does it extend to give our even-r-powers?
- [ ] The relationship between M[a,b] and immanants (character-weighted permanents)
- [ ] Whether the Cauchy-Binet decomposition M = E^T Lambda B can be symmetrized using the tournament constraint
- [ ] The critical group of tournaments and its connection to the skew-adjacency matrix
- [ ] Savchenko's exact cycle formulas for regular tournaments — do they imply transfer matrix symmetry?
- [ ] The cover polynomial refined for endpoints — does it give M[a,b]?
- [ ] El Sahili's type-specific converse invariance — does it extend to endpoint-conditioned counts?
- [ ] Chow's path-cycle symmetric function — extract M[a,b] from it?
- [ ] The "matrix cover polynomial" (Chung-Graham) — does it encode M[a,b]?

## Sources

- [Revisiting Rédei-Berge via Matrix Algebra](https://arxiv.org/html/2412.10572)
- [Grinberg-Stanley: Rédei-Berge symmetric function](https://arxiv.org/abs/2307.05569)
- [Rédei-Berge Hopf algebra](https://arxiv.org/abs/2402.07606)
- [Mitrovic: Noncommuting Rédei-Berge](https://arxiv.org/abs/2504.20968)
- [Feng: Dual Burnside Process](https://arxiv.org/abs/2510.25202)
- [Björklund-Husfeldt: Parity of Hamiltonian cycles](https://arxiv.org/abs/1301.7250)
- [Björklund: Determinant sums for Hamiltonicity](https://arxiv.org/abs/1008.0541)
- [Ai et al: Converse invariant digraph polynomials](https://arxiv.org/abs/2407.17051)
- [El Sahili: Oriented Hamiltonian paths in tournaments](https://arxiv.org/abs/2101.00713)
- [Rédei revisited](https://arxiv.org/abs/2510.10659)
- [Properties of Rédei-Berge and Hopf algebras](https://arxiv.org/html/2407.18608)
- [Skew-adjacency matrices bounded principal minors](https://www.sciencedirect.com/science/article/abs/pii/S0012365X23002388)
- [DRT = Skew Hadamard](https://arxiv.org/abs/1202.5374)
- [Hamilton transversals](https://arxiv.org/abs/2307.00912)
- [Chow: Path-cycle symmetric function](https://timothychow.net/pathcycle.pdf)
- [Chung-Graham: Matrix cover polynomial](https://mathweb.ucsd.edu/~ronspubs/16_04_cover.pdf)
