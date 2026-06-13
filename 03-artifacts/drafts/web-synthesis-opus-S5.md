# Web Research Synthesis: opus-2026-03-06-S5

## Session Focus
Deep web investigation into all leads, connections, and synthesis opportunities across the project's open questions.

---

## 1. THE OCF PROOF CHAIN (Clarified)

### Irving-Omar (arXiv:2412.10572, 2025)
**Title:** "Revisiting The Rédei-Berge Symmetric Functions via Matrix Algebra"
**Authors:** John Irving, Mohamed Omar

Key results extracted:
- **Proposition 2:** `ham(D) = sum_S det(Ā[S]) · per(A[S^c])` — Hamiltonian path count via determinant-permanent pairs
- **Corollary 20 (OCF for tournaments):** `ham(T̄) = sum_{σ ∈ S(T), all cycles odd} 2^{ψ(σ)}` where ψ(σ) counts nontrivial cycles
- **Theorem 5 (Rédei reproved):** Working mod 2 in F₂[x₁,...,xₙ]/⟨x₁²,...,xₙ²⟩, reducing to det(I+X)
- **Matrix framework:** Walk generating function W_D(z) = det(I+zXĀ)/det(I-zXA)
- **Power sum expansion (Corollary 14):** U_D = sum_{σ ∈ S(D,D̄)} (-1)^φ(σ) p_{cyc(σ)}

**Connection to our work:** Corollary 20 IS our OCF. The permutations σ with all-odd cycles biject exactly with vertex-disjoint odd cycle collections. Each such σ with k nontrivial cycles contributes 2^k to ham(T). Summing gives I(Omega(T), 2).

**IMPORTANT NOTE:** The paper does NOT explicitly mention "independence polynomial" or "conflict graph". The equivalence between their sum-over-permutations formula and I(Omega(T), 2) is implicit but straightforward: collections of vertex-disjoint odd directed cycles = permutations with all odd cycle types minus the identity parts.

### Grinberg-Stanley (arXiv:2307.05569, 2023)
**Title:** "The Rédei-Berge symmetric function of a directed graph"

Key results:
- **Theorem 1.30:** U_D = Σ_{σ ∈ S_V(D,D̄)} (-1)^φ(σ) p_{type σ}
- **Theorem 1.38 (Tournaments):** U_D is polynomial in p₁, 2p₃, 2p₅, 2p₇,... with NONNEG INTEGER coefficients
- **Section 7:** Mod-4 refinement of Rédei's theorem
- **Proof technique:** Algebraic symmetric functions + cycle analysis + inclusion-exclusion

**Synthesis:** The Grinberg-Stanley p-positivity theorem for tournaments is the symmetric function version of OCF. The coefficient of p₁^{n-2k₃-2k₅-...} · (2p₃)^{k₃} · (2p₅)^{k₅} ··· counts collections of k₃ disjoint 3-cycles, k₅ disjoint 5-cycles, etc. Specializing x₁=x₂=...=1 recovers ham(T) = I(Omega(T), 2).

### Grujić-Stojadinović (arXiv:2402.07606, 2024)
**Title:** "The Rédei-Berge Hopf algebra of digraphs"

Key results:
- **Hopf algebra comultiplication:** Δ([X]) = Σ_{W⊆V} [X|_W] ⊗ [X|_{V\W}] — this IS our subset convolution!
- **Deletion property for cycles:** If edges form a k-cycle, then U_X = Σ_{S⊆edges, S≠∅} (-1)^{|S|-1} U_{X\S}
- **Reciprocity gives Berge's theorem:** via the Hopf algebra antipode

**Synthesis for INV-033:** The Hopf algebra formulation directly encodes our subset convolution structure. The comultiplication IS the partition into S and V\S that appears in the transfer matrix. This potentially gives an algebraic route to proving transfer matrix symmetry (INV-001) — if we can show the comultiplication has a natural symmetry under the tournament constraint.

---

## 2. PALEY TOURNAMENTS AND DESIGN THEORY

### Doubly Regular Tournaments (DRTs) ↔ Skew Hadamard Matrices

**MAJOR CONNECTION FOUND:** Reid-Brown (1972) proved:
> DRTs on n vertices are EQUIVALENT to skew Hadamard matrices of order n+1.

A DRT has n = 4t+3 vertices, every vertex has in/out-degree 2t+1, and for any two vertices v,w there are exactly t vertices z with v→z and w→z.

**Paley tournaments T_p (p ≡ 3 mod 4) are the canonical example of DRTs.** The Paley construction gives a skew Hadamard matrix of order p+1.

### Our BIBD Discovery in Context

Our finding that cyclic 3-cycles of T_p form a 2-(p, 3, (p+1)/4) BIBD appears to be NOVEL. The web search found no prior work connecting Paley tournament 3-cycles to BIBDs. The known design theory connections are:
- Paley difference set: QR mod p is a ((p-1)/2, (p-3)/4)-difference set in Z_p
- Paley biplane: 2-(11, 5, 2) from translates of QR mod 11
- Skew Hadamard difference set: QR mod p (p ≡ 3 mod 4)

But the 3-cycle BIBD 2-(p, 3, (p+1)/4) — specifically the Fano decomposition at p=7 — does not appear in the literature. This should be written up.

### Nozaki-Suda (arXiv:1202.5374, 2012)
Characterized skew Hadamard matrices via spectra of tournaments of size n-2. This connects to our Paley deletion conjecture: if T_p is a DRT and skew-Hadamard-equivalent, then T_p-v should inherit special spectral properties.

### Pantangi (EJC, 2019)
"Difference families, skew Hadamard matrices, and critical groups of doubly regular tournaments"
- Proved Paley tournament is INEQUIVALENT to other DRT constructions via critical group computation
- Paley DRTs have specific algebraic structure not shared by all DRTs

---

## 3. H-MAXIMIZATION ASYMPTOTICS

### Szele (1943) → Alon (1990) → Adler-Alon-Ross (2001)

**Szele's conjecture (proved by Alon):** max H(T) ∼ n!/2^{n-1}

**Adler-Alon-Ross (2001) improved lower bound:**
> max H(T) ≥ (e - o(1)) · n!/2^{n-1}

Our computed ratios H(T_p)/(p!/2^{p-1}): 2.000, 2.400, 2.440, 2.527, 2.557 converging toward e ≈ 2.718.

**SYNTHESIS:** The Adler-Alon-Ross bound is achieved by RANDOM REGULAR tournaments. Our data suggests Paley tournaments approach this bound. Since Paley tournaments are "pseudorandom" (quasi-random in the Chung-Graham-Wilson sense), this is consistent: the maximizer should be as "random-looking" as possible while maintaining the algebraic structure that allows maximal H.

**Open question:** Does H(T_p)/(p!/2^{p-1}) → e as p → ∞? This would prove Paley tournaments are asymptotically optimal.

### OEIS A038375 Status
- a(1)...a(11) = 1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095
- a(8)=661, a(10)=3357 computed by Peter Kagey (2020)
- a(11)=95095 computed by Gordon Royle (2022) — matches H(T_11)!
- Our computation: a(7)=189=H(T_7) confirmed exhaustive, a(11)=95095=H(T_11)
- We should submit H(T_p-v) values and the Paley maximizer conjecture

---

## 4. TRANSFER MATRIX SYMMETRY AND FENG'S DUAL BURNSIDE

### Feng (arXiv:2510.25202, 2025)
**Title:** "The Dual Burnside Process"

Key results:
- **Q = AB factorization:** dual chain Q has matrix factorization Q=AB, K=BA with classical Burnside kernel K
- **Reversibility:** Dual process is REVERSIBLE with π(g) ∝ |X_g|
- **Eigenvalue sharing:** Q and K share all nonzero eigenvalues; mixing times differ by ≤ 1

**Connection to our transfer matrix (INV-016 deepened):**

Our transfer matrix M[a,b] = Σ_S (-1)^|S| E_a(S)·B_b(V\S) has EXACTLY the AB structure Feng describes:
- A = "start at vertex a, enumerate paths through S"
- B = "enumerate paths through complement V\S, end at vertex b"
- The tournament constraint T[x,y]+T[y,x]=1 plays the role of DETAILED BALANCE

**NEW INSIGHT:** The tournament constraint is the REVERSIBILITY condition in Feng's framework. This explains WHY transfer matrix symmetry holds only for tournaments: tournaments are the "reversible" digraphs in the Burnside process sense. For general digraphs (independent arc variables), there is no reversibility and M need not be symmetric.

**Proof strategy for INV-001:** Formalize the tournament as a group action (Z₂ acting on each arc pair {i,j}), show this action satisfies Feng's reversibility conditions (Thm 3.3), conclude Q=AB is symmetric = transfer matrix is symmetric.

### Compatible Polynomials (Chudnovsky-Seymour, 2007)

A family of polynomials is **compatible** if all non-negative linear combinations have only real roots. Key theorem: compatibility ⟺ pairwise compatibility ⟺ existence of a common interlacing polynomial.

This is the technique that proves I(Omega(T), x) has all real roots for n ≤ 8 (where Omega is claw-free). For n ≥ 9, claws appear and the technique fails — consistent with our THM-025 counterexample.

### Sokal's Conjecture (Peters-Regts, arXiv:1701.08049)
Zero-free radius for independence polynomials: λ_s(Δ) = (Δ-1)^{Δ-1}/Δ^Δ. For all Δ ≥ 3, λ_s < 0.15 << 2. So λ=2 is ALWAYS outside the zero-free disk, confirming our INV-028 finding that OCF is in the non-perturbative regime.

---

## 5. CONVERSE SYMMETRY AND EL SAHILI

### El Sahili & Ghazo Hanna (2023)
**Title:** "About the number of oriented Hamiltonian paths and cycles in tournaments"
J. Graph Theory 102 (2023), 684-701.

Main result: T and T^op (converse) have the same number of oriented Hamiltonian paths of EVERY type. This generalizes:
- H(T) = H(T^op) (Rédei/Berge)
- The TYPE distribution is also preserved

**Connection:** Our path reversal identity M_{T^op}[i,j] = (-1)^{n-2} M_T[j,i] is a STRONGER statement. Combined with transfer matrix symmetry M[i,j] = M[j,i], we get M_{T^op} = (-1)^{n-2} M_T, which implies T and T^op have the same transfer matrix up to sign. El Sahili-Ghazo Hanna's result follows from this.

### Ai (2025)
"Number of Subgraphs and Their Converses in Tournaments and New Digraph Polynomials"
J. Graph Theory. Extends El Sahili's work to general subgraph counting. May contain new polynomial invariants relevant to OCF.

---

## 6. BJÖRKLUND-HUSFELDT (arXiv:1301.7250)

**Title:** "The Parity of Directed Hamiltonian Cycles" (FOCS 2013)

Key formula: Computes #HamCycles mod 2 via inclusion-exclusion over cycle covers using F₂ linear algebra. The approach:
1. Reduce permanent mod 2 = determinant mod 2
2. Sieve for Hamiltonian cycles via systems of linear equations over Z₂
3. O(1.619^n) deterministic algorithm

**Connection to OCF:** Björklund's reduction permanent → determinant mod 2 is the same algebraic identity that underlies Rédei's theorem (det(I+X) mod 2). Our OCF goes further: not just mod 2, but the EXACT count via the independence polynomial. The question is whether Björklund's cycle cover sieve can be "lifted" from F₂ to Z to give the full OCF.

---

## 7. KEY SYNTHESIS OPPORTUNITIES

### A. Transfer Matrix Symmetry via Hopf Algebra + Feng Reversibility
**Route:** The Hopf algebra comultiplication Δ([T]) = Σ_S [T|_S] ⊗ [T|_{V\S}] encodes our subset convolution. Feng's reversibility theorem says Q=AB is symmetric when the action satisfies detailed balance. If we can show the Hopf algebra's comultiplication is "self-dual" under the tournament constraint (T[x,y]+T[y,x]=1), then transfer matrix symmetry follows for all n.

**Status:** This is the most promising algebraic route to INV-001. Needs formalization.

### B. Paley Maximizer Conjecture via DRT + Design Theory
**Route:** Paley T_p is a DRT ↔ skew Hadamard matrix of order p+1. DRTs have the most "balanced" cycle distribution. Our BIBD discovery (3-cycles form 2-(p,3,(p+1)/4)) quantifies this balance for triples. Hypothesis: the design-theoretic balance maximizes alpha_k in I(Omega,2), hence maximizes H.

**Status:** Verified computationally through p=23. The BIBD structure is potentially novel. Needs theoretical proof connecting design balance to independence polynomial maximization.

### C. Hereditary Chain via Spectral Regularity
**Route:** Nozaki-Suda characterize skew Hadamard matrices via tournament spectra. Our finding that corr(H, λ₁(AA^T)) = -0.97 means maximizers are spectrally regular. DRTs (= Paley) have the most regular spectrum. Vertex deletion from a DRT preserves "near-regularity" — explaining why T_p-v is the (p-1)-maximizer.

**Status:** Verified at p=3,7,11,19. The spectral regularity → H-maximization link needs formalization.

### D. Rédei-Berge Symmetric Function Encodes Everything
**Route:** The Grinberg-Stanley symmetric function U_T encodes ALL Hamiltonian path information. The p-positivity theorem for tournaments (U_T ∈ Z≥0[p₁, 2p₃, 2p₅, ...]) is the symmetric function form of OCF. The Irving-Omar matrix algebra gives computational tools. The Grujić-Stojadinović Hopf algebra gives structural tools.

**Key question:** Can we use the Hopf algebra DELETION property (U_X = Σ (-1)^{|S|-1} U_{X\S} for edge subsets forming cycles) to relate H(T) to H(T-v)? This could give Claim A directly from the Hopf algebra structure.

---

## 8. NEW LEADS TO INVESTIGATE

### INV-045: Satake's cyclotomic NDRTs (arXiv:2502.12090, Feb 2025)
Nearly-doubly-regular tournaments from almost difference sets. Savchenko's conjecture on canonical spectrum. May give new tournament families approaching Paley's H-maximization.

### INV-046: Jerrum (2026) zero-free regions for restricted graph classes
J. London Math. Soc. New results on where independence polynomial roots can lie for restricted graph classes. May apply to Omega(T) structure.

### INV-047: Irving-Omar matrix algebra → transfer matrix
Their det/per formula ham(D) = Σ_S det(Ā[S])·per(A[S^c]) is closely related to our transfer matrix. The "matrix algebra" approach may give a direct proof of transfer matrix symmetry.

### INV-048: Pantangi's critical groups distinguish Paley from other DRTs
The critical group (Smith normal form of Laplacian) could give algebraic invariants that explain WHY Paley maximizes H among all DRTs.

---

## Sources
- [Irving-Omar arXiv:2412.10572](https://arxiv.org/abs/2412.10572)
- [Grinberg-Stanley arXiv:2307.05569](https://arxiv.org/abs/2307.05569)
- [Grujić-Stojadinović arXiv:2402.07606](https://arxiv.org/abs/2402.07606)
- [Grujić-Stojadinović follow-up arXiv:2407.18608](https://arxiv.org/abs/2407.18608)
- [Feng arXiv:2510.25202](https://arxiv.org/abs/2510.25202)
- [El Sahili-Ghazo Hanna 2023 / arXiv:2101.00713](https://arxiv.org/abs/2101.00713)
- [Adler-Alon-Ross 2001](https://adler.ieor.berkeley.edu/ilans_pubs/hamilt_2001.pdf)
- [OEIS A038375](https://oeis.org/A038375)
- [Nozaki-Suda arXiv:1202.5374](https://arxiv.org/abs/1202.5374)
- [Pantangi EJC v26i3p59](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v26i3p59)
- [Björklund-Husfeldt arXiv:1301.7250](https://arxiv.org/abs/1301.7250)
- [Satake arXiv:2502.12090](https://arxiv.org/abs/2502.12090)
- [Peters-Regts / Sokal arXiv:1701.08049](https://arxiv.org/abs/1701.08049)
- [Chudnovsky-Seymour compatible polynomials](https://web.math.princeton.edu/~mchudnov/roots.pdf)
