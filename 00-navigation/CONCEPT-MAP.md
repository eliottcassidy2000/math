# Concept Map — Tournament Parity Research

**Purpose:** Complete, structured database of every mathematical concept, object, technique, and connection in this project. Organized for rapid lookup by future Claude instances. Created by kind-pasteur-2026-03-07-S34.

**Last updated:** kind-pasteur-2026-03-07-S34

---

## I. CORE MATHEMATICAL OBJECTS

### Tournaments
| Concept | Definition | Key Properties | Where Used |
|---------|-----------|---------------|------------|
| **Tournament** T | Complete directed graph on n vertices | Every pair has exactly one arc | Everywhere |
| **Opposite tournament** T^op | Reverse all arcs: T^op(u,v) = T(v,u) | H(T) = H(T^op) for Hamiltonian paths | THM-002, THM-030 |
| **Self-complementary (SC) tournament** | T ≅ T^op via some permutation | Exists only at n ≡ 0,1 (mod 4); |Aut| always odd (Moon) | THM-024, INV-043, T019 |
| **Self-converse tournament** | T has an anti-automorphism | SC ⊂ Self-converse; all circulant tournaments are SC | THM-052 |
| **Paley tournament** T_p | p ≡ 3 (mod 4) prime; i→j iff j−i is QR mod p | Vertex-transitive, self-complementary, max H(T) | T019, T053, INV-042 |
| **Regular tournament** | Every vertex has out-degree (n−1)/2 | Only at odd n; includes all Paley | THM-027, BIBD |
| **Doubly regular tournament (DRT)** | Regular + every pair of vertices has same # common out-neighbors | Unique at n=3,5,7; multiple at n≥11 | INV-068 |
| **Transitive tournament** | Acyclic; unique up to labeling | H = 1; base case for induction | OCF base case |
| **Cycle-rich tournament** | Every vertex in a directed 3-cycle | ⟹ no source/sink (Lemma Q); key for H=21 proof | THM-079, Part R |
| **R-cone** | One vertex beats/loses to all others | Rajkumar et al.; every tournament is flip-class of R-cone | T046, INV-015 |
| **Locally transitive tournament** | rank 2; sub-neighborhoods are transitive | Still has 5-cycles (67% at n=5) | T046 |
| **Vertex-transitive (VT) tournament** | Aut(T) acts transitively on vertices | Includes Paley; M = (H/n)I iff self-converse | THM-052 |
| **Frobenius tournament** F_21 | Non-circulant VT at n=21 via Frobenius group | First VT that is NOT self-converse | MISTAKE-013 |

### Paths and Cycles
| Concept | Definition | Key Properties | Where Used |
|---------|-----------|---------------|------------|
| **H(T)** | # directed Hamiltonian paths | Always odd (Rédei); H = I(Ω(T), 2) (OCF) | Central |
| **Directed k-cycle** | k vertices forming directed cycle | Odd cycles only matter for OCF | THM-002 |
| **c_k(T)** | # directed k-cycles in T | c_3 = sum C(s_i, 2) formula; divisibility by C(p,k) | T036, Savchenko |
| **3-cycle matching number** mm(T) | Max # pairwise vertex-disjoint 3-cycles | mm ≤ 2 ⟹ poisoning graph argument | THM-079 Part R |
| **Base path** P_0 | n → n−1 → ⋯ → 1 | Fixed reference path for tiling model | definitions.md |
| **Oriented Hamiltonian path** | Ham path with prescribed arc directions | El Sahili: resilient under arc deletion for n≥8 | INV-099 (new) |

### Conflict Graph and Independence Polynomial
| Concept | Definition | Key Properties | Where Used |
|---------|-----------|---------------|------------|
| **Conflict graph** Ω(T) | Vertices = odd cycles of T; edge iff share vertex | Central to OCF | THM-002 |
| **Independence polynomial** I(G,x) | Σ_k α_k x^k; α_k = # independent sets of size k | I(Ω(T), 2) = H(T) | THM-002 (OCF) |
| **α_k** | # independent sets of size k in Ω(T) | H = 1 + 2α₁ + 4α₂ + 8α₃ + ... | OCF decomposition |
| **μ(C)** | I(Ω(T−v)|_{avoid C\{v}}, 2) | Weights in Claim A RHS | definitions.md |
| **Ω_3(T)** | Subgraph of Ω restricted to 3-cycles only | Quasi-regular; degree ≈ Johnson graph J(n,3) | INV-041 |
| **Hard-core lattice gas** | I(G,x) = partition function at fugacity x | λ=2 outside uniqueness for Ω; special structure | T006, T050 |
| **Typed independence polynomial** | I_typed(Ω; y_3, y_5, ...) | Separates cycle lengths; specializes to I(Ω,x) at y_k=x | THM-068/069/070 |

### Tiling Model
| Concept | Definition | Key Properties | Where Used |
|---------|-----------|---------------|------------|
| **Tile** (a,b) | a ≥ b+2; bit 0 = forward arc, bit 1 = backward | m = C(n−1, 2) tiles total | definitions.md |
| **Pin grid** Grid(n) | (r,c): r≥1, c≥1, r+c≤n−1 | Isomorphic to δ_{n−2} (staircase diagram) | T009 |
| **Strip** Str(k) | {(r,c): r+c=k} | k−1 tiles per strip | definitions.md |
| **Hamming weight** | # backward arcs in tiling | Bell-curve correlation with H | S15 |

---

## II. KEY THEOREMS AND RESULTS

### Proved Theorems
| ID | Name | Statement (brief) | Proof Method |
|----|------|-------------------|-------------|
| **THM-002** | OCF | H(T) = I(Ω(T), 2) | Grinberg-Stanley (arXiv:2307.05569) + specialization |
| **THM-003** | Claim B | I(Ω(T),2) − I(Ω(T−v),2) = 2Σ_{C∋v} μ(C) | A-clique argument (internal) |
| **CONJ-001** | Claim A | H(T) − H(T−v) = 2Σ_{C∋v} μ(C) | OCF + Claim B |
| **THM-016/017** | Even-odd split | Σ_{|S| even} Δ(S,R) = Σ_{|S| odd} Δ(S,R) | Internal induction |
| **THM-020** | Real roots n≤8 | I(Ω(T),x) all real negative roots for n≤8 | Claw-free + Chudnovsky-Seymour |
| **THM-024** | SC involution anti-aut | Every SC tournament has involution anti-aut | Moon + Cauchy's theorem |
| **THM-025** | Real roots FAIL n=9 | Counterexample: score [1,1,3,4,4,4,6,6,7] | Explicit construction |
| **THM-027** | BIBD maximizes H | BIBD arrangement maximizes α₁ (total directed cycles) | Counting + Jensen |
| **THM-029** | H=7 impossible | α₁=3 with i₂=0 forces c₅≥1 ⟹ α₁≥4 | Combinatorial forcing |
| **THM-030** | Transfer matrix symmetry | M[a,b] = M[b,a] for all tournaments | Induction on |W| via Walsh analysis |
| **THM-052** | Scalar M for SC VT | M = (H/n)I for self-converse VT tournaments | Aut+Anti transitivity |
| **THM-079** | H=21 impossible | H(T) ≠ 21 for ALL tournaments on ALL n | Dichotomy + base cases |

### Key Lemmas in H=21 Proof
| Part | Name | Statement | Method |
|------|------|-----------|--------|
| **Part C** | 3-matching bound | 3 disjoint 3-cycles ⟹ α₁+2α₂+4α₃ ≥ 13 | Independence counting |
| **Part G** | Base case | H=21 absent for n≤8 | Exhaustive computation |
| **Part J** | Non-cyclic removal | v not in any 3-cycle ⟹ v not in any cycle | Structural (shortest cycle argument) |
| **Part Q** | Cycle-rich no source/sink | Every vertex in 3-cycle ⟹ no source/sink | Trivial contradiction |
| **Part R** | Dichotomy Theorem | Cycle-rich n≥9: mm≥3 OR safe deletion exists | Poisoning graph DAG |

### Disproved Conjectures
| ID | Name | What Failed | Where |
|----|------|------------|-------|
| **CONJ-002** | Paley Ham paths | H(T_p) ≠ expected at p=11 | H=95095 ≠ 4455 |
| **THM-025** | Real roots all n | Fails at n=9 | Score [1,1,3,4,4,4,6,6,7] |
| **MISTAKE-010** | Hereditary maximizer | Only regular maximizers at odd n | n=5: 24/64 hereditary |
| **MISTAKE-013** | VT ⟹ SC | False at n=21 (Frobenius F₂₁) | McKay database |

---

## III. PROOF TECHNIQUES

### Internal Techniques
| Technique | Description | Used In |
|-----------|------------|---------|
| **A-clique argument** | Through-v cycles form clique in Ω; deletion = clique removal | THM-003 (Claim B) |
| **Poisoning graph** | DAG on outer vertices R; w→v iff all w's 3-cycles contain v | THM-079 Part R |
| **Swap involution** | Swap i,j positions in Ham path; matches adj(i,j) with adj'(j,i) | T031, THM-014 |
| **Bracket structure** | 4-way vertex classification (M+, M−, Z0, Z1) under arc flip | T047 |
| **Fourier decomposition** | OCF decomposes into degree-homogeneous identities | INV-050 |
| **Even-odd split** | Alternating-sum vanishing of subset contributions | THM-016/017, T040 |
| **Transfer matrix** | M[a,b] = Σ_S (−1)^|S| E_a(S)·B_b(R\S) | THM-030 |
| **Walsh spectrum analysis** | Expand M[a,b] in Walsh basis; show only even r-powers | THM-030, THM-080 |
| **Pin grid / tiling encoding** | Tournament ↔ binary string via staircase Young diagram | definitions.md |

### External Techniques Referenced
| Technique | Source | How Used |
|-----------|--------|---------|
| **Chudnovsky-Seymour (2007)** | Real roots for claw-free I.P. | THM-020: proves real roots n≤8 |
| **Heilmann-Lieb (1972)** | Matching polynomial real roots | T054: Ω(T) line graph ⟹ real roots |
| **Lindström-Gessel-Viennot (LGV)** | Lattice path determinants | T046: potential bijection route |
| **RSK correspondence** | Tableaux ↔ permutations | T009: pin grid = Young diagram |
| **Stanley-Stembridge conjecture** | e-positivity of chromatic SF | T056: via Mitrovic-Stojadinovic bridge |
| **Lovász Local Lemma** | Independence polynomial zero-free ↔ LLL | Hard-core lattice gas connection |
| **Bethe ansatz** | Statistical mechanics transfer matrices | T006: potential import |
| **Frankl matching bound** | Max k-sets with no (s+1)-matching | INV-100: bounds on Ω with bounded mm |
| **Lichiardopol's conjecture** | min outdeg ≥ (q−1)k−1 ⟹ k disjoint q-cycles | INV-098: proved for q=3 |

---

## IV. NAMED MATHEMATICIANS AND THEIR CONTRIBUTIONS

| Person(s) | Contribution | Connection to Project |
|-----------|-------------|----------------------|
| **Rédei (1934)** | Every tournament has odd # Hamiltonian paths | The foundational theorem we refine |
| **Berge** | Extended Rédei to digraphs (via permanent) | Grinberg-Stanley generalizes |
| **Grinberg & Stanley (2023)** | Rédei-Berge symmetric function; OCF via specialization | Proves THM-002 (OCF) |
| **Forcade (1973)** | F₂-invariance of H for k-block decompositions | We give first combinatorial proof |
| **Moon** | |Aut(T)| is odd for all tournaments | Used in THM-024 |
| **Chudnovsky & Seymour (2007)** | I(G,x) real roots for claw-free G | Proves THM-020 |
| **Jerrum & Patel (JLMS 2026)** | Zero-free regions for H-free graph classes | Extends real roots to subdivided-claw-free |
| **Lichiardopol** | Conjecture: min outdeg condition for disjoint cycles | Proved for q=3; context for H=21 |
| **Bang-Jensen, Bessy, Thomassé** | Proved Lichiardopol's conjecture for q=3 | INV-098 |
| **Chen & Chang (2024)** | Disjoint cycles in tournaments (JGT) | INV-099 |
| **Frankl** | Erdős matching conjecture for k=3 | INV-100 |
| **El Sahili & Ghazo Hanna (2023)** | Oriented Ham path types in tournaments | T and T^op same distribution |
| **El Sahili & El Zein (2025)** | Ham paths stable under arc deletion for n≥8 | NEW: arXiv:2512.09332 |
| **Mitrovic (2025)** | Noncommuting Rédei-Berge; deletion-contraction | arXiv:2504.20968 |
| **Mitrovic & Stojadinovic (2025)** | Chromatic ↔ Rédei-Berge bridge | arXiv:2506.08841 |
| **Savchenko (2016-2024)** | Exact c_k formulas for regular tournaments | Phase transition at n=39 |
| **Herman (2024)** | Terwilliger algebras classify DRTs to n=23 | 237 non-iso DRTs at n=27 |
| **Feng (2025)** | Dual Burnside process; Q=AB factorization | arXiv:2510.25202 |
| **Schweser, Stiebitz, Toft (2025)** | Rédei's theorem revisited | arXiv:2510.10659 |
| **Irving & Omar (2024)** | Walk generating function; det/per formula | arXiv:2412.10572 |
| **Rajkumar et al.** | Flip classes; R-cones; sign rank | arXiv:2110.05188 |
| **Hetyei (2017)** | Alternation acyclic tournaments; Genocchi numbers | arXiv:1704.07245 |
| **Alon (1990)** | max H(T) = Θ(n!/2^{n−1}) via random regular | Upper bound theory |
| **Szele** | max H(T) ≥ n!/2^{n−1} | Lower bound on max H |
| **Guo, Gutin et al. (2026)** | Forward arc maximization in tournament generalizations | NEW: arXiv:2602.10713 |
| **Bencs & Buys (2025)** | Zero-free regions for hypergraph I.P. | Random Struct. Algorithms |
| **Galvin, McKinley, Perkins, Sarantis, Tetali** | Zeroes of hypergraph independence polynomials | CPC 2024 |
| **Stembridge** | Self-evacuating SYT count | T008 connection |
| **Chapman (2001)** | ASM ↔ monotone triangles | INV-020 |
| **Striker (2011)** | Unifying poset perspective (ASM, PP, Catalan) | INV-020 |
| **Leake & Ryder** | Same-phase stability for claw-free | THM-078 connection |

---

## V. CROSS-FIELD CONNECTIONS

### Established Connections
| Field | Connection | Strength | Reference |
|-------|-----------|----------|-----------|
| **Statistical mechanics** | I(G,x) = hard-core partition function at fugacity x | Strong | T006, T050 |
| **Representation theory** | Pin grid = staircase Young diagram δ_{n−2}; hook lengths all odd | Medium | T009 |
| **Algebraic combinatorics** | U_T = Rédei-Berge symmetric function | Strong | Grinberg-Stanley |
| **Chromatic polynomial theory** | Chromatic SF ≈ Rédei-Berge at poset level | Strong | Mitrovic-Stojadinovic |
| **Hopf algebras** | Deletion-contraction for noncommuting W_X | Medium | Mitrovic 2025 |
| **Number theory** | Paley tournaments via quadratic residues; 1729 appearance | Medium | T019, T025 |
| **2-adic analysis** | H(T) mod 2^k tower from I(Ω,2^k) | Speculative | T007 |
| **Plane partitions / ASMs** | 2^{m²} count; TSSCPP connection | Weak | T008 |
| **Matching theory** | 3-cycle matchings ↔ hypergraph matching | Strong | Frankl, Lichiardopol |
| **Spectral graph theory** | Gauss sums in Paley tournament eigenvalues | Medium | T024 |
| **Design theory (BIBD)** | Regular tournament 3-cycles as block designs | Strong | THM-027 |

### Speculative / Novel Connections (kind-pasteur-S34)

#### 1. GEOMETRIC MODEL: Tournament as Flow on Simplex
**Idea:** A tournament on n vertices defines an orientation of the complete graph K_n. The complete graph K_n is the 1-skeleton of the (n−1)-simplex Δ^{n−1}. An orientation of K_n is equivalent to a discrete vector field on the 1-skeleton. Hamiltonian paths are maximal gradient-like trajectories. The parity constraint (Rédei: H odd) becomes a topological statement about the Euler characteristic of the space of such trajectories.

**Connection to Morse theory:** A discrete Morse function on a simplicial complex assigns values to faces such that the resulting discrete gradient has special properties. A tournament is a discrete Morse function on the 0-skeleton extended to the 1-skeleton. H(T) counts the number of gradient paths of maximum length.

**Potential:** Could import discrete Morse theory (Forman) machinery. The independence polynomial I(Ω,x) might have a topological interpretation as a Poincaré polynomial of some associated complex.

#### 2. TROPICAL GEOMETRY: Permanent as Tropical Determinant
**Idea:** The permanent of a matrix is the tropical determinant (max-plus algebra determinant). For tournaments, H(T) involves counting paths, which is related to permanents of {0,1} matrices. In tropical geometry, the permanent computes the weight of the optimal assignment problem. The tournament arc matrix A defines a tropical hypersurface whose Newton polygon encodes the cycle structure.

**Connection:** Irving-Omar's formula ham(D) = Σ_S det(Ā[S])·per(A[S^c]) mixes determinants and permanents. This is exactly the structure that appears in tropical geometry when studying mixed volumes.

#### 3. HOMOLOGICAL ALGEBRA: Cycle Complex of a Tournament
**Idea:** Define a chain complex C_*(T) where C_k = free abelian group on directed k-cycles, with boundary ∂ = alternating sum of "face" maps (remove one vertex from cycle). The homology H_*(C_*(T)) could encode the independence structure of Ω(T).

**Specific conjecture:** H_0(C_*(T)) = Z (connected), and the Euler characteristic χ(C_*) = Σ (−1)^k rank(C_k) relates to H(T) mod 2.

**Connection to persistent homology:** As we vary a threshold on arc "strength" (e.g., score difference), the cycle complex changes. The persistence diagram could reveal the hierarchical structure of tournament cycles.

#### 4. QUANTUM GROUPS: Tournament as R-matrix
**Idea:** A tournament T defines a solution to the Yang-Baxter equation if we set R(i,j) = q when i→j and R(i,j) = q^{-1} when j→i. The trace of the resulting braid representation counts something related to Hamiltonian paths. The quantum group U_q(sl_2) representation theory at q = root of unity could connect to the 2-adic tower (T007).

**Connection to knot invariants:** Tournaments on n vertices define braids on n strands. The Jones polynomial of the resulting closure might encode H(T). Since the Jones polynomial has connections to statistical mechanics (Potts model), this could unify T006 (hard-core lattice gas) with the tournament structure.

#### 5. ERGODIC THEORY: Tournament Shift Dynamics
**Idea:** The arc-flip operation (reverse one arc) defines a random walk on the space of tournaments. The mixing time of this walk, the spectral gap of the transition matrix, and the stationary distribution all encode information about H(T). The OCF identity H(T) = I(Ω(T), 2) might be a fixed-point equation for some natural dynamics.

**Connection to Markov chain Monte Carlo:** The Jerrum-Patel results on zero-free regions of I(G,x) are motivated by MCMC sampling. Their algorithms use the polynomial method to sample independent sets. For tournament conflict graphs, the specific structure of Ω(T) might allow faster sampling.

#### 6. MODULAR FORMS: Paley Tournament L-functions
**Idea:** The Paley tournament T_p is defined via the Legendre symbol χ_p. The adjacency matrix eigenvalues involve Gauss sums (T024). Define L_T(s) = Σ_{n≥1} a_n n^{−s} where a_n encodes cycle counts. For Paley tournaments, L_{T_p}(s) might factor through Dirichlet L-functions L(s, χ_p).

**Connection:** The factorization H(T_p) = 55 × 1729 at p=11 could reflect the factorization of L(1, χ_11) or related special values. The sequence H(T_p)/|Aut| = 1, 9, 1729 should be checked against modular form coefficient tables.

#### 7. ALGEBRAIC K-THEORY: Tournament Grothendieck Group
**Idea:** Define K_0(Tour_n) as the free abelian group on isomorphism classes of n-tournaments modulo the relation [T] = [T'] + [T''] when T is obtained by "gluing" T' and T'' at a vertex. The OCF H(T) = I(Ω(T), 2) defines a ring homomorphism K_0(Tour_n) → Z.

**Connection:** The vertex deletion formula H(T) − H(T−v) = 2Σμ(C) is a Euler characteristic relation in this Grothendieck group.

#### 8. INDEPENDENCE COMPLEX TOPOLOGY: Ω(T) as Simplicial Complex
**Idea:** The independence complex Ind(Ω(T)) is the simplicial complex whose faces are independent sets of Ω(T). For claw-free graphs (n≤8), Ind(G) is homotopy equivalent to a wedge of spheres (Engström's work). The reduced Euler characteristic of Ind(Ω(T)) equals (−1)^d · I(Ω(T), −1) where d is the dimension. At x=2, I(Ω,2) = H(T) counts 2-colored independent sets, which has a topological interpretation as the number of simplices weighted by a 2-coloring.

**Connection to persistent homology:** Define a filtration on Ind(Ω(T)) by the size of cycles (3-cycles first, then 5-cycles, etc.). The resulting persistence diagram tracks how the topology changes as longer cycles are added. This could reveal WHY certain H values are impossible (gaps correspond to topological obstructions).

**Concrete test:** Compute the homology of Ind(Ω(T)) at n=5,6,7. Check if the Betti numbers relate to α_k coefficients. For the H=21 gap: does the corresponding independence complex have a topological obstruction?

#### 9. SIGNED HP PERMANENT AS PFAFFIAN VARIANT
**Idea (from opus S35c):** S(T) = Σ_P Π B[P_i][P_{i+1}] is a "path permanent" — distinct from the standard permanent (which sums over cycle covers). The standard perm(B) = 0 at odd n for skew-symmetric B, but S(T) ≠ 0. S(T) captures Hamiltonian PATH information invisible to both det(B) and perm(B).

**Key discovery:** S(T) mod 2^{n−1} depends ONLY on n (Universal Congruence Theorem, THM-H). This means the "residual" c_0 = S(T)/2^{n−1} has universal fractional part determined by n alone. At n=5: c_0 ∈ Z (integer). At n=7: c_0 has fractional part 3/4 always. This 2-adic structure connects to the 2-adic tower (T007).

**Connection to knot theory:** The signed adjacency matrix B is skew-symmetric, like the Seifert matrix of a knot. The Pfaffian of B gives the Alexander polynomial value. S(T) is a different "Pfaffian-like" invariant that uses paths instead of matchings. Could S(T) be a tournament analogue of a knot invariant?

#### 10. RANDOM MATRIX THEORY: Tournament Eigenvalue Distribution
**Idea:** The skew-adjacency matrix S_T (with S_{ij} = 1 if i→j, S_{ij} = −1 if j→i, S_{ii} = 0) has purely imaginary eigenvalues. For random tournaments, the empirical spectral distribution converges to the semicircle law scaled by √n. The Pfaffian Pf(S_T) relates to the number of perfect matchings of the underlying graph, which connects to the cycle structure.

**Connection to H(T):** H(T) = permanent of a certain {0,1} matrix. The permanent-determinant gap (permanent is #P-hard, determinant is polynomial) is the central problem in algebraic complexity. For tournament matrices, the permanent might have special structure due to the constraint T(i,j) + T(j,i) = 1.

#### 11. INFORMATION THEORY: Tournament Entropy
**Idea:** Define the entropy of a tournament as H_ent(T) = −Σ_P (1/H(T)) log(1/H(T)) = log(H(T)) (uniform distribution over Hamiltonian paths). The OCF H(T) = I(Ω(T), 2) then gives log H(T) = log I(Ω(T), 2), connecting tournament entropy to the free energy of the hard-core model. The 2-adic valuation v_2(H(T)) measures how "2-divisible" the path count is — a kind of 2-adic entropy.

**Connection to coding theory:** Tournament tilings (binary strings of length m = C(n−1,2)) form a code where each codeword represents a tournament. The Hamming weight structure of this code correlates with H(T) (bell curve, T053). Error-correcting properties of the "tournament code" could import coding theory machinery.

#### 12. CLUSTER ALGEBRAS: Tournament Mutations
**Idea:** In cluster algebra theory, mutations on quivers (directed graphs) produce new quivers. A tournament IS a quiver (complete, no 2-cycles). Cluster mutations preserve certain algebraic structure. The OCF identity H(T) = I(Ω(T), 2) could be a cluster variable identity, with arc-flip = cluster mutation.

**Connection:** Fomin-Zelevinsky cluster algebras on complete quivers have been studied. The exchange graph (quivers related by mutations) for tournaments could be the arc-flip graph. The Laurent phenomenon (cluster variables are Laurent polynomials) might constrain H(T) values, potentially explaining gaps.

**Concrete test:** Check if the cluster variable associated to a tournament quiver equals H(T) or I(Ω(T), x) for some specialization.

#### 13. MATROID THEORY: Tournament Matroid
**Idea:** Define a matroid on the arcs of a tournament where independent sets are "acyclic subsets" (subsets of arcs that contain no directed cycle). The rank function r(S) = |S| − (# cycles in S) could encode the cycle structure. The Tutte polynomial T_M(x,y) of this matroid at specific evaluations could give H(T) or cycle counts.

**Connection to Rédei-Berge:** The Rédei-Berge symmetric function U_T is a "symmetric function analogue" of the Tutte polynomial. The chromatic-Rédei-Berge bridge (Mitrovic-Stojadinovic) makes this explicit: the chromatic symmetric function IS the Tutte polynomial in disguise.

#### 14. OPERAD THEORY: Composition of Tournaments
**Idea:** Tournaments form an operad under substitution: given T on n vertices and T_1,...,T_n, the substitution T(T_1,...,T_n) replaces each vertex i by T_i. This operad structure governs how H(T) decomposes. The vertex-deletion formula H(T) − H(T−v) = 2Σμ(C) is a derivation in this operad.

**Connection:** The operad of associahedra (Stasheff) governs compositions of operations. Tournament substitution is the "directed" analogue. The free operad generated by tournaments might have a presentation that encodes the OCF identity.

#### 15. MIRROR SYMMETRY: Tournament Duality
**Idea:** H(T) = H(T^op) is a "mirror symmetry" for tournaments — the path count is invariant under arc reversal. In mirror symmetry, the Hodge numbers h^{p,q} of a Calabi-Yau manifold X equal h^{n−p,q} of the mirror X^∨. For tournaments, the "Hodge numbers" could be the coefficients of the W-polynomial, and the mirror is T^op.

**Connection to THM-030:** The transfer matrix symmetry M[a,b] = M[b,a] is analogous to the symmetry of the Hodge diamond. The W-polynomial W(r) having only even powers of r is analogous to the vanishing of odd Hodge numbers.

#### 16. PERSISTENT HOMOLOGY OF TOURNAMENT FILTRATIONS
**Idea:** Order tournaments by "complexity": T_0 (transitive) → T_1 (one arc flip) → ⋯ → T_m (all arcs flipped). Each step changes Ω(T), and the persistence barcode of this filtration tracks which cycles are born and die. The OCF identity H(T) = I(Ω(T), 2) says the partition function is preserved through the filtration (since H is determined by Ω).

**Novel observation:** The arc-flip graph on tournaments is the Boolean lattice {0,1}^m where m = C(n−1,2). The OCF is a function on this lattice. Its "topological complexity" (number of critical points in the Morse sense) could measure the difficulty of proving OCF.

#### 17. REPRESENTATION STABILITY: Tournament Invariants as n→∞
**Idea:** Church-Ellenberg-Farb representation stability: for families of S_n-representations {V_n}, the decomposition into irreducibles stabilizes. The space of tournament invariants (functions T ↦ f(T) invariant under relabeling) is an S_n-representation. Does it stabilize?

**Connection:** The Rédei-Berge symmetric function U_T lives in Λ_n (symmetric functions of degree n). As n grows, the relevant pieces of Λ_n could stabilize, giving eventually-polynomial formulas for H(T) statistics.

#### 18. PHYSICS: TOURNAMENT AS CAUSAL STRUCTURE
**Idea:** A tournament defines a total preorder on pairs — exactly one of i→j or j→i holds, like a causal relation. In physics, causal sets are partially ordered sets modeling spacetime. A tournament is a "maximally connected" causal set. Hamiltonian paths are "world lines" visiting every event exactly once. The parity of H(T) (always odd, Rédei) is a "topological charge" of the causal structure.

**Connection to loop quantum gravity:** The spin foam models of loop quantum gravity use labeled graphs and their amplitudes. A tournament with weighted arcs defines a spin foam amplitude. The independence polynomial I(Ω,x) at x=2 could be a partition function in this framework.

---

## VI. POLYNOMIAL OBJECTS

| Polynomial | Definition | Key Properties | Where |
|-----------|-----------|---------------|-------|
| **I(Ω(T), x)** | Independence polynomial of conflict graph | I(Ω,2) = H(T); real roots n≤8 | THM-002, THM-020 |
| **U_T** | Rédei-Berge symmetric function | U_T ∈ Λ_n; p-positive for tournaments | Grinberg-Stanley |
| **G_T(t, x)** | Two-variable generating function | t^m P(u,x), u = t+1/t; P(2,x) = n! | THM-064, S28 |
| **P(u, x)** | Reduced polynomial from G_T | p_m(x) = I(Ω(T), x) as leading coeff | S28 |
| **Q_T(w)** | Size-weighted I.P.: u_T(√w)/√w | All real roots n≤8; fails n≥9 with I(Ω) | THM-078 |
| **W(z)** | Walsh generating function | W(i/2) = (−2)^m p_0(2) | THM-064 |
| **M[a,b]** (transfer matrix entry) | Σ_S (−1)^|S| E_a(S)·B_b(R\S) | Symmetric; Walsh spectrum has 2^s factors | THM-030, THM-080 |
| **E_T^perm(t)** | All-permutation forward-edge polynomial | G_T(t; 2,2,...) specialization | S28 |
| **S(T)** (signed HP permanent) | Σ_P Π B[P_i][P_{i+1}], B=2A−1 | S=0 at even n; S mod 2^{n−1} universal | THM-A through THM-H |
| **W(r)** (W-polynomial) | tr(M(r)); only even r-powers | c_{n−1}=n!, c_{n−3} depends on t₃ | signed-hp-permanent-skeleton |
| **D_k** | Σ_P C(forward(P), k) | D_k mod 2^{n−1−k} is universal | THM-H |

---

## VII. GRAPH-THEORETIC HIERARCHY OF Ω(T)

| Property | Holds for n ≤ | Fails at n = | Reference |
|----------|-------------|-------------|-----------|
| Complete | 5 | 6 | T028 |
| Interval | 5 | 6 (13.9%) | T055 |
| Chordal | 5 | 6 | T055 |
| Line graph | 5 | 6 (K₅−e found, 45%) | INV-032 |
| Comparability | 6 | 7 (1%) | T049 |
| Quasi-line | 7 | 8 (49%) | INV-032 |
| Perfect | 7 | 8 (53.8%) | INV-032 |
| Claw-free | 8 | 9 (86%) | INV-032 |
| S_{1,1,1}-free | 11 | 12 | INV-032 |
| S_{2,1,1}-free | 9 | 10 (92%) | INV-032 |
| Real roots of I(Ω,x) | 8 (proved) | 9 (extremely rare failure) | THM-020, THM-025 |

---

## VIII. PERMANENT GAPS IN H-SPECTRUM

| Value | Status | Method |
|-------|--------|--------|
| H = 7 | IMPOSSIBLE for all n | THM-029: α₁=3 with i₂=0 forces extra cycles |
| H = 21 | IMPOSSIBLE for all n | THM-079: Dichotomy + base cases |
| H = 23 | OPEN — candidate gap | opus investigating (h23_achievability.c) |
| H even | IMPOSSIBLE at odd n | OCF: I(Ω,2) always odd when n odd |
| H = 2 (mod 4) | OPEN at even n | Need systematic enumeration |

---

## IX. OEIS CONNECTIONS

| Sequence | What | Values | Status |
|----------|------|--------|--------|
| **A038375** | Max H(T) for n-vertex tournaments | 1,1,3,5,45,189,661,... | Paley achieves max at prime n |
| **A000213** | Tribonacci numbers | Related to transitive tournaments | Connection unclear |
| **H(T_p)/|Aut|** | 1, 9, 1729 | NOT in OEIS | Candidate new sequence |
| **H(T_p)** | 3, 189, 95095 | NOT in OEIS | Candidate new sequence |
| **Tangent numbers** | P_n(0,0) = 2^{(n−1)/2} T_n | Connected via EGF | INV-093 |

---

## X. MISTAKES LOG (DO NOT REPEAT)

| ID | What Went Wrong | Correct Statement |
|----|----------------|-------------------|
| MISTAKE-001 | ind_poly_at_2_restricted() bug | Never use old scripts |
| MISTAKE-003 | Per-path identity extends to all n | FAILS at n≥6 |
| MISTAKE-008 | Even-odd split ≡ OCF | NOT equivalent |
| MISTAKE-010 | Hereditary maximizer for all maximizers | Only REGULAR at odd n |
| MISTAKE-011 | Paley at p ≡ 1 mod 4 | Only p ≡ 3 mod 4 |
| MISTAKE-013 | VT ⟹ SC | FALSE at n=21 (Frobenius) |
| MISTAKE-014 | Scalar M for all VT | Only SC VT at odd n |

---

## XI. SOFTWARE AND DATA

### Key Scripts
| Script | Purpose | Status |
|--------|---------|--------|
| `symbolic_proof.py` | Exhaustive OCF verification | n≤7 complete |
| `fourier_homogeneity.py` | Fourier decomposition of OCF | Proved for n=5,7 |
| `symmetry_check.py` | Transfer matrix M[a,b] = M[b,a] | Verified n=4,...,8 |
| `h21_dichotomy_proof.py` | Dichotomy verification at n=9 | 106,424 tests, 0 failures |
| `h21_poisoning_graph.py` | Poisoning graph structure analysis | 51,280 mm=2 cases |
| `interlacing_verify.py` | Clique-deletion interlacing | 0 failures n=5,...,8 |
| `paley_deletion_test.py` | Paley hereditary maximizer | Verified p=3,7,11 |

### External Databases
| Database | What | Used For |
|----------|------|---------|
| McKay tournament database | All tournaments to n=10; DRTs to n=23 | Verification |
| OEIS A038375 | Max Hamiltonian paths | Paley maximizer conjecture |
| arXiv | Papers | Literature connections |
