# Investigation Backlog

**Purpose:** Systematic catalog of every lead, reference, connection, and unexplored direction extracted from the repo. Claude agents should consult this before choosing what to work on, and add new leads as they emerge. Prioritized by potential impact on proving OCF (Claim A).

**Last full repo scour:** kind-pasteur-2026-03-06-S18e
**Last web research:** opus-2026-03-05-S9 (Paley maximizer, n=8 anomaly)

---

## Priority A: Key structural questions (OCF PROVED by Grinberg-Stanley)

### INV-032: Omega(T) structural properties — PARTIALLY DISPROVED
**Source:** Web research opus-S5, opus-S7 (disproof), opus-S9 (line graph disproof), opus-S10 (structure analysis)
**Status:** DISPROVED: Omega(T) is NOT always claw-free (fails n=9, 90%) or perfect (fails n=8, 53.8%). NOT a line graph (K_5-e found at n=6, 45%). S_{1,1,1}-free through n=11, fails n=12.
**What remained true:** All-real-roots of I(Omega(T), x) appears to hold even for imperfect/non-claw-free Omega (tested n<=10, 0 failures). This is a deep structural conjecture NOT explained by any forbidden subgraph property.
**Note:** OCF is now proved by Grinberg-Stanley, so this is no longer a proof strategy — it's a structural question. Real-rootedness explanation must be algebraic (Irving-Omar/Grinberg-Stanley symmetric function framework).
**Next step:** (1) Test all-real-roots at n=9 where Omega is genuinely non-claw-free. (2) If confirmed, formulate as conjecture and seek algebraic proof.

### INV-038: Clique-deletion interlacing for Omega(T)
**Source:** opus-2026-03-06-S17, T100, interlacing-clique-deletion.md
**Status:** STRUCTURAL INSIGHT. Proof sketch for n<=8.
**What:** Through-v cycles always form a CLIQUE in Omega(T) (proved: sharing vertex v). Deleting vertex v = deleting this clique from Omega(T). The sequential deletion recurrence I(G,x) = I(G-u,x) + x*I(G\N[u],x) can be applied step-by-step. At n=5, 100% of remaining cycles are adjacent to some through-v cycle (Omega is very dense).
**Key insight:** For n<=8 (claw-free), Chudnovsky-Seymour guarantees each step preserves real roots. For n>=9, the clique structure + high density may still force interlacing.
**Verification:** 0 failures: n=5 (5120 exhaustive), n=6 (196608 exhaustive), n=7-8 (random).
**Impact:** If provable for ALL n, gives inductive proof of real-rootedness of I(Omega(T),x).
**Next step:** (1) Verify n=7 exhaustive. (2) Check whether the "remainder" G\N[u] always has real-rooted I.P. (3) Look for common interlacing framework.
**Scripts:** `04-computation/interlacing_verify.py`, `04-computation/interlacing_structure.py`
**Writeup:** `03-artifacts/drafts/interlacing-clique-deletion.md`

### INV-039: Blueself odd-n obstruction — PROVED for ALL odd n
**Source:** opus-2026-03-06-S17, THM-022 Theorem 5 (upgraded)
**Status:** PROVED. Pure algebraic proof, no exhaustive search needed.
**What:** No blueself tilings exist at any odd n. Grid-symmetry forces k_0+k_{n-1}=n-2 (endpoint constraint). Flip changes endpoint multisets: {1+k_0, n-2-k_0} -> {n-1-k_0, k_0}. For these to be equal as multisets, need k_0=(n-2)/2 (non-integer at odd n) or 1=0 (impossible). Therefore sorted scores always differ, so flip(T) is never isomorphic to T.
**Script:** `04-computation/blueself_odd_n_proof.py`
**Impact:** Upgrades THM-022 Theorem 5 from "proved n<=7" to "proved all n". Completes the odd half of the blueself existence dichotomy.

### INV-040: Blueself vs SC maximizer — DISPROVED
**Source:** opus-2026-03-06-S17, T099
**Status:** DISPROVED at n=6.
**What:** At n=6, blueself class with H=41 is NOT the SC maximizer in score class (3,3,3,2,2,2) (SC max is H=45, also blueself). Blueself classes are always SC and have regular scores, but not always max-H. The blueself with higher disjoint pair count (alpha_2=4) beats the one with more total cycles (alpha_1=16).
**Script:** `04-computation/blueself_sc_maximizer_connection.py`

### INV-041: Quasi-regularity of Omega(T) — EXPLAINED
**Source:** opus-2026-03-06-S17 (T101), opus-2026-03-06-S18 (T103, proof)
**Status:** EXPLAINED. Theoretical argument + verified n=5-20.
**What:** Omega_3(T) is quasi-regular because adjacency depends on vertex-set intersection (sharing ≥1 vertex), not arc orientations. This makes Omega_3 an induced subgraph of J(n,3) (Johnson graph), inheriting its regularity. All 3-element subsets have identical intersection statistics, so degree of each 3-cycle concentrates around E[deg] = (C(n,3)−C(n−3,3))/4−1. The coefficient of variation CV = O(1/√m) → 0 as n→∞, giving λ_max/avg_deg ≈ 1+CV² → 1. Verified: CV drops from 0.05 (n=6) to 0.03 (n=20). This does not directly explain real-rootedness but constrains spectral structure.
**Scripts:** `04-computation/omega_spectral_fast.py`, `04-computation/omega_quasireg_proof.py`

### INV-042: Paley deletion maximizer — VERIFIED p=3,7,11
**Source:** kind-pasteur-2026-03-06-S18e (T097), opus-2026-03-06-S18
**Status:** VERIFIED at p=3,7,11. Conjecture for all Paley primes.
**What:** H(T_p − v) = a(p−1) (OEIS A038375 max H at n=p−1). Verified: T_3−v: H=1=a(2), T_7−v: H=45=a(6), T_11−v: H=15745=a(10). By vertex-transitivity all deletions equivalent. Combined with T053 (T_p achieves a(p)), the maximizer chain is "hereditary" via vertex deletion. Claim A decomposition for T_7: diff=144=2×72, sum_mu=6×3+30+24=72 (all 3-cycle complements have a 3-cycle in Paley).
**Next step:** Verify at p=19 (need H(T_19−v) = a(18)). Investigate: does the n=p−1 maximizer always come from Paley deletion, or can it be achieved by non-Paley tournaments too?
**Scripts:** `04-computation/paley_deletion_test.py`

### INV-043: Anti-aut involution existence — PROVED (THM-024)
**Source:** opus-2026-03-06-S18 (T102, THM-024), correcting kind-pasteur S18e (T095)
**Status:** PROVED. Clean group theory argument.
**What:** Every SC tournament has ≥1 involution anti-automorphism. Proof: (1) Moon's theorem: |Aut(T)| is odd. (2) H = ⟨Aut(T), σ₀⟩ has order 2|Aut(T)| (even). (3) By Cauchy, H has order-2 element. (4) Can't be in Aut(T) (odd order group). (5) Must be in σ₀·Aut(T) = set of anti-auts. NOT all anti-auts are involutions (counterexamples at n=6 with |Aut|>1), but at least one always is.
**Scripts:** `04-computation/anti_aut_involution_test.py`, THM-024

### INV-033: Redei-Berge Hopf algebra formalization of OCF
**Source:** Web research opus-S5, arXiv:2402.07606 (Grinberg)
**Status:** CONNECTION IDENTIFIED. NOT formalized.
**What:** The Redei-Berge symmetric function U_X for digraphs has comultiplication Delta([X]) = sum_S [X|S] tensor [X|V\S] — this IS our subset convolution. The character zeta counts Hamiltonian paths. The antipode S(U_X) = (-1)^|V| U(X-bar) encodes Berge's theorem. OCF could be a Hopf algebra identity relating zeta to the independence polynomial of Omega(T).
**Next step:** Read arXiv:2402.07606 in full. Express OCF in Hopf algebra language. Check if I(Omega, 2) has a natural coalgebra interpretation.

### INV-034: Björklund cycle cover reduction adapted for OCF
**Source:** Web research opus-S5, arXiv:1008.0541, arXiv:1301.7250
**Status:** CONNECTION IDENTIFIED. NOT attempted.
**What:** Björklund reduces Hamiltonian cycle counting to cycle cover counting via inclusion-exclusion and determinants. Could a directed version for Hamiltonian PATHS in tournaments reduce specifically to ODD cycle covers, yielding OCF? The characteristic-2 aspects are particularly relevant since Redei is a mod-2 statement.
**Next step:** Study whether Björklund's labeled cycle cover approach specializes to odd cycle covers for tournaments.

### INV-035: Tribonacci structure — OCF for T_full family via interval graphs
**Source:** opus-2026-03-05-S6 (Tribonacci web research), kind-pasteur-S11 (Tribonacci discovery)
**Status:** VERIFIED n=3,...,8. Both sides match Tribonacci(n) = A000213 independently.
**What:** T_full_n (full tiling tournament) has H(T_full_n) = Tribonacci(n) (proved via run decompositions). INDEPENDENTLY, Omega(T_full_n) is an INTERVAL GRAPH on odd-length consecutive intervals [k, k+2j], and I(Omega, 2) satisfies the same Tribonacci recurrence via a weighted interval packing DP that telescopes: f(n) = f(n-1) + 2f(n-3) + 2f(n-5) + ... = f(n-1) + f(n-2) + f(n-3).
**Key structural insight:** All directed odd cycles of T_full_n are consecutive intervals. The clique-cutset decomposition of this interval graph mirrors the DP structure computing H(T_full_n). Both sides produce Tribonacci by the same algebraic mechanism (telescoping) through different combinatorial objects.
**Why this matters:** Shows OCF's "both sides match" emerges from parallel decomposition structures. If this parallelism generalizes (clique-cutset of Omega mirrors Ham path DP), it could prove OCF.
**Extended results (opus-S13):**
- **Transitive+flip(i,j):** H = 1 + 2^(j-i-1). All odd cycles form a clique in Omega, so I(Omega,2) = 1 + 2·(#cycles) = 1 + 2^(j-i-1). Clean OCF-based proof.
- **Cone theorem:** H(source_cone(T')) = H(sink_cone(T')) = H(T') for ALL T'. Proved: source must be first in every Ham path. Verified exhaustively through n'=6.
- **Partial cones palindromic:** H(k) = H(n'-k) where k = out-degree of cone vertex. From self-converse symmetry.
- **Circulant S={1}:** H(T_{n,{1}}) = n, order-2 recurrence. No circulant with |S|>=2 has low-order recurrence.
- **Best circulant at n=9 gives H=3267 < 3357 = max.** H-maximizer at non-prime n is NOT circulant.
**Next step:** (1) Find direct bijection between run decompositions and weighted interval packings. (2) Check if the transfer matrix for T_full has Tribonacci characteristic polynomial factor.

### INV-001: Prove transfer matrix symmetry for all n
**Source:** T045 (tangents), symmetry_check.py, paper-connections.md
**Status:** Verified n=4,...,8 (7500+ tests). NOT yet proved.
**What:** M[a,b] = sum_S (-1)^|S| E_a(S)*B_b(M\S) is always symmetric. This is STRONGER than the even-odd split (MISTAKE-008 shows even-odd split != OCF). Proving this symmetry for all n would prove OCF.
**Key insight:** Connects to Feng's Dual Burnside (Q=AB symmetric under detailed balance). The "forward leg" (paths ending at vertex) and "backward leg" (paths starting at vertex) appear to satisfy a hidden reversibility.
**Next step:** Attempt algebraic proof at n=4,5 using the bracket structure (T047). Try to express M in terms of the bracket table {M-, M+, Z1, Z0} and show symmetry follows from the block structure.

### INV-002: Subset convolution identity — the core algebraic challenge
**Source:** proof-landscape-for-general-ocf.md (Approach B), T047
**Status:** Correct framework identified. No simplification found.
**What:** sum_S [f_i(S)*g_j(R) - f_j(S)*g_i(R)] = sum_{k>=1} 2^k Delta(alpha_k). Both sides are multilinear polynomials. The bracket B(u,w) has a 4-way type structure where Z0 rows and Z1 columns vanish.
**Next step:** Use the bracket table to decompose the convolution into 6 nonzero bracket types. Try to show the resulting expression telescopes via induction on the number of M-/M+ vertices.

### INV-003: Sign-reversing involution on the subset convolution
**Source:** proof-landscape-for-general-ocf.md (bottom), signed-adjacency-identity.md
**Status:** Idea only. Not attempted.
**What:** Find an involution on the terms of sum_S [f_i(S)*g_j(R) - f_j(S)*g_i(R)] that cancels everything except the cycle terms. This is the "bijective" approach to the algebraic identity.
**Connection:** The sigma-invariance (s -> -s, B is even in s-variables) reduces OCF to proving all s-degree-1 terms vanish: C_w + D_w = 0 for each w. This is a per-vertex condition that might be provable.

### INV-004: Flip-class proof strategy (prove for R-cones, extend via cut-flip)
**Source:** T046, paper-connections.md (CONNECTION 1, "Flip Class + OCF")
**Status:** Strategy identified. Not attempted.
**What:** Rajkumar et al. show every tournament is in the flip class of an R-cone. For R-cones (vertex beating/losing to everyone), Ham paths all start or end at universal vertex, simplifying both H(T) and I(Omega,2). Prove OCF for R-cones, then show cut-flip phi_S preserves E(T) = H(T) - I(Omega(T),2) = 0.
**Key gap:** Need to track how both H and I change under phi_S (reversing all arcs across a cut). This is a MULTI-arc flip, not a single arc flip.

### INV-005: Induction on mu(T) (flip-feedback dimension)
**Source:** paper-connections.md (CONNECTION 4)
**Status:** Idea only. Not formalized.
**What:** Rajkumar et al.'s mu(T) measures minimum flip-feedback node set. Dim bound: 2(mu(T)+1). For mu=0, T is transitive (OCF trivial). Could induct on mu(T), with each step being a cut-flip. Need: does cut-flip increase or preserve mu?

### INV-006: n=8 exhaustive proof completion
**Source:** OPEN-Q-009, ocf_n8_full.c
**Status:** PROVED by opus-S4 (2^27 configs, 57min, all passing). Independent verification by opus-S4b C implementation (3M+ configs, 0 fails through partial run).
**Next:** Close this out. Focus on n=9 strategy or general proof.

---

## Priority B: Important structural understanding

### INV-007: Odd-cycle bijection (Open Problem 3 in paper)
**Source:** oq:bijection in tex, bijection_search.py, T046
**Status:** Searched computationally at n=3,4. No natural bijection found.
**What:** Construct Phi: Ham(T) -> {(C, f) : C vertex-disjoint odd-cycle collection, f: C -> {0,1}}. The "2-colored cycles" interpretation means each independent set of k cycles contributes 2^k paths. At n=3 with one 3-cycle: 3 paths = 1 + 2. Need to identify which paths correspond to which colored cycle sets.
**Key obstacle:** The correspondence is NOT local/contiguous-block based (T035, confirmed dead end). Must be global.

### INV-008: Striker-Chapman S3-equivariance (Open Problem 5 in paper)
**Source:** oq:striker in tex, \cite{striker2011,chapman2001}
**Status:** NOT INVESTIGATED AT ALL.
**What:** Is the Striker-Chapman bijection between ASMs and tournaments S3-equivariant under the barycentric identification? Striker (2011) gives a unifying poset perspective connecting ASMs, plane partitions, Catalan objects, tournaments, and tableaux. Chapman (2001) connects ASMs to tournaments directly.
**Why it matters:** If S3-equivariant, then the S3 orbit counts in our Section 6 could be computed via ASM symmetries, potentially giving new structural constraints on H(T).
**Next step:** Read Striker (2011) and Chapman (2001) carefully. Test equivariance computationally at small n.

### INV-009: Self-evacuating SYT bijection (Open Problem 6 in paper)
**Source:** oq:se_bijection in tex, Section 8 (tetrahedral geometry)
**Status:** Count verified (2^{m^2} for n=5,7). Bijection NOT constructed.
**What:** Natural bijection between 2^{m^2} self-evacuating SYT of delta_{n-2} and 2^{m^2} sigma-fixed tilings of Grid(2m+1). Both count the same thing but via very different combinatorial objects.
**Connection:** TSSCPPs (Totally Symmetric Self-Complementary Plane Partitions) and the ASM conjecture. The TSSCPP count for order m is known to equal the ASM count.

### INV-010: Mixed graphs extension (Open Problem 4 in paper)
**Source:** oq:mixed in tex, \cite{schweser2025}
**Status:** NOT INVESTIGATED.
**What:** Extend the Q-Lemma to complete mixed graphs, recovering Schweser-Stiebitz-Toft (2025) strengthening of Redei.
**Next step:** Read arXiv:2510.10659 carefully. Understand what "complete mixed graph" means and how the Q-Lemma generalizes.

### INV-011: Mod-4 score-sequence criterion (Open Problem 9 in paper)
**Source:** oq:mod4_struct in tex
**Status:** Conjectural. NOT tested computationally beyond small cases.
**What:** Does alpha_1(T) (equiv H(T) mod 4) admit a score-sequence characterization? Moon's formula determines |C_3| from score sequence. Question: is alpha_1 = |C_3| (mod 2) in general?
**Next step:** Test computationally at n=5,6,7. Extract alpha_1 and |C_3| for all tournaments.

### INV-012: BlackSelf(8) exceptional class (Open Problem 7 in paper) — PARTIALLY RESOLVED
**Source:** oq:n8 in tex
**Status:** DEEP INVESTIGATION by opus-2026-03-05-S8. See `03-artifacts/drafts/n8-anomaly-deep-dive.md`.
**What:** Unique isomorphism class at n=8 that is self-converse, has |Aut|>1, |Fix(beta)| odd, but H(T)/|Fix(beta)| is even. Is it related to a Hadamard matrix of order 8 or a skew conference matrix?
**Key findings:**
- Exhaustive search over ALL 65536 SC tournaments (alpha = reversal): ZERO have Fix(beta) | H with H/Fix even.
- The definition likely means (H-Fix)/2 is even (number of beta-orbit pairs is even).
- Under this interpretation: T_657 (H=657, Fix=33, |Aut|=3) is the best candidate: (657-33)/2 = 312 (even).
- T_657 CONTAINS P(7) (Paley tournament on 7 vertices) as vertex-deletion. P(7) ↔ unique skew Hadamard matrix of order 8. THIS IS THE HADAMARD CONNECTION.
- T_657 has perfectly uniform D_v = 54 (mu-weighted 3-cycle count) for all 8 vertices.
- Full survey: 10 distinct (H, Fix, Aut) combinations among 2560 SC+Aut>1 tournaments.
**Next step:** Confirm T_657 is isomorphic to a known Paley extension. Resolve definition ambiguity with paper author.

### INV-013: Realizable odd-cycle conflict graphs (Open Problem 8 in paper) — INVESTIGATED
**Source:** oq:realizable in tex
**Status:** INVESTIGATED (opus-S13 background agent). Key structural findings.
**What:** Which graphs G arise as Omega_3(T) for some tournament T?
**Results:**
- Realizable isomorphism classes: 2 (n=3), 3 (n=4), 6 (n=5), 18 (n=6), ~97+ (n=7 sampled)
- n≤5: Omega_3 is ALWAYS a complete graph (two 3-cycles on ≤5 vertices must share a vertex)
- n=6: Always "complete minus matching" — complement is disjoint union of edges
- n=7: First non-perfect graphs (11/97 classes have chi=omega+1)
- alpha(Omega_3) ≤ floor(n/3) — proved by vertex counting (3k vertices for k disjoint 3-cycles)
- 100% real-rooted across all realizable classes (exhaustive n≤6, sampled n=7)
**Key insight:** The low independence number (alpha ≤ 2 for n≤8) means I(Omega_3, x) has degree ≤ 2, so real-rootedness is "easy" for Omega_3. The full Omega (including 5,7-cycles) has alpha ≤ floor(n/3) similarly, keeping degree low.
**Next step:** Characterize which "complete minus matching" graphs at n=6 are NOT realizable. Extend to n=8 where alpha=2 still but degree may increase.

### INV-014: 2-adic tower / higher Redei theorems — PARTIALLY RESOLVED
**Source:** OPEN-Q-008, T007, tex Section 5.5
**Status:** COMPUTED (opus-S13). v_2(H(T)) = 0 ALWAYS (= Redei's theorem).
**What:** I(Omega(T), x) at x=4,8,... gives mod-4, mod-8 invariants of H(T). v_2(H(T)) = 0 universally.
**Results:** H mod 4 ≡ 1+2*alpha_1 (mod 4) via OCF. At n=3,4 this equals 1+2*c3 (mod 4) exactly. At n≥5 the c3 formula breaks (5-cycles contribute to alpha_1). H mod 2^k approaches uniform on odd residues as n grows.
**Impact:** OPEN-Q-008 partially resolved. No deeper 2-adic structure at level of H(T). The mod-4 structure is fully explained by alpha_1 parity via OCF.

---

## Priority C: References to investigate

### INV-015: Rajkumar et al. (arXiv:2110.05188) — tournament representations
**Source:** paper-connections.md, T046, paper-deep-connections.md Section 2
**Status:** INVESTIGATED (opus-S4b, opus-S5). Key theorems extracted. Connections documented.
**What:** Flip classes, locally transitive = rank 2, R-cones, sign rank. Key results: every T in flip class of R-cone (Prop 1, distance exactly 1 cut-flip); mu(T) dimension bound <= 2(mu(T)+1) (Thm 11); sign-rank bound (Thm 12).
**Key finding:** Proposition 1 is constructive: T' = phi_{i ∪ T_i^-}(T) is R-coned by i. Cut-flip distance to R-cone is exactly 1. This directly enables INV-004 strategy.
**Assessment:** mu(T) induction (INV-005) less promising than FAS induction — mu may change unpredictably under cut-flips.
**Action needed:** ADD to bibliography. Concretely develop INV-004 (R-cone + cut-flip proof).
**Tested:** Locally transitive tournaments DO have 5/7-cycles (T046). OCF passes 100% for LT, R-cones, automorphism-symmetric tournaments.

### INV-016: Feng (arXiv:2510.25202) — dual Burnside process
**Source:** paper-connections.md, paper-deep-connections.md Section 3, tex line 1600/2151
**Status:** INVESTIGATED (opus-S4b, opus-S5). Key theorems extracted. Deep connection found.
**What:** Q=AB factorization, primal-dual spectral correspondence, lumping theory.
**Key findings:** (1) Q=AB is REVERSIBLE with pi(g)=|X_g|/(|G|*z) — detailed balance gives symmetry (Thm 3.3). (2) Block-flip M=[[0,A],[B,0]], M^2=[[Q,0],[0,K]] — bipartite structure with period 2. (3) Eigenvector intertwining (Thm 3.10): A maps K-eigenvectors to Q-eigenvectors, B maps back. (4) Our transfer matrix has EXACTLY this AB structure: A maps subsets to "path ends at vertex", B maps "path starts at vertex" to complement subsets. Transfer matrix symmetry (INV-001) = hidden detailed balance condition.
**Action needed:** Try to formalize the "hidden detailed balance" — identify the group action and show it satisfies Feng's reversibility conditions. This could prove INV-001.

### INV-017: El Sahili & Abi Aad (2020) — parity of paths in tournaments
**Source:** tex bibliography, \cite{elsahili2020}
**Status:** Referenced in tex for decisive/concordant classification. NOT deeply investigated for connections to OCF.
**What:** Discrete Math 343 (2020), Art. 111695. Mod-4 congruences.
**Action needed:** Read the paper. Check if their mod-4 results constrain or relate to our alpha_1 = |C_3| (mod 2) conjecture (INV-011).

### INV-018: El Sahili & Ghazo Hanna (2023) — number of Ham paths/cycles
**Source:** tex bibliography, \cite{elsahili2023}
**Status:** Referenced. NOT investigated for OCF connections.
**What:** J. Graph Theory 102 (2023), 684-701. About the number of oriented Hamiltonian paths and cycles in tournaments.
**Action needed:** Read the paper. They study H(T) directly — any bounds or structural results could inform OCF.

### INV-019: Schweser-Stiebitz-Toft (arXiv:2510.10659) — Redei revisited
**Source:** tex bibliography, \cite{schweser2025}, paper-deep-connections.md Section 1
**Status:** INVESTIGATED (opus-S5). Key theorems extracted.
**What:** Redei's Stronger Theorem (Thm 1.1): add non-oriented edges to tournament, #Ham paths beginning AND ending in tournament vertices is EVEN. Berge's Stronger Theorem (Thm 1.2): G and G-bar have same Ham path parity. Dirac's Stronger Theorem (Thm 2.1): inclusion-exclusion on edge subsets.
**Key findings:** (1) Direct connection to Open Problem 4 (mixed graphs). Non-oriented edges DOUBLE insertion opportunities in Q-Lemma. (2) Berge's theorem gives H(T) ≡ H(T^op) (mod 2) for tournaments, and constrains I(Omega(G),2) under complementation. (3) Strategy for extending Q-Lemma: verify computationally that inshat remains odd for mixed graphs.
**Action needed:** Test inshat parity for mixed graphs computationally. If confirmed, Q-Lemma proof extends directly.

### INV-020: Striker (2011) — unifying poset perspective
**Source:** tex bibliography, \cite{striker2011}
**Status:** Referenced in Open Problem 5. NOT investigated.
**What:** Adv. Appl. Math. 46 (2011), 583-609. Connects ASMs, plane partitions, Catalan objects, tournaments, tableaux via posets.
**Action needed:** Read the paper. Check S3-equivariance (INV-008). The poset perspective may give new structural insights for OCF.

### INV-021: Chapman (2001) — alternating sign matrices and tournaments
**Source:** tex bibliography, \cite{chapman2001}
**Status:** Referenced in Open Problem 5. NOT investigated.
**What:** Adv. Appl. Math. 27 (2001), 290-298. Direct connection between ASMs and tournaments.
**Action needed:** Read the paper. The ASM connection could provide algebraic tools (determinantal formulas, etc.) for H(T).

### INV-022: Eplett (1979) — self-converse tournaments
**Source:** tex bibliography, \cite{eplett1979}
**Status:** Referenced briefly. NOT investigated for OCF.
**What:** Canad. Math. Bull. 22 (1979), 23-27. Self-converse tournament counts.
**Action needed:** Check if self-converse tournaments have special OCF properties. The BlackSelf(8) class (INV-012) is self-converse.

### INV-023: Forcade (1973) — parity of paths and circuits
**Source:** tex bibliography, \cite{forcade1973}
**Status:** Referenced heavily. Our paper gives a NEW combinatorial proof of his F2-invariance.
**What:** Discrete Math 6 (1973), 115-118. Original F2-invariance proof via generating functions.
**Action needed:** Compare his generating function approach to our subset convolution. His GF machinery may contain seeds of an OCF proof.

---

## Priority D: Computational targets

### INV-024: H(T_p) for Paley primes — EXTENDED
**Source:** OPEN-Q-013, opus-S9
**Status:** COMPUTED through p=23.
**Results:**
- H(P(3))=3, H(P(7))=189, H(P(11))=95095, H(P(19))=1,172,695,746,915, H(P(23))=15,760,206,976,379,349
- All match OEIS A038375 where known (a(3)=3, a(7)=189, a(11)=95095)
- P(7) confirmed as GLOBAL maximizer by exhaustive check of all 2^21 n=7 tournaments
- |Aut(T_19)| = 171, H/|Aut| = 6,857,869,865
- Ratio H(P(p))/(p!/2^{p-1}): 2.000, 2.400, 2.440, 2.527, 2.557 — converging toward e=2.718
**Sequence:** H/|Aut| = 1 (p=3), 9 (p=7), 1729 (p=11), 6857869865 (p=19).
**Next step:** Submit H(P(p)) values and H(P(p))/|Aut| sequence to OEIS. Compute H(P(31)) if feasible (2^31*31 DP — ~66B ops, might take hours).

### INV-025: Integrality conjecture C(p,k) | c_k(T_p) for k >= (p+1)/2
**Source:** T036/T153 (tangents), OPEN-Q-013 table
**Status:** Observed for p=11. NOT tested for p=7 or p=19.
**What:** For Paley primes p = 3 mod 4, the cycle count c_k(T_p) is divisible by C(p,k) when k >= (p+1)/2.
**Next step:** Verify for p=7 (compute all c_k). Then prove using Aut(T_p) symmetry.

### INV-026: Alpha_1 vs |C_3| mod 2 — systematic test
**Source:** INV-011, oq:mod4_struct
**Status:** TESTED (opus-S5). CONJECTURE IS FALSE.
**Result:** Counterexamples at every n tested:
  - n=3: 2/8 counterexamples (the 3-cycle tournaments have c3=1 odd, alpha_1=0 even)
  - n=4: 16/64 counterexamples (R-cone and near-R-cone tournaments)
  - n=5: 384/1024 counterexamples
  All counterexamples have alpha_1=0 but c3 odd. The conjecture fails because alpha_1 counts independent PAIRS of odd cycles in Omega(T), which is 0 whenever #cycles <= 1.
**Impact:** Open Problem 9 needs reformulation. Alpha_1 ≠ c_3 mod 2 in general.

### INV-027: Realizable conflict graphs catalog
**Source:** INV-013, conflict_graph_catalog.py
**Status:** DONE (opus-S5). Major structural finding.
**Results:**
  - n=3: 2 distinct Omega structures. n=4: 3. n=5: 6. n=6: 24. n=7 (sampled): 172.
  - **Omega(T) is ALWAYS PERFECT** (exhaustive n<=6, 2000 random n=7). This is a significant constraint — independence number = clique cover number.
  - Omega(T) is NOT always chordal (14% non-chordal at n=6, 12% at n=7).
  - For n<=5, Omega is always complete (pigeonhole: any two subsets of size>=3 share a vertex in a 5-element set).
  - At n=6, non-edges correspond exclusively to complementary 3-cycles (vertex sets partition {0,...,5}).
  - Omega can be disconnected at n=6 (80/32768 tournaments, always exactly 2 complementary 3-cycles).
**Impact:** Perfectness of Omega(T) constrains the independence polynomial and could simplify OCF proof strategies. For perfect graphs, Lovasz theta = clique number, and the fractional chromatic number = chromatic number.

### INV-028b: Fix DR mod-4 proof (Thm 7.4 in tex)
**Source:** tex-deep-analysis.md (ISSUE-1)
**Status:** Proof is BROKEN (arithmetic produces v_2 = -2). Result verified for n=3,7,11 only.
**What:** The proof attempts Moon's formula arithmetic but fails. Need proper v_2 analysis using Kummer's theorem, or prove via alpha_1 parity directly (not just |C_3|).
**Next step:** Compute alpha_1 mod 2 for DR_n using OCF. Possibly downgrade to "Verified Conjecture" in tex.

### INV-029b: Fix SE-SYT formula (Thm 7.3 in tex)
**Source:** tex-deep-analysis.md (ISSUE-2)
**Status:** Classical formula cited gives non-integer (2^{3/2} for m=2). Result verified n=5,7.
**What:** Find correct classical reference for SE-SYT count on 2-core shapes. Likely Stembridge (1996) or similar.
**Next step:** Look up Stembridge's "Canonical bases and self-evacuating tableaux." Give clean proof or correct citation.

### INV-030b: Pin grid S_3 symmetry for OCF
**Source:** tex-deep-analysis.md (Section E)
**Status:** NOT explored.
**What:** The S_3 action on barycentric coordinates constrains the polynomial identity. Can it reduce the proof of delta_H = delta_I by exploiting the 6-fold symmetry? The subset convolution lives on Boolean lattice 2^{others} which is a sublattice of the pin grid.
**Next step:** Check if delta_H = delta_I as polynomial has S_3 symmetry. If so, proving it on a fundamental domain suffices.

---

## Priority E: Tangents needing investigation

### INV-028: Hard-core lattice gas at fugacity 2
**Source:** T006, hard_core_lattice_gas.py, hard_core_fast.py
**Status:** INVESTIGATED (opus-S5). Key finding: non-perturbative regime.
**What:** H(T) = I(Omega(T), 2) = Z(Omega(T), lambda=2). Lambda=2 is ABOVE all cluster expansion convergence thresholds for any max degree Delta >= 2:
  - Shearer bound: 1/(Delta-1) << 2 for Delta >= 2
  - LLL/tree bound: (Delta-1)^{Delta-1}/Delta^Delta << 2 for Delta >= 2
  - Kotecky-Preiss: 1/(e*(Delta+1)) << 2 for Delta >= 1
  This means OCF is a non-perturbative identity — standard polymer expansion / cluster expansion methods CANNOT prove it.
**Omega(T) structure (n=4,5):** #cycles dist ranges from 0 to 6 (n=5). Max degree of Omega grows with n. Density is moderate. Independence number = 1 for all n=4 tournaments with cycles (all cycles share vertices).
**Impact:** Rules out perturbative approaches. OCF requires exact cancellations, not convergence arguments.

### INV-029: Ballot sequence / Dyck path connection
**Source:** T001, OPEN-Q-005, ballot_sequence_test.py
**Status:** RESOLVED (opus-S5). Bijective proof FOUND.
**What:** C(L-2, 2k-1) counts signatures with exactly k Type-II positions in an L-cycle window.
**Bijective proof:** The L-cycle through v has L-1 non-v vertices, giving L-1 signature values (s_1=1 forced, s_{L-1}=0 forced, L-3 free). There are L-2 consecutive pairs. Define transition indicators t_j = (s_j != s_{j+1}). Since s starts at 1 and ends at 0, total transitions must be odd. Transitions alternate fall-rise-fall...fall, so k Type-II = (2k-1 transitions + 1)/2. Choosing which 2k-1 of the L-2 positions are transitions gives C(L-2, 2k-1). QED.
**Convention note:** Initial attempt with wrong convention (L-4 free vars, sig length L-2) gave C(L-3, 2k-1). Correct convention: L-1 non-v vertices, sig length L-1, L-3 free vars, L-2 pairs.

### INV-030: Tower hypothesis (L-cycle corrections from (L+2)-cycles)
**Source:** T012, OPEN-Q-012
**Status:** Hypothesis only. NOT tested.
**What:** At n=2k, the first cycle with mu>1 has length 2k-1. Excess from shorter cycles may be compensated by (L+2)-cycle contributions. Is there a recursive tower structure?

### INV-031: Lindstrom-Gessel-Viennot (LGV) approach to bijection
**Source:** T046
**Status:** Idea only. NOT attempted.
**What:** The bijection between Ham paths and 2-colored cycle sets, if it exists, might require a global construction like LGV lattice path counting. The non-local nature of the correspondence (T035 dead end) suggests a determinantal approach.

### INV-036: Tiling grid geometry and class structure
**Source:** opus-2026-03-06-S1 (deep tiling investigation)
**Status:** INVESTIGATED. Key structural findings.
**What:** How does the {0,1}^m tiling space geometry relate to tournament isomorphism classes?
**Results:**
- **Sigma (converse) acts cleanly on classes:** sigma permutes bits (no complement), preserves weight. Self-converse: 2,2,8,12 classes at n=3,4,5,6. Sigma-fixed tilings = 2^floor((n-1)^2/4).
- **Complement does NOT respect classes:** unlike sigma, flipping all non-path arcs does not map classes to classes.
- **Standard invariants almost distinguish:** At n=6, (score, c3, c5, omega_deg, H) fails only for sigma pairs (converse-paired classes) plus occasional self-converse coincidences.
- **Triangle 3-cycle probability:** P=1/2 for consecutive triples (path arcs), P=1/4 for all others. E[c3] = (C(n,3) + n-2)/4.
- **Strong H~c3 correlation:** r=0.956 at n=5,6. H = 1+2c3 exact at n<=4, breaks at n>=5.
- **Bit-position variance:** Longest arc (gap=n-1) most predictive of class. Middle arcs vary most.
- **Class transition graph:** Always connected. ΔH always even. E[ΔH]=0 for every arc position.
- **Weight distributions distinguish:** Can separate classes sharing all tournament invariants.
**Full writeup:** `03-artifacts/drafts/tiling-symmetry-analysis.md`
**Scripts:** `04-computation/tiling_*.py` (5 files)
**Next step:** (1) Investigate which tiling properties predict H beyond c3. (2) Connect sigma reduction to arc-flip proof strategy. (3) Look for grid-local rules that determine class.

### INV-037: Pin-grid sigma vs tournament sigma — two-sigma structure
**Source:** opus-2026-03-06-S2 (sigma structure investigation)
**Status:** INVESTIGATED. Clean structural results, but no proof path yet.
**What:** The pin-grid sigma (r,c)->(c,r) and tournament sigma (i,j)->(n-1-j,n-1-i) are DIFFERENT symmetries. Pin sigma acts within strips; tournament sigma acts across strips. They agree only on diagonal r=c.
**Key results:**
- **POS-free identity:** free(strip k) = cumul_POS(k) = floor(k/2). Growth rate: delta_free(k) = POS(k) = [k even].
- **n->n+2 structure:** Adds strips n and n+1 with exactly n sigma-free bits and exactly 1 POS (midpoint arc).
- **Tournament sigma always preserves H** (converse operation, verified n=3,...,7).
- **Pin-grid sigma does NOT preserve H** in general (only 5% at n=7).
- **Two sigmas don't commute;** composition has order 3; generate S_3-like group.
- **Mod-4 structure:** Neither sigma preserves H mod 4 reliably.
**Scripts:** `04-computation/sigma_structure.py`
**Next step:** (1) Understand algebraic significance of the S_3 group. (2) Can the n->n+2 POS structure be used differently (not through H preservation)? (3) Relate to transfer matrix symmetry (INV-001).

### INV-039: SC Maximizer Theorem and sigma* structure
**Source:** kind-pasteur-2026-03-06-S18/S18e, T091, T093, T095, OPEN-Q-016
**Status:** VERIFIED exhaustive n=4,5,6,7. Mechanism identified. NOT proved.
**What:** Within each self-complementary score class, max H is always achieved by SC tournament. The mechanism: involutory anti-automorphism sigma induces sigma* on directed odd cycles, which is an involutory automorphism of Omega(T). At even n, sigma* is fixed-point-free, pairing all cycles. Some pairs are vertex-disjoint (giving alpha_2 contributions). At even n, sigma is fixed-point-free on vertices (proved: fixed point implies score=(n-1)/2, non-integer).
**Key results:**
- n=7: all 15 self-comp score classes confirmed
- sigma* is ALWAYS involution and ALWAYS preserves Omega conflict graph
- At n=6: SC H=45 has 4 disjoint 3-cycle pairs; NSC H=43 has only 1
- Two routes to max H: more cycles (alpha_1) OR more disjoint pairs (alpha_2)
- Involutory anti-aut always exists (even when non-involutory ones also exist)
**Scripts:** sc_maximizer_n7_fast.py, anti_aut_analysis.py, sigma_star_analysis.py, anti_aut_involution_test.py
**Draft:** sc-maximizer-mechanism.md
**Next step:** (1) Test at n=8 (even n, computationally expensive but important). (2) Prove algebraically that sigma* structure maximizes I(Omega,2). (3) Connect to Kahn-Galvin-Tetali or Godsil's matching polynomial results for graphs with involutory automorphisms.

### INV-038: Blueself parity theorem and census structure
**Source:** opus-2026-03-06-S3 (deep census investigation)
**Status:** THM-023 PROVED. Census in progress through n=8.
**What:** Blueself (GS + self-flip) exists if and only if n is even. Proved algebraically: flip changes endpoint scores by score'(0) = n - score(0), so same-score requires score(0) = n/2 (integer only at even n).
**Census results (exhaustive n=3,...,6, in progress n=7,8):**
- POS orientation is perfectly UNIFORM: each pattern gets exactly 2^(m-#POS) tilings
- GS POS is also perfectly UNIFORM
- SC always maximizes H within each score sequence class (confirmed with kind-pasteur findings)
- Blueself at n=4: H=5 (rank 1/4), n=6: H=41,45 (ranks 5,1/56) — near or at global maximum
- Blackself at odd n is in SC classes; at even n exclusively in NSC (paired) classes
- SF tilings come in flip-pairs; SF count per class is 2 at n=6, 4 at n=5
- Self-flip fraction decreases: 25%, 12.5%, 1.56% at n=4,5,6
**Scripts:** `04-computation/deep_census_analysis.py`, `04-computation/pos_tiling_census.py`, `04-computation/census_n8.py`
**Theorem:** `01-canon/theorems/THM-023-blueself-parity.md`
**Next step:** (1) Complete n=7 and n=8 census. (2) Investigate why blueself achieves max H. (3) Count blueself at n=8 (1280 eligible GS tilings, need canonicalization).

---

## Completed / Closed investigations

- [DONE] OCF verified n<=8 exhaustive (opus-S3/S4)
- [DONE] Transfer matrix symmetry discovered and verified (opus-S4b)
- [DONE] Locally transitive tournaments tested — DO have 5/7-cycles (opus-S4b)
- [DONE] Feng + Rajkumar connections documented in paper-connections.md (opus-S4b)
- [DONE] T_11 cycle table complete, H(T_11)=95095 confirmed (kind-pasteur-S2/S5)
- [DONE] Per-path identity failure characterized (THM-009)
- [DONE] Even-odd split is consequence not equivalent to OCF (MISTAKE-008)
- [DONE] Bracket structure B(u,w) analyzed (T047, bracket_structure.py)
- [DONE] H(T_19) computed: 1,172,695,746,915; H/|Aut|=6,857,869,865 (opus-S5)
- [DONE] Deep paper analysis: SST, Rajkumar, Feng — all key theorems extracted (opus-S5)
- [DONE] Ballot sequence bijective proof for C(L-2, 2k-1) (opus-S5)
- [DONE] Hard-core lattice gas: lambda=2 is non-perturbative regime (opus-S5)
- [DONE] Alpha_1 ≡ c_3 (mod 2) conjecture DISPROVED (opus-S5)
- [DONE] Conflict graph catalog: Omega(T) is always PERFECT (opus-S5, exhaustive n<=6)
- [DONE] Omega(T) is always CLAW-FREE (opus-S5, exhaustive n<=6, sampled n=7,8)
- [DONE] Web research: 9 new connections documented in web-research-connections.md (opus-S5)
- [DEAD] Per-vertex decomposition of unmatched counts (T045)
- [DEAD] Cycle bijection under arc reversal (MISTAKE-005)
- [DEAD] Contiguous block decomposition (T035)
- [DEAD] Contraction approach (T017)
- [DEAD] Alpha_1 = c_3 (mod 2) conjecture — counterexamples at all n (opus-S5)
- [DONE] Paley maximizer conjecture verified at p=3,7,11 (exhaustive at p=7), extended with H(P(19)), H(P(23)) (opus-S9)
- [DONE] n=8 H-maximizer identified: H=661=a(8), self-converse, |Aut|=1, does NOT contain P(7) (opus-S9)
- [DONE] Full Omega structure at n=8: 76-78 vertices, density 0.98, 70-75% of H from 5/7-cycles (opus-S9)
- [DONE] Ratio H(P(p))/(p!/2^{p-1}) converges toward e: 2.00, 2.40, 2.44, 2.53, 2.56 for p=3,7,11,19,23 (opus-S9)
