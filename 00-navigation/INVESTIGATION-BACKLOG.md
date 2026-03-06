# Investigation Backlog

**Purpose:** Systematic catalog of every lead, reference, connection, and unexplored direction extracted from the repo. Claude agents should consult this before choosing what to work on, and add new leads as they emerge. Prioritized by potential impact on proving OCF (Claim A).

**Last full repo scour:** opus-2026-03-05-S4b
**Last web research:** opus-2026-03-05-S5

---

## Priority A: Directly blocks or could prove OCF

### INV-032: Omega(T) is always claw-free AND perfect — Dyer-Jerrum decomposition
**Source:** Web research opus-S5, arXiv:1909.03414 (Dyer-Jerrum-Müller-Vušković)
**Status:** CLAW-FREENESS VERIFIED computationally (exhaustive n<=6, sampled n=7,8). NOT proved.
**What:** Omega(T) is always (claw, odd-hole)-free. This is exactly the class studied by Dyer-Jerrum, who show the partition function of independent sets can be computed via clique-cutset decomposition into atoms. For claw-free perfect graphs, the structure theorem (Chvátal-Sbihi) decomposes into line graphs of bipartite graphs (computable via permanent) or "peculiar" graphs.
**Why this could prove OCF:** If the clique-cutset decomposition of Omega(T) mirrors the Hamiltonian path decomposition of the tournament, then I(Omega, 2) = H(T) follows from the decomposition structure. The claw-free perfect structure is highly constrained and may force the identity.
**Next step:** (1) PROVE Omega(T) is always claw-free. (2) Study the clique-cutset decomposition of Omega(T). (3) Check if atoms are line graphs of bipartite graphs.

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

### INV-012: BlackSelf(8) exceptional class (Open Problem 7 in paper)
**Source:** oq:n8 in tex
**Status:** Computationally identified. NOT classified.
**What:** Unique isomorphism class at n=8 that is self-converse, has |Aut|>1, |Fix(beta)| odd, but H(T)/|Fix(beta)| is even. Is it related to a Hadamard matrix of order 8 or a skew conference matrix?
**Next step:** Extract the specific tournament. Check against known Hadamard/conference matrices.

### INV-013: Realizable odd-cycle conflict graphs (Open Problem 8 in paper)
**Source:** oq:realizable in tex
**Status:** NOT INVESTIGATED.
**What:** Which graphs G arise as Omega(T) for some tournament T? This characterization could constrain the structure of the independence polynomial and potentially simplify OCF proofs.
**Next step:** Compute Omega(T) for all small tournaments. Catalog which graphs appear. Look for forbidden subgraph characterizations.

### INV-014: 2-adic tower / higher Redei theorems
**Source:** OPEN-Q-008, T007, tex Section 5.5
**Status:** Concept identified. NOT explored computationally.
**What:** I(Omega(T), x) at x=4,8,... gives mod-4, mod-8 invariants of H(T). What is v_2(H(T))? Is there a combinatorial characterization?
**Next step:** Compute v_2(H(T)) for all tournaments at n=5,6,7. Look for patterns related to cycle structure.

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

### INV-024: H(T_19) for the Paley prime p=19
**Source:** OPEN-Q-013
**Status:** COMPUTED (opus-S5). H(T_19) = 1,172,695,746,915.
**Result:** |Aut(T_19)| = 171. H/|Aut| = 6,857,869,865 (exact integer). H is odd. c_3=285, c_5=11628. Per-endpoint count = 61,720,828,785 (same for all 19 endpoints — Paley symmetry). Computed via C DP in 0.5s.
**Sequence:** H/|Aut| = 1 (p=3), 9 (p=7), 1729 (p=11), 6857869865 (p=19).
**Note:** c_3 formula gives 285 but expected 199 — discrepancy needs investigation (formula may be for non-isomorphic rather than labeled 3-cycles).

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
