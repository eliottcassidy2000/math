# Investigation Backlog

**Purpose:** Systematic catalog of every lead, reference, connection, and unexplored direction extracted from the repo. Claude agents should consult this before choosing what to work on, and add new leads as they emerge. Prioritized by potential impact on proving OCF (Claim A).

**Last full repo scour:** opus-2026-03-05-S4b

---

## Priority A: Directly blocks or could prove OCF

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
**Source:** paper-connections.md, T046
**Status:** PARTIALLY investigated (opus-S4b). Connections documented. Paper NOT in our bibliography.
**What:** Flip classes, locally transitive = rank 2, R-cones, sign rank. Key results: every T in flip class of R-cone; mu(T) dimension bound; forbidden configurations.
**Action needed:** ADD to bibliography. Develop INV-004 and INV-005 concretely.
**Tested:** Locally transitive tournaments DO have 5/7-cycles (T046). OCF passes 100% for LT, R-cones, automorphism-symmetric tournaments.

### INV-016: Feng (arXiv:2510.25202) — dual Burnside process
**Source:** paper-connections.md, tex line 1600/2151
**Status:** PARTIALLY investigated (opus-S4b). In bibliography but with minimal annotation.
**What:** Q=AB factorization, primal-dual spectral correspondence, lumping theory.
**Connections:** Transfer matrix symmetry (INV-001) resembles detailed balance. The Z_2^m group action on subsets with (-1)^|S| as character connects to Burnside orbit counting.
**Action needed:** Deepen connection. Can the primal-dual eigenvector correspondence (Thm 3.10) be applied to our transfer matrix?

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
**Source:** tex bibliography, \cite{schweser2025}
**Status:** Referenced. Connection to INV-010 (mixed graphs).
**What:** Stronger Redei for mixed graphs. Compatible with Route A (Q-Lemma).
**Action needed:** Read the paper. Check if their strengthening gives a Q-Lemma generalization that constrains OCF.

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
**Status:** NOT computed. Requires 19-vertex Ham path counting.
**What:** Next value in the sequence H/|Aut| = 1, 9, 1729, ??? Could reveal pattern.
**Feasibility:** 2^171 tournaments is impossible. Must compute H(T_19) directly using DP on the specific Paley tournament. n=19 DP is O(2^19 * 19^2) ~ 10^8, feasible.

### INV-025: Integrality conjecture C(p,k) | c_k(T_p) for k >= (p+1)/2
**Source:** T036/T153 (tangents), OPEN-Q-013 table
**Status:** Observed for p=11. NOT tested for p=7 or p=19.
**What:** For Paley primes p = 3 mod 4, the cycle count c_k(T_p) is divisible by C(p,k) when k >= (p+1)/2.
**Next step:** Verify for p=7 (compute all c_k). Then prove using Aut(T_p) symmetry.

### INV-026: Alpha_1 vs |C_3| mod 2 — systematic test
**Source:** INV-011, oq:mod4_struct
**Status:** NOT done.
**What:** For all tournaments at n=5,6,7: compute alpha_1 (number of independent pairs of odd cycles) and |C_3| (number of 3-cycles). Test if alpha_1 = |C_3| (mod 2).

### INV-027: Realizable conflict graphs catalog
**Source:** INV-013
**Status:** NOT done.
**What:** For n=3,...,7, compute Omega(T) for all non-isomorphic tournaments. Catalog which graphs appear. Look for forbidden subgraphs.

---

## Priority E: Tangents needing investigation

### INV-028: Hard-core lattice gas at fugacity 2
**Source:** T006
**Status:** Concept identified. NOT explored.
**What:** H(T) = I(Omega(T), 2) is the partition function of the hard-core model on Omega(T) at fugacity lambda=2. Statistical mechanics tools (transfer matrices, Bethe ansatz, correlation inequalities) could apply.

### INV-029: Ballot sequence / Dyck path connection
**Source:** T001, OPEN-Q-005
**Status:** Formula C(L-2, 2k-1) proved. Bijective proof MISSING.
**What:** The distribution of Type-II counts within an L-cycle window is given by binomial coefficients. This should have a clean bijective proof via ballot sequences.

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
- [DEAD] Per-vertex decomposition of unmatched counts (T045)
- [DEAD] Cycle bijection under arc reversal (MISTAKE-005)
- [DEAD] Contiguous block decomposition (T035)
- [DEAD] Contraction approach (T017)
