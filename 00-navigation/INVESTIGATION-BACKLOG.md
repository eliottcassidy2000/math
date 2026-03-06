# Investigation Backlog

**Purpose:** Systematic catalog of every lead, reference, connection, and unexplored direction extracted from the repo. Claude agents should consult this before choosing what to work on, and add new leads as they emerge. Prioritized by potential impact on proving OCF (Claim A).

**Last full repo scour:** opus-2026-03-06-S6
**Last web research:** opus-2026-03-06-S5 (deep web synthesis — Hopf algebra, Feng reversibility, DRT, Irving-Omar)

---

## Priority A: Key structural questions (OCF PROVED by Grinberg-Stanley)

### INV-032: Omega(T) structural properties — PARTIALLY DISPROVED
**Source:** Web research opus-S5, opus-S7 (disproof), opus-S9 (line graph disproof), opus-S10 (structure analysis)
**Status:** DISPROVED: Omega(T) is NOT always claw-free (fails n=9, 90%) or perfect (fails n=8, 53.8%). NOT a line graph (K_5-e found at n=6, 45%). S_{1,1,1}-free through n=11, fails n=12.
**What remained true:** All-real-roots of I(Omega(T), x) appears to hold even for imperfect/non-claw-free Omega (tested n<=10, 0 failures). This is a deep structural conjecture NOT explained by any forbidden subgraph property.
**Note:** OCF is now proved by Grinberg-Stanley, so this is no longer a proof strategy — it's a structural question. Real-rootedness explanation must be algebraic (Irving-Omar/Grinberg-Stanley symmetric function framework).
**Extended testing (opus-S18):** Real-rootedness tested for I(Omega_3(T), x) at n=9-21 with 0 failures across 1470+ samples (degrees up to 5). Log-concavity and Newton's inequalities hold in all cases. The "Omega_3 complement = matching" structure holds exhaustively at n≤6 (31088/31088) but fails at n≥7 (75.3%).
**Turán-based proof for n≤11:** At n=9-11, alpha(Omega_3) = 3, so the disjoint-pair graph is triangle-free. Turán gives a2 ≤ c3²/4, proving Newton's first inequality a1² ≥ 3a2. Combined with the degree-3 discriminant bound, this could give a complete proof at n≤11. For n≥12, Turán alone fails.
**Next step:** (1) Complete Turán+discriminant proof for n=9-11. (2) Find tournament-specific bounds on a2 for n≥12. (3) Investigate Irving-Omar determinantal formula for algebraic proof.

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

### INV-044: Hereditary Maximizer Chain — CORRECTED (regular-only at odd n)
**Source:** kind-pasteur-2026-03-06-S18f, S18g (correction), T104, T105, MISTAKE-010
**Status:** CORRECTED. Only REGULAR maximizers at odd n are hereditary.
**What:** Previous claim "all maximizers at odd n hereditary" was WRONG. Exhaustive check:
- n=3: 2/2 hereditary (all regular)
- n=4: 0/24 hereditary (non-regular)
- n=5: 24/64 hereditary (only 24 regular, NOT 40 with score (1,2,2,2,3))
- n=6: 0/480 hereditary (non-regular)
- n=7: 240/240 hereditary (all regular)

**R-minimization theorem (NEW, S18g):** The H-maximizer minimizes R(T) = sum_v H(T-v)/H(T) among all tournaments. Proved algebraically: R(T) = n - E_weighted[|U(S)|] where E is over independent sets of Omega(T) weighted by 2^{|S|}, and |U(S)| = total vertices covered. Confirmed exhaustively n=3-6, checking n=7.

**Key insight:** Being hereditary (R = n*H_{n-1}/H_n for regular maximizers) is NOT the same as minimizing R. The non-regular n=5 maximizers have LOWER R (1.4) than regular ones (5/3 ≈ 1.667) despite not being hereditary.
**Next step:** (1) Verify R-minimization at n=7 (running). (2) Prove R-minimization from OCF. (3) Test if regular n=9 maximizers are hereditary.
**Scripts:** `04-computation/hereditary_maximizer.py`, `04-computation/hereditary_correction.py`, `04-computation/R_minimization_proof.py`, `04-computation/R_min_n7_check.py`

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

### INV-001: Prove transfer matrix symmetry for all n — MAJOR PROGRESS
**Source:** T045, T103 (tangents), symmetry_check.py, symbolic_symmetry_proof.py
**Status:** PROVED SYMBOLICALLY at n=4,5,6,7 (exact polynomial identity). Verified numerically n=4,...,8 (7500+ tests).
**What:** M[a,b] = sum_S (-1)^|S| E_a(S)*B_b(M\S) is always symmetric. This is STRONGER than the even-odd split.
**BREAKTHROUGH (opus-S4):** M[a,b]-M[b,a] = 0 as a polynomial in the arc variables t_{ij} AFTER applying the tournament constraint T[j,i]=1-T[i,j]. With independent arc variables the difference is NONZERO (12 terms at n=4, 48 at n=5). The tournament constraint is essential and sufficient.
**Equivalent formulation:** M_{T^op} = (-1)^{n-2} M_T (converse identity). Combined with path reversal M_{T^op}[i,j]=(-1)^{n-2}M_T[j,i], gives symmetry.
**Key insight:** Connects to Feng's Dual Burnside (Q=AB symmetric under detailed balance). The tournament constraint T[x,y]+T[y,x]=1 plays the role of the "detailed balance" condition.
**New findings (opus-S6):**
- **THM-027 PROVED:** Trace formula tr(M) = H(T) for odd n, 0 for even n. Clean bijection proof via (-1)^{pos(a,P)} formula for diagonal entries.
- **MISTAKE-011:** Old claim M = [[1,0],[0,-1]] always is FALSE (2199/2500 failures at n=4). M entries range from -3 to +3.
- **Off-diagonal sum:** sum_{a≠b} M[a,b] = 0 (odd n), 2*H(T) (even n). Verified n=3,...,7 but NOT yet proved.
- **Complement pairing D(S)+D(U\S) is constant at n=4 but NOT at n=5**, ruling out the simplest telescoping argument.
- **Cauchy-Binet decomposition:** M = E^T * Lambda * B where E[S,v]=E_v(S), B[S,v]=B_v(U\S), Lambda=diag((-1)^|S|). Symmetry equivalent to E^T*Lambda*B = B^T*Lambda*E.
**Next step:** Find a CONCEPTUAL proof that works for all n. Possible approaches: (1) Sign-reversing involution on subsets using tournament constraint. (2) Determinantal identity. (3) Induction using the cancellation structure observed at n=4 (terms factor into X and -X via T[c,d]+T[d,c]=1). **(4) NEW (opus-S5): Feng's dual Burnside reversibility (see INV-045). (5) NEW: Irving-Omar det/per formula (see INV-046). (6) NEW: Hopf algebra comultiplication self-duality (see T114).**
**Scripts:** `04-computation/symbolic_symmetry_proof.py`, `04-computation/transfer_symmetry_analysis.py`

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

### INV-042: Fano-Paley design structure and alpha_2 — MAJOR PROGRESS
**Source:** T102, opus-2026-03-06-S4; T114, kind-pasteur-2026-03-06-S18h
**Status:** PROVED at n=7. BIBD arrangement MINIMIZES alpha_2 but MAXIMIZES H.
**What:** The cyclic triples of Paley T_p form a 2-(p, 3, (p+1)/4) BIBD. At p=7: lambda=2, 7 disjoint pairs = MINIMUM among all regular tournaments. The BIBD is Aut-transitive.
**CRITICAL CORRECTION (S18h):** Previous hypothesis was that BIBD maximizes alpha_2, driving H. This is WRONG. The BIBD actually MINIMIZES alpha_2 (7 vs 10 or 14 for other regular tournaments). But H-maximization is driven by alpha_1 (total DIRECTED odd cycles), not alpha_2. The BIBD forces every 5-vertex subtournament to be regular T_5 (2 directed Ham cycles each), giving 42 directed 5-cycles vs 28-36 for non-BIBD. Combined: alpha_1=80, alpha_2=7, H=189 vs alpha_1=59, alpha_2=14, H=175. Three rigid classes at n=7 (THM-027).
**Formula:** For regular tournaments, D = C(b,2) - p*C(r,2) + sum C(lambda_e, 2). BIBD minimizes the convex sum by Jensen's inequality.
**Next step:** (1) Verify at p=11: does BIBD also maximize directed 5-cycles? (2) Prove that BIBD forces subtournament regularity. (3) Can we prove alpha_1 maximization from BIBD structure at general p?

### INV-043: Paley deletion extended to p=19
**Source:** T104, opus-2026-03-06-S4
**Status:** COMPUTED. Consistent with conjecture.
**What:** H(T_19 - v) = 117,266,659,317 for all vertices v. Scores: (8^9, 9^9), self-complementary. H(T_19)-H(T_19-v) = 1,055,429,087,598 = 2*527,714,543,799.
**Conjecture:** a(18) = 117,266,659,317 in OEIS A038375. Verified chain: a(2)=1=H(T_3-v), a(6)=45=H(T_7-v), a(10)=15745=H(T_11-v).
**Cannot verify** against OEIS (only goes to a(11)=95095).
**Next step:** If someone computes a(18), compare. Or prove Paley deletion gives maximizer using design theory + SC maximizer mechanism.

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
**Source:** kind-pasteur-2026-03-06-S18/S18e, opus-2026-03-06-S4, OPEN-Q-016
**Status:** VERIFIED exhaustive n=4,5,6,7. Mechanism deeply analyzed. NOT proved.
**What:** Within each self-complementary score class, max H is always achieved by SC tournament. The mechanism: involutory anti-automorphism sigma induces sigma* on directed odd cycles, which is an involutory automorphism of Omega(T). At even n, sigma* is fixed-point-free, pairing all cycles. Some pairs are vertex-disjoint (giving alpha_2 contributions). At even n, sigma is fixed-point-free on vertices (proved: fixed point implies score=(n-1)/2, non-integer).
**Key results (opus-S4 deepened):**
- **PROVED: sigma always induces Omega automorphism** (clean proof: sigma maps directed C to reverse of sigma(C), preserving vertex-sharing)
- **PROVED: At even n, sigma* has NO fixed cycles** (3-cycles: can't fix set of 3 with fpf involution; 5-cycles at n=6: can't fix set of 5 with 3 two-cycles)
- SC and NSC have the SAME number of 3-cycles within score class — the difference is in ARRANGEMENT (disjoint pairs vs all overlapping)
- SC max alpha_2 >= NSC max alpha_2 within every score class at n=6
- **Path reversal identity:** M_{T^op}[i,j] = (-1)^{n-2} M_T[j,i] (proved)
- **At odd n=5:** alpha_2 = 0 for ALL SC tournaments (fixed point forces all cycles to overlap)
- **At odd n=7 (Paley):** 21 anti-auts, 7 involutions (one per fixed point), each finds 1 disjoint 3-cycle pair
**Scripts:** sc_maximizer_n7_fast.py, anti_aut_analysis.py, anti_aut_exhaustive.py, clique_antiaut_connection.py
**Draft:** sc-maximizer-mechanism.md
**Next step:** (1) Test at n=8 (even n). (2) Prove SC always achieves max alpha_2 within score class. (3) Formalize the "arrangement advantage" into an algebraic proof.

### INV-040: Paley deletion gives H-maximizer
**Source:** kind-pasteur-2026-03-06-S18e, opus-2026-03-06-S4
**Status:** VERIFIED at p=3,7,11. Conjecture.
**What:** Deleting any vertex from Paley tournament T_p gives a tournament with H = max H at n=p-1 (= OEIS A038375(p-1)).
**Results:**
- T_3 → T_2: H=3 → H=1 = a(2) ✓
- T_7 → T_6: H=189 → H=45 = a(6) ✓
- T_11 → T_10: H=95095 → H=15745 = a(10) ✓
- All vertex deletions give the same H (by Aut(T_p) transitivity)
- T_11 - v has self-complementary scores (4,4,4,4,4,5,5,5,5,5)
- H(T_p) - H(T_p-v) = 2 * (sum of mu-weighted cycles through v) (Claim A)
**Conjecture:** T_p - v is the GLOBAL H-maximizer at n = p-1 for all Paley primes p ≡ 3 mod 4.
**Next step:** (1) Verify at p=19 (need H(T_19 - v), n=18 DP ~2^18*18 ~ 5M, feasible). (2) Test whether T_p - v is SC. (3) Relate to lattice theory or QR structure.

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
- [DONE] Deep web synthesis: Irving-Omar, Grinberg-Stanley, Grujić-Stojadinović, Feng, DRT theory (opus-2026-03-06-S5)

---

## Priority F: New leads from web synthesis (opus-2026-03-06-S5)

### INV-045: Hopf algebra route to transfer matrix symmetry
**Source:** T114, Grujić-Stojadinović arXiv:2402.07606, Feng arXiv:2510.25202
**Status:** CONNECTION IDENTIFIED. Most promising algebraic route.
**What:** The Hopf comultiplication Δ([T]) = Σ_S [T|_S]⊗[T|_{V\S}] encodes our subset convolution. Feng's dual Burnside process proves Q=AB is symmetric under detailed balance. The tournament constraint T[x,y]+T[y,x]=1 is the detailed balance condition. If formalized, this proves INV-001 for all n.
**Why it matters:** Transfer matrix symmetry ⟹ OCF (via even-odd split + s-coefficient identity). A Hopf algebra proof would be clean and publication-worthy.
**Next step:** (1) Express our transfer matrix M in Feng's AB framework precisely. (2) Verify detailed balance condition algebraically. (3) Apply Feng Thm 3.3 to conclude symmetry.

### INV-046: Irving-Omar det/per formula → transfer matrix
**Source:** T118, Irving-Omar arXiv:2412.10572 Proposition 2
**Status:** CONNECTION IDENTIFIED. Not attempted.
**What:** ham(D) = Σ_S det(Ā[S])·per(A[S^c]). The walk generating function W_D(z)=det(I+zXĀ)/det(I-zXA) encodes Hamiltonian path structure as ratio of determinants. This may directly relate to our transfer matrix decomposition.
**Next step:** Express our M[a,b] in terms of Irving-Omar's matrix framework. Check if their det/per identity implies M symmetry.

### INV-047: Paley maximizer via DRT theory
**Source:** T116, Reid-Brown 1972, Nozaki-Suda arXiv:1202.5374
**Status:** CLASSICAL EQUIVALENCE + our computational evidence.
**What:** DRTs ↔ skew Hadamard matrices. Paley T_p is the canonical DRT. Nozaki-Suda characterize skew Hadamard via spectra of tournaments of size n-2. Our spectral regularity finding (corr(H,λ₁)=-0.97) + DRT theory could explain why Paley maximizes H.
**Key question:** Among all DRTs on p vertices, does Paley ALWAYS maximize H? Or is this specific to Paley among all tournaments?
**Next step:** (1) Check if non-Paley DRTs exist at small p and compare H values. (2) Relate DRT cycle balance to alpha_k maximization.

### INV-048: Asymptotic convergence H(T_p)/(p!/2^{p-1}) → e
**Source:** T117, Adler-Alon-Ross 2001
**Status:** COMPUTATIONAL EVIDENCE. Not proved.
**What:** Adler-Alon-Ross proved max H(T) ≥ (e-o(1))·n!/2^{n-1} using random regular tournaments. Our Paley ratios 2.00→2.56 suggest convergence to e. Paley tournaments are quasi-random in Chung-Graham-Wilson sense.
**Next step:** (1) Compute H(T_p)/(p!/2^{p-1}) for p=31 if feasible. (2) Check if quasi-randomness implies near-optimal H. (3) Try Stirling approximation on the cycle-count formula.

### INV-049: El Sahili-Ghazo Hanna type-preserving converse symmetry
**Source:** arXiv:2101.00713 (2023), J. Graph Theory 102
**Status:** PUBLISHED RESULT, connection to our work identified.
**What:** T and T^op have the same number of oriented Hamiltonian paths of EVERY type. Our transfer matrix identity M_{T^op} = (-1)^{n-2} M_T is a STRONGER result that implies this. The El Sahili result follows from transfer matrix symmetry + path reversal.
**Impact:** Our transfer matrix results strengthen known literature results.
**Next step:** Note in paper draft. Check if Ai (2025) "New Digraph Polynomials" extends further.

### INV-050: Satake's cyclotomic NDRTs (arXiv:2502.12090, Feb 2025)
**Source:** T116 area
**Status:** IDENTIFIED. Not investigated.
**What:** Nearly-doubly-regular tournaments from almost difference sets. Savchenko's conjecture on canonical spectrum. Under Hardy-Littlewood conjecture F, infinitely many NDRTs with canonical spectrum exist.
**Next step:** Read paper. Check if NDRTs approach Paley's H-maximization. Compare spectra.
