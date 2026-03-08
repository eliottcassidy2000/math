# Investigation Backlog

**Purpose:** Systematic catalog of every lead, reference, connection, and unexplored direction extracted from the repo. Claude agents should consult this before choosing what to work on, and add new leads as they emerge. Prioritized by potential impact on proving OCF (Claim A).

**Last full repo scour:** opus-2026-03-06-S10
**Last web research:** kind-pasteur-2026-03-07-S33 (Lichiardopol, Chen-Chang disjoint cycles, Frankl Erdos matching)

---

## Priority A: Key structural questions (OCF PROVED by Grinberg-Stanley)

### INV-050: Fourier Decomposition Proof of OCF — OCF PROVED AT n=5 AND n=7
**Source:** opus-2026-03-06-S11b (continued^7, ^8)
**Status:** OCF PROVED AT n=5 AND n=7 via Fourier decomposition. All identities at both n proved.
**What:** OCF decomposes into independent degree-homogeneous Fourier identities:
- **Fourier Homogeneity Theorem:** w_{n-1-2k} is a homogeneous polynomial of degree 2k in centered edge variables s_e = A_e - 1/2.
- **Degree 0:** Trivial (expected values). PROVED for all n.
- **Degree 2:** Proved via proportionality constants c_{2j+1} = C(n-3,2j-2)*(2j-2)!/2^{2j-2}. PROVED for all n.
- **Degree n-1:** Proved via path-cycle bijection: w_0 = 2*[deg-(n-1) of t_n]. PROVED for all n.
- **Degree 4 at n=7:** PROVED via counting lemmas. The degree-4 Fourier space is 2-dimensional:
  - Type P: 5-vertex spanning paths (coefficients ±1 in t5_d4)
  - Type Q: 6-vertex disjoint P₂ pairs (coefficients ±1 in α₂_d4)
  - [deg-4 of t₇] = (1/2)·[deg-4 of t₅] + [deg-4 of α₂] (EXACT, verified all 5985 monomials)
  - w₂/4 = 3·[deg-4 of t₅] + 6·[deg-4 of α₂] (EXACT, verified all 5985 monomials)
  - Counting: c₅=2, c₇=4, paths=12 (type P); c₇=8, paths=24, a₂=1 (type Q)
**Key insight:** Fourier supports of t₅_d4 and α₂_d4 are DISJOINT (spanning paths vs P₂ pairs), reducing the identity to independent counting arguments.
**Scripts:** `04-computation/fourier_homogeneity.py`, `fourier_degree2_identity.py`, `ocf_fourier_proof_framework.py`, `degree4_identity_n7.py`, `degree4_proof_n7.py`
**CRITICAL FINDING (opus-2026-03-07-S37):** At n=9, the degree-4 Fourier space has dimension **>> 200** (no saturation at 200 random tournaments with 300 probe monomials). Within the P4 type alone, coefficients are NOT proportional (117/126 vertex sets have non-constant |coeff|). The |coeff| values range from 1 to 27. This means the n=7 "two type" decomposition DOES NOT generalize. The Fourier approach is **infeasible at n=9** for middle degrees.
**Scripts:** `04-computation/degree4_n9_rank.py`, `degree4_n9_rank2.py`, `degree4_n9_saturation.py`
**Next step:** (1) The Fourier proof cannot extend to n>=9 for middle degrees. Focus on algebraic approaches (OCF already proved by Grinberg-Stanley). (2) The degree-0, degree-2, and degree-(n-1) identities still hold for all n and have clean proofs. Can they be combined differently?

### INV-123: THM-086 Universal Taylor Zeros mod 3 — PROOF SKETCH COMPLETE
**Source:** kind-pasteur-2026-03-07-S37
**Status:** PROOF SKETCH COMPLETE. Verified n=5-10, inductive structure identified.
**What:** c_j(T) = 0 mod 3 for all tournaments T on n vertices and all j < val(n), where val(n) = 2*floor((n-1)/2). This means (x-1)^{val(n)} | F(T,x) mod 3. For n odd, F(T,x) mod 3 is determined by a SINGLE parameter alpha = c_{n-1}(T) mod 3.
**Proved cases:** j=0,1,2 (THM-085, algebraic). j=3 (palindrome + THM-085). j>=4 (DC induction + palindrome, verified computationally).
**Key corollary:** Eulerian conjecture: 3|A(n,k) => 3|F_k(T) for all T. Follows from (x-1)-adic valuation of A_n(x) mod 3 being exactly val(n).
**What remains:** The "almost-tournament claim" — c_j(T\e) = 0 mod 3 for j < val(n)-1 — needs formal proof, likely via nested DC induction. Verified exhaustively at n=5, sampled at n=6-8.
**Scripts:** `04-computation/thm086_verify.py`, `dc_induction_proof.py`, `c4_induction_test.py`, `taylor_cj_mod3_analysis.py`, `eulerian_zeros_from_palindrome.py`
**Next step:** (1) Prove almost-tournament claim algebraically (N_uv formula reduces it to Taylor zeros of the "adjacent pair" polynomial). (2) Extend to mod 9. (3) Mod p for p>=5 INVESTIGATED (S38): universal zeros match Eulerian val for n >= p+2 but Eulerian conjecture FAILS for p>=5 (multiple free parameters).

### INV-124: THM-094 F_k mod 2 Tournament-Independent — PROOF SKETCH COMPLETE
**Source:** kind-pasteur-2026-03-07-S38
**Status:** PROOF SKETCH COMPLETE. Verified exhaustively n<=6, sampled n=7,8.
**What:** F_k(T) = A(n,k) = C(n-1, k) mod 2 for ALL tournaments T. F(T,x) = (1+x)^{n-1} mod 2 is COMPLETELY tournament-independent. Proof via universal Taylor zeros mod 2 (c_j = 0 for j < n-1) + Redei's theorem (F_{n-1} = Hamiltonian path count is always odd). The mod-2 result is the strongest possible: individual F_k are determined, not just linear combinations.
**Key insight:** p=2 is special because (1) val_2(A_n(x)) = n-1 (maximal), giving a single free parameter, and (2) Redei pins that parameter to 1.
**Mod-p generalization (S38):** For p >= 5, universal Taylor zeros match Eulerian valuation only for n >= p+2. The Eulerian conjecture (p|A(n,k) => p|F_k(T)) FAILS for p=5 at n=7 because multiple free parameters in F(T,x) mod 5 allow different zero patterns.
**Scripts:** `04-computation/fk_mod2_proof.py`, `taylor_zeros_mod_p.py`, `mod_p_general_conjecture.py`
**Next step:** (1) Prove universal Taylor zeros mod 2 algebraically (c_j = 0 for j < n-1). (2) Is there an elementary proof not using THM-086 machinery?

### INV-032: Omega(T) structural properties — PARTIALLY DISPROVED
**Source:** Web research opus-S5, opus-S7 (disproof), opus-S9 (line graph disproof), opus-S10 (structure analysis)
**Status:** DISPROVED: Omega(T) is NOT always claw-free (fails n=9, 90%) or perfect (fails n=8, 53.8%). NOT a line graph (K_5-e found at n=6, 45%). S_{1,1,1}-free through n=11, fails n=12.
**DISPROVED (THM-025, opus-S18):** Real-rootedness of I(Omega(T),x) FAILS at n=9. Counterexample: score [1,1,3,4,4,4,6,6,7], I=[1,94,10,1], Newton k=2 fails (100 < 141), complex roots confirmed.
**Structural characterization (opus-S19):** Failure requires (a) three vertex-disjoint 3-cycles partitioning V, AND (b) near-total inter-group domination (9-0, 9-0, 7-2 arc counts), creating a transitivity bottleneck with hub vertex in 92/94 cycles. The extreme tournament (full domination) gives I=(1+x)^3 with disc=0 exactly. One arc flip creates disc<0. Failure is MAXIMALLY RARE: 0 in 10000 random samples at n=9, 0 at n=10,11.
**What remains true:** Real-rootedness holds at n<=8 (Chudnovsky-Seymour, claw-free) and for "generic" tournaments at all n. But it is NOT a universal structural property.
**Next step:** (1) Characterize exact tournament class where failure occurs. (2) Check if the clique-deletion interlacing approach (INV-038) can prove real-rootedness under a claw-free assumption.

### INV-038: Clique-deletion interlacing for Omega(T)
**Source:** opus-2026-03-06-S17, T100, interlacing-clique-deletion.md
**Status:** STRUCTURAL INSIGHT. Proof sketch for n<=8.
**What:** Through-v cycles always form a CLIQUE in Omega(T) (proved: sharing vertex v). Deleting vertex v = deleting this clique from Omega(T). The sequential deletion recurrence I(G,x) = I(G-u,x) + x*I(G\N[u],x) can be applied step-by-step. At n=5, 100% of remaining cycles are adjacent to some through-v cycle (Omega is very dense).
**Key insight:** For n<=8 (claw-free), Chudnovsky-Seymour guarantees each step preserves real roots. For n>=9, real-rootedness can FAIL (THM-025), so the interlacing approach cannot extend to all tournaments. However, 84 claws at n=9 counterexample all share the same 3 leaves — the claw structure is very specific.
**Verification:** 0 failures: n=5 (5120 exhaustive), n=6 (196608 exhaustive), n=7-8 (random).
**Impact:** Proves real-rootedness for n<=8. For n>=9, would need a tournament-specific claw-free condition.
**Next step:** (1) Characterize when Omega(T) is claw-free at n>=9. (2) Check if "generically claw-free" suffices for applications.
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
**Source:** Web research opus-S5, arXiv:2402.07606 (Grinberg), arXiv:2506.08841 (Mitrovic-Stojadinovic)
**Status:** CONNECTION IDENTIFIED. Key bridge found (S36).
**What:** The Redei-Berge symmetric function U_X for digraphs has comultiplication Delta([X]) = sum_S [X|S] tensor [X|V\S] — this IS our subset convolution. The character zeta counts Hamiltonian paths. The antipode S(U_X) = (-1)^|V| U(X-bar) encodes Berge's theorem.
**NEW (Mitrovic-Stojadinovic, arXiv:2506.08841, June 2025):**
  - X_{inc(P)} = omega(U_P): Chromatic function of incomparability graph = omega of Redei-Berge
  - "Converse of Redei": if poset is not a chain, quasi-linear extensions are even
  - Bags-of-sticks decomposition: U_X = sum of simpler digraphs via inclusion-exclusion on edges
  - Stanley-Stembridge connection: e-positivity of X_{inc(P)} <=> h-positivity of U_P
  - Noncommutative deletion-contraction: W_X = W_{X\e} - W_{X/e}^up
  - Mitrovic-Stojadinovic phi(pi) = sum_{gamma X-cycle} (len(gamma)-1) is EXACTLY our S = sum(l_i-1)!
**Verified (S36):** OCF specialization p_1->1, p_{odd>=3}->2, p_{even}->0 gives H(T) from U_T.
**NEW (Mitrovic, arXiv:2504.20968, April 2025):** Noncommutative Redei-Berge function W_X has deletion-contraction: W_X = W_{X\e} - W_{X/e}↑. Thm 3.16: cycle decomposition via inclusion-exclusion over cycle edges. Cor 3.12: tournament formula W_X = Σ(2^{ψ(σ)} p_{Type(σ)}) for odd-cycle permutations = exactly OCF.
**Next step:** (1) Express OCF via bags-of-sticks decomposition. (2) Check if deletion-contraction on W_T gives a direct proof of Claim A. (3) Explore chromatic function connection for imperfect Omega(T). (4) Study Thm 3.16 cycle decomposition for odd cycles.

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

### INV-001: Prove transfer matrix symmetry for all n — PROVED (THM-030)
**Source:** T045, T103 (tangents), symmetry_check.py, symbolic_symmetry_proof.py
**Status:** PROVED FOR ALL n by induction (kind-pasteur-2026-03-06-S25). Verified computationally m=2,...,6 all (a,b) pairs.
**What:** M[a,b] = sum_S (-1)^|S| E_a(S)*B_b(M\S) is always symmetric. This is STRONGER than the even-odd split.
**BREAKTHROUGH (opus-S4):** M[a,b]-M[b,a] = 0 as a polynomial in the arc variables t_{ij} AFTER applying the tournament constraint T[j,i]=1-T[i,j]. With independent arc variables the difference is NONZERO (12 terms at n=4, 48 at n=5). The tournament constraint is essential and sufficient.
**Equivalent formulation:** M_{T^op} = (-1)^{n-2} M_T (converse identity). Combined with path reversal M_{T^op}[i,j]=(-1)^{n-2}M_T[j,i], gives symmetry.
**Key insight:** Connects to Feng's Dual Burnside (Q=AB symmetric under detailed balance). The tournament constraint T[x,y]+T[y,x]=1 plays the role of the "detailed balance" condition.
**c-TOURNAMENT GENERALIZATION (opus-S19b):** Symmetry holds for ALL c-tournaments where t_ij + t_ji = c (any constant c, not just c=1). Verified symbolically n=3,4,5; numerically n=6,7. The constraint must be UNIFORM across pairs (non-uniform c_ij gives 100% failure) and ALL pairs need it. In skew coordinates t_ij = c/2 + s_ij: M is EVEN in s for n even, ODD for n odd. The c^{n-2} coefficient is (n-2)!/2^{n-2} for even n, 0 for odd n. This c-generalization SIMPLIFIES the proof problem: we can work with skew-symmetric part S = A-A^T and ignore the specific value c=1.
**New findings (opus-S6):**
- **THM-027 PROVED:** Trace formula tr(M) = H(T) for odd n, 0 for even n. Clean bijection proof via (-1)^{pos(a,P)} formula for diagonal entries.
- **MISTAKE-011:** Old claim M = [[1,0],[0,-1]] always is FALSE (2199/2500 failures at n=4). M entries range from -3 to +3.
- **Off-diagonal sum:** sum_{a≠b} M[a,b] = 0 (odd n), 2*H(T) (even n). Verified n=3,...,7 but NOT yet proved.
- **Complement pairing D(S)+D(U\S) is constant at n=4 but NOT at n=5**, ruling out the simplest telescoping argument.
- **Cauchy-Binet decomposition:** M = E^T * Lambda * B where E[S,v]=E_v(S), B[S,v]=B_v(U\S), Lambda=diag((-1)^|S|). Symmetry equivalent to E^T*Lambda*B = B^T*Lambda*E.
**PATH REVERSAL PROOF AT c=0 (kind-pasteur-S23):** COMPLETE proof when c=0 (pure skew weights). Path reversal: B_v(S+v) = (-1)^|S| E_v(S+v). This gives M[a,b] = (-1)^{n-2} sum_S E_a(S+a) E_b(R+b) — unsigned, manifestly symmetric by S<->R relabeling. Verified n=3,4,5,6.
**EVEN r-POWERS CONJECTURE (kind-pasteur-S23):** At general c, M(r,s) where r=c/2 has ONLY even r-powers. Equivalent to symmetry. Verified n=3,4,5,6. Path reversal gives B_v(c,s) = E_v(c,-s), which yields M[b,a](c,s) = (-1)^{n-2} M[a,b](c,-s). So symmetry reduces to M having definite s-parity (-1)^{n-2}, i.e., only even r-powers.
**ALGEBRAIC PROOF (kind-pasteur-S23b): M[a,b](-r) = M[b,a](r)** — 5-step proof: T(-r)=-T^T, path reversal under negated transpose, sign bookkeeping, S↔R relabeling. Verified n=4,5,6. This proves the EQUIVALENCE between (i) symmetry, (ii) even-r-powers, (iii) s-parity.
**TOGGLE ANALYSIS (S23b):** At n=4, r^1 monomials cancel pairwise between different S-subsets. At n>=5, cancellation is multi-way (not simple pairwise). No clean single-vertex toggle involution found on whole subsets.
**H(U) MATRIX (S23b, from Kogan/Hamiltonian cycle polynomial):** H(U)_{i,j} = sum of Ham path weights from i to j. Identity: U*H(U)^T = H(U)*U^T. For c-tournaments U+U^T = c(J-I), this gives UH^T = H(cJ-cI-U), but does NOT directly imply H=H^T. Also note: M[a,b] is NOT the same as H(T)_{a,b} (M has inclusion-exclusion signs, H is a direct sum).
**PROVED INDEPENDENTLY by both opus-S25 and kind-pasteur-S25 (THM-030).** Key Identity: odd_r(B_b(W)) = r * col_sum_W(b), equivalently B_b(W)+(-1)^m E_b(W) = 2r*col_sum(b). Inductive proof using column recurrence + first-edge decomposition. The Sigma identity (r*Sigma = odd(T)) follows from summing the inductive hypothesis. The proof closes because odd(sum s*even(B_v)) = 0. Verified computationally m=2..6 all (a,b) pairs. See complete_even_r_proof.py, key_identity_complete_proof.py.
**Scripts:** `04-computation/symbolic_symmetry_proof.py`, `04-computation/transfer_symmetry_analysis.py`
**COEFFICIENT STRUCTURE (opus-S22 continuation):**
- **[r^{n-2}] = (n-2)!** when n even, **0** when n odd. Proof: counting argument, sum_k C(n-2,k)(-1)^k k!(n-2-k)! = (n-2)! * sum_k(-1)^k. Verified n=3,...,6.
- **[r^2] for n=5 = 2·sum_{u∈U}(s_{au}+s_{bu})**. For n=6: degree-2 in s with all coefficients ±2. For n=4: just 2 (constant).
- **[r^1] telescoping (n=4):** Each s_{uv} (u∈U, v∈{a,b}) appears exactly once with + and once with - across subsets. Moving vertex u between S and R flips the sign contribution.
- **M is NOT a cofactor** of A=rJ'+S (exhaustive test n=3,4). Cofactors have degree n-1; M has degree n-2.
- **M is NOT a permanent minor** of A either. The fundamental identity A(-r)=-A^T is clean but M does not decompose as a simple matrix function of A.
- **Key algebraic identity:** A(-r,s) = -A(r,s)^T (since J' symmetric, S skew). Any expression built from A·A^T+A^T·A (which is even in r) could explain the property, but no such expression matching M has been found.
- **Literature update:** Irving-Omar (arXiv:2412.10572), Mitrovic noncommuting Redei-Berge (arXiv:2504.20968) with deletion-contraction W_X = W_{X\e} - W_{X/e}^up, El Sahili-Ghazo Hanna proving T and T^op have same Hamiltonian path type distribution.
**APPROACH RULING (opus-S22 continuation):**
- ❌ Simple cofactor/minor of A (degree mismatch)
- ❌ Permanent of A minor (doesn't match)
- ❌ Adjugate entries of A, I±A, J-A (all fail)
- ❓ Deletion-contraction via Mitrovic noncommuting Redei-Berge (unexplored, most promising NEW lead)
- ❓ Irving-Omar walk generating function det(I+zXĀ)/det(I-zXA) (connection to M unclear)
- ❓ Direct r^1=0 proof via telescoping + induction on n (promising for base case)
**Next step:** (1) Try Mitrovic deletion-contraction approach — express M[a,b] recursively and prove even-r by induction. (2) Understand Irving-Omar matrix formula and whether it encodes M[a,b]. (3) Prove [r^1]=0 directly via the telescoping structure observed at n=4,5. (4) Previous approaches (Feng, Hopf, involution) remain viable but untested.
**Scripts:** `04-computation/symbolic_symmetry_proof.py`, `04-computation/transfer_symmetry_analysis.py`, `04-computation/determinantal_identity_test.py`, `04-computation/det_compare_explicit.py`, `04-computation/r1_coefficient_analysis.py`, `04-computation/r_coefficient_structure.py`

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

### INV-053: Even Cycle Vanishing Theorem — PROVED
**Source:** opus-2026-03-06-S10, T148
**Status:** PROVED. Clean involution argument.
**What:** For any tournament T on [n], p_mu(U_T) = 0 whenever mu has an even part. The proof pairs each permutation sigma with even k-cycle c with sigma' (c reversed); the sign flips because (-1)^{k-1} = -1 for even k. Verified computationally n=3 through n=7.
**Consequences:** (1) U_T lives in the subspace spanned by p_mu with all odd parts — drastically fewer terms. (2) At n=4, only types (1^4) and (3,1) contribute; at n=5, only (1^5), (3,1,1), (5); at n=6, only (1^6), (3,1,1,1), (3,3), (5,1). (3) The Schur expansion simplifies: [s_lambda]U_T = sum over odd-part mu only. (4) This is the SAME T<->T^op involution as in the path reversal proof (T147).
**Connection to INV-001:** The even-r-powers conjecture (kind-pasteur-S23) is the transfer matrix version of this same phenomenon. Both arise from the perpendicular grid symmetry.

### INV-054: Hook Schur Positivity for Tournaments — PARTIAL (fails at n=7)
**Source:** opus-2026-03-06-S10, T149
**Status:** PROVED at n=4 (clean sign argument). VERIFIED at n=5 (11/11), n=6 (40/40). **FAILS at n=7** (231/242, 11 failures all for middle hook (4,1,1,1)).
**What:** [s_{(k,1^{n-k})}]U_T >= 0 holds at n=4,5,6 but fails at n=7.

**n=4 proof:** Only p-types (1^4) and (3,1) matter (by INV-053). All hook characters non-negative at both → sum of non-negative terms.

**n=7 failure mechanism:** chi^{(4,1,1,1)}((7)) = -1 and chi^{(4,1,1,1)}((5,1,1)) = 0. Regular tournament T_7 has 48 directed 7-cycles contributing -48/7, overwhelming positive 3-cycle terms. Result: [s_{(4,1,1,1)}]U_{T_7} = -83/28 ≈ -2.96.

**Which hooks always hold:** Hooks (n) and (1^n) have all-positive characters at odd types → always positive. Hooks (n-2,1,1) and (3,1^{n-3}) also have all-positive chars (at n=7). Only hooks with j odd and j near n/2 can fail.

**Non-hook negativity:** Non-hook chars at (3,1,...,1) are always negative → non-hook coefficients negative for all non-transitive tournaments (verified n=4,5,6).

**Refined question:** For which hooks is positivity universal? Is there a simple characterization (e.g., j even, or |j - n/2| > threshold)?
**Scripts:** `04-computation/schur_hook_analysis.py`, `04-computation/tournament_cycle_structure.py`, `04-computation/even_cycle_vanishing_proof.py`, `04-computation/hook_positivity_n6.py`, `04-computation/hook_positivity_n7.py`

---

## Priority B: Important structural understanding

### INV-007: Odd-cycle bijection (Open Problem 3 in paper)
**Source:** oq:bijection in tex, bijection_search.py, T046
**Status:** Searched computationally at n=3,4. No natural bijection found.
**What:** Construct Phi: Ham(T) -> {(C, f) : C vertex-disjoint odd-cycle collection, f: C -> {0,1}}. The "2-colored cycles" interpretation means each independent set of k cycles contributes 2^k paths. At n=3 with one 3-cycle: 3 paths = 1 + 2. Need to identify which paths correspond to which colored cycle sets.
**Key obstacle:** The correspondence is NOT local/contiguous-block based (T035, confirmed dead end). Must be global.

### INV-008: Striker-Chapman S3-equivariance (Open Problem 5 in paper)
**Source:** oq:striker in tex, \cite{striker2011,chapman2001}
**Status:** INVESTIGATED (kind-pasteur-S22 agent). Question imprecise as stated.
**What:** Is the Striker-Chapman bijection between ASMs and tournaments S3-equivariant under the barycentric identification? Striker (2011) gives a unifying poset perspective connecting ASMs, plane partitions, Catalan objects, tournaments, and tableaux. Chapman (2001) connects ASMs to tournaments directly.
**Finding:** Chapman's bijection maps oriented monotone triangles (not ASMs directly) to tournaments. Striker's tournament bijection uses the disjoint three-color poset {b,r,(g)}, while ASMs use the four-color poset {b,y,o,g} — different subposets. No S3 action is defined in either paper. The question needs precise formulation: (a) define "the bijection" (Chapman's Phi? Striker's poset? composition?), (b) define how S3 acts on both sides. Only then can it be tested at n=3,4.
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

### INV-011: Mod-4 score-sequence criterion (Open Problem 9 in paper) — RESOLVED (NO)
**Source:** oq:mod4_struct in tex
**Status:** RESOLVED NEGATIVELY (kind-pasteur-2026-03-07-S39b).
**What:** Does score sequence determine H(T) mod 4?
**Answer:** NO for n >= 5. The formula H mod 4 = (1 + 2·c3) mod 4 holds only for n <= 4.
**Results:**
- n=3,4: YES (exhaustive, H mod 4 constant within every score class)
- n=5: NO — score (1,2,2,2,3) has H mod 4 in {1,3} (c5 varies within this class)
- n=6: NO — 5/22 score classes have varying H mod 4
- n=7: NO — 27/59 sampled score classes have varying H mod 4
**Key insight:** H mod 4 = (1 + 2·alpha_0) mod 4 where alpha_0 = total odd cycle count.
For n <= 4, alpha_0 = c3 which is score-determined. For n >= 5, 5-cycles contribute to alpha_0 but c5 is NOT determined by score sequence (confirmed independently).
**Also found:** c5 is NOT determined by (score, sum_d², edge_score_sum). Even (score, common_out_neighbor) pairs vary. c5 requires genuine graph structure beyond all local/pairwise statistics.
**Scripts:** `04-computation/mod4_score_test.py`, `04-computation/c5_score_determination.py`

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

### INV-045: Hopf algebra route to transfer matrix symmetry — SUPERSEDED by THM-030
**Source:** T114, Grujić-Stojadinović arXiv:2402.07606, Feng arXiv:2510.25202
**Status:** INVESTIGATED (kind-pasteur-S22 agent). No direct proof obtained.
**What:** The Hopf comultiplication Δ([T]) = Σ_S [T|_S]⊗[T|_{V\S}] encodes our subset convolution. Feng's dual Burnside process proves Q=AB is symmetric under detailed balance. The tournament constraint T[x,y]+T[y,x]=1 is the detailed balance condition.
**Finding:** Feng's framework FAILS because it requires positivity (probability transitions), but our M = E^T * Lambda * B has Lambda = diag((-1)^|S|) with negative entries. Grujic-Stojadinovic gives U_X = U_{X^op} (proved) but this is a GLOBAL identity, not per-entry M[a,b]=M[b,a]. The Hopf comultiplication IS cocommutative but this doesn't directly imply transfer matrix symmetry because E_a and B_b are different types of objects. **Most promising remaining paths:** (1) inductive s-variable approach, (2) Irving-Omar det/per formula (INV-046), (3) "signed Feng" extension (new theorem needed), (4) pointed Hopf algebra tracking distinguished vertices.
**Why it matters:** Transfer matrix symmetry ⟹ OCF (via even-odd split + s-coefficient identity). A Hopf algebra proof would be clean and publication-worthy.
**Next step:** (1) Express our transfer matrix M in Feng's AB framework precisely. (2) Verify detailed balance condition algebraically. (3) Apply Feng Thm 3.3 to conclude symmetry.

### INV-046: Irving-Omar det/per formula → transfer matrix — SUPERSEDED by THM-030
**Source:** T118, Irving-Omar arXiv:2412.10572 Proposition 2
**Status:** SUPERSEDED. Transfer matrix symmetry now proved directly (THM-030, opus-S25) without needing Irving-Omar det/per.
**What:** ham(D) = Σ_S det(Ā[S])·per(A[S^c]). The walk generating function W_D(z)=det(I+zXĀ)/det(I-zXA) encodes Hamiltonian path structure.
**Remaining interest:** Irving-Omar's framework may still provide insight into WHY the Key Identity works — e.g., is there a matrix-algebraic interpretation of B_b + (-1)^m E_b = 2r·col_sum?

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

---

## Priority G: New leads from web search (kind-pasteur-2026-03-06-S19)

### INV-051: Mitrovic noncommuting Rédei-Berge function (arXiv:2504.20968, Apr 2025)
**Source:** kind-pasteur-2026-03-06-S19 web search; opus-2026-03-06-S9 (detailed paper read)
**Status:** DEEPLY READ — HIGH PRIORITY
**What:** Stefan Mitrovic introduces the Rédei-Berge function in NONCOMMUTING variables, which satisfies deletion-contraction (W_X = W_{X\e} - W_{X/e}↑). The commutative version does NOT have deletion-contraction. Key properties: W_X = W_{X^op}, product rule W_{X·Y} = W_X·W_Y. For tournaments: W_X = sum over permutations with all odd cycles of 2^{psi(sigma)} p_{Type(sigma)} with positive integer coefficients.
**Why it matters:** Deletion-contraction enables INDUCTIVE PROOFS. This could provide an inductive framework for OCF or transfer matrix symmetry. The noncommutative structure preserves more information than the commutative version.
**TESTED (kind-pasteur S19):** Direct deletion-contraction does NOT preserve OCF. At n=4:
  - H(T) = H(T\e) - H(T/e): only 18.8% match (DC is for W_X not H)
  - OCF for T\e (non-tournament): only 39.3% hold
  - OCF for T/e (contracted): only 60.7% hold
  OCF is TOURNAMENT-SPECIFIC and does not hold for general digraphs from deletion/contraction.
  The noncommuting framework operates at a different level than H(T).
**DETAILED READING (opus-S9):**
  - Definition 3.1: W_X = sum_{f:V->P} sum_{sigma in Sigma_V(f,X)} delta_f(sigma) x_{f(v1)}...x_{f(vn)} where x_i are NONCOMMUTING. Depends on vertex labeling (unlike commutative U_X).
  - Theorem 3.7 (Deletion-Contraction): W_X = W_{X\e} - W_{X/e}↑ where e=(v_{n-1},v_n). The ↑ operation doubles the last variable: (x_{i1}...x_{in-1})↑ = x_{i1}...x_{in-2} x_{in-1}^2. This is the KEY technical device.
  - Theorem 3.10: W_X = sum_{sigma in S_V(X,Xbar)} (-1)^{phi(sigma)} p_{Type(sigma)} where Type is now a SET PARTITION (not integer partition). The noncommutative p_pi tracks which vertices are in which cycle.
  - Corollary 3.12: For tournaments, W_X = sum_{sigma in S_V(X), all odd cycles} 2^{psi(sigma)} p_{Type(sigma)}. This is the NONCOMMUTATIVE OCF — same sum but tracking set-partition cycle types.
  - Key insight: The ↑ operation has no obvious tiling interpretation, but it corresponds to "merging the last two vertices while remembering which was which." This is precisely what our contraction approach (T017) failed to handle in commutative setting.
**Next step:** (1) Instead of naive DC -> OCF induction, investigate whether DC can be used to relate the SYMMETRIC FUNCTION U_T across tournaments (e.g., U_T = U_{T'} + correction for single arc reversal). (2) Check if the noncommutative framework gives a new proof of even-odd split or transfer matrix symmetry. (3) CRITICAL: Investigate whether Theorem 3.10 (noncommutative power-sum expansion over SET partitions) implies transfer matrix symmetry when projected to integer partitions.

### INV-052: Mitrovic-Stojadinovic chromatic↔Rédei-Berge connection (arXiv:2506.08841, Jun 2025)
**Source:** kind-pasteur-2026-03-06-S19 web search; opus-2026-03-06-S9 (FULL PAPER READ)
**Status:** DEEPLY READ — HIGH PRIORITY. Contains multiple directly relevant results.
**What:** Proves the chromatic function of a graph and the Rédei-Berge function of a digraph are "almost identical" at the poset level. For any poset P: X_{inc(P)} = omega(U_P) (Theorem 3.2). This extends to noncommutative versions: Y_{inc(P)} = omega(W_P) (Theorem 5.7). Commutative diagram of four Hopf algebra morphisms (p.8, Remark 3.4).
**KEY RESULTS FROM FULL READ (opus-S9):**
  1. **Theorem 4.1 (Converse of Rédei):** If P is a poset that is NOT a chain, then the number of quasi-linear extensions is EVEN. Proof: chi_{inc(P)}(1) and (-1)^|P| u_P(-1) have the same parity. u_P(1) counts quasi-linear extensions; chi_{inc(P)}(1)=0 unless inc(P) is discrete (P is a chain). THIS GENERALIZES RÉDEI'S THEOREM TO POSETS.
  2. **Theorem 4.8:** For any poset P and partition lambda: #{broken-cycle-free subsets of inc(P) with component sizes lambda} = #{permutations in S_V(D_P-bar) with cycle type lambda}. This is a CYCLE-TYPE-PRESERVING bijection between broken cycles and permutations — potentially a route to bijective OCF proof.
  3. **Corollary 3.3:** chi_{inc(P)}(m) = (-1)^|P| u_P(-m). The chromatic polynomial at m colors = Rédei-Berge polynomial at -m. Our H(T) = u_T(1) = (-1)^n chi_{inc(P)}(-1) for the associated poset.
  4. **Theorem 4.6:** D_P-bar contains a Hamiltonian cycle iff P is irreducible. Combined with Corollary 4.5: U_{P_n} forms an algebraic basis of Sym for any sequence of irreducible posets.
  5. **Section 6 (Positivity):** Conjecture 6.3: If P is (3+1)-free, U_P is h-positive. Theorem 6.4: Already proved s-positive. Theorem 6.11: h-positivity propagates through deletion-contraction if there's a "sink-like" vertex.
  6. **Section 7 (Bags of sticks):** Decomposition into bags of sticks (disjoint unions of directed paths) gives explicit formulas. Triple deletion property generalized.
**Why it matters for us:**
  - The bridge X_{inc(P)} = omega(U_P) means we can use 30 years of chromatic symmetric function theory on tournament problems
  - Theorem 4.1 is a clean generalization of the parity theorem we're studying
  - Theorem 4.8 gives a TYPE-PRESERVING bijection — this is exactly what a bijective OCF proof would need
  - The h-positivity results (Section 6) may apply to tournament-associated posets
**Next step:** (1) For a tournament T, identify the associated poset P such that D_P = T. (2) Check if Theorem 4.8 gives a new proof of OCF when specialized to tournaments. (3) Investigate whether tournament posets are (3+1)-free (this would give h-positivity of U_T). (4) Test computationally: does the broken-cycle bijection from Thm 4.8 match our OCF terms?

### INV-053: Savchenko cycle counting formulas for regular tournaments — VERIFIED AT n=7
**Source:** kind-pasteur-2026-03-06-S19 web search; Savchenko J. Graph Theory 83 (2016), Discrete Math (2017), arXiv:2403.07629 (2024)
**Status:** VERIFIED at n=7. Cycle counts are class invariants. DRT vs LTT classification matches.
**What:** Savchenko has a series of papers giving EXACT polynomial formulas for c_k(T) (number of k-cycles) in regular tournaments:
- c5, c6 formulas (2016, J. Graph Theory 83)
- c7 formula (2017, Discrete Math)
- c8 for DRTs vs locally transitive tournaments (2024, arXiv:2403.07629)
Key finding: c8(DRT_n) is INDEPENDENT of which DRT is chosen. Phase transition at n=39: DRTs have more 8-cycles than locally-transitive for n≤35 but FEWER for n≥39.
**Why it matters:** These exact formulas could determine whether Paley tournaments maximize cycle counts at EVERY length, or only for short cycles. The phase transition at n=39 suggests our cycle-maximization mechanism may reverse at larger n. Also, the spectral methods used (eigenvalue-based cycle counting) could connect to our transfer matrix work.
**VERIFIED (kind-pasteur S19):** At n=7, the three regular tournament classes are EXACTLY:
  - DRT (Paley): 240 tours, dc={3:14, 5:42, 7:24}, H=189
  - Locally Transitive: 720 tours, dc={3:14, 5:28, 7:17}, H=175
  - Other Regular: 1680 tours, dc={3:14, 5:36, 7:15}, H=171
Cycle counts are CLASS INVARIANTS (exactly one vector per class). DRT maximizes directed 5-cycles and 7-cycles. LTT has "diametrically opposite" properties per Savchenko.
**EXTENDED (kind-pasteur S21):** Savchenko (2024) proves c_m(DR_n) > c_m(RLT_n) for ALL m = 1,2,3 mod 4 (including all odd m). Only m = 0 mod 4 has the phase transition. This directly explains DRT's H-maximization via OCF.
**DRT n=11 ANALYSIS (kind-pasteur S21):** Two DRT classes at n=11 (from (11,5,2)-difference sets): Paley QR={1,3,4,5,9} (H=95095, c3=55, c5=594, |Aut|=55) and non-Paley {1,2,3,5,8} (H=69311, c3=44, c5=407, |Aut|=11). Cycle counts are NOT invariant across different DRT classes — Paley strictly dominates. Savchenko's c_k invariance is within a single isomorphism class.
**Next step:** (1) Obtain Savchenko's exact polynomial c_k formulas. (2) Test at n=19 or n=23 (multiple DRT classes). (3) Prove Paley maximizes H among ALL DRTs.

### INV-054: Komarov-Mackey exact 5-cycle formula (arXiv:1410.6828, JGT 2017) — PARTIALLY INVESTIGATED
**Source:** kind-pasteur-2026-03-06-S19 web search
**Status:** PARTIALLY INVESTIGATED (kind-pasteur-S39b). c5 is NOT score-determined.
**What:** Exact formula for number of directed 5-cycles in any tournament in terms of edge score sequence. Maximum c5 ≈ (3/4)*C(n,5), achieved by almost all random tournaments. Lower bounds also proved.
**NEW FINDING (S39b):** c5 is NOT determined by score sequence, even combined with sum_d², edge_score, or common_out_neighbor statistics. At n=5, score (1,2,2,2,3) has c5 in {1,2,3}; at n=6, 9/22 score sequences have varying c5. The Komarov-Mackey formula likely involves CUBIC or higher-order graph statistics (e.g., directed walks of length 3+). For regular tournaments, c5 IS a class invariant (Savchenko, verified n=5,7).
**Why it matters:** This rules out O(n²) c5 computation from scores alone. Cycle enumeration (O(n^5) for 5-cycles, or O(2^n) bitmask DP) remains necessary.
**Next step:** (1) Read the exact Komarov-Mackey formula. (2) Identify which cubic invariant determines c5.

### INV-055: Linial-Morgenstern cycle density conjecture and extremal tournaments
**Source:** kind-pasteur-2026-03-06-S19 web search; arXiv:2011.14142 (Ma-Tang), arXiv:1902.00572
**Status:** NEW LEAD — MEDIUM PRIORITY
**What:** Linial-Morgenstern conjecture: among tournaments with fixed c3 density d, the c4 density is minimized by random blowups of transitive tournaments. Proved for d ≥ 1/36 using spectral methods. Ma-Tang extend to c_ℓ for ℓ ≢ 2 mod 4 when d is near 1.
**Why it matters:** This is the "dual" to our maximization question. We show Paley maximizes total directed cycles; this literature characterizes minimizers. The spectral methods used here (eigenvalue-based cycle density bounds) could provide tools for our Paley maximizer proof.
**Next step:** Check if the extremal results constrain H(T) via OCF.

### INV-056: Jerrum-Patel zero-free regions for H-free graphs (JLMS 2026)
**Source:** kind-pasteur-2026-03-06-S19 web search
**Status:** NEW LEAD — MEDIUM PRIORITY (for real-rootedness question)
**What:** Jerrum & Patel (2026, JLMS) prove zero-free regions for the independence polynomial of H-free graphs for various H. For claw-free: all zeros on negative real line (= Chudnovsky-Seymour). For subdivided claws: related zero-free regions. KEY: for H NOT a subdivided claw or path, there exist H-free graphs of max degree 3 with zeros NOT on the negative real line.
**Why it matters:** Our Omega_3(T) has all real roots for n≤20 but is NOT always claw-free (fails n≥9). Jerrum-Patel's results on subdivided claw avoidance may explain why real roots persist beyond n=8. The tournament-specific constraint on Omega_3 structure may ensure avoidance of exactly the "bad" subgraphs.
**Next step:** (1) Check what specific subdivided claws appear in Omega_3(T) at n≥9. (2) Apply Jerrum-Patel to determine if their zero-free regions explain our observations.

### INV-057: Herman's Terwilliger algebras of DRTs (arXiv:2404.11560, 2024)
**Source:** kind-pasteur-2026-03-06-S19 web search
**Status:** NEW LEAD — LOW-MEDIUM PRIORITY
**What:** Allen Herman computes Terwilliger algebras for DRTs (asymmetric rank-3 association schemes). Thin irreducible modules, dimension 4k+9. Key: Terwilliger algebras distinguish non-isomorphic DRTs up to n=23, but FAIL at n=27 (need rational Terwilliger algebras). There are 237 non-isomorphic DRTs at n=27.
**Why it matters:** (1) If all DRTs at small n have the same H(T), that would be a DRT invariant. (2) If Terwilliger algebra structure constrains H(T), this gives an algebraic route to Paley maximization. (3) The n=27 DRT classification gives test cases for our conjectures beyond Paley primes.
**Next step:** (1) Check if all DRTs at n=7 (there's only one: Paley) or n=11 have the same H. (2) At n=27, compare H across different DRT isomorphism classes.

### INV-058: Pantangi critical groups distinguish Paley from other DRTs
**Source:** kind-pasteur-2026-03-06-S19 web search; Pantangi arXiv:1905.08568 (2019)
**Status:** CONNECTION IDENTIFIED
**What:** Pantangi shows critical groups (sandpile groups) distinguish Paley from non-Paley DRTs. Chandler-Sin-Xiang computed Smith/critical groups of Paley GRAPHS. Different DRT constructions (Szekeres-Whiteman 2-block, Wallis-Whiteman 4-block) are distinguished by their critical groups.
**Why it matters:** If H(T) is a DRT invariant AND different DRTs have different critical groups, then H could be read off the critical group. This would give a purely algebraic characterization of the H-maximizer.
**Next step:** Compute critical groups for DRTs at n=11,19 and check correlation with H values.

### INV-062: Universal Master Polynomial and Central Factorial Numbers — VERIFIED (THM-059)
**Source:** opus-2026-03-06-S30
**Status:** VERIFIED computationally (22/22 cases, n=4..9). Algebraic proof pending.
**What:** The per-invariant r-polynomial C_I(r,n) = 2^{parts(I)} * F_f(r) where f = free position count and F_j is determined by the central factorial number triangle (OEIS A036969) via b_{k,j} = b_{k-1,j-1} + j^2 * b_{k-1,j}. This completely determines the entire W-coefficient hierarchy for all tournaments at all n. The shift principle is a corollary.
**Key findings:** (1) F_j(1/2) = 1 for all j. (2) Leading coefficient = (j+1)!. (3) Predictions made for n=11 without computation. (4) Complete n=8 even-n table computed.
**Next step:** (1) Algebraic proof of the central factorial recurrence from position pattern analysis. (2) Verify F_8 prediction at n=11 computationally. (3) Investigate C_0(r) (constant/background polynomial).
**Scripts:** `04-computation/universal_master_polynomial.py`, `04-computation/w1_n8_complete.py`

### INV-059: Cyclic subsets of tournaments (arXiv:2508.03634, Aug 2025)
**Source:** kind-pasteur-2026-03-06-S19 web search; Hunter-Liu-Milojević-Sudakov
**Status:** NEW LEAD — LOW PRIORITY
**What:** Optimal bounds on probability that random induced subtournament of high min-degree tournament is Hamiltonian. Extends to p-biased measure. Proves that high min-degree forces high Hamiltonicity probability.
**Why it matters:** Paley T_p has min-degree (p-1)/2 (doubly regular). This paper could give explicit bounds on the fraction of induced subtournaments that are Hamiltonian, which connects to our cycle counting.
**Next step:** Apply their bounds to Paley tournaments. Check if this gives lower bounds on c_k counts.

### INV-060: Eulerian cycle trace formula (arXiv:2502.02915, Feb 2025) — CLOSED
**Source:** kind-pasteur-2026-03-06-S19 web search; Ye Luo
**Status:** CLOSED (too remote, kind-pasteur-S22 agent investigation)
**What:** Trace formula counting Eulerian cycles via "twisted" vertex and edge adjacency matrices. Uses homological spectral graph theory.
**Finding:** Connection too remote to be actionable. Luo counts Eulerian cycles (edge traversals) via H_1 characters; we count Hamiltonian paths (vertex traversals) via inclusion-exclusion. Different domain, different algebraic structure. Closest classical analogue to our transfer matrix is Ryser's permanent formula, not Luo's trace formula.
**Why it matters:** Our transfer matrix tr(M) = H(T) (THM-027) is also a trace formula. This paper's approach—using twisted adjacency matrices with spectral antisymmetry—could provide a template for proving our trace formula properties (symmetry, off-diagonal sum) at general n.
**Next step:** Read the paper and check if "twisted adjacency" techniques apply to tournament transfer matrices.

### INV-061: Hamilton transversals in tournaments (Combinatorica 2024)
**Source:** kind-pasteur-2026-03-06-S19 web search; Chakraborti-Kim-Lee-Seo arXiv:2307.00912
**Status:** NEW LEAD — LOW PRIORITY
**What:** For collections of sufficiently many tournaments on the same vertex set, transversal Hamilton paths/cycles exist. For m=|V|-1 tournaments, there's a transversal Ham path; for m=|V| with m-1 strongly connected, transversal Ham cycle.
**Why it matters:** The "transversal" perspective could give a new way to relate Ham paths across different tournaments, potentially connecting to how H(T) changes under arc reversals.

### INV-062: Forward arc maximization in tournaments (arXiv:2602.10713, Feb 2026)
**Source:** kind-pasteur-2026-03-06-S19 web search; Guo-Gutin-Lan-Shao-Yeo-Zhou
**Status:** NEW LEAD — LOW PRIORITY
**What:** Characterizes maximum forward arcs in Hamilton cycles/paths for semicomplete and locally semicomplete digraphs. Polynomial-time algorithms.
**Why it matters:** Forward arcs in Hamilton paths relate to our "position-based" analysis (pos(a,P) in THM-027 trace formula). The maximum forward arc structure could inform transfer matrix properties.

### INV-063: Spectral pseudorandomness and Paley clique bounds (Exp. Math. 2024)
**Source:** kind-pasteur-2026-03-06-S19 web search; Kunisky-Yu arXiv:2303.16475
**Status:** NEW LEAD — LOW PRIORITY
**What:** Studies spectral pseudorandomness of Paley graphs via subgraph eigenvalue distributions. Conjecturally, minimum eigenvalue convergence would improve clique number bounds beyond √p.
**Why it matters:** Spectral properties of Paley graphs/tournaments are central to our theory. If Paley tournaments have stronger spectral pseudorandomness than other DRTs, this could explain H-maximization via eigenvalue-based cycle counting formulas.

### INV-064: Mitrovic Hopf algebra new bases (arXiv:2407.18608v3, Mar 2025)
**Source:** kind-pasteur-2026-03-06-S19 web search
**Status:** CONNECTION IDENTIFIED — supplements INV-033
**What:** Introduces two new combinatorial Hopf algebras of posets and permutations with Rédei-Berge functions. Constructs new bases for symmetric functions whose generators are Rédei-Berge functions. Investigates which digraph invariants are extractable from the Rédei-Berge function.
**Why it matters:** If H(T) can be expressed as a coefficient in one of these new bases, it gives an algebraic handle on Hamiltonian path counting.
**Next step:** Check which digraph invariants the paper extracts. Is H(T) among them?

### INV-065: Independence polynomial root gap (arXiv:2510.09197, FSTTCS 2025)
**Source:** kind-pasteur-2026-03-06-S19 web search; Om Prakash & Vikram Sharma
**Status:** NEW LEAD — LOW PRIORITY
**What:** Quantifies the gap between the smallest real root β(G) of I(G,x) and all other roots. For connected graphs, β(G) is a simple real root smaller than 1, but previous proofs gave no gap bound. This paper provides explicit bounds.
**Why it matters:** For our Omega(T) real-rootedness question, having a gap bound could help prove that all roots are real by showing they're well-separated from the complex plane.

### INV-066: Low-rank matrices from tournaments and symmetric designs (arXiv:2401.14015, 2024)
**Source:** kind-pasteur-2026-03-06-S19 web search; Balachandran-Sankarnarayanan
**Status:** NEW LEAD — LOW-MEDIUM PRIORITY
**What:** Constructs symmetric matrices from tournament structures where rank depends on design-theoretic properties. Symmetric designs (BIBDs) give matrices with rank near n/2. The rank-topology relationship involves bipartite graph eigenvalues.
**Why it matters:** Our transfer matrix M is constructed from a tournament and is symmetric. This paper's framework connecting tournament-derived matrices with design theory could explain structural properties of M (e.g., why symmetry holds, what the rank structure is).
**Next step:** Check if our M fits their M_T(f,a) framework.

### INV-067: Alpha_1 gap theorem and converse of Redei — CORRECTED S22
**Source:** kind-pasteur-2026-03-06-S21 (computation), CORRECTED S22
**Status:** PARTIALLY PROVED (THM-029 corrected)
**What:** alpha_1=3 is impossible at n<=6 but ACHIEVABLE at n>=7 (~9.2% of c3=3 at n=7). Common-vertex property fails at n>=7. HOWEVER, H=7 remains impossible for ALL n by refined argument: H=7 requires (alpha_1=3, i_2=0), but i_2=0 forces common vertex => c5>=1 => alpha_1>=4; while alpha_1=3 implies i_2>=1 => H>=11. H=21 absent through n=7.
**Achievable H values:** n=5: {1,3,5,9,11,13,15}. n=6: {1,3,5,9,11,...,45}\{7,21,35,39}. At n=7: 35 and 39 become achievable but 7 and 21 remain gaps.
**Connection:** Relates to Mitrovic-Stojadinovic "converse of Redei" (INV-052).
**Next step:** (1) Prove alpha_1=3 impossibility for ALL n (not just n<=6). (2) Find all permanent H-gaps. (3) Check OEIS for the sequence of achievable alpha_1 values.
**Scripts:** h7_impossibility.py, alpha1_gaps.py, alpha1_gap3_proof.py, c3_forces_c5.py, redei_converse_fast.py
**Theorem:** THM-029

### INV-068: DRT non-uniqueness and Paley dominance at n=11 — NEW FINDING
**Source:** kind-pasteur-2026-03-06-S21 (computation)
**Status:** COMPUTED
**What:** At n=11, exactly 2 DRT isomorphism classes exist (from (11,5,2)-difference sets). Paley (QR={1,3,4,5,9}) has H=95095, c3=55, c5=594, |Aut|=55. Non-Paley ({1,2,3,5,8}) has H=69311, c3=44, c5=407, |Aut|=11. Both are self-converse. Paley strictly dominates in ALL cycle counts.
**Key insight:** "DRT" is NOT a single tournament — different DRTs can have different cycle counts and H values. Paley is the BEST DRT. This strengthens the Paley maximizer conjecture from "Paley maximizes among all tournaments" to "Paley maximizes even among DRTs."
**Next step:** Test at n=19 and n=23 where more DRT classes exist.
**Scripts:** drt_n11_analysis.py, drt_n11_verify.py, drt_n11_deeper.py

### INV-069: Scalar M characterization — M=(H/n)*I ↔ H-maximizer at n=5, ↔ VT at n=7
**Source:** kind-pasteur-2026-03-06-S25c (T156), opus-2026-03-06-S26 (T158)
**Status:** COMPUTED, OPEN CONJECTURE
**What:** Transfer matrix M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b) is scalar (= (H/n)*I) for certain tournaments. At n=5: ALL 64 H-maximizers (H=15) have M=3*I, including 40 non-regular (scores 1,2,2,2,3) with |Aut|=3 (Z/3Z). These non-VT tournaments still have uniform position distribution. At n=7: M scalar ⟺ VT (100 sampled regular, 0 counterexamples). All circulant tournaments have scalar M. The characterization seems to be "uniform position distribution" (every vertex appears equally often at each position in Ham paths), which is weaker than VT at n=5 but equivalent at n=7.
**Key results:**
- n=5: M scalar ⟹ H=15 (max), but NOT ⟹ VT or regular
- n=7: M scalar ⟺ VT (conjecture, verified 100 samples)
- All circulant ⟹ M scalar (verified n=3,5,7)
- Paley T_7: |Aut|=21 (Frobenius Z/7Z⋊Z/3Z), sigma(x)=2x mod 7 is extra aut
- opus finding: non-regular n=5 scalar M has Z/3Z aut, uniform endpoint counts, uniform 3-cycle participation
**NONHAM VANISHING (kind-pasteur S25c):** The split pair decomposition M[a,b] = HAM(a,b) + NONHAM(a,b) shows:
- NONHAM(a,b) = 0 trivially when T[a,b]=1 (all split pairs are Ham paths)
- NONHAM(a,b) = M[a,b] when T[a,b]=0 (junction always fails)
- So NONHAM=0 for all (a,b) ⟺ M[a,b]=0 whenever T[a,b]=0
- VERIFIED: NONHAM=0 for ALL position-uniform n=5 (64/64 exhaustive) and ALL circulant n=7 (8/8)
- NONHAM ≠ 0 for general (non-uniform) tournaments at n=3,4,5
**PROOF CHAIN (verified n=3,5,7):**
1. Position-uniform ⟹ NONHAM=0 ⟹ M[a,b]=0 when T[a,b]=0
2. THM-030: M[a,b]=M[b,a]
3. For T[a,b]=1: T[b,a]=0 ⟹ M[b,a]=0 ⟹ M[a,b]=0 by symmetry
4. M[a,a]=H/n from position uniformity
5. M=(H/n)*I. QED (modulo proving step 1 at general n)
**Algebraic proof of step 1:** OPEN. The cancellation mechanism: for uniform T, nonzero E*B products at adjacent subset sizes pair up and cancel. For non-uniform T, "orphan" terms remain.
**Connection to Björklund et al.:** M[a,b] is a subset convolution in their "Fourier meets Möbius" framework.
**Next step:** (1) Prove NONHAM=0 for position-uniform tournaments at general n. (2) Test at n=9 (circulant only, n=9 exhaustive infeasible). (3) Find algebraic proof of the pairing mechanism.
**Scripts:** `04-computation/transfer_matrix_maximizers.py`, `04-computation/circulant_scalar_m_proof.py`, `04-computation/scalar_m_characterization.py`, `04-computation/scalar_m_n5_analysis.py`, `04-computation/circulant_m_scalar_proof.py`, `04-computation/nonham_vanish_general.py`, `04-computation/nonham_vanish_uniform.py`, `04-computation/nonham_proof_analysis.py`

### INV-070: Fibonacci determinant formula for transitive tournament
**Source:** kind-pasteur-2026-03-06-S25c (T157)
**Status:** VERIFIED n=2,...,11
**What:** For the transitive tournament T_n (i beats j iff i<j), det(M) = (-1)^{n(n-1)/2} * F(n+1) where F is Fibonacci. The matrix D*M = I+U-L (tridiagonal: 1s on diagonal/superdiagonal, -1s on subdiagonal) satisfies the Fibonacci recurrence. The Chebyshev eigenvalue conjecture (eigenvalues = 2cos(kπ/(n+1))) is FALSE.
**Connection:** The transitive tournament is the unique acyclic tournament, so Omega(T_n) is empty and I(Omega,2)=1=H(T_n). The transfer matrix has a clean tridiagonal structure reflecting the total order.
**Next step:** (1) Does the Fibonacci structure extend to near-transitive tournaments? (2) What is det(M) for other named tournament families?
**Scripts:** `04-computation/fibonacci_determinant_proof.py`

### INV-071: det(M) as tournament invariant — exhaustive n=5
**Source:** kind-pasteur-2026-03-06-S25c (T157)
**Status:** COMPUTED
**What:** At n=5, det(M) takes values {-27,-9,1,8,9,16,17,243} across 10 isomorphism classes (up to converse). The 10 eigenvalue patterns of M classify tournaments into exactly the 10 isomorphism classes. det(M) is NOT determined by (H,c3) — e.g., H=5,c3=2 gives both det=9 and det=17.
**Key finding:** det(M) is not related to det(A), per(A), or any simple matrix expression of the adjacency matrix A. Tested det(I±S), det(A+A^T), det(A*A^T), det(I+A*A^T) — none match. The transfer matrix is fundamentally different from the adjacency matrix.
**opus finding:** char poly det(λI-M) encodes hierarchy of path correlations. Scalar M ⟹ char poly = (λ-H/n)^n. PSD threshold at n=5: H≥13 ⟹ M is PSD.
**Next step:** (1) Compute det(M) exhaustively at n=7 (sample). (2) Find if any spectral graph invariant matches. (3) Investigate PSD threshold at general odd n.
**Scripts:** `04-computation/det_m_general_formula.py`, `04-computation/det_m_vs_adjacency_spectrum.py`, `04-computation/char_poly_M_analysis.py`, `04-computation/pos_skeleton_connection.py`

### INV-072: IO walk GF vs transfer matrix bridge
**Source:** opus-2026-03-06-S26, kind-pasteur-2026-03-06-S25c
**Status:** STRUCTURAL COMPARISON, NO BRIDGE FOUND
**What:** Irving-Omar walk GF W_D(z) = det(I+zXA^T)/det(I-zXA) uses cycle covers (det/per), while M[a,b] uses path decomposition (E_a*B_b). M[a,b] ≠ H(a→b) (direct endpoint-conditioned count). W(-z,-r) = W(z,r) at commutative level (opus finding). No simple matrix expression of A gives M.
**Key structural difference:** IO is multiplicative (det/per = products over cycles), M is additive (sum over subsets with inclusion-exclusion signs). Bridge might exist through Hopf algebra coproduct structure.
**Next step:** (1) Express M[a,b] in terms of IO's det/per framework. (2) Check if deletion-contraction on M matches Mitrovic's W_X = W_{X\e} - W_{X/e}^up.

### INV-073: Palindromic N => Scalar M for circulant tournaments — PROVED (THM-052)
**Source:** kind-pasteur-2026-03-06-S25d, building on THM-050 (opus-S26), THM-051 (opus-S26)
**Status:** PROVED for all circulant tournaments at odd n.
**What:** For circulant T on Z/nZ, the consecutive-position count N(a,b,j) = f(b-a mod n, j) is palindromic: f(d,j) = f(d,n-2-j). The proof uses three ingredients: (1) translation symmetry gives f depends only on d, (2) N symmetry gives f(d)=f(n-d), (3) self-complementarity via sigma: i->-i gives f(d,j)=f(n-d,n-2-j). Combining (2)+(3): f(d,j)=f(d,n-2-j). At odd n, palindromic N forces alternating sum = 0, so M[a,b]=0 for a!=b. Combined with M[a,a]=H/n, gives M=(H/n)*I.
**Verification:** n=5 (64/64 exhaustive pos-uniform), n=7 (8/8 circulant), n=9 (16/16), n=11 Paley, n=13 Paley.
**Key finding:** ALL position-uniform n=5 tournaments are self-complementary (64/64). No position-uniform tournaments exist at even n (n=4,6).
**Open extension:** Prove for general vertex-transitive (non-circulant) tournaments. The proof uses circulant-specific translation symmetry. At n=15, non-circulant VT tournaments exist (Babai-Kantor doubly-regular tournaments).
**Scripts:** `04-computation/palindromic_N_proof.py`, `04-computation/palindromic_N_posuniform.py`, `04-computation/palindromic_N_n9.py`, `04-computation/palindromic_N_n11.py`, `04-computation/selfcomp_posuniform_n7.py`

### INV-074: Diagonal signed position theorem — VERIFIED n=5
**Source:** opus-2026-03-06-S11b (continued³)
**Status:** VERIFIED computationally at n=5 (all 12 iso classes).
**What:** M[v,v] = sum_P (-1)^{pos(v,P)} where the sum is over all Hamiltonian paths P of T and pos(v,P) is the 0-indexed position of v in P. This means M[v,v] counts even-position appearances minus odd-position appearances. M[v,v] can be NEGATIVE (not a path count). For VT tournaments at odd n: M[v,v] = H/n. "Defect vertex" = vertex whose position distribution is biased relative to average.
**Connection to THM-027:** The trace formula tr(M) = H(T) = sum_v sum_P (-1)^{pos(v,P)} = sum_P sum_v (-1)^{pos(v,P)} = sum_P 1 (at odd n, since alternating sum of 0..n-1 positions = 1). This reproves the trace formula.
**Next step:** Prove from IE formula definition. Connect to defect vertex characterization.
**Scripts:** `04-computation/diagonal_signed_position_theorem.py`

### INV-075: Perpendicularity of M-directions across H-classes — CONFIRMED n=7
**Source:** opus-2026-03-06-S11b (continued³)
**Status:** CONFIRMED computationally at n=7 (790 iso classes).
**What:** The non-scalar part of M (i.e., M - (H/n)*I) has a "direction" in matrix space. Measuring cosine similarity between these directions across iso classes shows:
  - Low H (H<85): positive cosine (~0.5-0.8, aligned)
  - Mid H (H≈85-105): near zero (perpendicular!)
  - High H (H>105): negative cosine (anti-aligned)
  - Overall mean cosine = -0.0485 (near perpendicular)
The crossover (true perpendicularity) occurs near the MEDIAN H value. This means the non-scalar perturbation "rotates" continuously through eigenspace as H varies.
**Connection:** This is the "perpendicularity" the user hypothesized in earlier sessions. It connects to the inverted-U of position variance and the grid symmetry structure.
**Next step:** (1) Prove analytically why crossover occurs at median H. (2) Check at n=9. (3) Connect to the even cycle vanishing theorem (INV-053) which uses the same T↔T^op involution.
**Scripts:** `04-computation/perpendicularity_cosine_n7.py`

### INV-076: All H-maximizers have scalar M — VERIFIED n=5,7
**Source:** opus-2026-03-06-S11b (continued³)
**Status:** VERIFIED exhaustively at n=5 (both max-H classes) and n=7 (all 43 maximizers).
**What:** Every tournament achieving max H has M = (H/n)*I (scalar transfer matrix). At n=7: all 43 maximizers have M = 27*I. At n=5: both H=15 classes (VT circulant and non-regular VT) have M = 3*I.
**Conjecture:** For ALL odd n, the H-maximizer has scalar M. This is equivalent to saying H-maximizers are always vertex-transitive (or at least "position-uniform").
**Connection:** Combines with THM-052 (circulant => scalar) and the Paley maximizer conjecture. If Paley maximizes H and Paley is circulant, then Paley gives scalar M. The deeper question is whether scalar M is NECESSARY for H-maximization.
**Next step:** Verify at n=9 (need to find maximizers first).

### INV-077: VT tournament NOT self-converse at n=21 — THM-052 DISPROVED for non-SC VT
**Source:** kind-pasteur-2026-03-06-S25e
**Status:** RESOLVED. M is NOT scalar for non-SC VT tournaments.
**What:** ALL 22 non-circulant VT tournaments at n=21 (from McKay's database) are NOT self-converse. These are Cayley tournaments on F_21 = Z/7 x| Z/3. All 88 circulant VT tournaments at n=21 ARE self-converse.
**Computation (n=21, 1075s):**
- H(T) = 123,522,430,238,361 (divisible by 21)
- N(0,1,j) is NOT palindromic: N[0]=581,223,220,317 vs N[19]=581,314,958,778
- Alternating sum = M[0,1] = 45,478,409 != 0
- Therefore M != (H/n)*I for this VT tournament
**Conclusion:** THM-052 is PROVED for self-converse VT (including all circulants) but DISPROVED for non-SC VT. Self-converse is the exact boundary.
**Scripts:** `04-computation/frobenius21_palindromic_N.py`, `04-computation/mcKay_vt21_selfconverse.py`

### INV-078: Aut(T) union Anti(T) transitivity characterizes scalar M
**Source:** opus-2026-03-06-S26 (scalar_m_aut_anti_characterization.py)
**Status:** VERIFIED at n=5 (exhaustive). CONFIRMED at n=21: F_21 non-SC has Anti=empty, M not scalar.
**What:** Scalar M (M = (H/n)*I) holds iff Aut(T) union Anti(T) acts transitively on V. For the F_21 non-normal tournament: Aut=F_21 (transitive) but Anti=empty, so Aut+Anti = Aut alone. The conjecture predicts M not scalar, which is CONFIRMED by computation.
**Next step:** Verify at n=7 (exhaustive).

### INV-079: W(r) coefficient stratification by odd-cycle complexity
**Source:** kind-pasteur-2026-03-06-S25f, opus-2026-03-06-S27 (THM-055)
**Status:** PROVED (k=0,1), VERIFIED (k=2). Connected to Hopf algebra coproduct.
**What:** W(r) = sum_P prod(r + s_e) has coefficients w_{n-1-2k} = sum_P e_{2k}(s_P).
  - w_{n-1} = n! (universal)
  - w_{n-3} = 2*(n-2)!*t_3 - const (depends on t_3 only — PROVED)
  - w_{n-5} depends on t_3 AND the 4th moment of f_P (opus THM-055)
  - At n=5: w_0 = -t_3 + 2*t_5 + 1 (EXACT, kind-pasteur verified exhaustive)
**Key identity:** H = 1 + 2*(t_3 + t_5) at n=5 (OCF simplification since a_2=0)
**Recursive Hopf structure:** overlap=3 contribution to w_{n-5} at n=7 uses OCF at n=5 on each 5-element sub-tournament. This is the Hopf algebra coproduct Delta([T]) evaluated on fibers.
**Connection to THM-055:** e_{2k}(s_P) is a polynomial of degree 2k in f_P. All power sums p_{2l} are constant; only p_1 = f - (n-1)/2 varies. So everything reduces to moments of f_P.
**Next step:** (1) Find explicit formula for w_0 at n=7 in terms of tournament invariants. (2) Prove the Hopf algebra recursion algebraically. (3) Determine whether the 4th moment has a cycle-theoretic interpretation.

### INV-080: Pfaffian-path duality at even/odd n
**Source:** kind-pasteur-2026-03-06-S25f (pfaffian_path_duality.py)
**Status:** COMPUTATIONAL. Interesting correlations found.
**What:** At even n: det(S) = Pf(S)^2 is a nonzero odd square. At odd n: det(S) = 0 but tr(M) = H > 0.
  - n=4: det(S) is EXACTLY determined by t_3: det=1 if t_3 even, det=9 if t_3 odd
  - n=6: det(S) is NOT determined by t_3 alone (needs finer invariants)
  - |Pf| always odd (Fisher-Ryan), values in {1,3,5,...,n-1} at even n
  - The Pfaffian and H are NOT functionally related but are both determined by S
**Connections:** D_k classification (Zeng-You-Zhao 2025), Seidel tournament matrices (determinant gap phenomena)
**Speculative:** Is there a formal duality: paths (odd n) ↔ cycle covers (even n)?
**Next step:** (1) Check if det(S) at n=6 is determined by (t_3, t_5). (2) Investigate the eigenvalue connection between Pf and H.

### INV-081: Paley tournament W(r) structure
**Source:** kind-pasteur-2026-03-06-S25f (eigenvalue_W_connection.py)
**Status:** COMPUTED at p=7.
**What:** Paley T_7 has W(r)/7! = [1/320, 0, 1/80, 0, 1/4, 0, 1].
  - Eigenvalues of A_Paley: (p-1)/2 = 3 (mult 1), (-1±sqrt(p))/2 (mult (p-1)/2 each)
  - All non-trivial eigenvalues degenerate → scalar M → W coefficients have maximal symmetry
  - Paley requires p ≡ 3 (mod 4) (so -1 is not a QR)
**Next step:** Compute W(r) for T_11 and T_3. Check if W/p! ratios have a closed form involving p.

### INV-082: EXACT W-coefficient hierarchy — spectral decomposition of tournaments
**Source:** kind-pasteur-S25g
**Status:** VERIFIED at n=5 (exhaustive) and n=7 (20 random samples, 0 error)
**What:** W(r) coefficients form a hierarchy: w_{n-1-2k} depends on cycle data of complexity k.
  - w_{n-1} = n! (universal)
  - w_{n-3} = (n-2)! * [2*t_3 - C(n,3)/2] (depends on t_3 only; CENTERED: zero mean over random T)
  - w_{n-5} at n=7: -60*t_3 + 12*t_5 + 24*alpha_2 + 231
  - w_0 at n=7: 2*t_3 - t_5 + 2*t_7 - 2*alpha_2 - 17/4
  - Each level adds ONE new cycle complexity (t_{2k+1} or alpha_k)
  - H - w_0 penalty SHIFTS to higher-order cycles at larger n
  - n=5: H - w_0 = 3*t_3; n=7: H - w_0 = 3*t_5 + 6*alpha_2 + 21/4
  - Analogous to spectral decomposition / renormalization group / characteristic classes
**Next step:** Verify w_{n-3} = (n-2)! * [2*t_3 - C(n,3)/2] at n=9. Compute full hierarchy at n=9.

### INV-083: Rooted tournaments = OEIS A093934
**Source:** kind-pasteur-S25g
**Status:** VERIFIED n=2 through 6
**What:** P(n) = sum over iso classes of (# vertex orbits) = # rooted tournament iso classes.
  - P(2)=2, P(3)=4, P(4)=12, P(5)=48, P(6)=296
  - Matches OEIS A093934 (with offset)
  - P(n) = 2*(n-1)! for n<=5 (coincidence), FAILS at n=6
  - Orbit distributions: n=5 {1:1,3:4,5:7}, n=6 {2:5,4:10,6:41}
**Next step:** Check A093934 description more carefully. Compute P(7) if feasible.

### INV-084: W-coefficient hierarchy as Mobius inversion on cycle complex
**Source:** kind-pasteur-S25g (creative synthesis)
**Status:** CONJECTURED
**What:** The W-coefficients can be viewed as evaluations of I(Omega, x) at different points.
  - H = I(Omega, 2) (OCF, x=2)
  - chi = I(Omega, -1) (Euler characteristic, x=-1)
  - w_0 is an "intermediate" evaluation (not simply I(Omega, c) for some c)
  - The hierarchy parallels Fourier decomposition, but REVERSED: high "frequencies" are simple
  - Connects to renormalization: each level "integrates out" one cycle scale
**Next step:** Determine if w_0 = I(Omega, c) for some specific c, or if it's genuinely different.

### INV-085: Bipartite skeleton and t3 parity (THM-060)
**Source:** kind-pasteur-2026-03-06-S25h
**Status:** PROVED at n=3,5,7,9. Structural argument + exhaustive verification.
**What:** Blue line skeleton (GS flip graph on SC classes) is bipartite at odd n, with bipartition determined by t3 parity. At even n, skeleton is NOT bipartite (has 3-cycles).
**Key mechanism:** Consecutive triples each contribute 1 to t3(T)+t3(flip(T)), total n-2 (odd at odd n). Non-consecutive triples contribute even total for GS tilings.
**Open:** Algebraic proof of Type B evenness; spectral structure at large n; connection between skeleton eigenvalues and tournament invariants.
**Scripts:** 04-computation/t3_parity_proof_complete.py, bipartition_invariant.py, bipartition_n7_verify.py
**Writeup:** 01-canon/theorems/THM-060-bipartite-skeleton.md, 03-artifacts/drafts/bipartite-skeleton-synthesis-S25h.md

### INV-086: Silver ratio in skeleton eigenvalues
**Source:** kind-pasteur-2026-03-06-S25h, skeleton_spectral.py
**Status:** OBSERVED at n=5. Eigenvalues {±(1+√2), ±1, ±1, ±(√2-1)}.
**What:** The skeleton adjacency matrix at n=5 has eigenvalues involving the silver ratio 1+√2. K^2 diagonal = GS class sizes. Is this coincidence or does it persist at n=7?
**Next step:** Compute skeleton eigenvalues at n=7 (88×88 matrix). Check if silver ratio generalizes.

### INV-087: Antiferromagnetic interpretation of skeleton
**Source:** kind-pasteur-2026-03-06-S25h
**Status:** CONCEPTUAL. Skeleton = Ising model with antiferromagnetic coupling.
**What:** SC classes have "spin" = (-1)^{t3}. GS flip edges connect opposite spins. At odd n: perfect Neel order (unfrustrated). At even n: frustrated (odd cycles). Connects tournament theory to statistical mechanics.
**Next step:** Compute partition function Z(beta) = sum over SC of H(T)^beta. Check for phase transitions.

### INV-088: Schweser-Stiebitz-Toft — Rédei revisited (expository)
**Source:** arXiv:2510.10659 (Oct 2025, revised Feb 2026), found by web search opus-2026-03-07-S36
**Status:** CATALOGED. Expository paper, likely low priority.
**What:** Revisits the classical theorems of Rédei, Dirac, and Berge on Hamiltonian paths in tournaments. Exhibits the stronger theorems and explains connections between them. Does NOT mention independence polynomials, odd cycles, or conflict graphs.
**Next step:** Skim for any novel structural insight about H(T) parity not already in our framework. Low priority.

### INV-089: Irving-Omar authorship correction
**Source:** opus-2026-03-07-S37 (this session)
**Status:** CORRECTED in THM-002, CONJ-001, THM-070.
**What:** arXiv:2412.10572 ("Revisiting The Rédei-Berge Symmetric Functions via Matrix Algebra") is by **Irving & Omar**, NOT Grinberg & Stanley. Their Corollary 20 restates Grinberg-Stanley's Theorem 1.39 + Lemma 6.5 from arXiv:2307.05569. The OCF result itself is correctly attributed to Grinberg-Stanley; only the paper authorship was wrong.
**Remaining:** Some computation scripts and broadcast messages still reference "Grinberg-Stanley" for arXiv:2412.10572. These are historical and low priority to fix.

### INV-090: Three equivalent H(T) formulas and the even-cycle cancellation
**Source:** opus-2026-03-07-S37
**Status:** VERIFIED computationally. Structural explanation needed.
**What:** Three equivalent ways to compute H(T):
1. **Direct**: count Hamiltonian paths (Held-Karp)
2. **OCF = I(Ω(T), 2)**: sum over T-only disjoint odd cycle collections, weight 2^k
3. **ps_1(U_T)(1)**: sum of ALL Rédei-Berge p-coefficients (uses T + T^op cycles, signed by (-1)^φ)

KEY FINDINGS at n=7 (Paley):
- All p-coefficients of U_T have ALL-ODD partitions (no even-part lambdas appear!)
- Sum of coefficients = 189 = H(T) [ps_1 at m=1]
- OCF specialization of FULL U_T = 433 ≠ H(T)
- OCF specialization of T-ONLY part = 189 = H(T)
- The "OCF specialization" (p_1→1, p_odd→2, p_even→0) is NOT how GS proves OCF
- GS uses ps_1(U_T)(1) which uses ALL cycles, not just T-direction ones

OPEN: Why are the p-coefficients supported only on all-odd partitions at n=7? Is this true for all tournaments? Omega constraint: (-1)^{n-l(λ)} symmetry forces sum over even-l terms = 0, but here they're individually zero. Is the all-odd support a coincidence of n=7 (Mersenne prime) or universal?

**Mixed-direction findings at n=9:** Mixed T/T^op cycle pairs DO exist at n=9 (100+ per tournament). So U_T at n=9 WILL have even-part lambda contributions, but they must cancel in ps_1(1).
**RESOLVED (opus-S38 agent):** The all-odd support IS universal for tournaments (proved by Grinberg-Stanley Theorem 1.39). The factor 2^k in OCF is the **orientation multiplicity**: each odd cycle has exactly 2 directed orientations (T and T^op), both contributing sign +1 (since (-1)^{k-1}=+1 for odd k). An independent set of k cycles thus contributes 2^k copies. Even-part lambdas vanish because even cycles contribute opposite signs from T and T^op directions.
**Scripts:** `04-computation/omega_constraint.py`, `mixed_sum_n9.py`, `ut_specialization_n9.py`

### INV-091: H=21 permanent gap — PROVED through n=7
**Source:** opus-2026-03-07-S38, kind-pasteur-2026-03-07-S28
**Status:** PROVED through n=7 (exhaustive). Strong conjecture for all n.
**What:** H(T)=21 is never achieved by any tournament on n<=7 vertices. Exhaustive computation: 2,097,152 tournaments at n=7, zero with H=21. The gap 19→23 appears at both n=6 and n=7. OCF analysis: none of the valid (alpha_1,alpha_2) decompositions for H=21 are achievable.
**Key structural insight:** At n=6, the achievable alpha_1 values jump in a way that skips all decompositions summing to H=21. The constraint is that certain cycle counts force additional cycles or disjoint pairs, pushing H past 21.
**Scripts:** `04-computation/h21_n7_fast.py`, `04-computation/h21_theory_fixed.py`
**Theorem:** THM-075
**Next step:** (1) Prove H=21 impossibility at general n. (2) Characterize ALL permanent gaps. (3) Is {7, 21} the complete list, or are there more?

### INV-092: Type-count sequence = A000009 (partitions into distinct parts)
**Source:** opus-2026-03-07-S38 agent, kind-pasteur-S29
**Status:** PROVED (bijective argument).
**What:** The number of OCF cycle types at size n (multisets of odd parts >=3 summing to <=n, plus empty) equals A000009(n) = number of partitions into distinct parts. Bijection: remove all 1's from partition into odd parts. Null-dim sequence 0,0,1,3,6,11,19,29,44,65 is NOT in OEIS (likely novel).

### INV-093: Tangent number connection proved
**Source:** opus-2026-03-07-S38 agent, kind-pasteur-S29
**Status:** PROVED.
**What:** P_n(0,0) = 2^{(n-1)/2} * T_n where T_n is the n-th tangent number. Proof: P_n(u,0) = A_n(t)/t^{(n-1)/2} with u=t+1/t. Setting u=0 gives t=i, and A_n(i)/i^{(n-1)/2} = 2^{(n-1)/2} * T_n by the EGF of Eulerian polynomials evaluated at t=i.
**Connection:** Hetyei (2017, arXiv:1704.07245) on "alternation acyclic tournaments" connects tournament counts to median Genocchi numbers (same family as tangent numbers).

### INV-094: Mitrovic noncommutative deletion-contraction explored
**Source:** opus-2026-03-07-S38, arXiv:2504.20968 (Mitrovic, Apr 2025)
**Status:** EXPLORED. Not useful for OCF.
**What:** W_X = W_{X\e} - W_{X/e}^up for the noncommutative Redei-Berge function. Edge contraction T/e is NOT a tournament (bidirectional edges possible). OCF fails for T\e and T/e. The deletion-contraction approach is useful for algebraic properties of U_T but not for proving OCF.
**Scripts:** `04-computation/edge_deletion_contraction.py`

### INV-095: Bags-of-sticks for OCF — DEAD END
**Source:** opus-2026-03-07-S38 agent
**Status:** CLOSED (dead end).
**What:** The bags-of-sticks decomposition (Mitrovic-Stojadinovic Theorem 4.2) expresses U_X via inclusion-exclusion over edge deletions. Under OCF specialization, every bag of sticks contributes 1 (acyclic digraphs have empty Omega). So the decomposition reduces to: H(T) = sum of inclusion-exclusion coefficients, which is trivially true. No new information for OCF.

### INV-096: H=21 Component Reduction (THM-079) — PROVED FOR ALL n
**Source:** opus-2026-03-07-S39 (partial), kind-pasteur-2026-03-07-S33 (completion)
**Status:** PROVED. H(T) ≠ 21 for ALL tournaments on ALL n vertices.
**What:** For H(T)=21, the OCF requires I(Omega(T),2)=21. Component factorization gives:
  - Disconnected case: IMPOSSIBLE. 21=3*7, but I(component)=7 impossible by THM-029 argument.
  - P_4 case: IMPOSSIBLE. P_4 realization blocked because sharing 3-cycles force extra cycles.
  - K_6-2e case: SUPERSEDED by Dichotomy Theorem proof.
**Dichotomy Theorem (Part R):** For cycle-rich T on n≥9, either (a) 3 disjoint 3-cycles exist (⟹ H≥27), or (b) safe deletion to cycle-rich T−v exists (⟹ H≥H(T−v)+2≥27). Proved via poisoning graph DAG argument. Combined with base case n≤8 (exhaustive) and Part J (non-cyclic vertex removal), gives H(T)≠21 for ALL n.
**Key ingredients:** Lemma Q (cycle-rich ⟹ no source/sink), poisoning graph has out-degree ≤1 and is acyclic, DAG source deletion preserves cycle-rich.
**Scripts:** `04-computation/h21_gap_mechanism.py`, `h21_dichotomy_proof.py`, `h21_poisoning_graph.py`, `h21_cycle_rich_auto_no_ss.py`
**Writeup:** `03-artifacts/drafts/dichotomy-proof-formal.md`

### INV-097: u_T Size-Weighted Independence Polynomial (THM-078) — PROVED
**Source:** opus-2026-03-07-S39
**Status:** PROVED. u_T(m) = sum_j sw(j)*m^{n-2j} where sw(j) = sum over j-element independent sets of 2^|S|.
**What:** The size-weighted independence polynomial identity connects u_T(m) to the Omega(T) independence structure. Q_T(w) = u_T(sqrt(w))/sqrt(w) has all real non-positive roots for n<=8 (Leake-Ryder/Chudnovsky-Seymour stability for claw-free graphs). Fails at n>=9 when claws appear.
**Next step:** Check if Q_T root structure has implications for achievable H values.

### INV-098: Lichiardopol's Conjecture and Disjoint Cycle Forcing — EXPLORED
**Source:** kind-pasteur-2026-03-07-S33 (web research)
**Status:** EXPLORED. Used in H=21 proof context but not directly needed.
**What:** Lichiardopol's conjecture (proved for q=3 by Bang-Jensen, Bessy, Thomassé): tournaments with min out-degree ≥ (q-1)k-1 have k vertex-disjoint q-cycles. For 3-cycles with k=3: min outdeg ≥ 5. However, at n=9 cycle-rich, min outdeg is ALWAYS ≤ 4 (100% of 106,424 tested), so Lichiardopol doesn't directly fire. The poisoning graph argument covers ALL cases including those below the Lichiardopol threshold.
**Papers:** [Bang-Jensen-Bessy-Thomassé](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v27i2p52)
**Next step:** Could be useful for other permanent H-gap proofs where k≥4 disjoint cycles are needed.

### INV-099: Chen-Chang 2024 Disjoint Cycles in Tournaments — CATALOGED
**Source:** kind-pasteur-2026-03-07-S33 (web research)
**Status:** CATALOGED. Not yet deeply investigated.
**What:** Chen-Chang (2024, J. Graph Theory) prove results on disjoint cycles in tournaments. Extends cycle-matching theory. Could provide tools for proving existence of k≥4 disjoint 3-cycles under weaker conditions than Lichiardopol.
**Paper:** [Chen-Chang 2024](https://onlinelibrary.wiley.com/doi/10.1002/jgt.23038)
**Next step:** Read the paper for potentially stronger theorems applicable to H-gap proofs.

### INV-100: Frankl's Proof of Erdos Matching Conjecture (k=3) — CATALOGED
**Source:** kind-pasteur-2026-03-07-S33 (web research)
**Status:** CATALOGED. Provides context for 3-uniform hypergraph matching.
**What:** Frankl proved the Erdos matching conjecture for k=3: bounds the max number of 3-element sets with no matching of size s+1. The cycle vertex sets in Omega(T) form a 3-uniform hypergraph (for 3-cycles). Frankl's bound could constrain the maximum number of 3-cycles with bounded matching number.
**Connection:** If mm(T) ≤ 2, Frankl's bound limits |Omega_3(T)| ≤ max(C(5,3), 3*3-3+1) depending on exact formulation. This could give an independent route to the dichotomy.
**Next step:** Check exact Frankl bound for our setting (n vertices, 3-uniform, matching ≤ 2).

### INV-101: Other Permanent H-Gaps Beyond 7 and 21 — CONFIRMED THROUGH n=8 EXHAUSTIVE
**Source:** kind-pasteur-2026-03-07-S33, opus-2026-03-07-S43
**Status:** STRONG CONJECTURE that H=7 and H=21 are the ONLY permanent gaps.
**What:** With H=7 and H=21 both proved as permanent gaps (never achieved for ANY n), the natural question is: are there other permanent gaps?
**Computational evidence (opus-S43):**
  - ALL n=7 gaps (63, 107, 119, 149, 161-169, 173) fill at n=8 (sampling, very quickly)
  - n=8 exhaustive computation running (268M tournaments)
  - For H≥27 (w≥13): 20+ graph-feasible decompositions, blocking all seems impossible
  - Decomposition analysis: for w≥4, the (w-2,1,0,...) decomposition is available; it fails at w=10 due to cascade forcing (Part N) but works at all other w
**Algebraic argument (opus-S43):**
  - For w≥13: alpha_3≥1 decompositions become available (3 disjoint cycles feasible)
  - The number of decompositions grows rapidly with w (14 at w=10, 20 at w=13, 60 at w=20)
  - Each decomposition needs an independent tournament obstruction to block it
  - Only w=3 (1 feasible decomp) and w=10 (4 feasible decomps, all blocked) have this property
**Mod-4 result (Grinberg-Stanley Theorem 7.1):** H(T) ≡ 1 + 2·(# nontrivial odd cycles) mod 4. Does not directly rule out any odd H.
**Conjecture: H=7 and H=21 are the ONLY permanent gaps in the H-spectrum.**
**n=8 EXHAUSTIVE RESULT (opus-S45):** All 268,435,456 tournaments enumerated. Max H=661. Only missing odd values in [1,300]: H=7 and H=21. This CONFIRMS the conjecture through n=8. No new gaps appear. All n=7 gaps (63, 107, 119, etc.) fill at n=8.
**Status:** STRONG EVIDENCE. Conjecture holds through n=8 exhaustive enumeration.

### INV-102: Grinberg-Stanley Mod-4 Theorem (Theorem 7.1) — CATALOGED
**Source:** opus-2026-03-07-S43 (web research)
**Status:** CATALOGED. Read and extracted from arXiv:2307.05569.
**What:** Theorem 7.1: H(T) ≡ 1 + 2·(# nontrivial odd D-cycles) (mod 4). This is the OCF mod-4 reduction: since H = 1 + 2·alpha_1 + 4·(...), H mod 4 = 1 + 2·alpha_1 mod 4. The proof uses the power-sum expansion and the specialization map zeta. Not directly useful for gaps but confirms the algebraic structure.
**Next step:** Check if higher modular refinements (mod 8, mod 16) exist in the Grinberg-Stanley framework.

### INV-103: Non-Separating Vertices in Tournaments — CATALOGED
**Source:** opus-2026-03-07-S43 (web research)
**Status:** CATALOGED. Related but weaker than cycle-rich deletion.
**What:** A vertex in a strongly connected tournament is "non-separating" if its removal preserves strong connectivity. For min in/out-degree ≥ p, at least min{|V|, 4p-2} non-separating vertices exist. Our "good deletion" requirement is stronger: preserve cycle-richness (every vertex in a 3-cycle), not just strong connectivity.
**Next step:** Could the non-separating vertex techniques be adapted to our stronger requirement?

### INV-104: "Cycle-Rich" as Novel Concept — NOTED
**Source:** opus-2026-03-07-S43 (web research)
**Status:** The term "cycle-rich" (every vertex in a directed 3-cycle, no source/sink) does not appear in the literature. This is a novel concept from our project. The poisoning graph argument (Part R) may be publishable as a standalone result about cycle-rich tournaments.

### INV-105: Deletion-Contraction for H(T) — VERIFIED
**Source:** kind-pasteur-2026-03-07-S34
**Status:** VERIFIED COMPUTATIONALLY. H(T) = H(T\e) + H(T/e) holds 100% at n=4,5. Commutative specialization of Mitrovic's W_X = W_{X\e} - W_{X/e}↑ (arXiv:2504.20968).
**Convention:** Contraction merges tail/head: w inherits IN from tail, OUT from head.
**Next step:** Prove algebraically (should follow from Mitrovic by specialization). Use for inductive proof of Redei/OCF/H-gaps.

### INV-106: GLMY Path Homology of Tournaments — IDENTIFIED
**Source:** kind-pasteur-2026-03-07-S34 (web research)
**Status:** IDENTIFIED. Tang-Yau (2026, arXiv:2602.04140) compute path homology of circulant digraphs using Fourier decomposition — directly applicable to Paley tournaments.
**Key connection:** H_1(T) should relate to cycle space of T and thus to Omega(T) and alpha_1.
**Next step:** Compute GLMY path homology for small tournaments (n=4,5). Check if rank(H_1) = alpha_1 or c_3.

### INV-107: Extended Root Polytope Deletion-Contraction — IDENTIFIED
**Source:** kind-pasteur-2026-03-07-S34 (web research). Kalman-Tothmeresz 2024, arXiv:2409.18902.
**Status:** IDENTIFIED. h*-polynomial of extended root polytope is monotone under deletion-contraction.
**Next step:** Compute root polytope for small tournament digraphs. Check if h* relates to I(Omega,x).

### INV-108: Lee-Yang Zeros of F(T,x) — DISCOVERED (opus-S44)
**Source:** opus-2026-03-07-S44
**Status:** MAJOR DISCOVERY. F(T,x) zeros come in reciprocal pairs (palindrome). Cluster at ±2pi/3 on unit circle. H=9 at n=5: ALL zeros on unit circle.
**Key findings:** F(T,omega) real at n=7, F(T,i) pure imaginary, universal divisibilities F(T,omega) = 0 mod 9 and F(T,i) = 0 mod 16i.
**Next step:** Prove Lee-Yang property for specific tournament classes. Connect to phase transitions in statistical mechanics.

### INV-109: Walsh/Fourier Spectral OCF — DISCOVERED (opus-S35c)
**Source:** opus-2026-03-07-S35c9
**Status:** MAJOR DISCOVERY. THM-081: hat{t_k}[S] = (1/2^k) sum (-1)^{asc(S,C)}. Counting identity provides new proof path for OCF.
**Next step:** Prove counting identity algebraically for d=1 (single edge). Extend to n=7 (need hat{t7} and hat{bc35}).

### INV-110: Ihara Zeta Function of Tournaments — TESTED
**Source:** kind-pasteur-2026-03-07-S34
**Status:** TESTED. z_inv(1/2) = det(I-A/2+(D-I)/4) strongly correlated with H (r=-0.95 at n=5) but NOT uniquely determined.
**Conclusion:** Ihara zeta constrains H but doesn't determine it. Consistent with cycles constraining but not determining independence structure.

### INV-111: p-adic Structure Beyond p=2 — TESTED
**Source:** kind-pasteur-2026-03-07-S34
**Status:** TESTED. H mod 3 = (1+2*alpha_1+alpha_2) mod 3 from OCF. At n=4: H mod 3 uniquely determined by c3 via (1+2c3) mod 3. H mod 7 = 0 impossible at n<=6, first achievable at n=7.
**Next step:** Investigate H mod p for larger p. Is there a p-adic tower for p=3?

### INV-112: Converse Invariant Digraph Polynomials — IDENTIFIED
**Source:** kind-pasteur-2026-03-07-S34 (web research). Ai-Gutin-Lei-Yeo-Zhou 2024, arXiv:2407.17051.
**Status:** IDENTIFIED. New digraph polynomial for converse invariance testing. H(T)=H(T^op) means Ham paths are converse invariant.
**Next step:** Read the paper. Check if their polynomial gives new information about H(T).

### INV-113: Stanley-Stembridge Resolution Implications — IDENTIFIED
**Source:** kind-pasteur-2026-03-07-S34 (web research). Hikita 2024.
**Status:** IDENTIFIED. Stanley-Stembridge (e-positivity of chromatic SF for 3+1-free posets) proved by Hikita 2024.
**Connection:** Via Mitrovic-Stojadinovic, Redei-Berge U_T connects to chromatic SF. If tournament poset is (3+1)-free, U_T inherits e-positivity.
**Next step:** Check if tournament arc ordering posets are (3+1)-free. Investigate Hessenberg varieties approach.

### INV-114: Flip Formula F(T,x) - F(T',x) = (x-1)*D(x) — PROVED (THM-083)
**Source:** opus-2026-03-07-S45 (computational discovery), kind-pasteur-2026-03-07-S35 (algebraic proof)
**Status:** PROVED algebraically (THM-083). Verified at n=4,5 exhaustive.
**What:** For arc u->v in T, flip gives T'. The difference F(T)-F(T') factors as (x-1)*D(x), where D = F(T/e) - F(T'/e') is anti-palindromic.
**Key identities (THM-083):**
  - F_T(x) = F_{T\e}(x) + (x-1) * F(T/e, x)  (polynomial deletion-contraction)
  - G_{u,v}(x) = F(T/e, x)  (contraction = conditional path polynomial)
  - D(x) = -x^{n-2} D(1/x)  (anti-palindromicity from tournament palindrome)
  - H(T) - H(T') = D_{n-2}  (leading coefficient of D)
**CORRECTION (kind-pasteur-S35):** H(T) ≠ H(T') under arc flip in general (deltas up to ±12 at n=5). The opus claim "H(T)=H(T')" was WRONG. The correct statement: F(T,1)=n!=F(T',1) (total permutation count).
**CORRECTION (kind-pasteur-S35):** G_uv + G_vu = 2*F(T/e) only when T/e is a tournament (requires u,v to have identical profiles to other vertices). In general, G_uv + G_vu = F(T/e) + F(T'/e') which is palindromic but ≠ 2*F(T/e).
**Scripts:** `04-computation/f_poly_flip_formula.py`, `04-computation/flip_formula_D_analysis.py`, `04-computation/poly_deletion_contraction.py`, `04-computation/flip_reduction_via_contraction.py`

### INV-115: Matroid Structure of Vertex-Disjoint Odd Cycles — BOUNDARY at n=5
**Source:** opus-2026-03-07-S45 (computational discovery)
**Status:** VERIFIED. Exchange axiom holds at n=5 (1024/1024) but FAILS at n>=6 (15360/32768 at n=6).
**What:** Collections of vertex-disjoint odd directed cycles in T form the independent sets of a matroid if and only if n<=5. At n>=6, maximal independent sets can have different sizes.
**Script:** `04-computation/gammoid_matroid_test.py`
**Next step:** (1) Prove n<=5 case (small enough for case analysis). (2) Characterize failure at n=6 — which exchange pairs fail? (3) Relationship to Omega(T) perfectness boundary (perfect through n=7, fails n=8).

### INV-116: Transfer Matrix W(x) and per(W) — EXPLORED
**Source:** opus-2026-03-07-S45
**Status:** EXPLORED. per(W(1)) = D_n (subfactorial) universally. F(T,x) = Hamiltonian path sum over W entries. per(W(x)) palindromic for certain tournament classes.
**Scripts:** `04-computation/transfer_matrix_F_connection.py`, `04-computation/per_W_analysis.py`
**Next step:** (1) Explore eigenvalues of W(x) at specific x values. (2) Connection to Irving-Omar det formula (INV-046). (3) Can det(I-zW) generating function extract F?

### INV-117: Archer-Gessel-Graves-Liang Strong Tournament Descent Polynomial — RESEARCHED
**Source:** opus-2026-03-07-S45 (background agent), Discrete Math 343 (2020)
**Status:** RESEARCHED. Paper fully reviewed. t_n(u) = descent poly for strong tournaments. Palindromic, divisible by (1+u)^{floor(n/2)}. GF: U(x) = 1/(1-T(x)).
**Connection:** Different statistic from F(T,x) (global descent vs path-local forward edges), but same palindromic structure. The (1+u)^{floor(n/2)} divisibility may connect to OCF factor-2 structures. The "Eulerian graphic GF" framework (q=(1+uy)/(1+y)) is the natural algebraic home for descent statistics on tournaments.
**Notes:** `04-computation/gessel_strong_tournament_notes.md`
**Next step:** (1) Compute t_n(u) at small n and compare to F(T,x) aggregates. (2) Check if the strong component decomposition gives new structural insights for H(T). (3) Test the (1+u)^{floor(n/2)} divisibility analogue for F(T,x).

### INV-118: F(T,omega) mod 9 Universality — CONFIRMED NOVEL
**Source:** opus-2026-03-07-S44/S45 (computational discovery + background agent literature search)
**Status:** CONFIRMED NOVEL. Extensive literature search found NO prior work on roots-of-unity evaluations of F(T,x). Not a consequence of Grinberg-Stanley mod-4. Chebikin et al. studied cyclotomic factors of descent set polynomial Q_n but Phi_3 doesn't appear for n<=23.
**What:** F(T,omega) ≡ 0 mod 9 at n=7 (all 5040 tournaments). Equivalently S_0 = sum_{k≡0 mod 3} F_k ≡ 0 mod 6.
**Next step:** (1) Prove algebraically using OCF or Fourier decomposition. (2) Check at n=9,10. (3) Generalize: are there universal congruences for F(T,zeta_k) mod k^2?

### INV-119: Deletion-Contraction for Hamiltonian Paths — PROVED (THM-082)
**Source:** kind-pasteur-2026-03-07-S35
**Status:** PROVED by clean bijection argument. Verified exhaustive n=4,5.
**What:** For any digraph D with directed edge e=(u→v):
  H(D) = H(D\e) + H(D/e)
where D\e = deletion, D/e = contraction (w inherits IN from u, OUT from v).
**Proof:** Ham paths not using e = H(D\e). Ham paths using e biject with Ham paths of D/e via collapsing ...→u→v→... to ...→w→...
**Corollary:** Arc-flip H-difference reduces to contraction: H(T)-H(T') = H(T/e) - H(T'/e').
**Key structural insight:** T/e and T'/e' differ ONLY in how w connects to other vertices (profile swap). If u,v have identical profiles, T/e = T'/e' and H(T) = H(T').
**Connection to Mitrovic:** Commutative specialization of W_X = W_{X\e} - W_{X/e}↑ (arXiv:2504.20968).
**Scripts:** `04-computation/deletion_contraction_test.py`, `04-computation/flip_reduction_via_contraction.py`

### INV-120: Polynomial Deletion-Contraction for F(T,x) — PROVED (THM-083)
**Source:** kind-pasteur-2026-03-07-S35
**Status:** PROVED algebraically. Verified exhaustive n=4,5.
**What:** F_T(x) = F_{T\e}(x) + (x-1) * F(T/e, x). Generalizes THM-082 to polynomial level.
**Key identification:** G_{u,v}(x) = F(T/e, x) — the "conditional path polynomial" summing over permutations with u immediately before v equals the forward-edge polynomial of the contraction.
**Flip formula as corollary:** F_T - F_{T'} = (x-1) * [F(T/e) - F(T'/e')], with D anti-palindromic.
**Anti-palindromicity proof:** D(x) = -x^{n-2}D(1/x) follows from palindromicity of F_T, F_{T'}.
**Scripts:** `04-computation/poly_deletion_contraction.py`

### INV-121: F(T,omega) mod 9 universality — PROVED (THM-085)
**Source:** kind-pasteur-2026-03-07-S36 (extending S35 analysis)
**Status:** PROVED algebraically. Complete proof via Taylor expansion.
**What:** 9 | F(T,omega) for ALL tournaments on n >= 6 vertices. Proof:
1. Taylor expansion F(T,x) = sum c_k (x-1)^k. Over F_3: x^3-1 = (x-1)^3.
2. c_0 = n! (tournament-indep), c_1 = n!(n-1)/2 (tournament-INDEPENDENT!), both divisible by 3.
3. c_2 = A_non + (n-2)!*dp(T), where A_non is tournament-independent and dp(T) = directed 2-path count. Both A_non and (n-2)! divisible by 3 for n >= 5, so c_2 = 0 mod 3 regardless of T.
4. Therefore (x-1)^3 | F(T,x) mod 3 for n >= 5, giving S_r = 0 mod 3.
5. Combined with v_3(n!) >= 2 for n >= 6: 9 | F(T,omega).
**Additional:** Eulerian conjecture: 3|A(n,k) => 3|F_k(T) for all T (verified n=5-8). But this alone doesn't explain n=9,10 (where all A(n,k) = 1 mod 3); the Taylor proof covers all n.
**Sharp:** n=5 has S_r=0 mod 3 but v_3(5!)=1 blocks mod 9. n=4 has c_2 NOT forced 0.
**Scripts:** `04-computation/c2_mod3_proof.py`, `fk_mod3_conjecture.py`, `sr_mod3_n9_check.py`, `f_omega_mod27_analysis.py`

### INV-122: THM-084 naming fix + Corollary 2 error
**Source:** kind-pasteur-2026-03-07-S36
**Status:** FIXED.
**What:** opus-S46 created THM-082-flip-factorization-anti-palindrome.md, colliding with kind-pasteur's THM-082-deletion-contraction-ham-paths.md. Renamed opus's to THM-084. Also fixed Corollary 2 which incorrectly claimed H(T)=H(T') under arc flip (FALSE: H(T) != H(T') in general, deltas up to +-12 at n=5). Correct: F(T,1)=n!=F(T',1) trivially.

### INV-123: Worpitzky Expansion of F(T,x) — PROVED (THM-084)
**Source:** opus-2026-03-07-S46b
**Status:** PROVED algebraically, verified n=3..7
**What:** F(T,x)/(1-x)^n = sum a_m x^m where a_m is polynomial in m of degree n-1.
  - Top 2 coefficients are UNIVERSAL: n and C(n,2)
  - For transitive tournament: a_m = (m+1)^n - m^n (binomial coefficients)
  - Deviation from binomial: delta_2 = 2(n-2)*t3, delta_3 = (n-2)(n-3)*t3
  - At n=6, deeper coefficients need invariants beyond t3
  - Spectral connection: delta_2 = 2(n-2)/3 * tr(A^3)
**Analogy:** F(T,x) is an h*-vector; a_m is an Ehrhart-like polynomial. Transitive tournament corresponds to unit cube h*-vector.
**Scripts:** `04-computation/worpitzky_coefficients.py`, `04-computation/worpitzky_deeper.py`

### INV-124: Signed Forward-Edge Polynomial SF(T,x) — PROVED (THM-085b)
**Source:** opus-2026-03-07-S46b
**Status:** PROVED algebraically
**What:** SF(T,x) = sum sgn(sigma) x^{fwd_T(sigma)} is palindromic with parity (-1)^{C(n,2)}.
  - SF(T,1) = 0 always, so (x-1) | SF(T,x)
  - Quotient SF/(x-1) is anti-palindromic
  - At n=4: SF = c(T) * (x-1)^2(x+1) (since anti-palindromic of even degree has (x+1) factor)
  - SF determines F at n<=5 but NOT at n>=6 (coarser invariant)
**Connection:** SF is a "path immanant" for the sign character. F is the "path permanent."
**Scripts:** `04-computation/signed_F_analysis.py`

### INV-125: Forward-Edge Variance Formula — PROVED (THM-086)
**Source:** opus-2026-03-07-S46b
**Status:** PROVED
**What:** Var[fwd] = (n+1)/12 + 4*t3/(n(n-1)). Exact formula.
  - Non-adjacent forward indicators are UNCORRELATED (tournament completeness: C(n-2,2) edges always)
  - Adjacent covariance = -1/12 + 2*t3/(n(n-1)(n-2))
  - Directed 2-path count = C(n,3) + 2*t3
  - At t3=0 (transitive): Var = (n+1)/12 = Eulerian variance
**Scripts:** `04-computation/worpitzky_w_connection.py`

### INV-126: Cross-Domain Connections for F(T,x)
**Source:** opus-2026-03-07-S46b
**Status:** EXPLORED
**What:** Multiple connections between F(T,x) and other mathematical structures:
  1. **q-analogue F(T,x,q):** q-marginal (inv distribution) is UNIVERSAL = [n]_q! for all T
  2. **det(W(x)):** det((J-I)+(x-1)A) at x=1 = (-1)^{n-1}(n-1) for all T
  3. **Descent algebra:** SF is the sign-character evaluation of the "tournament descent" element
  4. **Worpitzky = Ehrhart:** a_m analogous to Ehrhart polynomial, F to h*-vector
**Next step:** (1) Explore F(T,x,q) as bivariate polynomial. (2) Find polytope whose h*-vector is F(T,x). (3) Connect Worpitzky coefficients to W-polynomial hierarchy (INV-082).

### INV-127: GLMY Path Homology of Tournaments — EVEN BETTI VANISHING
**Source:** opus-2026-03-07-S46e (path_homology_phase2.py output + new analysis)
**Status:** CONJECTURE (strong numerical evidence)
**What:** Tournaments have β_{2k} = 0 for all k >= 1 in GLMY path homology. Only odd Betti numbers can be nonzero. β₁ and β₃ are mutually exclusive. Euler characteristic χ ∈ {0, 1}.
**Evidence:** Exhaustive n=3,4,5 (all 1024+64+8 tournaments); sampled n=6 (500), n=7 (600+). 0 failures.
**Why:** The completeness of tournaments (every pair connected) makes Ω_p ⊂ A_p with ratio ~0.5 at p=2, constraining the chain complex. Transitive tournaments have Ω_p = A_p (ratio 1.0).
**Not in literature:** Checked Tang-Yau (arXiv:2602.04140), GLMY papers — no prior result on tournament even Betti vanishing.
**Connection to OCF:** Both H (via OCF) and χ (via path homology) are expressed as alternating sums over odd-indexed quantities.
**Scripts:** `path_homology_tournament_betti.py`, `even_betti_chain_complex.py`, `even_betti_quick.py`
**Next step:** (1) Prove algebraically for β₂ = 0. (2) Investigate β₁/β₃ mutual exclusivity. (3) Check whether the Ω_p/A_p ratio has a closed form. (4) Try n=8 sampling.

### INV-128: Universal Coefficient Theorem — PROVED (THM-095)
**Source:** opus-2026-03-07-S46e
**Status:** PROVED
**What:** coeff(t_{2k+1} in κ_{2k}) = 2/C(n, 2k) for all k. Proved via forward path formula + OCF + multinomial expansion. Resolves OPEN-Q-023.
**Scripts:** `universal_coeff_proof.py`

### INV-129: Celano-Sieger-Spiro A_T(t) — NOT same as F(T,x)
**Source:** opus-2026-03-07-S46e (web research)
**Status:** CLARIFIED (dead end for direct application)
**What:** arXiv:2309.07240 defines A_T(t) = sum over labelings of t^{des_T(sigma)} where des counts descents across ALL arcs. This has degree C(n,2), not n-1. The (1+t)^{floor(n/2)} divisibility applies to A_T(t), not to our F(T,x). The two polynomials encode different statistics (all-arc descents vs Hamiltonian-path forward edges).
**Impact:** The Celano-Sieger-Spiro result cannot be directly applied to F(T,x). However, it establishes that tournaments have a universal structural constraint on A_T(t) depending only on n, which is analogous to our universal constraint on F(T,x) mod 3 (THM-086).

### INV-130: Pfaffian-Betti Connection — EXHAUSTIVE at n=6
**Source:** opus-2026-03-07-S46e
**Status:** VERIFIED EXHAUSTIVE (32768 tournaments)
**What:** The Pfaffian of the skew-adjacency matrix constrains path homology Betti numbers. At n=6 (exhaustive): β₁>0 ⟹ |Pf(S)| ∈ {1,3}; β₃>0 ⟹ |Pf(S)| ∈ {7,9}. The Pfaffian completely separates β₁ from β₃. Since Pf(S)² = det(S) counts signed perfect matchings, this connects directed topology to algebraic cycle invariants.
**Scripts:** `pfaffian_betti_check.py`, `spectral_betti_connection.py`
**Next step:** (1) Verify at n=7 sampling. (2) Prove the implication β₁>0 ⟹ |Pf| small algebraically. (3) Find exact threshold.

### INV-131: Path Homology Hidden Invariant — Cycle Overlap Pattern
**Source:** opus-2026-03-07-S46e
**Status:** VERIFIED EXHAUSTIVE n=5
**What:** Path homology β₁ is NOT determined by (F-polynomial, t₃, score sequence, strong connectivity). At n=5, the distinguishing invariant is the 3-cycle overlap pattern: β₁=1 iff all 3-cycles form a "star" (share a common edge); β₁=0 when cycle overlaps are heterogeneous. Two non-isomorphic tournament types with identical F-poly and t₃ have different β₁.
**Scripts:** `betti_hidden_invariant.py`, `f_poly_betti_deep.py`
**Next step:** (1) Formalize the overlap pattern criterion. (2) Check at n=6,7.

### INV-132: Mod-2 Cumulant Collapse via OCF
**Source:** opus-2026-03-07-S46e
**Status:** VERIFIED
**What:** THM-094 (F mod 2 universal) implies all integer moment sums n!·E[fwd^r] ≡ 0 mod 2 for n≥3. Root cause: OCF evaluates at x=2 ≡ 0 mod 2, so I(Ω,2) mod 2 = I(Ω,0) = 1, giving H ≡ 1 = Rédei. Combined with THM-086 (mod 3 zeros), all moments vanish mod 6 for n≥3.
**Scripts:** `mod2_cumulant_connection.py`

### INV-133: Tang-Yau Circulant Path Homology (arXiv:2602.04140)
**Source:** opus-2026-03-07-S46e (web research)
**Status:** NEW LEAD
**What:** Tang-Yau (Feb 2026) compute GLMY Betti numbers of circulant digraphs using Fourier decomposition via shift automorphism. Additive structure of connection set S governs low-dimensional chains. Could be adapted to circulant tournaments (odd n, S={1,...,(n-1)/2}).
**Next step:** Read paper. Check if their Fourier method can prove β₂=0 for tournament subclass.

### INV-134: GLMY Torsion-Zeta Bridge
**Source:** opus-2026-03-07-S46e (web research)
**Status:** NEW LEAD
**What:** Grigoryan-Muranov define Reidemeister/analytic torsion on digraphs via path homology Hodge Laplacian, proving Cheeger-Müller analogue. Torsion invariants are classically expressible through zeta-type regularized determinants. This could connect the Ihara/Bartholdi zeta of a tournament to its path homology.
**References:** arXiv:2012.07302, grigor/torsion.pdf
**Next step:** Investigate whether Bartholdi zeta of a tournament encodes path homology info.

### INV-135: Asao — Transitive Tournaments Have Trivial Path Homology
**Source:** opus-2026-03-07-S46e (web research, arXiv:2503.06722)
**Status:** CONFIRMED (known)
**What:** Asao (March 2025) proves transitive tournaments are "regularly diagonal" digraphs with vanishing reduced regular path homology. Confirms our computational finding that transitive tournaments have β=(1,0,...,0).
**Impact:** Provides a published reference for one endpoint of the topological trichotomy.

### INV-136: Chaplin — Random Digraph β₁ Phase Transitions
**Source:** opus-2026-03-07-S46e (web research, arXiv:2111.13493)
**Status:** NEW LEAD
**What:** Chaplin (2022) shows β₁ of random Erdős-Rényi digraphs has two phase transitions. Since tournaments are "density 1/2" digraphs, this places them in a specific regime. Could explain why ~30% of tournaments at n=5 have β₁>0.
**Next step:** Check if their density threshold matches tournament β₁ fraction.
