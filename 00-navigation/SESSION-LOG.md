# Session Log

Chronological record of all sessions. Every new Claude instance adds an entry at the **top** of this file before doing any work.

Entry format:
```
## [INSTANCE-ID] — [DATE]
**Account:** [A/B/C/...]
**Continuation of:** [previous instance ID, or "fresh start"]
**Files read:** [list of files read at session start]
**Summary of work:** [brief description]
**New contributions:** [theorem IDs, court cases, tangents added]
**Unresolved threads:** [things left open for next session]
```

## kind-pasteur-2026-03-06-S25f — 2026-03-06 (Grand Synthesis + W(r) Stratification + Pfaffian Duality)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-06-S25e (context limit)
**Files read:** grand-synthesis-S25f.md (prior context), SESSION-LOG.md, INVESTIGATION-BACKLOG.md, opus S27 scripts
**Summary of work:**
  (1) **Grand Synthesis document** (`03-artifacts/drafts/grand-synthesis-S25f.md`): Comprehensive map of all mathematical structures — 5 equivalent algebraic perspectives (independence polynomial, transfer matrix, symmetric functions, Hopf algebra, tiling geometry), 6 novel creative connections (Mobius inversion, perpendicular plane/Mobius strip, Steiner systems, Cayley transform, groups/SC boundary, 1729 mystery).
  (2) **W(r) coefficient stratification** (`04-computation/W_coefficient_stratification.py`): Verified w_0 = -t_3 + 2*t_5 + 1 at n=5 (EXACT for all 11 iso classes). All odd-indexed W coefficients are exactly 0. W(r) stratifies tournament invariants by odd-cycle complexity.
  (3) **W(r)/OCF connection** (`04-computation/W_ocf_connection_v4.py`): At n=5, H = 1 + 2*(total directed odd cycles) since a_2 = 0 (no room for disjoint cycle pairs). Fixed cycle enumeration bug (canonical form filter was wrong).
  (4) **Recursive Hopf structure**: Discovered that overlap=3 contribution to w_{n-5} at n=7 uses OCF at n=5 as subroutine via the Hopf algebra coproduct. This creates a recursive hierarchy: W(r) coefficients at n use OCF at smaller n as building blocks.
  (5) **Pfaffian-path duality** (`04-computation/pfaffian_path_duality.py`): At n=4, det(S) is exactly determined by t_3 (det=9 iff t_3 odd). At n=6, needs finer invariants. Path-cycle duality: odd n has H (paths survive), even n has Pf (cycles survive).
  (6) **Paley tournament eigenvalue analysis**: Fixed Paley construction (requires p ≡ 3 mod 4). Paley T_7 has W(r)/7! = [1/320, 0, 1/80, 0, 1/4, 0, 1]. All non-trivial eigenvalues degenerate → scalar M.
  (7) **Deep connections document** (`03-artifacts/drafts/deep-connections-S25f.md`): Extended analysis of W(r)=tr(M(r)), Cayley transform, Pfaffian structure, BIBD embedding, transfer matrix as random walk, Hopf algebra recursion, 1729 number theory.
  (8) **Integrated opus S27 results**: THM-055 (coefficient hierarchy theorem) connects perfectly — e_{2k}(s_P) is polynomial in f_P, moments of f_P determine W coefficients.
**New contributions:** INV-079 (W(r) stratification), INV-080 (Pfaffian-path duality), INV-081 (Paley W(r) structure), grand-synthesis-S25f.md, deep-connections-S25f.md
**Unresolved threads:**
  - Explicit formula for w_0 at n=7 (depends on finer invariants than (t_3, t_5))
  - Prove Hopf algebra recursion algebraically
  - Does the 4th moment of f_P have a cycle-theoretic interpretation?
  - Compute W(r) for Paley T_11 and T_3
  - Formal path-cycle duality between H and Pf(S)?

## opus-2026-03-06-S11b (continued^4) — 2026-03-06 (Coefficient Hierarchy Theorem)
**Account:** opus
**Continuation of:** opus-2026-03-06-S11b (continued^3)
**Summary of work:**
  (1) **THM-055: Coefficient Hierarchy Theorem** — Proved that tr(c_{n-1-2k}) = sum_P e_{2k}(s_P) where e_{2k} is a polynomial of degree 2k in f_P (forward arc count). Via Newton's identity reduction, all power sums p_j of s_i = ±1/2 collapse to functions of p_1 = f - (n-1)/2.
  (2) **Explicit formulas computed** via sympy for e_0, e_2, e_4, e_6 at n=5,7,9.
  (3) **Moment hierarchy discovered:** sum_P f^j for j≤3 depends only on t_3; sum_P f^4 depends on MORE than t_3. This is why tr(c_{n-5}) cannot be expressed via simple cycle counts.
  (4) **Position decomposition at n=7:** f_(4) splits into 4 types by position adjacency pattern: "4 consec" (3 subsets, 5 vertices = subtournament H sums), "3+1" (6 subsets, 6 vertices), "2+2" (3 subsets, 6 vertices), "2+1+1" (3 subsets, 7 vertices = irreducible whole-tournament invariant).
  (5) **At n=5:** f_(4) = H directly (only 1 subset of 4 positions), giving tr(c_0) = H - 3*t_3.
  (6) Read algebraic proof of c_{n-3} from opus-S27. Read grand synthesis from S25f.
**New contributions:** THM-055, coefficient_hierarchy_proof.py
**Unresolved threads:**
  - Prove sum_P f^3 depends only on t_3 algebraically (currently only verified computationally)
  - What is the explicit invariant that sum_P f^4 captures at n≥7?
  - Can the hierarchy be used to bound or relate c_0 to H?

## kind-pasteur-2026-03-06-S25e — 2026-03-06 (THM-052 DISPROVED for non-SC VT; McKay database)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-06-S25d (context limit, then continuation)
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, INVESTIGATION-BACKLOG.md, frobenius21_palindromic_N.py
**Summary of work:**
  (1) **THM-052 DISPROVED for non-SC VT tournaments:** Computed N(0,1,j) for the F_21 non-normal Cayley tournament at n=21 via prefix+suffix bitmask DP (1075s total). N is NOT palindromic (all 20 values differ from mirror), alternating sum = M[0,1] = 45,478,409 != 0. M is NOT scalar.
  (2) **McKay database validation:** Downloaded McKay's VT tournament data for n=21 (users.cecs.anu.edu.au/~bdm/data/). 88 circulant + 22 non-circulant = 110 total VT tournaments. ALL 88 circulant are self-converse, ALL 22 non-circulant are NOT self-converse. n=21 is the SMALLEST order with non-circulant VT tournaments.
  (3) **digraph6 decoder bug fixed:** Initial decoder skipped diagonal entries in digraph6 format, corrupting adjacency matrices. Fixed by advancing bit index on diagonal. Verified: all circulant tournaments now correctly decode as regular (all out-degree 10).
  (4) **3-cycle counts:** All 22 non-circulant VT tournaments have identical 3-cycle count (385), same as our F_21 construction. All are regular.
  (5) **Literature search:** El Sahili-Ghazo Hanna (2023): T and T^op have same oriented Hamiltonian path type distribution. Ai et al. (2025): converse-invariant digraph polynomial. Neither directly addresses position distributions.
  (6) **MISTAKE-013, MISTAKE-014 logged:** Self-converse assumption false for non-abelian VT; THM-052 scope must be restricted.
**New contributions:** mcKay_vt21_selfconverse.py, MISTAKE-013, MISTAKE-014, INV-077, INV-078
**Unresolved threads:**
  - Which of McKay's 22 non-circulant VT tournaments corresponds to our F_21 construction?
  - What is the exact M matrix structure for the non-SC tournament? (only M[0,1] computed)
  - Verify Aut+Anti characterization at n=7 exhaustively
  - Does any structure theorem hold for non-SC VT M matrices?

## opus-2026-03-06-S26 (continued²) — 2026-03-06 (position-uniform, converse disproof, even-r hierarchy)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-06-S26 (continued) (context limit)
**Summary of work:**
  (1) **THM-052 CONVERSE DISPROVED at n=5:** Scalar M does NOT imply vertex-transitive. Counterexample: "cone over C_3" with scores (1,2,2,2,3), |Aut|=3, but M=3I. Has anti-automorphisms swapping poles.
  (2) **POSITION-UNIFORM <=> SCALAR M at n=3,5 (exhaustive):** N[v,j] = H/n for all v,j is exactly equivalent to M=(H/n)*I. 64/1024 labeled tournaments on 5 vertices satisfy this.
  (3) **Aut+Anti transitivity analysis:** Aut(T) union Anti(T) transitive => scalar M (= VT for tournaments), but NOT conversely. Non-VT scalar-M tournaments have 2 orbits under Aut+Anti.
  (4) **Even-r polynomial for non-VT PU:** VT tournaments have EACH c_k individually scalar. Non-VT PU tournaments have non-scalar c_0 and c_2 that cancel at r=1/2. Key: c_0 off-diagonal = 0.5 cancelled by c_2*0.25 = -0.5.
  (5) **Cone construction does NOT generalize to n=7:** H=135, H%7=2. No non-circulant PU found at n=7 (1000 random + 84 flips tested).
  (6) **FORMULA: tr(c_{n-3}) = 2*(n-2)!*t_3 + const(n).** Verified at n=5 (12*t_3-30) and n=7 (240*t_3-2100). Constants: -0.5, -30, -2100.
  (7) **tr(c_2) at n=7 depends on full iso class,** not just cycle counts. Varies from -57 to +63 within same score sequence.
  (8) **n=9 even-r: top 3 coefficients universal for regular.** tr(c_8)=362880=9!, tr(c_6)=90720, tr(c_4)=6480 all universal. Only tr(c_2), tr(c_0) vary.
  (9) **Even-r vs OCF connection:** H = tr(c_0) + tr(c_2)/4 + ... refines the OCF H = 1 + 2*alpha_1 + 4*alpha_2. At n=5: tr(c_0) = H - 3*t_3.
**New contributions:** THM-052 converse disproof, position-uniform characterization, cone construction analysis, c_{n-3} formula, even-r hierarchy
**Unresolved threads:**
  - Prove tr(c_{n-3}) = 2*(n-2)!*t_3 + const algebraically
  - What is const(n)? Ratios are 60, 70 — pattern unclear
  - At n=7 (prime), is PU = VT? (evidence: no non-VT PU found)
  - What determines tr(c_2) at n>=7? (finer than score sequence)
  - Extend OCF connection: express tr(c_k) in terms of independence polynomial coefficients

## opus-2026-03-06-S11b (continued³) — 2026-03-06 (diagonal signed position + perpendicularity)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-06-S11b (continued²) (context limit)
**Summary of work:**
  (1) **DIAGONAL SIGNED POSITION THEOREM:** M[v,v] = sum_P (-1)^{pos(v,P)} where sum is over all Ham paths. Verified all 12 iso classes at n=5. M[v,v] can be NEGATIVE (not a path count). "Defect vertex" = vertex with position bias far from H/n.
  (2) **c_0 RANK STRUCTURE across H-classes at n=7:**
    - H=189 (Paley): c_0 = 2.25·I (SCALAR, rank 0)
    - H=175: c_0 = 0.25·I (SCALAR, rank 0)
    - H=171: c_0 has rank-2 perturbation, 1 defect vertex with diag=-3.75 vs 0.25
  (3) **PERPENDICULARITY CONFIRMED at n=7 (790 iso classes):**
    Mean cosine of M-directions = -0.0485 (near 0 = perpendicular).
    Low-H classes have positive cosine (aligned), high-H negative (anti-aligned).
    Crossover (true perpendicularity) at H≈95-105 (near median).
  (4) **ALL 43 H-maximizers at n=7 have M = 27·I** (scalar). Supports conjecture: max-H => VT => scalar M.
  (5) **Defect count decreases with H** within score class (1,2,2,2,3) at n=5:
    H=11: 2 defects, H=13: 1 defect, H=15: 0 defects (scalar).
  (6) **M(r) symmetry for ALL r** verified: each c_{2k} individually symmetric (polynomial identity). Already proved by THM-030.
  (7) **MISTAKE-012 corrected:** Blue pair ≠ tournament complement. Complement formula: M(T^c) = diag(M(T)) - offdiag(M(T)).
**New contributions:** diagonal_signed_position_theorem.py, perpendicularity_cosine_n7.py, c0_concentration_theorem.py, rank2_signed_adjacency_n7.py
**Unresolved threads:**
  - Prove diagonal signed position theorem from IE formula
  - Why does perpendicularity occur near median H? Is there a spectral explanation?
  - Extend rank-2 structure to n=9

## opus-2026-03-06-S11b (continued²) — 2026-03-06 (rank-2 signed-adjacency theorem)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-06-S11b (continued) (context limit)
**Summary of work:**
  (1) **RANK-2 SIGNED-ADJACENCY THEOREM at n=7 (H=171):** For ALL 57 regular tilings with H=171, M = 25I + e_v·s_v^T + s_v·e_v^T - 4·e_v·e_v^T where v is the defect vertex (unique fixed point of Aut(T), |Aut|=3) and s_v[w]=±1 is the signed adjacency. Eigenvalues: {23-√10, 25(×5), 23+√10}. Verified 57/57 with zero error.
  (2) **Defect vertex identification:** Fixed point of automorphism group. Has 30 five-cycles vs 25 for normals (invisible at 3-cycle level). The defect is a GLOBAL property — not detectable from local graph statistics (degree, 3-cycles all uniform).
  (3) **Mechanism:** All 3 out-edges from defect have C(v,b,j) = [8,7,8,11,8,7] with alt_sum=-1. Normal vertices have 1 edge to defect with alt_sum=+1, 2 inter-normal edges with alt_sum=0. Total: M[defect]=−3+24=21, M[normal]=1+24=25, difference=4=n−3.
  (4) **M is NOT diagonal** (correcting earlier claim): off-diagonal entries are ±1 in the defect row/column, matching the signed adjacency vector exactly.
**New contributions:** rank2_signed_adjacency_n7.py
**Unresolved threads:**
  - Prove the rank-2 formula algebraically (not just computationally)
  - Check if analogous structure exists at n=9
  - Connect rank-2 perturbation to the polynomial c_0 coefficient

## kind-pasteur-2026-03-06-S25d — 2026-03-06 (THM-052: scalar M PROVED for circulants)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-06-S25c (context limit)
**Summary of work:**
  (1) **THM-052 PROVED: M=(H/n)*I for all circulant tournaments at odd n.** Clean algebraic proof using three ingredients: translation symmetry (N depends on d=b-a mod n), N-symmetry (f(d)=f(n-d)), self-complementarity via sigma: i->-i (f(d,j)=f(n-d,n-2-j)). Combining gives palindromic f(d,j)=f(d,n-2-j), which forces alternating sum=0 at odd n.
  (2) **Exhaustive verification:** Position-uniform => palindromic N confirmed at n=5 (64/64), n=7 (8/8 circulant), n=9 (16/16 circulant), n=11 Paley, n=13 Paley.
  (3) **Position-uniform => self-complementary at n=5:** ALL 64 position-uniform n=5 are SC. Not all are vertex-transitive (40/64 have |Aut|=3, not VT).
  (4) **All circulant tournaments are SC:** Verified n=7 (8/8), n=9 (16/16). The map i->-i always gives self-complementarity for circulants.
  (5) **Processed opus messages:** THM-050 (corrected consecutive formula), THM-051 (reversal identity), eigenvalue formula, F-C decomposition, palindromic landscape analysis.
  (6) **H(T_13) = 1,579,968** computed. f(d,j) constant for Paley T_11 and T_13 (super-palindromic). f(2)=0 for T_13 (vertices 0,2 never adjacent in any path).
**New contributions:** THM-052, palindromic_N_proof.py, palindromic_N_posuniform.py, palindromic_N_n9.py, palindromic_N_n11.py, palindromic_N_proof_attempt.py, selfcomp_posuniform_n7.py, INV-073
**Unresolved threads:**
  - Extend THM-052 to non-circulant vertex-transitive tournaments (need pair-orbit argument)
  - The proof uses circulant-specific translation symmetry — what replaces this for general VT?
  - At n=15, non-circulant VT tournaments exist — these need testing
  - Web search for Babai-Kantor results on VT tournament automorphism groups

## opus-2026-03-06-S26 (continued) — 2026-03-06 (even-r polynomial, THM-052 circulant scalar)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-06-S26 (context limit)
**Summary of work:**
  (1) **Even-r polynomial decomposition M(r) = c_0 + c_2*r^2 + ... + c_{n-1}*r^{n-1}·I:** Proved c_{n-1} = (n-1)!·I universally at odd n. c_4 = 24I at n=5 (verified), c_6 = 720I at n=7 (verified).
  (2) **H = tr(c_0) + 3·(#3-cycles)** at n=5 — verified ALL 12 iso classes. Derived from tr(c_2) = 12·t_3 - 30 and tr(c_4)/16 = 7.5.
  (3) **c_2 eigenvalues determined by score sequence** (not full iso class). Classes within same score group share c_2 spectrum up to permutation. Only c_0 distinguishes iso classes within score groups.
  (4) **THM-052 PROVED: Vertex-transitive tournaments have scalar M = (H/n)·I at odd n.** Proof via reflection-reversal bijection φ giving palindromic N(d,j). Extended from circulant to all vertex-transitive (verified Z/3×Z/3 at n=9).
  (5) **IO reciprocity W(z)·W(-z) = 1 confirmed** independent from M(r) = M(-r) — different symmetries in different variables.
  (6) **Blue pair analysis:** complement does NOT preserve M (path edges fixed). Self-paired classes at n=5: classes 8 and 10.
**New contributions:** THM-052, c_2 spectrum analysis, H formula, even_parity_unification.py, blue_skeleton_even_r_synthesis.py, c2_spectrum_sharing.py, even_r_polynomial_full.py, c4_universal_proof.py, circulant_scalar_m_conjecture.py, circulant_scalar_proof.py, scalar_m_beyond_circulant.py, even_r_n7_circulant.py
**Unresolved threads:**
  - Prove H = tr(c_0) + f(t_3, t_5, ...) at general n
  - Off-diagonal c_2 formula (not just score differences)
  - Does scalar M imply vertex-transitive (converse of THM-052)?
  - Extend even-r polynomial analysis to n=9

## opus-2026-03-06-S11b (continued) — 2026-03-06 (eigenvalue formula, spectral skeleton, perpendicularity)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-06-S11b (context limit)
**Summary of work:**
  (1) **UNIFIED eigenvalue formula for M(T_full_n) at ALL n:** mu_j = 3+2cos(j*pi/K), K=(n+1)/2. Even n: K=half-integer, full ±pairing. Odd n: K=integer, plus lambda_0=1. Verified n=2..25. Spectral radius approaches sqrt(5) as O(pi^2/(2n)).
  (2) **Position variance as perpendicularity mechanism:** pvar(T) has inverted-U shape vs H. At position-uniform tournaments (odd n), C(a,b,j) is constant, so alternating sum 1-1+1-1...=0 forces M off-diag=0. This is WHY position-uniform => M=(H/n)I.
  (3) **SPECTRAL BIPARTITE STRUCTURE discovered at n=5:** H=9 class (18 tilings) splits into TWO spectral sub-classes of 9 each: A={2-sqrt7,1,1,3,2+sqrt7} and B={2-sqrt5,2-sqrt5,1,2+sqrt5,2+sqrt5}. Cross-spectral flip graph is BIPARTITE. Sub-class B has block-diagonal M with eigenvalues phi^3 (golden ratio cube!). Blocks satisfy A+B=4I, AB=-I (mutual negative inverses).
  (4) **n=7 skeleton analysis:** 242 distinct (score,H) classes. Regular score (3,...,3) splits into H=171 (NOT pos-uniform, M has defect 4 at one vertex), H=175 (pos-uniform, M=25I), H=189 (Paley, M=27I). The defect at H=171 is exactly n-3=4.
  (5) **Cross-scale eigenvalue flow:** dH=0 flips preserve lambda=1 eigenvalue. dH>0 compresses eigenvalues toward H/n (scalar). dH<0 spreads them. dH=0 cross-spectral flips rotate within spectral fiber — source of perpendicularity.
**New contributions:** unified_eigenvalue_formula.py, pvar_perpendicularity_mechanism.py, spectral_bipartite_skeleton.py
**Unresolved threads:**
  - Prove perpendicularity analytically (inverted-U of pvar)
  - Does bipartite spectral structure extend to n=7?
  - Prove defect = n-3 for near-uniform regular tournaments
