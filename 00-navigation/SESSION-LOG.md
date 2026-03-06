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

## opus-2026-03-05-S4c — 2026-03-05 (Claim B discovery — key reduction for OCF proof)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S4
**Files read:** All warm-up files, previous session's exploration scripts
**Summary of work:** Major progress on OPEN-Q-009. (1) Corrected the s-evenness reduction (C_w+D_w=0 was wrong; correct statement is s-linear coefficients vanish). (2) Verified B(Li,Rj)=B(Lj,Ri) as polynomial identity at n=3,...,7 with real-valued arcs in [-2,3]. (3) Discovered the PAIRING STRUCTURE: for each u!=v, left(u)=right(u) exactly, meaning the d-dependent part of alpha_v vanishes by the INDUCTIVE HYPOTHESIS. (4) Isolated CLAIM (B): a standalone identity about tournament Hamiltonian path weights: sum_{S⊆W\{v}} (-1)^|S| H(S) h_start(W\S,v) = (-1)^{m+1} h_end(W,v). Verified at m=1,...,8. Proved algebraically at m=1,2,3. (5) Established full proof chain: Claim B → B(Li,Rj)=B(Lj,Ri) → OCF → Claim A.
**New contributions:**
- 04-computation/q009_cw_dw_proof.py (showed C_w+D_w!=0, corrected reduction)
- 04-computation/q009_algebra_n4.py (exact n=4 algebra, s-coordinate formula)
- 04-computation/q009_s_quadratic.py (quadratic s-structure analysis)
- 04-computation/q009_direct_proof.py (alpha_v decomposition, pairing discovery)
- 04-computation/q009_claim_b.py (Claim B verification m=1,...,8)
- Updated signed-adjacency-identity.md with proof structure and Claim B
**Unresolved threads:**
- Prove Claim (B) for all m — the ONLY remaining step for OCF proof
- Claim (B) is an alternating-sum identity about internal tournament Hamiltonian paths
- Proved at m=1,2,3 by direct algebra; need general proof (possibly by induction on m)

## kind-pasteur-2026-03-05-S9 — 2026-03-05 (tex deep analysis)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S8 (context continuation)
**Files read:** Full parity_tournaments_fixed.tex (2189 lines), all opus-S4b new files (paper-connections.md, INVESTIGATION-BACKLOG.md, 4 new computation scripts), inbox messages MSG-004 and MSG-005
**Summary of work:** Pulled opus-S4b contributions (transfer matrix symmetry discovery, investigation backlog, paper connections). Then performed line-by-line analysis of entire tex file. Found 5 issues: (1) DR mod-4 proof (Thm 7.4) has broken arithmetic — v_2 analysis produces impossible value, falls back to 3 examples. (2) SE-SYT formula (Thm 7.3) produces non-integer 2^{3/2} for m=2. (3) Transitive uniqueness proof (Prop 2.1) is incomplete/hand-wavy. (4) Verification record outdated (missing n≤8 results). (5) Rajkumar et al. missing from bibliography. Extracted 5 geometric insights, cataloged all 12 references with investigation priority, identified 5 most promising unexploited directions, and mapped paper concepts to recent discoveries. Created comprehensive tex-deep-analysis.md. Added 3 new investigation leads to backlog (INV-028b, INV-029b, INV-030b).
**New contributions:**
- 03-artifacts/drafts/tex-deep-analysis.md (comprehensive analysis report)
- INV-028b: Fix DR mod-4 proof
- INV-029b: Fix SE-SYT formula
- INV-030b: Pin grid S_3 symmetry for OCF
- Investigation backlog updated
**Unresolved threads:**
- Forcade 1973 GF approach (INV-023) — highest priority unexplored lead
- Chapman 2001 ASM bijection (INV-021) — could give determinantal formula for H(T)
- Striker 2011 S_3-equivariance (INV-008) — completely unexplored open problem
- Transfer matrix symmetry proof (INV-001) — opus-S4b's discovery, not yet proved
- Run corrected sympy_proof_n8.py overnight

## opus-2026-03-05-S4b — 2026-03-05 (Signed Position Identity + C verifier)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S8 (context continuation)
**Files read:** Full parity_tournaments_fixed.tex (2189 lines), all opus-S4b new files (paper-connections.md, INVESTIGATION-BACKLOG.md, 4 new computation scripts), inbox messages MSG-004 and MSG-005
**Summary of work:** Pulled opus-S4b contributions (transfer matrix symmetry discovery, investigation backlog, paper connections). Then performed line-by-line analysis of entire tex file. Found 5 issues: (1) DR mod-4 proof (Thm 7.4) has broken arithmetic — v_2 analysis produces impossible value, falls back to 3 examples. (2) SE-SYT formula (Thm 7.3) produces non-integer 2^{3/2} for m=2. (3) Transitive uniqueness proof (Prop 2.1) is incomplete/hand-wavy. (4) Verification record outdated (missing n≤8 results). (5) Rajkumar et al. missing from bibliography. Extracted 5 geometric insights, cataloged all 12 references with investigation priority, identified 5 most promising unexploited directions, and mapped paper concepts to recent discoveries. Created comprehensive tex-deep-analysis.md. Added 3 new investigation leads to backlog (INV-028b, INV-029b, INV-030b).
**New contributions:**
- 03-artifacts/drafts/tex-deep-analysis.md (comprehensive analysis report)
- INV-028b: Fix DR mod-4 proof
- INV-029b: Fix SE-SYT formula
- INV-030b: Pin grid S_3 symmetry for OCF
- Investigation backlog updated
**Unresolved threads:**
- Forcade 1973 GF approach (INV-023) — highest priority unexplored lead
- Chapman 2001 ASM bijection (INV-021) — could give determinantal formula for H(T)
- Striker 2011 S_3-equivariance (INV-008) — completely unexplored open problem
- Transfer matrix symmetry proof (INV-001) — opus-S4b's discovery, not yet proved
- Run corrected sympy_proof_n8.py overnight

---

## opus-2026-03-05-S4b — 2026-03-05 (repo scour + backlog creation)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S4b (context continuation — second part)
**Files read:** Full tex paper (open problems, bibliography, literature), TANGENTS, OPEN-QUESTIONS, MISTAKES, MASTER_FINDINGS, proof-landscape, signed-adjacency-identity, all code file headers, paper-connections.md
**Summary of work:** Completed the user's directive to establish a systematic practice of scouring the repo for leads. Read every significant file in the repo. Created INVESTIGATION-BACKLOG.md with 31 cataloged leads across 5 priority tiers (A: blocks OCF proof, B: structural understanding, C: references to investigate, D: computational targets, E: tangent investigations). Updated CLAUDE.md with mandatory Step 5b (repo scouring every session). Notable findings: Rajkumar paper (2110.05188) is NOT in the tex bibliography despite being connected; 6 bibliography entries have NOT been investigated for OCF connections (Striker, Chapman, Eplett, El Sahili x2, Forcade GF approach); 3 open problems from the paper have zero investigation (mixed graphs, Striker-Chapman equivariance, realizable conflict graphs).
**New contributions:**
- `00-navigation/INVESTIGATION-BACKLOG.md` — comprehensive lead catalog (31 items)
- `CLAUDE.md` Step 5b — mandatory repo-scouring practice for all agents
- Memory file created at `.claude/projects/.../memory/MEMORY.md`
**Unresolved threads:**
- INV-001 (transfer matrix symmetry proof) — highest priority
- INV-008 (Striker-Chapman equivariance) — completely uninvestigated
- INV-015 (add Rajkumar to bibliography)
- INV-024 (compute H(T_19)) — feasible, not done
- n=8 C verifier may still be running from earlier in this session

## kind-pasteur-2026-03-05-S8 — 2026-03-05 (error audit)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S7 (context continuation)
**Files read:** All canon files, all opus-S3/S4 new files, all proof scripts, TANGENTS, OPEN-QUESTIONS
**Summary of work:** Pulled opus-S3 and opus-S4 contributions, then audited the entire codebase for errors. Found and fixed 4 issues: (1) BUG in sympy_proof_n8.py — used simplified n<=7 formula which THM-013 says FAILS at n=8; rewrote to use full A-clique formula. (2) LOGICAL ERROR in even-odd split lemma — claimed "equivalent to OCF" but it's only a consequence; the odd-S sum of Delta(S,R) differs from the cycle formula [g-l]*H(R) at the per-subset level. (3) Stale claim in proof-landscape saying n<=7 instead of n<=8. (4) Duplicate tangent number T040. Verified all SymPy proofs (n=4,5,6) and opus n=7 proof pass correctly. NOTE: opus-S4b ran in parallel and found additional structure (Signed Position Identity, tournament-specificity, bracket analysis) — see their session entry.
**New contributions:**
- MISTAKE-008: even-odd split equivalence claim corrected
- MISTAKE-009: sympy_proof_n8.py formula bug documented and fixed
- sympy_proof_n8.py rewritten with correct full A-clique formula
- Tangent numbering fixed (T040 duplicate -> T044-T047)
- Even-odd-split-lemma.md corrected, OPEN-Q-009 corrected, TANGENTS corrected
**Unresolved threads:**
- Run sympy_proof_n8.py overnight (corrected version, full A-clique formula)
- Prove OCF for all n — central open problem
- Reconcile even-odd split equivalence claim with opus-S4b's Signed Position Identity analysis

---

## opus-2026-03-05-S4b — 2026-03-05 (Transfer matrix symmetry + paper connections)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S3 (parallel session to S4, then extended)
**Files read:** All warm-up files, inbox snippet from S4 diff, all n=8 computation files, arXiv:2110.05188, arXiv:2510.25202
**Summary of work:** Three phases: (1) Verified Even-Odd Split Lemma, proved tournament-specific/polynomial, built C verifier for n=8. (2) Analyzed two external papers: "Tournament Representations" (flip classes, locally transitive, R-cones) and "Dual Burnside Process" (Q=AB factorization, spectral correspondence). (3) MAJOR DISCOVERY: the transfer matrix M[a,b] = sum_S (-1)^|S| E_a(S)*B_b(M\S) is ALWAYS SYMMETRIC (verified 7500+ tests, n=4..8). This is STRONGER than Even-Odd Split — each cross-term individually agrees. Connected to Burnside detailed balance. Also found locally transitive tournaments DO have 5-cycles and 7-cycles (conjecture "only 3-cycles" fails). All OCF tests pass for R-cones, automorphism-symmetric tournaments, locally transitive tournaments.
**New contributions:**
- 03-artifacts/drafts/paper-connections.md (10 connections to external papers)
- 04-computation/locally_transitive_test.py (cycle structure of rank-2 tournaments)
- 04-computation/burnside_connection_test.py (character/Fourier decomposition)
- 04-computation/transfer_matrix_test.py (transfer matrix structure)
- 04-computation/symmetry_check.py (VERIFIES transfer matrix symmetry, 7500+ tests)
- T045 (Transfer matrix symmetry discovery)
- T046 (Tournament representations / flip class connection)
**Unresolved threads:**
- PROVE transfer matrix symmetry for all n (would prove OCF)
- Explore detailed-balance / reversibility interpretation
- R-cone simplification strategy (prove OCF for R-cones, extend via flip)

---

## opus-2026-03-05-S4 — 2026-03-05 (n=7 and n=8 PROVED, even-odd split discovery)
**Account:** Eliott (opus machine)
**Continuation of:** opus-2026-03-05-S3 (context continuation)
**Files read:** TANGENTS.md (resolved merge conflict), THM-015, PROP-001, OPEN-QUESTIONS.md, SESSION-LOG.md, symbolic_proof.py, symbolic_proof_fast.py, kind-pasteur inbox messages
**Summary of work:** Resolved git merge conflict from S3 rebase (renumbered my T032 to T039). Then wrote optimized n=7 exhaustive verifier using numpy-vectorized bitmask approach — PROVED OCF at n=7 (2^20 = 1,048,576 configs in 4 seconds). Extended to n=8 with chunked processing — PROVED OCF at n=8 (2^27 = 134,217,728 configs, 57 minutes). Discovered the Even-Odd Split Lemma: the adj decomposition delta = sum_S Delta(S,R) splits equally between even-|S| and odd-|S| terms, so delta = 2*(odd-S sum). This connects directly to the cycle formula (only odd cycles). The alternating sum sum(-1)^|S| Delta(S,R) = 0 is equivalent to OCF but provides a clean algebraic reformulation. Verified n=5,...,8.
**New contributions:**
- THM-015 updated: n=7 and n=8 PROVED
- 04-computation/q009_prove_n7.py (numpy-vectorized exhaustive verifier)
- 04-computation/q009_prove_n8.py (chunked n=8 verifier)
- 04-computation/q009_even_odd_split.py (even-odd split discovery)
- 04-computation/q009_alternating_sum.py (alternating sum analysis)
- 03-artifacts/drafts/even-odd-split-lemma.md (new proof angle documentation)
- Tangent T040 (even-odd split)
- OPEN-Q-009 updated with proof frontier and even-odd split
**Unresolved threads:**
- n=8 exhaustive verification COMPLETE (134M/134M, 57 min)
- Prove the alternating sum identity for ALL n (equivalent to OCF)
- The even-odd split is a clean algebraic reformulation but doesn't simplify the proof

---

## kind-pasteur-2026-03-05-S7 — 2026-03-05 (n=7 proof + structural analysis)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S6 (context continuation)
**Files read:** SESSION-LOG.md, OPEN-QUESTIONS.md, TANGENTS.md, MISTAKES.md, definitions.md, THM-013, THM-014, THM-015, all S6 code files
**Summary of work:** Extended OCF proof to n=7 via two independent methods: (1) SymPy polynomial identity verification (20 arc variables, 1858 monomials, 77s) and (2) exhaustive {0,1} enumeration (1,048,576 cases, 775s). Both confirm delta_H = delta_I as polynomial identity at n=7. Also explored structural decomposition for general proof: proved delta_H = -sum s_x*(R1+R2) algebraically, showed R1+R2 = 2*H(B_x) at n=4, found excess terms at n=5 satisfy sum s_x*exc(x) = -2*(D5-C5) on {0,1} but NOT as polynomial identity (non-multilinear). Investigated f(S)+f(S^c) pairing: universal at n=4, fails at n>=5 due to 5-cycle corrections. Created complete hand proof for n=4 via f(S) decomposition.
**New contributions:**
- THM-015 updated to PROVED at n<=7 (was n<=6)
- 03-artifacts/code/sympy_proof_n5.py (SymPy polynomial identity proof at n=5)
- 03-artifacts/code/sympy_proof_n6.py (SymPy polynomial identity proof at n=6)
- 03-artifacts/code/sympy_proof_n7.py (SymPy polynomial identity proof at n=7)
- 03-artifacts/code/algebraic_proof_n4.py (complete hand proof at n=4)
- 03-artifacts/code/pairing_proof.py (f(S)+f(S^c) analysis, shows n>=5 failure)
- 03-artifacts/code/proof_structure_analysis.py (pred/succ decomposition, s-signature analysis)
**Unresolved threads:**
- Prove polynomial identity for ALL n — central remaining open problem
- n=8 SymPy verification feasible (user offered overnight compute)
- Structural proof approaches: transfer matrix, induction on n, permanent expansion
- Non-multilinear excess identity blocks naive decomposition approach

---

## opus-2026-03-05-S3 — 2026-03-05
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S2 (context overflow, resumed)
**Files read:** All warm-up + THM-012, THM-013, TANGENTS, OPEN-QUESTIONS
**Summary of work:** Extended tournament_lib.py with arc-flip functions. Independently verified THM-013 adjacency identity. Extended OCF verification to n=10 (first time n>=8 verified).
1. Added adj_count(), flip_arc(), verify_thm013(), verify_ocf(), independence_poly_at_fast() to tournament_lib.py
2. THM-013 adjacency identity: VERIFIED n=5 exhaustive (10240/10240), n=6 sampled (3000/3000)
3. OCF H(T)=I(Omega(T),2): PROVED n<=7 exhaustive (1,048,576 arc assignments at n=7, 27min), n=8 (500 random), n=9 (100), n=10 (30) -- all 0 failures
4. Algebraic analysis: decomposed adj(i,j)-adj'(j,i) as subset convolution. Per-vertex LHS/RHS decomposition does NOT match -- proof must work globally.
5. Key insight: at n<=8, max independent set size in Omega(T) is 2; at n=9, it's 3. Fast I computation exploits this.
6. Explored proof strategies: per-vertex nadj decomposition (dead end), transfer matrix/permanent (no clean formula found), bijection search (suggestive structure but no proof), inclusion-exclusion over blockers (correct framework but doesn't simplify).
7. Key structural finding: s=0 vertices never block swaps (cleanly separates U_T from U_T'), but nadj_sum/H(B_x) ratio is not constant -- the identity is fundamentally global.
**New contributions:**
- tournament_lib.py: 6 new functions (adj_count, flip_arc, compute_s_x, count_directed_5_cycles_through_arc, verify_thm013, verify_ocf, independence_poly_at_fast)
- verify_ocf_sweep.py: clean verification script
- symbolic_proof_n7.py: exhaustive n=7 OCF proof (1,048,576/1,048,576 PASS)
- OCF PROVED through n=7, verified through n=10
- Updated THM-013, THM-015, OPEN-QUESTIONS.md with n=7 exhaustive results
- 4 exploratory scripts: unmatched_decomposition_by_vertex.py, inductive_structure.py, transfer_matrix_approach.py, bijection_search.py, proof_strategy_analysis.py
**Unresolved threads:**
- PROVE the identity for all n (not just n<=7): the central open problem
- Subset convolution framework is the right language but per-vertex doesn't match; need global argument
- The "2-colored independent cycle set" interpretation of I(Omega,2) may lead to a bijective proof
- n=8 proof feasible with optimization (2^27 ~134M cases) but would take hours

---

## kind-pasteur-2026-03-05-S6 — 2026-03-05 (polynomial identity proof)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S5 (context continuation)
**Files read:** SESSION-LOG.md, OPEN-QUESTIONS.md, TANGENTS.md, THM-013, THM-014, all S5 code files
**Summary of work:** Resolved 6-file git merge conflict from S5 rebase. Then discovered and proved a NEW PROOF METHOD for OCF: the swap involution's unmatched path count U_T'-U_T can be expressed as a polynomial identity in arc variables, and this polynomial equals delta_I from THM-013. Proved by hand at n=4 (clean formula: U_T'-U_T = 4-2(p+q+r+t) = 2*sum(s_x)). Verified exhaustively at n=5 (512 cases) and n=6 (16384 cases). This constitutes a PROOF of OCF for n<=6. Found and fixed a sign convention error in THM-013's D/C labels.
**New contributions:**
- THM-015 in 01-canon/theorems/ (swap polynomial identity, proves OCF at n<=6)
- 03-artifacts/code/symbolic_proof.py (polynomial identity verifier for n=4,5,6)
- 03-artifacts/code/symbolic_proof_fast.py (optimized n=6 version)
- 03-artifacts/code/unmatched_decomposition.py (structural analysis of blocking)
- 03-artifacts/code/unmatched_vs_formula.py (aggregate identity tests)
- Tangents T037-T038 (polynomial identity proof, sign convention warning)
- Resolved S5 merge conflicts (TANGENTS renumbered T028-T036, all files clean)
**Unresolved threads:**
- Prove the polynomial identity for ALL n (not just n<=6) — this would prove OCF/Claim A
- The n=4 hand proof suggests a transfer matrix / generating function approach
- n=7 verification feasible (2^19 = 524K cases) but needs optimization
- OPEN-Q-009 partially resolved: proof method identified, works at n<=6

---

## kind-pasteur-2026-03-05-S5 — 2026-03-05 (arc-flip session)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S4
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, CONJ-001, THM-003 (Claim B proof), LaTeX paper claim_strategies section
**Summary of work:** Attacked OPEN-Q-009 (arc-reversal invariance) with new angle. Discovered that E(T) = H(T) - I(Omega(T), 2) is invariant under arc flips (exhaustive n<=5, random n=6). Derived algebraic formula for delta_I via A-clique argument adapted to arc flips. Showed the simple formula delta=2*(gained-lost) works at n<=5 but fails at n>=6. Tested and ruled out contiguous-block decomposition. Formulated PROP-001 (arc-flip identity) as a new proof strategy for OCF/Claim A equivalent to proving a combinatorial identity about Ham path changes under arc flips weighted by complement H-counts.
**New contributions:**
- PROP-001 in 01-canon/theorems/ (arc-flip identity, new proof strategy)
- 03-artifacts/code/arc_reversal_study.py (delta_H decomposition under arc flips)
- 03-artifacts/code/ocf_arc_flip_study.py (E(T) invariance, joint delta_H/delta_I distribution)
- 03-artifacts/code/delta_I_correct_formula.py (algebraic delta_I formula verification)
- 03-artifacts/code/delta_I_fast_test.py, paths_via_cycles_test.py, contiguous_block_test.py
- OPEN-Q-009 reformulated with cleaner E(T) version
- Tangents T028-T031 (arc-flip approach, dead ends)
**Unresolved threads:**
- PROP-001: prove the arc-flip Ham path identity combinatorially
- OPEN-Q-009: the E(T) formulation is cleaner but still needs proof
- Three possible proof approaches: transfer matrix, generating function, involution

---

## kind-pasteur-2026-03-05-S3 — 2026-03-05
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S2
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, DISC-001, CONJ-001, theorem files
**Summary of work:** Rebased kind-pasteur onto main (resolved 6-file conflict). Closed DISC-001 (mu bug does not contaminate verification). Investigated OPEN-Q-010/011 at n=7: wrote test_n7_ABD.py computing A (TypeII sum), B (mu-weighted), D (cycle mu sum) for 1050 pairs. Finding: A=/=D at n=7 in general (only 5.9% exact match). Near-cancellation is statistical (mean A-D~0.1) not algebraic.
**New contributions:**
- 03-artifacts/code/test_n7_ABD.py (correct A-B-D computation at n=7)
- 03-artifacts/code/test_perpath_n7.py (initial wrong test -- documents why naive approach fails)
- DISC-001 moved to resolved/; OPEN-Q-010/011 updated with n=7 findings
- TANGENT T027 (n=7 A-B-D near-cancellation is statistical not algebraic)
**Unresolved threads:**
- OPEN-Q-009: arc-reversal invariance -- the key step, still untouched
- OPEN-Q-012: tower hypothesis -- not yet investigated
- OPEN-Q-013: H(T_p) formula for Paley primes

---

## kind-pasteur-2026-03-05-S2 — 2026-03-05
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S1
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, DISC-001, CONJ-001, CONJ-002, LEM-001, LEM-002, THM-002
**Summary of work:** Resolved OPEN-Q-009 (kind-pasteur numbering) by direct computation. h_QR=h_NQR=201, c_9(T_11)=11055. Direct enumeration gives H(T_11)=95095=55*1729, refuting CONJ-002 for p=11. Added finish_session.py enforcement tooling.
**New contributions:**
- H(T_11)=95095 computed directly (03-artifacts/code/compute_H_T11.py)
- CONJ-002 REFUTED for p=11; OPEN-Q-013 opened (correct formula for H(T_p))
- MISTAKE-006: ratio-coincidence c_k/C(11,k) has no basis
- TANGENTS T022-T026 (Paley structure, 1729, symmetry)
- agents/finish_session.py and agents/check_session_closed.py (end-of-session enforcement)
**Unresolved threads:**
- OPEN-Q-013: What is the correct formula for H(T_p)?
- OPEN-Q-009 (opus): arc-reversal invariance — the key proof step
- DISC-001: needs formal close

---

## kind-pasteur-2026-03-05-S1 — 2026-03-05
**Account:** Eliott (primary)
**Continuation of:** fresh start — first session on this machine
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, inbox/PROCESSING-REPORT.md, agents/REGISTRY.md, all theorem files, DISC-001
**Summary of work:** Registered as new agent. Processed PALEY_T11_c9_ANALYSIS.md inbox doc: extracted LEM-001, LEM-002, CONJ-002. Argued Position A in DISC-001.
**New contributions:** LEM-001, LEM-002, CONJ-002, MISTAKE-006 (ratio-coincidence), DISC-001 Letter 2
**Unresolved threads:** DISC-001 response needed; h_QR/h_NQR computation

---

## kind-pasteur-2026-03-05-S4 — 2026-03-05
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S3
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, opus broadcast messages
**Summary of work:** Structural fix session. kind-pasteur was pushing to claude/kind-pasteur branch (invisible to opus and CI). Merged claude/kind-pasteur into main. Fixed Windows encoding issues in processor.py (CP1252 fails on Unicode). Updated finish_session.py to push HEAD:main always.
**New contributions:**
- All kind-pasteur S1-S3 research now on main (Paley computation, n=7 A-B-D, DISC-001 resolved)
- agents/processor.py: Windows UTF-8 fix (stdout.reconfigure + write_text encoding)
- agents/finish_session.py: push to origin HEAD:main (not current branch)
**Unresolved threads:**
- OPEN-Q-009 (arc-reversal invariance): top priority, untouched
- OPEN-Q-012 (tower hypothesis): untouched

---

## opus-2026-03-05-S2 — 2026-03-05
**Account:** Claude Opus 4.6 (e's MacBook)
**Continuation of:** opus-2026-03-05-S1 (parallel session)
**Files read:** All navigation, canon, court files; LaTeX paper; FINAL_FINDINGS.md; file.txt; tournament_lib.py
**Summary of work:** Proved THM-008b (general mu triviality), THM-012 (partial mu invariance under arc flips). Disproved MISTAKE-004 (OCF IS a valid closed form — verified H(T)=I(Omega(T),2) for all n<=6). Resolved OPEN-Q-011 (near-cancellation is statistical, not structural). Verified adjacency formula H(T)-H(T')=adj(i,j)-adj'(j,i). Established arc-reversal decomposition framework for Q-009.
**New contributions:**
- THM-008b: general mu triviality bound L >= n-2
- THM-012: mu invariant under arc flips when at least one endpoint in V(C)\{v}
- DISC-002: MISTAKE-004 is wrong — OCF is a valid closed form (verified 33,864 tournaments)
- OPEN-Q-011 resolved (near-cancellation is statistical artifact)
- Adjacency formula verified; arc-reversal decomposition framework built
- 6 computation scripts in 04-computation/
**Unresolved threads:**
- DISC-002 needs formal resolution (retract MISTAKE-004)
- Arc-reversal proof strategy (Q-009) partially developed — key obstacle is tracking cycle creation/destruction and mu changes simultaneously
- Claim A proof for general n remains open (OPEN-Q-002)

---

## opus-2026-03-05-S1 --- 2026-03-05
**Account:** opus (Claude Opus 4.6)
**Continuation of:** SYSTEM-2026-03-05-S1 (initial setup)
**Files read:** All navigation files, all canon files, DISC-001, agents/processor.py, inbox/processor.py, README.md
**Summary of work:** Built the complete tournament computation library (tournament_lib.py) implementing all core objects from definitions.md. Independently verified Claim A at n<=6 (196,608 pairs, 0 failures) with a MISTAKE-001-compliant implementation. Also verified Claim B at n<=5, Redei at n<=6, and Claim A at n=7 by random sampling (100 tournaments, 0 failures). Ingested two contributed files (FINAL_FINDINGS.md, file.txt) containing substantial new results. Extracted 4 new theorems (THM-008 through THM-011), 2 new mistakes (MISTAKE-004, MISTAKE-005), resolved OPEN-Q-001 and OPEN-Q-003, added 4 new open questions (Q-009 through Q-012), and documented 3 dead ends (T016-T018).
**New contributions:**
- 03-artifacts/code/tournament_lib.py (CODE-000): complete computation library
- 03-artifacts/code/verify.py (CODE-000v): CLI verification runner
- THM-008: mu=1 trivially for n<=5 (resolves OPEN-Q-001)
- THM-009: per-path failure characterization at n=6 (resolves OPEN-Q-003)
- THM-010: n=4 block-counting theorem (exactly 3 Ham paths)
- THM-011: general block-counting formula H_C^+(T)
- MISTAKE-004: OCF is recursive, not closed-form over all cycles
- MISTAKE-005: cycle bijection under arc reversal fails
- OPEN-Q-009: arc-reversal invariance (key unproved step)
- OPEN-Q-010 through Q-012: new questions from ingested files
- Dead ends T016-T018 documented
- Registered agent: opus
**Unresolved threads:**
- Claim A proof for general n (CONJ-001) -- arc-reversal invariance (OPEN-Q-009) is the key step
- DISC-001 can be moved to resolved
- The near-cancellation phenomenon (OPEN-Q-011) is a promising proof strategy lead
- Per-path formula with 3+5 cycles at n=7 (OPEN-Q-010) needs computational testing
- Computation library could be extended: inshat, per-path identity, arc-reversal testing

---

## SYSTEM-2026-03-05-S1 — 2026-03-05
**Account:** System (initial setup)
**Continuation of:** fresh start — first Cowork session
**Files read:** MASTER_FINDINGS.md (uploaded), parity_tournaments_fixed.tex (uploaded)
**Summary of work:** Built the full research directory system. Processed two uploaded source files as the first contributions. Extracted theorems F1–F5 and Claim A to canon. Logged MISTAKE-001. Opened DISC-001 as a potential court case seed (μ computation bug vs. paper's 0-failure verification claim).
**New contributions:**
- THM-001 through THM-007, CONJ-001 in 01-canon/theorems/
- MISTAKES.md: MISTAKE-001 (μ computation bug in scripts 6-9)
- 00-navigation/TANGENTS.md, OPEN-QUESTIONS.md populated from source files
- 03-artifacts/drafts/parity_tournaments_fixed.tex archived
- 02-court/active/DISC-001-mu-bug-vs-verification.md opened
**Unresolved threads:**
- Claim A proof (the central open problem)
- Why does per-path identity hold at n=5 despite 5-cycles? (OPEN-Q-001)
- μ computation bug in scripts 6-9 needs investigation (DISC-001)
- Formal proof of Fix(σ) = 2^{m²} for self-evacuating SYT (OPEN-Q-007)
