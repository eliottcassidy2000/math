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
3. OCF H(T)=I(Omega(T),2): VERIFIED n<=6 exhaustive (33,864), n=7 (5000), n=8 (500), n=9 (100), n=10 (30) -- all 0 failures
4. Algebraic analysis: decomposed adj(i,j)-adj'(j,i) as subset convolution. Per-vertex LHS/RHS decomposition does NOT match -- proof must work globally.
5. Key insight: at n<=8, max independent set size in Omega(T) is 2; at n=9, it's 3. Fast I computation exploits this.
**New contributions:**
- tournament_lib.py: 6 new functions (adj_count, flip_arc, compute_s_x, count_directed_5_cycles_through_arc, verify_thm013, verify_ocf, independence_poly_at_fast)
- verify_ocf_sweep.py: clean verification script
- OCF verified through n=10 for first time
- Updated OPEN-QUESTIONS.md with verification records
**Unresolved threads:**
- PROVE the adjacency identity (algebraic approach via subset convolution may be tractable)
- The per-vertex decomposition doesn't match between LHS and RHS -- global approach needed
- Consider transfer matrix / generating function methods

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
