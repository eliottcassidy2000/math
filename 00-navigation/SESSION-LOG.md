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

## opus-2026-03-05-S2 — 2026-03-05
**Account:** Claude Opus 4.6 (e's MacBook)
**Continuation of:** opus-2026-03-05-S1 (parallel session)
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, CONJ-001, THM-002, THM-003, DISC-001, LaTeX paper (claim_strategies section)
**Summary of work:** Proved THM-008 (mu triviality bound), resolving OPEN-Q-001 (n=5 mystery). Independently verified Claim A at n=4,5 with correct mu implementation. Characterized mu distribution at n=6: mu(3-cycle) in {1,3}, determined by cyclic/transitive structure of available vertices.
**New contributions:**
- THM-008 in 01-canon/theorems/ (mu triviality bound — resolves OPEN-Q-001)
- 04-computation/verify_claim_a.py (independent Claim A verification)
- 04-computation/n5_mu_analysis.py (n=5 mystery investigation)
- 04-computation/n5_mystery_deep.py (exploratory analysis)
- 04-computation/n6_mu_distribution.py (mu distribution at n=6)
- OPEN-Q-001 resolved, OPEN-Q-009 opened and resolved
**Unresolved threads:**
- Claim A proof for n>=6 still open (OPEN-Q-002)
- Need to investigate correlation between mu=3 and per-path failures at n=6
- Proof strategies for Claim A not yet extracted to standalone document

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
