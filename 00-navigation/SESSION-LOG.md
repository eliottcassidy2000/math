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
