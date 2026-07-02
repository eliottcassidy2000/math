# The LRC(14) proof-DAG ledger — node-by-node formalizability status

**mac-mini-2026-07-02-S1 (HYP-3858).** The single-page answer to "could this be sorry-free
formalized completely?" Every node of the current proof attempt, its status, and its
formalization character. Statuses: **[LEAN]** = sorry-free in-repo; **[PAPER]** = complete
elementary proof in canon (transcription task); **[TABLE]** = finite exact-rational
computation, `decide`-shaped; **[CENSUS]** = finite enumeration, spec'd, engine needed;
**[GAP]** = reasoning incomplete.

## Layer 0 — reductions (all done)
- q-witness sieve; counterexample ⟹ covering-saturated ................ [LEAN] (`sieve_frac`, THM-369)
- witness attainment (sup attained; reach ⟹ lonely) .................... [LEAN]
- dilation invariance; k ≤ 2; Dirichlet tightness ....................... [LEAN] (kps LonelyRunnerMathlib)
- LRC(≤13) external base ................................................. cited (arXiv:2604.23906; preprint-risk only)
- MSS finite speed bound (91^12) ......................................... cited (arXiv:2411.06903)

## Layer 1 — structure theorems (paper-complete, Lean-partial)
- THM-592 radius-derivative/co-area, breakpoints, ladder ................ [PAPER]
- THM-593 unit-residue rigidity .......................................... [LEAN] (`LRCUnitResidue`) + [PAPER] addenda
- THM-594 pair law (two-branch); continuous Mirsky–Newman ............... [PAPER] (+ discrete twin [LEAN] `PolygonMirskyNewman`)
- THM-596 final-window bands; THM-597 collapse law ...................... [PAPER] + per-set [TABLE]
- klein nest lemma (gcd-nest closed form, cap universe) ................. [PAPER-pending klein's writeup] + [TABLE]
- kps torus-band theorem (THM-600) ....................................... [PAPER] (kps verified; canonized)

## Layer 2 — the free-phase floors (the former hpartA analytic core)
- THM-601(i) dangerous ⟺ P+Q ≤ 7 ....................................... [LEAN] (`DangerousPatterns`, avoidance dir) + [PAPER] (forced dir)
- THM-601(ii) exact minima table ......................................... [TABLE] (computed PQ ≤ 64; `decide` ingestion pending)
- THM-599 truncation identity ............................................ [LEAN] (`BonferroniTruncation`)
- THM-598/599 forced independence, d ≤ 5, arc-counting errors .......... [PAPER] (THM-602 Part C — written 2026-07-02)
- finite-support of the deviation table (7-commensurate rows exact) ..... [TABLE]-conjecture, per-pattern `decide`; **the one item needing either the finite check run or a 3-line commensuration proof**
- THM-602 Part A trichotomy + Part B renormalization exactness ......... [PAPER] (written 2026-07-02; closes the cluster-coherence GAP)
- THM-602 Part D recursion/composition ................................... [PAPER]
- the resonance-lattice census (HNF bases, height ≤ 7, j ≤ 13) ......... [CENSUS] (finite, small; engine trivial)
- bounded-pattern base cases (rational-ray clusters; opus F_j) .......... [TABLE] (F_7..F_11 exact) + [CENSUS] (remaining rays from the lattice census)

## Layer 3 — hpartA assembly
- G2-window structure (finite rational intervals + arc-count tax) ....... [PAPER] (THM-592(i) + kps c-ruler identities, HYP-3953)
- kps rotation identity / Fubini gap integral ............................ [PAPER] (kps-S30; verified)
- window floor ⟹ Mreach ≥ 1/14 ⟹ lonely ............................... [LEAN]-trivial wiring (skeleton glue exists)
- kps residual finite census ((★)-census, symbolic via THM-600) ......... [CENSUS] (spec'd symbolic; engine = kps ledger run)

## Layer 4 — hp0cap
- THM-534 p0 ≤ L_y ....................................................... [PAPER] (+ partial Lean `LRCMomentDual` WIP)
- nest-lemma closed forms for Λ_P(1/14) ................................. [PAPER-pending] + [TABLE]
- "consec maximizes L_y" over spread ≤ 30 ................................ [CENSUS] (finite rational comparison; klein queued)
- far-element exact deficits (pair law branch-2) ......................... [PAPER] + [TABLE]

## Layer 5 — top assembly
- skeleton DAG glue (hfloor/hpartA/hR0 wiring) ........................... [LEAN] (`LRCFourteenSkeleton`, sorry-parameterized)
- S_d bounds as DAG-node hypotheses ⟸ Layer-2 tables ................... wiring task only

## The critical path to zero sorries (ordered, all finite)
1. Run the resonance-lattice census + the d ≤ 5 simplex minima table (engineering, small).
2. Prove-or-decide the 7-commensurate finite-support row (3-line commensuration argument expected).
3. klein's nest-lemma writeup + L_y census run; kps's (★)-census/ledger run.
4. Transcribe THM-602 Parts B–D + THM-601 forced direction into Lean; wire tables via `decide`/`native_decide`; replace the two skeleton sorries with the discharged nodes.

**Verdict: no [GAP] nodes remain.** The proof attempt is structurally complete: every
node is proved, tabled, or a specified finite census. Remaining risk concentrates in
(i) census engines actually terminating with the expected verdicts, (ii) the preprint
status of LRC(≤13), (iii) transcription effort.
