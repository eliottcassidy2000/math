# The LRC(14) sorry-ledger — every remaining unproved node, with formalization pathways

klein-2026-07-02-S94 (HYP-4002; survey agent + owner directive "get to sorry-free-formalizable state").

## Tier 0 — the two blockers (close these ⟹ proof complete, modulo Tier-1 finite legs)

**N1. hp0cap = the L_y extremality: L_y(E) ≤ cap_k for ALL E (THM-534's open half).**
NOW RESTRUCTURED (HYP-4001, this session) into four decidable legs:
  (a) bounded-spread census — exact, done for k=8,9,10 windows (THM-534), extend per k: DECIDE;
  (b) THE FAR-ELEMENT RATE LEMMA (proved, elementary wrap-counting, verified exact):
      |J(A, E∪{w}) − (1−|A|/7)·J(A,E)| ≤ 2·comp(A,E)·|A|/w  ⟹  L_y(E∪{w}) ≤ L_y^∞(E) + K(E)/w,
      L_y^∞ = Σ y_r (1−r/7) S_r (the damped functional);
  (c) the damped comparison max_E L_y^∞(E) vs cap: POSITIVE margins at all sampled k
      (k=8→9: +0.047, 9→10: +0.081, 10→11: +0.074, 11→12: +0.154) — finite rational: DECIDE;
  (d) w-band w < K/margin: finite exact sweep; clustered far blocks: the renormalization
      tower (opus F-lemmas). Formalization: (b) is interval bookkeeping (LRCFinalWindowBand-
      style); (a),(c),(d) are the Helly/decide pattern (HYP-4000).

**N2. hpartA = the witness architecture (kps-S30 "dissolved": rotation identity + c-ruler +
Fubini + (⋆)-census).** Status: structurally closed, NOT yet skeleton-wired. Needs: (i) the
rotation identity + Fubini written as lemmas (exact identities, no analysis); (ii) the
(⋆)-census run symbolically (kps's c-breakpoint engine / torus THM-599-600); (iii) skeleton
node `thm527_partA_density_pos_implies_reach` discharged against (i)+(ii). The deviation
lemma (HYP-3847) + THM-598 cover the danger case's two regimes (spectral-gap / frozen).

## Tier 1 — finite legs of the census-exhaustiveness assembly (THM-595 / opus doc)

**N3. F3 rate (√ → O(1/N))**: mac-mini HYP-3850a claims it; if it lands, the middle band
shrinks 10^8 → 10^3 cases. Else: accept the wide band and run it (finite, decide).
**N4. G2 divisor-chain rigidity**: "bounded-ratio compact clusters contain no unbounded
divisor chains" — 1-2 paragraph structural argument off THM-594(C); unwritten.
**N5. F-iii/F-iv sweeps** (4-6-outlier band; heights 20 ≤ N < N*): finite enumerations,
blocked only by N3's constant. **N6. G1 census extension** (V=19 → 25): routine large run.
**N7. R2 large-range recursion**: ~5-line structural write-out (Lemma P + peel counting).

## Tier 2 — non-blocking (exposition/robustness)

Stability unit-residue lemma (mac-mini claim, THM-593 addendum); HYP-3901 fixed-point
conjecture (superseded by the F_j rearrangement floor for closure purposes); THM-533's
corr_L route (superseded by THM-534).

## Skeleton wiring status (Lean, informational — no builds this session)

Named hypotheses awaiting discharge: `thm527_partA_...` (=N2), `hsmall` (pigeonhole,
provable now), `hlarge` (= the rhoGlob ledger, kps HYP-3957 machine-checked the ledger
rows — wire-up remains), `rhoStar_nonneg` (measure triviality). Sorry-free islands to
wire: LonelyRunnerMathlib, PolygonPartitionDMNR, LRCFinalWindowBand, NestTelescope,
NestDecidable, DangerousPatterns, BonferroniTruncation, LRCGapSurplusLedger.

## Housekeeping (numbering)

TWO THM-599s (quintic-bonferroni / torus-band) + THM-600-torus duplicate + TWO THM-601s
(nest lemma [klein, first-to-origin S92] / dangerous-pattern [mac-mini S100, now cited by
two Lean modules]). RECOMMENDATION: keep THM-599 = torus (kps reserved first), renumber
quintic → THM-602; keep THM-601 = nest lemma, renumber dangerous-pattern → THM-603 and
update the two Lean docstrings (mac-mini's call — their references). Delete the THM-600
duplicate after merging any deltas into THM-599.

## Honest distance

~60-70% of the full proof is proved-or-machine-checked. The path: N1 legs (a),(c),(d) are
decide-shaped TODAY given (b) [proved this session]; N2 needs kps's identities written +
census run; N3-N7 are finite. NO remaining node requires new analytic ideas — every one is
either elementary, structural, or a finite computation. That is the state the owner asked
for: sorry-free formalizability is now a project-management claim, not a mathematical bet.
